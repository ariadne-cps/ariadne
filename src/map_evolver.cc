/***************************************************************************
 *            map_evolver.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
#include "config.h"

#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "vector.h"
#include "function.h"
#include "enclosure.h"
#include "orbit.h"

#include "logging.h"

#include "map.h"
#include "map_evolver.h"

namespace {

using namespace Ariadne;

template<class ES> List<ES> subdivide(const ES& enclosure) {
    List<ES> result;
    Pair<ES,ES> split=enclosure.split();
    result.append(split.first);
    result.append(split.second);
    return result;
}

} // namespace


namespace Ariadne {

// Allow subdivisions in upper evolution
const bool ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution
const bool ENABLE_PREMATURE_TERMINATION = false;

static const int BLOCKING_EVENT = -2;
using std::shared_ptr;

class DegenerateCrossingException { };


MapEvolver::MapEvolver(const SystemType& system)
    : _sys_ptr(system.clone())
    , _configuration(new ConfigurationType())
{
}


Orbit<MapEvolver::EnclosureType>
MapEvolver::
orbit(const EnclosureType& initial_set,
      const TerminationType& termination,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     initial_set,termination,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}




enum PredicateKind { INVARIANT, ACTIVATION, GUARD, TIME, MIXED };
enum CrossingKind { TRANSVERSE, TOUCHING, NONE, UNKNOWN };



void
MapEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const EnclosureType& initial_set,
           const TerminationType& maximum_time,
           Semantics semantics,
           bool reach) const
{
    verbosity=0;

    typedef RealVectorFunction FunctionType;

    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    List< TimedEnclosureType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
        EnclosureType initial_enclosure(initial_set);
        ARIADNE_LOG(6,"initial_enclosure = "<<initial_enclosure<<"\n");
        TimeType initial_time = 0;
        ARIADNE_LOG(6,"initial_time = "<<initial_time<<"\n");
        working_sets.push_back(make_pair(initial_time,initial_enclosure));
    }


    while(!working_sets.empty()) {
        TimedEnclosureType current_set=working_sets.back();
        working_sets.pop_back();
        EnclosureType initial_enclosure=current_set.second;
        TimeType initial_time=current_set.first;
        Float initial_set_radius=radius(initial_enclosure.bounding_box());
        if(initial_time>=maximum_time) {
            final_sets.adjoin(EnclosureType(initial_enclosure));
        } else if(UPPER_SEMANTICS && ENABLE_SUBDIVISIONS
                  && (initial_set_radius>this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List<EnclosureType> subdivisions=subdivide(initial_enclosure);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_enclosure=subdivisions[i];
                working_sets.push_back(make_pair(initial_time,subdivided_enclosure));
            }
        } else if(LOWER_SEMANTICS && ENABLE_PREMATURE_TERMINATION && initial_set_radius>this->_configuration->maximum_enclosure_radius()) {
            ARIADNE_WARN("Terminating lower evolution at time " << initial_time
                         << " and set " << initial_enclosure << " due to maximum radius being exceeded.");
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  current_set,maximum_time,
                                  semantics,reach);
        }
    }

}





void
MapEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& current_set,
                const TimeType& maximum_time,
                Semantics semantics,
                bool reach) const
{
    EnclosureType initial_enclosure;
    TimeType initial_time;

    ARIADNE_LOG(9,"working_set = "<<current_set<<"\n");
    make_lpair(initial_time,initial_enclosure)=current_set;

    ARIADNE_LOG(6,"initial_time = "<<initial_time<<"");
    ARIADNE_LOG(6,"initial_enclosure = "<<initial_enclosure<<"\n");

    ARIADNE_LOG(2,"box = "<<initial_enclosure.bounding_box()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(initial_enclosure.bounding_box())<<"\n\n");
    //const uint nd=initial_enclosure.result_size();
    //const uint ng=initial_enclosure.argument_size();


    /////////////// Main Evolution ////////////////////////////////


    // Compute the map model
    EnclosureType final_enclosure=Ariadne::apply(_sys_ptr->function(),initial_enclosure);
    TimeType final_time=initial_time+1;
    ARIADNE_LOG(6,"final_enclosure = "<<final_enclosure<<"\n");
    ARIADNE_LOG(4,"Done computing evolution\n");

    reach_sets.adjoin(final_enclosure);

    intermediate_sets.adjoin(final_enclosure);
    working_sets.push_back(make_pair(final_time,final_enclosure));

    ARIADNE_LOG(2,"Done evolution_step.\n\n");

}


MapEvolverConfiguration::MapEvolverConfiguration()
{
    maximum_enclosure_radius(100.0);
}


std::ostream&
MapEvolverConfiguration::write(std::ostream& os) const
{
    os << "MapEvolverSettings"
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << "\n)\n";
    return os;
}


}  // namespace Ariadne

