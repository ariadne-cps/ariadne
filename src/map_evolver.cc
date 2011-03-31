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

#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "vector.h"
#include "function.h"
#include "taylor_model.h"
#include "taylor_function.h"
#include "taylor_set.h"
#include "orbit.h"
#include "taylor_calculus.h"
#include "settings.h"

#include "logging.h"

#include "map.h"
#include "map_evolver.h"

namespace {

using namespace Ariadne;

template<class V, class Iter> inline
void append(V& v, Iter begin, Iter end)
{
    for(;begin!=end;++begin) {
        v.push_back(*begin);
    }
}

template<class V, class C> inline
void append(V& v, const C& c)
{
    for(typename C::const_iterator iter=c.begin(); iter!=c.end(); ++iter) {
        v.push_back(*iter);
    }
}


}



namespace Ariadne {

// Allow subdivisions in upper evolution
const bool ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution
const bool ENABLE_PREMATURE_TERMINATION = false;

static const int BLOCKING_EVENT = -2;
using boost::shared_ptr;

class DegenerateCrossingException { };


MapEvolver::MapEvolver()
    : _settings(new EvolutionSettingsType()),
      _toolbox(new TaylorCalculus())
{
}



MapEvolver::MapEvolver(const EvolutionSettingsType& p)
    : _settings(new EvolutionSettingsType(p)),
      _toolbox(new TaylorCalculus())
{
}


Orbit<MapEvolver::EnclosureType>
MapEvolver::
orbit(const SystemType& system,
      const EnclosureType& initial_set,
      const TimeType& time,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics);
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
           const SystemType& system,
           const EnclosureType& initial_set,
           const TimeType& maximum_time,
           Semantics semantics) const
{
    verbosity=0;

    typedef VectorFunction FunctionType;
    typedef EnclosureType SetModelType;

    typedef boost::shared_ptr< const VectorFunction > FunctionConstPointer;

    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    typedef tuple<TimeType, SetModelType> TimedSetType;
    typedef Float RealType;

    std::vector< TimedSetType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
        SetModelType initial_set_model=this->_toolbox->set_model(initial_set);
        ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
        TimeType initial_time = 0;
        ARIADNE_LOG(6,"initial_time = "<<initial_time<<"\n");
        working_sets.push_back(make_tuple(initial_time,initial_set_model));

		// Checks for match between the enclosure cell size and the set size
		ARIADNE_ASSERT_MSG(this->_settings->maximum_enclosure_cell.size() == initial_set_model.size(), "Error: mismatch between the maximum_enclosure_cell size and the set size.");
    }

    while(!working_sets.empty()) {
        TimedSetType current_set=working_sets.back();
        working_sets.pop_back();
        SetModelType initial_set_model=current_set.second;
        TimeType initial_time=current_set.first;
		Vector<Interval> initial_set_model_range=initial_set_model.range();

		// Checks whether the range can be included into the maximum_enclosure_cell
		bool maximum_enclosure_reached = false;
		for (uint i=0;i<initial_set_model_range.size();++i) {
			if (initial_set_model_range[i].width() > this->_settings->maximum_enclosure_cell[i]) {
				maximum_enclosure_reached = true;
				break; 
			}
		}

        if(initial_time>=maximum_time) {
            final_sets.adjoin(this->_toolbox->enclosure(initial_set_model));
        } else if(UPPER_SEMANTICS && ENABLE_SUBDIVISIONS && maximum_enclosure_reached) {
            // Subdivide
            array< SetModelType > subdivisions=this->_toolbox->subdivide(initial_set_model);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                SetModelType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_tuple(initial_time,subdivided_set_model));
            }
        } else if(LOWER_SEMANTICS && ENABLE_PREMATURE_TERMINATION && maximum_enclosure_reached) {
            std::cerr << "WARNING: Terminating lower evolution at time " << initial_time
                      << " and set " << initial_set_model << " due to maximum radius being exceeded.";
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  system,current_set,maximum_time,
                                  semantics);
        }
    }

}





void
MapEvolver::
_evolution_step(std::vector< TimedSetType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const SystemType& system,
                const TimedSetType& current_set,
                const TimeType& maximum_time,
                Semantics semantics) const
{
    SetModelType initial_set_model;
    TimeType initial_time;

    ARIADNE_LOG(9,"working_set = "<<current_set<<"\n");
    make_ltuple(initial_time,initial_set_model)=current_set;

    ARIADNE_LOG(6,"initial_time = "<<initial_time<<"");
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");

    ARIADNE_LOG(2,"box = "<<initial_set_model.range()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(initial_set_model.range())<<"\n\n");
    //const uint nd=initial_set_model.result_size();
    //const uint ng=initial_set_model.argument_size();


    /////////////// Main Evolution ////////////////////////////////


    // Compute the map model
    SetModelType final_set_model=this->_toolbox->reset_step(system.function(),initial_set_model);
    TimeType final_time=initial_time+1;
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
    ARIADNE_LOG(4,"Done computing evolution\n");

    reach_sets.adjoin(final_set_model);

    intermediate_sets.adjoin(final_set_model);
    working_sets.push_back(make_tuple(final_time,final_set_model));

    ARIADNE_LOG(2,"Done evolution_step.\n\n");

}




}  // namespace Ariadne

