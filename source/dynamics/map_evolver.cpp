/***************************************************************************
 *            dynamics/map_evolver.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../utility/array.hpp"
#include "../utility/tuple.hpp"
#include "../utility/stlio.hpp"
#include "../algebra/vector.hpp"
#include "../function/function.hpp"
#include "../function/constraint.hpp"
#include "../dynamics/enclosure.hpp"
#include "../dynamics/orbit.hpp"

#include "../output/logging.hpp"

#include "../dynamics/map.hpp"
#include "../dynamics/map_evolver.hpp"

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

// The maximum allowable radius of a basic set.
MapEvolverConfiguration::RealType DEFAULT_MAXIMUM_ENCLOSURE_RADIUS = 100;
// Allow subdivisions in upper evolution.
const Bool DEFAULT_ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution.
const Bool DEFAULT_ENABLE_PREMATURE_TERMINATION = false;

using std::shared_ptr;

FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory();

class DegenerateCrossingException { };


MapEvolver::MapEvolver(const SystemType& system)
    : _sys_ptr(system.clone())
    , _configuration(new ConfigurationType())
{
    this->charcode = "m";
}

typename MapEvolver::FunctionFactoryType const MapEvolver::function_factory() const {
    return ValidatedFunctionModelDPFactory(make_taylor_function_factory());
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const ExactBoxType& box) const {
    return Enclosure(box,this->function_factory());
}

typename MapEvolver::EnclosureType MapEvolver::enclosure(const RealBox& box) const {
    return Enclosure(box,this->function_factory());
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

Void
MapEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const EnclosureType& initial_set,
           const TerminationType& maximum_time,
           Semantics semantics,
           Bool reach) const
{
    verbosity=0;

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
        FloatDPUpperBound initial_set_radius=initial_enclosure.bounding_box().radius();
        if(initial_time>=maximum_time) {
            final_sets.adjoin(EnclosureType(initial_enclosure));
        } else if(semantics == Semantics::UPPER && this->_configuration->enable_subdivisions()
                  && decide(initial_set_radius>this->_configuration->maximum_enclosure_radius())) {
            // Subdivide
            List<EnclosureType> subdivisions=subdivide(initial_enclosure);
            for(Nat i=0; i!=subdivisions.size(); ++i) {
                EnclosureType const& subdivided_enclosure=subdivisions[i];
                working_sets.push_back(make_pair(initial_time,subdivided_enclosure));
            }
        } else if(semantics == Semantics::LOWER && this->_configuration->enable_premature_termination() && decide(initial_set_radius>this->_configuration->maximum_enclosure_radius())) {
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





Void
MapEvolver::
_evolution_step(List< TimedEnclosureType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const TimedEnclosureType& current_set,
                const TimeType& maximum_time,
                Semantics semantics,
                Bool reach) const
{
    EnclosureType initial_enclosure;
    TimeType initial_time;

    ARIADNE_LOG(9,"working_set = "<<current_set<<"\n");
    make_lpair(initial_time,initial_enclosure)=current_set;

    ARIADNE_LOG(6,"initial_time = "<<initial_time<<"");
    ARIADNE_LOG(6,"initial_enclosure = "<<initial_enclosure<<"\n");

    ARIADNE_LOG(2,"box = "<<initial_enclosure.bounding_box()<<" ");
    ARIADNE_LOG(2,"radius = "<<initial_enclosure.bounding_box().radius()<<"\n\n");
    //const Nat nd=initial_enclosure.result_size();
    //const Nat ng=initial_enclosure.argument_size();


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
    set_maximum_enclosure_radius(DEFAULT_MAXIMUM_ENCLOSURE_RADIUS);
    set_enable_subdivisions(DEFAULT_ENABLE_SUBDIVISIONS);
    set_enable_premature_termination(DEFAULT_ENABLE_PREMATURE_TERMINATION);
}


OutputStream&
MapEvolverConfiguration::_write(OutputStream& os) const
{
    os << "MapEvolverSettings"
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << "\n)\n";
    return os;
}


}  // namespace Ariadne

