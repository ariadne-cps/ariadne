/***************************************************************************
 *            vector_field_evolver.cc
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
#include "expression_interface.h"
#include "function_interface.h"
#include "taylor_model.h"
#include "taylor_variable.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "orbit.h"
#include "taylor_calculus.h"
#include "evolution_parameters.h"

#include "logging.h"

#include "vector_field.h"
#include "vector_field_evolver.h"

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


VectorFieldEvolver::VectorFieldEvolver()
    : _parameters(new EvolutionParametersType()),
      _toolbox(new TaylorCalculus())
{
}



VectorFieldEvolver::VectorFieldEvolver(const EvolutionParametersType& p)
    : _parameters(new EvolutionParametersType(p)),
      _toolbox(new TaylorCalculus())
{
}


Orbit<VectorFieldEvolver::EnclosureType> 
VectorFieldEvolver::
orbit(const SystemType& system, 
      const EnclosureType& initial_set, 
      const TimeType& time,
      Semantics semantics) const 
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final; 
    EnclosureListType reachable; 
    EnclosureListType intermediate; 
    this->_evolution(final,reachable,intermediate,
                     system,initial_set,time,semantics,false); 
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}




enum PredicateKind { INVARIANT, ACTIVATION, GUARD, TIME, MIXED };
enum CrossingKind { TRANSVERSE, TOUCHING, NONE, UNKNOWN };



void
VectorFieldEvolver::
_evolution(EnclosureListType& final_sets, 
           EnclosureListType& reach_sets, 
           EnclosureListType& intermediate_sets, 
           const SystemType& system, 
           const EnclosureType& initial_set, 
           const TimeType& maximum_time, 
           Semantics semantics, 
           bool reach) const
{

  
    typedef FunctionInterface FunctionType;
    typedef Vector<Interval> BoxType;
    typedef TaylorFunction FunctionModelType;
    typedef TaylorFunction FlowModelType;
    typedef TaylorSet SetModelType;


    typedef boost::shared_ptr< const FunctionInterface > FunctionConstPointer;

    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    typedef tuple<TimeType,SetModelType> TimedSetType;
    typedef Float RealType;

    std::vector< TimedSetType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
        TimeType initial_time = 0.0;
        ARIADNE_LOG(6,"initial_time = "<<initial_time<<"\n");
        SetModelType initial_set_model=this->_toolbox->set_model(initial_set);
        ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
        working_sets.push_back(make_tuple(initial_time,initial_set_model));
    }


    while(!working_sets.empty()) {
        TimedSetType current_timed_set=working_sets.back();
        working_sets.pop_back();
        TimeType current_time=current_timed_set.first;
        SetModelType current_set_model=current_timed_set.second;
        RealType current_set_radius=radius(current_set_model.range());
        if(current_time>=maximum_time) {
            final_sets.adjoin(this->_toolbox->enclosure(current_set_model));
        } else if(UPPER_SEMANTICS && ENABLE_SUBDIVISIONS
                  && (current_set_radius>this->_parameters->maximum_enclosure_radius)) {
            // Subdivide
            array< SetModelType > subdivisions=this->_toolbox->subdivide(current_set_model);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                SetModelType const& subdivided_set_model=subdivisions[i];
                working_sets.push_back(make_tuple(current_time,subdivided_set_model));
            }
        } else if(LOWER_SEMANTICS && ENABLE_PREMATURE_TERMINATION && current_set_radius>this->_parameters->maximum_enclosure_radius) {
            std::cerr << "WARNING: Terminating lower evolution at time " << current_time
                      << " and set " << current_set_model << " due to maximum radius being exceeded.";
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  system,current_timed_set,maximum_time,
                                  semantics,reach);
        }
    }

}





void
VectorFieldEvolver::
_evolution_step(std::vector< TimedSetType >& working_sets,
                EnclosureListType& final_sets, 
                EnclosureListType& reach_sets, 
                EnclosureListType& intermediate_sets, 
                const SystemType& system, 
                const TimedSetType& working_timed_set_model,
                const TimeType& maximum_time, 
                Semantics semantics, 
                bool reach) const
{
    typedef FunctionInterface FunctionType;
    typedef Vector<Interval> BoxType;
    typedef TaylorFunction MapModelType;
    typedef TaylorFunction FlowModelType;
    typedef TaylorSet SetModelType;
  
    SetModelType current_set_model;
    TimeType current_time;
    ARIADNE_LOG(9,"working_timed_set_model = "<<working_timed_set_model<<"\n");
    make_ltuple(current_time,current_set_model)=working_timed_set_model;

    ARIADNE_LOG(6,"current_time = "<<current_time<<"");
    ARIADNE_LOG(6,"current_set_model = "<<current_set_model<<"\n");
  
    ARIADNE_LOG(2,"box = "<<current_set_model.range()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(current_set_model.range())<<"\n\n");
    //const uint nd=initial_set_model.result_size();
    //const uint ng=initial_set_model.argument_size();
  

    /////////////// Main Evolution ////////////////////////////////
    const FunctionType& dynamic=system.function();
  
    // Set evolution parameters
    const Float maximum_step_size=this->_parameters->maximum_step_size;
    const Float maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
    const Float zero_time=0.0;
  
    // Get bounding boxes for time and space range
    Vector<Interval> current_set_bounds=current_set_model.range();
    ARIADNE_LOG(4,"current_set_bounds = "<<current_set_bounds<<"\n");
  
    //ARIADNE_ASSERT(initial_time_range.width() <= maximum_step_size);
  
    // Compute flow bounds and find flow bounding box
    Vector<Interval> flow_bounds; 
    Float step_size;
    make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(dynamic,current_set_bounds,maximum_step_size,maximum_bounds_diameter);
    ARIADNE_LOG(4,"step_size = "<<step_size<<"\n");
    ARIADNE_LOG(4,"flow_bounds = "<<flow_bounds<<"\n");
  
    // Compute the flow model
    FlowModelType flow_model=this->_toolbox->flow_model(dynamic,current_set_bounds,step_size,flow_bounds);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
  
    // Compute the integration time model
    TimeType next_time=current_time+step_size;
    ARIADNE_LOG(6,"next_time = "<<next_time<<"\n");
    TimeType integration_time=next_time-current_time;
  
    // Compute the flow tube (reachable set) model and the final set
    SetModelType reach_set_model=this->_toolbox->reachability_step(flow_model,current_set_model,zero_time,integration_time);         
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    SetModelType next_set_model=this->_toolbox->integration_step(flow_model,current_set_model,integration_time);
    ARIADNE_LOG(6,"next_set_model = "<<next_set_model<<"\n");
    ARIADNE_LOG(4,"Done computing evolution\n");
  
    reach_sets.adjoin(reach_set_model);

    intermediate_sets.adjoin(EnclosureType(next_set_model));
    working_sets.push_back(make_tuple(next_time,next_set_model));

  



    ARIADNE_LOG(2,"Done evolution_step.\n\n");

}




}  // namespace Ariadne

