/***************************************************************************
 *            constraint_based_hybrid_evolver.code.h
 *
 *  Copyright  2008  Pieter Collins
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
 
#include <vector>

#include "base/array.h"
#include "base/tuple.h"
#include "base/stack.h"
#include "function/function_interface.h"
#include "geometry/set_interface.h"
#include "geometry/taylor_set.h"
#include "geometry/hybrid_basic_set.h"
#include "geometry/hybrid_set.h"
#include "geometry/timed_set.h"
#include "geometry/zonotope.h"
#include "system/hybrid_automaton.h"
#include "evaluation/hybrid_time.h"
#include "evaluation/evolution_profiler.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/dynamical_toolbox.h"

#include "output/epsstream.h"
#include "output/logging.h"

#include "constraint_based_hybrid_evolver.h"


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
 
class DegenerateCrossingException { };

template<class R>
HybridEvolver<R>::~HybridEvolver()
{
}


template<class R>
HybridEvolver<R>::HybridEvolver(const EvolutionParameters<R>& p)
  : _parameters(new EvolutionParameters<R>(p)),
    _toolbox(new DynamicalToolbox< ApproximateTaylorModel<R> >())
{
}


template<class R>
HybridEvolver<R>::HybridEvolver(const HybridEvolver<R>& evolver)
  : _parameters(new EvolutionParameters<R>(*evolver._parameters)),
    _toolbox(new DynamicalToolbox< ApproximateTaylorModel<R> >(*evolver._toolbox))
{
}




template<class R>
EvolutionParameters<R>&
HybridEvolver<R>::parameters() 
{
  return *this->_parameters;
}


template<class R>
const EvolutionParameters<R>&
HybridEvolver<R>::parameters() const
{
  return *this->_parameters;
}


template<class R>
Rational
HybridEvolver<R>::maximum_step_size() const
{
  return this->_parameters->maximum_step_size();
}


template<class R>
R
HybridEvolver<R>::maximum_enclosure_radius() const
{
  return this->_parameters->maximum_enclosure_radius();
}


template<class R>
Rational
HybridEvolver<R>::lock_to_grid_time() const
{
  return this->_parameters->lock_to_grid_time();
}




/*


template<class R>
tuple< CrossingKind, ApproximateTaylorModel<R> >
HybridEvolver<R>::_crossing(const ApproximateTaylorModel<R>& guard_model, const ApproximateTaylorModel<R>& flow_model)
{
  ARIADNE_ASSERT(guard_model.result_size()=1u);
  ARIADNE_ASSERT(guard_model.argument_size()=flow_model.result_size());
  ARIADNE_ASSERT(flow_model.argument_size()=flow_model.result_size()+1u);
  uint dimension=flow_model.result_size();

  CrossingKind crossing_kind;
  Interval<R> touching_time_interval;
  ApproximateTaylorModel<R> hitting_time_model;
  ApproximateTaylorModel<R> initial_set_model;
  ApproximateTaylorModel<R> final_set_model;

  R zero_time=0;
  R step_time=flow_model.domain()[dimension].upper();


  try {
    hitting_time_model=this->_toolbox->crossing_time(flow_model,guard_model,initial_set_model,zero_time,step_time);
    crossing_kind=TRANSVERSE;
    ARIADNE_LOG(4,"\nTRANSVERSE\n\n"); 
  } 
  catch(DegenerateCrossingException) { 
    final_set_model = this->_toolbox->integration_step(flow_model,initial_set_model,step_time);
    integration_time_model=ModelType::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),Vector<R>(1u,step_time),spacial_order,smoothness);
    final_time_model = integration_time_model;
    if(definitely(this->_toolbox->active(invariant,final_set_model))) {
      flow_type=CROSSING;
      ARIADNE_LOG(4,"\nCROSSING\n\n"); 
      make_lpair(lower_touching_time,upper_touching_time)=this->_toolbox->touching_time_interval(flow_model, invariant_model,initial_set_model, zero_time, step_time);
      ARIADNE_LOG(4,"crossing_time_interval="<<lower_touching_time<<","<<upper_touching_time<<"\n\n"); 
      active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,lower_touching_time,upper_touching_time);
      // TODO: put this in!
      //jump_time_model=this->_toolbox->reachability_time(initial_time_model,crossing_time_interval.lower(),crossing_time_interval.upper());
      reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,upper_touching_time);
    } else {
      flow_type=GRAZING;
      ARIADNE_LOG(4,"\nGRAZING \n\n"); 
      make_lpair(lower_touching_time,upper_touching_time)=this->_toolbox->touching_time_interval(flow_model,invariant_model,initial_set_model,zero_time,step_time);
      ARIADNE_LOG(4,"grazing_time_interval="<<lower_touching_time<<","<<upper_touching_time<<"\n\n"); 
      active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,lower_touching_time,upper_touching_time);
      ARIADNE_LOG(6,"active_set_model="<<active_set_model<<"\n\n"); 
      // TODO: put this in!
      //jump_time_model=this->_toolbox->reachability_time(initial_time_model,grazing_time_interval.lower(),grazing_time_interval.upper());
      ARIADNE_LOG(6,"jump_time_model="<<jump_time_model<<"\n\n"); 
      reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,step_time);
      ARIADNE_LOG(6,"reach_set_model="<<reach_set_model<<"\n\n"); 
    } 
  }
}


}

*/



template<class R>
void
HybridEvolver<R>::
_evolution(EnclosureListType& final_sets, 
           EnclosureListType& reach_sets, 
           EnclosureListType& intermediate_sets, 
           const SystemType& system, 
           const EnclosureType& initial_set, 
           const TimeType& maximum_hybrid_time, 
           Semantics semantics, 
           bool reach) const
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef typename traits<R>::interval_type I;
  
  typedef HybridTime HybridTimeType;
  typedef DiscreteMode<R> ModeType;
  
  typedef DiscreteTransition<R> TransitionType;
  typedef TaylorSet<R> ContinuousEnclosureType;
  typedef HybridBasicSet<ContinuousEnclosureType> HybridEnclosureType;
  typedef TimeModelHybridBasicSet<ContinuousEnclosureType> TimedSetType;
  typedef HybridListSet<ContinuousEnclosureType> HybridListSetType;

  typedef FunctionInterface<R> FunctionType;
  typedef Box<R> BoxType;
  typedef ModelType MapModelType;
  typedef ModelType FlowModelType;
  typedef ModelType ConstraintModelType; 
  typedef ModelType SetModelType;
  typedef ModelType TimeModelType;
  typedef ModelType TimedSetModelType;


  typedef boost::shared_ptr< const FunctionInterface<R> > FunctionConstPointer;


  ARIADNE_LOG(5,__PRETTY_FUNCTION__);
  assert(semantics==upper_semantics);

  const Rational maximum_time=maximum_hybrid_time.time();
  const Integer maximum_steps=maximum_hybrid_time.steps();

  const R max_time=R(A(maximum_time));
  const R max_step_size=R(A(this->maximum_step_size()));

  const uint spacial_order=2;
  const uint temporal_order=4;
  const uint order=spacial_order+temporal_order;
  const uint smoothness=1;



  typedef tuple<DiscreteState, Integer, SetModelType, TimeModelType> HybridTimedSetType;

  stack< HybridTimedSetType > working_sets;

  {
    // Set up initial timed set models
    DiscreteState initial_discrete_state=initial_set.state();
    ModelType initial_set_model=this->_toolbox->model(initial_set.set());
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    ModelType initial_time_model
      =ModelType::constant(initial_set_model.domain(),initial_set_model.centre(),
                                           Vector<A>(1,A(0)),order,smoothness);
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push(make_tuple(initial_discrete_state,Integer(0),initial_set_model,initial_time_model));
  }


  while(!working_sets.empty()) {
    DiscreteState initial_discrete_state(0);
    Integer initial_steps;
    SetModelType initial_set_model;
    TimeModelType initial_time_model;
    make_ltuple(initial_discrete_state,initial_steps,initial_set_model,initial_time_model)=working_sets.pop();
    ARIADNE_LOG(4,"initial_steps = "<<initial_steps<<"\n");
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ARIADNE_LOG(4,"initial_discrete_state = "<<initial_discrete_state<<"\n");
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");

    const uint nd=initial_set_model.result_size();
    const uint ng=initial_set_model.argument_size();

    const DiscreteMode<R>& initial_mode=system.mode(initial_discrete_state);
    const FunctionInterface<R>& dynamic=initial_mode.dynamic();
    const FunctionInterface<R>& invariant=initial_mode.invariant();

    reference_vector< const DiscreteTransition<R> > transitions = system.transitions(initial_discrete_state);
    
    if(Box<R>(initial_set_model.range()).radius()>this->_parameters->maximum_enclosure_radius()) {
      // Subdivide
      TimedSetModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
      array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
      for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        SetModelType subdivided_set_model=project(subdivided_timed_set_model,range(0,nd));
        TimeModelType subdivided_time_model=project(subdivided_timed_set_model,range(nd,nd+1));
        working_sets.push(make_tuple(initial_discrete_state,initial_steps,subdivided_time_model,subdivided_set_model));
      }
    } else {
      // Main evolution
      Vector<I> initial_set_bounds=initial_set_model.range();
      Interval<R> initial_time_range=initial_time_model.range()[0];
      ARIADNE_LOG(4,"initial_time_range = "<<initial_time_range<<"\n");
      ARIADNE_ASSERT(initial_time_range.width() <= max_step_size);
      
      ++this->_profiler->time_steps;
      Vector<I> flow_bounds; 
      R maximum_step_size=R(A(this->maximum_step_size()));
      R maximum_bounds_diameter=mul_approx(4,maximum_step_size);
      R step_size;
      I invariant_value(-1,1);
      
      // Try to unwind time interval
      if(initial_time_range.width()>0) {
        make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(dynamic,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
        invariant_value=invariant.evaluate(flow_bounds)[0];
      }
      if(invariant_value.upper()>0) {
        make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(dynamic,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
        invariant_value=invariant.evaluate(flow_bounds)[0];
      }
      
      this->_profiler->total_stepping_time+=step_size;
      this->_profiler->minimum_time_step=std::min(Rational(R(step_size)),this->_profiler->minimum_time_step);

      ARIADNE_LOG(6,"set_bounding_box = "<<initial_set_model.range()<<"\n");
      ARIADNE_LOG(6,"flow_bounds = "<<flow_bounds<<"\n");
      ARIADNE_LOG(6,"step_size = "<<step_size<<"\n");
      ARIADNE_LOG(6,"invariant_value = "<<invariant_value<<"\n\n");
      
      //MapModelType dynamic_model(flow_bounds,dynamic,order,smoothness);
      //ARIADNE_LOG(6,"dynamic_model = "<<dynamic_model<<"\n");
      //FlowModelType flow_model=this->_toolbox->flow(dynamic_model);
      FlowModelType flow_model=this->_toolbox->flow_model(dynamic,initial_set_bounds,step_size,flow_bounds);
      ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
      
      R zero_time=0;
      R step_time=step_size;
      R final_time;
      
      TimeModelType lower_touching_time_model;
      TimeModelType upper_touching_time_model;

      ModelType guard_model;
      
      ModelType crossing_time_model;
      ModelType integration_time_model;
      
      ModelType active_set_model;
      ModelType jump_set_model;
      ModelType jump_time_model;
      ModelType final_set_model;
      ModelType final_time_model;
      ModelType reach_set_model;
      ModelType reach_time_model;
      
      const FunctionType* tmp_reset_map_ptr=0;
      const FunctionType& reset_map = *tmp_reset_map_ptr;
      
      // Compute the discrete transitions and their activation times
      reference_vector< const DiscreteTransition<R> > transitions=system.transitions(initial_discrete_state);
      reference_vector< const DiscreteTransition<R> > possibly_active_transitions;
      reference_vector< const DiscreteTransition<R> > active_transitions;

      for(uint i=0; i!=transitions.size(); ++i) {
        const DiscreteTransition<R>& transition=transitions[i];
        if(possibly(this->_toolbox->active(transition.activation(),flow_bounds))) {
          possibly_active_transitions.push_back(transition);
        }
      }

      ARIADNE_LOG(6,"possibly_active_transitions="<<possibly_active_transitions);

      //*******************************************************
      // Compute the type (TRANSVERSE, TOUCHING, CROSSING or 
      // GRAZING) of each transition.
      //*******************************************************
      CrossingKind crossing_kind;
      for(uint i=0; i!=possibly_active_transitions.size(); ++i) {
        try {
          crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);
          crossing_kind=TRANSVERSE;
        }
        catch(DegenerateCrossingException) {
          make_lpair(lower_touching_time_model,upper_touching_time_model)
            =this->_toolbox->touching_time_interval(guard_model,flow_model,initial_set_model);
          crossing_kind=TOUCHING;
        }
      }



      //********************************************************
      // OK UP TO HERE
      //********************************************************

          /*
      R reduced_step_time=step_time;
      for(uint i=0; i!=4; ++i) {
        reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,reduced_step_time);
        if(definitely(!this->_toolbox->active(invariant,reach_set_model))) {
          flow_type=MISSING;
          step_time=reduced_step_time;
            break;
          } else {
            reduced_step_time=div_approx(reduced_step_time,2);
        }
        step_time=reduced_step_time;
      }
      
      

      if(flow_type==MISSING) {
        ARIADNE_LOG(4,"\nMISSING\n\n"); 
        final_time=add_approx(initial_time_range.lower(),step_time);
        final_time=std::min(max_time,final_time);
        final_time_model=ModelType::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),Vector<R>(1u,final_time),spacial_order,smoothness);
        integration_time_model=final_time_model-initial_time_model;
        final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
        reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);
      } else {
        invariant_model=ModelType(flow_bounds,invariant,order,smoothness);
        ARIADNE_LOG(6,"invariant_model = "<<invariant_model<<"\n");   
        try {
          hitting_time_model=this->_toolbox->crossing_time(flow_model,invariant_model,initial_set_model,zero_time,step_time);
          jump_time_model=initial_time_model+hitting_time_model;
          active_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,hitting_time_model);
          reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,hitting_time_model);
          flow_type=TRANSVERSE;
          ARIADNE_LOG(4,"\nTRANSVERSE\n\n"); 
        } 
        catch(DegenerateCrossingException) { 
          final_set_model = this->_toolbox->integration_step(flow_model,initial_set_model,step_size);
          integration_time_model=ModelType::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),Vector<R>(1u,step_size),spacial_order,smoothness);
          final_time_model = initial_time_model + integration_time_model;
          if(definitely(this->_toolbox->active(invariant,final_set_model))) {
            flow_type=CROSSING;
            ARIADNE_LOG(4,"\nCROSSING\n\n"); 
            make_lpair(lower_touching_time,upper_touching_time)=this->_toolbox->touching_time_interval(flow_model, invariant_model,initial_set_model, zero_time, step_time);
            ARIADNE_LOG(4,"crossing_time_interval="<<lower_touching_time<<","<<upper_touching_time<<"\n\n"); 
            active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,lower_touching_time,upper_touching_time);
            // TODO: put this in!
            //jump_time_model=this->_toolbox->reachability_time(initial_time_model,crossing_time_interval.lower(),crossing_time_interval.upper());
            reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,upper_touching_time);
          } else {
            flow_type=GRAZING;
            ARIADNE_LOG(4,"\nGRAZING \n\n"); 
            make_lpair(lower_touching_time,upper_touching_time)=this->_toolbox->touching_time_interval(flow_model,invariant_model,initial_set_model,zero_time,step_time);
            ARIADNE_LOG(4,"grazing_time_interval="<<lower_touching_time<<","<<upper_touching_time<<"\n\n"); 
            active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,lower_touching_time,upper_touching_time);
            ARIADNE_LOG(6,"active_set_model="<<active_set_model<<"\n\n"); 
            // TODO: put this in!
            //jump_time_model=this->_toolbox->reachability_time(initial_time_model,grazing_time_interval.lower(),grazing_time_interval.upper());
            ARIADNE_LOG(6,"jump_time_model="<<jump_time_model<<"\n\n"); 
            reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,step_time);
            ARIADNE_LOG(6,"reach_set_model="<<reach_set_model<<"\n\n"); 
          } 
        }
      }
      
      
      // Compute the jump set
      if(flow_type!=MISSING) {
        Vector<A> reset_model_centre(active_set_model.evaluate(Vector<A>(active_set_model.centre())));
        ARIADNE_LOG(6,"active_set_model = "<<active_set_model<<"\n");
        ModelType reset_model(active_set_model.range(),reset_model_centre,reset_map,spacial_order,smoothness);
        ARIADNE_LOG(6,"reset_model = "<<reset_model<<"\n");
        jump_set_model=compose(reset_model,active_set_model);
        ARIADNE_LOG(6,"jump_set_model = "<<jump_set_model<<"\n");
        ARIADNE_LOG(6,"jump_time_model = "<<jump_time_model<<"\n");
      }
      
      if(flow_type==GRAZING) {
        ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
        ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
      }
      
      reach_sets.adjoin(EnclosureType(initial_discrete_state,this->_toolbox->set(reach_set_model)));
      if(flow_type!=MISSING) {
        // FIXME: Use correct discrete state
        DiscreteState jump_discrete_state=initial_discrete_state;
        EnclosureType active_set=EnclosureType(initial_discrete_state,this->_toolbox->set(active_set_model));
        EnclosureType jump_set=EnclosureType(jump_discrete_state,this->_toolbox->set(jump_set_model));
        Integer jump_steps=initial_steps+1;
        intermediate_sets.adjoin(active_set);
        intermediate_sets.adjoin(jump_set);
        working_sets.push(make_tuple(jump_discrete_state,jump_steps,jump_set_model,jump_time_model));
      } 
      if(flow_type==MISSING || flow_type==GRAZING) {
        EnclosureType final_set=EnclosureType(initial_discrete_state,this->_toolbox->set(final_set_model));
        Rational final_time_midpoint=Rational(final_time_model.range()[0].midpoint());
        if(final_time_midpoint >= maximum_time) {
          final_sets.adjoin(final_set);
        } else {
          intermediate_sets.adjoin(final_set);
          working_sets.push(make_tuple(initial_discrete_state,initial_steps,final_set_model,final_time_model));
        }
      }
          */
    }
  }
}



}  // namespace Ariadne

