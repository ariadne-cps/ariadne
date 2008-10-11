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
#include "base/stlio.h"
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
 
const bool ENABLE_SUBDIVISIONS = false;

class DegenerateCrossingException { };

template<class R>
HybridEvolver<R>::~HybridEvolver()
{
}


template<class R>
HybridEvolver<R>*
HybridEvolver<R>::clone() const
{
  return new HybridEvolver<R>(*this);
}


template<class R>
HybridEvolver<R>::HybridEvolver()
  : _parameters(new EvolutionParameters<R>()),
    _toolbox(new DynamicalToolbox< ApproximateTaylorModel<R> >()),
    _profiler(new EvolutionProfiler)
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








template<class R>
void
HybridEvolver<R>::
_evolution(EnclosureListType& final_sets, 
           EnclosureListType& reach_sets, 
           EnclosureListType& intermediate_sets, 
           const SystemType& system, 
           const EnclosureType& initial_set, 
           const TimeType& maximum_real_time, 
           Semantics semantics, 
           bool reach) const
{
  uint verbosity=0;
  
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

  ARIADNE_LOG(5,__PRETTY_FUNCTION__<<"\n");
  assert(semantics==upper_semantics);

  const Rational maximum_time=maximum_real_time;
  const Integer maximum_steps=1u;

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
    DiscreteState initial_location=initial_set.state();
    ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
    ModelType initial_set_model=this->_toolbox->model(initial_set.set());
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    ModelType initial_time_model
      =ModelType::constant(initial_set_model.domain(),initial_set_model.centre(),
                                           Vector<A>(1,A(0)),order,smoothness);
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push(make_tuple(initial_location,Integer(0),initial_set_model,initial_time_model));
  }


  while(!working_sets.empty()) {
    DiscreteState initial_location(0);
    Integer initial_steps;
    SetModelType initial_set_model;
    TimeModelType initial_time_model;
    make_ltuple(initial_location,initial_steps,initial_set_model,initial_time_model)=working_sets.pop();

    ARIADNE_LOG(4,"initial_steps = "<<initial_steps<<"\n");
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ARIADNE_LOG(4,"initial_location = "<<initial_location<<"\n");
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");

    ARIADNE_LOG(2,"steps = "<<initial_steps<<" ");
    ARIADNE_LOG(2,"time = "<<initial_time_model.range()<<" ");
    ARIADNE_LOG(2,"location = "<<initial_location<<" ");
    ARIADNE_LOG(2,"box = "<<initial_set_model.range()<<" ");
    ARIADNE_LOG(2,"radius = "<<Box<R>(initial_set_model.range()).radius()<<"\n\n");
    const uint nd=initial_set_model.result_size();
    const uint ng=initial_set_model.argument_size();

    ARIADNE_LOG(6,"modes="<<system.modes()<<"\n");
    ARIADNE_LOG(6,"transitions="<<system.transitions()<<"\n");
    reference_vector< const DiscreteTransition<R> > transitions = system.transitions(initial_location);
    ARIADNE_LOG(6,"transitions="<<transitions<<"\n");

    if(initial_time_model.range()[0].lower()>=maximum_time) {
      final_sets.adjoin(HybridEnclosureType(initial_location,this->_toolbox->set(initial_set_model)));
    } else if(ENABLE_SUBDIVISIONS
              && (Box<R>(initial_set_model.range()).radius()>this->_parameters->maximum_enclosure_radius())) 
    {
      // Subdivide
      TimedSetModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
      array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
      for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        SetModelType subdivided_set_model=project(subdivided_timed_set_model,range(0,nd));
        TimeModelType subdivided_time_model=project(subdivided_timed_set_model,range(nd,nd+1));
        working_sets.push(make_tuple(initial_location,initial_steps,subdivided_time_model,subdivided_set_model));
      }
    } else {
      ++this->_profiler->time_steps;

      const DiscreteMode<R>& initial_mode=system.mode(initial_location);
      const FunctionType* dynamic_ptr=&initial_mode.dynamic();
      const FunctionType* invariant_ptr=&initial_mode.invariant();
      const FunctionType* activation_ptr=0;
      const FunctionType* reset_ptr=0;

      
      // Main evolution
      Vector<I> initial_set_bounds=initial_set_model.range();
      Interval<R> initial_time_range=initial_time_model.range()[0];
      ARIADNE_LOG(4,"initial_time_range = "<<initial_time_range<<"\n");
      ARIADNE_ASSERT(initial_time_range.width() <= max_step_size);
      
      Vector<I> flow_bounds; 
      R maximum_step_size=R(A(this->maximum_step_size()));
      R maximum_bounds_diameter=mul_approx(4,maximum_step_size);
      R step_size;
      I invariant_value(-1,1);
     

      // Compute flow bounds and find flow bounding box
      if(initial_time_range.width()>0) {
        make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
        invariant_value=invariant_ptr->evaluate(flow_bounds)[0];
      }
      if(invariant_value.upper()>0) {
        make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
        invariant_value=invariant_ptr->evaluate(flow_bounds)[0];
      }
      
      this->_profiler->total_stepping_time+=step_size;
      this->_profiler->minimum_time_step=std::min(Rational(R(step_size)),this->_profiler->minimum_time_step);

      ARIADNE_LOG(6,"set_bounding_box = "<<initial_set_model.range()<<"\n");
      ARIADNE_LOG(6,"flow_bounds = "<<flow_bounds<<"\n");
      ARIADNE_LOG(6,"step_size = "<<step_size<<"\n");
      ARIADNE_LOG(6,"invariant_value = "<<invariant_value<<"\n\n");


      // Compute the discrete transitions and their activation times
      reference_vector< const DiscreteTransition<R> > transitions=system.transitions(initial_location);
      reference_vector< const DiscreteTransition<R> > possibly_active_transitions;
      reference_vector< const DiscreteTransition<R> > active_transitions;

      for(uint i=0; i!=transitions.size(); ++i) {
        const DiscreteTransition<R>& transition=transitions[i];
        if(possibly(this->_toolbox->active(transition.activation(),flow_bounds))) {
          possibly_active_transitions.push_back(transition);
        }
      }

      ARIADNE_LOG(6,"possibly_active_transitions="<<possibly_active_transitions);
      // For now, only consider at most one possibly active transition
      if(possibly_active_transitions.size()>1) {
        throw std::runtime_error("Current version cannot handle more than one possibly active transition at a time.");
      }

      // Compute the flow model
      FlowModelType flow_model=this->_toolbox->flow_model(*dynamic_ptr,initial_set_bounds,step_size,flow_bounds);
      ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
      
      R zero_time=0;
      R step_time=step_size;
      R final_time;
      
      ModelType integration_time_model;
      ModelType final_time_model;
      
      // Compute the crossing times of the invariant
      bool blocking=false;
      if(this->_toolbox->active(*invariant_ptr,flow_model)) {
        ARIADNE_LOG(4,"\nMISSING\n\n"); 
        R final_time=add_approx(initial_time_range.lower(),step_time);
        final_time=std::min(max_time,final_time);
        final_time_model=ModelType::constant(Vector<I>(ng,I(-1,1)),Vector<R>(ng,R(0)),Vector<A>(1u,A(final_time)),spacial_order,smoothness);
        integration_time_model=final_time_model-initial_time_model;
      } else {
        ModelType invariant_model=ModelType(flow_bounds,*invariant_ptr,order,smoothness);
        try {
          ModelType crossing_time_model=this->_toolbox->crossing_time(invariant_model,flow_model,initial_set_model);
          integration_time_model=crossing_time_model;
          final_time_model=integration_time_model+initial_time_model;
          blocking=true;
        }
        catch(DegenerateCrossingException) {
          ModelType lower_touching_time_model,upper_touching_time_model;
          make_lpair(lower_touching_time_model,upper_touching_time_model)
            =this->_toolbox->touching_time_interval(invariant_model,flow_model,initial_set_model);
          integration_time_model=upper_touching_time_model;
          final_time_model=integration_time_model+initial_time_model;
        }
      }

      ModelType final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
      ModelType reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);         
      ARIADNE_LOG(4,"Done computing continuous evolution\n");

      if(!blocking) { working_sets.push(make_tuple(initial_location,initial_steps,final_set_model,final_time_model)); }
      intermediate_sets.adjoin(HybridEnclosureType(initial_location,this->_toolbox->set(final_set_model)));

      ARIADNE_LOG(4,"intermediate_set="<<intermediate_sets[intermediate_sets.size()-1]<<"\n\n");

      // Set up the data for the events
      std::map< id_type,ModelType > transverse_event_times;
      std::map< id_type,pair<ModelType,ModelType> > touching_event_times;
      
      typedef typename std::map<DiscreteState,ModelType>::const_iterator transverse_event_iterator;
      typedef typename std::map<DiscreteState,pair<ModelType,ModelType> >::const_iterator touching_event_iterator;

      // Compute the events and their types
      for(uint i=0; i!=possibly_active_transitions.size(); ++i) {
        DiscreteEvent event(possibly_active_transitions[i].id());
        ARIADNE_LOG(6,"Testing for event "<<event<<" ");
        DiscreteState jump_location(possibly_active_transitions[i].destination().id());
        ARIADNE_LOG(6,"with target "<<jump_location<<"\n");
        Integer jump_steps = initial_steps + 1;
        activation_ptr = &possibly_active_transitions[i].activation();
        ARIADNE_LOG(6,"  activation="<<*activation_ptr<<"\n");
        reset_ptr = &possibly_active_transitions[i].reset();
        ARIADNE_LOG(6,"  reset="<<*reset_ptr<<"\n");
        ModelType activation_model=ModelType(flow_bounds,*activation_ptr,order,smoothness);
        ARIADNE_LOG(6,"\nactivation_model="<<activation_model<<"\n");
        try {
          ModelType crossing_time_model=this->_toolbox->crossing_time(flow_model,activation_model,initial_set_model);
          ARIADNE_LOG(6,"TRANSVERSE\n");
          ModelType active_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,crossing_time_model);
          Vector<A> reset_model_centre(active_set_model.evaluate(Vector<A>(active_set_model.centre())));
          ModelType reset_model(active_set_model.range(),reset_model_centre,*reset_ptr,spacial_order,smoothness);
          ARIADNE_LOG(6,"\nreset_model = "<<reset_model<<"\n");
          ModelType jump_set_model=compose(reset_model,active_set_model);
          ARIADNE_LOG(6,"\njump_set_model = "<<jump_set_model<<"\n");
          ModelType jump_time_model=initial_time_model+crossing_time_model;
          ARIADNE_LOG(6,"jump_time_model = "<<jump_time_model<<"\n");
          working_sets.push(make_tuple(jump_location,jump_steps,jump_set_model,jump_time_model));
          intermediate_sets.adjoin(HybridEnclosureType(initial_location,this->_toolbox->set(active_set_model)));
          intermediate_sets.adjoin(HybridEnclosureType(jump_location,this->_toolbox->set(jump_set_model)));
        }
        catch(DegenerateCrossingException) {
          ARIADNE_LOG(6,"TOUCHING\n");
          ModelType lower_touching_time_model, upper_touching_time_model;
          make_lpair(lower_touching_time_model, upper_touching_time_model)
            =this->_toolbox->touching_time_interval(flow_model,activation_model,initial_set_model);
          ModelType active_set_model
            =this->_toolbox->reachability_step(flow_model,initial_set_model,lower_touching_time_model,upper_touching_time_model);
          ModelType active_time_model=this->_toolbox->reachability_time(lower_touching_time_model+initial_time_model,upper_touching_time_model+initial_time_model);
          Vector<A> reset_model_centre(active_set_model.evaluate(Vector<A>(active_set_model.centre())));
          ModelType reset_model(active_set_model.range(),reset_model_centre,*reset_ptr,spacial_order,smoothness);
          ARIADNE_LOG(6,"reset_model = "<<reset_model<<"\n");
          ModelType jump_set_model=compose(reset_model,active_set_model);
          ARIADNE_LOG(6,"jump_set_model = "<<jump_set_model<<"\n");
          ModelType jump_time_model=this->_toolbox->reachability_time(initial_time_model+lower_touching_time_model,initial_time_model+upper_touching_time_model);
          ARIADNE_LOG(6,"jump_time_model = "<<jump_time_model<<"\n");
          working_sets.push(make_tuple(jump_location,jump_steps,jump_set_model,jump_time_model));
          intermediate_sets.adjoin(HybridEnclosureType(initial_location,active_set_model));
          intermediate_sets.adjoin(HybridEnclosureType(jump_location,jump_set_model));
        }
      }
      

    }
  }
}



}  // namespace Ariadne

