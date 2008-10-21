/***************************************************************************
 *            hybrid_evolver.cc
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
#include "function_interface.h"
#include "approximate_taylor_model.h"
#include "dynamical_toolbox.h"
#include "evolution_parameters.h"

#include "logging.h"

#include "hybrid_automaton.h"
#include "hybrid_evolver.h"

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
static const int BLOCKING_EVENT = -2;
using boost::shared_ptr;

class DegenerateCrossingException { };


HybridEvolver::HybridEvolver()
  : _parameters(new EvolutionParameters()),
    _toolbox(new DynamicalToolbox<ModelType>())
{
}



HybridEvolver::HybridEvolver(const EvolutionParameters& p)
  : _parameters(new EvolutionParameters(p)),
    _toolbox(new DynamicalToolbox<ModelType>())
{
}








struct DetectionData {
  int id;
  DiscreteEvent event;
  shared_ptr<const FunctionInterface> guard;
  int kind;
  bool initially_active;
  bool finally_active;
  Interval touching_times;
  ApproximateTaylorModel hitting_time;
};

inline bool operator<(const DetectionData& d1, const DetectionData& d2) {
  return d1.id < d2.id;
}

void
HybridEvolver::
_evolution(EnclosureListType& final_sets, 
           EnclosureListType& reach_sets, 
           EnclosureListType& intermediate_sets, 
           const SystemType& system, 
           const EnclosureType& initial_set, 
           const TimeType& maximum_hybrid_time, 
           Semantics semantics, 
           bool reach) const
{
  uint verbosity=0;

  
  typedef FunctionInterface FunctionType;
  typedef Vector<Interval> BoxType;
  typedef ModelType MapModelType;
  typedef ModelType FlowModelType;
  typedef ModelType ConstraintModelType; 
  typedef ModelType SetModelType;
  typedef ModelType TimeModelType;
  typedef ModelType TimedSetModelType;

  typedef uint Integer;

  typedef boost::shared_ptr< const FunctionInterface > FunctionConstPointer;

  ARIADNE_LOG(5,__PRETTY_FUNCTION__<<"\n");
  assert(semantics==upper_semantics);

  const uint maximum_steps=maximum_hybrid_time.second;
  const Float maximum_time=maximum_hybrid_time.first;
  const Float maximum_step_size=this->_parameters->maximum_step_size;

  const uint spacial_order=2;
  const uint temporal_order=4;
  const uint order=spacial_order+temporal_order;
  const uint smoothness=1;



  typedef tuple<DiscreteState, Integer, SetModelType, TimeModelType> HybridTimedSetType;

  std::vector< HybridTimedSetType > working_sets;

  {
    // Set up initial timed set models
    ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
    DiscreteState initial_location;
    ContinuousEnclosureType initial_continuous_set;
    make_lpair(initial_location,initial_continuous_set)=initial_set;
    ARIADNE_LOG(6,"initial_location = "<<initial_location<<"\n");
    ModelType initial_set_model=this->_toolbox->model(initial_continuous_set);
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    ModelType initial_time_model
      =ModelType::constant(initial_set_model.domain(),initial_set_model.centre(),
                           Vector<Float>(1,Float(0)),order,smoothness);
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push_back(make_tuple(initial_location,Integer(0),initial_set_model,initial_time_model));
  }


  while(!working_sets.empty()) {
    DiscreteState initial_location(0);
    Integer initial_steps;
    SetModelType initial_set_model;
    TimeModelType initial_time_model;
    ARIADNE_LOG(9,"working_set = "<<working_sets.back()<<"\n");
    make_ltuple(initial_location,initial_steps,initial_set_model,initial_time_model)=working_sets.back();
    working_sets.pop_back();

    ARIADNE_LOG(4,"initial_steps = "<<initial_steps<<"\n");
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"");
    ARIADNE_LOG(4,"initial_location = "<<initial_location<<"\n");
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");

    ARIADNE_LOG(2,"steps = "<<initial_steps<<" ");
    ARIADNE_LOG(2,"time_range = "<<initial_time_model.range()<<" ");
    ARIADNE_LOG(2,"location = "<<initial_location<<" ");
    ARIADNE_LOG(2,"box = "<<initial_set_model.range()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(initial_set_model.range())<<"\n\n");
    const uint nd=initial_set_model.result_size();
    const uint ng=initial_set_model.argument_size();

    ARIADNE_LOG(6,"modes="<<system.modes()<<"\n");
    ARIADNE_LOG(6,"transitions="<<system.transitions()<<"\n");
    const std::set< DiscreteTransition > transitions = system.transitions(initial_location);
    ARIADNE_LOG(6,"transitions="<<transitions<<"\n");

    if(initial_time_model.range()[0].lower()>=maximum_time) {
      final_sets.adjoin(EnclosureType(initial_location,this->_toolbox->set(initial_set_model)));
    } else if(ENABLE_SUBDIVISIONS
              && (radius(initial_set_model.range())>this->_parameters->maximum_enclosure_radius)) 
    {
      // Subdivide
      TimedSetModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
      array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
      for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        SetModelType subdivided_set_model=Ariadne::project(subdivided_timed_set_model,range(0,nd));
        TimeModelType subdivided_time_model=Ariadne::project(subdivided_timed_set_model,range(nd,nd+1));
        working_sets.push_back(make_tuple(initial_location,initial_steps,subdivided_time_model,subdivided_set_model));
      }
    } else {

      /////////////// Main Evolution ////////////////////////////////
      const DiscreteMode& initial_mode=system.mode(initial_location);
      const FunctionType* dynamic_ptr=&initial_mode.dynamic();
      

      std::set<DetectionData> transition_data;
      std::vector< shared_ptr<const FunctionInterface> > invariants
        =initial_mode.invariants();
      ARIADNE_LOG(9,"invariants="<<invariants<<"\n");
      for(std::vector< shared_ptr<const FunctionInterface> >::const_iterator
            iter=invariants.begin(); iter!=invariants.end(); ++iter)
      {
        DetectionData new_data;
        new_data.id=transition_data.size();
        new_data.event=BLOCKING_EVENT;
        new_data.guard=(*iter);
        transition_data.insert(new_data);
      }

      std::set< DiscreteTransition > transitions
        =system.transitions(initial_mode.location());
      ARIADNE_LOG(9,"transitions="<<transitions<<"\n");
      for(std::set< DiscreteTransition >::const_iterator 
            iter=transitions.begin(); iter!=transitions.end(); ++iter)
      {
        DetectionData new_data;
        new_data.id=transition_data.size();
        new_data.event=iter->event();
        new_data.guard=iter->activation_ptr();
        transition_data.insert(new_data);
      }
      


      // Set evolution parameters
      const Float maximum_step_size=this->_parameters->maximum_step_size;
      const Float maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
      const Float zero_time=0.0;

      // Get bounding boxes for time and space range
      Vector<Interval> initial_set_bounds=initial_set_model.range();
      ARIADNE_LOG(4,"initial_set_range = "<<initial_set_bounds<<"\n");
      Interval initial_time_range=initial_time_model.range()[0];
      ARIADNE_LOG(4,"initial_time_range = "<<initial_time_range<<"\n");
      
      ARIADNE_ASSERT(initial_time_range.width() <= maximum_step_size);
      
      // Compute flow bounds and find flow bounding box
      Vector<Interval> flow_bounds; 
      Float step_size;
      make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
      ARIADNE_LOG(4,"step_size = "<<step_size<<"\n");
      ARIADNE_LOG(4,"flow_bounds = "<<flow_bounds<<"\n");
   
      // Compute the flow model
      FlowModelType flow_model=this->_toolbox->flow_model(*dynamic_ptr,initial_set_bounds,step_size,flow_bounds);
      ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");

      // Compute the integration time model
      TimeModelType final_time_model=initial_time_model+Vector<Float>(1u,step_size);
      ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
      TimeModelType integration_time_model=final_time_model-initial_time_model;

      // Compute the flow tube (reachable set) model and the final set
      ModelType final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
      ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
      ModelType reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);         
      ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
      ARIADNE_LOG(4,"Done computing continuous evolution\n");






      static const bool NO_DISCRETE_EVOLUTION=true;
      
      if(NO_DISCRETE_EVOLUTION) {
        working_sets.push_back(make_tuple(initial_location,initial_steps,final_set_model,final_time_model));
        intermediate_sets.adjoin(EnclosureType(initial_location,reach_set_model));
        final_sets.adjoin(EnclosureType(initial_location,final_set_model));
      }





      /*
      Interval invariant_value;
      const FunctionInterface* invariant_ptr;
      const FunctionInterface* activation_ptr;
      const FunctionInterface* reset_ptr;
      
    
     

      // Compute flow bounds and find flow bounding box
      if(initial_time_range.width()>0) {
        make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
      }
      
      // Compute the values of the invariants over the flow bo
      invariant_value=invariant_ptr->evaluate(flow_bounds)[0];
      
      if(invariant_value.upper()>0) {
        make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
        invariant_value=invariant_ptr->evaluate(flow_bounds)[0];
      }
      
      ARIADNE_LOG(6,"set_bounding_box = "<<initial_set_model.range()<<"\n");
      ARIADNE_LOG(6,"flow_bounds = "<<flow_bounds<<"\n");
      ARIADNE_LOG(6,"step_size = "<<step_size<<"\n");
      ARIADNE_LOG(6,"invariant_value = "<<invariant_value<<"\n\n");


      // Compute the discrete transitions and their activation times
      std::set< DiscreteTransition > transition_set=system.transitions(initial_location);
      std::vector< DiscreteTransition > all_transitions(transition_set.begin(),transition_set.end());
      std::vector< DiscreteTransition > possibly_active_transitions;
      std::vector< DiscreteTransition > active_transitions;

      for(uint i=0; i!=transitions.size(); ++i) {
        const DiscreteTransition& transition=all_transitions[i];
        if(possibly(this->_toolbox->active(transition.activation(),flow_bounds))) {
          possibly_active_transitions.push_back(transition);
        }
      }

      ARIADNE_LOG(6,"possibly_active_transitions="<<possibly_active_transitions);
      // For now, only consider at most one possibly active transition
      if(possibly_active_transitions.size()>1) {
        throw std::runtime_error("Current version cannot handle more than one possibly active transition at a time.");
      }

      Float step_time=step_size;
      Float final_time;
      
      // Compute the crossing times of the invariant
      bool blocking=false;
      if(this->_toolbox->active(*invariant_ptr,flow_model)) {
        ARIADNE_LOG(4,"\nMISSING\n\n"); 
        Float final_time=add_approx(initial_time_range.lower(),step_time);
        final_time=std::min(maximum_time,final_time);
        Vector<Interval> final_time_domain(ng,Interval(-1,1));
        Vector<Float> final_time_centre(ng,Float(0));
        Vector<Float> final_time_value(1u,final_time);
        final_time_model=ModelType::constant(final_time_domain,final_time_centre,final_time_value,spacial_order,smoothness);
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

      final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
      reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);         
      ARIADNE_LOG(4,"Done computing continuous evolution\n");

      if(!blocking) { working_sets.push_back(make_tuple(initial_location,initial_steps,final_set_model,final_time_model)); }
      intermediate_sets.adjoin(EnclosureType(initial_location,this->_toolbox->set(final_set_model)));

      ARIADNE_LOG(4,"intermediate_set="<<intermediate_sets[intermediate_sets.size()-1]<<"\n\n");

      // Set up the data for the events
      std::map< DiscreteEvent,ModelType > transverse_event_times;
      std::map< DiscreteEvent,pair<ModelType,ModelType> > touching_event_times;
      
      typedef std::map<DiscreteEvent,ModelType>::const_iterator transverse_event_iterator;
      typedef std::map<DiscreteEvent,pair<ModelType,ModelType> >::const_iterator touching_event_iterator;

      // Compute the events and their types
      for(uint i=0; i!=possibly_active_transitions.size(); ++i) {
        DiscreteEvent event(possibly_active_transitions[i].event());
        ARIADNE_LOG(6,"Testing for event "<<event<<" ");
        DiscreteState jump_location(possibly_active_transitions[i].target().location());
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
          Vector<Float> reset_model_centre(active_set_model.evaluate(Vector<Float>(active_set_model.centre())));
          ModelType reset_model(active_set_model.range(),reset_model_centre,*reset_ptr,spacial_order,smoothness);
          ARIADNE_LOG(6,"\nreset_model = "<<reset_model<<"\n");
          ModelType jump_set_model=compose(reset_model,active_set_model);
          ARIADNE_LOG(6,"\njump_set_model = "<<jump_set_model<<"\n");
          ModelType jump_time_model=initial_time_model+crossing_time_model;
          ARIADNE_LOG(6,"jump_time_model = "<<jump_time_model<<"\n");
          working_sets.push_back(make_tuple(jump_location,jump_steps,jump_set_model,jump_time_model));
          intermediate_sets.adjoin(EnclosureType(initial_location,this->_toolbox->set(active_set_model)));
          intermediate_sets.adjoin(EnclosureType(jump_location,this->_toolbox->set(jump_set_model)));
        }
        catch(DegenerateCrossingException) {
          ARIADNE_LOG(6,"TOUCHING\n");
          ModelType lower_touching_time_model, upper_touching_time_model;
          make_lpair(lower_touching_time_model, upper_touching_time_model)
            =this->_toolbox->touching_time_interval(flow_model,activation_model,initial_set_model);
          ModelType active_set_model
            =this->_toolbox->reachability_step(flow_model,initial_set_model,lower_touching_time_model,upper_touching_time_model);
          ModelType active_time_model=this->_toolbox->reachability_time(lower_touching_time_model+initial_time_model,upper_touching_time_model+initial_time_model);
          Vector<Float> reset_model_centre(active_set_model.evaluate(Vector<Float>(active_set_model.centre())));
          ModelType reset_model(active_set_model.range(),reset_model_centre,*reset_ptr,spacial_order,smoothness);
          ARIADNE_LOG(6,"reset_model = "<<reset_model<<"\n");
          ModelType jump_set_model=compose(reset_model,active_set_model);
          ARIADNE_LOG(6,"jump_set_model = "<<jump_set_model<<"\n");
          ModelType jump_time_model=this->_toolbox->reachability_time(initial_time_model+lower_touching_time_model,initial_time_model+upper_touching_time_model);
          ARIADNE_LOG(6,"jump_time_model = "<<jump_time_model<<"\n");
          working_sets.push_back(make_tuple(jump_location,jump_steps,jump_set_model,jump_time_model));
          intermediate_sets.adjoin(EnclosureType(initial_location,active_set_model));
          intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
        }
      }
      */  

    }
  }
}



}  // namespace Ariadne

