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
#include "orbit.h"
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


Orbit<HybridEvolver::EnclosureType> 
HybridEvolver::
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


struct DetectionData {
  int id;
  DiscreteEvent event;
  shared_ptr<const FunctionInterface> guard_ptr;
  PredicateKind predicate_kind;
  CrossingKind crossing_kind;
  tribool active;
  tribool initially_active;
  tribool finally_active;
  Interval touching_time_interval;
  ApproximateTaylorModel crossing_time_model;

  DetectionData()
    : id(0), event(0), guard_ptr()
    , predicate_kind(), crossing_kind(UNKNOWN)
    , active(indeterminate), initially_active(indeterminate), finally_active(indeterminate)
    , touching_time_interval(+1,-1), crossing_time_model() { }
};


std::ostream& operator<<(std::ostream& os, const PredicateKind& kind) {
  switch(kind) {
    case INVARIANT: return os << "INVARIANT"; 
    case ACTIVATION: return os << "ACTIVATION"; 
    case GUARD: return os << "GUARD"; 
    case TIME: return os << "TIME"; 
    case MIXED: return os << "MIXED"; 
  }
  return os << "INVALID"; 
}

std::ostream& operator<<(std::ostream& os, const CrossingKind& kind) {
  switch(kind) {
    case TRANSVERSE: return os << "TRANSVERSE"; 
      //case CROSSING: return os << "CROSSING"; 
      //case GRAZING: return os << "GRAZING"; 
    case TOUCHING: return os << "TOUCHING"; 
      //case ACTIVE: return os << "ACTIVE"; 
      //case MISSING: return os << "MISSING"; 
    case NONE: return os << "NONE"; 
    case UNKNOWN: return os << "UNKNOWN"; 
  }
  return os << "INVALID"; 
}

std::ostream& operator<<(std::ostream& os, const DetectionData& data) {
  os << std::boolalpha;
  os << " { id="<<data.id;
  os << ", event="<<data.event;
  os << ", kind="<<data.predicate_kind;
  os << ", crossing="<<data.crossing_kind;
  os << ", active="<<data.active;
  if(data.crossing_kind!=NONE) {
    os << ", initially_active="<<data.initially_active;
    os << ", finally_active="<<data.finally_active;
    if(data.touching_time_interval.lower()<=data.touching_time_interval.upper()) { 
      os << ", touching_times="<<data.touching_time_interval; }
  }
  return os <<" } "; 
}

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
  verbosity=0;

  
  typedef FunctionInterface FunctionType;
  typedef Vector<Interval> BoxType;
  typedef ModelType MapModelType;
  typedef ModelType FlowModelType;
  typedef ModelType ConstraintModelType; 
  typedef ModelType SetModelType;
  typedef ModelType TimeModelType;
  typedef ModelType TimedSetModelType;


  typedef boost::shared_ptr< const FunctionInterface > FunctionConstPointer;

  ARIADNE_LOG(5,__PRETTY_FUNCTION__<<"\n");
  assert(semantics==UPPER_SEMANTICS);

  const IntegerType maximum_steps=maximum_hybrid_time.second;
  const Float maximum_time=maximum_hybrid_time.first;

  const uint spacial_order=2;
  const uint temporal_order=4;
  const uint order=spacial_order+temporal_order;
  const uint smoothness=1;



  typedef tuple<DiscreteState, IntegerType, SetModelType, TimeModelType> HybridTimedSetType;

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
    working_sets.push_back(make_tuple(initial_location,IntegerType(0),initial_set_model,initial_time_model));
  }


  while(!working_sets.empty()) {
    HybridTimedSetType current_set=working_sets.back();
    working_sets.pop_back();
    DiscreteState initial_location=current_set.first;
    IntegerType initial_steps=current_set.second;
    SetModelType initial_set_model=current_set.third;
    TimeModelType initial_time_model=current_set.fourth;
    if(initial_time_model.range()[0].lower()>=maximum_time || initial_steps>=maximum_steps) {
      final_sets.adjoin(EnclosureType(initial_location,this->_toolbox->set(initial_set_model)));
    } else if(ENABLE_SUBDIVISIONS
              && (radius(initial_set_model.range())>this->_parameters->maximum_enclosure_radius)) 
    {
      // Subdivide
      uint nd=initial_set_model.result_size();
      TimedSetModelType initial_timed_set_model=join(initial_set_model,initial_time_model);
      array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
      for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        SetModelType subdivided_set_model=Ariadne::project(subdivided_timed_set_model,range(0,nd));
        TimeModelType subdivided_time_model=Ariadne::project(subdivided_timed_set_model,range(nd,nd+1));
        working_sets.push_back(make_tuple(initial_location,initial_steps,subdivided_time_model,subdivided_set_model));
      }
    } else {
      this->_evolution_step(working_sets,
                            final_sets,reach_sets,intermediate_sets,
                            system,current_set,maximum_hybrid_time,
                            semantics,reach);
    }
  }

}





void
HybridEvolver::
_evolution_step(std::vector< HybridTimedSetType >& working_sets,
                EnclosureListType& final_sets, 
                EnclosureListType& reach_sets, 
                EnclosureListType& intermediate_sets, 
                const SystemType& system, 
                const HybridTimedSetType& current_set,
                const TimeType& maximum_hybrid_time, 
                Semantics semantics, 
                bool reach) const
{
  typedef FunctionInterface FunctionType;
  typedef Vector<Interval> BoxType;
  typedef ModelType MapModelType;
  typedef ModelType FlowModelType;
  typedef ModelType ConstraintModelType; 
  typedef ModelType SetModelType;
  typedef ModelType TimeModelType;
  typedef ModelType TimedSetModelType;
  
  DiscreteState initial_location(0);
  IntegerType initial_steps;
  SetModelType initial_set_model;
  TimeModelType initial_time_model;
  ARIADNE_LOG(9,"working_set = "<<current_set<<"\n");
  make_ltuple(initial_location,initial_steps,initial_set_model,initial_time_model)=current_set;

  ARIADNE_LOG(4,"initial_steps = "<<initial_steps<<"\n");
  ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"");
  ARIADNE_LOG(4,"initial_location = "<<initial_location<<"\n");
  ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
  
  ARIADNE_LOG(2,"steps = "<<initial_steps<<" ");
  ARIADNE_LOG(2,"time_range = "<<initial_time_model.range()<<" ");
  ARIADNE_LOG(2,"location = "<<initial_location<<" ");
  ARIADNE_LOG(2,"box = "<<initial_set_model.range()<<" ");
  ARIADNE_LOG(2,"radius = "<<radius(initial_set_model.range())<<"\n\n");
  //const uint nd=initial_set_model.result_size();
  //const uint ng=initial_set_model.argument_size();
  

  /////////////// Main Evolution ////////////////////////////////

  const DiscreteMode& initial_mode=system.mode(initial_location);
  const FunctionType* dynamic_ptr=&initial_mode.dynamic();
  
  ARIADNE_LOG(6,"mode="<<initial_mode<<"\n");
  std::vector< shared_ptr<const FunctionInterface> > invariants
    =initial_mode.invariants();
  ARIADNE_LOG(7,"invariants="<<invariants<<"\n");
  const std::set< DiscreteTransition > transitions = system.transitions(initial_location);
  ARIADNE_LOG(7,"transitions="<<transitions<<"\n");
  
  
  
  // Set evolution parameters
  const Float maximum_step_size=this->_parameters->maximum_step_size;
  const Float maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
  const Float zero_time=0.0;
  
  // Get bounding boxes for time and space range
  Vector<Interval> initial_set_bounds=initial_set_model.range();
  ARIADNE_LOG(4,"initial_set_range = "<<initial_set_bounds<<"\n");
  Interval initial_time_range=initial_time_model.range()[0];
  ARIADNE_LOG(4,"initial_time_range = "<<initial_time_range<<"\n");
  
  //ARIADNE_ASSERT(initial_time_range.width() <= maximum_step_size);
  
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
  
  
  





      
  // Compute transitions
  
  ARIADNE_LOG(2,"Evolution step.\n");

  ARIADNE_LOG(7,"initial_time_bounds="<<initial_time_model.range()<<"\n");
  ARIADNE_LOG(7,"flow_bounds="<<flow_bounds<<"\n");
  ARIADNE_LOG(7,"intermediate_flow_bounds="<<reach_set_model.range()<<"\n");
  ARIADNE_LOG(7,"initial_flow_bounds="<<initial_set_model.range()<<"\n");
  ARIADNE_LOG(7,"final_flow_bounds="<<final_set_model.range()<<"\n");
  ARIADNE_LOG(7,"time_interval="<<Interval(0,step_size)<<"\n");
  
  std::map<int,DetectionData> detection_data;
  ARIADNE_LOG(9,"invariants="<<invariants<<"\n");
  for(std::vector< shared_ptr<const FunctionInterface> >::const_iterator
        iter=invariants.begin(); iter!=invariants.end(); ++iter)
  {
    DetectionData new_data;
    new_data.id=detection_data.size();
    new_data.event=BLOCKING_EVENT;
    new_data.predicate_kind=INVARIANT;
    new_data.crossing_kind=UNKNOWN;
    new_data.guard_ptr=(*iter);
    detection_data.insert(make_pair(new_data.id,new_data));
  }

  ARIADNE_LOG(9,"transitions="<<transitions<<"\n");
  for(std::set< DiscreteTransition >::const_iterator 
        iter=transitions.begin(); iter!=transitions.end(); ++iter)
  {
    DetectionData new_data;
    new_data.id=detection_data.size();
    new_data.event=iter->event();
    new_data.predicate_kind=iter->forced() ? GUARD : ACTIVATION;
    new_data.crossing_kind=UNKNOWN;
    new_data.guard_ptr=iter->activation_ptr();
        detection_data.insert(make_pair(new_data.id,new_data));
  }
  
  ARIADNE_LOG(8,"detection_data="<<detection_data<<"\n");
  
  typedef std::map<int,DetectionData>::iterator predicate_iterator;
  for(predicate_iterator iter=detection_data.begin(); 
      iter!=detection_data.end(); ++iter)
  {
    DetectionData& data=iter->second;
    if(possibly(this->_toolbox->active(*data.guard_ptr,flow_bounds))) {
      data.active=this->_toolbox->active(*data.guard_ptr,reach_set_model);
      data.initially_active=data.active;
      data.finally_active=data.active;
      if(data.active) {
        data.crossing_kind=NONE;
      } else if(!data.active) {
        data.crossing_kind=NONE;
      } else {
        data.initially_active=this->_toolbox->active(*data.guard_ptr,initial_set_model);
        data.finally_active=this->_toolbox->active(*data.guard_ptr,final_set_model);
        ModelType guard_model=this->_toolbox->predicate_model(*data.guard_ptr,flow_bounds);
        if(data.initially_active xor data.finally_active) {
          ARIADNE_LOG(7," Testing crossing time for: "<<data<<"\n");
          try {
            data.crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);
            data.touching_time_interval=data.crossing_time_model.range()[0];
            data.crossing_kind=TRANSVERSE;
          }
          catch(DegenerateCrossingException) { ARIADNE_LOG(7," DegenerateCrossing\n"); }
        }
        if(data.crossing_kind!=TRANSVERSE) {
          data.touching_time_interval=this->_toolbox->touching_time_interval(guard_model,flow_model,initial_set_model);
          data.crossing_kind=TOUCHING;
        }
      }
    } else {
      data.active=false;
      data.crossing_kind=NONE;                
    }
  }
  
  ARIADNE_LOG(6,"detection_data="<<detection_data<<"\n");
  
  // Remove inactive transitions
  for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ) {
    if(!iter->second.active) {
         detection_data.erase(iter++);
    } else {
      ++iter;
    }
  }

  ARIADNE_LOG(6,"active_detection_data="<<detection_data<<"\n");
  
  
  



  

  // Compute the evolution time for the blocking evolution
  ARIADNE_LOG(6,"Testing for blocking predicate\n");
  DetectionData blocking_data;
  blocking_data.predicate_kind=TIME;
  blocking_data.crossing_kind=NONE;
  blocking_data.touching_time_interval=step_size;
  blocking_data.active=false;
  blocking_data.initially_active=false;
  blocking_data.finally_active=false;
  for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ++iter) {
    DetectionData& data=iter->second;
    if(data.predicate_kind==INVARIANT || data.predicate_kind==GUARD) {
      if(data.initially_active) {
        blocking_data.predicate_kind=data.predicate_kind;
        blocking_data.active=indeterminate;
        blocking_data.initially_active=true;
        blocking_data.finally_active=indeterminate;
        blocking_data.crossing_kind=NONE;
        ARIADNE_LOG(6,"No continuous evolution possible\n");
        break;
      } else if(data.finally_active) {
        blocking_data.predicate_kind=data.predicate_kind;
        if(data.touching_time_interval.lower()<blocking_data.touching_time_interval.upper()) {
          blocking_data.active=indeterminate;
          if(data.touching_time_interval.upper()<blocking_data.touching_time_interval.lower() && data.crossing_kind==TRANSVERSE) {
            blocking_data.predicate_kind=data.predicate_kind;
            blocking_data.crossing_kind=TRANSVERSE;
            blocking_data.finally_active=true;
            blocking_data.touching_time_interval=data.touching_time_interval;
            blocking_data.crossing_time_model=data.crossing_time_model;
            blocking_data.event=data.event;
          } else {
            blocking_data.predicate_kind=MIXED;
            blocking_data.crossing_kind=TOUCHING;
            blocking_data.initially_active=blocking_data.initially_active or data.initially_active;
            blocking_data.finally_active=blocking_data.initially_active or data.finally_active;
            blocking_data.touching_time_interval=min(data.touching_time_interval,blocking_data.touching_time_interval);
          }
        }
      }
    }
  }
  
  ARIADNE_LOG(6,"blocking_data="<<blocking_data<<"\n");

  
  Float maximum_evolution_time;
  if(semantics==UPPER_SEMANTICS) {
    maximum_evolution_time=blocking_data.touching_time_interval.upper();
  } else {
    maximum_evolution_time=blocking_data.touching_time_interval.lower();
  }

    
  if(maximum_evolution_time>=step_size) {
    ARIADNE_LOG(6,"No restrictions on continuous evolution\n");
    maximum_evolution_time=step_size;
  }
  
  ARIADNE_LOG(6,"maximum_evolution_time="<<maximum_evolution_time<<"\n");
  
  
  // Remove transitions which only occur after the blocking time
  for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ) {
    if(!iter->second.initially_active && iter->second.touching_time_interval.lower()>maximum_evolution_time) {
      detection_data.erase(iter++);
    } else {
      ++iter;
    }
  }
  
 
  ARIADNE_LOG(6,"active_detection_data="<<detection_data<<"\n");

  // Compute continuous evolution
  ARIADNE_LOG(6,"Computing continuous evolution\n");
  if(blocking_data.initially_active) {
    ARIADNE_LOG(6,"  No continuous evolution\n");
    reach_set_model=initial_set_model;
    reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
  } else if(blocking_data.crossing_kind==TRANSVERSE) {
    ARIADNE_LOG(6,"  Continuous evolution to transverse invariant/guard\n");
    reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,blocking_data.crossing_time_model);
    reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
  } else if(blocking_data.finally_active) {
    ARIADNE_LOG(6,"  Continuous evolution for finite time to invariant/guard\n");
    reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,blocking_data.touching_time_interval.upper());
    final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,blocking_data.touching_time_interval.upper());
    reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
  } else {
    ARIADNE_LOG(6,"  Continuous evolution for maximum time; unwinding time differentce\n");
    ARIADNE_ASSERT(blocking_data.predicate_kind==TIME);
    // Try to unwind any differences in time
    Interval initial_time_range=initial_time_model.range()[0];
    Float final_time=initial_time_range.lower()+step_size;
    ARIADNE_LOG(8,"    initial_time_model="<<initial_time_model<<"\n");
    final_time_model=this->_toolbox->constant_time_model(final_time,initial_time_model.domain());
    ARIADNE_LOG(8,"    final_time_model="<<final_time_model<<"\n");
    TimeModelType integration_time_model=final_time_model-initial_time_model;
    ARIADNE_LOG(8,"    integration_time_model="<<integration_time_model<<"\n");
    reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);
    final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
    reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
    intermediate_sets.adjoin(EnclosureType(initial_location,final_set_model));
    working_sets.push_back(make_tuple(initial_location,initial_steps,final_set_model,final_time_model));
  }

  // Process discrete transitions
  ARIADNE_LOG(6,"Computing discrete transitions\n");
  for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ++iter) {
    DetectionData& data=iter->second;
    ARIADNE_LOG(6,"  transition"<<data<<"\n");
    if(data.predicate_kind==ACTIVATION || data.predicate_kind==GUARD) {
      const DiscreteTransition& transition=system.transition(data.event,initial_location);
      DiscreteState jump_location=transition.target().location();
      const FunctionInterface* reset_ptr=&transition.reset();
      TimeModelType active_time_model;
      TimeModelType jump_time_model;
      if(blocking_data.initially_active) {
        SetModelType jump_set_model
          =this->_toolbox->compose(*reset_ptr,initial_set_model);
        
        working_sets.push_back(make_tuple(jump_location,initial_steps+1,jump_set_model,initial_time_model));
        reach_sets.adjoin(EnclosureType(initial_location,initial_set_model));
        intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
      } else {
        if(data.crossing_kind==TRANSVERSE) {
          if(data.initially_active) {
            active_time_model=this->_toolbox->reachability_time(zero_time,data.crossing_time_model);
            jump_time_model=this->_toolbox->reachability_time(initial_time_model,initial_time_model+data.crossing_time_model);
          } else {
            active_time_model=this->_toolbox->reachability_time(data.crossing_time_model,maximum_evolution_time);
            jump_time_model=this->_toolbox->
              reachability_time(initial_time_model+data.crossing_time_model,
                                initial_time_model+Vector<Float>(1u,maximum_evolution_time));
          } 
        } else {
          Float lower_active_time
            = data.initially_active ? zero_time : data.touching_time_interval.lower();
          Float upper_active_time
            = data.finally_active ? step_size : data.touching_time_interval.upper();
          TimeModelType lower_active_time_model=this->_toolbox->constant_time_model(lower_active_time,initial_time_model.domain());
          TimeModelType upper_active_time_model=this->_toolbox->constant_time_model(upper_active_time,initial_time_model.domain());
          active_time_model=this->_toolbox->reachability_time(lower_active_time_model,upper_active_time_model);
          jump_time_model=this->_toolbox->reachability_time(initial_time_model+lower_active_time_model,initial_time_model+upper_active_time_model);
        }
        SetModelType active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,active_time_model);
        SetModelType jump_set_model=this->_toolbox->compose(*reset_ptr,active_set_model);
        
        working_sets.push_back(make_tuple(jump_location,initial_steps+1,jump_set_model,jump_time_model));
        intermediate_sets.adjoin(EnclosureType(initial_location,active_set_model));
        intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
      }
    }
  }
  
  ARIADNE_LOG(2,"Done evolution_step.\n\n");

}




}  // namespace Ariadne

