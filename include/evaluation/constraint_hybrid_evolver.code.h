/***************************************************************************
 *            constraint_hybrid_evolver.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#include "../geometry/rectangle_expression.h"
#include "../geometry/set_interface.h"
#include "../geometry/hybrid_set.h"
#include "../geometry/timed_set.h"
#include "../system/map.h"
#include "../system/vector_field.h"
#include "../system/constraint_hybrid_automaton.h"
#include "../evaluation/applicator.h"
#include "../evaluation/integrator.h"

#include "../evaluation/lohner_integrator.h"

#include "../output/epsfstream.h"
#include "../output/logging.h"

#include "constraint_hybrid_evolver.h"


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
 
namespace Evaluation { static int& verbosity = hybrid_evolver_verbosity; }

using namespace Numeric;
using namespace LinearAlgebra;
using namespace Geometry;
using namespace System;

template<class R>
Evaluation::ConstraintHybridEvolver<R>::~ConstraintHybridEvolver()
{
  delete this->_applicator;
  delete this->_integrator;
}


template<class R>
Evaluation::ConstraintHybridEvolver<R>::ConstraintHybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _applicator(a.clone()), _integrator(dynamic_cast<LohnerIntegrator<R>*>(i.clone()))
{
  if(!this->_integrator) {
    throw std::runtime_error("ConstraintHybridEvolver::ConstraintHybridEvolver(Applicator a, Integrator i): Invalid integrator");
  }
}

template<class R>
Evaluation::ConstraintHybridEvolver<R>::ConstraintHybridEvolver(const ConstraintHybridEvolver<R>& evolver)
  :_applicator(evolver._applicator), _integrator(evolver._integrator)
{
}



template<class R>
const std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>&
Evaluation::ConstraintHybridEvolver<R>::trace() const
{
  return this->_trace;
}


template<class R>
time_type
Evaluation::ConstraintHybridEvolver<R>::maximum_step_size() const
{
  return this->_integrator->maximum_step_size();
}


template<class R>
R
Evaluation::ConstraintHybridEvolver<R>::maximum_basic_set_radius() const
{
  return this->_integrator->maximum_basic_set_radius();
}


template<class R> inline
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::subdivide(const timed_set_type& timed_set)
{
  ARIADNE_LOG(8,"HybridEvolver:subdivide(...)\n");
  ARIADNE_LOG(9,"  radius="<<timed_set.continuous_state_set().radius()<<"\n");
  ARIADNE_LOG(9,"  bounding_box="<<timed_set.continuous_state_set().bounding_box()<<"\n");
  ARIADNE_LOG(9,"  set="<<timed_set<<"\n");
  std::vector<timed_set_type> result;
  Zonotope<R> z=orthogonal_over_approximation(timed_set.continuous_state_set());
  ListSet< Zonotope<R> > subdivisions=Geometry::subdivide(z);
  for(typename ListSet< Zonotope<R> >::const_iterator iter=subdivisions.begin();
      iter!=subdivisions.end(); ++iter)
  {
    result.push_back(timed_set_type(timed_set.time(),timed_set.steps(),timed_set.discrete_state(),*iter));
  }
  return result;
}


template<class R> 
typename Evaluation::ConstraintHybridEvolver<R>::working_sets_type
Evaluation::ConstraintHybridEvolver<R>::_compute_working_sets(const hybrid_list_set_type& set) const
{
  working_sets_type working_sets;
  for(typename hybrid_list_set_type::const_iterator loc_iter=set.begin();
      loc_iter!=set.end(); ++loc_iter)
  {
    id_type loc_id=loc_iter->first;
    const Geometry::ListSet<continuous_basic_set_type>& list_set=*loc_iter->second;
    for(typename Geometry::ListSet<continuous_basic_set_type>::const_iterator bs_iter=list_set.begin();
        bs_iter!=list_set.end(); ++bs_iter)
    {
      const continuous_basic_set_type& bs=*bs_iter;
      working_sets.push_back(timed_set_type(0,0,loc_id,bs));
    }
  }
  return working_sets;
}


template<class R> inline
typename Evaluation::ConstraintHybridEvolver<R>::hybrid_list_set_type
Evaluation::ConstraintHybridEvolver<R>::_compute_list_set(const working_sets_type& working_sets, const HybridSpace& locations) const
{
  hybrid_list_set_type result_set(locations);
  for(typename working_sets_type::const_iterator iter=working_sets.begin();
      iter!=working_sets.end(); ++iter)
  {
    result_set.adjoin(hybrid_basic_set_type(iter->discrete_state(),iter->continuous_state_set()));
  }

    return result_set;
}


template<class R>
typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type
Evaluation::ConstraintHybridEvolver<R>::discrete_step(const transition_type& transition,
                                                      const timed_set_type& initial_set) const
{
  assert(transition.source_id()==initial_set.discrete_state());
  continuous_basic_set_type continuous_state_set = 
    this->_applicator->evaluate(transition.reset(),initial_set.continuous_state_set());
  return timed_set_type(initial_set.time(),initial_set.steps()+1,transition.destination_id(),continuous_state_set);
}


template<class R>
typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type
Evaluation::ConstraintHybridEvolver<R>::integration_step(const mode_type& mode,
                                                         const timed_set_type& initial_set,
                                                         const bounding_box_type& bounding_box,
                                                         const time_type& time_step) const
{
  assert(mode.id()==initial_set.discrete_state());
  continuous_basic_set_type continuous_state_set = 
    this->_integrator->bounded_integration_step(mode.dynamic(),initial_set.continuous_state_set(),bounding_box,time_step);
  return timed_set_type(time_type(initial_set.time()+time_step),initial_set.steps(),initial_set.discrete_state(),continuous_state_set);
}



/*! \brief Compute the possible states reached by an fored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::forced_jump(const transition_type& transition,
                                                    const timed_set_type& initial_set,
                                                    const bounding_box_type& source_bounding_box,
                                                    const time_type& required_time_step) const
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

  ARIADNE_LOG(2,"HybridEvolver::forced_jump(...)\n");
  ARIADNE_LOG(3,"  transition_set="<<transition<<"\n");
  ARIADNE_LOG(3,"  initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(3,"  source_bounding_box="<<source_bounding_box<<"\n");
  ARIADNE_LOG(3,"  time_step="<<required_time_step<<"\n");
  
  const System::VectorFieldInterface<R>& dynamic1=transition.source().dynamic();
  const System::VectorFieldInterface<R>& dynamic2=transition.destination().dynamic();
  const System::MapInterface<R>& reset=transition.reset();
  const Geometry::ConstraintInterface<R>& guard=transition.activation();
  Geometry::Zonotope<I,I> continuous_state_set=initial_set.continuous_state_set();
  time_type time_step1=required_time_step;
  time_type time_step2=required_time_step;

  continuous_state_set = orthogonal_over_approximation(continuous_state_set);

  Rectangle<R> bounding_box = source_bounding_box;

  // Estimate the time at which the centre crosses the guard set
  Point<I> centre = continuous_state_set.centre();
  Rectangle<R> centre_bounding_box = this->_integrator->refine_flow_bounds(dynamic1,Rectangle<R>(centre),bounding_box,required_time_step);
  I centre_normal_derivative = inner_product(guard.gradient(centre_bounding_box),dynamic1(centre_bounding_box));
  assert((bool)(centre_normal_derivative>0));
  I centre_crossing_time_interval = -guard.value(centre)/centre_normal_derivative;
  
  // Refine the estimate of the centre crossing time
  time_type minimum_centre_crossing_time = centre_crossing_time_interval.lower();
  time_type maximum_centre_crossing_time = centre_crossing_time_interval.upper();
  time_type estimated_centre_crossing_time = centre_crossing_time_interval.midpoint();
  centre_bounding_box = this->_integrator->refine_flow_bounds(dynamic1,Rectangle<R>(centre),centre_bounding_box,maximum_centre_crossing_time);
  Point<I> flowed_centre = this->_integrator->bounded_flow(dynamic1,centre,centre_bounding_box,minimum_centre_crossing_time);
  Rectangle<R> flowed_centre_bounding_box = this->_integrator->refine_flow_bounds(dynamic1,Rectangle<R>(flowed_centre),centre_bounding_box,time_type(maximum_centre_crossing_time-minimum_centre_crossing_time));
  I flowed_normal_derivative = inner_product(guard.gradient(flowed_centre_bounding_box),dynamic1(flowed_centre_bounding_box));
  I flowed_centre_crossing_time_interval = -guard.value(flowed_centre)/centre_normal_derivative;
  
  // Flow the set to the estimated crossing time interval
  estimated_centre_crossing_time = time_type(minimum_centre_crossing_time)+time_type(flowed_centre_crossing_time_interval.midpoint());
  maximum_centre_crossing_time = time_type(minimum_centre_crossing_time)+time_type(flowed_centre_crossing_time_interval.upper());
  bounding_box = this->_integrator->refine_flow_bounds(dynamic1,continuous_state_set.bounding_box(),bounding_box,estimated_centre_crossing_time);
  timed_set_type flowed_set = this->integration_step(transition.source(),initial_set,bounding_box,estimated_centre_crossing_time);

  this->_trace.push_back(flowed_set);

  // Estimate range of hitting times and the bounding box for flowed set
  I normal_derivative = inner_product(guard.gradient(bounding_box),dynamic1(bounding_box));
  assert((bool)(normal_derivative>0));
  I crossing_time_interval = -guard.value(bounding_box)/normal_derivative;
  Interval<Rational> flowed_crossing_time_interval = Interval<Rational>(crossing_time_interval) - estimated_centre_crossing_time;

  Rectangle<R> flowed_crossing_bounding_box = this->_integrator->refine_flow_bounds(dynamic1,flowed_set.bounding_box(),bounding_box,flowed_crossing_time_interval);
  flowed_crossing_time_interval = -guard.value(flowed_crossing_bounding_box)/inner_product(guard.gradient(flowed_crossing_bounding_box),dynamic1(flowed_crossing_bounding_box));
  flowed_crossing_bounding_box = this->_integrator->refine_flow_bounds(dynamic1,flowed_set.bounding_box(),bounding_box,flowed_crossing_time_interval);
  
  Rectangle<R> jump_crossing_bounding_box = this->_applicator->evaluate(reset,flowed_crossing_bounding_box);
  
  // Compute the gradient of the crossing times
  Point<I> source_centre = flowed_set.continuous_state_set().centre();
  Matrix<I> generators = flowed_set.continuous_state_set().generators();
  bounding_box = flowed_crossing_bounding_box;
  Rectangle<R> destination_bounding_box = flowed_crossing_bounding_box;
  Vector<I> guard_gradient = guard.gradient(source_bounding_box);
  Vector<I> time_gradient = (-guard_gradient)/inner_product(guard_gradient,dynamic1(source_bounding_box));
  
  // Compute the jump set
  Matrix<I> reset_jacobian = reset.jacobian(bounding_box);
  Matrix<I> jump_derivative = reset_jacobian + outer_product(reset_jacobian*dynamic1(bounding_box)-dynamic2(destination_bounding_box),time_gradient);
  Point<I> jump_centre = reset(source_centre);
  timed_set_type jump_set(flowed_set.time(),flowed_set.steps()+1,transition.destination_id(),Zonotope<I>(jump_centre,jump_derivative*generators));
  this->_trace.push_back(jump_set);

  Zonotope<I> integral_set = this->_integrator->integrate(dynamic2,jump_set.continuous_state_set(),required_time_step-estimated_centre_crossing_time);

  timed_set_type result_set(initial_set.time()+required_time_step,initial_set.steps()+1,
                            transition.destination_id(),integral_set);

  return std::vector<timed_set_type>(1,result_set);
}


/*! \brief Compute the possible states reached by an unfored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::unforced_jump(const transition_type& transition,
                                                      const timed_set_type& initial_set,
                                                      const bounding_box_type& source_flow_bounding_box,
                                                      const time_type& required_time_step) const
{
  ARIADNE_LOG(2,"HybridEvolver::unforced_jump(...)\n");
  ARIADNE_LOG(3,"  transition="<<transition<<"\n  initial_set="<<initial_set<<"\n  required_time_step="<<required_time_step<<"\n");  

  const System::VectorFieldInterface<R>& dynamic1=transition.source().dynamic();
  const System::VectorFieldInterface<R>& dynamic2=transition.destination().dynamic();
  const System::MapInterface<R>& reset=transition.reset();
  const Geometry::ConstraintInterface<R>& activation=transition.activation();
  Geometry::Zonotope<I,I> continuous_state_set=initial_set.continuous_state_set();

  std::vector<timed_set_type> result;

  const time_type finish_time = initial_set.time()+required_time_step;

  Rectangle<R> source_bounding_box,destination_bounding_box;

  // Find a bounding box and time step size for the flow in the destination mode
  
  timed_set_type current_set = initial_set;
  while(current_set.time() < finish_time) {
    time_type time_step = finish_time - current_set.time();
    
    destination_bounding_box = this->_integrator->estimate_flow_bounds(dynamic2,this->_applicator->evaluate(reset,source_flow_bounding_box),time_step);
    source_bounding_box = this->_integrator->refine_flow_bounds(dynamic1,current_set.bounding_box(),source_bounding_box,time_step);

    // Test if jump is activated
    if( possibly(activation.value(source_bounding_box)>=0) ) {
      // Jump  activated
      destination_bounding_box = this->_integrator->refine_flow_bounds(dynamic2,this->_applicator->evaluate(reset,source_flow_bounding_box),destination_bounding_box,time_step);
      timed_set_type jump_set = this->integration_step(transition.source(),current_set,source_bounding_box,time_step/2);
      jump_set = this->discrete_step(transition,jump_set);
      jump_set = this->integration_step(transition.destination(),jump_set,destination_bounding_box,time_step/2);
      
      Matrix<I> DPhi2 = this->_integrator->estimate_flow_jacobian_bounds(dynamic2,destination_bounding_box,time_step/2);
      Vector<I> F1 = dynamic1(source_bounding_box);
      Vector<I> F2 = dynamic2(destination_bounding_box);
      Matrix<I> Dr = reset.jacobian(source_bounding_box);

      I hh(time_type(time_step/2));
      Vector<I> correction =  hh * ( (DPhi2*(Dr*F1)) - F2 );

      const Point<I>& centre = jump_set.continuous_state_set().centre();
      const Matrix<I>& generators = jump_set.continuous_state_set().generators();
      
      Zonotope<I> reached_set(centre,generators,correction);
      
      assert(reached_set.number_of_generators()==current_set.continuous_state_set().number_of_generators() + 1);
      result.push_back(timed_set_type(current_set.time()+time_step,current_set.steps()+1,transition.destination_id(),reached_set));

    }
       
    current_set = this->integration_step(transition.source(),current_set,source_bounding_box,time_step);

  }
  return result;
  
}


template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                             const timed_set_type& initial_set,
                                                             time_type& time_step)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                             const timed_set_type& initial_set,
                                                             time_type& time_step)
{
  ARIADNE_LOG(2,"HybridEvolver::upper_evolution_step(...)\n");
  ARIADNE_LOG(3,"  initial_set="<<initial_set<<"\n");  

  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

      
  std::vector<timed_set_type> result;

  const id_type& discrete_state = initial_set.discrete_state();
  const continuous_basic_set_type& basic_set = initial_set.continuous_state_set();
  const time_type& time = initial_set.time();
  const discrete_time_type& steps = initial_set.steps();

  const mode_type& mode = automaton.mode(discrete_state);
  const System::VectorFieldInterface<R>& dynamic=mode.dynamic();

  const std::vector< constraint_const_pointer >& invariants=automaton.invariants(discrete_state);
  const std::map< id_type, constraint_const_pointer >& activations=automaton.activations(discrete_state);
  const std::map< id_type, constraint_const_pointer >& guards=automaton.guards(discrete_state);

  ARIADNE_LOG(4,"  time_step="<<time_step<<", time="<<time<<"\n");
  Rectangle<R> bounding_box=this->_integrator->estimate_flow_bounds(dynamic,basic_set.bounding_box(),time_step);
  ARIADNE_LOG(4,"  bounding_box="<<bounding_box<<"\n");
  uint old_verbosity=verbosity;
  verbosity=0;
  continuous_basic_set_type reachable_set = this->_integrator->bounded_reachability_step(dynamic,basic_set,bounding_box,time_step);
  verbosity=old_verbosity;
  ARIADNE_LOG(4,"  reachable_set="<<reachable_set<<"\n");
  continuous_basic_set_type integral_set = this->_integrator->bounded_integration_step(dynamic,basic_set,bounding_box,time_step);
  ARIADNE_LOG(4,"  integral_set="<<integral_set<<"\n");

  if(time_step < this->_integrator->minimum_step_size()) {
    throw std::runtime_error("Minimum step size reached -- aborting.");
  }

  bool block = false;
  std::vector<id_type> events;

  // Check for blocking by invariants
  ARIADNE_LOG(7,"invariants.size()="<<invariants.size()<<"\n");
  for(typename std::vector< constraint_const_pointer >::const_iterator constraint_ptr_iter = invariants.begin();
      constraint_ptr_iter != invariants.end(); ++constraint_ptr_iter)
  {
    const ConstraintInterface<R>& constraint=**constraint_ptr_iter;
    ARIADNE_LOG(7,"constraint_ptr="<<*constraint_ptr_iter<<"\n");
    ARIADNE_LOG(7,"constraint="<<constraint<<"\n");
    if(!satisfies(integral_set,constraint)) {
      ARIADNE_LOG(7,"  constraint not satisfied\n");
      block=true;
      break;
    }
  }

  // Check for blocking by guards (blocking if constraint is satisfied)
  ARIADNE_LOG(7,"guards.size()="<<guards.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    const id_type event_id=constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(possibly(satisfies(reachable_set,constraint))) {
      events.push_back(event_id);
      if(satisfies(integral_set,constraint)) {
        block=true;
      } 
    }
  }

  // Check for activation 
  ARIADNE_LOG(7,"activations.size()="<<activations.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = activations.begin();
      constraint_iter != activations.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(possibly(satisfies(reachable_set,constraint))) {
      events.push_back(event_id);
    }
  }

  ARIADNE_LOG(5,"block="<<block<<", events="<<events<<"\n");


  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = activations.begin();
      constraint_iter != activations.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(possibly(satisfies(reachable_set,constraint))) {
      break;
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      std::vector<timed_set_type> jump_set = this->unforced_jump(transition,initial_set,bounding_box,time_step);
      ::append(result,jump_set);
    }
  }

  // Run forced transitions
  for(uint i=0; i!=events.size(); ++i) {
    id_type event_id = events[i];
    const transition_type& transition = automaton.transition(event_id,discrete_state);
    if(transition.forced()) {
      std::vector<timed_set_type> jump_sets = this->forced_jump(transition,initial_set,bounding_box,time_step);
      ::append(result,jump_sets);
    } else {
      std::vector<timed_set_type> jump_sets = this->unforced_jump(transition,initial_set,bounding_box,time_step);
      ::append(result,jump_sets);
    }
  }

  if(block==false) {
    ARIADNE_LOG(7,"continuous_evolution="<<guards.size()<<"\n");
    result.push_back(timed_set_type(time+time_step,steps,discrete_state,integral_set));
    ARIADNE_LOG(7,"continuous_evolution="<<result.back()<<"\n");
  }

  if(!events.empty()) {
    ARIADNE_LOG(7,"result="<<result<<"\n");
    //assert(false);
  }
  return result;
}



template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                                const timed_set_type& initial_set,
                                                                time_type& time_step)
{
  ARIADNE_LOG(2,"HybridEvolver::upper_reachability_step(...)\n");
  ARIADNE_LOG(3,"  initial_set="<<initial_set<<" time_step="<<time_step<<"\n");  

  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;

      
  std::vector<timed_set_type> result;

  const id_type& discrete_state = initial_set.discrete_state();
  const continuous_basic_set_type& basic_set = initial_set.continuous_state_set();
  const time_type& time = initial_set.time();
  const discrete_time_type& steps = initial_set.steps();

  const mode_type& mode = automaton.mode(discrete_state);
  const System::VectorFieldInterface<R>& dynamic=mode.dynamic();

  const std::vector< constraint_const_pointer >& invariants=automaton.invariants(discrete_state);
  const std::map< id_type, constraint_const_pointer >& activations=automaton.activations(discrete_state);
  const std::map< id_type, constraint_const_pointer >& guards=automaton.guards(discrete_state);

  Rectangle<R> bounding_box=this->_integrator->estimate_flow_bounds(dynamic,basic_set.bounding_box(),time_step);
  bounding_box=this->_integrator->refine_flow_bounds(dynamic,basic_set.bounding_box(),bounding_box(),time_step);

  continuous_basic_set_type reach_set = this->_integrator->bounded_reachability_step(dynamic,basic_set,bounding_box,time_step);
  continuous_basic_set_type integral_set = this->_integrator->bounded_integration_step(dynamic,basic_set,bounding_box,time_step);
  ARIADNE_LOG(4,"time_step="<<time_step<<", time="<<time<<"\n");
  ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
  ARIADNE_LOG(4,"integral_set="<<integral_set<<"\n");

  if(time_step < this->_integrator->minimum_step_size()) {
    throw std::runtime_error("Minimum step size reached -- aborting.");
  }

  bool block = false;

  // Check for blocking by invariants
  ARIADNE_LOG(7,"invariants.size()="<<invariants.size()<<"\n");
  for(typename std::vector< constraint_const_pointer >::const_iterator constraint_ptr_iter = invariants.begin();
      constraint_ptr_iter != invariants.end(); ++constraint_ptr_iter)
  {
    const ConstraintInterface<R>& constraint=**constraint_ptr_iter;
    ARIADNE_LOG(7,"constraint="<<constraint<<"\n");
    if(!satisfies(integral_set,constraint)) {
      ARIADNE_LOG(7,"  invariant not satisfied\n");
      block=true;
      break;
    }
  }

  // Check for blocking by guards (blocking if constraint is satisfied)
  ARIADNE_LOG(7,"guards.size()="<<guards.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    ARIADNE_LOG(7,"constraint="<<constraint<<"\n");
    if(satisfies(integral_set,constraint)) {
      ARIADNE_LOG(7,"  guard not satisfied\n");
      block=true;
      break;
    }
  }

  if(block==false) {
    result.push_back(timed_set_type(time+time_step,steps,discrete_state,integral_set));
  }

  return result;

  // Check for activation 
  ARIADNE_LOG(7,"activations.size()="<<activations.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = activations.begin();
      constraint_iter != activations.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      std::vector<timed_set_type> jump_set = this->unforced_jump(transition,initial_set,bounding_box,time_step);
      ::append(result,jump_set);
    }

  }

  // Check for guards 
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      std::vector<timed_set_type> jump_sets = this->forced_jump(transition,initial_set,bounding_box,time_step);
      ::append(result,jump_sets);
    }

  }

  return result;
}



template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                     time_type evolution_time,
                                                     size_type maximum_number_of_events)
{
  verbosity=8;
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<", evolution_time="<<evolution_time<<"\n");
  ARIADNE_LOG(7,"maximum_step_size="<<this->maximum_step_size()<<", maximum_basic_set_radius="<<this->maximum_basic_set_radius()<<"\n");
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets=this->_compute_working_sets(initial_set);
  std::vector< timed_set_type > final_sets;
  std::vector< timed_set_type >& trace_sets=this->_trace;
  trace_sets=working_sets;

  while(!working_sets.empty()) {
    ARIADNE_LOG(9,"working_sets.size()="<<working_sets.size()<<"\n");
    timed_set_type working_set=working_sets.back(); working_sets.pop_back();
    time_type time_step=Numeric::min(this->maximum_step_size(),time_type(evolution_time-working_set.time()));
    ARIADNE_LOG(9,"working_set="<<working_set<<"\n");
    assert((bool)(working_set.time()<=evolution_time));
    
    if(working_set.time()==evolution_time) {
      ARIADNE_LOG(9,"  reached evolution time\n");
      final_sets.push_back(working_set);
    } else if(working_set.continuous_state_set().radius()>this->maximum_basic_set_radius()) {
      ::append(working_sets,this->subdivide(working_set));
    } else {
      time_step=Numeric::min(time_step,time_type(evolution_time-working_set.time()));
      ::append(working_sets,this->upper_evolution_step(automaton,working_set,time_step));
      ::append(this->_trace,this->upper_evolution_step(automaton,working_set,time_step));
    }
  }

  return this->_compute_list_set(final_sets,automaton.locations());
}


template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events)
{
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<", evolution_time="<<evolution_time<<"\n");
  ARIADNE_LOG(7,"maximum_step_size="<<this->maximum_step_size()<<", maximum_basic_set_radius="<<this->maximum_basic_set_radius()<<"\n");

  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets=this->_compute_working_sets(initial_set);
  std::vector< timed_set_type > reached_sets;

  while(!working_sets.empty()) {
    timed_set_type working_set=working_sets.back(); working_sets.pop_back();
    time_type time_step=Numeric::min(this->maximum_step_size(),time_type(evolution_time-working_set.time()));
    if(working_set.time()==evolution_time) {
    } else if(working_set.continuous_state_set().radius()>this->maximum_basic_set_radius()) {
      ::append(working_sets,this->subdivide(working_set));
    } else {
      ::append(working_sets,this->upper_evolution_step(automaton,working_set,time_step));
      ::append(reached_sets,this->upper_reachability_step(automaton,working_set,time_step));
      if(!working_sets.empty()) { std::cout << "evolution step: " << working_sets.back() << std::endl; }
      if(!reached_sets.empty()) { std::cout << "reach step: " << reached_sets.back() << std::endl; }
    }
  }

  return this->_compute_list_set(reached_sets,automaton.locations());
}



      

} // namespace Ariadne
