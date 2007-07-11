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

namespace Ariadne {
 
namespace Evaluation { static int& verbosity = hybrid_evolver_verbosity; }



template<class R>
Evaluation::ConstraintHybridEvolver<R>::~ConstraintHybridEvolver()
{
}

template<class R>
Evaluation::ConstraintHybridEvolver<R>::ConstraintHybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _applicator(&a), _integrator(dynamic_cast<LohnerIntegrator<R>*>(&i))
{
  if(!this->_integrator) {
    throw std::runtime_error("ConstraintHybridEvolver::ConstraintHybridEvolver(Applicator a, Integrator i): Invalid integrator");
  }
}



/*! \brief Compute the possible states reached by an unfored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
typename Evaluation::ConstraintHybridEvolver<R>::basic_set_type
Evaluation::ConstraintHybridEvolver<R>::unforced_jump(const System::VectorFieldInterface<R>& dynamic1,
                                                      const System::VectorFieldInterface<R>& dynamic2,
                                                      const System::MapInterface<R>& reset,
                                                      const basic_set_type& initial_set,
                                                      const Geometry::ConstraintInterface<R>& activation,
                                                      const time_type& step_size)
{
  typedef Numeric::Interval<R> I;

  time_type required_step_size=step_size;

  Geometry::Rectangle<R> bounding_box = initial_set.bounding_box();
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic1,bounding_box,required_step_size);
  Geometry::Point<I> mode1_bounds = bounding_box;
  bounding_box = this->_applicator->evaluate(reset,bounding_box);
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic2,bounding_box,required_step_size);
  Geometry::Point<I> mode2_bounds = bounding_box;
  assert(required_step_size==step_size);
  
  // Compute evolution assuming transition at time step_size/2
  time_type half_step_size=step_size/2;
  basic_set_type integral_set = initial_set;
  integral_set = this->_integrator->integration_step(dynamic1,integral_set,half_step_size);
  integral_set = this->_applicator->evaluate(reset,integral_set);
  integral_set = this->_integrator->integration_step(dynamic2,integral_set,half_step_size);
  
  // Compute extra generator associated to transition
  // FIXME: I don't think formula for jacobian is correct
  I hh=half_step_size;
  LinearAlgebra::Vector<I> v=(dynamic2.jacobian(mode2_bounds)*(reset.jacobian(mode1_bounds)*dynamic1.image(mode1_bounds))-dynamic2.image(mode2_bounds))*hh;

  return basic_set_type(integral_set.centre(),integral_set.generators(),v);
}


/*! \brief Compute the possible states reached by an unfored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
typename Evaluation::ConstraintHybridEvolver<R>::basic_set_type
Evaluation::ConstraintHybridEvolver<R>::forced_jump(const System::VectorFieldInterface<R>& dynamic1,
                                                    const System::VectorFieldInterface<R>& dynamic2,
                                                    const System::MapInterface<R>& reset,
                                                    const basic_set_type& initial_set,
                                                    const Geometry::ConstraintInterface<R>& guard,
                                                    const time_type& step_size)
{
  typedef Numeric::Interval<R> I;

  time_type required_step_size=step_size;

  Geometry::Rectangle<R> bounding_box = initial_set.bounding_box();
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic1,bounding_box,required_step_size);
  Geometry::Point<I> mode1_bounds = bounding_box;
  bounding_box = this->_applicator->evaluate(reset,bounding_box);
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic2,bounding_box,required_step_size);
  Geometry::Point<I> mode2_bounds = bounding_box;
  assert(required_step_size==step_size);
  
  // Compute time to jump as an interval affine map on the state

  // Reference centre and generators of zonotope
  const Geometry::Point<I>& c=initial_set.centre();
  const LinearAlgebra::Matrix<I>& G=initial_set.generators();
  
  Numeric::Interval<R> cval=guard.value(c);
  LinearAlgebra::Vector<I> cgrad=guard.gradient(mode1_bounds);
  
  // Direction of the vector field
  LinearAlgebra::Vector<I> dir=dynamic1(mode1_bounds);

  // Normal direction of the vector field
  Numeric::Interval<R> ndir=LinearAlgebra::inner_product(cgrad,dir);

  // Time to jump for centre of zonoope
  Numeric::Interval<R> tcentre=(cval-LinearAlgebra::inner_product(cgrad,c.position_vector()))/ndir;
  // Modification to jump time for generators of zonotope
  LinearAlgebra::Vector<I> tgrad=(cgrad*G)/(-ndir);
  
  // FIXME: need to use time interval for centre
  time_type centre_step_size1=(time_type(tcentre.upper())-time_type(tcentre.lower()))/2;
  time_type centre_step_size2=step_size-centre_step_size1;

  basic_set_type integral_set = initial_set;
  integral_set = this->_integrator->integration_step(dynamic1,integral_set,centre_step_size1);
  integral_set = this->_applicator->evaluate(reset,integral_set);
  integral_set = this->_integrator->integration_step(dynamic2,integral_set,centre_step_size2);
  
  // Compute extra generator associated to transition
  // FIXME: I don't think formula for jacobian is correct
  LinearAlgebra::Vector<I> v(dynamic2.jacobian(mode2_bounds)*(reset.jacobian(mode1_bounds)*dynamic1.image(mode1_bounds))-dynamic2.image(mode2_bounds));
   

  return basic_set_type(integral_set.centre(),integral_set.generators()+LinearAlgebra::outer_product(v,tgrad));
}


template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                             const timed_set_type& initial_set,
                                                             time_type& step_size)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                             const timed_set_type& initial_set,
                                                             time_type& step_size)
{
      
  std::vector<timed_set_type> result;

  const id_type& discrete_state = initial_set.set().discrete_state();
  const basic_set_type& basic_set = initial_set.set().continuous_state_set();
  const time_type& time = initial_set.time();

  const mode_type& mode = automaton.mode(discrete_state);
  const System::VectorFieldInterface<R>& dynamic=mode.dynamic();

  const std::vector< constraint_pointer >& invariants=automaton.invariants(discrete_state);
  const std::map< id_type, constraint_pointer >& activations=automaton.activations(discrete_state);
  const std::map< id_type, constraint_pointer >& guards=automaton.guards(discrete_state);


  basic_set_type reach_set = this->_integrator->reachability_step(dynamic,basic_set,step_size);
  time_type reach_step_size=step_size;
  basic_set_type integral_set = this->_integrator->integration_step(dynamic,basic_set,step_size);
  // FIXME: The step size for a reach step should always be less than that of an integration step.
  while(reach_step_size!=step_size) {
    ARIADNE_LOG(4,"reach_step_size="<<reach_step_size<<", integrate_step_size="<<step_size<<"\n");
    basic_set_type reach_set = this->_integrator->reachability_step(dynamic,basic_set,step_size);
    reach_step_size=step_size;
    basic_set_type integral_set = this->_integrator->integration_step(dynamic,basic_set,step_size);
  }
  ARIADNE_LOG(4,"reach_step_size="<<reach_step_size<<", integrate_step_size="<<step_size<<"\n");
  ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
  ARIADNE_LOG(4,"integral_set="<<integral_set<<"\n");
  assert(reach_step_size==step_size);

  bool block = false;

  // Check for blocking by invariants
  ARIADNE_LOG(7,"invariants.size()="<<invariants.size()<<"\n");
  for(typename std::vector< constraint_pointer >::const_iterator constraint_ptr_iter = invariants.begin();
      constraint_ptr_iter != invariants.end(); ++constraint_ptr_iter)
  {
    const Geometry::ConstraintInterface<R>& constraint=**constraint_ptr_iter;
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
  for(typename std::map< id_type, constraint_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    const Geometry::ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(satisfies(integral_set,constraint)) {
      block=true;
      break;
    }
  }

  if(block==false) {
    result.push_back(timed_set_type(time+step_size,hybrid_basic_set_type(discrete_state,integral_set)));
  }


  // Check for activation 
  ARIADNE_LOG(7,"activations.size()="<<activations.size()<<"\n");
  for(typename std::map< id_type, constraint_pointer >::const_iterator constraint_iter = activations.begin();
      constraint_iter != activations.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const Geometry::ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      const mode_type& destination = transition.destination();
      id_type destination_state = destination.id();
      const System::VectorFieldInterface<R>& destination_dynamic = destination.dynamic();
      const System::MapInterface<R>& reset = transition.reset();
      basic_set_type jump_set = this->unforced_jump(dynamic,destination_dynamic,reset,basic_set,constraint,step_size);
      result.push_back(timed_set_type(time+step_size,destination_state,jump_set));
      assert(false); // FINISH
    }

  }

  // Check for guards 
  for(typename std::map< id_type, constraint_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const Geometry::ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      const mode_type& destination = transition.destination();
      id_type destination_state = destination.id();
      const System::VectorFieldInterface<R>& destination_dynamic = destination.dynamic();
      const System::MapInterface<R>& reset = transition.reset();
      basic_set_type jump_set = this->forced_jump(dynamic,destination_dynamic,reset,basic_set,constraint,step_size);
      result.push_back(timed_set_type(time+step_size,destination_state,jump_set));
    }

  }

  return result;
}



template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                     const hybrid_basic_set_type& initial_set, 
                                                     time_type evolution_time,
                                                     size_type maximum_number_of_events)
{
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set.continuous_state_set()<<"\n");
  using namespace Numeric;
  
  ARIADNE_LOG(7,"integrator="<<this->_integrator);
  time_type step_size=this->_integrator->maximum_step_size();
  ARIADNE_LOG(7,"step_size=");
  R maximum_set_radius=this->_integrator->maximum_basic_set_radius();
  ARIADNE_LOG(7,"step_size="<<step_size<<", maximum_set_radius="<<maximum_set_radius<<"\n");
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  Geometry::HybridListSet< basic_set_type > final_set(automaton.locations());

  working_sets.push_back(timed_set_type(t,initial_set));
  
  typedef typename std::vector< timed_set_type >::const_iterator list_set_const_iterator;
  while(!working_sets.empty()) {
    timed_set_type timed_set=working_sets.back();
    working_sets.pop_back();
    
    if(verbosity>5) { 
    }

    const hybrid_basic_set_type& hbs=timed_set.set();
    if(timed_set.time()==evolution_time) {
      final_set.adjoin(hbs.discrete_state(),hbs.continuous_state_set());
    } else if(hbs.radius()>this->_integrator->maximum_basic_set_radius()) {
      assert(false);
    } else {
      std::vector<timed_set_type> new_sets=this->upper_evolution_step(automaton,timed_set,h);
      for(typename std::vector<timed_set_type>::const_iterator set_iterator=new_sets.begin();
          set_iterator!=new_sets.end(); ++set_iterator)
      {
        working_sets.push_back(*set_iterator);
      } 
      new_sets.clear();
    }

  }

  //if(verbosity>6) { std::clog << "  final_set=" << final_set << std::endl; }
  return final_set;
}



      

} // namespace Ariadne
