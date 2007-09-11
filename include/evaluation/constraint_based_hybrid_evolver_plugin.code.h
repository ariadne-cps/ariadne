/***************************************************************************
 *            constraint_based_hybrid_evolver_plugin.code.h
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
#include <cassert>

#include "../base/reference_container.h"
#include "../linear_algebra/vector.h"
#include "../geometry/rectangle_expression.h"
#include "../geometry/set_interface.h"
#include "../geometry/hybrid_set.h"
#include "../geometry/timed_set.h"
#include "../system/map.h"
#include "../system/vector_field.h"
#include "../system/constraint_based_hybrid_automaton.h"
#include "../evaluation/applicator_plugin_interface.h"
#include "../evaluation/bounder_plugin_interface.h"
#include "../evaluation/integrator_plugin_interface.h"
#include "../evaluation/detector_plugin_interface.h"

#include "../evaluation/bounder_plugin.h"
#include "../evaluation/detector_plugin.h"
#include "../evaluation/lohner_integrator.h"

#include "../output/epsstream.h"
#include "../output/logging.h"

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

template<class T>
std::set<id_type> ids(const reference_set< T >& transitions) {
  std::set<id_type> result;
  for(typename reference_set< T >::const_iterator transition_iter=transitions.begin();
      transition_iter!=transitions.end(); ++transition_iter)
  {
    result.insert(transition_iter->id()); 
  }
  return result;
}

}


namespace Ariadne {

namespace Evaluation { 
static int& verbosity = hybrid_evolver_verbosity; 
static const uint time_step_event=uint(-1);
static const uint final_time_event=uint(-2);
static const uint regularizing_time_event=uint(-3);
}
  
using namespace Numeric;
using namespace LinearAlgebra;
using namespace Geometry;
using namespace System;


template<class R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::~ConstraintBasedHybridEvolverPlugin()
{
  delete _applicator;
  delete _bounder;
  delete _integrator;
  delete _detector;
}


template<class R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::ConstraintBasedHybridEvolverPlugin(const ApplicatorPluginInterface<BS>& a, const IntegratorPluginInterface<BS>& i, const DetectorPluginInterface<R>& d)
  : _applicator(a.clone()), 
    _bounder(new BounderPlugin<R>()),
    _integrator(dynamic_cast<C1LohnerIntegrator<R>*>(i.clone())),
    _detector(d.clone())
{
  if(!this->_integrator) {
    throw std::runtime_error("ConstraintBasedHybridEvolverPlugin::ConstraintBasedHybridEvolverPlugin(MapEvolver a, VectorFieldEvolver i, Detector d): Invalid integrator");
  }
}


template<class R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::ConstraintBasedHybridEvolverPlugin(const ConstraintBasedHybridEvolverPlugin<R>& plugin)
  : _applicator(plugin._applicator->clone()), 
    _bounder(plugin._bounder->clone()), 
    _integrator(plugin._integrator->clone()),
    _detector(plugin._detector->clone())
{
}




template<class R> inline
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::regularize(const timed_set_type& timed_set) const
{
  ARIADNE_LOG(8,"ConstraintBasedHybridEvolverPlugin::regularize(...)\n");
  ARIADNE_LOG(9,"  radius="<<timed_set.continuous_state_set().radius()<<"\n");
  ARIADNE_LOG(9,"  bounding_box="<<timed_set.continuous_state_set().bounding_box()<<"\n");
  ARIADNE_LOG(9,"  set="<<timed_set<<"\n");

  if(possibly(timed_set.time().gradient()!=Vector<I>(timed_set.time().gradient().size()))) {
    ARIADNE_LOG(9,"  time="<<timed_set.time()<<"\n");
    ARIADNE_LOG(9,"  cannot regularize if time is not constant\n");
    return timed_set;
  }
  return timed_set;
  //if(Geometry::Rectangle<R>(timed_set.centre()).radius()*2>timed_set.radius()) {
  if(true) {
    Zonotope<R> z=orthogonal_over_approximation(timed_set.continuous_state_set());
    time_model_type t(timed_set.time().average(),LinearAlgebra::Vector<I>(z.number_of_generators()));
    std::cerr << t << "\n" << z << std::endl;
    assert(t.number_of_generators()==z.number_of_generators());
    return timed_set_type(t,timed_set.steps(),timed_set.discrete_state(),Zonotope<I>(z));
  } else {
    return timed_set;
  }

}

template<class R> inline
std::vector<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::subdivide(const timed_set_type& timed_set) const
{
  using namespace LinearAlgebra;
  
  ARIADNE_LOG(8,"ConstraintBasedHybridEvolverPlugin:subdivide(...)\n");
  ARIADNE_LOG(9,"  radius="<<timed_set.continuous_state_set().radius()<<"\n");
  ARIADNE_LOG(9,"  bounding_box="<<timed_set.continuous_state_set().bounding_box()<<"\n");
  ARIADNE_LOG(9,"  set="<<timed_set<<"\n");
  std::vector<timed_set_type> result;
  const Zonotope<I,I>& z=timed_set.continuous_state_set();
  const TimeModel<R>& t=timed_set.time();
  assert(t.number_of_generators()==z.number_of_generators());

  dimension_type dimension=z.dimension();
  size_type number_of_generators=z.number_of_generators();
  
  // Combine generators and time
  
  Zonotope<I,I> combined_set(Point<I>(concatenate(z.centre().position_vector(),t.average())),concatenate_rows(z.generators(),t.gradient()));
  
  ListSet< Zonotope<I,I> > subdivisions=Geometry::divide(combined_set);
  //ListSet< Zonotope<I,R> > subdivisions=Geometry::subdivide(orthogonal_approximation);
  for(typename ListSet< Zonotope<I,I> >::const_iterator iter=subdivisions.begin();
      iter!=subdivisions.end(); ++iter)
  {
    Zonotope<I,I> combined_set=*iter;
    Zonotope<I,I> continuous_state_set(Point<I>(VectorSlice<const I>(dimension,combined_set.centre().position_vector().begin())),
                                       Matrix<I>(MatrixSlice<const I>(dimension,number_of_generators,combined_set.generators().begin())));
    TimeModel<R> time_model(combined_set.centre()[dimension],combined_set.generators().row(dimension));
    
    if(time_model.number_of_generators()!=continuous_state_set.number_of_generators()) {
      ARIADNE_LOG(2,"ConstraintBasedHybridEvolverPlugin:subdivide(...):\n"
                    "combined_set="<<*iter<<", continuous_state_set="<<continuous_state_set<<", time_model="<<time_model<<"\n");
      assert(time_model.number_of_generators()==continuous_state_set.number_of_generators());
    }
    result.push_back(timed_set_type(time_model,timed_set.steps(),timed_set.discrete_state(),continuous_state_set));
  }
  return result;
}


template<class R>
Geometry::Rectangle<R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::estimate_flow_bounds(const mode_type& mode, 
                                                                   const timed_set_type& initial_set,
                                                                   time_type& maximum_step_size) const
{
  return this->_bounder->estimate_flow_bounds(mode.dynamic(),initial_set.bounding_box(),maximum_step_size);
}

template<class R>
Geometry::Rectangle<R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::refine_flow_bounds(const mode_type& mode, 
                                                                 const timed_set_type& initial_set, 
                                                                 const bounding_box_type& bounding_set,
                                                                 const time_type& maximum_step_size) const
{
  return this->_bounder->refine_flow_bounds(mode.dynamic(),initial_set.bounding_box(),bounding_set,maximum_step_size);
}






template<class R>
tribool
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::enabled(const transition_type& transition,
                                                           const bounding_box_type& bounding_box) const
{
  return this->_detector->value(transition.constraint(),bounding_box) >= 0;
}


template<class R>
tribool
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::enabled(const transition_type& transition,
                                                           const timed_set_type& set) const
{
  DetectorPlugin<R>* detector=dynamic_cast<DetectorPlugin<R>*>(this->_detector);
  assert(detector);
  return detector->value(transition.constraint(),set.continuous_state_set()) >= 0;
}


/*! Tests if constraint2 is always unsatisfied if constraint1 is unsatisfied. */
template<class R>
tribool
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::forces(const constraint_type& constraint1,
                                                          const constraint_type& constraint2,
                                                          const bounding_box_type& bounding_box) const
{
  return this->_detector->forces(constraint1,constraint2,bounding_box);
}


/*! Computes an estimation of the \em absolute crossing time step time, using the same base variables as the initial set current time.
 */
template<class R>
Numeric::Interval<R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::estimate_normal_derivative(const transition_type& transition,
                                                                              const bounding_box_type& bounding_box) const
{
  const vector_field_type& dynamic=transition.source().dynamic();
  const differentiable_constraint_type& constraint=dynamic_cast<const differentiable_constraint_type&>( transition.constraint());
  const Geometry::Point<I> bounds=bounding_box;
  assert(&constraint);

  return this->_detector->normal_derivative(dynamic,constraint,bounds);
}


/*! Computes an estimation of the \em absolute crossing time step time, using the same base variables as the initial set current time.
 */
template<class R>
Numeric::Interval<R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::estimate_crossing_time_step(const transition_type& transition, 
                                                                               const timed_set_type& initial_set,
                                                                               const bounding_box_type& bounding_box) const
{
  const vector_field_type& dynamic=transition.source().dynamic();
  const constraint_type& constraint=transition.constraint();
  const Geometry::Point<I> initial_point=initial_set.bounding_box();

  return this->_detector->crossing_time(dynamic,constraint,initial_point,bounding_box);
}
  

/*! Computes the \em absolute crossing time step time, using the same base variables as the initial set current time.
 *  Detects non-transverse crossings and provides first-order estimates
 */
template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::time_model_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::compute_crossing_time_step(const transition_type& transition, 
                                                               const timed_set_type& initial_set,
                                                               const bounding_box_type& bounding_box) const
{
  ARIADNE_LOG(8,"    ConstraintBasedHybridEvolverPlugin::compute_crossing_time_step(...)\n");
  const System::VectorFieldInterface<R>& dynamic=transition.source().dynamic();
  const Geometry::ConstraintInterface<R>& constraint=transition.constraint();
  const Geometry::Rectangle<R> domain=initial_set.bounding_box();

  TimeModel<R> spacial_crossing_time_step=this->_detector->crossing_time(dynamic,constraint,domain,bounding_box);
  ARIADNE_LOG(9,"    spacial_crossing_time_step="<<spacial_crossing_time_step<<"\n");

  // FIXME: Is this formula correct? Maybe the average crossing time is too small for the centre...
  TimeModel<R> parameter_crossing_time_step(spacial_crossing_time_step.average(),spacial_crossing_time_step.gradient()*initial_set.generators());
  ARIADNE_LOG(9,"    parameter_crossing_time_step="<<parameter_crossing_time_step<<"\n\n");

  return parameter_crossing_time_step;
}

  


                                                            
template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::crossing_data_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::compute_crossing_data(const transition_type& transition,
                                                                    const timed_set_type& initial_set,
                                                                    const timed_set_type& final_set,
                                                                    const bounding_box_type& bounding_box) const
{
  DetectorPlugin<R>* detector=dynamic_cast<DetectorPlugin<R>*>(this->_detector);
  const constraint_type& constraint = transition.constraint();

  // constraint might not be satisfied; compute further
  I value_range = this->_detector->value(constraint,bounding_box);
  I initial_constraint_value = detector->value(constraint,initial_set.continuous_state_set());
  I final_constraint_value = detector->value(constraint,final_set);
  I normal_derivative = this->estimate_normal_derivative(transition,bounding_box);
  I crossing_time_bounds = this->estimate_crossing_time_step(transition,initial_set,bounding_box);
  TimeModel<R> crossing_time_model=this->compute_crossing_time_step(transition,initial_set,bounding_box);
  int direction = ( normal_derivative>0 ? 1 : (normal_derivative < 0 ? -1 : 0) );
  ARIADNE_LOG(7,"     initial_constraint_value="<<initial_constraint_value<<"\n");
  ARIADNE_LOG(7,"     final_constraint_value="<<final_constraint_value<<"\n");
  ARIADNE_LOG(7,"     normal_derivative="<<normal_derivative<<"\n");
  ARIADNE_LOG(7,"     crossing_time_bounds="<<crossing_time_bounds<<"\n");
  ARIADNE_LOG(7,"     crossing_time_model="<<crossing_time_model<<"\n");
  
  CrossingData<R> crossing_data = { initial_constraint_value, final_constraint_value, 
                                    normal_derivative, crossing_time_bounds, 
                                    crossing_time_model, direction };
  return crossing_data;
}



template<class R>
std::map< id_type, typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::crossing_data_type >
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::compute_crossing_data(const reference_set<const transition_type>& transitions,
                                                                    const timed_set_type& initial_set,
                                                                    const timed_set_type& final_set,
                                                                    const bounding_box_type& bounding_box,
                                                                    const time_type& maximum_time_step) const
{
  ARIADNE_LOG(6," ConstraintBasedHybridEvolverPlugin::compute_crossing_data(...)\n");
  ARIADNE_LOG(7,"   initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(7,"   final_set="<<final_set<<"\n");
  ARIADNE_LOG(7,"   bounding_box="<<bounding_box<<"\n\n");
  
  std::map<id_type, crossing_data_type> result;
  
  for(typename reference_set<const transition_type>::const_iterator transition_iter=transitions.begin();
      transition_iter != transitions.end(); ++transition_iter)
  {
    const transition_type& transition = *transition_iter;
    const id_type& id = transition_iter->id();
    const constraint_type& constraint = transition_iter->constraint();
    
    I value_range=this->_detector->value(constraint,bounding_box);
    ARIADNE_LOG(7,"   constraint["<<id<<"]={"<<constraint<<"}\n");
    if(value_range<0) {
      // constraint is satisfied throughout evolution
      ARIADNE_LOG(7,"     value_range="<<value_range<<"\n");
    } else {
      crossing_data_type crossing_data = this->compute_crossing_data(transition,initial_set,final_set,bounding_box);
      // Remove events which are obviously not possible
      if( (crossing_data.crossing_time_bounds<0 || crossing_data.crossing_time_bounds>maximum_time_step)
          && (crossing_data.initial_constraint_value<0 || crossing_data.final_constraint_value<0) )
      {
        // Don't insert event as it cannot occur
      } else {
        result.insert(std::make_pair(transition.id(),crossing_data));
      }
    }
  }
  ARIADNE_LOG(7,"\n");

  return result;
}
    

template<class R>
Numeric::Interval<R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::
compute_evolution_time_bounds(const std::map<id_type, crossing_data_type>& crossing_data) const
{
  typedef typename std::map<id_type,crossing_data_type>::const_iterator crossing_data_const_iterator;
  Numeric::Interval<R> time_bounds = crossing_data.begin()->second.crossing_time_bounds;
  for(crossing_data_const_iterator cd_iter=crossing_data.begin(); cd_iter!=crossing_data.end(); ++cd_iter) {
    time_bounds=Numeric::min(time_bounds,cd_iter->second.crossing_time_bounds);
  }
  return time_bounds;
}


template<class R>
Numeric::Interval<R>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::
compute_evolution_time_bounds(const std::map<id_type, time_model_type>& event_times) const
{
  typedef typename std::map<id_type,time_model_type>::const_iterator const_iterator;
  Numeric::Interval<R> time_bounds = event_times.begin()->second.bound();
  for(const_iterator iter=event_times.begin(); iter!=event_times.end(); ++iter) {
    time_bounds=Numeric::min(time_bounds,iter->second.bound());
  }
  return time_bounds;
}


template<class R>
std::map< id_type, typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::time_model_type>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::compute_terminating_event_times(const mode_type& mode,
                                                                              const timed_set_type& initial_set,
                                                                              const timed_set_type& final_set,
                                                                              const time_type& final_time,
                                                                              const bounding_box_type& bounding_box,
                                                                              const time_type& time_step_size) const
{
  typedef boost::shared_ptr< const Geometry::ConstraintInterface<R> > constraint_const_pointer;
  typedef typename reference_set<const transition_type>::const_iterator transitions_const_iterator;
  typedef typename std::set< id_type >::const_iterator events_const_iterator;
  typedef typename std::map< id_type, time_model_type >::const_iterator event_times_const_iterator;
  typedef typename std::map< id_type, crossing_data_type>::const_iterator crossings_const_iterator;
  typedef typename std::map< id_type, time_model_type>::const_iterator crossing_times_const_iterator;
  typedef Numeric::Interval<R> time_interval_type;

  const reference_set<const transition_type>& transitions=mode.transitions();

  time_model_type zero_time_step(I(0),Vector<I>(initial_set.time().gradient().size()));
  
  // Compute the integration time step
  time_type initial_time_lower_bound=initial_set.time().bound().lower();
  time_type initial_time_upper_bound=initial_set.time().bound().upper();
  time_type initial_time_variation=initial_set.time().bound().width();
  time_type regularizing_time=initial_time_lower_bound+time_step_size;
  time_model_type integration_time_step=time_step_size+zero_time_step;
  time_model_type final_integration_time_step=final_time-initial_set.time();
  time_model_type regularizing_integration_time_step=regularizing_time-initial_set.time();
  

  // Compute events which are possibly enabled based on bounding box
  reference_set<const transition_type> possibly_enabled_transitions;
  for(transitions_const_iterator transition_iter=transitions.begin();
      transition_iter!=transitions.end(); ++transition_iter)
  {
    if(transition_iter->kind()==System::invariant_tag || transition_iter->kind()==System::guard_tag) {
      if(possibly(this->enabled(*transition_iter,bounding_box))) {
        possibly_enabled_transitions.insert(*transition_iter);
        ARIADNE_LOG(6,"  event "<<transition_iter->id()<<" with constraint "<<transition_iter->constraint()<<" possibly not satisfied by "<<bounding_box<<"\n");
      }
    }
  }
  
  // Compute full crossing data
  std::map<id_type,crossing_data_type> crossing_data
    = this->compute_crossing_data(possibly_enabled_transitions,initial_set,final_set,bounding_box,time_step_size);
  ARIADNE_LOG(6,"  crossing_data="<<crossing_data<<"\n");


  // Compute the times of the events
  std::map<id_type,time_model_type> event_time_steps;
  for(crossings_const_iterator crossing_iter = crossing_data.begin();
      crossing_iter != crossing_data.end(); ++crossing_iter)
  {
    const id_type& id = crossing_iter->first;
    const time_model_type& crossing_time = crossing_iter->second.crossing_time_model;
    event_time_steps.insert(std::make_pair(id,crossing_time));
  }
  // Add the integration termination times
  event_time_steps.insert(std::make_pair(time_step_event,integration_time_step));
  event_time_steps.insert(std::make_pair(final_time_event,final_integration_time_step));
  ARIADNE_LOG(6,"  event_time_steps="<<event_time_steps<<"\n");
 



  // Find minimal blocking times
  std::set<id_type> blocking_events;
  for(event_times_const_iterator event_iter = event_time_steps.begin();
      event_iter != event_time_steps.end(); ++event_iter)
  {
    const time_model_type& crossing_time = event_iter->second;
    bool possible_blocking_transition=true;
    for(events_const_iterator comparison_iter = blocking_events.begin();
        comparison_iter != blocking_events.end(); ++comparison_iter)
    {
      ARIADNE_LOG(8,"event=" << event_iter->first << " compare=" << *comparison_iter);
      const time_model_type& comparison_time = event_time_steps.find(*comparison_iter)->second;
      ARIADNE_LOG(8,"comparison_time=" << comparison_time << "\n");
      if(crossing_time<=comparison_time) {
        blocking_events.erase(comparison_iter);
      } else if(comparison_time<=crossing_time) {
        possible_blocking_transition=false;
        break;
      }
    }
    if(possible_blocking_transition) {
      blocking_events.insert(event_iter->first);
    }
  }
  ARIADNE_LOG(6,"  blocking_events=" << blocking_events << "\n");

  std::map<id_type,time_model_type> terminating_event_times;
  for(events_const_iterator event_iter=blocking_events.begin();
      event_iter!=blocking_events.end(); ++event_iter)
  {
    ARIADNE_LOG(8,"event=" << *event_iter << "\n");
    terminating_event_times.insert(std::make_pair(*event_iter,event_time_steps.find(*event_iter)->second));
  }

  ARIADNE_LOG(6,"  terminating_event_time_steps=" << terminating_event_times<<"\n");
  
  if(terminating_event_times.size()==1 && terminating_event_times.begin()->first==time_step_event) {
    terminating_event_times.clear();
    terminating_event_times.insert(std::make_pair(regularizing_time_event,regularizing_integration_time_step));
    ARIADNE_LOG(6,"  terminating_event_time_steps="<<terminating_event_times<<"\n");
  }
  return terminating_event_times;
}
      




template<class R>
reference_set<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::transition_type>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::compute_possibly_enabled_transitions(const mode_type& mode,
                                                                                   const bounding_box_type& bounding_box) const
{
  reference_set<transition_type> possibly_enabled_transitions;
  for(typename reference_set<const transition_type>::const_iterator transition_iter=mode.transitions().begin();
      transition_iter != mode.transitions().end(); ++transition_iter)
  {
    I value_range=this->_detector->value(transition_iter->constraint(),bounding_box);
    ARIADNE_LOG(7,"   constraint["<<transition_iter->id()<<"]={"<<transition_iter->constraint()<<"}\n");
    ARIADNE_LOG(7,"     value_range="<<value_range<<"\n");
    if(possibly(value_range>=0)) {
      possibly_enabled_transitions.insert(*transition_iter);
    }
  }
  ARIADNE_LOG(7,"\n");

  return possibly_enabled_transitions;
}








template<class R>
std::map<id_type, typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::time_model_pair_type>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::compute_enabled_activation_times(const mode_type& mode,
                                                                               const timed_set_type& initial_set,
                                                                               const timed_set_type& final_set,
                                                                               const time_model_type& maximum_time_step,
                                                                               const bounding_box_type& bounding_box) const
{
  typedef typename reference_set<const transition_type>::const_iterator transitions_const_iterator;
  typedef typename std::map<id_type, time_model_pair_type>::const_iterator crossings_const_iterator;
  typedef typename std::set<id_type>::const_iterator events_const_iterator;

  // Aliases for commonly used objects
  const reference_set<const transition_type>& transitions=mode.transitions();

  const time_model_type minimum_time_step(I(0),Vector<I>(maximum_time_step.gradient().size()));

  // Compute events which are possibly enabled based on bounding box
  reference_set<const transition_type> possibly_enabled_activations;
  for(transitions_const_iterator transition_iter=transitions.begin();
      transition_iter!=transitions.end(); ++transition_iter)
  {
    if(transition_iter->kind()==System::activation_tag) {
      if(possibly(this->enabled(*transition_iter,bounding_box))) {
        possibly_enabled_activations.insert(*transition_iter);
        ARIADNE_LOG(6,"  event "<<transition_iter->id()<<" with constraint "<<transition_iter->constraint()<<" possibly not satisfied by "<<bounding_box<<"\n");
      }
    }
  }

  std::map< id_type, std::pair<time_model_type,time_model_type> > enabled_activations;
  for(transitions_const_iterator transition_iter = possibly_enabled_activations.begin();
      transition_iter!=possibly_enabled_activations.end(); ++transition_iter)
  {
    const transition_type& transition=*transition_iter;
    const id_type event_id = transition.id();
    const differentiable_constraint_type& constraint=dynamic_cast<const differentiable_constraint_type&>(transition.constraint());
    assert(&constraint);
    const vector_field_type& dynamic=transition.source().dynamic();
    time_model_type crossing_time_step=this->compute_crossing_time_step(transition,initial_set,bounding_box);
    I normal_derivative = LinearAlgebra::inner_product(dynamic(bounding_box),constraint.gradient(bounding_box));
    ARIADNE_LOG(7,"    event "<<transition_iter->id()<<" crossing time step is "<<crossing_time_step<<"\n");
  
    tribool initially_enabled=this->enabled(transition,initial_set);
    tribool finally_enabled=this->enabled(transition,final_set);
    if(minimum_time_step >= crossing_time_step || crossing_time_step >= maximum_time_step) {
      if(definitely(initially_enabled)) {
        ARIADNE_LOG(7,"      activated entire evolution time\n");
        enabled_activations.insert(std::make_pair(event_id,std::make_pair(minimum_time_step,maximum_time_step)));
      } else if(possibly(initially_enabled)) {
        ARIADNE_LOG(7,"      appears to be activated entire evolution time\n");
        enabled_activations.insert(std::make_pair(event_id,std::make_pair(minimum_time_step,maximum_time_step)));
      } else {
        ARIADNE_LOG(7,"      inactivated entire evolution time\n");
      }
    } else if(not initially_enabled) {
      ARIADNE_LOG(7,"      event inactive at beginning of evolution step\n")
      enabled_activations.insert(std::make_pair(event_id,std::make_pair(crossing_time_step,maximum_time_step)));
    } else if(not finally_enabled) { 
      ARIADNE_LOG(7,"      event inactive at end end of evolution step\n");
      enabled_activations.insert(std::make_pair(event_id,std::make_pair(minimum_time_step,crossing_time_step)));
    } else {
      if(normal_derivative>0) {
        ARIADNE_LOG(7,"      event boundary traversed in positive direction; event becomes activated\n");
        enabled_activations.insert(std::make_pair(event_id,std::make_pair(crossing_time_step,maximum_time_step)));
      } else if(normal_derivative<0) {
        ARIADNE_LOG(7,"      event boundary traversed in negative direction; event becomes deactivated\n");
        enabled_activations.insert(std::make_pair(event_id,std::make_pair(minimum_time_step,crossing_time_step)));
      } else {
        ARIADNE_LOG(7,"  Warning: normal derivative to event boundary has indeterminate sign; cannot determine activation time\n");
        throw std::runtime_error("Normal derivative to event boundary has indeterminate sign; cannot determine activation time");
      }

    }
  }
  ARIADNE_LOG(6,"  activation_times"<<enabled_activations<<"\n\n");
  return enabled_activations;
}



template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::final_continuous_evolution_step(const mode_type& mode,
                                                                              const timed_set_type& initial_set,
                                                                              const time_type& final_time,
                                                                              const bounding_box_type& bounding_box) const
{
  ARIADNE_LOG(7," ConstraintBasedHybridEvolverPlugin::final_continuous_evolution_step(...)\n");
  ARIADNE_LOG(8,"    dynamic="<<mode.dynamic()<<"\n    initial_set="<<initial_set<<"\n    final_time="<<final_time<<"\n    bounding_box="<<bounding_box<<"\n");
  time_model_type zero_time_model(I(0),LinearAlgebra::Vector<I>(initial_set.number_of_generators()));
  time_model_type final_time_model=final_time+zero_time_model;
  time_model_type time_step=final_time-initial_set.time();
  timed_set_type final_set=this->continuous_evolution_step(mode,initial_set,time_step,bounding_box);
  return timed_set_type(final_time_model,final_set.steps(),final_set.discrete_state(),final_set.continuous_state_set());
}



template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::continuous_evolution_step(const mode_type& mode,
                                                                        const timed_set_type& initial_set,
                                                                        const time_model_type& time_step,
                                                                        const bounding_box_type& bounding_box) const
{
  ARIADNE_LOG(7," ConstraintBasedHybridEvolverPlugin::continuous_evolution_step(...)\n");
  ARIADNE_LOG(8,"    dynamic="<<mode.dynamic()<<"\n    initial_set="<<initial_set<<"\n    time_step="<<time_step<<"\n    bounding_box="<<bounding_box<<"\n");
  assert(mode.id()==initial_set.discrete_state());
  I average_time_step = time_step.average();
  continuous_basic_set_type continuous_state_set = 
      this->_integrator->integration_step(mode.dynamic(),initial_set.continuous_state_set(),average_time_step,bounding_box);
  Matrix<I>& generators = const_cast<Matrix<I>&>(continuous_state_set.generators());
  generators += outer_product(mode.dynamic()(bounding_box),time_step.gradient());
  timed_set_type integrated_set(initial_set.time()+time_step,initial_set.steps(),initial_set.discrete_state(),continuous_state_set);
  ARIADNE_LOG(8,"    integrated_set="<<integrated_set<<"\n\n");
  return integrated_set;
}



template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::continuous_reachability_step(const mode_type& mode,
                                                                           const timed_set_type& initial_set,
                                                                           const time_model_type& lower_time_step,
                                                                           const time_model_type& upper_time_step,
                                                                           const bounding_box_type& bounding_box) const
{
  ARIADNE_LOG(7," ConstraintBasedHybridEvolverPlugin::continuous_reachability_step(...)\n");
  ARIADNE_LOG(8,"    dynamic="<<mode.dynamic()<<"\n    initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(8,"    lower_time_step="<<lower_time_step<<"\n    upper_time_step="<<upper_time_step<<"\n    bounding_box="<<bounding_box<<"\n");
  assert(mode.id()==initial_set.discrete_state());

  time_model_type average_time_step=(lower_time_step+upper_time_step)/2;
  time_type integration_time_step=average_time_step.average().midpoint();
  continuous_basic_set_type continuous_state_set = 
      this->_integrator->integration_step(mode.dynamic(),initial_set.continuous_state_set(),integration_time_step,bounding_box);

  R lower_time_bound=lower_time_step.bound().lower();
  R upper_time_bound=upper_time_step.bound().upper();
  ARIADNE_LOG(8,"    lower_time_bound="<<lower_time_bound<<", upper_time_bound="<<upper_time_bound<<"\n");
  I time_interval=Interval<R>(lower_time_bound,upper_time_bound);
  R time_interval_radius=time_interval.radius();

  Vector<I> new_generator=time_interval_radius*(mode.dynamic().image(bounding_box));
  Geometry::Zonotope<I> reach_continuous_state_set(continuous_state_set.centre(),continuous_state_set.generators(),new_generator);

  time_model_type reach_time(initial_set.time()+average_time_step);
  reach_time=time_model_type(reach_time.average(),reach_time.gradient(),time_interval_radius);
  
  timed_set_type reached_set(reach_time,initial_set.steps(),initial_set.discrete_state(),reach_continuous_state_set);
  ARIADNE_LOG(8,"    reached_set="<<reached_set<<"\n\n");
  return reached_set;
}


template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::discrete_event_step(const transition_type& transition,
                                                                  const timed_set_type& initial_set) const
{
  ARIADNE_LOG(7," ConstraintBasedHybridEvolverPlugin::discrete_event_step(...)\n");
  assert(transition.source().id()==initial_set.discrete_state());
  continuous_basic_set_type continuous_state_set = 
    this->_applicator->evaluate(transition.reset(),initial_set.continuous_state_set());
  return timed_set_type(initial_set.time(),initial_set.steps()+1,transition.destination().id(),continuous_state_set);
}



/*! \brief Compute the possible states reached by an fored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::forced_jump_step(const transition_type& transition,
                                                               const timed_set_type& initial_set,
                                                               const time_model_type& crossing_time,
                                                               const bounding_box_type& bounding_box) const
{
  ARIADNE_LOG(7,"ConstraintBasedHybridEvolverPlugin::forced_jump_step(...)\n");
  ARIADNE_LOG(8,"  transition="<<transition<<"\n");
  ARIADNE_LOG(8,"  initial_set="<<initial_set.continuous_state_set()<<"\n");
  ARIADNE_LOG(8,"  initial_time="<<initial_set.time()<<"\n");
  ARIADNE_LOG(8,"  crossing_time="<<crossing_time<<"\n");
  ARIADNE_LOG(8,"  bounding_box="<<bounding_box<<"\n");
  
  timed_set_type flowed_set = this->continuous_evolution_step(transition.source(),initial_set,crossing_time,bounding_box);
  timed_set_type mapped_set = this->discrete_event_step(transition,flowed_set);
  return mapped_set;
}



/*! \brief Compute the possible states reached by an unfored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::unforced_jump_step(const transition_type& transition,
                                                                 const timed_set_type& initial_set,
                                                                 const time_model_type& minimum_time,
                                                                 const time_model_type& maximum_time,
                                                                 const bounding_box_type& bounding_box) const
{
  ARIADNE_LOG(7," ConstraintBasedHybridEvolverPlugin::unforced_jump_step(...)\n");
  ARIADNE_LOG(8,"  transition="<<transition<<"\n");
  ARIADNE_LOG(8,"  initial_set="<<initial_set.continuous_state_set()<<"\n");
  ARIADNE_LOG(8,"  initial_time="<<initial_set.time()<<"\n");
  ARIADNE_LOG(8,"  minimum_time="<<minimum_time<<"\n");
  ARIADNE_LOG(8,"  maximum_time="<<maximum_time<<"\n");
  ARIADNE_LOG(8,"  bounding_box="<<bounding_box<<"\n");
  
  time_model_type zero_time(maximum_time.gradient().size());
  timed_set_type flowed_set = this->continuous_reachability_step(transition.source(),initial_set,zero_time,maximum_time,bounding_box);
  timed_set_type mapped_set = this->discrete_event_step(transition,flowed_set);
  return mapped_set;
}







template<class R>
std::vector<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type>
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::evolution_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                  const timed_set_type& initial_set,
                                                                  const time_type& final_time,
                                                                  EvolutionSemantics evolution_semantics,
                                                                  EvolutionKind evolution_kind) const
{
  // FIXME: this is a hack!
  R maximum_splitting_set_radius = 0.0625;

  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

  typedef typename reference_set<const transition_type>::const_iterator transitions_const_iterator;
  typedef typename std::map<id_type,time_model_type>::const_iterator crossing_times_const_iterator;

  this->trace.clear();

  ARIADNE_LOG(2,"\nConstraintBasedHybridEvolverPlugin::evolution_step(...)\n");
  ARIADNE_LOG(3,"  evolution_semantics="<<(evolution_semantics==lower_semantics?"lower":"upper"));
  ARIADNE_LOG(3,", evolution_kind="<<(evolution_kind==compute_evolved_set?"evolve":"reach")<<"\n");  
  ARIADNE_LOG(3,"  initial_steps="<<initial_set.steps()<<"  initial_mode="<<initial_set.discrete_state()<<"\n")
  ARIADNE_LOG(3,"  initial_set="<<Geometry::over_approximation(initial_set.continuous_state_set())<<"\n");  
  ARIADNE_LOG(3,"  initial_time="<<initial_set.time()<<"\n\n");  

  std::vector<timed_set_type> result;

  const id_type& discrete_state = initial_set.discrete_state();
  const continuous_basic_set_type& basic_set = initial_set.continuous_state_set();
  //const time_model_type& initial_time = initial_set.time();
  //const discrete_time_type& initial_steps = initial_set.steps();

  const mode_type& mode = automaton.mode(discrete_state);
  const System::VectorFieldInterface<R>& dynamic=mode.dynamic();
  const reference_set<const transition_type>& transitions=mode.transitions();

  ARIADNE_LOG(6,"  transitions="<<transitions<<"\n\n");
  
  // Easy references to the spacial time and the time relative to the set generators
  time_model_type zero_spacial_time_step(I(0),LinearAlgebra::Vector<I>(initial_set.dimension()));
  time_model_type zero_time_step=initial_set.time()*0;

  // Compute rough bounding box
  time_type time_step_size=(final_time-initial_set.time()).bound().upper();
  Rectangle<R> bounding_box=this->_bounder->estimate_flow_bounds(dynamic,basic_set.bounding_box(),time_step_size);
  ARIADNE_LOG(6,"  maximum_time_step_size="<<time_step_size<<"\n");
  ARIADNE_LOG(6,"  rough_bounding_box="<<bounding_box<<"\n");
  bounding_box=this->_bounder->refine_flow_bounds(dynamic,basic_set.bounding_box(),bounding_box,time_step_size);
  ARIADNE_LOG(6,"  improved_bounding_box="<<bounding_box<<"\n");
  bounding_box=this->_bounder->refine_flow_bounds(dynamic,basic_set.bounding_box(),bounding_box,time_step_size);
  ARIADNE_LOG(6,"  further_improved_bounding_box="<<bounding_box<<"\n\n");

  // Log some facts about the flow
  ARIADNE_LOG(6,"  dynamic(initial_set.centre())="<<dynamic(midpoint(initial_set.centre()))<<"\n");
  ARIADNE_LOG(6,"  dynamic(bounding_box)="<<dynamic(Geometry::Point<I>(bounding_box))<<"\n");
  ARIADNE_LOG(6,"  dynamic.jacobian(bounding_box)="<<dynamic.jacobian(Geometry::Point<I>(bounding_box))<<"\n\n");

  // Set the preferred step for integration
  time_model_type integration_time_step = zero_time_step + time_step_size;
  time_type integration_time_step_size = time_step_size;
  // Compute the evolved set after the maximum evolution time
  timed_set_type integrated_set = this->continuous_evolution_step(mode,initial_set,integration_time_step,bounding_box);
  ARIADNE_LOG(6,"  integrated_set="<<integrated_set<<"\n\n");

  // Schedule events bounding the evolution time
  std::map<id_type,time_model_type> terminating_event_times = this->compute_terminating_event_times(mode,initial_set,integrated_set,final_time,bounding_box,time_step_size);
  assert(terminating_event_times.size()>=1);

  // Perform switching logic if there is more than one terminating event, and compute the terminating event time step for upper evolution
  time_model_type terminating_time_step;
  if(terminating_event_times.size()==1) {
    terminating_time_step=terminating_event_times.begin()->second;
  } else {
    // Compute the absolute maximim and minimum time steps
    I terminating_time_bound=this->compute_evolution_time_bounds(terminating_event_times);
    if(terminating_time_bound>integration_time_step_size) {
      // events are never triggered
      terminating_time_step=integration_time_step;
      terminating_event_times.clear();
    } else if(terminating_time_bound<integration_time_step_size) {
      if (terminating_event_times.size()==1) {
        terminating_time_step=terminating_event_times.begin()->second;
      } else {
        terminating_time_step=terminating_time_bound.upper()+zero_time_step;
      }
    } else {
      // See if we can reduce the step size to prevent hanging on the final event 
      if(terminating_time_bound.lower() > time_step_size) { 
        // Make the (unique) terminating event the flow event
        integration_time_step=time_type(terminating_time_bound.lower())+zero_time_step;
        terminating_time_step=integration_time_step;
        terminating_event_times.clear();
      } else {
        // Make the terminating time the integration time for safety
        terminating_time_step=integration_time_step;
      }
    }
  }
  ARIADNE_LOG(6,"  terminating_time_step="<<terminating_time_step<<"\n\n");
  
  if (terminating_event_times.size()>1) {
    // Subdivide if upper semantics and set is not too small
    if(evolution_semantics==upper_semantics && initial_set.radius()>maximum_splitting_set_radius) {
      ARIADNE_LOG(3," Multiply enabled events: subdividing\n");
      return this->subdivide(initial_set);
    }
  }


  // Schedule events between current time and integration time step
  std::map<id_type,time_model_pair_type> enabled_activations=this->compute_enabled_activation_times(mode,initial_set,integrated_set,terminating_time_step,bounding_box);
  
  


  // Add the sets reached by unforced steps (same for integration and reachability)
  for(typename reference_set<transition_type>::const_iterator transition_iter = transitions.begin();
      transition_iter!=transitions.end(); ++transition_iter)
  {
    if(enabled_activations.count(transition_iter->id())) {
      const transition_type& transition = *transition_iter;
      const id_type& event_id = transition_iter->id();
      const time_model_pair_type& switching_times = enabled_activations.find(event_id)->second;
      const time_model_type& lower_switching_time = switching_times.first;
      const time_model_type& upper_switching_time = switching_times.second;
      
      timed_set_type jump_set=this->unforced_jump_step(transition,initial_set,lower_switching_time,upper_switching_time,bounding_box);
      result.push_back(jump_set);
      ARIADNE_LOG(6,"  unforced_jump_set["<<event_id<<"]="<<jump_set<<"\n");
      this->trace.push_back(jump_set);
    }
  }
  ARIADNE_LOG(6,"\n");


  // Add the final sets after the evolution for both integration and reachability 
  if(lower_semantics && terminating_event_times.size()>1) {
    ARIADNE_LOG(6," Multiple terminating events for lower evolution; cannot decide which to take.\n");
    // If lower semantics and terminating_event_times.size()>1, then standard evolution blocks
  } else {
    // Add all possible sets reached by continuous evolution possibly followed by a jump step
    for(crossing_times_const_iterator crossing_time_iter=terminating_event_times.begin();
        crossing_time_iter!=terminating_event_times.end(); ++crossing_time_iter)
    {
      id_type terminating_event=crossing_time_iter->first;
      const time_model_type& terminating_time_step=crossing_time_iter->second;

      if(terminating_event==time_step_event) {
        // The terminating event is time passing
        timed_set_type continuous_evolved_set=this->continuous_evolution_step(mode,initial_set,terminating_time_step,bounding_box);
        result.push_back(continuous_evolved_set);
        this->trace.push_back(continuous_evolved_set);
        ARIADNE_LOG(6,"  continuous_evolved_set="<<continuous_evolved_set<<"\n");
      } else if(terminating_event==final_time_event) {
        // The terminating event is reaching the requested integration time
        timed_set_type continuous_evolved_set=this->final_continuous_evolution_step(mode,initial_set,final_time,bounding_box);
        result.push_back(continuous_evolved_set);
        this->trace.push_back(continuous_evolved_set);
        ARIADNE_LOG(6,"  continuous_evolved_set="<<continuous_evolved_set<<"\n");
      } else if(terminating_event==regularizing_time_event) {
        // The terminating event is reaching the requested integration time
        time_type regularizing_time=(initial_set.time()+terminating_time_step).average().midpoint();
        timed_set_type continuous_evolved_set=this->final_continuous_evolution_step(mode,initial_set,regularizing_time,bounding_box);
        result.push_back(continuous_evolved_set);
        this->trace.push_back(continuous_evolved_set);
        ARIADNE_LOG(6,"  continuous_evolved_set="<<continuous_evolved_set<<"\n");
      } else {
        // The terminating event is a transition, which may be a blocking (invariant) transition, or a forced jump
        const transition_type& blocking_transition = automaton.transition(terminating_event,discrete_state);
        if(blocking_transition.kind()==guard_tag) {
          // Note that we do not need to check if the terminating time step here exceeds the maximum time step, since either:
          //  - we are computing over approximations, and the badly behaved part of the flow is covered by the integration step, or
          //  - we are computing lower approximations, and the terminating time step used cannot we cannot exceed the maximum time step
          timed_set_type jump_set=forced_jump_step(blocking_transition,initial_set,terminating_time_step,bounding_box);
          result.push_back(jump_set);
          this->trace.push_back(jump_set);
          ARIADNE_LOG(6,"  forced_jump_set["<<blocking_transition.id()<<"]="<<jump_set<<"\n");
        } else if(blocking_transition.kind()==invariant_tag) {
          ARIADNE_LOG(6,"  invariant["<<blocking_transition.id()<<"] reached\n");
        } else {
          throw std::runtime_error("Terminating event should block evolution, but unrecognised event type");
        }

      }
    }
  }

  // Compute the continuous reachable set
  timed_set_type continuous_reachable_set = this->continuous_reachability_step(mode,initial_set,zero_time_step,terminating_time_step,bounding_box);
  ARIADNE_LOG(6,"  continuous_reachable_set="<<continuous_reachable_set<<"\n");


  if(evolution_kind==compute_reachable_set) {
    result.push_back(continuous_reachable_set);
    //this->trace.push_back(continuous_reachable_set);
  } else {
  }

  ARIADNE_LOG(6,"\n  evolution_step::result="<<result<<"\n\n");
  return result;
}













template<class R>
std::vector<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type> 
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::lower_evolution_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                   const timed_set_type& initial_set,
                                                                   const time_type& maximum_time) const
{
  return this->evolution_step(automaton, initial_set, maximum_time, lower_semantics, compute_evolved_set);
}


template<class R>
std::vector<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type> 
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::upper_evolution_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                   const timed_set_type& initial_set,
                                                                   const time_type& maximum_time) const
{
  ARIADNE_LOG(2,"\nConstraintBasedHybridEvolverPlugin::upper_evolution_step(...)\n");
  return this->evolution_step(automaton, initial_set, maximum_time, upper_semantics, compute_evolved_set);
}


template<class R>
std::vector<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type> 
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::lower_reachability_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const timed_set_type& initial_set,
                                                                      const time_type& maximum_time) const
{
  return this->evolution_step(automaton, initial_set, maximum_time, lower_semantics, compute_reachable_set);
}


template<class R>
std::vector<typename Evaluation::ConstraintBasedHybridEvolverPlugin<R>::timed_set_type> 
Evaluation::ConstraintBasedHybridEvolverPlugin<R>::upper_reachability_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const timed_set_type& initial_set,
                                                                      const time_type& maximum_time) const
{
  ARIADNE_LOG(2,"\nConstraintBasedHybridEvolverPlugin::upper_reachability_step(...)\n");
  return this->evolution_step(automaton, initial_set, maximum_time, upper_semantics, compute_reachable_set);
}

} // namespace Ariadne
