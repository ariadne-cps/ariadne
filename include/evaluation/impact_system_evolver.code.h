/***************************************************************************
 *            impact_system_evolver.code.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
#include "impact_system_evolver.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "base/stack.h"
#include "base/tuple.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/exceptions.h"
#include "function/function_interface.h"
#include "function/affine_function.h"
#include "function/constant_function.h"
#include "function/approximate_taylor_model.h"

#include "geometry/box.h"
#include "geometry/list_set.h"
#include "geometry/timed_list_set.h"
#include "geometry/taylor_set.h"

#include "system/impact_system.h"

#include "evaluation/evolution_profiler.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/standard_satisfier.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/standard_reducer.h"
#include "evaluation/cascade_reducer.h"

#include "output/logging.h"

namespace Ariadne {
  
class DegenerateCrossingException { };

const uint spacial_order=2;
const uint temporal_order=4;
const uint order=spacial_order+temporal_order;
const uint smoothness=1;



template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver()
{
  EvolutionParameters<R> parameters;
  StandardApplicator< ES > applicator;
  StandardIntegrator< ES > integrator;
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  *this = Evolver(parameters,applicator,integrator,satisfier,subdivider,reducer);
}


template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver(const EvolutionParameters<R>& parameters)
{
  StandardApplicator< ES > applicator;
  StandardIntegrator< ES > integrator;
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  *this = Evolver(parameters,applicator,integrator,satisfier,subdivider,reducer);
}

template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver(const EvolutionParameters<R>& parameters,
        const ApplicatorInterface<ES>& applicator, 
        const IntegratorInterface<ES>& integrator)
{
  StandardSatisfier< ES > satisfier;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  *this = Evolver(parameters,applicator,integrator,satisfier,subdivider,reducer);
}


template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver(const EvolutionParameters<R>& parameters,
        const ApplicatorInterface<ES>& applicator, 
        const IntegratorInterface<ES>& integrator, 
        const SatisfierInterface<ES>& satisfier, 
        const SubdividerInterface<ES>& subdivider, 
        const ReducerInterface<ES>& reducer)
  : _parameters(parameters.clone()),
    _integrator(integrator.clone()),
    _subdivider(subdivider.clone()),
    _reducer(reducer.clone()),
    _profiler(new EvolutionProfiler),
    verbosity(0)
{ }





template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
model(const ES& set, ushort d) const
{
  ushort order=d;
  ushort smoothness=d;
  Vector<I> domain=set.domain().position_vectors();
  Vector<R> midpoint=Ariadne::midpoint(domain);
  Vector<A> centre=Vector<A>(set.centre().position_vector());
  Matrix<A> const& generators=reinterpret_cast<Matrix<A>const&>(set.generators());
  return ApproximateTaylorModel<R>::affine(domain,midpoint,centre,generators,order,smoothness);
}

  
template<class ES>
ES
Evolver<ImpactSystem<typename ES::real_type>,ES>::
set(const ApproximateTaylorModel<R>& model) const
{
  uint n=model.result_size();
  uint m=model.argument_size();
  Vector<I> set_centre(n);
  Matrix<R> set_generators(n,m);
  make_lpair(set_centre,set_generators) = affine_model(model);
  ES enclosure_set(Point<I>(set_centre),set_generators);
  return enclosure_set;
}



  
template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_integration_step(const ATM& flow_model, const ATM& initial_set_model, const A& integration_time) const
{
  uint n=flow_model.result_size();
  ApproximateTaylorModel<R> integration_time_model
    =ApproximateTaylorModel<R>::constant(Vector<I>(n,I(-1,1)),Vector<R>(n,R(0)),Vector<A>(1,integration_time),order,smoothness);
  ARIADNE_LOG(6,"integration_time_step_model = "<<integration_time_model<<"\n");

  return this->_integration_step(flow_model,initial_set_model,integration_time_model);
}

template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_integration_step(const ATM& flow_model, const ATM& initial_set_model, const ATM& integration_time_model) const
{
    ApproximateTaylorModel<R> set_step_model=join(initial_set_model, integration_time_model);
    ARIADNE_LOG(6,"set_step_model = "<<set_step_model<<"\n");
    ApproximateTaylorModel<R> final_set_model=compose(flow_model,set_step_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");

    return final_set_model;
}



template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_reachability_step(const ATM& flow_model, const ATM& initial_set_model, const A& initial_time, const A& final_time) const
{
  A time_midpoint=(A(initial_time)+A(final_time))/2;
  A time_radius=(A(final_time)-A(initial_time))/2;
  ApproximateTaylorModel<R> time_interval_model=ApproximateTaylorModel<R>::affine(I(-1,1),R(0),time_midpoint,time_radius,order,smoothness);
  ARIADNE_LOG(6,"time_interval_time_model="<<time_interval_model<<"\n");
  ApproximateTaylorModel<R> expanded_timed_set_model=combine(initial_set_model,time_interval_model);
  ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
  ApproximateTaylorModel<R> reach_set_model=compose(flow_model,expanded_timed_set_model);
  ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
  
  return reach_set_model;
}


template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_reachability_step(const ATM& flow_model, const ATM& initial_set_model, const A& initial_time, const ATM& final_time_model) const
{
  ApproximateTaylorModel<R> initial_time_model
    =ApproximateTaylorModel<R>::constant(final_time_model.domain(),final_time_model.centre(),Vector<A>(1u,initial_time),
                                         final_time_model.order(),final_time_model.smoothness());
  
  return this->_reachability_step(flow_model, initial_set_model, initial_time_model, final_time_model);
}


// Compute the grazing time using bisections
template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_hitting_time_model(const ATM& flow_model, const ATM& guard_model, const ATM& initial_set_model, const A& initial_time, const A& final_time) const
{
  ApproximateTaylorModel<R> hitting_model=compose(guard_model,flow_model);
  ARIADNE_LOG(6,"hitting_model = "<<hitting_model<<"\n");
  ApproximateTaylorModel<R> free_hitting_time_model;
  try {
    free_hitting_time_model=implicit(hitting_model); 
  } catch(NonInvertibleFunctionException) {
    throw DegenerateCrossingException();
  }
  ARIADNE_LOG(6,"free_hitting_time_model = "<<free_hitting_time_model<<"\n");
  ApproximateTaylorModel<R> hitting_time_model=compose(free_hitting_time_model,initial_set_model);
  ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
  Interval<R> hitting_time_range=hitting_time_model.range()[0];
  ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
  if(hitting_time_range.lower()<R(initial_time) || hitting_time_range.upper()>R(final_time)) {
    throw DegenerateCrossingException();
  }
 
  return hitting_time_model;
}


// Compute the grazing time using bisections
template<class ES>
Interval<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_crossing_time_interval(const ATM& flow_model, const ATM& guard_model, const ATM& initial_set_model, 
                        const A& initial_time, const A& final_time) const
{
  ApproximateTaylorModel<R> final_set_model=this->_integration_step(flow_model,initial_set_model,final_time);
  assert(definitely(this->_active(guard_model,final_set_model)));

  uint refinements=5;

  A lower_time=final_time;
  A upper_time=initial_time;
  
  if(possibly(this->_active(guard_model,initial_set_model))) {
    lower_time=initial_time;
  } else {
    A min_lower_time=initial_time;
    A max_lower_time=final_time;
    for(uint i=0; i!=refinements; ++i) {
      A new_lower_time=(min_lower_time+max_lower_time)/2;
      if(possibly(this->_active(guard_model,this->_reachability_step(flow_model,initial_set_model,initial_time,new_lower_time)))) {
        max_lower_time=new_lower_time;
      } else {
        min_lower_time=new_lower_time;
        lower_time=new_lower_time;
      }
    }
  }

  if(definitely(this->_active(guard_model,initial_set_model))) {
    upper_time=initial_time;
  } else {
    A min_upper_time=initial_time;
    A max_upper_time=final_time;
    for(uint i=0; i!=refinements; ++i) {
      A new_upper_time=(min_upper_time+max_upper_time)/2;
      if(possibly(this->_active(guard_model,this->_integration_step(flow_model,initial_set_model,new_upper_time)))) {
        min_upper_time=new_upper_time;
      } else {
        max_upper_time=new_upper_time;
        upper_time=new_upper_time;
      }
    }
  }

  return Interval<R>(R(lower_time),R(upper_time));
}



// Compute the grazing time using bisections
template<class ES>
Interval<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_grazing_time_interval(const ATM& flow_model, const ATM& guard_model, const ATM& initial_set_model, const A& initial_time, const A& final_time) const
{
  uint refinements=5;

  ATM final_set_model=this->_integration_step(flow_model,initial_set_model,final_time);
  ATM reach_set_model=this->_reachability_step(flow_model,initial_set_model,initial_time,final_time);
  
  A lower_time=final_time;
  A upper_time=initial_time;
  
  if(definitely(!this->_active(guard_model,reach_set_model))) {
    return Interval<R>(R(lower_time),R(upper_time));
  }

  if(possibly(this->_active(guard_model,initial_set_model))) {
    lower_time=initial_time;
  } else {
    A min_lower_time=initial_time;
    A max_lower_time=final_time;
    for(uint i=0; i!=refinements; ++i) {
      A new_lower_time=(min_lower_time+max_lower_time)/2;
      if(possibly(this->_active(guard_model,this->_reachability_step(flow_model,initial_set_model,initial_time,new_lower_time)))) {
        max_lower_time=new_lower_time;
      } else {
        min_lower_time=new_lower_time;
        lower_time=new_lower_time;
      }
    }
  }

  if(possibly(this->_active(guard_model,final_set_model))) {
    upper_time=final_time;
  } else {
    A min_upper_time=initial_time;
    A max_upper_time=final_time;
    for(uint i=0; i!=refinements; ++i) {
      A new_upper_time=(min_upper_time+max_upper_time)/2;
      if(possibly(this->_active(guard_model,this->_reachability_step(flow_model,initial_set_model,new_upper_time,final_time)))) {
        min_upper_time=new_upper_time;
      } else {
        max_upper_time=new_upper_time;
        upper_time=new_upper_time;
      }
    }
  }

  return Interval<R>(R(lower_time),R(upper_time));
}
  
  


    
template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_reachability_step(const ATM& flow_model, const ATM& initial_set_model, const ATM& initial_time_model, const ATM& final_time_model) const
{

    // Compute the reachable set
    // Need an extra independent variable to represent time
    uint n=flow_model.result_size();
    
    ApproximateTaylorModel<R> expanded_initial_set_model=embed(initial_set_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
    ApproximateTaylorModel<R> expanded_initial_time_model=embed(initial_time_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_initial_time_model="<<expanded_initial_time_model<<"\n");
    ApproximateTaylorModel<R> expanded_final_time_model=embed(final_time_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_final_time_model="<<expanded_final_time_model<<"\n");
    
    ApproximateTaylorModel<R> time_interval_model=ApproximateTaylorModel<R>::affine(I(-1,1),R(0),A(0.5),A(0.5),order,smoothness);
    ARIADNE_LOG(6,"time_interval_time_model="<<time_interval_model<<"\n");
    ApproximateTaylorModel<R> expanded_time_interval_model=embed(time_interval_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),n);
    ARIADNE_LOG(6,"expanded_time_interval_model="<<expanded_time_interval_model<<"\n");
    ApproximateTaylorModel<R> expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);
    ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
    ApproximateTaylorModel<R> expanded_timed_set_model=join(expanded_initial_set_model,expanded_reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    ApproximateTaylorModel<R> reach_set_model=compose(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    
    return reach_set_model;
}

template<class ES>
tribool
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_active(const ATM& guard_model, const ATM& set_model) const
{
  I range=compose(guard_model,set_model).range()[0];
  if(range.upper()<0) { 
    return false; 
  } else if(range.lower()>0) {
    return true;
  } else {
    return indeterminate;
  }
}

template<class ES>
tribool
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_active(const FN& guard_function, const ATM& set_model) const
{
  ApproximateTaylorModel<R> guard_model(set_model.range(),guard_function,spacial_order,smoothness);
  return this->_active(guard_model,set_model);
}

template<class ES>
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_evolution(ESL& final_sets, 
           ESL& reach_sets, 
           ESL& intermediate_sets, 
           const Sys& system, 
           const ES& initial_set, 
           const T& maximum_hybrid_time, 
           Semantics semantics, 
           bool reach) const
{
  enum FlowType { TRANSVERSE, CROSSING, GRAZING, MISSING, UNKNOWN };

  ARIADNE_LOG(5,__PRETTY_FUNCTION__);
  assert(semantics==upper_semantics);

  const Rational maximum_time=maximum_hybrid_time.time();
  const Integer maximum_steps=maximum_hybrid_time.steps();

  const A max_time=A(maximum_time);
  const A max_step_size=A(this->maximum_step_size());

  const uint n=system.state_space().dimension();
  const uint spacial_order=2;
  const uint temporal_order=4;
  const uint order=spacial_order+temporal_order;
  const uint smoothness=1;


  const FunctionInterface<R>& reset_map=system.impact_map();
  const FunctionInterface<R>& vector_field=system.vector_field();
  const FunctionInterface<R>& guard_function=system.guard_condition();

  typedef tuple<Integer, ApproximateTaylorModel<R> > hybrid_timed_set_type;

  stack< hybrid_timed_set_type > working_sets;

  {
    // Set up initial timed set models
    ApproximateTaylorModel<R> initial_set_model=this->model(initial_set,spacial_order);
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    ApproximateTaylorModel<R> initial_time_model
      =ApproximateTaylorModel<R>::constant(initial_set_model.domain(),initial_set_model.centre(),
                                           Vector<A>(1,A(0)),order,smoothness);
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ApproximateTaylorModel<R> initial_timed_set_model=join(initial_set_model,initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push(make_tuple(Integer(0),initial_timed_set_model));
  }

  while(!working_sets.empty()) {
    Integer initial_steps;
    ApproximateTaylorModel<R> initial_timed_set_model;
    make_ltuple(initial_steps,initial_timed_set_model)=working_sets.pop();
    ARIADNE_LOG(4,"initial_steps = "<<initial_steps<<"\n");
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");

    const ApproximateTaylorModel<R> initial_set_model=project(initial_timed_set_model,range(0,n));
    ARIADNE_LOG(4,"initial_set_model = "<<initial_set_model<<"\n");
    const ApproximateTaylorModel<R> initial_time_model=project(initial_timed_set_model,range(n,n+1));
    ARIADNE_LOG(4,"initial_time_model = "<<initial_time_model<<"\n");
    I initial_time_range=initial_time_model.range()[0];
    ARIADNE_LOG(4,"initial_time_range = "<<initial_time_range<<"\n");
    ARIADNE_ASSERT(A(initial_time_range.width()) <= max_step_size);

    ++this->_profiler->time_steps;
    IVec flow_bounds; 
    Rational h;
    I guard_value(-1,1);

    // Try to unwind time interval
    if(initial_time_range.width()>0) {
      h=Rational(initial_time_range.width());
      make_lpair(h,flow_bounds)=this->flow_bounds(vector_field,initial_set_model.range(),h);
      guard_value=guard_function.evaluate(flow_bounds)[0];
    }
    if(guard_value.upper()>0) {
      h=this->maximum_step_size();
      make_lpair(h,flow_bounds)=this->flow_bounds(vector_field,initial_set_model.range(),h);
      guard_value=guard_function.evaluate(flow_bounds)[0];
    }

    this->_profiler->total_stepping_time+=h;
    this->_profiler->minimum_time_step=std::min(h,this->_profiler->minimum_time_step);
    A step_size=A(h);
    ARIADNE_LOG(6,"set_bounding_box = "<<initial_set_model.range()<<"\n");
    ARIADNE_LOG(6,"flow_bounds = "<<flow_bounds<<"\n");
    ARIADNE_LOG(6,"step_size = "<<step_size<<"\n");
    ARIADNE_LOG(6,"guard_value = "<<guard_value<<"\n\n");
  
    FlowType flow_type=UNKNOWN;
  
    if(guard_value <= 0) { flow_type = MISSING; }

    
    ApproximateTaylorModel<R> vector_field_model(flow_bounds,vector_field,order,smoothness);
    ARIADNE_LOG(6,"vector_field_model = "<<vector_field_model<<"\n");
    ApproximateTaylorModel<R> flow_model=Ariadne::flow(vector_field_model);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
    
    A zero_time=0;
    A step_time=step_size;
    A final_time;

    Interval<R> crossing_time_interval;
    Interval<R> grazing_time_interval;

    ApproximateTaylorModel<R> guard_model;

    ApproximateTaylorModel<R> hitting_time_model;
    ApproximateTaylorModel<R> integration_time_model;

    ApproximateTaylorModel<R> active_set_model;
    ApproximateTaylorModel<R> jump_set_model;
    ApproximateTaylorModel<R> jump_time_model;
    ApproximateTaylorModel<R> final_set_model;
    ApproximateTaylorModel<R> final_time_model;
    ApproximateTaylorModel<R> reach_set_model;
    ApproximateTaylorModel<R> reach_time_model;
    
    if(flow_type==UNKNOWN) {
      A reduced_step_time=step_time;
      for(uint i=0; i!=4; ++i) {
        reach_set_model=this->_reachability_step(flow_model,initial_set_model,zero_time,reduced_step_time);
        if(definitely(!this->_active(guard_function,reach_set_model))) {
          flow_type=MISSING;
          step_time=reduced_step_time;
          break;
        } else {
          reduced_step_time/=2;
        }
      }
    }

    if(flow_type==MISSING) {
      final_time=A(initial_time_range.lower())+step_time;
      final_time_model=ApproximateTaylorModel<R>::constant(Vector<I>(n,I(-1,1)),Vector<R>(n,R(0)),Vector<A>(1u,final_time),spacial_order,smoothness);
      integration_time_model=final_time_model-initial_time_model;
      final_set_model=this->_integration_step(flow_model,initial_set_model,integration_time_model);
      reach_set_model=this->_reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);
    } else {
      guard_model=ApproximateTaylorModel<R>(flow_bounds,guard_function,order,smoothness);
      ARIADNE_LOG(6,"guard_model = "<<guard_model<<"\n");   
      try {
        hitting_time_model=this->_hitting_time_model(flow_model,guard_model,initial_set_model,zero_time,step_time);
        jump_time_model=initial_time_model+hitting_time_model;
        active_set_model=this->_integration_step(flow_model,initial_set_model,hitting_time_model);
        reach_set_model=this->_reachability_step(flow_model,initial_set_model,zero_time,hitting_time_model);
        flow_type=TRANSVERSE;
      } 
      catch(DegenerateCrossingException) { 
        final_set_model = this->_integration_step(flow_model,initial_set_model,step_size);
        integration_time_model=ApproximateTaylorModel<R>::constant(Vector<I>(n,I(-1,1)),Vector<R>(n,R(0)),Vector<A>(1u,step_size),spacial_order,smoothness);
        final_time_model = initial_time_model + integration_time_model;
        if(definitely(this->_active(guard_function,final_set_model))) {
          flow_type=CROSSING;
          crossing_time_interval=this->_crossing_time_interval(flow_model, guard_model,initial_set_model, zero_time, step_time);
          A crossing_time_midpoint=A(crossing_time_interval.midpoint());
          //active_set_model=this->_reachable_set(flow_model, guard_model,initial_set_model,A(crossing_time_interval.lower()),A(crossing_time_interval.upper()));
          active_set_model=this->_integration_step(flow_model,initial_set_model,A(crossing_time_interval.midpoint()));
          jump_time_model=ApproximateTaylorModel<R>::constant(Vector<I>(n,I(-1,1)),Vector<R>(n,R(0)),Vector<A>(1u,crossing_time_midpoint),spacial_order,smoothness);
          reach_set_model=this->_reachability_step(flow_model,initial_set_model,zero_time,A(crossing_time_interval.upper()));
        } else {
          flow_type=GRAZING;
          grazing_time_interval=this->_grazing_time_interval(flow_model, guard_model,initial_set_model, zero_time, step_time);
          A grazing_time_midpoint=A(grazing_time_interval.midpoint());
          //active_set_model=this->_reachability_step(flow_model, guard_model,initial_set_model,A(grazing_time_interval.lower()),A(grazing_time_interval.upper()));
          active_set_model=this->_integration_step(flow_model,initial_set_model,A(grazing_time_midpoint));
          jump_time_model=ApproximateTaylorModel<R>::constant(Vector<I>(n,I(-1,1)),Vector<R>(n,R(0)),Vector<A>(1u,grazing_time_midpoint),spacial_order,smoothness);
          reach_set_model=this->_reachability_step(flow_model,initial_set_model,A(crossing_time_interval.lower()),A(grazing_time_interval.upper()));
        } 
      }
    }
        
    
    if(flow_type == TRANSVERSE) { ARIADNE_LOG(4,"\nTRANSVERSE\n\n"); }
    if(flow_type == CROSSING ) { ARIADNE_LOG(4,"\nCROSSING\n\n"); }
    if(flow_type == GRAZING ) { ARIADNE_LOG(4,"\nGRAZING \n\n"); }
    if(flow_type == MISSING) { ARIADNE_LOG(4,"\nMISSING\n\n"); }
    

    // Compute the jump set
    if(flow_type!=MISSING) {
      Vector<A> reset_model_centre(active_set_model.evaluate(Vector<A>(active_set_model.centre())));
      ARIADNE_LOG(6,"active_set_model = "<<active_set_model<<"\n");
      ApproximateTaylorModel<R> reset_model(active_set_model.range(),reset_model_centre,reset_map,spacial_order,smoothness);
      ARIADNE_LOG(6,"reset_model = "<<reset_model<<"\n");
      jump_set_model=compose(reset_model,active_set_model);
      ARIADNE_LOG(6,"jump_set_model = "<<jump_set_model<<"\n");
      ARIADNE_LOG(6,"jump_time_model = "<<jump_time_model<<"\n");
    }
  
    if(flow_type==GRAZING) {
      ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
      ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
    }

    reach_sets.adjoin(this->set(reach_set_model));
    if(flow_type!=MISSING) {
      ES active_set=this->set(active_set_model);
      ES jump_set=this->set(jump_set_model);
      Integer jump_steps=initial_steps+1;
      intermediate_sets.adjoin(active_set);
      intermediate_sets.adjoin(jump_set);
      working_sets.push(make_tuple(jump_steps,join(jump_set_model,jump_time_model)));
    } 
    if(flow_type==MISSING || flow_type==GRAZING) {
      ES final_set=this->set(final_set_model);
      Rational final_time_midpoint=Rational(final_time_model.range()[0].midpoint());
      if(final_time_midpoint >= maximum_time) {
        final_sets.adjoin(final_set);
      } else {
        intermediate_sets.adjoin(final_set);
        working_sets.push(make_tuple(initial_steps,join(final_set_model,final_time_model)));
      }
    }

  }
}


/*
template<class ES>
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_evolution(ESL& final,
           ESL& reachable,
           ESL& intermediate, 
           const Sys& sys,
           const ES& initial,
           const T& time,
           Semantics semantics,
           bool reach) const
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef typename traits<R>::interval_type I;

  typedef Box<R> Bx;

  const uint n=sys.state_space().dimension();
  const uint spacial_order=2;
  const uint temporal_order=4;
  const uint order=spacial_order+temporal_order;
  const uint smoothness=1;

  const FunctionInterface<R>& h=sys.impact_map();
  const FunctionInterface<R>& vf=sys.vector_field();
  const FunctionInterface<R>& guard_function=sys.guard_condition();

  
  uint verbosity=this->verbosity();

  TESL working;
  working.adjoin(TES(T(0),initial)); 

  while(working.size()!=0) {
    TES ts=working.pop();
    const ES& ws=ts.set();
    ARIADNE_LOG(6,"  ts="<<ts<<"\n"<<"\n");
    ARIADNE_ASSERT(ts.time()<=time);
    if(ts.time()==time) {
      final.adjoin(ts.set());
    } else if(radius(ts.set()) > maximum_enclosure_radius()) {
      if(semantics==upper_semantics) {
        this->adjoin_subdivision(working,ts);
        ++this->_profiler->subdivisions;
      } 
    } else {
      //working.adjoin(evolution_step(sys,ts));
    }
  }
}
*/


}  // namespace Ariadne
