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
_reachability_step(const ATM& flow_model, const ATM& initial_set_model, const ATM& integration_time_model) const
{

    // Compute the reachable set
    // Need an extra independent variable to represent time
    uint n=flow_model.result_size();
    
    ApproximateTaylorModel<R> expanded_initial_set_model=embed(initial_set_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
    ApproximateTaylorModel<R> expanded_integration_time_model=embed(integration_time_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_integration_time_model="<<expanded_integration_time_model<<"\n");
    
    ApproximateTaylorModel<R> reach_time_model=ApproximateTaylorModel<R>::affine(I(-1,1),R(0),A(0.5),A(0.5),order,smoothness);
    ARIADNE_LOG(6,"reach_time_model="<<reach_time_model<<"\n");
    ApproximateTaylorModel<R> expanded_reach_time_model=embed(reach_time_model,Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),n);
    ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
    expanded_reach_time_model=expanded_integration_time_model*expanded_reach_time_model;
    ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
    ApproximateTaylorModel<R> expanded_timed_set_model=join(expanded_initial_set_model,expanded_reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    ApproximateTaylorModel<R> reach_set_model=compose(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    
    return reach_set_model;
}


template<class ES>
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_evolution(ESL& final_sets, 
           ESL& reach_sets, 
           ESL& intermediate_sets, 
           const Sys& system, 
           const ES& initial_set, 
           const T& maximum_time, 
           Semantics semantics, 
           bool reach) const
{
  ARIADNE_LOG(5,__PRETTY_FUNCTION__);
  assert(semantics==upper_semantics);

  Integer max_steps=maximum_time.steps();
  A max_time=A(maximum_time.time());

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
    ARIADNE_LOG(6,"initial_steps = "<<initial_steps<<"\n");
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");

    const ApproximateTaylorModel<R> initial_set_model=project(initial_timed_set_model,range(0,n));
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    const ApproximateTaylorModel<R> initial_time_model=project(initial_timed_set_model,range(n,n+1));
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");

    ++this->_profiler->time_steps;
    IVec flow_bounds; Rational h;
    make_lpair(h,flow_bounds)=this->flow_bounds(vector_field,initial_set_model.range(),Rational(0.25));
    this->_profiler->total_stepping_time+=h;
    this->_profiler->minimum_time_step=std::min(h,this->_profiler->minimum_time_step);
    A step_size=A(h);
    ARIADNE_LOG(6,"flow_bounds = "<<flow_bounds<<"\nstep_size = "<<step_size<<"\n"<<"\n");
    
    ApproximateTaylorModel<R> vector_field_model(flow_bounds,vector_field,order,smoothness);
    ARIADNE_LOG(6,"vector_field_model = "<<vector_field_model<<"\n");
    ApproximateTaylorModel<R> flow_model=Ariadne::flow(vector_field_model);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");
    
    IVec guard_value = guard_function.evaluate(flow_bounds);
    ARIADNE_LOG(6,"guard_value = "<<guard_value<<"\n");

    tribool hitting;
    if(guard_value <= 0) {
      hitting=false;
    } else if(guard_value >= 0) {
      hitting=true;
    } else {
      hitting=indeterminate;
    }
    ARIADNE_LOG(6,std::boolalpha<<"hitting = "<<hitting<<std::endl<<"\n");

    ApproximateTaylorModel<R> integration_time_model;
    ApproximateTaylorModel<R> guard_model;
    ApproximateTaylorModel<R> hitting_time_model;
    ApproximateTaylorModel<R> final_set_model;
    ApproximateTaylorModel<R> final_time_model;
    ApproximateTaylorModel<R> reach_set_model;
    ApproximateTaylorModel<R> reach_time_model;
    
    if(possibly(hitting)) { // hitting
      ARIADNE_LOG(6,"guard_function = "<<guard_function<<"\n");

      guard_model=ApproximateTaylorModel<R>(flow_bounds,guard_function,order,smoothness);
      ARIADNE_LOG(6,"guard_model = "<<guard_model<<"\n");
      ApproximateTaylorModel<R> hitting_model=compose(guard_model,flow_model);
      ARIADNE_LOG(6,"hitting_model = "<<hitting_model<<"\n");
      try {
        // Transverse crossing
        hitting_time_model=implicit(hitting_model); 
        ARIADNE_LOG(6,"hitting_time_model = "<<integration_time_model<<"\n");
        integration_time_model=compose(hitting_time_model,initial_set_model);
        ARIADNE_LOG(6,"integration_time_model = "<<integration_time_model<<"\n");
        final_time_model=initial_time_model+integration_time_model;
        ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
      } catch(NonInvertibleFunctionException) {
        // Non-transverse crossing
        ARIADNE_LOG(6,"\nnon-transverse crossing\n\n");
        A minimum_first_hitting_time = step_size;
        
        
        
      }
    } else { // no hitting
      I initial_time_range=initial_time_model.range()[0];
      ARIADNE_ASSERT(A(initial_time_range.lower())+step_size > A(initial_time_range.upper()));
      A final_time=A(initial_time_range.lower())+step_size;
      if(final_time>max_time) { final_time=max_time; }
      final_time_model=ApproximateTaylorModel<R>::constant(initial_set_model.domain(),initial_set_model.centre(),Vector<A>(1u,final_time),order,smoothness);
      ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
      integration_time_model=final_time_model-initial_time_model;
      ARIADNE_LOG(6,"integration_time_model = "<<integration_time_model<<"\n");
    }


    // Compute the evolved set
    final_set_model=this->_integration_step(flow_model,initial_set_model,integration_time_model);
   
    // Compute the reachable set
    reach_set_model=this->_reachability_step(flow_model,initial_set_model,integration_time_model);
    

    // Compute the jump set
    ApproximateTaylorModel<R> jump_set_model;
    if(possibly(hitting)) {
      Vector<A> reset_model_centre(final_set_model.evaluate(Vector<A>(final_set_model.centre())));
      ApproximateTaylorModel<R> reset_model(final_set_model.range(),reset_model_centre,reset_map,spacial_order,smoothness);
      ARIADNE_LOG(6,"reset_model = "<<reset_model<<"\n");
      jump_set_model=compose(reset_model,final_set_model);
      ARIADNE_LOG(6,"jump_set_model = "<<jump_set_model<<"\n");
    }
  

    // Compute reach set from model
    ES reach_set=this->set(reach_set_model);
    // Compute final set from model
    ES final_set=this->set(final_set_model);

    // Estimate the time
    Rational t=Rational(midpoint(final_time_model.range()[0]));
    ARIADNE_LOG(6,"\ntime = "<<A(t)<<"\n\n");
    
    reach_sets.adjoin(reach_set);
    if(possibly(hitting) && initial_steps<max_steps) {
      ES jump_set=this->set(jump_set_model);
      Integer jump_steps=initial_steps+1;
      intermediate_sets.adjoin(final_set);
      intermediate_sets.adjoin(jump_set);
      working_sets.push(make_tuple(jump_steps,join(jump_set_model,final_time_model)));
    } else if(t>=maximum_time.time()) {
      final_sets.adjoin(final_set);
    } else {
      intermediate_sets.adjoin(final_set);
      working_sets.push(make_tuple(initial_steps,join(final_set_model,final_time_model)));
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
