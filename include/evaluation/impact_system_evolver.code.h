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

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "differentiation/sparse_differential.h"

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
    _profiler(new EvolutionProfiler)
{ }





template<class ES>
ApproximateTaylorModel<typename ES::real_type>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
model(const ES& set, ushort d) const
{
  uint n=set.dimension();
  uint m=set.number_of_generators();
  Vector<I> set_domain=set.domain().position_vectors();
  SparseDifferentialVector<A> set_expansion(n,m,d,
                                            Vector<A>(set.centre().position_vector()),
                                            reinterpret_cast<Matrix<A>const&>(set.generators()));
  ApproximateTaylorModel<R> set_model(set_domain,midpoint(set_domain),set_expansion);
  return set_model;
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
  make_lpair(set_centre,set_generators) = affine(model);
  ES enclosure_set(Point<I>(set_centre),set_generators);
  return enclosure_set;
}



  
template<class ES>
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
evolution(ESL& final_sets, ESL& reach_sets, ESL& intermediate_sets, const Sys& system, const ES& initial_set, const T& time) const
{
  uint verbosity=0;
  ARIADNE_LOG(5,__PRETTY_FUNCTION__);

  const uint n=system.state_space().dimension();
  const uint spacial_order=2;
  const uint temporal_order=4;
  const uint order=spacial_order+temporal_order;
  const uint smoothness=1;

  const FunctionInterface<R>& reset_map=system.impact_map();
  const FunctionInterface<R>& vector_field=system.vector_field();
  const FunctionInterface<R>& guard_function=system.guard_condition();

  stack< ApproximateTaylorModel<R> > working_sets;

  {
    // Set up initial timed set models
    ApproximateTaylorModel<R> initial_set_model=this->model(initial_set,spacial_order);
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    
    ApproximateTaylorModel<R> initial_time_model(initial_set_model.domain(),initial_set_model.centre(),
                                                 SparseDifferentialVector<A>::constant(1,n,order,Vector<A>(1,A(0))));
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    ApproximateTaylorModel<R> initial_timed_set_model=join(initial_set_model,initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push(initial_timed_set_model);
  }

  while(!working_sets.empty()) {
    const ApproximateTaylorModel<R> initial_timed_set_model=working_sets.pop();
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");

    const ApproximateTaylorModel<R> initial_set_model=project(initial_timed_set_model,range(0,n));
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    const ApproximateTaylorModel<R> initial_time_model=project(initial_timed_set_model,range(n,n+1));
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");

    ++this->_profiler->time_steps;
    IVec flow_bounds; T h;
    make_lpair(h,flow_bounds)=this->flow_bounds(vector_field,initial_set_model.range(),T(0.25));
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
    ApproximateTaylorModel<R> final_set_model;
    ApproximateTaylorModel<R> final_time_model;
    ApproximateTaylorModel<R> reach_set_model;
    ApproximateTaylorModel<R> reach_time_model;
    
    if(possibly(hitting)) { // hitting
      SparseDifferentialVector<A> guard_expansion=guard_function.expansion(midpoint(flow_bounds),order);
      ARIADNE_LOG(6,"guard_function = "<<guard_function<<"\n");
      ARIADNE_LOG(6,"guard_centre = "<<midpoint(flow_bounds)<<"\n");
      ARIADNE_LOG(6,"guard_expansion = "<<guard_expansion<<"\n");

      guard_model=ApproximateTaylorModel<R>(flow_bounds,guard_function,order,smoothness);
      ARIADNE_LOG(6,"guard_model = "<<guard_model<<"\n");
      ApproximateTaylorModel<R> hitting_model=compose(guard_model,flow_model);
      ARIADNE_LOG(6,"hitting_model = "<<hitting_model<<"\n");
      ApproximateTaylorModel<R> hitting_time_model=implicit(hitting_model);
      ARIADNE_LOG(6,"hitting_time_model = "<<integration_time_model<<"\n");
      integration_time_model=compose(hitting_time_model,initial_set_model);
      ARIADNE_LOG(6,"integration_time_model = "<<integration_time_model<<"\n");
      final_time_model=initial_time_model+integration_time_model;
      ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
    } else { // no hitting
      I initial_time_range=initial_time_model.range()[0];
      ARIADNE_ASSERT(A(initial_time_range.lower())+step_size > A(initial_time_range.upper()));
      A final_time=A(initial_time_range.lower())+step_size;
      SparseDifferentialVector<A> final_time_expansion=SparseDifferentialVector<A>::constant(1u,n,order,Vector<A>(1u,final_time));
      final_time_model=ApproximateTaylorModel<R>(initial_set_model.domain(),initial_set_model.centre(),final_time_expansion);
      ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
      integration_time_model=final_time_model-initial_time_model;
      ARIADNE_LOG(6,"integration_time_model = "<<integration_time_model<<"\n");
    }

    //SparseDifferentialVector<A> identity=SparseDifferentialVector<A>::variable(n,n,order,Vector<A>(n,A(0)));
    //ApproximateTaylorModel<R> hitting_model(hitting_time_model.domain(),hitting_time_model.centre(),
    //                                          join(identity,hitting_time_model.expansion()));
    //ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<std::cerr;
      
    //final_set_model=compose(flow_model,hitting_model<<"\n");
    //final_time_model=initial_time_model+hitting_time_model;

    SparseDifferentialVector<A> step_expansion(n+1,n,order);
    step_expansion[n]=integration_time_model.expansion()[0];
    for(uint i=0; i!=n; ++i) {
      step_expansion[i].gradient(i)=1.0;
    }
    ARIADNE_LOG(6,"step_expansion = "<<step_expansion<<"\n");
    ApproximateTaylorModel<R> time_flow_model(Vector<I>(n,I(-1,1)),Vector<R>(n,R(0)),step_expansion);
    ARIADNE_LOG(6,"time_flow_model = "<<time_flow_model<<"\n");
    ApproximateTaylorModel<R> integration_step_model=compose(flow_model,time_flow_model);
    ARIADNE_LOG(6,"integration_step_model = "<<integration_step_model<<"\n");
    final_set_model=compose(integration_step_model,initial_set_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
  
    ApproximateTaylorModel<R> set_step_model=join(initial_set_model, integration_time_model);
    ARIADNE_LOG(6,"set_step_model = "<<set_step_model<<"\n");
    final_set_model=compose(flow_model,set_step_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");

    final_time_model=initial_time_model+integration_time_model;
    ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
    
    // Compute the reachable set
    // Need an extra independent variable to represent time
    SparseDifferentialVector<A> space_expansion(n,n+1,order);
    for(uint i=0; i!=n; ++i) { space_expansion[i][i]=1; }
    //ARIADNE_LOG(6,"space_expansion = "<<space_expansion<<"\n");
    SparseDifferentialVector<A> time_expansion(1,n+1,order);
    time_expansion[0].value()=A(h)/2;
    time_expansion[0].gradient(n)=A(h)/2;
    //ARIADNE_LOG(6,"time_expansion = "<<time_expansion<<"\n");
    ApproximateTaylorModel<R> space_expansion_model(Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),space_expansion);
    ARIADNE_LOG(6,"space_expansion_model = "<<space_expansion_model<<"\n");
    ARIADNE_LOG(6,"intial_set_model = "<<initial_set_model<<"\n");
    ApproximateTaylorModel<R> expanded_set_model=compose(initial_set_model,space_expansion_model);
    ARIADNE_LOG(6,"expanded_set_model = "<<expanded_set_model<<"\n");
    ApproximateTaylorModel<R> expanded_time_model(Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),time_expansion);
    ARIADNE_LOG(6,"expanded_time_model = "<<expanded_time_model<<"\n");
    ApproximateTaylorModel<R> expanded_timed_set_model=join(expanded_set_model,expanded_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model = "<<expanded_timed_set_model<<"\n");
    
    reach_set_model=compose(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    

    ApproximateTaylorModel<R> jump_set_model;
    if(possibly(hitting)) {
      ApproximateTaylorModel<R> reset_model(final_set_model.range(),Vector<R>(final_set_model.expansion().value()),reset_map,spacial_order,smoothness);
      ARIADNE_LOG(6,"reset_model = "<<reset_model<<"\n");
      jump_set_model=compose(reset_model,final_set_model);
      ARIADNE_LOG(6,"jump_set_model = "<<jump_set_model<<"\n");
    }
  
    // Compute reach set from model
    ES reach_set=this->set(reach_set_model);
  
    // Compute final set from model
    ES final_set=this->set(final_set_model);
    T t=T(final_time_model.expansion().value()[0]._value);
    
    reach_sets.adjoin(reach_set);
    if(possibly(hitting)) {
      ES jump_set=this->set(jump_set_model);
      intermediate_sets.adjoin(final_set);
      intermediate_sets.adjoin(jump_set);
      working_sets.push(join(jump_set_model,final_time_model));
    } else if(t>=time) {
      final_sets.adjoin(final_set);
    } else {
      intermediate_sets.adjoin(final_set);
      working_sets.push(join(final_set_model,final_time_model));
    }
 
  }

}


template<class ES>
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_evolution(ESL& final,
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



}  // namespace Ariadne
