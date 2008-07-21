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
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_evolution(ESL& final_sets, 
           ESL& reach_sets, 
           ESL& intermediate_sets, 
           const Sys& system, 
           const ES& initial_set, 
           const T& time, 
           Semantics semantics, 
           bool reach) const
{
  uint verbosity=0;
  ARIADNE_LOG(5,__PRETTY_FUNCTION__);
  assert(semantics==upper_semantics);

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
    
    ApproximateTaylorModel<R> initial_time_model
      =ApproximateTaylorModel<R>::constant(initial_set_model.domain(),initial_set_model.centre(),
                                           Vector<A>(1,A(0)),order,smoothness);
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
      ARIADNE_LOG(6,"guard_function = "<<guard_function<<"\n");

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
      final_time_model=ApproximateTaylorModel<R>::constant(initial_set_model.domain(),initial_set_model.centre(),Vector<A>(1u,final_time),order,smoothness);
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

    ApproximateTaylorModel<R> identity_model=ApproximateTaylorModel<R>::identity(integration_time_model.domain(),integration_time_model.centre(),order,smoothness);
    ARIADNE_LOG(6,"identity_model = "<<identity_model<<"\n");
    ApproximateTaylorModel<R> time_flow_model=join(identity_model,integration_time_model);
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
    
    /*
    SparseDifferentialVector<A> space_expansion(n,n+1,order);
    for(uint i=0; i!=n; ++i) { space_expansion[i][i]=1; }
    ApproximateTaylorModel<R> space_expansion_model(Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),space_expansion);
    ARIADNE_LOG(6,"space_expansion_model = "<<space_expansion_model<<"\n");
    ARIADNE_LOG(6,"intial_set_model = "<<initial_set_model<<"\n");
    ApproximateTaylorModel<R> expanded_set_model=compose(initial_set_model,space_expansion_model);
    ARIADNE_LOG(6,"expanded_set_model = "<<expanded_set_model<<"\n");
    SparseDifferentialVector<A> time_expansion(1u,n+1,order);
    A hh=A(h)/2; time_expansion[0]=hh; time_expansion[0][n]=hh;
    ApproximateTaylorModel<R> expanded_time_model(Vector<I>(n+1,I(-1,1)),Vector<R>(n+1,R(0)),time_expansion);
    ARIADNE_LOG(6,"expanded_time_model = "<<expanded_time_model<<"\n");
    ApproximateTaylorModel<R> expanded_timed_set_model=join(expanded_set_model,expanded_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model = "<<expanded_timed_set_model<<"\n");
    */
    
    ApproximateTaylorModel<R> time_model;
    ApproximateTaylorModel<R> expanded_timed_set_model;

    
    //Currently, can only convert model with domain the unit box to an affine model,
    //so can't use real time domain
    //Vector<I> time_domain(1,I(0u,h));
    //time_model=ApproximateTaylorModel<R>::identity(time_domain,midpoint(time_domain),order,smoothness);

    A hh=A(h)/2; 
    Vector<I> unit_interval(1,I(-1,1));
    Vector<A> time_centre(1); time_centre[0]=hh;
    Matrix<A> time_generator(1,1); time_generator[0][0]=hh;
    time_model=ApproximateTaylorModel<R>::affine(unit_interval,midpoint(unit_interval),time_centre,time_generator,order,smoothness);
    
    ARIADNE_LOG(6,"time_model = "<<time_model<<"\n");
    expanded_timed_set_model=combine(initial_set_model,time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model = "<<expanded_timed_set_model<<"\n");

    reach_set_model=compose(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    

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
    //T t=T(final_time_model.expansion().value()[0]._value);
    T t=T(midpoint(final_time_model.range()[0]));

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
