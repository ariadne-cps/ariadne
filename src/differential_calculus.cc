/***************************************************************************
 *            differential_calculus.cc
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
 
#include <iomanip>

#include "macros.h"
#include "logging.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "sparse_differential.h"
#include "function.h"
#include "box.h"
#include "differential_calculus.h"

#include "approximate_taylor_model.h"



namespace Ariadne {

class DegenerateCrossingException { };
class NonInvertibleFunctionException { };



  
template<class Mdl>
DifferentialCalculus<Mdl>::
DifferentialCalculus() 
    : _spacial_order(2),
      _temporal_order(4),
      _order(6),
      _smoothness(1)
{ 
}
  

template<class Mdl>
typename DifferentialCalculus<Mdl>::SetModelType
DifferentialCalculus<Mdl>::
reset_step(const FlowModelType& function_model, 
           const SetModelType& set_model) const
{
    return Ariadne::compose(function_model,set_model);
}


template<class Mdl>
typename DifferentialCalculus<Mdl>::SetModelType
DifferentialCalculus<Mdl>::
integration_step(const FlowModelType& flow_model, 
                 const SetModelType& initial_set_model, 
                 const TimeModelType& integration_time_model) const
{
    Mdl set_step_model=join(initial_set_model, integration_time_model);
    ARIADNE_LOG(6,"set_step_model = "<<set_step_model<<"\n");
    Mdl final_set_model=compose(flow_model,set_step_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");

    return final_set_model;
}







template<class Mdl>
typename DifferentialCalculus<Mdl>::SetModelType
DifferentialCalculus<Mdl>::
reachability_step(const FlowModelType& flow_model, 
                  const SetModelType& initial_set_model, 
                  const TimeModelType& expanded_reach_time_model) const
{
    ARIADNE_ASSERT(initial_set_model.argument_size()+1==expanded_reach_time_model.argument_size());

    // Compute the reachable set
    // Need an extra independent variable to represent time
    uint ng=initial_set_model.argument_size();
  
    ModelType expanded_initial_set_model=embed(initial_set_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
    ModelType expanded_timed_set_model=join(expanded_initial_set_model,expanded_reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    ModelType reach_set_model=compose(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
  
    return reach_set_model;
}


template<class Mdl>
typename DifferentialCalculus<Mdl>::SetModelType
DifferentialCalculus<Mdl>::
reachability_step(const FlowModelType& flow_model, 
                  const SetModelType& initial_set_model, 
                  const TimeModelType& initial_time_model, 
                  const TimeModelType& final_time_model) const
{
    // Compute the reachable set
    // Need an extra independent variable to represent time
    uint ng=initial_set_model.argument_size();
  
    ModelType expanded_initial_set_model=embed(initial_set_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
    ModelType expanded_initial_time_model=embed(initial_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_initial_time_model="<<expanded_initial_time_model<<"\n");
    ModelType expanded_final_time_model=embed(final_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
    ARIADNE_LOG(6,"expanded_final_time_model="<<expanded_final_time_model<<"\n");
  
    ModelType time_interval_model=ModelType::affine(I(-1,1),R(0),A(0.5),A(0.5),_order,_smoothness);
    ARIADNE_LOG(6,"time_interval_time_model="<<time_interval_model<<"\n");
    ModelType expanded_time_interval_model=embed(time_interval_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),ng);
    ARIADNE_LOG(6,"expanded_time_interval_model="<<expanded_time_interval_model<<"\n");
    ModelType expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);
    ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
    ModelType expanded_timed_set_model=join(expanded_initial_set_model,expanded_reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    ModelType reach_set_model=compose(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
  
    return reach_set_model;
}





// Compute the grazing time using bisections
template<class Mdl>
typename DifferentialCalculus<Mdl>::ModelType
DifferentialCalculus<Mdl>::
crossing_time(const ModelType& guard_model,
              const ModelType& flow_model, 
              const ModelType& initial_set_model) const
{
    uint dimension=flow_model.result_size();
    RealType minimum_time=flow_model.domain()[dimension].lower(); 
    RealType maximum_time=flow_model.domain()[dimension].upper(); 

    //ARIADNE_ASSERT(minimum_time<=0);
    //ARIADNE_ASSERT(maximum_time>=0);
    ARIADNE_ASSERT(flow_model.argument_size()==flow_model.result_size()+1);
    ARIADNE_ASSERT(guard_model.argument_size()==flow_model.result_size());
    ARIADNE_ASSERT(initial_set_model.result_size()==flow_model.result_size());

    ModelType hitting_model=compose(guard_model,flow_model);
    ARIADNE_LOG(6,"hitting_model = "<<hitting_model<<"\n");
    ModelType free_hitting_time_model;
    try {
        free_hitting_time_model=implicit(hitting_model); 
    } catch(NonInvertibleFunctionException) {
        throw DegenerateCrossingException();
    }
    ARIADNE_LOG(6,"free_hitting_time_model = "<<free_hitting_time_model<<"\n");
    ModelType hitting_time_model=compose(free_hitting_time_model,initial_set_model);
    ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
    Interval hitting_time_range=hitting_time_model.range()[0];
    ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
    if(hitting_time_range.lower()<R(minimum_time) || hitting_time_range.upper()>R(maximum_time)) {
        throw DegenerateCrossingException();
    }

    return hitting_time_model;
}


// Compute the grazing time using bisections
template<class Mdl>
Interval
DifferentialCalculus<Mdl>::
touching_time_interval(const ModelType& guard_model, 
                       const ModelType& flow_model, 
                       const ModelType& initial_set_model) const
{
    //uint old_verbosity = verbosity;
    //verbosity = 4;
    ARIADNE_LOG(3,"DifferentialCalculus<Mdl>::touching_time_interval(...)\n");

    ARIADNE_ASSERT(flow_model.result_size()+1==flow_model.argument_size());
    ARIADNE_ASSERT(guard_model.argument_size()==flow_model.result_size());
    ARIADNE_ASSERT(guard_model.result_size()==1u);

    uint dimension=guard_model.argument_size();
    RealType minimum_time=flow_model.domain()[dimension].lower(); 
    RealType maximum_time=flow_model.domain()[dimension].upper(); 
  
    ARIADNE_LOG(4,"\nminimum_time="<<minimum_time<<" maximum_time="<<maximum_time<<"\n");
    ARIADNE_ASSERT(minimum_time<=0);
    ARIADNE_ASSERT(maximum_time>=0);

    if(minimum_time<=0.0) minimum_time = 0.0;

    ARIADNE_LOG(4,"flow_model="<<flow_model<<"\n");
    ARIADNE_LOG(4,"guard_model="<<guard_model<<"\n");
    ARIADNE_LOG(4,"initial_set_model="<<initial_set_model<<"\n");

    ModelType maximum_time_model=this->time_model(maximum_time,initial_set_model.domain());
    ModelType final_set_model=this->integration_step(flow_model,initial_set_model,maximum_time_model);

    ARIADNE_LOG(4,"final_set_model="<<final_set_model<<"\n");    

    uint refinements=5;

    RealType lower_time=minimum_time;
    RealType upper_time=maximum_time;
  
    if(!possibly(this->active(guard_model,initial_set_model))) {
        RealType min_lower_time=minimum_time;
        RealType max_lower_time=maximum_time;
        for(uint i=0; i!=refinements; ++i) {
            RealType new_lower_time=med_approx(min_lower_time,max_lower_time);
            ARIADNE_LOG(4,"step "<<i<<", testing time "<<new_lower_time<<"\n");
            if(possibly(this->active(guard_model,this->reachability_step(flow_model,initial_set_model,minimum_time,new_lower_time)))) {
                max_lower_time=new_lower_time;
            } else {
                min_lower_time=new_lower_time;
                lower_time=new_lower_time;
            }
        }
        
    }
    ARIADNE_LOG(4,"lower_time="<<lower_time<<" upper_time="<<upper_time<<"\n");


    if(!possibly(this->active(guard_model,final_set_model))) {
        RealType min_upper_time=lower_time;
        RealType max_upper_time=maximum_time;
        for(uint i=0; i!=refinements; ++i) {
            RealType new_upper_time=med_approx(min_upper_time,max_upper_time);
            ARIADNE_LOG(4,"step "<<i<<", testing time "<<new_upper_time<<"\n");
            if(possibly(this->active(guard_model,this->integration_step(flow_model,initial_set_model,new_upper_time)))) {
                min_upper_time=new_upper_time;
            } else {
                max_upper_time=new_upper_time;
                upper_time=new_upper_time;
            }
        }
    }
    ARIADNE_LOG(4,"lower_time="<<lower_time<<" upper_time="<<upper_time<<"\n");

    //verbosity=old_verbosity;
    //ModelType lower_time_model=ModelType::constant(initial_set_model.domain(),initial_set_model.centre(),Vector<Float>(1u,lower_time),_order,_smoothness);
    //ModelType upper_time_model=ModelType::constant(initial_set_model.domain(),initial_set_model.centre(),Vector<Float>(1u,upper_time),_order,_smoothness);
    return Interval(lower_time,upper_time);
}




template<class Mdl>
std::pair<Float, Vector<Interval> >
DifferentialCalculus<Mdl>::flow_bounds(FunctionInterface const& vf, 
                                   Vector<Interval> const& r, 
                                   Float const& hmax, 
                                   Float const& dmax) const
{ 
    // Try to find a time h and a set b such that inside(r+Interval<R>(0,h)*vf(b),b) holds
    ARIADNE_LOG(6,"flow_bounds(Function,Box,Time hmax)\n");
    ARIADNE_LOG(7,"  r="<<r<<" hmax="<<hmax<<"\n");
  
    ARIADNE_ASSERT(vf.argument_size()==r.size());

    // Set up constants of the method.
    // TODO: Better estimates of constants
    const Float INITIAL_MULTIPLIER=2;
    const Float MULTIPLIER=1.125;
    //const Float BOX_RADIUS_MULTIPLIER=1.03125;
    const uint EXPANSION_STEPS=8;
    const uint REDUCTION_STEPS=8;
    const uint REFINEMENT_STEPS=4;
  
    Vector<Interval> delta=r-midpoint(r);
  
    Float h=hmax;
    Float hmin=hmax/(1<<REDUCTION_STEPS);
    bool success=false;
    Vector<Interval> b,nb,df;
    Interval ih(0,h);
    while(!success) {
        ARIADNE_ASSERT(h>hmin);
        b=r+INITIAL_MULTIPLIER*ih*vf.evaluate(r)+delta;
        for(uint i=0; i!=EXPANSION_STEPS; ++i) {
            df=vf.evaluate(b);
            nb=r+ih*df;
            ARIADNE_LOG(9,"  h="<<h<<" b="<<b<<" vf="<<vf.evaluate(b)<<" nb="<<nb<<"\n");
            if(subset(nb,b)) {
                success=true;
                break;
            } else {
                b=r+MULTIPLIER*ih*df+delta;
            }
        }
        if(!success) {
            h/=2;
            ih=Interval(0,h);
        }
    }

    ARIADNE_ASSERT(subset(nb,b));
  
    Vector<Interval> vfb;
    vfb=vf.evaluate(b);
    ARIADNE_LOG(9,std::setprecision(20)<<"\n\n   b="<<b<<" vf="<<vfb<<"\n  nb="<<nb<<"\n");

    for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
        b=nb;
        vfb=vf.evaluate(b);
        nb=r+ih*vfb;
        ARIADNE_LOG(9,std::setprecision(20)<<"   b="<<b<<" vf="<<vfb<<"\n  nb="<<nb<<"\n");
        ARIADNE_ASSERT_MSG(subset(nb,b),std::setprecision(20)<<"refinement "<<i<<": "<<nb<<" is not a inside of "<<b);
    }
  
    // Check result of operation
    // We use "possibly" here since the bound may touch 
    ARIADNE_ASSERT(subset(nb,b));
  
    ARIADNE_LOG(7,"  h="<<h<<" b="<<nb<<"\n");

    return std::make_pair(h,nb);
}



template<class Mdl>
tribool
DifferentialCalculus<Mdl>::
active(const PredicateModelType& guard_model, const SetModelType& set_model) const
{
    BoxType range=compose(guard_model,set_model).range();
 /*
    uint dimension=range.size();
    for(uint i = 0; i < dimension ; i++) {
      tribool test = this->_tribool(range[i]);
      if(!possibly(test)) return false;   // at least one component is negative
      if(indeterminate(test)) return indeterminate;   // one component is indeterminate
    }
    // if the for loop is completed, all components are positive and the guard is true
    return true;
*/
    return this->_tribool(range[0]);
}




template<class Mdl>
Mdl
DifferentialCalculus<Mdl>::map_model(FunctionInterface const& f, Vector<Interval> const& bx) const
{ 
    ARIADNE_ASSERT(f.argument_size()==bx.size());

    Mdl map_model(bx,f,_order,_smoothness);
    ARIADNE_LOG(6,"map_model = "<<map_model<<"\n");

    return map_model;
}


template<class Mdl>
Mdl
DifferentialCalculus<Mdl>::flow_model(FunctionInterface const& vf, Vector<Interval> const& bx, Float const& h, Vector<Interval> const& bb) const
{ 
    Mdl vector_field_model(bb,vf,_order,_smoothness);
    ARIADNE_LOG(6,"vector_field_model = "<<vector_field_model<<"\n");
  
    // Use flow function on model type
    Mdl flow_model=Ariadne::flow(vector_field_model);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");

    return flow_model;
}


template<class Mdl>
Mdl
DifferentialCalculus<Mdl>::predicate_model(FunctionInterface const& g, Vector<Interval> const& bx) const
{ 
    ARIADNE_ASSERT(g.result_size()==1);
    ARIADNE_ASSERT(g.argument_size()==bx.size());

    Mdl predicate_model(bx,g,_order,_smoothness);
    ARIADNE_LOG(6,"predicate_model = "<<predicate_model<<"\n");

    return predicate_model;
}

//! \brief A model for the constant time function \a t over the domain \a d.
template<class Mdl>
Mdl
DifferentialCalculus<Mdl>::
time_model(const Float& t, 
           const BoxType& bx) const
{
    Mdl time_model=Mdl::constant(bx,midpoint(bx),Vector<Float>(1u,t),_order,_smoothness);
    ARIADNE_LOG(6,"time_model = "<<time_model<<"\n");

    return time_model;
}


//! \brief A model for the set f\a bx.
template<class Mdl>
Mdl
DifferentialCalculus<Mdl>::
set_model(const BoxType& bx) const
{
    SetModelType set_model(bx,IdentityFunction(bx.size()),this->_order,this->_smoothness);
    return set_model;
 }


//! \brief A model for the enclosure \a es.
template<class Mdl>
Mdl
DifferentialCalculus<Mdl>::
set_model(const EnclosureType& es) const
{
    return es;
}


//! \brief An enclosure for the set model \a s.
template<class Mdl>
typename DifferentialCalculus<Mdl>::EnclosureType
DifferentialCalculus<Mdl>::
enclosure(const SetModelType& s) const
{
    return s;
}





template<class Mdl>
tribool
DifferentialCalculus<Mdl>::disjoint(Mdl const&, Vector<Interval> const&) const
{ 
    ARIADNE_NOT_IMPLEMENTED;
}

template<class Mdl>
Vector<Interval>
DifferentialCalculus<Mdl>::bounding_box(Mdl const&) const
{ 
    ARIADNE_NOT_IMPLEMENTED;
}

template<class Mdl>
array<Mdl> 
DifferentialCalculus<Mdl>::subdivide(Mdl const&) const
{ 
    ARIADNE_NOT_IMPLEMENTED;
}

template<class Mdl>
Mdl 
DifferentialCalculus<Mdl>::simplify(Mdl const&) const
{ 
    ARIADNE_NOT_IMPLEMENTED;
}



template class DifferentialCalculus<ApproximateTaylorModel>;

}  // namespace Ariadne
