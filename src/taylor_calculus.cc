/***************************************************************************
 *            taylor_calculus.cc
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
#include "pointer.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "taylor_model.h"
#include "taylor_variable.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "function.h"
#include "box.h"
#include "taylor_calculus.h"




namespace Ariadne {

class DegenerateCrossingException : public std::runtime_error {
  public:
    DegenerateCrossingException(const std::string& msg) : std::runtime_error(msg) { }
};

class NonInvertibleFunctionException  : public std::runtime_error {
  public:
    NonInvertibleFunctionException(const std::string& msg) : std::runtime_error(msg) { }
};


std::pair<Float, Vector<Interval> >
flow_bounds(FunctionInterface const& vf,
            Vector<Interval> const& d,
            Float const& hmax)
{

    ARIADNE_ASSERT(vf.argument_size()==d.size());

    // Set up constants of the method.
    // TODO: Better estimates of constants
    const double INITIAL_MULTIPLIER=2;
    const double MULTIPLIER=1.125;
    const double BOX_RADIUS_MULTIPLIER=1.03125;
    const uint EXPANSION_STEPS=8;
    const uint REDUCTION_STEPS=8;
    const uint REFINEMENT_STEPS=4;

    Vector<Interval> delta=(d-midpoint(d))*(BOX_RADIUS_MULTIPLIER-1);

    Float h=hmax;
    Float hmin=hmax/(1<<REDUCTION_STEPS);
    bool success=false;
    Vector<Interval> b,nb,df;
    Interval ih(0,h);
    while(!success) {
        ARIADNE_ASSERT(h>hmin);
        b=d+INITIAL_MULTIPLIER*ih*vf.evaluate(d)+delta;
        for(uint i=0; i!=EXPANSION_STEPS; ++i) {
            df=vf.evaluate(b);
            nb=d+delta+ih*df;
            if(subset(nb,b)) {
                success=true;
                break;
            } else {
                b=d+delta+MULTIPLIER*ih*df;
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

    for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
        b=nb;
        vfb=vf.evaluate(b);
        nb=d+delta+ih*vfb;
        ARIADNE_ASSERT_MSG(subset(nb,b),std::setprecision(20)<<"refinement "<<i<<": "<<nb<<" is not a inside of "<<b);
    }


    // Check result of operation
    // We use subset rather than inner subset here since the bound may touch
    ARIADNE_ASSERT(subset(nb,b));

    b=nb;

/*
    // Expand the resulting bound by a small constant for robustness
    volatile float new_bound;
    set_rounding_mode(upward);
    for(uint i=0; i!=b.size(); ++i) {
        new_bound=b[i].upper();
        new_bound+=1e-8;
        new_bound+=1e-8;
        b[i].u=new_bound;
    }
    set_rounding_mode(downward);
    for(uint i=0; i!=b.size(); ++i) {
        new_bound=b[i].lower();
        new_bound-=1e-8;
        new_bound-=1e-8;
        b[i].l=new_bound;
    }
    set_rounding_mode(to_nearest);
*/

    ARIADNE_ASSERT(subset(d,b));

    ARIADNE_ASSERT_MSG(subset(d+h*vf.evaluate(b),b),
        "d="<<d<<"\nh="<<h<<"\nf(b)="<<vf.evaluate(b)<<"\nd+hf(b)="<<Vector<Interval>(d+h*vf.evaluate(b))<<"\nb="<<b<<"\n");

    return std::make_pair(h,b);
}





TaylorCalculus::
TaylorCalculus()
    : _spacial_order(2),
      _temporal_order(4),
      _order(6),
      _smoothness(1),
      _maximum_step_size(1.0)
{
}




TaylorCalculus::SetModelType
TaylorCalculus::
reset_step(const FunctionType& map,
           const SetModelType& set_model) const
{
    // Direct computation from function
    return map.evaluate(set_model.models());
    // Indirect computation via model
    BoxType range=set_model.range();
    return Ariadne::apply(FunctionModelType(range,map),set_model);
}


TaylorCalculus::SetModelType
TaylorCalculus::
reset_step(const FlowModelType& function_model,
           const SetModelType& set_model) const
{
    ARIADNE_ASSERT(subset(set_model.range(),function_model.domain()));
    return Ariadne::unchecked_apply(function_model,set_model);
}



TaylorCalculus::SetModelType
TaylorCalculus::
integration_step(const FlowModelType& flow_model,
                 const SetModelType& initial_set_model,
                 const TimeModelType& integration_time_model) const
{
    //SetModelType set_step_model=join(initial_set_model, integration_time_model);
    SetModelType set_step_model(join(initial_set_model.models(),integration_time_model));
    ARIADNE_LOG(6,"set_step_model = "<<set_step_model<<"\n");
    //FIXME: Make sure that we never need to check subset of domain in flow range.
    //ARIADNE_ASSERT(subset(flow_model.domain()[flow_model.size()],integration_time_model.range()));
    //ARIADNE_ASSERT(subset(flow_model.domain(),set_step_model.range()));
    SetModelType final_set_model=unchecked_apply(flow_model,set_step_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
    return final_set_model;
}










TaylorCalculus::TimeModelType
TaylorCalculus::
reachability_time(const TimeModelType& initial_time_model,
                  const TimeModelType& final_time_model) const
{
    ARIADNE_ASSERT(initial_time_model.argument_size()==final_time_model.argument_size());
    uint ng=initial_time_model.argument_size();
    Interval unit_ivl(-1,1);

    TimeModelType expanded_initial_time_model=embed(initial_time_model,1u);
    TimeModelType expanded_final_time_model=embed(final_time_model,1u);

    TimeModelType time_interval_model=TimeModelType::scaling(1u,0u,Interval(0,1));
    TimeModelType expanded_time_interval_model=embed(ng,time_interval_model);
    TimeModelType expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);

    return expanded_reach_time_model;
}


TaylorCalculus::SetModelType
TaylorCalculus::
reachability_step(const FlowModelType& flow_model,
                  const SetModelType& initial_set_model,
                  const TimeType& initial_time,
                  const TimeType& final_time) const
{
    ARIADNE_ASSERT(initial_set_model.dimension()+1==flow_model.argument_size());

    TimeModelType reach_time_model(1);
    reach_time_model.set_value((final_time+initial_time)/2);
    reach_time_model.set_gradient(0u,(final_time-initial_time)/2);

    // FIXME: Embed set model correctly
    SetModelType expanded_timed_set_model=combine(initial_set_model.models(),reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    SetModelType reach_set_model=apply(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");

    return reach_set_model;
}

TaylorCalculus::SetModelType
TaylorCalculus::
reachability_step(const FlowModelType& flow_model,
                  const SetModelType& initial_set_model,
                  const TimeModelType& expanded_reach_time_model) const
{
    ARIADNE_ASSERT(initial_set_model.generators_size()+1==expanded_reach_time_model.argument_size());


    // FIXME: Embed set model correctly
    SetModelType expanded_initial_set_model=embed(initial_set_model.models(),1u);
    ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
    SetModelType expanded_timed_set_model=join(expanded_initial_set_model.models(),expanded_reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    SetModelType reach_set_model=apply(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");

    return reach_set_model;
}


TaylorCalculus::SetModelType
TaylorCalculus::
reachability_step(const FlowModelType& flow_model,
                  const SetModelType& initial_set_model,
                  const TimeModelType& initial_time_model,
                  const TimeModelType& final_time_model) const
{
    // Compute the reachable set
    // Need an extra independent variable to represent time
    uint ng=initial_set_model.generators_size();
    Interval ui(-1,+1);
    SetModelType expanded_initial_set_model=embed(initial_set_model.models(),1u);
    ARIADNE_LOG(6,"expanded_initial_set_model="<<expanded_initial_set_model<<"\n");
    TimeModelType expanded_initial_time_model=embed(initial_time_model,1u);
    ARIADNE_LOG(6,"expanded_initial_time_model="<<expanded_initial_time_model<<"\n");
    TimeModelType expanded_final_time_model=embed(final_time_model,1u);
    ARIADNE_LOG(6,"expanded_final_time_model="<<expanded_final_time_model<<"\n");

    TimeModelType expanded_time_interval_model=TimeModelType::scaling(ng+1,ng,Interval(0,1));
    ARIADNE_LOG(6,"expanded_time_interval_model="<<expanded_time_interval_model<<"\n");
    TimeModelType expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);

    ARIADNE_LOG(6,"expanded_reach_time_model="<<expanded_reach_time_model<<"\n");
    SetModelType expanded_timed_set_model=join(expanded_initial_set_model.models(),expanded_reach_time_model);
    ARIADNE_LOG(6,"expanded_timed_set_model="<<expanded_timed_set_model<<"\n");
    SetModelType reach_set_model=Ariadne::unchecked_apply(flow_model,expanded_timed_set_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");

    return reach_set_model;
}





// Compute the grazing time using bisections

TaylorCalculus::TimeModelType
TaylorCalculus::
crossing_time(const PredicateModelType& guard_model,
              const FlowModelType& flow_model,
              const SetModelType& initial_set_model) const
{
    uint dimension=flow_model.result_size();
    RealType minimum_time=flow_model.domain()[dimension].lower();
    RealType maximum_time=flow_model.domain()[dimension].upper();

    Interval flow_time_interval=flow_model.domain()[dimension];
    ARIADNE_ASSERT(subset(join(initial_set_model.range(),flow_time_interval),flow_model.domain()));
    
    //ARIADNE_ASSERT(minimum_time<=0);
    //ARIADNE_ASSERT(maximum_time>=0);
    ARIADNE_ASSERT(flow_model.argument_size()==flow_model.result_size()+1);
    ARIADNE_ASSERT(guard_model.argument_size()==flow_model.result_size());
    ARIADNE_ASSERT(initial_set_model.dimension()==flow_model.result_size());

    PredicateModelType hitting_model=unchecked_compose(guard_model,flow_model);
    ARIADNE_LOG(6,"hitting_model = "<<hitting_model<<"\n");
    PredicateModelType free_hitting_time_model;
    try {
        free_hitting_time_model=Ariadne::implicit(hitting_model);
    } catch(const ImplicitFunctionException& e) {
        throw DegenerateCrossingException(e.what());
    } catch(const std::runtime_error& e) {
        std::cerr<<e.what();
        throw e;
    }

    ARIADNE_LOG(6,"free_hitting_time_model = "<<free_hitting_time_model<<"\n");
    ARIADNE_ASSERT_MSG(free_hitting_time_model.argument_size()==initial_set_model.dimension(),free_hitting_time_model<<initial_set_model);
    TimeModelType hitting_time_model=Ariadne::apply(free_hitting_time_model,initial_set_model);
    ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
    Interval hitting_time_range=hitting_time_model.range();
    ARIADNE_LOG(6,"hitting_time_model = "<<hitting_time_model<<"\n");
    if(hitting_time_range.lower()<R(minimum_time) || hitting_time_range.upper()>R(maximum_time)) {
        ARIADNE_THROW(DegenerateCrossingException,"TaylorCalculus::crossing_time","Hitting time model "<<hitting_time_model<<" range does not lie in integration time range");
    }
    ARIADNE_ASSERT_MSG(hitting_time_model.argument_size()==initial_set_model.argument_size(),hitting_time_model<<initial_set_model);
    return hitting_time_model;
}


// Compute the grazing time using bisections

TaylorCalculus::IntervalType
TaylorCalculus::
touching_time_interval(const PredicateModelType& guard_model,
                       const FlowModelType& flow_model,
                       const SetModelType& initial_set_model) const
{
    //std::cerr<<"touching_time_interval"<<std::endl;
    ARIADNE_ASSERT(flow_model.result_size()+1==flow_model.argument_size());
    ARIADNE_ASSERT(guard_model.argument_size()==flow_model.result_size());

    uint dimension=guard_model.argument_size();
    ARIADNE_ASSERT(subset(initial_set_model.range(),guard_model.domain()));
    ARIADNE_ASSERT(subset(initial_set_model.range(),project(flow_model.domain(),range(0,dimension))));
    //ARIADNE_ASSERT_MSG(subset(flow_model.range(),guard_model.domain()),
    //    flow_model.range()<<"\n"<<guard_model.domain()<<"\n");

    RealType minimum_time=flow_model.domain()[dimension].lower();
    RealType maximum_time=flow_model.domain()[dimension].upper();

    ARIADNE_LOG(6,"\nminimum_time="<<minimum_time<<" maximum_time="<<maximum_time<<"\n");
    ARIADNE_ASSERT(minimum_time<=0);
    ARIADNE_ASSERT(maximum_time>=0);

    SetModelType final_set_model=this->integration_step(flow_model,initial_set_model,maximum_time);

    uint refinements=5;

    RealType lower_time=minimum_time;
    RealType upper_time=maximum_time;

    if(possibly(this->active(guard_model,initial_set_model))) {
        lower_time=minimum_time;
    } else {
        RealType min_lower_time=minimum_time;
        RealType max_lower_time=maximum_time;
        for(uint i=0; i!=refinements; ++i) {
            RealType new_lower_time=med_approx(min_lower_time,max_lower_time);
            if(possibly(this->active(guard_model,this->reachability_step(flow_model,initial_set_model,minimum_time,new_lower_time)))) {
                max_lower_time=new_lower_time;
            } else {
                min_lower_time=new_lower_time;
                lower_time=new_lower_time;
            }
        }
    }

    if(definitely(this->active(guard_model,initial_set_model)) || definitely(!this->active(guard_model,initial_set_model))) {
        upper_time=minimum_time;
    } else {
        RealType min_upper_time=minimum_time;
        RealType max_upper_time=maximum_time;
        for(uint i=0; i!=refinements; ++i) {
            RealType new_upper_time=med_approx(min_upper_time,max_upper_time);
            if(definitely(this->active(guard_model,this->integration_step(flow_model,initial_set_model,new_upper_time)))) {
                max_upper_time=new_upper_time;
                upper_time=new_upper_time;
            } else {
                min_upper_time=new_upper_time;
            }
        }
    }

    //FunctionModelType lower_time_model=FunctionModelType::constant(initial_set_model.domain(),initial_set_model.centre(),Vector<Float>(1u,lower_time),_order,_smoothness);
    //FunctionModelType upper_time_model=FunctionModelType::constant(initial_set_model.domain(),initial_set_model.centre(),Vector<Float>(1u,upper_time),_order,_smoothness);
    return Interval(lower_time,upper_time);

}





std::pair<Float, Vector<Interval> >
TaylorCalculus::flow_bounds(FunctionInterface const& vf,
                            Vector<Interval> const& r,
                            Float const& hmax,
                            Float const& dmax) const
{
    // Try to find a time h and a set b such that inside(r+Interval<R>(0,h)*vf(b),b) holds
    ARIADNE_LOG(6,"flow_bounds(Function,Box,Time hmax)\n");
    ARIADNE_LOG(7,"  r="<<r<<" hmax="<<hmax<<"\n");

    // Expand initial domain slightly if interior is nonempty
    
    static const Float EPS=1e-8;
    Vector<Interval> bx=r;
    for(uint i=0; i!=bx.size(); ++i) { if(bx[i].lower()==bx[i].upper()) { bx[i]+=Interval(-EPS,+EPS); } }

    return Ariadne::flow_bounds(vf,bx,hmax);
}


tribool
TaylorCalculus::
active(const PredicateModelType& guard_model, const SetModelType& set_model) const
{
    //ARIADNE_ASSERT_MSG(subset(set_model.range(),guard_model.domain()),
    //    "\nset_model.range()="<<set_model.range()<<"\nguard_model="<<guard_model<<"\nset_model="<<set_model);
    IntervalType range=Ariadne::unchecked_apply(guard_model,set_model).range();
    return this->_tribool(range);
}





TaylorCalculus::FunctionModelType
TaylorCalculus::map_model(FunctionInterface const& f, Vector<Interval> const& bx) const
{
    ARIADNE_ASSERT(f.argument_size()==bx.size());

    FunctionModelType map_model(bx,f);
    ARIADNE_LOG(6,"map_model = "<<map_model<<"\n");

    return map_model;
}



TaylorCalculus::FlowModelType
TaylorCalculus::flow_model(FunctionInterface const& vf, Vector<Interval> const& ibx, Float const& h, Vector<Interval> const& bb) const
{
    Vector<Interval> bx=ibx;

    ARIADNE_ASSERT(subset(bx+Interval(0,h)*vf.evaluate(bb),bb));

    // We need the initial box to have nonempty interior to be a valid domain for a function model,
    // so we expand slightly if one of the components fails this test
    for(uint i=0; i!=bx.size(); ++i) {
        static const Float EPS=1e-8;
        if(bx[i].lower()==bx[i].upper()) {
            //std::cerr<<"Warning: expanding initial box "<<bx<<std::endl;
            bx[i]+=Interval(-EPS,+EPS); }
    }

    ARIADNE_ASSERT(subset(bx+Interval(0,h)*vf.evaluate(bb),bb));

    TaylorModel::Accuracy* accuracy_raw_ptr=new TaylorModel::Accuracy;
    boost::shared_ptr<TaylorModel::Accuracy> accuracy(accuracy_raw_ptr);
    accuracy->_maximum_degree=this->_temporal_order;
    accuracy->_sweep_threshold=1e-10;
    FunctionModelType vector_field_model(bb,vf,accuracy);
    ARIADNE_LOG(6,"vector_field_model = "<<vector_field_model<<"\n");
    for(uint i=0; i!=vector_field_model.size(); ++i) {
        vector_field_model.models()[i].accuracy_ptr()->_maximum_degree+=this->_spacial_order;
    }

    const uint temporal_order=this->_temporal_order;
    // Use flow function on model type
    FlowModelType flow_model=Ariadne::unchecked_flow(vector_field_model,bx,Interval(0,h),temporal_order);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");

    return flow_model;
}



TaylorCalculus::PredicateModelType
TaylorCalculus::predicate_model(FunctionInterface const& g, Vector<Interval> const& bx) const
{
    //ARIADNE_DEPRECATED("TaylorCalculus::predicate_model(FunctionInterface,Vector<Interval>","Use ExpressionInterface instead");
    ARIADNE_ASSERT(g.argument_size()==bx.size());

    FunctionModelType predicate_model(bx,g);
    ARIADNE_LOG(6,"predicate_model = "<<predicate_model<<"\n");

    return predicate_model[0];
}

TaylorCalculus::PredicateModelType
TaylorCalculus::predicate_model(ExpressionInterface const& g, Vector<Interval> const& bx) const
{
    ARIADNE_ASSERT(g.argument_size()==bx.size());

    PredicateModelType predicate_model(bx,g);
    ARIADNE_LOG(6,"predicate_model = "<<predicate_model<<"\n");

    return predicate_model;
}

//! \brief A model for the constant time function \a t over the domain \a d.

TaylorCalculus::TimeModelType
TaylorCalculus::
time_model(const Float& t,
           const BoxType& d) const
{
    TimeModelType time_model=TimeModelType::constant(d.size(),t);
    ARIADNE_LOG(6,"time_model = "<<time_model<<"\n");

    return time_model;
}


//! \brief A model for the set f\a bx.

TaylorCalculus::SetModelType
TaylorCalculus::
set_model(const BoxType& bx) const
{
    SetModelType set_model(bx);
    return set_model;
 }


//! \brief A model for the enclosure \a es.

TaylorCalculus::SetModelType
TaylorCalculus::
set_model(const EnclosureType& es) const
{
    return es;
}


//! \brief An enclosure for the set model \a s.

TaylorCalculus::EnclosureType
TaylorCalculus::
enclosure(const SetModelType& s) const
{
    return s;
}






tribool
TaylorCalculus::disjoint(SetModelType const& set, BoxType const& bx) const
{
    return set.disjoint(bx);
}


TaylorCalculus::BoxType
TaylorCalculus::bounding_box(SetModelType const& set) const
{
    return set.bounding_box();
}


array<TaylorCalculus::SetModelType>
TaylorCalculus::subdivide(SetModelType const& set) const
{
    std::pair<SetModelType,SetModelType> split=set.split();
    array<SetModelType> result(2);
    result[0]=split.first;
    result[1]=split.second;
    return result;
}


TaylorCalculus::SetModelType
TaylorCalculus::simplify(SetModelType const&) const
{
    ARIADNE_NOT_IMPLEMENTED;
}




}  // namespace Ariadne
