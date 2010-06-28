/***************************************************************************
 *            hybrid_enclosure.cc
 *
 *  Copyright  2009-10  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "operators.h"
#include "function.h"
#include "constraint.h"
#include "propagator.h"
#include "taylor_model.h"
#include "taylor_function.h"
#include "box.h"
#include "grid_set.h"
#include "hybrid_time.h"
#include "discrete_event.h"
#include "discrete_location.h"

#include "linear_programming.h"
#include "nonlinear_programming.h"
#include "constraint_solver.h"
#include "taylor_set.h"
#include "affine_set.h"
#include "polytope.h"
#include "polyhedron.h"

#include "graphics_interface.h"
#include "hybrid_enclosure.h"
#include "hybrid_set.h"
#include <expression.h>


namespace Ariadne {



template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }





//-------------- HybridEnclosure -----------------------------------------//

HybridEnclosure::~HybridEnclosure() {
}

HybridEnclosure::HybridEnclosure()
    : _location("")
{
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const Box& box)
    : _location(location), _events(), _domain(box), _state(VectorIntervalFunction::identity(box)),
      _time(ScalarIntervalFunction::constant(box,0.0)), _reduced_domain(_domain)
{
}

/*
HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const ContinuousStateSetType& set)
    : _location(location), _events(), _set(set),
      _time(ScalarIntervalFunction::constant(set.domain(),0.0))
{
}
*/

HybridEnclosure::HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>& hpair)
    : _location(hpair.first), _events(), _domain(hpair.second.domain()), _state(_domain,hpair.second.function()),
      _time(ScalarIntervalFunction::constant(_domain,0.0)), _reduced_domain(_domain)
{
    _negative_constraints=hpair.second.negative_constraints();
    _zero_constraints=hpair.second.zero_constraints();
}


HybridEnclosure* HybridEnclosure::clone() const {
    return new HybridEnclosure(*this);
}

List<DiscreteEvent> const&
HybridEnclosure::previous_events() const
{
    return this->_events;
}

VectorIntervalFunction const&
HybridEnclosure::space_function() const
{
    return this->_state;
}

ScalarIntervalFunction const&
HybridEnclosure::time_function() const
{
    return this->_time;
}

IntervalVector
HybridEnclosure::space_bounding_box() const
{
    return this->_state.codomain();
}

Interval
HybridEnclosure::time_range() const
{
    return this->_time.codomain();
}

uint
HybridEnclosure::number_of_constraints() const
{
    return this->_negative_constraints.size()+this->_zero_constraints.size();
}

uint
HybridEnclosure::number_of_parameters() const
{
    return this->_domain.size();
}

const IntervalVector&
HybridEnclosure::parameter_domain() const
{
    return this->_domain;
}

HybridBox
HybridEnclosure::bounding_box() const
{
    return HybridBox(this->_location,this->space_bounding_box());
}


void HybridEnclosure::new_invariant(DiscreteEvent event, ScalarFunction constraint) {
    ScalarIntervalFunction constraint_wrt_params=compose(constraint,this->_state);
    Interval range=constraint_wrt_params.evaluate(this->_domain);
    if(range.upper()>=0.0) {
        //this->_constraint_events.push_back((this->_events,event));
        this->_negative_constraints.append(constraint_wrt_params);
    }
}

void HybridEnclosure::new_invariant(DiscreteEvent event, ScalarIntervalFunction constraint) {
    ScalarIntervalFunction constraint_wrt_params=compose(constraint,this->_state);
    Interval range=constraint_wrt_params.evaluate(this->_domain);
    if(range.upper()>=0.0) {
        //this->_constraint_events.push_back((this->_events,event));
        this->_negative_constraints.append(constraint_wrt_params);
    }
}

void HybridEnclosure::new_activation(DiscreteEvent event, ScalarFunction constraint) {
    //this->_constraint_events.push_back((this->_events,event));
    this->_negative_constraints.append(compose(-constraint,this->_state));
}

void HybridEnclosure::new_guard(DiscreteEvent event, ScalarFunction constraint) {
    //this->_constraint_events.push_back((this->_events,event));
    this->_zero_constraints.append(compose(constraint,this->_state));
}

void HybridEnclosure::new_guard(DiscreteEvent event, ScalarFunction constraint, ScalarIntervalFunction crossing_time) {
    ARIADNE_NOT_IMPLEMENTED;
}




void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, VectorFunction map)
{
    ARIADNE_ASSERT(map.argument_size()==this->dimension());
    this->_events.append(event);
    this->_location=target;
    this->_state=compose(map,this->_state);
}

void HybridEnclosure::_apply_flow(VectorIntervalFunction phi, Float h, ScalarIntervalFunction elps)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    Interval time_domain=Interval(0,h);
    IntervalVector new_domain=join(this->_domain,time_domain);
    this->_domain=new_domain;
    this->_reduced_domain=join(this->_reduced_domain,time_domain);
    ScalarIntervalFunction time_step_function=ScalarIntervalFunction::coordinate(new_domain,new_domain.size()-1u);
    this->_time=embed(this->_time,time_domain)+time_step_function;
    this->_state=unchecked_compose(phi,join(embed(this->_state,time_domain),time_step_function));
    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        this->_negative_constraints[i]=embed(this->_negative_constraints[i],time_domain);
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        this->_zero_constraints[i]=embed(this->_zero_constraints[i],time_domain);
    }
    this->_negative_constraints.append(time_step_function-embed(elps,time_domain));
}

void HybridEnclosure::apply_flow(VectorIntervalFunction phi, Float h)
{
    this->_apply_flow(phi,h,ScalarIntervalFunction::constant(this->_domain,h));
}

void HybridEnclosure::apply_flow(VectorIntervalFunction phi, ScalarIntervalFunction eps)
{
    Float h=phi.domain()[phi.domain().size()-1].upper();
    this->_apply_flow(phi,h,unchecked_compose(eps,this->_state));
}

void HybridEnclosure::apply_flow_and_bound_time(VectorIntervalFunction phi, ScalarIntervalFunction omega)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    ARIADNE_ASSERT_MSG(omega.argument_size()==this->number_of_parameters(),
                       "Final time "<<omega<<" must be a function of "<<this->number_of_parameters()<<" parameterisation variables");
    Float h=phi.domain()[phi.domain().size()-1].upper();
    this->_apply_flow(phi,h,omega-this->_time);
    this->bound_time(embed(omega,Interval(0,h)));
}

void HybridEnclosure::apply_flow_and_bound_time(VectorIntervalFunction phi, Float tmax)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    Float h=phi.domain()[phi.domain().size()-1].upper();
    this->_apply_flow(phi,h,tmax-this->_time);
    this->bound_time(ScalarIntervalFunction::constant(this->_domain,tmax));
}

void HybridEnclosure::_apply_flow_step(VectorIntervalFunction phi, ScalarIntervalFunction elps)
{
    // xi'(s) = phi(xi(s),eps(s)) where range(eps) in [0,h]
    // tau'(s) = tau(s)+eps(s)
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    this->_time=this->_time+elps;
    this->_state=unchecked_compose(phi,join(this->_state,elps));
}

void HybridEnclosure::apply_flow_step(VectorIntervalFunction phi, Float h)
{
    // xi'(s) = phi(xi(s),eps(s)) where range(eps) in [0,h]
    // tau'(s) = tau(s)+eps(s)
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    this->_apply_flow_step(phi,ScalarIntervalFunction::constant(this->_domain,h));
}

void HybridEnclosure::apply_flow_step(VectorIntervalFunction phi, ScalarIntervalFunction eps)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    ARIADNE_ASSERT_MSG(eps.argument_size()==this->dimension(),
                       "Evolved time "<<eps<<" must be a function of "<<this->dimension()<<" variables");
    this->_apply_flow_step(phi,unchecked_compose(eps,this->_state));
}

void HybridEnclosure::apply_flow_and_set_time(VectorIntervalFunction phi, ScalarIntervalFunction omega)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    ARIADNE_ASSERT_MSG(omega.argument_size()==this->number_of_parameters(),
                       "Final time "<<omega<<" must be a function of "<<this->number_of_parameters()<<" parameterisation variables");
    this->_apply_flow_step(phi,omega-this->_time);
    this->_time=omega;
}

void HybridEnclosure::apply_flow_and_set_time(VectorIntervalFunction phi, Float tmax)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    this->_apply_flow_step(phi,tmax-this->_time);
    this->_time=tmax;
}


void HybridEnclosure::bound_time(Float tmax) {
    this->_negative_constraints.append(this->_time-tmax);
}

void HybridEnclosure::bound_time(ScalarFunction tmax) {
    this->_negative_constraints.append(this->_time-ScalarIntervalFunction(this->_domain,tmax));
}

void HybridEnclosure::bound_time(ScalarIntervalFunction tmax) {
    this->_negative_constraints.append(this->_time-tmax);
}

void HybridEnclosure::set_time(Float time)
{
    this->_zero_constraints.append(this->_time-time);
}

void HybridEnclosure::set_time(ScalarFunction time)
{
    this->_zero_constraints.append(this->_time-ScalarIntervalFunction(this->_domain,time));
}


void HybridEnclosure::set_maximum_time(DiscreteEvent event, Float final_time)
{
    this->_negative_constraints.append(this->_time-final_time); // Deprecated
}

void HybridEnclosure::new_time_step_bound(DiscreteEvent event, ScalarIntervalFunction constraint) {
    ARIADNE_NOT_IMPLEMENTED; // Deprecated
}

void HybridEnclosure::set_step_time(Float time)
{
    ARIADNE_NOT_IMPLEMENTED; // Deprecated
}



const DiscreteLocation& HybridEnclosure::location() const {
    return this->_location;
}

/*
ConstrainedImageSet HybridEnclosure::continuous_state_set() const {
    return ConstrainedImageSet(this->_domain,this->_state,this->_constraints);
}
*/


TaylorConstrainedImageSet HybridEnclosure::continuous_state_set() const {
    TaylorConstrainedImageSet set(this->_state);
    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        set.new_negative_constraint(this->_negative_constraints[i]);
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        set.new_zero_constraint(this->_zero_constraints[i]);
    }
    return set;
}


uint HybridEnclosure::dimension() const {
    return this->_state.result_size();
}

tribool HybridEnclosure::empty() const {
    List<NonlinearConstraint> constraints=this->constraints();
    if(constraints.empty()) { return this->parameter_domain().empty(); }
    for(uint i=0; i!=constraints.size(); ++i) {
        if(Ariadne::disjoint(constraints[i].function().evaluate(this->_reduced_domain),constraints[i].bounds())) {
            return true;
        }
    }
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(constraints,this->_reduced_domain);
    if(this->_reduced_domain.empty()) { return true; }
    else { return indeterminate; }
}

tribool HybridEnclosure::subset(const HybridBox& hbx) const {
    const Box& bx = hbx.continuous_state_set();
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"HybridEnclosure::subset(HybridBox): self="<<*this<<", box="<<hbx);
    List<NonlinearConstraint> constraints=this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(constraints,this->_reduced_domain);

    if(_reduced_domain.empty()) { return true; }
    if(this->_location!=hbx.location()) { return indeterminate; }

    for(uint i=0; i!=bx.dimension(); ++i) {
        const Box test_domain=this->_reduced_domain;
        constraints.append(ScalarIntervalFunction(this->_state[i]).function() <= bx[i].lower());
        if(possibly(contractor.feasible(constraints,test_domain).first)) { return indeterminate; }
        constraints.back()=(ScalarIntervalFunction(this->_state[i]).function() >= bx[i].upper());
        if(possibly(contractor.feasible(constraints,test_domain).first)) { return indeterminate; }
        constraints.pop_back();
    }
    return true;
}

tribool HybridEnclosure::disjoint(const HybridBox& hbx) const {
    const Box& bx = hbx.continuous_state_set();
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"HybridEnclosure::subset(HybridBox): self="<<*this<<", box="<<hbx);
    List<NonlinearConstraint> constraints=this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(constraints,this->_reduced_domain);

    if(_reduced_domain.empty()) { return true; }
    if(this->_location!=hbx.location()) { return indeterminate; }

    const Box test_domain=this->_reduced_domain;
    for(uint i=0; i!=bx.dimension(); ++i) {
        constraints.append(ScalarIntervalFunction(this->_state[i]).function() >= bx[i].lower());
        constraints.append(ScalarIntervalFunction(this->_state[i]).function() <= bx[i].upper());
    }
    return !contractor.feasible(constraints,test_domain).first;
}

tribool HybridEnclosure::satisfies(NonlinearConstraint c) const
{
    return this->continuous_state_set().satisfies(c);
}


List<NonlinearConstraint> HybridEnclosure::constraints() const {
    List<NonlinearConstraint> result;
    for(List<ScalarIntervalFunction>::const_iterator iter=_negative_constraints.begin(); iter!=_negative_constraints.end(); ++iter) {
        result.append(iter->function()<=0);
    }
    for(List<ScalarIntervalFunction>::const_iterator iter=_zero_constraints.begin(); iter!=_zero_constraints.end(); ++iter) {
        result.append(iter->function()==0);
    }
    return result;
}


void HybridEnclosure::draw(CanvasInterface& canvas) const
{
    this->continuous_state_set().draw(canvas);
}


std::ostream& HybridEnclosure::write(std::ostream& os) const
{
    return os << "HybridEnclosure"
              << "( events=" << this->_events
              << ", range=" << this->_state(this->_domain)
              << ", location=" << this->_location
              << ", domain=" << this->_domain
              << ", range=" << this->_state(this->_domain)
              << ", subdomain=" << this->_reduced_domain
              << ", empty=" << this->empty()
              << ", state=" << polynomial(this->_state)
              << ", negative=" << polynomials(this->_negative_constraints)
              << ", zero=" << polynomials(this->_zero_constraints)
              << ", time="<<polynomial(this->_time) << ")";
}









} // namespace Ariadne
