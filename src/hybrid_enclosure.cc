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
    : _location(""), _events(), _set(), _time(ScalarIntervalFunction::constant(_set._domain,0.0))
{
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const Box& box)
    : _location(location), _events(), _set(box), _time(ScalarIntervalFunction::constant(_set._domain,0.0))
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
    : _location(hpair.first), _events(), _set(hpair.second), _time(ScalarIntervalFunction::constant(_set._domain,0.0))
{
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
    return this->_set._function;
}

ScalarIntervalFunction const&
HybridEnclosure::time_function() const
{
    return this->_time;
}

IntervalVector
HybridEnclosure::space_bounding_box() const
{
    return this->_set._function.codomain();
}

Interval
HybridEnclosure::time_range() const
{
    return this->_time.codomain();
}

uint
HybridEnclosure::number_of_constraints() const
{
    return this->_set.number_of_constraints();
}

uint
HybridEnclosure::number_of_parameters() const
{
    return this->_set.number_of_parameters();
}

const IntervalVector&
HybridEnclosure::parameter_domain() const
{
    return this->_set._domain;
}

HybridBox
HybridEnclosure::bounding_box() const
{
    return HybridBox(this->_location,this->space_bounding_box());
}


void HybridEnclosure::new_invariant(DiscreteEvent event, ScalarFunction constraint) {
    ScalarIntervalFunction constraint_wrt_params=compose(constraint,this->_set._function);
    Interval range=constraint_wrt_params.evaluate(this->_set._domain);
    if(range.upper()>=0.0) {
        //this->_constraint_events.push_back((this->_events,event));
        this->_set._constraints.append(constraint_wrt_params);
    }
}

void HybridEnclosure::new_invariant(DiscreteEvent event, ScalarIntervalFunction constraint) {
    ScalarIntervalFunction constraint_wrt_params=compose(constraint,this->_set._function);
    Interval range=constraint_wrt_params.evaluate(this->_set._domain);
    if(range.upper()>=0.0) {
        //this->_constraint_events.push_back((this->_events,event));
        this->_set._constraints.append(constraint_wrt_params);
    }
}

void HybridEnclosure::new_activation(DiscreteEvent event, ScalarFunction constraint) {
    //this->_constraint_events.push_back((this->_events,event));
    this->_set._constraints.append(compose(-constraint,this->_set._function));
}

void HybridEnclosure::new_guard(DiscreteEvent event, ScalarFunction constraint) {
    //this->_constraint_events.push_back((this->_events,event));
    this->_set._equations.append(compose(constraint,this->_set._function));
}

void HybridEnclosure::new_guard(DiscreteEvent event, ScalarFunction constraint, ScalarIntervalFunction crossing_time) {
    ARIADNE_NOT_IMPLEMENTED;
}

void HybridEnclosure::new_parameter_constraint(DiscreteEvent event, NonlinearConstraint constraint) {
    ARIADNE_ASSERT_MSG(constraint.function().argument_size()==parameter_domain().size(),
                       "constraint "<<constraint<<" is incompatible with parameter domain "<<parameter_domain());
    ScalarIntervalFunction function(this->_set._domain,constraint.function());
    const Interval bounds=constraint.bounds();
    if(bounds.singleton()) {
        this->_set._equations.append(function-bounds.upper());
    } else {
        if(bounds.upper()!=+inf<Float>()) {
            this->_set._constraints.append(function-bounds.upper());
        }
        if(bounds.lower()!=-inf<Float>()) {
            this->_set._constraints.append(bounds.lower()-function);
        }
    }
}

void HybridEnclosure::new_state_constraint(DiscreteEvent event, NonlinearConstraint constraint) {
    ARIADNE_ASSERT_MSG(constraint.function().argument_size()==dimension(),
                       "constraint "<<constraint<<" is incompatible with dimension "<<dimension());
    NonlinearConstraint parameter_constraint(compose(constraint.function(),this->_set._function).function(),constraint.bounds());
    this->new_parameter_constraint(event,parameter_constraint);
}




void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, VectorFunction map)
{
    ARIADNE_ASSERT(map.argument_size()==this->dimension());
    this->_events.append(event);
    this->_location=target;
    this->_set._function=compose(map,this->_set._function);
}

void HybridEnclosure::_apply_flow(VectorIntervalFunction phi, Float h, ScalarIntervalFunction elps)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    Interval time_domain=Interval(0,h);
    IntervalVector new_domain=join(this->_set._domain,time_domain);
    this->_set._domain=new_domain;
    this->_set._reduced_domain=join(this->_set._reduced_domain,time_domain);
    ScalarIntervalFunction time_step_function=ScalarIntervalFunction::coordinate(new_domain,new_domain.size()-1u);
    this->_time=embed(this->_time,time_domain)+time_step_function;
    this->_set._function=unchecked_compose(phi,join(embed(this->_set._function,time_domain),time_step_function));
    for(uint i=0; i!=this->_set._constraints.size(); ++i) {
        this->_set._constraints[i]=embed(this->_set._constraints[i],time_domain);
    }
    for(uint i=0; i!=this->_set._equations.size(); ++i) {
        this->_set._equations[i]=embed(this->_set._equations[i],time_domain);
    }
    this->_set._constraints.append(time_step_function-embed(elps,time_domain));
}

void HybridEnclosure::_apply_flow(VectorIntervalFunction phi, Float h)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    Interval time_domain=Interval(0,h);
    IntervalVector new_domain=join(this->_set._domain,time_domain);
    this->_set._domain=new_domain;
    this->_set._reduced_domain=join(this->_set._reduced_domain,time_domain);
    ScalarIntervalFunction time_step_function=ScalarIntervalFunction::coordinate(new_domain,new_domain.size()-1u);
    this->_time=embed(this->_time,time_domain)+time_step_function;
    this->_set._function=unchecked_compose(phi,join(embed(this->_set._function,time_domain),time_step_function));
    for(uint i=0; i!=this->_set._constraints.size(); ++i) {
        this->_set._constraints[i]=embed(this->_set._constraints[i],time_domain);
    }
    for(uint i=0; i!=this->_set._equations.size(); ++i) {
        this->_set._equations[i]=embed(this->_set._equations[i],time_domain);
    }
}

void HybridEnclosure::apply_flow_for(VectorIntervalFunction phi, Float h)
{
    this->_apply_flow(phi,h,ScalarIntervalFunction::constant(this->_set._domain,h));
}

void HybridEnclosure::apply_flow_for(VectorIntervalFunction phi, ScalarIntervalFunction eps)
{
    Float h=phi.domain()[phi.domain().size()-1].upper();
    // Pre-compute the evolved time in the new domain
    ScalarIntervalFunction evolve_time=embed(compose(eps,this->_set._function),Interval(0,h));
    ScalarIntervalFunction time_step_function=embed(this->_set._domain,ScalarIntervalFunction::identity(Interval(0,h)));
    this->_apply_flow(phi,h,unchecked_compose(eps,this->_set._function));
}

void HybridEnclosure::apply_flow_to(VectorIntervalFunction phi, ScalarIntervalFunction omega)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    ARIADNE_ASSERT_MSG(omega.argument_size()==this->number_of_parameters(),
                       "Final time "<<omega<<" must be a function of "<<this->number_of_parameters()<<" parameterisation variables");
    Float h=phi.domain()[phi.domain().size()-1].upper();
    this->_apply_flow(phi,h);
    this->bound_time(embed(omega,Interval(0,h)));
}

void HybridEnclosure::apply_flow_to(VectorIntervalFunction phi, Float tmax)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    Float h=phi.domain()[phi.domain().size()-1].upper();
    this->_apply_flow(phi,h);
    this->bound_time(ScalarIntervalFunction::constant(this->_set._domain,tmax));
}

void HybridEnclosure::_apply_flow_step(VectorIntervalFunction phi, ScalarIntervalFunction elps)
{
    // xi'(s) = phi(xi(s),eps(s)) where range(eps) in [0,h]
    // tau'(s) = tau(s)+eps(s)
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    this->_time=this->_time+elps;
    this->_set._function=unchecked_compose(phi,join(this->_set._function,elps));
}

void HybridEnclosure::apply_flow_step_for(VectorIntervalFunction phi, Float h)
{
    // xi'(s) = phi(xi(s),eps(s)) where range(eps) in [0,h]
    // tau'(s) = tau(s)+eps(s)
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    this->_apply_flow_step(phi,ScalarIntervalFunction::constant(this->_set._domain,h));
}

void HybridEnclosure::apply_flow_step_for(VectorIntervalFunction phi, ScalarIntervalFunction eps)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    ARIADNE_ASSERT_MSG(eps.argument_size()==this->dimension(),
                       "Evolved time "<<eps<<" must be a function of "<<this->dimension()<<" variables");
    this->_apply_flow_step(phi,unchecked_compose(eps,this->_set._function));
}

void HybridEnclosure::apply_flow_step_to(VectorIntervalFunction phi, ScalarIntervalFunction omega)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    ARIADNE_ASSERT_MSG(omega.argument_size()==this->number_of_parameters(),
                       "Final time "<<omega<<" must be a function of "<<this->number_of_parameters()<<" parameterisation variables");
    this->_set._function=unchecked_compose(phi,join(this->_set._function,omega-this->_time));
    this->_time=omega;
}

void HybridEnclosure::apply_flow_step_to(VectorIntervalFunction phi, Float tmax)
{
    ARIADNE_ASSERT_MSG(phi.argument_size()==this->dimension()+1,
                       "Flow "<<phi<<" must be a function of "<<this->dimension()<<"+1 variables");
    this->_set._function=unchecked_compose(phi,join(this->_set._function,tmax-this->_time));
    this->_time=tmax;
}


void HybridEnclosure::bound_time(Float tmax) {
    this->_set._constraints.append(this->_time-tmax);
}

void HybridEnclosure::bound_time(ScalarFunction tmax) {
    this->_set._constraints.append(this->_time-ScalarIntervalFunction(this->_set._domain,tmax));
}

void HybridEnclosure::bound_time(ScalarIntervalFunction tmax) {
    this->_set._constraints.append(this->_time-tmax);
}

void HybridEnclosure::set_time(Float time)
{
    this->_set._equations.append(this->_time-time);
}

void HybridEnclosure::set_time(ScalarFunction time)
{
    this->_set._equations.append(this->_time-ScalarIntervalFunction(this->_set._domain,time));
}


void HybridEnclosure::set_maximum_time(DiscreteEvent event, Float final_time)
{
    this->_set._constraints.append(this->_time-final_time); // Deprecated
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
    return ConstrainedImageSet(this->_set._domain,this->_set._function,this->_constraints);
}
*/


const TaylorConstrainedImageSet&
HybridEnclosure::continuous_state_set() const {
    return this->_set;
}


uint HybridEnclosure::dimension() const {
    return this->_set._function.result_size();
}

tribool HybridEnclosure::empty() const {
    return this->_set.empty();
}

tribool HybridEnclosure::subset(const HybridBox& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_state_set().subset(hbx.continuous_state_set()); }
    else { return this->continuous_state_set().empty(); }
}

tribool HybridEnclosure::disjoint(const HybridBox& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_state_set().disjoint(hbx.continuous_state_set()); }
    else { return true; }
}

tribool HybridEnclosure::satisfies(NonlinearConstraint c) const
{
    return this->continuous_state_set().satisfies(c);
}


List<NonlinearConstraint> HybridEnclosure::constraints() const {
    return this->continuous_state_set().constraints();
}


void HybridEnclosure::draw(CanvasInterface& canvas) const
{
    this->continuous_state_set().draw(canvas);
}


std::ostream& HybridEnclosure::write(std::ostream& os) const
{
    return os << "HybridEnclosure"
              << "( events=" << this->_events
              << ", range=" << this->_set._function(this->_set._domain)
              << ", location=" << this->_location
              << ", domain=" << this->_set._domain
              << ", range=" << this->_set._function(this->_set._domain)
              << ", subdomain=" << this->_set._reduced_domain
              << ", empty=" << this->empty()
              << ", state=" << polynomial(this->_set._function)
              << ", negative=" << polynomials(this->_set._constraints)
              << ", zero=" << polynomials(this->_set._equations)
              << ", time="<<polynomial(this->_time) << ")";
}









} // namespace Ariadne
