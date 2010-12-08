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
#include "polynomial.h"
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

#include "graphics_interface.h"
#include "hybrid_enclosure.h"
#include "hybrid_set.h"


namespace Ariadne {



template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }





//-------------- HybridEnclosure -----------------------------------------//

HybridEnclosure::~HybridEnclosure() {
}

HybridEnclosure::HybridEnclosure()
    : _location(""), _events(), _set(), _time(ScalarIntervalFunction::constant(_set._domain,0.0)), _dwell_time(_time)
{
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const Box& box)
    : _location(location), _events(), _set(box), _time(ScalarIntervalFunction::constant(_set._domain,0.0)), _dwell_time(_time)
{
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const ContinuousStateSetType& set, const ScalarTaylorFunction& time)
    : _location(location), _events(), _set(set),
      _time(Ariadne::restrict(time,set.domain())), _dwell_time(set.domain())
{
}

HybridEnclosure::HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>& hpair)
    : _location(hpair.first), _events(), _set(hpair.second), _time(ScalarIntervalFunction::constant(_set._domain,0.0)), _dwell_time(_time)
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

ScalarIntervalFunction const&
HybridEnclosure::dwell_time_function() const
{
    return this->_dwell_time;
}

void HybridEnclosure::set_time_function(const ScalarIntervalFunction& time_function)
{
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->parameter_domain(),time_function.domain()),
                       "Domain of "<<time_function<<" does not contain parameter domain "<<this->parameter_domain());
    ARIADNE_ASSERT_MSG(this->parameter_domain()==time_function.domain(),
                       "Domain of "<<time_function<<" does not equal parameter domain "<<this->parameter_domain());
    this->_time=time_function;
}

IntervalVector
HybridEnclosure::space_bounding_box() const
{
    ARIADNE_LOG(8,"space_codomain="<<this->_set._function.codomain()<<" space_range="<<this->_set._function(this->_set._reduced_domain)<<"\n");
    //return this->_set._function(this->_set._reduced_domain);
    return this->_set.bounding_box();
}

Interval
HybridEnclosure::time_range() const
{
    ARIADNE_LOG(8,"time_codomain="<<this->_time.codomain()<<" time_range="<<this->_time(this->_set._reduced_domain)<<"\n");
    return this->_time.codomain();
    //return this->_time(this->_set._reduced_domain);
}

Interval
HybridEnclosure::dwell_time_range() const
{
    ARIADNE_LOG(8,"dwell_time_codomain="<<this->_dwell_time.codomain()<<" dwell_time_range="<<this->_dwell_time(this->_set._reduced_domain)<<"\n");
    return this->_dwell_time.codomain();
    //return this->_time(this->_set._reduced_domain);
}

uint
HybridEnclosure::number_of_constraints() const
{
    return this->_set.number_of_constraints();
}

uint
HybridEnclosure::number_of_inequality_constraints() const
{
    return this->_set._constraints.size();
}

uint
HybridEnclosure::number_of_equality_constraints() const
{
    return this->_set._equations.size();
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


void HybridEnclosure::new_parameter(Interval ivl)
{
    IntervalVector new_domain=join(this->_set._domain,ivl);
    this->_set._domain=new_domain;
    this->_set._reduced_domain=join(this->_set._reduced_domain,ivl);
    this->_time=embed(this->_time,ivl);
    this->_dwell_time=embed(this->_dwell_time,ivl);
    this->_set._function=embed(this->_set._function,ivl);
    for(uint i=0; i!=this->_set._constraints.size(); ++i) {
        this->_set._constraints[i]=embed(this->_set._constraints[i],ivl);
    }
    for(uint i=0; i!=this->_set._equations.size(); ++i) {
        this->_set._equations[i]=embed(this->_set._equations[i],ivl);
    }
}

void HybridEnclosure::new_invariant(DiscreteEvent event, RealScalarFunction constraint) {
    ScalarIntervalFunction constraint_wrt_params=compose(constraint,this->_set._function);
    Interval range=constraint_wrt_params.evaluate(this->_set._domain);
    if(range.upper()>=0.0) {
        //this->_constraint_events.push_back((this->_events,event));
        this->_set._constraints.append(constraint_wrt_params);
    }
}

void HybridEnclosure::new_invariant(DiscreteEvent event, ScalarIntervalFunction constraint) {
    ScalarIntervalFunction constraint_wrt_params=unchecked_compose(constraint,this->_set._function);
    Interval range=constraint_wrt_params.evaluate(this->_set._domain);
    if(range.upper()>=0.0) {
        //this->_constraint_events.push_back((this->_events,event));
        this->_set._constraints.append(constraint_wrt_params);
    }
}

void HybridEnclosure::new_activation(DiscreteEvent event, RealScalarFunction constraint) {
    //this->_constraint_events.push_back((this->_events,event));
    this->_set._constraints.append(compose(-constraint,this->_set._function));
}

void HybridEnclosure::new_guard(DiscreteEvent event, RealScalarFunction constraint) {
    //this->_constraint_events.push_back((this->_events,event));
    this->_set._equations.append(compose(constraint,this->_set._function));
}

void HybridEnclosure::new_guard(DiscreteEvent event, RealScalarFunction constraint, ScalarIntervalFunction crossing_time) {
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
    NonlinearConstraint parameter_constraint(compose(constraint.function(),this->_set._function).real_function(),constraint.bounds());
    this->new_parameter_constraint(event,parameter_constraint);
}




void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, RealVectorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"*this="<<*this<<", event="<<event<<", target="<<target<<", map="<<map);
    this->_events.append(event);
    this->_location=target;
    this->_set._function=compose(map,this->_set._function);
    this->_dwell_time=0.0;
}

void HybridEnclosure::apply_evolve_step(const VectorIntervalFunction& phi, const ScalarIntervalFunction& elps)
{
    // xi'(s) = phi(xi(s),eps(s)) where range(eps) in [0,h]
    // tau'(s) = tau(s)+eps(s)
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    this->_time=this->_time+elps;
    this->_dwell_time=this->_dwell_time+elps;
    this->_set._function=unchecked_compose(phi,join(this->_set._function,elps));
}

void HybridEnclosure::apply_finishing_evolve_step(const VectorIntervalFunction& phi, const ScalarIntervalFunction& omega)
{
    // xi'(s) = phi(xi(s),omega(s)-tau(s)) where range(omega-tau) in [0,h]
    // tau'(s) = tau(s)+eps(s)
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(omega.argument_size()==this->number_of_parameters());
    this->_set._function=unchecked_compose(phi,join(this->_set._function,omega-this->_time));
    this->_dwell_time=this->_dwell_time+(omega-this->_time);
    this->_time=omega;
}


void HybridEnclosure::apply_reach_step(const VectorIntervalFunction& phi, const ScalarIntervalFunction& elps)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    Interval time_domain=phi.domain()[phi.domain().size()-1];
    this->new_parameter(time_domain);
    const IntervalVector& new_domain=this->parameter_domain();
    ScalarIntervalFunction time_step_function=ScalarIntervalFunction::coordinate(new_domain,new_domain.size()-1u);
    this->_time=this->_time+time_step_function;
    this->_dwell_time=this->_dwell_time+time_step_function;
    this->_set._function=unchecked_compose(phi,join(this->_set._function,time_step_function));
    if(elps.range().lower()<time_domain.upper()) {
        this->_set._constraints.append(time_step_function-embed(elps,time_domain));
    }
}

void HybridEnclosure::apply_full_reach_step(const VectorIntervalFunction& phi)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    Interval time_domain=phi.domain()[phi.domain().size()-1];
    this->new_parameter(time_domain);
    const IntervalVector& new_domain=this->parameter_domain();
    ScalarIntervalFunction time_step_function=ScalarIntervalFunction::coordinate(new_domain,new_domain.size()-1u);
    this->_time=this->_time+time_step_function;
    this->_dwell_time=this->_dwell_time+time_step_function;
    this->_set._function=unchecked_compose(phi,join(this->_set._function,time_step_function));
}



void HybridEnclosure::bound_time(Real tmax) {
    if(this->time_range().upper()>Interval(tmax).lower()) {
        this->_set._constraints.append(this->_time-Interval(tmax));
    }
}

void HybridEnclosure::bound_time(RealScalarFunction tmax) {
    this->_set._constraints.append(this->_time-ScalarIntervalFunction(this->_set._domain,tmax));
}

void HybridEnclosure::bound_time(ScalarIntervalFunction tmax) {
    this->_set._constraints.append(this->_time-tmax);
}

void HybridEnclosure::set_time(Real time)
{
    this->_set._equations.append(this->_time-Interval(time));
}

void HybridEnclosure::set_time(RealScalarFunction time)
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

void
HybridEnclosure::restrict(const IntervalVector& subdomain)
{
    this->_set.restrict(subdomain);
    this->_time.restrict(this->parameter_domain());
    this->_dwell_time.restrict(this->parameter_domain());
}

List<HybridEnclosure>
HybridEnclosure::split() const
{
    List<IntervalVector> subdomains = this->continuous_state_set().splitting_subdomains_zeroth_order();
    List<HybridEnclosure> result(subdomains.size(),*this);
    for(uint i=0; i!=result.size(); ++i) {
        result[i].restrict(subdomains[i]);
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


HybridGridTreeSet outer_approximation(const ListSet<HybridEnclosure>& hls, const HybridGrid& g, uint d) {
    HybridGridTreeSet result(g);
    for(ListSet<HybridEnclosure>::const_iterator iter=hls.begin(); iter!=hls.end(); ++iter) {
        result[iter->first].adjoin_outer_approximation(iter->second,d);
    }
    return result;
}





} // namespace Ariadne
