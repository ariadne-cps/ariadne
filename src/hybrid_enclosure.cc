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

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;






//-------------- HybridEnclosure -----------------------------------------//

HybridEnclosure::~HybridEnclosure() {
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const Box& box)
    : _location(location), _events(), _set(box),
      _time(ScalarTaylorFunction::constant(box,0.0))
{
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const ContinuousStateSetType& set)
    : _location(location), _events(), _set(set),
      _time(ScalarTaylorFunction::constant(set.domain(),0.0))
{
}

HybridEnclosure::HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>& hpair)
    : _location(hpair.first), _events(), _set(hpair.second),
      _time(ScalarTaylorFunction::constant(hpair.second.domain(),0.0))
{
}

HybridEnclosure* HybridEnclosure::clone() const {
    return new HybridEnclosure(*this);
}

VectorTaylorFunction const&
HybridEnclosure::space_function() const
{
    return this->_set.taylor_function();
}

ScalarTaylorFunction const&
HybridEnclosure::time_function() const
{
    return this->_time;
}


void HybridEnclosure::new_invariant(DiscreteEvent event, ScalarFunction constraint) {
    Interval range=compose(constraint,this->space_function()).evaluate(this->continuous_state_set().domain());
    if(range.upper()>=0.0) {
        this->_constraint_events.push_back((this->_events,event));
        this->_set.new_negative_constraint(compose(constraint,this->_set.function()));
    }
}

void HybridEnclosure::new_activation(DiscreteEvent event, ScalarFunction constraint) {
    this->_constraint_events.push_back((this->_events,event));
    this->_set.new_negative_constraint(compose(-constraint,this->_set.function()));
}

void HybridEnclosure::new_guard(DiscreteEvent event, ScalarFunction constraint) {
    this->_constraint_events.push_back((this->_events,event));
    this->_set.new_negative_constraint(compose(constraint,this->_set.function()));
    this->_constraint_events.push_back((this->_events,event));
    this->_set.new_negative_constraint(compose(-constraint,this->_set.function()));
}

void HybridEnclosure::new_guard(DiscreteEvent event, ScalarFunction constraint, ScalarTaylorFunction crossing_time) {
    // Remove the invariant corresponding to the event
    ARIADNE_ASSERT(crossing_time.argument_size()+1u==this->continuous_state_set ().number_of_parameters());
    List<DiscreteEvent> events=(this->_events,event);
    for(uint i=this->_constraint_events.size()-1u; i!=uint(-1); --i) {
        if(this->_constraint_events[i]==events) {
            this->_constraint_events[i]=this->_constraint_events.back();
            this->_set._constraints[i]=this->_set._constraints.back();
            this->_constraint_events.pop_back();
            this->_set._constraints.pop_back();
        }
    }
    this->set_dwell_time(event,crossing_time);
}

void HybridEnclosure::new_time_bound(DiscreteEvent event, ScalarFunction constraint) {
    this->_set.new_negative_constraint(constraint);
}


void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, VectorFunction map)
{
    this->_events.append(event);
    this->_location=target;
    this->_set.apply_map(map);
}

void HybridEnclosure::apply_flow(VectorTaylorFunction phi, Float time)
{
    this->_set.apply_map(partial_evaluate(phi,phi.argument_size()-1,time));
    this->_time+=time;
}

void HybridEnclosure::apply_flow(VectorFunction phi, Interval time_domain)
{
    this->_set.apply_flow(phi,time_domain);
    Vector<Interval> const& new_domain=this->_set.domain();
    this->_time=embed(this->_time,time_domain)+ScalarTaylorFunction::coordinate(new_domain,new_domain.size()-1u);
}

void HybridEnclosure::apply_flow(VectorTaylorFunction phi, Interval time_domain)
{
    this->_set.apply_flow(phi,time_domain);
    Vector<Interval> const& new_domain=this->_set.domain();
    this->_time=embed(this->_time,time_domain)+ScalarTaylorFunction::coordinate(new_domain,new_domain.size()-1u);
}



void HybridEnclosure::set_maximum_time(DiscreteEvent event, Float final_time)
{
    this->_set.new_negative_constraint(this->_time-final_time);
}

void HybridEnclosure::set_step_time(Float step_time)
{
    uint n=this->continuous_state_set().number_of_parameters()-1u;
    this->_set.substitute(n,step_time);
    this->_time=partial_evaluate(this->_time,n,step_time);
}

void HybridEnclosure::set_time(DiscreteEvent event, Float final_time)
{
    ScalarTaylorFunction dwell_time=implicit(this->_time-final_time);
    this->set_dwell_time(event,dwell_time);
    this->_time=ScalarTaylorFunction::constant(this->_set.domain(),final_time);
}

void HybridEnclosure::set_dwell_time(DiscreteEvent event, Float time_step)
{
    const uint m=this->_time.argument_size()-1;
    IntervalVector new_domain=project(this->_time.domain(),range(0,m));
    ScalarTaylorFunction dwell_time=ScalarTaylorFunction::constant(new_domain,time_step);
    this->set_dwell_time(event,dwell_time);
}

void HybridEnclosure::set_dwell_time(DiscreteEvent event, ScalarFunction time)
{
    ARIADNE_ASSERT_MSG(time.argument_size()+1==this->_set.number_of_parameters(),*this<<" "<<time);
    const uint n=time.argument_size();
    IntervalVector new_domain=project(this->_time.domain(),range(0,n));
    ScalarTaylorFunction time_model(new_domain,time);
    this->set_dwell_time(event,time_model);
}

void HybridEnclosure::set_dwell_time(DiscreteEvent event, ScalarTaylorFunction dwell_time)
{
    const uint n=dwell_time.argument_size();
    IntervalVector old_domain=this->_set.domain();
    IntervalVector new_domain=project(old_domain,range(0,n));
    ARIADNE_ASSERT_MSG(n+1==this->_set.number_of_parameters(),*this<<" "<<dwell_time);
    ARIADNE_ASSERT(Box(new_domain)==Box(dwell_time.domain()));

    this->_set.substitute(n,dwell_time);
    this->_time=substitute(this->_time,n,dwell_time);
}

const DiscreteLocation& HybridEnclosure::location() const {
    return this->_location;
}

const TaylorConstrainedImageSet& HybridEnclosure::continuous_state_set() const {
    return this->_set;
}

uint HybridEnclosure::dimension() const {
    return this->_set.dimension();
}

tribool HybridEnclosure::empty() const {
    return this->_set.empty();
}


void HybridEnclosure::draw(CanvasInterface& canvas) const
{
    this->_set.draw(canvas);
}

List< Polynomial<Interval> > polynomials(const List<ScalarTaylorFunction>& tf) {
    List< Polynomial<Interval> > p;
    for(uint i=0; i!=tf.size(); ++i) {
        p.append(tf[i].polynomial());
    }
    return p;
}

std::ostream& HybridEnclosure::write(std::ostream& os) const
{
    return os << "HybridEnclosure"
              << "( events=" << this->_events
              << ", location=" << this->_location
              << ", domain=" << this->_set.domain()
              << ", range=" << this->_set.bounding_box()
              << ", function=" << this->_set.function()
              << ", negative=" << polynomials(this->_set.negative_constraints())
              << ", zero=" << polynomials(this->_set.zero_constraints())
              << ", time="<<this->_time.polynomial() << ")";
}









} // namespace Ariadne
