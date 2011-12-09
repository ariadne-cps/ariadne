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
#include "config.h"

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
#include "affine_set.h"

#include "graphics_interface.h"
#include "hybrid_enclosure.h"
#include "hybrid_set.h"


namespace Ariadne {


std::ostream& operator<<(std::ostream& os, const EnclosureVariableType& evt) {
    switch (evt) {
        case INITIAL: return os << "x";
        case TEMPORAL: return os << "t";
        //case CROSSING: return os << "t";
        //case STEP: return os << "h";
        case PARAMETER: return os << "a";
        case INPUT: return os << "u";
        case NOISE: return os << "v";
        case ERROR: return os << "e";
        case UNKNOWN: default: return os << "s";
    }
}

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

List<String> variable_names(const List<EnclosureVariableType>& vt) {
    std::map<EnclosureVariableType,uint> counts;
    List<String> result;
    for(uint i=0; i!=vt.size(); ++i) {
        result.append( str(vt[i]) + str(counts[vt[i]]++) );
    }
    return result;
}

struct Variables2d {
    RealVariable _x,_y;
    Variables2d(const RealVariable& x, const RealVariable& y) : _x(x), _y(y) { }
    RealVariable const& x_variable() const { return this->_x; };
    RealVariable const& y_variable() const { return this->_y; };
};


//-------------- HybridEnclosure -----------------------------------------//

HybridEnclosure::~HybridEnclosure() {
}

HybridEnclosure::HybridEnclosure()
    : _location(""), _events(), _space(), _set(), _variables()
{
}

HybridEnclosure::HybridEnclosure(const HybridRealExpressionBoundedConstraintSet& hybrid_set,
                                 const RealSpace& space,
                                 const IntervalFunctionModelFactoryInterface& factory)
    : _location(hybrid_set.location()), _events(), _space(space.variable_names()), _set(),
      _variables(space.dimension(),INITIAL)
{
    BoundedConstraintSet euclidean_set=hybrid_set.continuous_state_set(space);
    this->_set=Enclosure(euclidean_set,factory);
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const RealSpace& space,
                                 const Box& box, const IntervalFunctionModelFactoryInterface& factory)
    : _location(location), _events(), _space(space.variable_names()), _set(box,factory),
      _variables(box.dimension(),INITIAL)
{
}

HybridEnclosure::HybridEnclosure(const HybridBox& hbox, const IntervalFunctionModelFactoryInterface& factory)
    : _location(hbox.location()), _events(), _space(hbox.space().variable_names()), _set(hbox.continuous_state_set(),factory),
      _variables(hbox.continuous_state_set().dimension(),INITIAL)
{
}



HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const RealSpace& spc, const Enclosure& set)
    : _location(location), _events(), _space(spc.variable_names()), _set(set),
      _variables(catenate(List<EnclosureVariableType>(set.dimension(),INITIAL),List<EnclosureVariableType>(set.number_of_parameters()-set.dimension(),UNKNOWN)))
{
}


HybridEnclosure* HybridEnclosure::clone() const {
    return new HybridEnclosure(*this);
}

const RealSpace
HybridEnclosure::space() const
{
    return RealSpace(this->_space);
}

List<DiscreteEvent> const&
HybridEnclosure::previous_events() const
{
    return this->_events;
}

IntervalFunctionModelFactoryInterface const&
HybridEnclosure::function_factory() const
{
    return this->_set.function_factory();
}

IntervalVectorFunctionModel const&
HybridEnclosure::space_function() const
{
    return this->_set.space_function();
}

IntervalScalarFunctionModel const&
HybridEnclosure::time_function() const
{
    return this->_set.time_function();
}

IntervalScalarFunctionModel const&
HybridEnclosure::dwell_time_function() const
{
    return this->_set.dwell_time_function();
}

void HybridEnclosure::set_time_function(const IntervalScalarFunctionModel& time_function)
{
    ARIADNE_NOT_IMPLEMENTED;
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->parameter_domain(),time_function.domain()),
                       "Domain of "<<time_function<<" does not contain parameter domain "<<this->parameter_domain());
    ARIADNE_ASSERT_MSG(this->parameter_domain()==time_function.domain(),
                       "Domain of "<<time_function<<" does not equal parameter domain "<<this->parameter_domain());
    //this->_set._time_function=time_function;
}

IntervalVector
HybridEnclosure::space_bounding_box() const
{
    ARIADNE_LOG(8,"space_codomain="<<this->space_function().codomain()<<" space_range="<<this->space_function()(this->_set.reduced_domain())<<"\n");
    //return this->space_function()(this->_set.reduced_domain());
    return this->_set.bounding_box();
}

Interval
HybridEnclosure::time_range() const
{
    ARIADNE_LOG(8,"time_codomain="<<this->time_function().codomain()<<" time_range="<<this->time_function()(this->_set.reduced_domain())<<"\n");
    //return this->time_function().codomain();
    return this->time_function().evaluate(this->_set.reduced_domain());
}

Interval
HybridEnclosure::dwell_time_range() const
{
    ARIADNE_LOG(8,"dwell_time_codomain="<<this->dwell_time_function().codomain()<<
                  " dwell_time_range="<<this->dwell_time_function().evaluate(this->_set.reduced_domain())<<"\n");
    return this->dwell_time_function().evaluate(this->_set.reduced_domain());
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

const IntervalVector
HybridEnclosure::parameter_domain() const
{
    return this->_set.domain();
}

HybridBox
HybridEnclosure::bounding_box() const
{
    return HybridBox(this->_location,this->_space,this->space_bounding_box());
}


void HybridEnclosure::new_parameter(Interval ivl, EnclosureVariableType vt)
{
    this->_set.new_parameter(ivl);
    this->_variables.append(vt);
}

void HybridEnclosure::new_variable(Interval ivl, EnclosureVariableType vt)
{
    this->_set.new_variable(ivl);
    this->_variables.append(vt);
}

void HybridEnclosure::new_invariant(DiscreteEvent event, IntervalScalarFunction constraint_function) {
    this->_set.new_negative_state_constraint(constraint_function);
}


void HybridEnclosure::new_activation(DiscreteEvent event, IntervalScalarFunction constraint_function) {
    this->_set.new_positive_state_constraint(constraint_function);
}

void HybridEnclosure::new_guard(DiscreteEvent event, IntervalScalarFunction constraint_function) {
    this->_set.new_zero_state_constraint(constraint_function);
}

void HybridEnclosure::new_parameter_constraint(DiscreteEvent event, IntervalConstraint constraint) {
    this->_set.new_parameter_constraint(constraint);
}

void HybridEnclosure::new_state_constraint(DiscreteEvent event, IntervalConstraint constraint) {
    this->_set.new_state_constraint(constraint);
}


void HybridEnclosure::new_constraint(DiscreteEvent event, IntervalConstraint constraint) {
    this->new_state_constraint(event,constraint);
}


void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, RealSpace space, const IntervalVectorFunction& map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"*this="<<*this<<", event="<<event<<", space="<<space<<", target="<<target<<", map="<<map);
    this->_events.append(event);
    this->_location=target;
    this->_space=space.variable_names();
    this->_set.apply_map(map);
}

void HybridEnclosure::apply_spacetime_evolve_step(const IntervalVectorFunctionModel& phi, const IntervalScalarFunctionModel& elps)
{
    this->_set.apply_spacetime_evolve_step(phi,elps);
}

void HybridEnclosure::apply_spacetime_reach_step(const IntervalVectorFunctionModel& phi, const IntervalScalarFunctionModel& elps)
{
    this->_set.apply_spacetime_reach_step(phi,elps);
}


void HybridEnclosure::apply_evolve_step(const IntervalVectorFunctionModel& phi, const IntervalScalarFunctionModel& elps)
{
    this->_set.apply_parameter_evolve_step(phi,elps);
}

void HybridEnclosure::apply_finishing_evolve_step(const IntervalVectorFunctionModel& phi, const IntervalScalarFunctionModel& omega)
{
    this->_set.apply_finishing_parameter_evolve_step(phi,omega);
}


void HybridEnclosure::apply_reach_step(const IntervalVectorFunctionModel& phi, const IntervalScalarFunctionModel& elps)
{
    this->_set.apply_parameter_reach_step(phi,elps);
}

void HybridEnclosure::apply_full_reach_step(const IntervalVectorFunctionModel& phi)
{
    this->_set.apply_full_reach_step(phi);
}



void HybridEnclosure::bound_time(Real tmax) {
    if(this->time_range().upper()>Interval(tmax).lower()) {
        this->_set.new_negative_parameter_constraint(this->time_function()-Interval(tmax));
    }
}

void HybridEnclosure::bound_time(IntervalScalarFunction tmax) {
    this->_set.new_negative_parameter_constraint(this->time_function()-this->_set.function_factory().create(this->_set.domain(),tmax));
}

void HybridEnclosure::set_time(Real time)
{
    this->_set.new_zero_parameter_constraint(this->time_function()-Interval(time));
}

void HybridEnclosure::set_time(IntervalScalarFunction time)
{
    this->_set.new_zero_parameter_constraint(this->time_function()-this->function_factory().create(this->_set.domain(),time));
}


void HybridEnclosure::set_maximum_time(DiscreteEvent event, Float final_time)
{
    this->_set.new_negative_parameter_constraint(this->time_function()-final_time); // Deprecated
}

void HybridEnclosure::new_time_step_bound(DiscreteEvent event, IntervalScalarFunction constraint) {
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
ValidatedConstrainedImageSet HybridEnclosure::continuous_state_set() const {
    return ValidatedConstrainedImageSet(this->_set.domain(),this->space_function(),this->_constraints);
}
*/


const Enclosure&
HybridEnclosure::continuous_state_set() const {
    return this->_set;
}


uint HybridEnclosure::dimension() const {
    return this->space_function().result_size();
}

tribool HybridEnclosure::empty() const {
    return this->_set.empty();
}

tribool HybridEnclosure::inside(const HybridBox& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_state_set().inside(hbx.continuous_state_set()); }
    else { return this->continuous_state_set().empty(); }
}

tribool HybridEnclosure::separated(const HybridBox& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_state_set().separated(hbx.continuous_state_set()); }
    else { return true; }
}

tribool HybridEnclosure::satisfies(RealConstraint c) const
{
    return this->continuous_state_set().satisfies(c);
}


List<IntervalConstraint> HybridEnclosure::constraints() const {
    return this->continuous_state_set().constraints();
}

void
HybridEnclosure::restrict(const IntervalVector& subdomain)
{
    this->_set.restrict(subdomain);
}

void
HybridEnclosure::recondition()
{
    this->uniform_error_recondition();
    this->kuhn_recondition();
}

void
HybridEnclosure::uniform_error_recondition()
{
    uint old_number_of_parameters = this->number_of_parameters();
    this->_set.uniform_error_recondition();
    IntervalVector new_variables = project(this->parameter_domain(),range(old_number_of_parameters,this->number_of_parameters()));
    this->_variables.concatenate(List<EnclosureVariableType>(new_variables.size(),ERROR));
    this->_check();
}


void
HybridEnclosure::kuhn_recondition()
{
    this->_set.kuhn_recondition();
    this->_check();
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


void check_subset(const IntervalVector& dom1, const IntervalVector& dom2, const char* msg)
{
    if(dom1.size()!=dom2.size()) {
        ARIADNE_FAIL_MSG(msg<<" size("<<dom1<<")!=size("<<dom2<<")");
    } else if(!subset(dom1,dom2)) {
        ARIADNE_FAIL_MSG(msg<<" !subset("<<dom1<<","<<dom2<<")");
    }
}

void
HybridEnclosure::_check() const
{
    //this->_set.check();
    const IntervalVector& reduced_domain = this->_set.reduced_domain();
    check_subset(reduced_domain,this->_set.domain(),"domain");
    check_subset(reduced_domain,this->space_function().domain(),"function domain");
    for(uint i=0; i!=this->_set.constraints().size(); ++i) {
        check_subset(reduced_domain,this->_set.constraint(i).function().domain(),"constraint");
    }
    check_subset(reduced_domain,this->time_function().domain(),"time");
}

void HybridEnclosure::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const
{
    Projection2d proj=Ariadne::projection(this->space(),axes);
    this->continuous_state_set().draw(canvas,proj);
}

std::ostream& HybridEnclosure::write(std::ostream& os) const
{
    return os << "HybridEnclosure"
              << "( variables = " << variable_names(this->_variables)
              << ", events=" << this->_events
              << ", location=" << this->_location
              << ", space=" << this->space()
              << ", range=" << this->space_function()(this->_set.domain())
              << ", domain=" << this->_set.domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", empty=" << Ariadne::empty(this->_set.reduced_domain())
              << ", state=" << (this->_set.space_function())
              << ", constraints=" << (this->_set.constraints())
              << ", time="<< (this->_set.time_function())
              << ")";
}

std::ostream& HybridEnclosure::print(std::ostream& os) const
{
    return os << "HybridEnclosure"
              << "( events=" << this->_events
              << ", location=" << this->_location
              << ", range=" << this->space_function()(this->_set.domain())
              << ", domain=" << this->_set.domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", #constraints=" << this->_set.constraints().size()
              << ", time_range="<<this->time_range() << ")";
}

std::ostream& HybridEnclosure::repr(std::ostream& os) const
{
    return os << "HybridEnclosure"
              << "( " << this->_location
              << ", " << this->_set.domain()
              << ", " << this->_set.reduced_domain()
              << ", " << this->_set.space_function()
              << ", " << this->_set.constraints()
              << ", " << this->time_function() << ")";
}


void HybridEnclosure::adjoin_outer_approximation_to(HybridGridTreeSet& hgts, int depth) const {
    const Enclosure& set = this->continuous_state_set();
    GridTreeSet& paving = hgts[this->location()];
    set.adjoin_outer_approximation_to(paving,depth);
}

HybridGridTreeSet outer_approximation(const ListSet<HybridEnclosure>& hls, const HybridGrid& g, int depth) {
    HybridGridTreeSet result(g);
    for(ListSet<HybridEnclosure>::const_iterator iter=hls.begin(); iter!=hls.end(); ++iter) {
        result[iter->location()].adjoin_outer_approximation(iter->continuous_state_set(),depth);
    }
    return result;
}

void
draw(FigureInterface& figure, const ListSet<HybridEnclosure>& hels) {
    for(ListSet<HybridEnclosure>::const_iterator iter=hels.begin(); iter!=hels.end(); ++iter) {
        draw(figure,iter->continuous_state_set());
    }
}

ListSet<HybridEnclosure::ContinuousStateSetType>
ListSet<HybridEnclosure>::operator[](const DiscreteLocation& loc) const
{
    ListSet<HybridEnclosure::ContinuousStateSetType> result;
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->location()==loc) {
            result.adjoin(iter->continuous_state_set());
        }
    }
    return result;
}




} // namespace Ariadne
