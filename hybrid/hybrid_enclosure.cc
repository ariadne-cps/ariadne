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

#include "function/functional.h"
#include "config.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "numeric/operators.h"
#include "function/function.h"
#include "function/constraint.h"
#include "function/polynomial.h"
#include "function/taylor_function.h"
#include "geometry/box.h"
#include "geometry/grid_set.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/discrete_event.h"
#include "hybrid/discrete_location.h"

#include "solvers/linear_programming.h"
#include "solvers/nonlinear_programming.h"
#include "solvers/constraint_solver.h"
#include "geometry/affine_set.h"

#include "output/graphics_interface.h"
#include "hybrid/hybrid_enclosure.h"
#include "hybrid/hybrid_set.h"


namespace Ariadne {


OutputStream& operator<<(OutputStream& os, const EnclosureVariableType& evt) {
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

template<class T> StringType str(const T& t) { StringStream ss; ss<<t; return ss.str(); }

List<String> variable_names(const List<EnclosureVariableType>& vt) {
    std::map<EnclosureVariableType,Nat> counts;
    List<String> result;
    for(Nat i=0; i!=vt.size(); ++i) {
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
    : _location(), _events(), _state_space(), _set(), _variables()
{
}

HybridEnclosure::HybridEnclosure(const HybridBoundedConstraintSet& hybrid_set,
                                 const RealSpace& state_space,
                                 const IntervalFunctionModelFactoryInterface& factory)
    : _location(hybrid_set.location()), _events(), _state_space(state_space.variable_names()), _set(),
      _variables(state_space.dimension(),INITIAL)
{
    BoundedConstraintSet euclidean_set=hybrid_set.euclidean_set(this->_location,state_space);
    this->_set=Enclosure(euclidean_set,factory);
}

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const RealSpace& state_space,
                                 const ExactBoxType& box, const IntervalFunctionModelFactoryInterface& factory)
    : _location(location), _events(), _state_space(state_space.variable_names()), _set(box,factory),
      _variables(box.dimension(),INITIAL)
{
}

HybridEnclosure::HybridEnclosure(const HybridBoxType& hbox, const IntervalFunctionModelFactoryInterface& factory)
    : _location(hbox.location()), _events(), _state_space(hbox.space().variable_names()), _set(hbox.continuous_set(),factory),
      _variables(hbox.continuous_set().dimension(),INITIAL)
{
}

template<class T> List<T> catenate(List<T> lst1, T const& val2, List<T> const& lst3) {
    lst1.append(val2); return catenate(lst1,lst3); }

HybridEnclosure::HybridEnclosure(const DiscreteLocation& location, const RealSpace& spc, const Enclosure& set)
    : _location(location), _events(), _state_space(spc.variable_names()), _set(set),
      _variables(catenate(List<EnclosureVariableType>(set.state_dimension(),INITIAL),TEMPORAL,List<EnclosureVariableType>(set.number_of_parameters()-set.state_dimension()-1u,UNKNOWN)))
{
}


HybridEnclosure* HybridEnclosure::clone() const {
    return new HybridEnclosure(*this);
}

const RealSpace
HybridEnclosure::state_time_auxiliary_space() const
{
    return join(join(this->state_space(),this->time_variable()),this->auxiliary_space());
}

const RealSpace
HybridEnclosure::state_space() const
{
    return RealSpace(this->_state_space);
}

const RealVariable
HybridEnclosure::time_variable() const
{
    return TimeVariable();
}

const RealSpace
HybridEnclosure::auxiliary_space() const
{
    return RealSpace(this->_auxiliary_space);
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

ValidatedVectorFunctionModel const&
HybridEnclosure::state_function() const
{
    return this->_set.state_function();
}

ValidatedScalarFunctionModel const&
HybridEnclosure::time_function() const
{
    return this->_set.time_function();
}

ValidatedScalarFunctionModel const&
HybridEnclosure::dwell_time_function() const
{
    return this->_set.dwell_time_function();
}

ValidatedVectorFunctionModel const
HybridEnclosure::state_time_auxiliary_function() const
{
    return this->_set.state_time_auxiliary_function();
}


Void HybridEnclosure::set_time_function(const ValidatedScalarFunctionModel& time_function)
{
    ARIADNE_NOT_IMPLEMENTED;
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->parameter_domain(),time_function.domain()),
                       "Domain of "<<time_function<<" does not contain parameter domain "<<this->parameter_domain());
    ARIADNE_ASSERT_MSG(this->parameter_domain()==time_function.domain(),
                       "Domain of "<<time_function<<" does not equal parameter domain "<<this->parameter_domain());
    //this->_set._time_function=time_function;
}

UpperBoxType
HybridEnclosure::state_bounding_box() const
{
    ARIADNE_LOG(8,"space_codomain="<<this->state_function().codomain()<<" space_range="<<apply(this->state_function(),this->_set.reduced_domain())<<"\n");
    //return this->state_function()(this->_set.reduced_domain());
    return cast_exact_box(this->_set.bounding_box());
}

UpperIntervalType
HybridEnclosure::range_of(EffectiveScalarFunction const& g) const
{
    return apply(compose(g,this->state_function()),this->_set.reduced_domain());
}

UpperIntervalType
HybridEnclosure::time_range() const
{
    ARIADNE_LOG(8,"time_codomain="<<this->time_function().codomain()<<" time_range="<<apply(this->time_function(),this->_set.reduced_domain())<<"\n");
    //return this->time_function().codomain();
    return apply(this->time_function(),this->_set.reduced_domain());
}

UpperIntervalType
HybridEnclosure::dwell_time_range() const
{
    ARIADNE_LOG(8,"dwell_time_codomain="<<this->dwell_time_function().codomain()<<
                  " dwell_time_range="<<apply(this->dwell_time_function(),this->_set.reduced_domain())<<"\n");
    return apply(this->dwell_time_function(),this->_set.reduced_domain());
}

SizeType
HybridEnclosure::number_of_constraints() const
{
    return this->_set.number_of_constraints();
}

SizeType
HybridEnclosure::number_of_parameters() const
{
    return this->_set.number_of_parameters();
}

const ExactBoxType
HybridEnclosure::parameter_domain() const
{
    return this->_set.domain();
}

HybridUpperBoxType
HybridEnclosure::bounding_box() const
{
    return HybridUpperBoxType(this->_location,this->_state_space,cast_exact_box(this->state_bounding_box()));
}


Void HybridEnclosure::set_auxiliary(List<RealVariable> vars, EffectiveVectorFunction aux)
{
    if(vars.size()!=aux.result_size()) { std::cerr<<vars<<" "<<aux<<"\n"; }
    ARIADNE_ASSERT(this->_state_space.size()==aux.argument_size());
    ARIADNE_ASSERT(vars.size()==aux.result_size());
    this->_auxiliary_space=vars;
    this->_set.set_auxiliary(aux);
}

Void HybridEnclosure::new_parameter(ExactIntervalType ivl, EnclosureVariableType vt)
{
    this->_set.new_parameter(ivl);
    this->_variables.append(vt);
}

Void HybridEnclosure::new_variable(ExactIntervalType ivl, EnclosureVariableType vt)
{
    this->_set.new_variable(ivl);
    this->_variables.append(vt);
}

Void HybridEnclosure::new_invariant(DiscreteEvent event, ValidatedScalarFunction constraint_function) {
    this->_set.new_negative_state_constraint(constraint_function);
}


Void HybridEnclosure::new_activation(DiscreteEvent event, ValidatedScalarFunction constraint_function) {
    this->_set.new_positive_state_constraint(constraint_function);
}

Void HybridEnclosure::new_guard(DiscreteEvent event, ValidatedScalarFunction constraint_function) {
    this->_set.new_zero_state_constraint(constraint_function);
}

Void HybridEnclosure::new_parameter_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->_set.new_parameter_constraint(constraint);
}

Void HybridEnclosure::new_state_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->_set.new_state_constraint(constraint);
}

Void HybridEnclosure::new_state_time_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->_set.new_state_time_constraint(constraint);
}


Void HybridEnclosure::new_constraint(DiscreteEvent event, ValidatedConstraint constraint) {
    this->new_state_constraint(event,constraint);
}


Void HybridEnclosure::clear_time()
{
    this->_set.clear_time();
}

Void HybridEnclosure::clear_events()
{
    this->_events.clear();
}

Void HybridEnclosure::apply_reset(DiscreteEvent event, DiscreteLocation target, RealSpace state_space, const ValidatedVectorFunction& map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->state_dimension(),"*this="<<*this<<", event="<<event<<", state_space="<<state_space<<", target="<<target<<", map="<<map);
    this->_events.append(event);
    this->_location=target;
    this->_state_space=state_space.variables();
    this->_set.apply_map(map);
}

Void HybridEnclosure::apply_fixed_evolve_step(const ValidatedVectorFunctionModel& phi, const Float64Value& elps)
{
    this->_set.apply_fixed_evolve_step(phi,elps);
}

Void HybridEnclosure::apply_spacetime_evolve_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_spacetime_evolve_step(phi,elps);
}

Void HybridEnclosure::apply_spacetime_reach_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_spacetime_reach_step(phi,elps);
}


Void HybridEnclosure::apply_evolve_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_parameter_evolve_step(phi,elps);
}

Void HybridEnclosure::apply_finishing_evolve_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& omega)
{
    this->_set.apply_finishing_parameter_evolve_step(phi,omega);
}


Void HybridEnclosure::apply_reach_step(const ValidatedVectorFunctionModel& phi, const ValidatedScalarFunctionModel& elps)
{
    this->_set.apply_parameter_reach_step(phi,elps);
}

Void HybridEnclosure::apply_full_reach_step(const ValidatedVectorFunctionModel& phi)
{
    this->_set.apply_full_reach_step(phi);
}



Void HybridEnclosure::bound_time(Real tmax) {
    if(possibly(this->time_range().upper()>tmax)) {
       this->_set.new_negative_parameter_constraint(this->time_function()-tmax);
    }
}

Void HybridEnclosure::bound_time(ValidatedScalarFunction tmax) {
    this->_set.new_negative_parameter_constraint(this->time_function()-this->_set.function_factory().create(this->_set.domain(),tmax));
}

Void HybridEnclosure::set_time(Real time)
{
    this->_set.new_zero_parameter_constraint(this->time_function()-time);
}

Void HybridEnclosure::set_time(ValidatedScalarFunction time)
{
    this->_set.new_zero_parameter_constraint(this->time_function()-this->function_factory().create(this->_set.domain(),time));
}


Void HybridEnclosure::set_maximum_time(DiscreteEvent event, Float64 final_time)
{
    this->_set.new_negative_parameter_constraint(this->time_function()-ValidatedNumericType(final_time)); // Deprecated
}

Void HybridEnclosure::new_time_step_bound(DiscreteEvent event, ValidatedScalarFunction constraint) {
    ARIADNE_NOT_IMPLEMENTED; // Deprecated
}

Void HybridEnclosure::set_step_time(Float64Value time)
{
    ARIADNE_NOT_IMPLEMENTED; // Deprecated
}



const DiscreteLocation& HybridEnclosure::location() const {
    return this->_location;
}

/*
ValidatedConstrainedImageSet HybridEnclosure::continuous_set() const {
    return ValidatedConstrainedImageSet(this->_set.domain(),this->state_function(),this->_constraints);
}
*/


const Enclosure&
HybridEnclosure::continuous_set() const {
    return this->_set;
}


DimensionType HybridEnclosure::dimension() const {
    return this->_set.dimension();
}

DimensionType HybridEnclosure::state_dimension() const {
    return this->_set.state_dimension();
}

ValidatedSierpinskian HybridEnclosure::is_empty() const {
    return this->_set.is_empty();
}

ValidatedSierpinskian HybridEnclosure::inside(const HybridBoxType& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_set().inside(hbx.continuous_set()); }
    else { return this->continuous_set().is_empty(); }
}

ValidatedSierpinskian HybridEnclosure::separated(const HybridBoxType& hbx) const {
    if(this->_location==hbx.location()) { return this->continuous_set().separated(hbx.continuous_set()); }
    else { return true; }
}

ValidatedKleenean HybridEnclosure::satisfies(EffectiveConstraint c) const
{
    return this->continuous_set().satisfies(c);
}


List<ValidatedConstraint> HybridEnclosure::constraints() const {
    return this->continuous_set().constraints();
}

Void
HybridEnclosure::restrict(const ExactBoxType& subdomain)
{
    this->_set.restrict(subdomain);
}

Void
HybridEnclosure::recondition()
{
    this->uniform_error_recondition();
    this->kuhn_recondition();
}

Void
HybridEnclosure::uniform_error_recondition()
{
    Nat old_number_of_parameters = this->number_of_parameters();
    this->_set.uniform_error_recondition();
    ExactIntervalVectorType new_variables = project(this->parameter_domain(),range(old_number_of_parameters,this->number_of_parameters()));
    this->_variables.concatenate(List<EnclosureVariableType>(new_variables.size(),ERROR));
    this->_check();
}


Void
HybridEnclosure::kuhn_recondition()
{
    this->_set.kuhn_recondition();
    this->_check();
}


List<HybridEnclosure>
HybridEnclosure::split() const
{
    List<ExactBoxType> subdomains = this->continuous_set().splitting_subdomains_zeroth_order();
    List<HybridEnclosure> result(subdomains.size(),*this);
    for(Nat i=0; i!=result.size(); ++i) {
        result[i].restrict(subdomains[i]);
    }
    return result;
}


Void check_subset(const ExactBoxType& dom1, const ExactBoxType& dom2, const char* msg)
{
    if(dom1.size()!=dom2.size()) {
        ARIADNE_FAIL_MSG(msg<<" size("<<dom1<<")!=size("<<dom2<<")");
    } else if(!subset(dom1,dom2)) {
        ARIADNE_FAIL_MSG(msg<<" !subset("<<dom1<<","<<dom2<<")");
    }
}

Void
HybridEnclosure::_check() const
{
    //this->_set.check();
    const ExactIntervalVectorType& reduced_domain = this->_set.reduced_domain();
    check_subset(reduced_domain,this->_set.domain(),"domain");
    check_subset(reduced_domain,this->state_function().domain(),"function domain");
    for(Nat i=0; i!=this->_set.constraints().size(); ++i) {
        check_subset(reduced_domain,this->_set.constraint(i).function().domain(),"constraint");
    }
    check_subset(reduced_domain,this->time_function().domain(),"time");
}

Void HybridEnclosure::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const
{
    Projection2d proj=Ariadne::projection(this->state_time_auxiliary_space(),axes);
    this->continuous_set().draw(canvas,proj);
}


OutputStream& HybridEnclosure::write(OutputStream& os) const
{
    return os << "HybridEnclosure"
              << "( variables = " << variable_names(this->_variables)
              << ", events=" << this->_events
              << ", location=" << this->_location
              << ", state_space=" << this->state_space()
              << ", range=" << apply(this->state_function(),this->_set.domain())
              << ", domain=" << this->_set.domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", empty=" << this->_set.reduced_domain().is_empty()
              << ", state=" << this->_set.state_function()
              << ", constraints=" << this->_set.constraints()
              << ", time="<< this->_set.time_function()
              << ")";
}

OutputStream& HybridEnclosure::print(OutputStream& os) const
{
    return os << "HybridEnclosure"
              << "( events=" << this->_events
              << ", location=" << this->_location
              << ", range=" << apply(this->state_function(),this->_set.domain())
              << ", domain=" << this->_set.domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", subdomain=" << this->_set.reduced_domain()
              << ", #constraints=" << this->_set.constraints().size()
              << ", time_range="<<this->time_range() << ")";
}

OutputStream& HybridEnclosure::repr(OutputStream& os) const
{
    return os << "HybridEnclosure"
              << "( " << this->_location
              << ", " << this->_set.domain()
              << ", " << this->_set.reduced_domain()
              << ", " << this->_set.state_function()
              << ", " << this->_set.constraints()
              << ", " << this->time_function() << ")";
}


Void HybridEnclosure::adjoin_outer_approximation_to(HybridGridTreeSet& hgts, Int depth) const {
    const Enclosure& set = this->continuous_set();
    GridTreeSet& paving = hgts[this->location()];
    set.adjoin_outer_approximation_to(paving,depth);
}

HybridGridTreeSet outer_approximation(const ListSet<HybridEnclosure>& hls, const HybridGrid& g, Int depth) {
    HybridGridTreeSet result(g);
    for(ListSet<HybridEnclosure>::ConstIterator iter=hls.begin(); iter!=hls.end(); ++iter) {
        result[iter->location()].adjoin_outer_approximation(iter->continuous_set(),depth);
    }
    return result;
}

Void
draw(FigureInterface& figure, const ListSet<HybridEnclosure>& hels) {
    for(ListSet<HybridEnclosure>::ConstIterator iter=hels.begin(); iter!=hels.end(); ++iter) {
        draw(figure,iter->continuous_set());
    }
}

ListSet<HybridEnclosure::ContinuousStateSetType>
ListSet<HybridEnclosure>::operator[](const DiscreteLocation& loc) const
{
    ListSet<HybridEnclosure::ContinuousStateSetType> result;
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->location()==loc) {
            result.adjoin(iter->continuous_set());
        }
    }
    return result;
}




} // namespace Ariadne
