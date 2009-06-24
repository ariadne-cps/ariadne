/***************************************************************************
 *            system_submodule.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include <iostream>
#include <iomanip>

#include "expression_interface.h"
#include "function_interface.h"
#include "formula.h"
#include "hybrid_automaton.h"
#include "hybrid_time.h"
#include "hybrid_set.h"

#include "utilities.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;



Space*
make_space(const boost::python::object& obj)
{
    Space* spcptr=new Space();
    boost::python::list elements=boost::python::extract<boost::python::list>(obj);
    int m=boost::python::len(elements);
    for(int i=0; i!=m; ++i) {
        RealVariable v=boost::python::extract<RealVariable>(elements[i]);
        *spcptr,v;
    }
    return spcptr;
}

Formula exp(RealVariable v) { return exp(Formula(v)); }
Formula log(RealVariable v) { return log(Formula(v)); }
Formula sin(RealVariable v) { return sin(Formula(v)); }
Formula cos(RealVariable v) { return cos(Formula(v)); }
Formula tan(RealVariable v) { return tan(Formula(v)); }

// Need to wrap pure virtual function explicitly for some reason...
ExpressionInterface* expression(const FormulaInterface& f, const Space& s) { return f.expression(s); }

void export_formula()
{
    class_<RealVariable> variable_class("RealVariable", init<std::string>());
    variable_class.def(init<RealVariable>());
    variable_class.def("name", &RealVariable::name, return_value_policy<reference_existing_object>());
    variable_class.def("expression", &RealVariable::expression, return_value_policy<manage_new_object>());
    variable_class.def("__eq__", &RealVariable::operator==);
    variable_class.def("__neq__", &RealVariable::operator!=);
    variable_class.def("__add__", &__add__<Formula,RealVariable,RealVariable>);
    variable_class.def("__sub__", &__sub__<Formula,RealVariable,RealVariable>);
    variable_class.def("__mul__", &__mul__<Formula,RealVariable,RealVariable>);
    variable_class.def("__div__", &__div__<Formula,RealVariable,RealVariable>);
    variable_class.def("__add__", &__add__<Formula,RealVariable,double>);
    variable_class.def("__sub__", &__sub__<Formula,RealVariable,double>);
    variable_class.def("__mul__", &__mul__<Formula,RealVariable,double>);
    variable_class.def("__div__", &__div__<Formula,RealVariable,double>);
    variable_class.def("__add__", &__add__<Formula,RealVariable,Interval>);
    variable_class.def("__sub__", &__sub__<Formula,RealVariable,Interval>);
    variable_class.def("__mul__", &__mul__<Formula,RealVariable,Interval>);
    variable_class.def("__div__", &__div__<Formula,RealVariable,Interval>);
    variable_class.def("__radd__", &__radd__<Formula,RealVariable,double>);
    variable_class.def("__rsub__", &__rsub__<Formula,RealVariable,double>);
    variable_class.def("__rmul__", &__rmul__<Formula,RealVariable,double>);
    variable_class.def("__rdiv__", &__rdiv__<Formula,RealVariable,double>);
    variable_class.def("__radd__", &__radd__<Formula,RealVariable,Interval>);
    variable_class.def("__rsub__", &__rsub__<Formula,RealVariable,Interval>);
    variable_class.def("__rmul__", &__rmul__<Formula,RealVariable,Interval>);
    variable_class.def("__rdiv__", &__rdiv__<Formula,RealVariable,Interval>);
    variable_class.def(self_ns::str(self));

    def("exp", (Formula(*)(RealVariable)) &exp);
    def("log", (Formula(*)(RealVariable)) &log);
    def("sin", (Formula(*)(RealVariable)) &sin);
    def("cos", (Formula(*)(RealVariable)) &cos);
    def("tan", (Formula(*)(RealVariable)) &tan);


    class_<Space> space_class("Space");
    space_class.def("__init__", make_constructor(&make_space) );
    space_class.def("dimension", &Space::dimension);
    space_class.def("variable", &Space::variable, return_value_policy<reference_existing_object>());
    space_class.def("index", &Space::index);
    space_class.def(self_ns::str(self));


    class_<Formula> formula_class("Formula", no_init);
    formula_class.def("expression", &Formula::expression, return_value_policy<manage_new_object>());
    formula_class.def("__add__", &__add__<Formula,Formula,Formula>);
    formula_class.def("__sub__", &__sub__<Formula,Formula,Formula>);
    formula_class.def("__mul__", &__mul__<Formula,Formula,Formula>);
    formula_class.def("__div__", &__div__<Formula,Formula,Formula>);
    formula_class.def("__add__", &__add__<Formula,Formula,double>);
    formula_class.def("__sub__", &__sub__<Formula,Formula,double>);
    formula_class.def("__mul__", &__mul__<Formula,Formula,double>);
    formula_class.def("__div__", &__div__<Formula,Formula,double>);
    formula_class.def("__add__", &__add__<Formula,Formula,Interval>);
    formula_class.def("__sub__", &__sub__<Formula,Formula,Interval>);
    formula_class.def("__mul__", &__mul__<Formula,Formula,Interval>);
    formula_class.def("__div__", &__div__<Formula,Formula,Interval>);
    formula_class.def("__add__", &__add__<Formula,Formula,RealVariable>);
    formula_class.def("__sub__", &__sub__<Formula,Formula,RealVariable>);
    formula_class.def("__mul__", &__mul__<Formula,Formula,RealVariable>);
    formula_class.def("__div__", &__div__<Formula,Formula,RealVariable>);
    formula_class.def("__radd__", &__radd__<Formula,Formula,double>);
    formula_class.def("__rsub__", &__rsub__<Formula,Formula,double>);
    formula_class.def("__rmul__", &__rmul__<Formula,Formula,double>);
    formula_class.def("__rdiv__", &__rdiv__<Formula,Formula,double>);
    formula_class.def("__radd__", &__radd__<Formula,Formula,Interval>);
    formula_class.def("__rsub__", &__rsub__<Formula,Formula,Interval>);
    formula_class.def("__rmul__", &__rmul__<Formula,Formula,Interval>);
    formula_class.def("__rdiv__", &__rdiv__<Formula,Formula,Interval>);
    formula_class.def("__radd__", &__radd__<Formula,Formula,RealVariable>);
    formula_class.def("__rsub__", &__rsub__<Formula,Formula,RealVariable>);
    formula_class.def("__rmul__", &__rmul__<Formula,Formula,RealVariable>);
    formula_class.def("__rdiv__", &__rdiv__<Formula,Formula,RealVariable>);
    formula_class.def(self_ns::str(self));

    def("exp", (Formula(*)(Formula)) &exp);
    def("log", (Formula(*)(Formula)) &log);
    def("sin", (Formula(*)(Formula)) &sin);
    def("cos", (Formula(*)(Formula)) &cos);
    def("tan", (Formula(*)(Formula)) &tan);

    def("expression", &Formula::expression, return_value_policy<manage_new_object>());

}



void export_hybrid_automaton()
{
    typedef boost::shared_ptr<const ExpressionInterface> ExpressionPtr;
    typedef boost::shared_ptr<const FunctionInterface> FunctionPtr;

/*
    class_<DiscreteState> discrete_state_class("DiscreteState",no_init);
    discrete_state_class.def(self_ns::str(self));

    class_<DiscreteEvent> discrete_event_class("DiscreteEvent",no_init);
    discrete_event_class.def(self_ns::str(self));
*/

    class_<DiscreteMode, shared_ptr<DiscreteMode> > discrete_mode_class("DiscreteMode",no_init);
    discrete_mode_class.def("location",&DiscreteMode::location);
    discrete_mode_class.def("dynamic",&DiscreteMode::dynamic_ptr);
    discrete_mode_class.def("invariants",&DiscreteMode::invariants,return_value_policy<reference_existing_object>());

    class_<DiscreteTransition, shared_ptr<DiscreteTransition> > discrete_transition_class("DiscreteTransition",no_init);
    discrete_transition_class.def("event",&DiscreteTransition::event);
    discrete_transition_class.def("source",&DiscreteTransition::source,return_value_policy<reference_existing_object>());
    discrete_transition_class.def("target",&DiscreteTransition::target,return_value_policy<reference_existing_object>());
    discrete_transition_class.def("reset",&DiscreteTransition::reset_ptr);
    discrete_transition_class.def("guard",&DiscreteTransition::activation_ptr);
    discrete_transition_class.def("urgency",&DiscreteTransition::forced);

    class_<HybridTime> hybrid_time_class("HybridTime",init<double,int>());
    //hybrid_time_class.def("continuous_time",&HybridTime::continuous_time);
    //hybrid_time_class.def("discrete_time",&HybridTime::discrete_time);

    class_<HybridAutomaton> hybrid_automaton_class("HybridAutomaton",no_init);
    hybrid_automaton_class.def("mode",&HybridAutomaton::mode,return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("transition",&HybridAutomaton::transition,return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("modes",&HybridAutomaton::modes,return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("transitions",(const std::set<DiscreteTransition>&(HybridAutomaton::*)()const) &HybridAutomaton::transitions,return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("transitions",(std::set<DiscreteTransition>(HybridAutomaton::*)(DiscreteState)const) &HybridAutomaton::transitions);
    hybrid_automaton_class.def("blocking_guards",(std::map<DiscreteEvent,FunctionPtr>(HybridAutomaton::*)(DiscreteState)const) &HybridAutomaton::blocking_guards);

}



void system_submodule() {
    export_formula();
    export_hybrid_automaton();
}

