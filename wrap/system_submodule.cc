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
        spcptr->append(v.name());
    }
    return spcptr;
}


void export_formula()
{
    implicitly_convertible<RealVariable,RealExpression>();
    implicitly_convertible<double,RealExpression>();
    implicitly_convertible<Interval,RealExpression>();

    class_<Space> space_class("Space");
    space_class.def("__init__", make_constructor(&make_space) );
    space_class.def("dimension", &Space::dimension);
    space_class.def("variable", &Space::variable, return_value_policy<reference_existing_object>());
    space_class.def("index", &Space::index);
    space_class.def(self_ns::str(self));


    class_<RealExpression> expression_class("RealExpression", no_init);
    expression_class.def("__neg__", &__neg__<RealExpression,RealExpression>);
    expression_class.def("__add__", &__add__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__sub__", &__sub__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__mul__", &__mul__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__div__", &__div__<RealExpression,RealExpression,RealExpression>);
    expression_class.def(self_ns::str(self));

    def("neg", (RealExpression(*)(RealExpression)) &neg);
    def("rec", (RealExpression(*)(RealExpression)) &rec);
    def("sqr", (RealExpression(*)(RealExpression)) &sqr);
    def("pow", (RealExpression(*)(RealExpression,int)) &pow);
    def("sqrt", (RealExpression(*)(RealExpression)) &sqrt);
    def("exp", (RealExpression(*)(RealExpression)) &exp);
    def("log", (RealExpression(*)(RealExpression)) &log);
    def("sin", (RealExpression(*)(RealExpression)) &sin);
    def("cos", (RealExpression(*)(RealExpression)) &cos);
    def("tan", (RealExpression(*)(RealExpression)) &tan);


}



void export_hybrid_automaton()
{
    typedef boost::shared_ptr<const ScalarFunctionInterface> ExpressionPtr;
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

