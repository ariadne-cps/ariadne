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



template<class T>
Space<T>*
make_space(const boost::python::object& obj)
{
    Space<T>* spcptr=new Space<T>();
    boost::python::list elements=boost::python::extract<boost::python::list>(obj);
    int m=boost::python::len(elements);
    for(int i=0; i!=m; ++i) {
        boost::python::extract<String> xs(elements[i]);
        if(xs.check()) { spcptr->append(Variable<T>(xs())); }
        else { Variable<T> v=extract< Variable<T> >(elements[i]); spcptr->append(v); }
    }
    return spcptr;
}

template<class T>
struct space_from_python_list {
    space_from_python_list() { converter::registry::push_back(&convertible,&construct,type_id< Space<T> >()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        storage = make_space<T>(extract<object>(obj_ptr));
        data->convertible = storage;
    }
};



namespace Ariadne {
    RealExpression var(const std::string& s) { return RealExpression(RealVariable(s)); }
    RealExpression operator+(const RealVariable& v) { return +RealExpression(v); }
    RealExpression operator-(const RealVariable& v) { return -RealExpression(v); }
    RealExpression operator+(const RealVariable& v, const RealExpression& e) { return RealExpression(v)+e; }
    RealExpression operator-(const RealVariable& v, const RealExpression& e) { return RealExpression(v)-e; }
    RealExpression operator*(const RealVariable& v, const RealExpression& e) { return RealExpression(v)*e; }
    RealExpression operator/(const RealVariable& v, const RealExpression& e) { return RealExpression(v)/e; }
    RealExpression operator+(const RealExpression& e, const RealVariable& v) { return e+RealExpression(v); }
    RealExpression operator-(const RealExpression& e, const RealVariable& v) { return e-RealExpression(v); }
    RealExpression operator*(const RealExpression& e, const RealVariable& v) { return e*RealExpression(v); }
    RealExpression operator/(const RealExpression& e, const RealVariable& v) { return e/RealExpression(v); }
}

namespace Ariadne { int length(const array<std::string>& a) { return a.size(); } }

void export_formula()
{
    implicitly_convertible<RealVariable,RealExpression>();
    implicitly_convertible<double,RealExpression>();
    implicitly_convertible<Interval,RealExpression>();

    class_<RealVariable> variable_class("RealVariable", init<std::string>());
    variable_class.def("__pos__", &__pos__<RealExpression,RealVariable>);
    variable_class.def("__neg__", &__neg__<RealExpression,RealVariable>);
    variable_class.def("__add__", &__add__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__sub__", &__sub__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__mul__", &__mul__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__div__", &__div__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__radd__", &__radd__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__rsub__", &__rsub__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__rmul__", &__rmul__<RealExpression,RealVariable,RealExpression>);
    variable_class.def("__rdiv__", &__rdiv__<RealExpression,RealVariable,RealExpression>);
    variable_class.def(self_ns::str(self));

    class_<RealSpace> space_class("RealSpace");
    space_class.def("__init__", make_constructor(&make_space<Real>) );
    space_class.def("dimension", &RealSpace::dimension);
    space_class.def("variable", &RealSpace::variable, return_value_policy<reference_existing_object>());
    space_class.def("index", &RealSpace::index);
    space_class.def(self_ns::str(self));

    def("variable",(RealVariable(*)(const String& s)) &variable);
    def("variables",&make_space<Real>,return_value_policy<manage_new_object>());

    space_from_python_list<Real>();

    class_<RealExpression> expression_class("RealExpression", no_init);
    expression_class.def("__pos__", &__pos__<RealExpression,RealExpression>);
    expression_class.def("__neg__", &__neg__<RealExpression,RealExpression>);
    expression_class.def("__add__", &__add__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__sub__", &__sub__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__mul__", &__mul__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__div__", &__div__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__radd__", &__radd__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__rsub__", &__rsub__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__rmul__", &__rmul__<RealExpression,RealExpression,RealExpression>);
    expression_class.def("__rdiv__", &__rdiv__<RealExpression,RealExpression,RealExpression>);
    expression_class.def(self_ns::str(self));

    def("neg", (RealExpression(*)(RealExpression)) &neg);
    def("rec", (RealExpression(*)(RealExpression)) &rec);
    def("sqr", (RealExpression(*)(RealExpression)) &sqr);
    //def("pow", (RealExpression(*)(RealExpression,int)) &pow);
    def("sqrt", (RealExpression(*)(RealExpression)) &sqrt);
    def("exp", (RealExpression(*)(RealExpression)) &exp);
    def("log", (RealExpression(*)(RealExpression)) &log);
    def("sin", (RealExpression(*)(RealExpression)) &sin);
    def("cos", (RealExpression(*)(RealExpression)) &cos);
    def("tan", (RealExpression(*)(RealExpression)) &tan);

    class_<RealAssignment> assignment_class("RealAssignment",no_init);
    assignment_class.def(self_ns::str(self));

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

