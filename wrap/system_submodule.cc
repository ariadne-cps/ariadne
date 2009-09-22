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
#include "real.h"
#include "formula.h"
#include "hybrid_automaton.h"
#include "hybrid_time.h"
#include "hybrid_set.h"

#include "utilities.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {


template<class T>
struct from_python< Space<T> > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Space<T> >()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        boost::python::list elements=extract<boost::python::list>(obj_ptr);
        Space<T>* spc_ptr = new (storage) Space<T>();
        for(int i=0; i!=len(elements); ++i) {
            extract<String> xs(elements[i]);
            if(xs.check()) { spc_ptr->append(Variable<T>(xs())); }
            else { Variable<T> v=extract< Variable<T> >(elements[i]); spc_ptr->append(v); }
        }
        data->convertible = storage;
    }
};


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

} // namespace Ariadne


namespace Ariadne { int length(const array<std::string>& a) { return a.size(); } }

void export_formula()
{
    implicitly_convertible<RealVariable,RealExpression>();

    to_python< List<RealExpression> >();


    // TODO: These interval conversions are dangerous since they are applied when they sometimes should not be.
    //implicitly_convertible<double,RealExpression>();
    //implicitly_convertible<Interval,RealExpression>();

    class_<RealVariable> real_variable_class("RealVariable", init<std::string>());
    real_variable_class.def("__pos__", &__pos__<RealExpression,RealVariable>);
    real_variable_class.def("__neg__", &__neg__<RealExpression,RealVariable>);
    real_variable_class.def("__add__", &__add__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__sub__", &__sub__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__mul__", &__mul__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__div__", &__div__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__radd__", &__radd__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__rsub__", &__rsub__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__rmul__", &__rmul__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def("__rdiv__", &__rdiv__<RealExpression,RealVariable,RealExpression>);
    real_variable_class.def(self_ns::str(self));

    class_<RealSpace> real_space_class("RealSpace", init<RealSpace>());
    real_space_class.def("dimension", &RealSpace::dimension);
    real_space_class.def("variable", &RealSpace::variable, return_value_policy<reference_existing_object>());
    real_space_class.def("index", &RealSpace::index);
    real_space_class.def(self_ns::str(self));

    def("variable",(RealVariable(*)(const String& s)) &variable<Real>);
    def("variables",(RealSpace(*)(const List<String>& s)) &variables<Real>);

    from_python<RealSpace>();

    class_<RealExpression> real_expression_class("RealExpression", init<RealExpression>());
    real_expression_class.def("name", &RealExpression::operator_name);
    real_expression_class.def("subexpressions", &RealExpression::subexpressions);
    real_expression_class.def("substitute", &RealExpression::substitute<Real>);
    real_expression_class.def("simplify", &RealExpression::simplify);
    real_expression_class.def("__pos__", &__pos__<RealExpression,RealExpression>);
    real_expression_class.def("__neg__", &__neg__<RealExpression,RealExpression>);
    real_expression_class.def("__add__", &__add__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__sub__", &__sub__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__mul__", &__mul__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__div__", &__div__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__radd__", &__radd__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__rsub__", &__rsub__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__rmul__", &__rmul__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def("__rdiv__", &__rdiv__<RealExpression,RealExpression,RealExpression>);
    real_expression_class.def(self_ns::str(self));

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

    class_<RealAssignment> real_assignment_class("RealAssignment",no_init);
    real_assignment_class.def(self_ns::str(self));

    typedef Variable<Tribool> TriboolVariable;
    typedef Expression<Tribool> TriboolExpression;

    to_python< List<TriboolExpression> >();

    class_<TriboolVariable> tribool_variable_class("TriboolVariable", init<std::string>());
    tribool_variable_class.def(self_ns::str(self));

    class_<TriboolExpression> tribool_expression_class("TriboolExpression",init<TriboolExpression>());
    tribool_expression_class.def("name", &TriboolExpression::operator_name);
    tribool_expression_class.def("subexpressions", &TriboolExpression::subexpressions);
    tribool_expression_class.def("substitute", &TriboolExpression::substitute<Real>);
    tribool_expression_class.def("substitute", &TriboolExpression::substitute<Tribool>);
    tribool_expression_class.def("simplify", &TriboolExpression::simplify);
    tribool_expression_class.def("__and__", &__and__<TriboolExpression,TriboolExpression,TriboolExpression>);
    tribool_expression_class.def("__or__", &__or__<TriboolExpression,TriboolExpression,TriboolExpression>);
    tribool_expression_class.def("__neg__", &__not__<TriboolExpression,TriboolExpression>);
    tribool_expression_class.def(self_ns::str(self));

    implicitly_convertible<TriboolVariable,TriboolExpression>();

    def("sgn", (TriboolExpression(*)(RealExpression)) &sgn);

    real_expression_class.def("__less__", &__gte__<TriboolExpression,RealExpression,RealExpression>);
    real_expression_class.def("__gtr__", &__lte__<TriboolExpression,RealExpression,RealExpression>);
    real_expression_class.def("__gt__", &__gte__<TriboolExpression,RealExpression,double>);
    real_expression_class.def("__lt__", &__lte__<TriboolExpression,RealExpression,double>);


    //class_<RealVariable> real_variable_class("RealVariable", init<std::string>());

}



void export_hybrid_automaton()
{
    typedef boost::shared_ptr<const ScalarFunctionInterface> ExpressionPtr;
    typedef boost::shared_ptr<const VectorFunctionInterface> FunctionPtr;

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

