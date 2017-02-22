/***************************************************************************
 *            system_submodule.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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
#include <functional>

#include <boost/python.hpp>

#include "utility/tribool.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "expression/expression.hpp"
#include "expression/space.hpp"
#include "hybrid/discrete_event.hpp"
#include "hybrid/hybrid_automata.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_set.hpp"

#include "boost_python.hpp"
#include "utilities.hpp"

using namespace boost::python;
using namespace Ariadne;


namespace Ariadne {


template<class T>
struct from_python< Space<T> > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Space<T> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((converter::rvalue_from_python_storage<ExactIntervalType>*)data)->storage.bytes;
        boost::python::list elements=boost::python::extract<boost::python::list>(obj_ptr);
        Space<T>* spc_ptr = new (storage) Space<T>();
        for(Int i=0; i!=len(elements); ++i) {
            boost::python::extract<String> xs(elements[i]);
            if(xs.check()) { spc_ptr->append(Variable<T>(xs())); }
            else { Variable<T> v=boost::python::extract< Variable<T> >(elements[i]); spc_ptr->append(v); }
        }
        data->convertible = storage;
    }
};

template<>
struct from_python< Set<DiscreteEvent> > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Set<DiscreteEvent> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((converter::rvalue_from_python_storage< Set<DiscreteEvent> >*)data)->storage.bytes;
        boost::python::list elements=boost::python::extract<boost::python::list>(obj_ptr);
        Set<DiscreteEvent>* evnts_ptr = new (storage) Set<DiscreteEvent>();
        for(Int i=0; i!=len(elements); ++i) {
            boost::python::extract<String> xs(elements[i]);
            if(xs.check()) { evnts_ptr->insert(DiscreteEvent(xs())); }
            else { DiscreteEvent e=boost::python::extract< DiscreteEvent >(elements[i]); evnts_ptr->insert(e); }
        }
        data->convertible = storage;
    }
};

template<class X>
struct from_python< Vector<X>> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Vector<X> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*) data)->storage.bytes;
        Vector<X> res(len(lst));
        for(Nat i=0; i!=res.size(); ++i) { res[i]=boost::python::extract<X>(lst[i]); }
        new (storage) Vector<X>(res);
        data->convertible = storage;
    }
};

template<class T> Nat __hash__(const T&);
template<> Nat __hash__<DiscreteEvent>(const DiscreteEvent& e) {
    return reinterpret_cast<const ushort&>(e.name().c_str()[0]); }
template<> Nat __hash__<DiscreteLocation>(const DiscreteLocation& q) {
    return reinterpret_cast<const ushort&>(to_string(q).c_str()[0]); }

DottedRealVariable dot(const RealVariable&);

RealExpression var(const StringType& s) { return RealExpression(RealVariable(s)); }
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
RealExpression operator+(const RealVariable& v, const Real& x) { return RealExpression(v)+RealExpression(x); }
RealExpression operator-(const RealVariable& v, const Real& x) { return RealExpression(v)-RealExpression(x); }
RealExpression operator*(const RealVariable& v, const Real& x) { return RealExpression(v)*RealExpression(x); }
RealExpression operator/(const RealVariable& v, const Real& x) { return RealExpression(v)/RealExpression(x); }
RealExpression operator+(const Real& x, const RealVariable& v) { return RealExpression(x)+RealExpression(v); }
RealExpression operator-(const Real& x, const RealVariable& v) { return RealExpression(x)-RealExpression(v); }
RealExpression operator*(const Real& x, const RealVariable& v) { return RealExpression(x)*RealExpression(v); }
RealExpression operator/(const Real& x, const RealVariable& v) { return RealExpression(x)/RealExpression(v); }
RealExpression neg(const RealExpression&);
RealExpression rec(const RealExpression&);
RealExpression sqr(const RealExpression&);
RealExpression pow(const RealExpression&,Int);
RealExpression sqrt(const RealExpression&);
RealExpression exp(const RealExpression&);
RealExpression log(const RealExpression&);
RealExpression sin(const RealExpression&);
RealExpression cos(const RealExpression&);
RealExpression tan(const RealExpression&);
KleeneanExpression sgn(const RealExpression&);

Int length(const Array<StringType>& a) { return a.size(); }

} // namespace Ariadne



Void export_formula()
{
    implicitly_convertible<String,StringExpression>();
    implicitly_convertible<StringVariable,StringExpression>();

    implicitly_convertible<Int,IntegerExpression>();
    implicitly_convertible<Integer,IntegerExpression>();
    implicitly_convertible<IntegerVariable,IntegerExpression>();

    implicitly_convertible<RealVariable,RealExpression>();

    to_python< List<RealExpression> >();

    from_python< List<DiscreteEvent> >();

/*
    implicitly_convertible<Event,EventSet>();

    class_<Event> event_class("Event", init<StringType>());
    event_class.def(self_ns::str(self));

    class_<EventSet> event_set_class("EventSet", init<EventSet>());
    event_set_class.def(init<>());
    event_set_class.def("__invert__", &__not__<EventSet,EventSet>);
    event_set_class.def(self_ns::str(self));

    from_python<EventSet>();
*/

    // TODO: These interval conversions are dangerous since they are applied when they sometimes should not be.
    //implicitly_convertible<double,RealExpression>();
    //implicitly_convertible<ExactIntervalType,RealExpression>();

    class_<StringVariable> string_variable_class("StringVariable", init<StringType>());
    string_variable_class.def("__eq__", &__eq__<Expression<Boolean>,StringVariable,StringType>);
    string_variable_class.def("__ne__", &__ne__<Expression<Boolean>,StringVariable,StringType>);
    string_variable_class.def(self_ns::str(self));

    class_<StringExpression> string_expression_class("StringExpression", init<StringExpression>());
    string_expression_class.def(self_ns::str(self));

    class_<PrimedStringVariable> string_next_variable_class("PrimedStringVariable", no_init);
    string_next_variable_class.def("__lshift__", (PrimedStringAssignment(PrimedStringVariable::*)(const StringExpression&)const) &PrimedStringVariable::operator=);
    string_next_variable_class.def(self_ns::str(self));

    def("next", (PrimedStringVariable(*)(const StringVariable&)) &next);


    class_<IntegerVariable> integer_variable_class("IntegerVariable", init<StringType>());
    integer_variable_class.def("__lshift__", (IntegerAssignment(LetIntegerVariable::*)(const IntegerExpression&)const) &LetIntegerVariable::operator=);
    integer_variable_class.def("__pos__", &__pos__<IntegerExpression,IntegerVariable>);
    integer_variable_class.def("__neg__", &__neg__<IntegerExpression,IntegerVariable>);
    integer_variable_class.def("__add__", &__add__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__sub__", &__sub__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__mul__", &__mul__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__radd__", &__radd__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__rsub__", &__rsub__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__rmul__", &__rmul__<IntegerExpression,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__eq__", &__eq__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__ne__", &__ne__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__le__", &__le__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__ge__", &__ge__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__lt__", &__lt__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def("__gt__", &__gt__<DiscretePredicate,IntegerVariable,IntegerExpression>);
    integer_variable_class.def(self_ns::str(self));

    class_<PrimedIntegerVariable> integer_next_variable_class("PrimedIntegerVariable", no_init);
    integer_next_variable_class.def("__lshift__", (PrimedIntegerAssignment(PrimedIntegerVariable::*)(const IntegerExpression&)const) &PrimedIntegerVariable::operator=);
    integer_next_variable_class.def(self_ns::str(self));

    def("next", (PrimedIntegerVariable(*)(const IntegerVariable&)) &next);

    class_<IntegerExpression> integer_expression_class("IntegerExpression", init<IntegerExpression>());
    integer_expression_class.def("__pos__", &__pos__<IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__neg__", &__neg__<IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__add__", &__add__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__sub__", &__sub__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__mul__", &__mul__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__radd__", &__radd__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__rsub__", &__rsub__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__rmul__", &__rmul__<IntegerExpression,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__eq__", &__eq__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__ne__", &__ne__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__le__", &__le__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__ge__", &__ge__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__lt__", &__lt__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def("__gt__", &__gt__<DiscretePredicate,IntegerExpression,IntegerExpression>);
    integer_expression_class.def(self_ns::str(self));


    class_<RealVariable> real_variable_class("RealVariable", init<StringType>());
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
    real_variable_class.def("__add__", &__add__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__sub__", &__sub__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__mul__", &__mul__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__div__", &__div__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__radd__", &__radd__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__rsub__", &__rsub__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__rmul__", &__rmul__<RealExpression,RealVariable,Real>);
    real_variable_class.def("__rdiv__", &__rdiv__<RealExpression,RealVariable,Real>);

    real_variable_class.def("__le__", &__le__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__ge__", &__ge__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__lt__", &__lt__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__gt__", &__gt__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__lshift__", (RealAssignment(LetRealVariable::*)(const RealExpression&)const) &LetRealVariable::operator=);
    real_variable_class.def(self_ns::str(self));

    class_<RealVariables> real_variables_class("RealVariables", init<StringType,SizeType>());
    real_variables_class.def("__getitem__", &__getitem__<RealVariables,SizeType,RealVariable>);
    real_variables_class.def(self_ns::str(self));

    class_<DottedRealVariable> real_dotted_variable_class("DottedRealVariable", no_init);
    real_dotted_variable_class.def("__lshift__", (DottedRealAssignment(DottedRealVariable::*)(const RealExpression&)const) &DottedRealVariable::operator=);
    real_dotted_variable_class.def(self_ns::str(self));

    def("dot", (DottedRealVariable(*)(const RealVariable&)) &dot);

    class_<PrimedRealVariable> real_next_variable_class("PrimedRealVariable", no_init);
    real_next_variable_class.def("__lshift__", (PrimedRealAssignment(PrimedRealVariable::*)(const RealExpression&)const) &PrimedRealVariable::operator=);
    real_next_variable_class.def(self_ns::str(self));

    def("next", (PrimedRealVariable(*)(const RealVariable&)) &next);

    class_<RealSpace> real_space_class("RealSpace", init<RealSpace>());
    real_space_class.def("dimension", &RealSpace::dimension);
    real_space_class.def("variable", &RealSpace::variable);
    real_space_class.def("index", (SizeType(RealSpace::*)(const Identifier&)const) &RealSpace::index);
    real_space_class.def("index", (SizeType(RealSpace::*)(const RealVariable&)const) &RealSpace::index);
    real_space_class.def(self_ns::str(self));

    from_python<RealSpace>();

    class_<RealExpression> real_expression_class("RealExpression", init<RealExpression>());
    real_expression_class.def("simplify", (RealExpression(*)(const RealExpression&)) &simplify);
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
    real_expression_class.def("__le__", &__le__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def("__ge__", &__ge__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def("__lt__", &__lt__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def("__gt__", &__gt__<ContinuousPredicate,RealExpression,RealExpression>);
    //real_expression_class.def("__cmp__", &__cmp__<ContinuousPredicate,RealExpression,RealExpression>);
    real_expression_class.def(self_ns::str(self));

    real_expression_class.def("__add__", &__add__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__sub__", &__sub__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__mul__", &__mul__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__div__", &__div__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__radd__", &__radd__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__rsub__", &__rsub__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__rmul__", &__rmul__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__rdiv__", &__rdiv__<RealExpression,RealExpression,Real>);

    def("neg", (RealExpression(*)(RealExpression const&)) &neg);
    def("rec", (RealExpression(*)(RealExpression const&)) &rec);
    def("sqr", (RealExpression(*)(RealExpression const&)) &sqr);
    def("pow", (RealExpression(*)(RealExpression const&,Int)) &pow);
    def("sqrt", (RealExpression(*)(RealExpression const&)) &sqrt);
    def("exp", (RealExpression(*)(RealExpression const&)) &exp);
    def("log", (RealExpression(*)(RealExpression const&)) &log);
    def("sin", (RealExpression(*)(RealExpression const&)) &sin);
    def("cos", (RealExpression(*)(RealExpression const&)) &cos);
    def("tan", (RealExpression(*)(RealExpression const&)) &tan);

    class_<RealAssignment> real_assignment_class("RealAssignment",no_init);
    real_assignment_class.def(self_ns::str(self));
    class_<DottedRealAssignment> dotted_real_assignment_class("DottedRealAssignment",no_init);
    dotted_real_assignment_class.def(self_ns::str(self));
    class_<PrimedRealAssignment> primed_real_assignment_class("PrimedRealAssignment",no_init);
    primed_real_assignment_class.def(self_ns::str(self));
    class_<PrimedStringAssignment> primed_string_assignment_class("PrimedStringAssignment",no_init);
    primed_string_assignment_class.def(self_ns::str(self));
    class_<IntegerAssignment> integer_assignment_class("IntegerAssignment",no_init);
    integer_assignment_class.def(self_ns::str(self));
    class_<PrimedIntegerAssignment> primed_integer_assignment_class("PrimedIntegerAssignment",no_init);
    primed_integer_assignment_class.def(self_ns::str(self));

    to_python< List<KleeneanExpression> >();

    class_<DiscretePredicate> discrete_predicate_class("DiscretePredicate", init<DiscretePredicate>());
    discrete_predicate_class.def(init<Bool>());
    discrete_predicate_class.def("__and__", &__and__<DiscretePredicate,DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__or__", &__or__<DiscretePredicate,DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__invert__", &__not__<DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def(self_ns::str(self));

    class_<ContinuousPredicate> continuous_predicate_class("ContinuousPredicate", init<ContinuousPredicate>());
    continuous_predicate_class.def(init<Kleenean>());
    continuous_predicate_class.def("__and__", &__and__<ContinuousPredicate,ContinuousPredicate,ContinuousPredicate>);
    continuous_predicate_class.def("__or__", &__or__<ContinuousPredicate,ContinuousPredicate,ContinuousPredicate>);
    continuous_predicate_class.def("__invert__", &__not__<ContinuousPredicate,ContinuousPredicate>);
    continuous_predicate_class.def(self_ns::str(self));

    class_<KleeneanVariable> tribool_variable_class("KleeneanVariable", init<StringType>());
    tribool_variable_class.def(self_ns::str(self));

    from_python<Vector<RealExpression>>();
    def("make_function", (RealScalarUnivariateFunction(*)(RealVariable const&, RealExpression const&)) &make_function);
    def("make_function", (RealScalarFunction(*)(RealSpace const&, RealExpression const&)) &make_function);
    def("make_function", (RealVectorFunction(*)(RealSpace const&, Vector<RealExpression> const&)) &make_function);

    /*
    class_<KleeneanExpression> tribool_expression_class("KleeneanExpression",init<KleeneanExpression>());
    tribool_expression_class.def("name", &KleeneanExpression::operator_name);
    tribool_expression_class.def("subexpressions", &KleeneanExpression::subexpressions);
    tribool_expression_class.def("substitute", &KleeneanExpression::substitute<Real>);
    tribool_expression_class.def("substitute", &KleeneanExpression::substitute<Kleenean>);
    tribool_expression_class.def("simplify", &KleeneanExpression::simplify);
    tribool_expression_class.def("__and__", &__and__<KleeneanExpression,KleeneanExpression,KleeneanExpression>);
    tribool_expression_class.def("__or__", &__or__<KleeneanExpression,KleeneanExpression,KleeneanExpression>);
    tribool_expression_class.def("__neg__", &__not__<KleeneanExpression,KleeneanExpression>);
    tribool_expression_class.def(self_ns::str(self));
    */

    implicitly_convertible<KleeneanVariable,ContinuousPredicate>();

    def("sgn", (KleeneanExpression(*)(RealExpression)) &sgn);



    //class_<RealVariable> real_variable_class("RealVariable", init<StringType>());

}


Void export_hybrid_automaton()
{
    // Don't use return_value_policy<copy_const_reference> since reference lifetime should not exceed automaton lifetime

    to_python< Set<DiscreteEvent> >();
    to_python< Set<DiscreteLocation> >();

    class_<DiscreteLocation> discrete_state_class("DiscreteLocation",init<DiscreteLocation>());
    discrete_state_class.def("__eq__", &__eq__<Bool,DiscreteLocation,DiscreteLocation>);
    discrete_state_class.def("__ne__", &__ne__<Bool,DiscreteLocation,DiscreteLocation>);
    discrete_state_class.def("__hash__", &__hash__<DiscreteLocation>);
    discrete_state_class.def(self_ns::str(self));

    class_<DiscreteEvent> discrete_event_class("DiscreteEvent",init<DiscreteEvent>());
    discrete_event_class.def("__eq__", &__eq__<Bool,DiscreteEvent,DiscreteEvent>);
    discrete_event_class.def("__ne__", &__ne__<Bool,DiscreteEvent,DiscreteEvent>);
    discrete_event_class.def("__hash__", &__hash__<DiscreteEvent>);
    discrete_event_class.def(self_ns::str(self));
    implicitly_convertible<Int,DiscreteEvent>();
    implicitly_convertible<StringType,DiscreteEvent>();


    class_<HybridTime> hybrid_time_class("HybridTime",init<double,Int>());
    hybrid_time_class.def("continuous_time",&HybridTime::continuous_time,return_value_policy<copy_const_reference>());
    hybrid_time_class.def("discrete_time",&HybridTime::discrete_time,return_value_policy<copy_const_reference>());

    class_<HybridAutomaton> hybrid_automaton_class("HybridAutomaton",init<>());
    hybrid_automaton_class.def("locations", &HybridAutomaton::locations);
    hybrid_automaton_class.def("events", &HybridAutomaton::events);
    hybrid_automaton_class.def("event_kind", &HybridAutomaton::event_kind);
    hybrid_automaton_class.def("dynamic_function", &HybridAutomaton::dynamic_function);
    hybrid_automaton_class.def("guard_function", &HybridAutomaton::guard_function);
    hybrid_automaton_class.def("reset_function", &HybridAutomaton::reset_function);
    hybrid_automaton_class.def("new_mode",(Void(HybridAutomaton::*)(DiscreteLocation,List<DottedRealAssignment> const&)) &HybridAutomaton::new_mode, return_value_policy<reference_existing_object>());
    hybrid_automaton_class.def("new_invariant",(Void(HybridAutomaton::*)(DiscreteLocation,ContinuousPredicate const&,DiscreteEvent)) &HybridAutomaton::new_invariant);
    hybrid_automaton_class.def("new_transition",(Void(HybridAutomaton::*)(DiscreteLocation,DiscreteEvent,DiscreteLocation,List<PrimedRealAssignment> const&,ContinuousPredicate const&,EventKind)) &HybridAutomaton::new_transition);
    hybrid_automaton_class.def(self_ns::str(self));

}



Void system_submodule() {
    export_formula();
    export_hybrid_automaton();
}

