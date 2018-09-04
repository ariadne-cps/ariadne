/***************************************************************************
 *            symbolic_submodule.cpp
 *
 *  Copyright 2009--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
#include "utilities.hpp"

#include <iostream>
#include <iomanip>
#include <functional>

#include "utility/tribool.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "geometry/function_set.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression_set.hpp"

using namespace Ariadne;


namespace Ariadne {

template<class T> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Variable<T>>& repr) {
    Variable<T> const& var=repr.reference();
    return os << class_name<T>() << "Variable(\"" << var.name() << "\")";
}
template<class T> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Expression<T>>& repr) {
    Expression<T> const& expr=repr.reference();
    return os << class_name<T>() << "Expression(" << expr << ")";
}

template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const Map<Variable<Y>,Expression<Y>>& a) {
    List<Assignment<Variable<Y>,Expression<Y>>> lst_a;
    for(auto key_val : a) { lst_a.append(let(key_val.first)=key_val.second); }
    return substitute(e,lst_a);
}

template<class T, class Y> Expression<T> substitute(const Expression<T>& e, const Map<Variable<Y>,Y>& a) {
    List<Assignment<Variable<Y>,Expression<Y>>> lst_a;
    for(auto key_val : a) { lst_a.append(let(key_val.first)=Expression<Y>(key_val.second)); }
    return substitute(e,lst_a);
}

/*
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

template<class X>
struct from_python< Vector<X>> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Vector<X> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*) data)->storage.bytes;
        Vector<X> res(static_cast<SizeType>(len(lst)));
        for(Nat i=0; i!=res.size(); ++i) { res[i]=boost::python::extract<X>(lst[i]); }
        new (storage) Vector<X>(res);
        data->convertible = storage;
    }
};
*/

LetIntegerVariable let(const IntegerVariable&);
LetRealVariable let(const RealVariable&);
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
RealExpression operator+(const RealVariable& v1, const RealVariable& v2) { return RealExpression(v1)+RealExpression(v2); }
RealExpression operator-(const RealVariable& v1, const RealVariable& v2) { return RealExpression(v1)-RealExpression(v2); }
RealExpression operator*(const RealVariable& v1, const RealVariable& v2) { return RealExpression(v1)*RealExpression(v2); }
RealExpression operator/(const RealVariable& v1, const RealVariable& v2) { return RealExpression(v1)/RealExpression(v2); }
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

//RealVariableInterval operator&&(RealVariableLowerInterval elivl, RealVariableUpperInterval euivl) {
//    assert(elivl.variable()==euivl.variable());
//    return RealVariableInterval(elivl.lower(),elivl.variable(),euivl.upper());
//}

RealVariableInterval operator|(RealVariable v, RealInterval ivl) {
    return RealVariableInterval(ivl.lower(),v,ivl.upper());
}

void foo(Map<int,double>const&) { }
} // namespace Ariadne



Void export_formula(pybind11::module& module)
{

//    to_python< List<RealExpression> >();

//    from_python< List<RealVariableInterval> >();
//    from_python< List<ContinuousPredicate> >();

    // TODO: These interval conversions are dangerous since they are applied when they sometimes should not be.
    //pybind11::implicitly_convertible<double,RealExpression>();
    //pybind11::implicitly_convertible<ExactIntervalType,RealExpression>();

    pybind11::class_<StringVariable> string_variable_class(module,"StringVariable");
    string_variable_class.def(pybind11::init<StringType>());
    string_variable_class.def("__eq__", &__eq__<Expression<Boolean>,StringVariable,StringType>);
    string_variable_class.def("__ne__", &__ne__<Expression<Boolean>,StringVariable,StringType>);
    string_variable_class.def("__str__", &__cstr__<StringVariable>);
    string_variable_class.def("__repr__", &__repr__<StringVariable>);

    pybind11::class_<StringExpression> string_expression_class(module,"StringExpression");
    string_expression_class.def(pybind11::init<StringExpression>());
    string_expression_class.def("__str__", &__cstr__<StringExpression>);
    string_expression_class.def("__repr__", &__repr__<StringExpression>);

    pybind11::class_<PrimedStringVariable> string_next_variable_class(module,"PrimedStringVariable");
    string_next_variable_class.def("__lshift__", (PrimedStringAssignment(PrimedStringVariable::*)(const StringExpression&)const) &PrimedStringVariable::operator=);
    string_next_variable_class.def("__str__", &__cstr__<PrimedStringVariable>);

    module.def("next", (PrimedStringVariable(*)(const StringVariable&)) &next);


    pybind11::class_<IntegerVariable> integer_variable_class(module,"IntegerVariable");
    integer_variable_class.def(pybind11::init<StringType>());
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
    integer_variable_class.def("__str__", &__cstr__<IntegerVariable>);
    integer_variable_class.def("__repr__", &__repr__<IntegerVariable>);

    pybind11::class_<LetIntegerVariable> let_integer_variable_class(module,"LetIntegerVariable");
    let_integer_variable_class.def("__lshift__", (IntegerAssignment(LetIntegerVariable::*)(const IntegerExpression&)const) &LetIntegerVariable::operator=);
    module.def("let", (LetIntegerVariable(*)(const IntegerVariable&)) &let);


    pybind11::class_<PrimedIntegerVariable> integer_next_variable_class(module,"PrimedIntegerVariable");
    integer_next_variable_class.def("__lshift__", (PrimedIntegerAssignment(PrimedIntegerVariable::*)(const IntegerExpression&)const) &PrimedIntegerVariable::operator=);
    integer_next_variable_class.def("__str__", &__cstr__<PrimedIntegerVariable>);
    module.def("next", (PrimedIntegerVariable(*)(const IntegerVariable&)) &next);

    pybind11::class_<IntegerExpression> integer_expression_class(module,"IntegerExpression");
    integer_expression_class.def(pybind11::init<IntegerExpression>());
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
    integer_expression_class.def("__str__", &__cstr__<IntegerExpression>);
    integer_expression_class.def("__repr__", &__repr__<IntegerExpression>);


    pybind11::class_<RealVariable> real_variable_class(module,"RealVariable");
    real_variable_class.def(pybind11::init<StringType>());
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
    real_variable_class.def("__pow__", &__pow__<RealExpression,RealVariable,Int>);

    real_variable_class.def("__le__", &__le__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__ge__", &__ge__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__lt__", &__lt__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__gt__", &__gt__<ContinuousPredicate,RealVariable,RealExpression>);
    real_variable_class.def("__str__", &__cstr__<RealVariable>);
    real_variable_class.def("__repr__", &__repr__<RealVariable>);

    pybind11::class_<RealVariables> real_variables_class(module,"RealVariables");
    real_variables_class.def(pybind11::init<StringType,SizeType>());
    real_variables_class.def("__getitem__", &__getitem__<RealVariables,SizeType,RealVariable>);
    real_variables_class.def("__str__", &__cstr__<RealVariables>);

    pybind11::class_<LetRealVariable> real_let_variable_class(module,"LetRealVariable");
    real_let_variable_class.def("__lshift__", (RealAssignment(LetRealVariable::*)(const RealExpression&)const) &LetRealVariable::operator=);
    module.def("let", (LetRealVariable(*)(const RealVariable&)) &let);

    pybind11::class_<DottedRealVariable> real_dotted_variable_class(module,"DottedRealVariable");
    real_dotted_variable_class.def("__lshift__", (DottedRealAssignment(DottedRealVariable::*)(const RealExpression&)const) &DottedRealVariable::operator=);
    real_dotted_variable_class.def("__str__", &__cstr__<DottedRealVariable>);
    module.def("dot", (DottedRealVariable(*)(const RealVariable&)) &dot);

    pybind11::class_<PrimedRealVariable> real_next_variable_class(module,"PrimedRealVariable");
    real_next_variable_class.def("__lshift__", (PrimedRealAssignment(PrimedRealVariable::*)(const RealExpression&)const) &PrimedRealVariable::operator=);
    real_next_variable_class.def("__str__", &__cstr__<PrimedRealVariable>);
    module.def("next", (PrimedRealVariable(*)(const RealVariable&)) &next);

    pybind11::class_<RealSpace> real_space_class(module,"RealSpace");
    real_space_class.def(pybind11::init<RealSpace>());
    real_space_class.def("dimension", &RealSpace::dimension);
    real_space_class.def("variable", &RealSpace::variable);
    real_space_class.def("index", (SizeType(RealSpace::*)(const Identifier&)const) &RealSpace::index);
    real_space_class.def("index", (SizeType(RealSpace::*)(const RealVariable&)const) &RealSpace::index);
    real_next_variable_class.def("__str__", &__cstr__<RealSpace>);

//    from_python<RealSpace>();

    pybind11::class_<RealExpression> real_expression_class(module,"RealExpression");
    real_expression_class.def(pybind11::init<RealExpression>());
    real_expression_class.def(pybind11::init<Real>());
    real_expression_class.def(pybind11::init<RealVariable>());
    module.def("simplify", (RealExpression(*)(const RealExpression&)) &simplify);
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
    real_expression_class.def("__str__", &__cstr__<RealExpression>);
    real_expression_class.def("__repr__", &__repr__<RealExpression>);

    real_expression_class.def("__add__", &__add__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__sub__", &__sub__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__mul__", &__mul__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__div__", &__div__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__radd__", &__radd__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__rsub__", &__rsub__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__rmul__", &__rmul__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__rdiv__", &__rdiv__<RealExpression,RealExpression,Real>);
    real_expression_class.def("__le__", &__le__<ContinuousPredicate,RealExpression,Real>);
    real_expression_class.def("__ge__", &__ge__<ContinuousPredicate,RealExpression,Real>);
    real_expression_class.def("__le__", &__le__<ContinuousPredicate,Real,RealExpression>);
    real_expression_class.def("__ge__", &__ge__<ContinuousPredicate,Real,RealExpression>);

    module.def("neg", (RealExpression(*)(RealExpression const&)) &neg);
    module.def("rec", (RealExpression(*)(RealExpression const&)) &rec);
    module.def("sqr", (RealExpression(*)(RealExpression const&)) &sqr);
    module.def("pow", (RealExpression(*)(RealExpression const&,Int)) &pow);
    module.def("sqrt", (RealExpression(*)(RealExpression const&)) &sqrt);
    module.def("exp", (RealExpression(*)(RealExpression const&)) &exp);
    module.def("log", (RealExpression(*)(RealExpression const&)) &log);
    module.def("sin", (RealExpression(*)(RealExpression const&)) &sin);
    module.def("cos", (RealExpression(*)(RealExpression const&)) &cos);
    module.def("tan", (RealExpression(*)(RealExpression const&)) &tan);

    pybind11::class_<RealAssignment> real_assignment_class(module,"RealAssignment");
    real_assignment_class.def("__str__",&__cstr__<RealAssignment>);
    pybind11::class_<DottedRealAssignment> dotted_real_assignment_class(module,"DottedRealAssignment");
    dotted_real_assignment_class.def("__str__",&__cstr__<DottedRealAssignment>);
    pybind11::class_<PrimedRealAssignment> primed_real_assignment_class(module,"PrimedRealAssignment");
    primed_real_assignment_class.def("__str__",&__cstr__<PrimedRealAssignment>);
    pybind11::class_<PrimedStringAssignment> primed_string_assignment_class(module,"PrimedStringAssignment");
    primed_string_assignment_class.def("__str__",&__cstr__<PrimedStringAssignment>);
    pybind11::class_<IntegerAssignment> integer_assignment_class(module,"IntegerAssignment");
    integer_assignment_class.def("__str__",&__cstr__<IntegerAssignment>);
    pybind11::class_<PrimedIntegerAssignment> primed_integer_assignment_class(module,"PrimedIntegerAssignment");
    primed_integer_assignment_class.def("__str__",&__cstr__<PrimedIntegerAssignment>);
    
//    from_python<Pair<RealVariable,RealExpression>>();
//    from_python<Map<RealVariable,RealExpression>>();

//    from_python<Pair<RealVariable,Real>>();
//    from_python<Map<RealVariable,Real>>();

//    module.def("substitute",(RealExpression(*)(RealExpression const&, const List<Assignment<RealVariable,RealExpression>>&)) &substitute);
    module.def("substitute",(RealExpression(*)(RealExpression const&, const Map<RealVariable,RealExpression>&)) &substitute);
    module.def("substitute",(RealExpression(*)(RealExpression const&, const Map<RealVariable,Real>&)) &substitute);

//    to_python< List<KleeneanExpression> >();

    pybind11::class_<DiscretePredicate> discrete_predicate_class(module,"DiscretePredicate");
    discrete_predicate_class.def(pybind11::init<DiscretePredicate>());
    discrete_predicate_class.def(pybind11::init<Bool>());
    discrete_predicate_class.def("__and__", &__and__<DiscretePredicate,DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__or__", &__or__<DiscretePredicate,DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__invert__", &__not__<DiscretePredicate,DiscretePredicate>);
    discrete_predicate_class.def("__str__",&__cstr__<DiscretePredicate>);

    pybind11::class_<ContinuousPredicate> continuous_predicate_class(module,"ContinuousPredicate");
    continuous_predicate_class.def(pybind11::init<ContinuousPredicate>());
    continuous_predicate_class.def(pybind11::init<Kleenean>());
    continuous_predicate_class.def("__and__", &__and__<ContinuousPredicate,ContinuousPredicate,ContinuousPredicate>);
    continuous_predicate_class.def("__or__", &__or__<ContinuousPredicate,ContinuousPredicate,ContinuousPredicate>);
    continuous_predicate_class.def("__invert__", &__not__<ContinuousPredicate,ContinuousPredicate>);
    continuous_predicate_class.def("__str__",&__cstr__<ContinuousPredicate>);

    pybind11::class_<KleeneanVariable> tribool_variable_class(module,"KleeneanVariable");
    tribool_variable_class.def(pybind11::init<StringType>());
    tribool_variable_class.def("__str__",&__cstr__<KleeneanVariable>);

//    from_python<Vector<RealExpression>>();
    module.def("make_function", (RealScalarUnivariateFunction(*)(RealVariable const&, RealExpression const&)) &make_function);
    module.def("make_function", (RealScalarFunction(*)(RealSpace const&, RealExpression const&)) &make_function);
    module.def("make_function", (RealVectorFunction(*)(RealSpace const&, Vector<RealExpression> const&)) &make_function);

    /*
    pybind11::class_<KleeneanExpression> tribool_expression_class(module,"KleeneanExpression",pybind11::init<KleeneanExpression>());
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

    pybind11::implicitly_convertible<KleeneanVariable,ContinuousPredicate>();

    module.def("sgn", (KleeneanExpression(*)(RealExpression const&)) &sgn);



    //pybind11::class_<RealVariable> real_variable_class(module,"RealVariable", pybind11::init<StringType>());

    pybind11::class_<RealVariableLowerInterval> real_variable_lower_interval_class(module,"RealVariableLowerInterval");
    real_variable_lower_interval_class.def(pybind11::init<Real,RealVariable>());
    real_variable_lower_interval_class.def("__str__",&__cstr__<RealVariableLowerInterval>);

    pybind11::class_<RealVariableUpperInterval> real_variable_upper_interval_class(module,"RealVariableUpperInterval");
    real_variable_upper_interval_class.def(pybind11::init<RealVariable,Real>());
    real_variable_upper_interval_class.def("__str__",&__cstr__<RealVariableUpperInterval>);

    real_variable_class.def("__le__", &__le__<RealVariableLowerInterval,Real,RealVariable>);
    real_variable_class.def("__ge__", &__ge__<RealVariableLowerInterval,RealVariable,Real>);
    real_variable_class.def("__le__", &__le__<RealVariableUpperInterval,RealVariable,Real>);
    real_variable_lower_interval_class.def("__le__", &__le__<RealVariableInterval,RealVariableLowerInterval,Real>);
    real_variable_upper_interval_class.def("__le__", &__le__<RealVariableInterval,Real,RealVariableUpperInterval>);
    real_variable_upper_interval_class.def("__ge__", &__ge__<RealVariableInterval,RealVariableUpperInterval,Real>);

//    real_variable_lower_interval_class.def("__and__", &__and__<RealVariableInterval,RealVariableLowerInterval,RealVariableUpperInterval>);

    real_variable_class.def("__or__", &__bitor__<RealVariableInterval,RealVariable,RealInterval>);

    pybind11::class_<RealVariableInterval> real_variable_interval_class(module,"RealVariableInterval");
    real_variable_interval_class.def(pybind11::init<Real,RealVariable,Real>());
    real_variable_interval_class.def("variable", &RealVariableInterval::variable);
    real_variable_interval_class.def("interval", &RealVariableInterval::interval);
    real_variable_interval_class.def("lower", &RealVariableInterval::lower);
    real_variable_interval_class.def("upper", &RealVariableInterval::upper);
    real_variable_interval_class.def("__str__",&__cstr__<RealVariableInterval>);

    pybind11::class_<RealVariablesBox> real_variables_box_class(module,"RealVariablesBox");
    real_variables_box_class.def(pybind11::init<List<RealVariableInterval>>());
    real_variables_box_class.def("__getitem__", &RealVariablesBox::operator[]);
    real_variables_box_class.def("__str__",&__cstr__<RealVariablesBox>);


    pybind11::class_<RealExpressionConstraintSet> real_expression_constraint_set_class(module,"RealExpressionConstraintSet");
    real_expression_constraint_set_class.def(pybind11::init<List<ContinuousPredicate>>());
    real_expression_constraint_set_class.def("__str__",&__cstr__<RealExpressionConstraintSet>);

    pybind11::class_<RealExpressionBoundedConstraintSet> real_expression_bounded_constraint_set_class(module,"RealExpressionBoundedConstraintSet");
    real_expression_bounded_constraint_set_class.def(pybind11::init<List<RealVariableInterval>,List<ContinuousPredicate>>());
    real_expression_bounded_constraint_set_class.def(pybind11::init<RealVariablesBox,RealExpressionConstraintSet>());
    real_expression_bounded_constraint_set_class.def("__str__",&__cstr__<RealExpressionBoundedConstraintSet>);
    
    module.def("make_box", (RealBox(*)(RealSpace const&, RealVariablesBox const&)) &make_box);
    module.def("make_set", (RealBox(*)(RealSpace const&, RealVariablesBox const&)) &make_box);
    module.def("make_set", (ConstraintSet(*)(RealSpace const&, RealExpressionConstraintSet const&)) &make_set);
    module.def("make_set", (BoundedConstraintSet(*)(RealSpace const&, RealExpressionBoundedConstraintSet const&)) &make_set);
    module.def("make_set", (BoundedConstraintSet(*)(RealSpace const&, RealVariablesBox const&, RealExpressionConstraintSet const&)) &make_set);


//    pybind11::implicitly_convertible<String,StringExpression>();
//    pybind11::implicitly_convertible<StringVariable,StringExpression>();

    pybind11::implicitly_convertible<Int,IntegerExpression>();
    pybind11::implicitly_convertible<Integer,IntegerExpression>();
    pybind11::implicitly_convertible<IntegerVariable,IntegerExpression>();

    // FIXME: Allowing this conversion over-eagerly prevent Real in RealVector([...])
    // pybind11::implicitly_convertible<Real,RealExpression>();
    pybind11::implicitly_convertible<RealVariable,RealExpression>();

    pybind11::implicitly_convertible<List<RealVariableInterval>,RealVariablesBox>();
    pybind11::implicitly_convertible<List<ContinuousPredicate>,RealExpressionConstraintSet>();

}



Void symbolic_submodule(pybind11::module& module) {
    export_formula(module);
}

