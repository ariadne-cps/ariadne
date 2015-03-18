/***************************************************************************
 *            function_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "boost_python.h"
#include "utilities.h"

#include <boost/python.hpp>

#include <iostream>
#include <iomanip>

#include "utility/array.h"
#include "utility/container.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/expansion.h"
#include "algebra/multi_index.h"
#include "function/taylor_model.h"
#include "algebra/differential.h"
#include "function/polynomial.h"
#include "function/affine.h"
#include "function/taylor_function.h"
#include "function/constraint.h"
#include "function/function.h"
#include "function/function_mixin.h"
#include "expression/expression.h"
#include "expression/space.h"
#include "expression/assignment.h"
#include "expression/formula.h"
#include "expression/function_expression.h"
#include "solvers/constraint_solver.h"

#include "function/function_mixin.tcc"

using namespace boost::python;

namespace Ariadne {

template<>
struct from_python< MultiIndex >
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<MultiIndex>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::tuple tup=boost::python::extract<boost::python::tuple>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage<MultiIndex>*)   data)->storage.bytes;
        MultiIndex res(len(tup));
        for(Nat i=0; i!=res.size(); ++i) { res.set(i,boost::python::extract<Nat>(tup[i])); }
        new (storage) MultiIndex(res);
        data->convertible = storage;
    }
};




template<>
struct from_python<EffectiveVectorFunction>
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<EffectiveVectorFunction>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=boost::python::extract<list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< EffectiveVectorFunction >*)   data)->storage.bytes;
        EffectiveVectorFunction res(len(lst),0);
        for(Nat i=0; i!=res.result_size(); ++i) { res.set(i,boost::python::extract<EffectiveScalarFunction>(lst[i])); }
        new (storage) EffectiveVectorFunction(res);
        data->convertible = storage;
    }
};




template<class X, class D>
inline Matrix<X> get_jacobian(const Vector<D>& d) {
    const Nat rs=d.size(); const Nat as=d[0].argument_size();
    Matrix<X> J(rs,as);
    for(Nat i=0; i!=rs; ++i) {
        for(Nat j=0; j!=as; ++j) {
            J[i][j]=d[i][j];
        }
    }
    return J;
}

template<class X> OutputStream& operator<<(OutputStream& os, const Representation< ScalarFunction<X> >& frepr) {
    static_cast<const ScalarFunctionInterface<X>&>(frepr.reference()).repr(os); return os;
}

template<class X> OutputStream& operator<<(OutputStream& os, const Representation< VectorFunction<X> >& frepr) {
    static_cast<const VectorFunctionInterface<X>&>(frepr.reference()).repr(os); return os;
}

class ScalarPythonFunction
    : public ScalarFunctionMixin<ScalarPythonFunction,EffectiveTag>
{
  public:
    friend class ScalarFunctionMixin<ScalarPythonFunction,EffectiveTag>;
    template<class T> Void _compute(T& r, const Vector<T>& a) const {
        r=boost::python::extract<T>(this->_pyf(a)); }
  public:
    ScalarPythonFunction(StringType& nm, SizeType as, const boost::python::object& pyf) : _name(nm), _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(SizeType as, const boost::python::object& pyf) : _name(),  _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(const boost::python::object& pyf)
        : _name(),
          _argument_size(boost::python::extract<SizeType>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    ScalarPythonFunction* clone() const { return new ScalarPythonFunction(*this); }
    virtual SizeType argument_size() const { return this->_argument_size; }

    virtual EffectiveScalarFunctionInterface* _derivative (SizeType j) const {
        ARIADNE_FAIL_MSG("Cannot symbolically differentiate a Python function"); }
    virtual OutputStream& repr(OutputStream& os) const { return os; }
    virtual OutputStream& write(OutputStream& os) const {
        os << "ScalarUserFunction( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "argument_size="<<this->_argument_size;
        return os << " )"; }
    EffectiveScalarFunction derivative(SizeType j) const { return EffectiveScalarFunction(this->_derivative(j)); }
  private:
    StringType _name;
    SizeType _argument_size;
    boost::python::object _pyf;
};


class VectorPythonFunction
    : public VectorFunctionMixin<VectorPythonFunction,EffectiveTag>
{
  public:
    friend class VectorFunctionMixin<VectorPythonFunction,EffectiveTag>;
    template<class T> Void _compute(Vector<T>& r, const Vector<T>& a) const {
        r=boost::python::extract< Vector<T> >(this->_pyf(a)); }
  public:
    VectorPythonFunction(StringType& nm, SizeType rs, SizeType as, const boost::python::object& pyf) : _name(nm), _result_size(rs), _argument_size(as), _pyf(pyf) { }
    VectorPythonFunction(SizeType rs, SizeType as, const object& pyf) : _name(), _result_size(rs), _argument_size(as), _pyf(pyf) { }
    VectorPythonFunction(const object& pyf)
        : _name(),
          _result_size(boost::python::extract<Int>(pyf.attr("result_size"))),
          _argument_size(boost::python::extract<Int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    VectorPythonFunction* clone() const { return new VectorPythonFunction(*this); }
    virtual SizeType result_size() const { return this->_result_size; }
    virtual SizeType argument_size() const { return this->_argument_size; }

    virtual EffectiveVectorFunctionInterface* _derivative (SizeType j) const {
        ARIADNE_FAIL_MSG("Cannot symbolically differentiate a Python function"); }
    virtual EffectiveScalarFunctionInterface* _get(SizeType i) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }
    virtual EffectiveScalarFunction operator[](SizeType i) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }

    virtual OutputStream& write(OutputStream& os) const {
        os << "VectorUserFunction( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "result_size="<<this->_result_size;
        os << ", argument_size="<<this->_argument_size;
        return os << " )"; }
  private:
    StringType _name;
    SizeType _result_size;
    SizeType _argument_size;
    boost::python::object _pyf;
};



}

using namespace Ariadne;

typedef ApproximateFloat64 F;
typedef ExactIntervalType I;
typedef Vector<ApproximateFloat64> FV;
typedef Vector<ExactIntervalType> IV;
typedef Matrix<ApproximateFloat64> FMx;
typedef Matrix<ExactIntervalType> IMx;
typedef Vector< Differential<ApproximateFloat64> > FSDV;
typedef Vector< Differential<ExactIntervalType> > ISDV;
typedef Vector<ValidatedTaylorModel> TMV;
typedef VectorTaylorFunction TFM;
typedef ValidatedTaylorModel TM;

template<class X> using Monomial = ExpansionValue<X>;





Void export_multi_index()
{
    class_< MultiIndex > multi_index_class("MultiIndex", init<Nat>());
    multi_index_class.def(init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def(self_ns::str(self));
    multi_index_class.def(self_ns::repr(self));

    from_python<MultiIndex>();
    to_python< List<MultiIndex> >();
}

template<class X>
Void export_monomial()
{
    class_< ExpansionValue<X> > monomial_class(python_name<X>("Monomial"), init<MultiIndex,X>());
    monomial_class.def("key",(const MultiIndex&(ExpansionValue<X>::*)()const)&ExpansionValue<X>::key,return_value_policy<copy_const_reference>());
    monomial_class.def("data",(const X&(ExpansionValue<X>::*)()const) &ExpansionValue<X>::data,return_value_policy<copy_const_reference>());
    monomial_class.def(self_ns::str(self));
}

template<class X>
Void export_polynomial()
{
    X real;

    class_< Polynomial<X> > polynomial_class(python_name<X>("Polynomial"), init< Polynomial<X> >());
    polynomial_class.def(init<Nat>());
    polynomial_class.def("constant", (Polynomial<X>(*)(Nat,double)) &Polynomial<X>::constant);
    polynomial_class.staticmethod("constant");
    polynomial_class.def("variable", (Polynomial<X>(*)(Nat,Nat)) &Polynomial<X>::variable);
    polynomial_class.staticmethod("variable");
    polynomial_class.def("coordinate", (Polynomial<X>(*)(Nat,Nat)) &Polynomial<X>::variable);
    polynomial_class.staticmethod("coordinate");
    polynomial_class.def("variables", (Vector< Polynomial<X> >(*)(Nat)) &Polynomial<X>::variables);
    polynomial_class.staticmethod("variables");

    polynomial_class.def("argument_size", &Polynomial<X>::argument_size);
    polynomial_class.def("insert", &Polynomial<X>::insert);
    polynomial_class.def(+self);
    polynomial_class.def(-self);
    polynomial_class.def(self+self);
    polynomial_class.def(self-self);
    polynomial_class.def(self*self);
    polynomial_class.def(self+real);
    polynomial_class.def(self-real);
    polynomial_class.def(self*real);
    polynomial_class.def(self/real);
    polynomial_class.def(real+self);
    polynomial_class.def(real-self);
    polynomial_class.def(real*self);
    polynomial_class.def(self_ns::str(self));
    //polynomial_class.def(self_ns::repr(self));

    to_python< Vector< Polynomial<X> > >();
}

Void export_scalar_function()
{
    class_<EffectiveScalarFunction>
        scalar_function_class("EffectiveScalarFunction", init<EffectiveScalarFunction>());
    scalar_function_class.def(init<Nat>());
    scalar_function_class.def("argument_size", &EffectiveScalarFunction::argument_size);
    scalar_function_class.def("derivative", &EffectiveScalarFunction::derivative);
    scalar_function_class.def("__call__", (BoundedFloat64(EffectiveScalarFunction::*)(const Vector<BoundedFloat64>&)const)&EffectiveScalarFunction::operator() );
    scalar_function_class.def("__call__", (ApproximateFloat64(EffectiveScalarFunction::*)(const Vector<ApproximateFloat64>&)const)&EffectiveScalarFunction::operator() );
    scalar_function_class.def("__call__", (Differential<BoundedFloat64>(EffectiveScalarFunction::*)(const Vector<Differential<BoundedFloat64>>&)const)&EffectiveScalarFunction::evaluate );
    scalar_function_class.def("__call__", (Differential<ApproximateFloat64>(EffectiveScalarFunction::*)(const Vector<Differential<ApproximateFloat64>>&)const)&EffectiveScalarFunction::evaluate );
    scalar_function_class.def("gradient", (Covector<BoundedFloat64>(EffectiveScalarFunction::*)(const Vector<BoundedFloat64>&)const)&EffectiveScalarFunction::gradient );
    scalar_function_class.def("gradient", (Covector<BoundedFloat64>(EffectiveScalarFunction::*)(const Vector<BoundedFloat64>&)const)&EffectiveScalarFunction::gradient );
    scalar_function_class.def("differential", (Differential<BoundedFloat64>(EffectiveScalarFunction::*)(const Vector<BoundedFloat64>&,DegreeType)const) &EffectiveScalarFunction::differential);
    scalar_function_class.def("differential", (Differential<ApproximateFloat64>(EffectiveScalarFunction::*)(const Vector<ApproximateFloat64>&,DegreeType)const) &EffectiveScalarFunction::differential);
    scalar_function_class.def("__pos__", &__pos__<EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__neg__", &__neg__<EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__add__", &__add__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__sub__", &__sub__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__mul__", &__mul__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__div__", &__div__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__add__", &__add__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__sub__", &__sub__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__mul__", &__mul__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__div__", &__div__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__radd__", &__radd__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__rsub__", &__rsub__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__rmul__", &__rmul__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__rdiv__", &__rdiv__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__eq__", &__eq__<EffectiveConstraint,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__le__", &__le__<EffectiveConstraint,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__ge__", &__ge__<EffectiveConstraint,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__str__", &__cstr__<EffectiveScalarFunction>);
    scalar_function_class.def("__repr__", &__crepr__<EffectiveScalarFunction>);

    scalar_function_class.def("constant", (EffectiveScalarFunction(*)(SizeType,Real)) &EffectiveScalarFunction::constant);
    scalar_function_class.def("coordinate", (EffectiveScalarFunction(*)(SizeType,SizeType)) &EffectiveScalarFunction::coordinate);
    scalar_function_class.staticmethod("constant");
    scalar_function_class.staticmethod("coordinate");

    def("evaluate", (ApproximateFloat64(*)(const EffectiveScalarFunction&,const Vector<ApproximateFloat64>&)) &evaluate);
    def("evaluate", (BoundedFloat64(*)(const EffectiveScalarFunction&,const Vector<BoundedFloat64>&)) &evaluate);

    def("derivative", (EffectiveScalarFunction(EffectiveScalarFunction::*)(SizeType)const) &EffectiveScalarFunction::derivative);

    def("pow", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&,Int)) &pow);
    def("rec", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &rec);
    def("sqr", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &sqr);
    def("sqrt", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &sqrt);
    def("exp", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &exp);
    def("log", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &log);
    def("sin", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &sin);
    def("cos", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &cos);
    def("tan", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &tan);

    def("lie_derivative", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&,const EffectiveVectorFunction&)) &lie_derivative);

    //scalar_function_class.def("__call__", (RealScalarFunctionExpression(EffectiveScalarFunction::*)(const Vector<RealVariable>&)const)&EffectiveScalarFunction::operator() );

    scalar_function_class.def("__call__", (RealExpression(*)(EffectiveScalarFunction const&, Vector<RealVariable> const&)) &evaluate);
    def("evaluate", (RealExpression(*)(EffectiveScalarFunction const&, Vector<RealVariable> const&)) &evaluate);


    class_<ValidatedScalarFunction> interval_scalar_function_class("ValidatedScalarFunction", init<ValidatedScalarFunction>());
    interval_scalar_function_class.def(init<EffectiveScalarFunction>());
    interval_scalar_function_class.def(init<Nat>());
    interval_scalar_function_class.def("argument_size", &ValidatedScalarFunction::argument_size);
    interval_scalar_function_class.def("__str__", &__cstr__<ValidatedScalarFunction>);
    interval_scalar_function_class.def("__repr__", &__crepr__<ValidatedScalarFunction>);

    implicitly_convertible<EffectiveScalarFunction,ValidatedScalarFunction>();
}


Void export_vector_function()
{

    class_<EffectiveVectorFunction>
        vector_function_class("EffectiveVectorFunction", init<EffectiveVectorFunction>());
    vector_function_class.def(init<Nat,Nat>());

    vector_function_class.def("result_size", &EffectiveVectorFunction::result_size);
    vector_function_class.def("argument_size", &EffectiveVectorFunction::argument_size);
    vector_function_class.def("__getitem__", &EffectiveVectorFunction::get);
    vector_function_class.def("__setitem__", &EffectiveVectorFunction::set);
    vector_function_class.def("__call__", (Vector<BoundedFloat64>(EffectiveVectorFunction::*)(const Vector<BoundedFloat64>&)const)&EffectiveVectorFunction::operator() );
    vector_function_class.def("__call__", (Vector<ApproximateFloat64>(EffectiveVectorFunction::*)(const Vector<ApproximateFloat64>&)const)&EffectiveVectorFunction::operator() );
    vector_function_class.def("__call__", (Vector<Differential<BoundedFloat64>>(EffectiveVectorFunction::*)(const Vector<Differential<BoundedFloat64>>&)const)&EffectiveVectorFunction::evaluate );
    vector_function_class.def("__call__", (Vector<Differential<ApproximateFloat64>>(EffectiveVectorFunction::*)(const Vector<Differential<ApproximateFloat64>>&)const)&EffectiveVectorFunction::evaluate );
    vector_function_class.def("jacobian", (Matrix<BoundedFloat64>(EffectiveVectorFunction::*)(const Vector<BoundedFloat64>&)const) &EffectiveVectorFunction::jacobian);
    vector_function_class.def("jacobian", (Matrix<ApproximateFloat64>(EffectiveVectorFunction::*)(const Vector<ApproximateFloat64>&)const) &EffectiveVectorFunction::jacobian);
    vector_function_class.def("differential", (Vector<Differential<BoundedFloat64> >(EffectiveVectorFunction::*)(const Vector<BoundedFloat64>&,DegreeType)const) &EffectiveVectorFunction::differential);
    vector_function_class.def("differential", (Vector<Differential<ApproximateFloat64> >(EffectiveVectorFunction::*)(const Vector<ApproximateFloat64>&,DegreeType)const) &EffectiveVectorFunction::differential);
    vector_function_class.def("__str__", &__cstr__<EffectiveVectorFunction>);
    vector_function_class.def("__repr__", &__crepr__<EffectiveVectorFunction>);

    vector_function_class.def("identity", (EffectiveVectorFunction(*)(SizeType)) &EffectiveVectorFunction::identity);
    vector_function_class.staticmethod("identity");

    def("evaluate", (Vector<ApproximateFloat64>(*)(const EffectiveVectorFunction&,const Vector<ApproximateFloat64>&)) &evaluate);
    def("evaluate", (Vector<BoundedFloat64>(*)(const EffectiveVectorFunction&,const Vector<BoundedFloat64>&)) &evaluate);

    def("join", (EffectiveVectorFunction(*)(const EffectiveScalarFunction&, const EffectiveScalarFunction&)) &join);
    def("join", (EffectiveVectorFunction(*)(const EffectiveVectorFunction&, const EffectiveScalarFunction&)) &join);
    def("join", (EffectiveVectorFunction(*)(const EffectiveScalarFunction&, const EffectiveVectorFunction&)) &join);
    def("join", (EffectiveVectorFunction(*)(const EffectiveVectorFunction&, const EffectiveVectorFunction&)) &join);

    def("compose", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&,const EffectiveVectorFunction&)) &compose);
    def("compose", (EffectiveVectorFunction(*)(const EffectiveVectorFunction&,const EffectiveVectorFunction&)) &compose);

    from_python<EffectiveVectorFunction>();

    class_<ValidatedVectorFunction> interval_vector_function_class("ValidatedVectorFunction", init<ValidatedVectorFunction>());
    interval_vector_function_class.def(init<EffectiveVectorFunction>());
    interval_vector_function_class.def(init<Nat,Nat>());
    interval_vector_function_class.def("result_size", &ValidatedVectorFunction::result_size);
    interval_vector_function_class.def("argument_size", &ValidatedVectorFunction::argument_size);
    interval_vector_function_class.def("__getitem__", &ValidatedVectorFunction::get);
    interval_vector_function_class.def("__setitem__", &ValidatedVectorFunction::set);
    interval_vector_function_class.def("__str__", &__cstr__<ValidatedVectorFunction>);
    interval_vector_function_class.def("__repr__", &__crepr__<ValidatedVectorFunction>);

    implicitly_convertible<EffectiveVectorFunction,ValidatedVectorFunction>();
}


Void export_scalar_python_function()
{
    class_<ScalarPythonFunction, bases< EffectiveScalarFunctionInterface > > scalar_python_function_class("ScalarUserFunction", init<object>());
    scalar_python_function_class.def(init<Nat,object>());
}

Void export_vector_python_function()
{
    class_<VectorPythonFunction, bases< EffectiveVectorFunctionInterface > > vector_python_function_class("VectorUserFunction", init<object>());
    vector_python_function_class.def(init<Nat,Nat,object>());
}


Void function_submodule() {
    to_python< Array<StringType> >();
    from_python< Array<StringType> >();

    export_multi_index();

    export_scalar_function();
    export_vector_function();

    //export_scalar_python_function();
    //export_vector_python_function();
}

