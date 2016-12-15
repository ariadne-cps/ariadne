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

typedef ScalarFunction<EffectiveTag> ESF;
ESF pow(ESF const&, Int);
ESF pos(ESF const&); ESF neg(ESF const&); ESF sqr(ESF const&); ESF rec(ESF const&);
ESF sqrt(ESF const&); ESF exp(ESF const&); ESF log(ESF const&); ESF atan(ESF const&);
ESF sin(ESF const&); ESF cos(ESF const&); ESF tan(ESF const&);

typedef ScalarFunction<ValidatedTag> VSF;
VSF pow(VSF const&, Int);
VSF pos(VSF const&); VSF neg(VSF const&); VSF sqr(VSF const&); VSF rec(VSF const&);
VSF sqrt(VSF const&); VSF exp(VSF const&); VSF log(VSF const&); VSF atan(VSF const&);
VSF sin(VSF const&); VSF cos(VSF const&); VSF tan(VSF const&);

/*
RealExpression evaluate(EffectiveScalarFunction const& f, Vector<RealVariable> const& vars);
Float64Approximation evaluate(const ScalarFunction<ValidatedTag>&, const Vector<Float64Approximation>&);
Float64Approximation evaluate(const ScalarFunction<EffectiveTag>&, const Vector<Float64Approximation>&);
Float64Bounds evaluate(const ScalarFunction<ValidatedTag>&, const Vector<Float64Bounds>&);
Float64Bounds evaluate(const ScalarFunction<EffectiveTag>&, const Vector<Float64Bounds>&);
Vector<Float64Approximation> evaluate(const VectorFunction<ValidatedTag>&, const Vector<Float64Approximation>&);
Vector<Float64Approximation> evaluate(const VectorFunction<EffectiveTag>&, const Vector<Float64Approximation>&);
Vector<Float64Bounds> evaluate(const VectorFunction<ValidatedTag>&, const Vector<Float64Bounds>&);
Vector<Float64Bounds> evaluate(const VectorFunction<EffectiveTag>&, const Vector<Float64Bounds>&);
*/

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

typedef Float64Approximation F;
typedef ExactIntervalType I;
typedef Vector<Float64Approximation> FV;
typedef Vector<ExactIntervalType> IV;
typedef Matrix<Float64Approximation> FMx;
typedef Matrix<ExactIntervalType> IMx;
typedef Vector< Differential<Float64Approximation> > FSDV;
typedef Vector< Differential<ExactIntervalType> > ISDV;
typedef Vector<ValidatedTaylorModel64> TMV;
typedef VectorTaylorFunction TFM;
typedef ValidatedTaylorModel64 TM;

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


Void export_univariate_function()
{
    class_<EffectiveScalarUnivariateFunction>
        function_class("EffectiveScalarUnivariateFunction", init<EffectiveScalarUnivariateFunction>());
    function_class.def("__call__", (Float64Bounds(EffectiveScalarUnivariateFunction::*)(const Float64Bounds&)const)&EffectiveScalarUnivariateFunction::operator() );
    function_class.def("__call__", (Differential<Float64Bounds>(EffectiveScalarUnivariateFunction::*)(const Differential<Float64Bounds>&)const)&EffectiveScalarUnivariateFunction::operator() );

    function_class.def("constant", (EffectiveScalarUnivariateFunction(*)(IntervalDomain,EffectiveNumber)) &EffectiveScalarUnivariateFunction::constant);
    function_class.def("coordinate", (EffectiveScalarUnivariateFunction(*)()) &EffectiveScalarUnivariateFunction::coordinate);
    function_class.staticmethod("constant");
    function_class.staticmethod("coordinate");

}


template<class P> Void export_scalar_function()
{
    class_<ScalarFunction<P>>
        scalar_function_class((class_name<P>()+"ScalarFunction").c_str(), init<ScalarFunction<P>>());
    scalar_function_class.def(init<SizeType>());
    scalar_function_class.def("argument_size", &ScalarFunction<P>::argument_size);
    scalar_function_class.def("derivative", &ScalarFunction<P>::derivative);
    scalar_function_class.def("__pos__", &__pos__<ScalarFunction<P>,ScalarFunction<P>>);
    scalar_function_class.def("__neg__", &__neg__<ScalarFunction<P>,ScalarFunction<P>>);
    scalar_function_class.def("__add__", &__add__<ScalarFunction<P>,ScalarFunction<P>,ScalarFunction<P>>);
    scalar_function_class.def("__sub__", &__sub__<ScalarFunction<P>,ScalarFunction<P>,ScalarFunction<P>>);
    scalar_function_class.def("__mul__", &__mul__<ScalarFunction<P>,ScalarFunction<P>,ScalarFunction<P>>);
    scalar_function_class.def("__div__", &__div__<ScalarFunction<P>,ScalarFunction<P>,ScalarFunction<P>>);
    scalar_function_class.def("__add__", &__add__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__sub__", &__sub__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__mul__", &__mul__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__div__", &__div__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__radd__", &__radd__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__rsub__", &__rsub__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__rmul__", &__rmul__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);
    scalar_function_class.def("__rdiv__", &__rdiv__<ScalarFunction<P>,ScalarFunction<P>,Number<P>>);

//FIXME
//    scalar_function_class.def("__eq__", &__eq__<Constraint<ScalarFunction<P>,Number<P>>,ScalarFunction<P>,Number<P>>);
//    scalar_function_class.def("__le__", &__le__<Constraint<ScalarFunction<P>,Number<P>>,ScalarFunction<P>,Number<P>>);
//    scalar_function_class.def("__ge__", &__ge__<Constraint<ScalarFunction<P>,Number<P>>,ScalarFunction<P>,Number<P>>);

    scalar_function_class.def("__call__", (Float64Bounds(ScalarFunction<P>::*)(const Vector<Float64Bounds>&)const)&ScalarFunction<P>::operator() );
    scalar_function_class.def("__call__", (Float64Approximation(ScalarFunction<P>::*)(const Vector<Float64Approximation>&)const)&ScalarFunction<P>::operator() );
    scalar_function_class.def("__call__", (Differential<Float64Bounds>(ScalarFunction<P>::*)(const Vector<Differential<Float64Bounds>>&)const)&ScalarFunction<P>::evaluate );
    scalar_function_class.def("__call__", (Differential<Float64Approximation>(ScalarFunction<P>::*)(const Vector<Differential<Float64Approximation>>&)const)&ScalarFunction<P>::evaluate );
    scalar_function_class.def("gradient", (Covector<Float64Bounds>(ScalarFunction<P>::*)(const Vector<Float64Bounds>&)const)&ScalarFunction<P>::gradient );
    scalar_function_class.def("gradient", (Covector<Float64Bounds>(ScalarFunction<P>::*)(const Vector<Float64Bounds>&)const)&ScalarFunction<P>::gradient );
    scalar_function_class.def("differential", (Differential<Float64Bounds>(ScalarFunction<P>::*)(const Vector<Float64Bounds>&,DegreeType)const) &ScalarFunction<P>::differential);
    scalar_function_class.def("differential", (Differential<Float64Approximation>(ScalarFunction<P>::*)(const Vector<Float64Approximation>&,DegreeType)const) &ScalarFunction<P>::differential);

    scalar_function_class.def("__str__", &__cstr__<ScalarFunction<P>>);
    scalar_function_class.def("__repr__", &__crepr__<ScalarFunction<P>>);

    scalar_function_class.def("constant", (ScalarFunction<P>(*)(SizeType,Number<P>)) &ScalarFunction<P>::constant);
    scalar_function_class.def("coordinate", (ScalarFunction<P>(*)(SizeType,SizeType)) &ScalarFunction<P>::coordinate);
    scalar_function_class.staticmethod("constant");
    scalar_function_class.staticmethod("coordinate");

    def("pow", (ScalarFunction<P>(*)(const ScalarFunction<P>&,Int)) &pow);
    def("rec", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &rec);
    def("sqr", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &sqr);
    def("sqrt", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &sqrt);
    def("exp", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &exp);
    def("log", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &log);
    def("sin", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &sin);
    def("cos", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &cos);
    def("tan", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &tan);

    def("evaluate", (Float64Approximation(*)(const ScalarFunction<P>&,const Vector<Float64Approximation>&)) &evaluate<P,IntervalDomain,Float64Approximation>);
    def("evaluate", (Float64Bounds(*)(const ScalarFunction<P>&,const Vector<Float64Bounds>&)) &evaluate<P,IntervalDomain,Float64Bounds>);

    def("derivative", (ScalarFunction<P>(ScalarFunction<P>::*)(SizeType)const) &ScalarFunction<P>::derivative);



}


template<class P> Void export_vector_function()
{
    typedef VectorFunction<P> VF;

    class_<VectorFunction<P>>
        vector_function_class("VectorFunction<P>", init<VectorFunction<P>>());
    vector_function_class.def(init<Nat,Nat>());

    vector_function_class.def("result_size", &VectorFunction<P>::result_size);
    vector_function_class.def("argument_size", &VectorFunction<P>::argument_size);
    vector_function_class.def("__getitem__", &VectorFunction<P>::get);
    vector_function_class.def("__setitem__", &VectorFunction<P>::set);
    vector_function_class.def("__call__", (Vector<Float64Bounds>(VectorFunction<P>::*)(const Vector<Float64Bounds>&)const)&VectorFunction<P>::operator() );
    vector_function_class.def("__call__", (Vector<Float64Approximation>(VectorFunction<P>::*)(const Vector<Float64Approximation>&)const)&VectorFunction<P>::operator() );
    vector_function_class.def("__call__", (Vector<Differential<Float64Bounds>>(VectorFunction<P>::*)(const Vector<Differential<Float64Bounds>>&)const)&VectorFunction<P>::evaluate );
    vector_function_class.def("__call__", (Vector<Differential<Float64Approximation>>(VectorFunction<P>::*)(const Vector<Differential<Float64Approximation>>&)const)&VectorFunction<P>::evaluate );
    vector_function_class.def("jacobian", (Matrix<Float64Bounds>(VectorFunction<P>::*)(const Vector<Float64Bounds>&)const) &VectorFunction<P>::jacobian);
    vector_function_class.def("jacobian", (Matrix<Float64Approximation>(VectorFunction<P>::*)(const Vector<Float64Approximation>&)const) &VectorFunction<P>::jacobian);
    vector_function_class.def("differential", (Vector<Differential<Float64Bounds> >(VectorFunction<P>::*)(const Vector<Float64Bounds>&,DegreeType)const) &VectorFunction<P>::differential);
    vector_function_class.def("differential", (Vector<Differential<Float64Approximation> >(VectorFunction<P>::*)(const Vector<Float64Approximation>&,DegreeType)const) &VectorFunction<P>::differential);
    vector_function_class.def("__str__", &__cstr__<VectorFunction<P>>);
    vector_function_class.def("__repr__", &__crepr__<VectorFunction<P>>);

    vector_function_class.def("identity", (VectorFunction<P>(*)(SizeType)) &VectorFunction<P>::identity);
    vector_function_class.staticmethod("identity");

    def("evaluate", (Vector<Float64Approximation>(*)(const VectorFunction<P>&,const Vector<Float64Approximation>&)) &evaluate<P,BoxDomain,Float64Approximation>);
    def("evaluate", (Vector<Float64Bounds>(*)(const VectorFunction<P>&,const Vector<Float64Bounds>&)) &evaluate<P,BoxDomain,Float64Bounds>);

    def("join", (VectorFunction<P>(*)(const ScalarFunction<P>&, const ScalarFunction<P>&)) &join);
    def("join", (VectorFunction<P>(*)(const VectorFunction<P>&, const ScalarFunction<P>&)) &join);
    def("join", (VectorFunction<P>(*)(const ScalarFunction<P>&, const VectorFunction<P>&)) &join);
    def("join", (VectorFunction<P>(*)(const VectorFunction<P>&, const VectorFunction<P>&)) &join);

    def("compose", (ScalarFunction<P>(*)(const ScalarFunction<P>&,const VectorFunction<P>&)) &compose);
    def("compose", (VectorFunction<P>(*)(const VectorFunction<P>&,const VectorFunction<P>&)) &compose);

}

Void export_scalar_functions() {
    export_scalar_function<EffectiveTag>();
    export_scalar_function<ValidatedTag>();
    implicitly_convertible<ScalarFunction<EffectiveTag>,ScalarFunction<ValidatedTag>>();
    def("lie_derivative", (ScalarFunction<EffectiveTag>(*)(const ScalarFunction<EffectiveTag>&,const VectorFunction<EffectiveTag>&)) &lie_derivative);
};

Void export_vector_functions() {
    export_vector_function<EffectiveTag>();
    export_vector_function<ValidatedTag>();
    implicitly_convertible<VectorFunction<EffectiveTag>,VectorFunction<ValidatedTag>>();
    from_python<VectorFunction<EffectiveTag>>();
};


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

    export_univariate_function();
    export_scalar_functions();
    export_vector_functions();

    //export_scalar_python_function();
    //export_vector_python_function();
}

