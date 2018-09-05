/***************************************************************************
 *            function_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "utility/array.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/multi_index.hpp"
#include "function/taylor_model.hpp"
#include "algebra/differential.hpp"
#include "function/formula.hpp"
#include "function/polynomial.hpp"
#include "function/affine.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/space.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/function_expression.hpp"
#include "solvers/constraint_solver.hpp"

#include "function/function_mixin.tpl.hpp"

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

typedef ScalarFunction<ApproximateTag> ASF;
ASF pow(ASF const&, Int);
ASF pos(ASF const&); ASF neg(ASF const&); ASF sqr(ASF const&); ASF rec(ASF const&);
ASF sqrt(ASF const&); ASF exp(ASF const&); ASF log(ASF const&); ASF atan(ASF const&);
ASF sin(ASF const&); ASF cos(ASF const&); ASF tan(ASF const&);

/*
RealExpression evaluate(EffectiveScalarFunction const& f, Vector<RealVariable> const& vars);
FloatDPApproximation evaluate(const ScalarFunction<ValidatedTag>&, const Vector<FloatDPApproximation>&);
FloatDPApproximation evaluate(const ScalarFunction<EffectiveTag>&, const Vector<FloatDPApproximation>&);
FloatDPBounds evaluate(const ScalarFunction<ValidatedTag>&, const Vector<FloatDPBounds>&);
FloatDPBounds evaluate(const ScalarFunction<EffectiveTag>&, const Vector<FloatDPBounds>&);
Vector<FloatDPApproximation> evaluate(const VectorFunction<ValidatedTag>&, const Vector<FloatDPApproximation>&);
Vector<FloatDPApproximation> evaluate(const VectorFunction<EffectiveTag>&, const Vector<FloatDPApproximation>&);
Vector<FloatDPBounds> evaluate(const VectorFunction<ValidatedTag>&, const Vector<FloatDPBounds>&);
Vector<FloatDPBounds> evaluate(const VectorFunction<EffectiveTag>&, const Vector<FloatDPBounds>&);
*/

/*
template<>
struct from_python< MultiIndex >
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<MultiIndex>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::tuple tup=boost::python::extract<boost::python::tuple>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage<MultiIndex>*)   data)->storage.bytes;
        MultiIndex res(static_cast<SizeType>(len(tup)));
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
        assert(len(lst)!=0);
        EffectiveScalarFunction sf0 = boost::python::extract<EffectiveScalarFunction>(lst[0]);
        EffectiveVectorFunction res(static_cast<SizeType>(len(lst)),sf0.argument_size());
        for(Nat i=0; i!=res.result_size(); ++i) { res.set(i,boost::python::extract<EffectiveScalarFunction>(lst[i])); }
        new (storage) EffectiveVectorFunction(res);
        data->convertible = storage;
    }
};
*/



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

/*
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
          _result_size(boost::python::extract<Nat>(pyf.attr("result_size"))),
          _argument_size(boost::python::extract<Nat>(pyf.attr("argument_size"))),
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
*/


}

using namespace Ariadne;

typedef FloatDPApproximation F;
typedef ExactIntervalType I;
typedef Vector<FloatDPApproximation> FV;
typedef Vector<ExactIntervalType> IV;
typedef Matrix<FloatDPApproximation> FMx;
typedef Matrix<ExactIntervalType> IMx;
typedef Vector< Differential<FloatDPApproximation> > FSDV;
typedef Vector< Differential<ExactIntervalType> > ISDV;
typedef Vector<ValidatedTaylorModelDP> TMV;
typedef ValidatedVectorTaylorFunctionModelDP TFM;
typedef ValidatedTaylorModelDP TM;

template<class X> using Monomial = ExpansionValue<MultiIndex,X>;

static constexpr auto self = pybind11::detail::self;



Void export_multi_index(pybind11::module& module)
{
    pybind11::class_< MultiIndex > multi_index_class(module,"MultiIndex");
    multi_index_class.def(pybind11::init<Nat>());
    multi_index_class.def(pybind11::init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def("__str__", &__cstr__<MultiIndex>);
//    multi_index_class.def("__repr__", &__repr__<MultiIndex>);

//    from_python<MultiIndex>();
//    to_python< List<MultiIndex> >();
}

template<class X>
Void export_monomial(pybind11::module& module)
{
    pybind11::class_< ExpansionValue<MultiIndex,X> > monomial_class(python_name<X>("Monomial"));
    monomial_class.def(pybind11::init<MultiIndex,X>());
    monomial_class.def("key",(const MultiIndex&(ExpansionValue<MultiIndex,X>::*)()const)&ExpansionValue<MultiIndex,X>::key);
    monomial_class.def("data",(const X&(ExpansionValue<MultiIndex,X>::*)()const) &ExpansionValue<MultiIndex,X>::data);
    monomial_class.def("__str__", &__cstr__< ExpansionValue<MultiIndex,X> >);
}

template<class X>
Void export_polynomial(pybind11::module& module)
{
    X real;

    pybind11::class_< Polynomial<X> > polynomial_class(python_name<X>("Polynomial"));
    polynomial_class.def(pybind11::init< Polynomial<X> >());
    polynomial_class.def(pybind11::init<Nat>());
    polynomial_class.def_static("constant", (Polynomial<X>(*)(Nat,double)) &Polynomial<X>::constant);
    polynomial_class.def_static("variable", (Polynomial<X>(*)(Nat,Nat)) &Polynomial<X>::variable);
    polynomial_class.def_static("coordinate", (Polynomial<X>(*)(Nat,Nat)) &Polynomial<X>::variable);
    polynomial_class.def_static("variables", (Vector< Polynomial<X> >(*)(Nat)) &Polynomial<X>::variables);

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
    polynomial_class.def("__str__",&__cstr__<Polynomial<X>>);
    //polynomial_class.def(self_ns::repr(self));

    to_python< Vector< Polynomial<X> > >();
}


Void export_univariate_function(pybind11::module& module)
{
    pybind11::class_<EffectiveScalarUnivariateFunction> function_class(module,"EffectiveScalarUnivariateFunction");
    function_class.def(pybind11::init<EffectiveScalarUnivariateFunction>());
    function_class.def("__call__", (FloatDPBounds(EffectiveScalarUnivariateFunction::*)(const FloatDPBounds&)const)&EffectiveScalarUnivariateFunction::operator() );
    function_class.def("__call__", (Differential<FloatDPBounds>(EffectiveScalarUnivariateFunction::*)(const Differential<FloatDPBounds>&)const)&EffectiveScalarUnivariateFunction::operator() );

    function_class.def_static("constant", (EffectiveScalarUnivariateFunction(*)(IntervalDomainType,EffectiveNumber)) &EffectiveScalarUnivariateFunction::constant);
    function_class.def_static("coordinate", (EffectiveScalarUnivariateFunction(*)()) &EffectiveScalarUnivariateFunction::coordinate);

}

template<class C, class F> void def(pybind11::class_<C>& c, const char* n, F const& f) {
    c.def(n,f); }
template<class A, class F> void def_call(pybind11::class_<F>& f) {
    typedef ResultOf<F(A)> R; f.def("__call__", (R(F::*)(A const&)const) &F::operator()); }
template<class A, class F> void def_evaluate(pybind11::class_<F>& f) {
    typedef ResultOf<F(A)> R; f.def("__call__", (R(*)(F const&, A const&)) &evaluate); }
template<class A, class F> void def_gradient(pybind11::class_<F>& f) {
    typedef decltype(declval<F>().gradient(declval<A>())) R; f.def("gradient", (R(F::*)(A const&)const) &F::gradient); }
template<class A, class F> void def_differential(pybind11::class_<F>& f) {
    typedef DegreeType D; typedef decltype(declval<F>().differential(declval<A>(),declval<D>())) R;
    f.def("differential", (R(F::*)(A const&,D)const) &F::differential); }

template<class F, class T> using ArgumentType = typename F::template Argument<T>;

template<class F> Void export_function_evaluation(pybind11::class_<F>& function_class, ApproximateTag) {
    def_call<ArgumentType<F,FloatDPApproximation>>(function_class);
    def_call<ArgumentType<F,FloatMPApproximation>>(function_class);
    def_call<ArgumentType<F,Differential<FloatDPApproximation>>>(function_class);
    def_call<ArgumentType<F,Differential<FloatMPApproximation>>>(function_class);

    def_evaluate<ArgumentType<F,FloatDPApproximation>>(function_class);
    def_evaluate<ArgumentType<F,FloatMPApproximation>>(function_class);

    def_differential<ArgumentType<F,FloatDPApproximation>>(function_class);
    def_differential<ArgumentType<F,FloatMPApproximation>>(function_class);
}

template<class F> Void export_function_evaluation(pybind11::class_<F>& function_class, ValidatedTag) {
    export_function_evaluation(function_class,ApproximateTag());
    def_call<ArgumentType<F,FloatDPBounds>>(function_class);
    def_call<ArgumentType<F,FloatMPBounds>>(function_class);
    def_call<ArgumentType<F,Differential<FloatDPBounds>>>(function_class);
    def_call<ArgumentType<F,Differential<FloatMPBounds>>>(function_class);
    def_evaluate<ArgumentType<F,FloatDPBounds>>(function_class);
    def_evaluate<ArgumentType<F,FloatMPBounds>>(function_class);
    def_differential<ArgumentType<F,FloatDPBounds>>(function_class);
    def_differential<ArgumentType<F,FloatMPBounds>>(function_class);
}

template<class F> Void export_function_evaluation(pybind11::class_<F>& function_class, EffectiveTag) {
    export_function_evaluation(function_class,ValidatedTag());
}

template<class F> Void export_function_evaluation(pybind11::class_<F>& function_class)
{
    export_function_evaluation(function_class, Paradigm<F>());
}


/*
    function_class.def("gradient", (Covector<FloatDPApproximation>(ScalarFunction<P>::*)(const Argument<FloatDPApproximation>&)const)&ScalarFunction<P>::gradient );
    function_class.def("gradient", (Covector<FloatDPBounds>(ScalarFunction<P>::*)(const Argument<FloatDPBounds>&)const)&ScalarFunction<P>::gradient );
    function_class.def("differential", (Differential<FloatDPApproximation>(ScalarFunction<P>::*)(const Argument<FloatDPApproximation>&,DegreeType)const) &ScalarFunction<P>::differential);
    function_class.def("differential", (Differential<FloatDPBounds>(ScalarFunction<P>::*)(const Argument<FloatDPBounds>&,DegreeType)const) &ScalarFunction<P>::differential);
*/


template<class P> Void export_scalar_function_evaluation(pybind11::class_<ScalarFunction<P>>& scalar_function_class) {
    using FP=ScalarFunction<P>;
    export_function_evaluation(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatDPApproximation>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatMPApproximation>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatDPBounds>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatMPBounds>>(scalar_function_class);
}

template<> Void export_scalar_function_evaluation<ApproximateTag>(pybind11::class_<ScalarFunction<ApproximateTag>>& scalar_function_class) {
    using P=ApproximateTag;
    using FP=ScalarFunction<P>;
    export_function_evaluation(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatDPApproximation>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatMPApproximation>>(scalar_function_class);
}

Void export_function_evaluation(pybind11::class_<ScalarFunction<ApproximateTag>>& scalar_function_class)
{
    def_call<Vector<FloatDPApproximation>>(scalar_function_class);
    def_call<Vector<FloatMPApproximation>>(scalar_function_class);
    def_call<Vector<Differential<FloatDPApproximation>>>(scalar_function_class);
    def_call<Vector<Differential<FloatMPApproximation>>>(scalar_function_class);

    def_evaluate<Vector<FloatDPApproximation>>(scalar_function_class);
    def_evaluate<Vector<FloatMPApproximation>>(scalar_function_class);
}

template<class P> Void export_vector_function_evaluation(pybind11::class_<VectorFunction<P>>& vector_function_class){
    export_function_evaluation(vector_function_class);
}

template<class P> Void export_scalar_function(pybind11::module& module)
{
    pybind11::class_<ScalarFunction<P>> scalar_function_class(module,(class_name<P>()+"ScalarFunction").c_str());
    scalar_function_class.def(pybind11::init<ScalarFunction<P>>());
    scalar_function_class.def(pybind11::init<SizeType>());
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

    scalar_function_class.def("__str__", &__cstr__<ScalarFunction<P>>);
    scalar_function_class.def("__repr__", &__crepr__<ScalarFunction<P>>);

    scalar_function_class.def_static("constant", (ScalarFunction<P>(*)(SizeType,Number<P>)) &ScalarFunction<P>::constant);
    scalar_function_class.def_static("coordinate", (ScalarFunction<P>(*)(SizeType,SizeType)) &ScalarFunction<P>::coordinate);

    module.def("pow", (ScalarFunction<P>(*)(const ScalarFunction<P>&,Int)) &pow);
    module.def("rec", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &rec);
    module.def("sqr", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &sqr);
    module.def("sqrt", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &sqrt);
    module.def("exp", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &exp);
    module.def("log", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &log);
    module.def("sin", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &sin);
    module.def("cos", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &cos);
    module.def("tan", (ScalarFunction<P>(*)(const ScalarFunction<P>&)) &tan);

    module.def("derivative", (ScalarFunction<P>(ScalarFunction<P>::*)(SizeType)const) &ScalarFunction<P>::derivative);

    export_scalar_function_evaluation(scalar_function_class);
}


template<class P> Void export_vector_function(pybind11::module& module)
{
    pybind11::class_<VectorFunction<P>> vector_function_class(module,(class_name<P>()+"VectorFunction").c_str());
    vector_function_class.def(pybind11::init<VectorFunction<P>>());
    vector_function_class.def(pybind11::init<Nat,Nat>());

    vector_function_class.def("result_size", &VectorFunction<P>::result_size);
    vector_function_class.def("argument_size", &VectorFunction<P>::argument_size);
    vector_function_class.def("__getitem__", &VectorFunction<P>::get);
    vector_function_class.def("__setitem__", &VectorFunction<P>::set);
    export_vector_function_evaluation(vector_function_class);
//    vector_function_class.def("__call__", (Vector<FloatDPBounds>(VectorFunction<P>::*)(const Vector<FloatDPBounds>&)const)&VectorFunction<P>::operator() );
//    vector_function_class.def("__call__", (Vector<FloatDPApproximation>(VectorFunction<P>::*)(const Vector<FloatDPApproximation>&)const)&VectorFunction<P>::operator() );
//    vector_function_class.def("__call__", (Vector<Differential<FloatDPBounds>>(VectorFunction<P>::*)(const Vector<Differential<FloatDPBounds>>&)const)&VectorFunction<P>::evaluate );
//    vector_function_class.def("__call__", (Vector<Differential<FloatDPApproximation>>(VectorFunction<P>::*)(const Vector<Differential<FloatDPApproximation>>&)const)&VectorFunction<P>::evaluate );
//    vector_function_class.def("jacobian", (Matrix<FloatDPBounds>(VectorFunction<P>::*)(const Vector<FloatDPBounds>&)const) &VectorFunction<P>::jacobian);
//    vector_function_class.def("jacobian", (Matrix<FloatDPApproximation>(VectorFunction<P>::*)(const Vector<FloatDPApproximation>&)const) &VectorFunction<P>::jacobian);
//    vector_function_class.def("differential", (Vector<Differential<FloatDPBounds> >(VectorFunction<P>::*)(const Vector<FloatDPBounds>&,DegreeType)const) &VectorFunction<P>::differential);
//    vector_function_class.def("differential", (Vector<Differential<FloatDPApproximation> >(VectorFunction<P>::*)(const Vector<FloatDPApproximation>&,DegreeType)const) &VectorFunction<P>::differential);
    vector_function_class.def("__str__", &__cstr__<VectorFunction<P>>);
    vector_function_class.def("__repr__", &__crepr__<VectorFunction<P>>);

    vector_function_class.def_static("identity", (VectorFunction<P>(*)(SizeType)) &VectorFunction<P>::identity);

    module.def("evaluate", (Vector<FloatDPApproximation>(*)(const VectorFunction<P>&,const Vector<FloatDPApproximation>&)) &evaluate);
//    module.def("evaluate", (Vector<FloatDPBounds>(*)(const VectorFunction<P>&,const Vector<FloatDPBounds>&)) &evaluate);

    module.def("join", (VectorFunction<P>(*)(const ScalarFunction<P>&, const ScalarFunction<P>&)) &join);
    module.def("join", (VectorFunction<P>(*)(const VectorFunction<P>&, const ScalarFunction<P>&)) &join);
    module.def("join", (VectorFunction<P>(*)(const ScalarFunction<P>&, const VectorFunction<P>&)) &join);
    module.def("join", (VectorFunction<P>(*)(const VectorFunction<P>&, const VectorFunction<P>&)) &join);

    module.def("compose", (ScalarFunction<P>(*)(const ScalarFunction<P>&,const VectorFunction<P>&)) &compose);
    module.def("compose", (VectorFunction<P>(*)(const VectorFunction<P>&,const VectorFunction<P>&)) &compose);

    export_vector_function_evaluation(vector_function_class);
}

template<class Y, class X> Void export_procedure(pybind11::module& module) {
    typedef Paradigm<Y> P;
    pybind11::class_<Procedure<Y>> procedure_class(module,(class_name<P>()+"Procedure").c_str());
    procedure_class.def("__str__", &__cstr__<Procedure<Y>>);
    module.def("make_procedure", (Procedure<Y>(*)(ScalarFunction<P> const&)) &make_procedure);
    module.def("evaluate", (X(*)(Procedure<Y> const&, Vector<X> const&)) &evaluate);
    module.def("gradient", (Covector<X>(*)(Procedure<Y> const&, Vector<X> const&)) &gradient);
    module.def("hessian", (X(*)(Procedure<Y> const&, Vector<X> const&, Vector<X> const&)) &hessian);
}

Void export_scalar_functions(pybind11::module& module) {
    export_scalar_function<ApproximateTag>(module);
    export_scalar_function<ValidatedTag>(module);
    export_scalar_function<EffectiveTag>(module);
    pybind11::implicitly_convertible<ScalarFunction<EffectiveTag>,ScalarFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<ScalarFunction<EffectiveTag>,ScalarFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<ScalarFunction<ValidatedTag>,ScalarFunction<ApproximateTag>>();
    module.def("lie_derivative", (ScalarFunction<EffectiveTag>(*)(const ScalarFunction<EffectiveTag>&,const VectorFunction<EffectiveTag>&)) &lie_derivative);
}

Void export_vector_functions(pybind11::module& module) {
    export_vector_function<ApproximateTag>(module);
    export_vector_function<ValidatedTag>(module);
    export_vector_function<EffectiveTag>(module);
    pybind11::implicitly_convertible<VectorFunction<EffectiveTag>,VectorFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<VectorFunction<EffectiveTag>,VectorFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<VectorFunction<ValidatedTag>,VectorFunction<ApproximateTag>>();
//    from_python<VectorFunction<EffectiveTag>>();
}


/*
Void export_scalar_python_function(pybind11::module& module)
{
    pybind11::class_<ScalarPythonFunction, bases< EffectiveScalarFunctionInterface > > scalar_python_function_class(module,"ScalarUserFunction", pybind11::init<object>());
    scalar_python_function_class.def(pybind11::init<Nat,object>());
}

Void export_vector_python_function(pybind11::module& module)
{
    pybind11::class_<VectorPythonFunction, bases< EffectiveVectorFunctionInterface > > vector_python_function_class(module,"VectorUserFunction", pybind11::init<object>());
    vector_python_function_class.def(pybind11::init<Nat,Nat,object>());
}
*/


Void function_submodule(pybind11::module& module) {
//    to_python< Array<StringType> >();
//    from_python< Array<StringType> >();

    export_multi_index(module);

    export_univariate_function(module);
    export_scalar_functions(module);
    export_vector_functions(module);

    export_procedure<ApproximateNumber, FloatDPApproximation>(module);
    export_procedure<ValidatedNumber, FloatDPBounds>(module);


    //export_scalar_python_function(module);
    //export_vector_python_function(module);
}

