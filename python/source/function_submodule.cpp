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




MultiIndex multi_index_from_python(pybind11::tuple pytup) {
    MultiIndex res(static_cast<SizeType>(len(pytup)));
    for(SizeType i=0; i!=res.size(); ++i) { res.set(i,pybind11::cast<Nat>(pytup[i])); }
    return res;
}

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

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<MultiIndex>& arepr) {
    MultiIndex const& a=arepr.reference(); os << "("; for(SizeType i=0; i!=a.size(); ++i) { if(i!=0) { os << ','; } os << a[i]; } os << ")"; return os;
}

template<class P, class D, class C> OutputStream& operator<<(OutputStream& os, const Representation< Function<P,D,C> >& frepr) {
    static_cast<const FunctionInterface<P,D,C>&>(frepr.reference()).repr(os); return os;
}

}// namespace Ariadne

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
typedef ValidatedVectorMultivariateTaylorFunctionModelDP TFM;
typedef ValidatedTaylorModelDP TM;



Void export_multi_index(pybind11::module& module)
{
    pybind11::class_< MultiIndex > multi_index_class(module,"MultiIndex");
    multi_index_class.def(pybind11::init<Nat>());
    multi_index_class.def(pybind11::init(&multi_index_from_python));
    multi_index_class.def(pybind11::init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def("__str__", &__cstr__<MultiIndex>);
    multi_index_class.def("__repr__", &__repr__<MultiIndex>);

    pybind11::implicitly_convertible<pybind11::tuple,MultiIndex>();
}


template<class X>
Void export_polynomial(pybind11::module& module)
{
    pybind11::class_< Monomial<X> > monomial_class(module,python_name<X>("Monomial").c_str());
    monomial_class.def(pybind11::init<MultiIndex,X>());
    monomial_class.def("index", (MultiIndex const&(Monomial<X>::*)()const) &Monomial<X>::index);
    monomial_class.def("coefficient", (X const&(Monomial<X>::*)()const) &Monomial<X>::coefficient);
    monomial_class.def("__str__", &__cstr__<Monomial<X>>);


    pybind11::class_< Polynomial<X> > polynomial_class(module,python_name<X>("Polynomial").c_str());
    polynomial_class.def(pybind11::init< Polynomial<X> >());
    polynomial_class.def(pybind11::init<Nat>());
    polynomial_class.def_static("constant", (Polynomial<X>(*)(SizeType,X const&)) &Polynomial<X>::constant);
    polynomial_class.def_static("variable", (Polynomial<X>(*)(SizeType,SizeType)) &Polynomial<X>::variable);
    polynomial_class.def_static("coordinate", (Polynomial<X>(*)(SizeType,SizeType)) &Polynomial<X>::variable);

    polynomial_class.def_static("variables", [](Nat as){return Polynomial<X>::variables(as).array();});

    polynomial_class.def("argument_size", &Polynomial<X>::argument_size);
    polynomial_class.def("insert", &Polynomial<X>::insert);

    define_algebra(module,polynomial_class);
    polynomial_class.def("__str__",&__cstr__<Polynomial<X>>);

    export_vector<Polynomial<X>>(module, (python_name<X>("PolynomialVector")).c_str());
}

Void export_polynomials(pybind11::module& module)
{
    export_polynomial<FloatDPBounds>(module);
    export_polynomial<FloatDPApproximation>(module);
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



template<class P> Void export_scalar_function_evaluation(pybind11::class_<ScalarMultivariateFunction<P>>& scalar_function_class) {
    using FP=ScalarMultivariateFunction<P>;
    export_function_evaluation(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatDPApproximation>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatMPApproximation>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatDPBounds>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatMPBounds>>(scalar_function_class);
}

template<> Void export_scalar_function_evaluation<ApproximateTag>(pybind11::class_<ScalarMultivariateFunction<ApproximateTag>>& scalar_function_class) {
    using P=ApproximateTag;
    using FP=ScalarMultivariateFunction<P>;
    export_function_evaluation(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatDPApproximation>>(scalar_function_class);
    def_gradient<ArgumentType<FP,FloatMPApproximation>>(scalar_function_class);
}

Void export_function_evaluation(pybind11::class_<ScalarMultivariateFunction<ApproximateTag>>& scalar_function_class)
{
    def_call<Vector<FloatDPApproximation>>(scalar_function_class);
    def_call<Vector<FloatMPApproximation>>(scalar_function_class);
    def_call<Vector<Differential<FloatDPApproximation>>>(scalar_function_class);
    def_call<Vector<Differential<FloatMPApproximation>>>(scalar_function_class);

    def_evaluate<Vector<FloatDPApproximation>>(scalar_function_class);
    def_evaluate<Vector<FloatMPApproximation>>(scalar_function_class);
}

template<class P> Void export_vector_function_evaluation(pybind11::class_<VectorMultivariateFunction<P>>& vector_function_class){
    export_function_evaluation(vector_function_class);
}

template<class P> Void export_scalar_function(pybind11::module& module)
{
    pybind11::class_<ScalarMultivariateFunction<P>> scalar_function_class(module,(class_name<P>()+"ScalarMultivariateFunction").c_str());
    scalar_function_class.def(pybind11::init<ScalarMultivariateFunction<P>>());
    scalar_function_class.def(pybind11::init<SizeType>());
    scalar_function_class.def("argument_size", &ScalarMultivariateFunction<P>::argument_size);
    scalar_function_class.def("derivative", &ScalarMultivariateFunction<P>::derivative);
    define_elementary_algebra<ScalarMultivariateFunction<P>,Number<P>>(module,scalar_function_class);

    if constexpr (IsSame<P,ValidatedTag>::value) {
        scalar_function_class.def(pybind11::init<ScalarMultivariateFunction<EffectiveTag>>());
        pybind11::implicitly_convertible<ScalarMultivariateFunction<EffectiveTag>,ScalarMultivariateFunction<ValidatedTag>>();
    }
    
    
//FIXME
//    scalar_function_class.def("__eq__", &__eq__<Constraint<ScalarMultivariateFunction<P>,Number<P>>,ScalarMultivariateFunction<P>,Number<P>>);
//    scalar_function_class.def("__le__", &__le__<Constraint<ScalarMultivariateFunction<P>,Number<P>>,ScalarMultivariateFunction<P>,Number<P>>);
//    scalar_function_class.def("__ge__", &__ge__<Constraint<ScalarMultivariateFunction<P>,Number<P>>,ScalarMultivariateFunction<P>,Number<P>>);

    scalar_function_class.def("__str__", &__cstr__<ScalarMultivariateFunction<P>>);
    scalar_function_class.def("__repr__", &__crepr__<ScalarMultivariateFunction<P>>);

    scalar_function_class.def_static("constant", (ScalarMultivariateFunction<P>(*)(SizeType,Number<P>)) &ScalarMultivariateFunction<P>::constant);
    scalar_function_class.def_static("coordinate", (ScalarMultivariateFunction<P>(*)(SizeType,SizeType)) &ScalarMultivariateFunction<P>::coordinate);

    scalar_function_class.def("gradient", (Covector<FloatDPApproximation>(ScalarMultivariateFunction<P>::*)(const Vector<FloatDPApproximation>&)const) &ScalarMultivariateFunction<P>::gradient);
    if constexpr (not IsSame<P,ApproximateTag>::value) {
        scalar_function_class.def("gradient", (Covector<FloatDPBounds>(ScalarMultivariateFunction<P>::*)(const Vector<FloatDPBounds>&)const) &ScalarMultivariateFunction<P>::gradient);
    }

    module.def("derivative", (ScalarMultivariateFunction<P>(ScalarMultivariateFunction<P>::*)(SizeType)const) &ScalarMultivariateFunction<P>::derivative);

    module.def("evaluate", (Scalar<FloatDPApproximation>(*)(const ScalarMultivariateFunction<P>&,const Vector<FloatDPApproximation>&)) &evaluate);
    if constexpr (not IsSame<P,ApproximateTag>::value) {
        module.def("evaluate", (Scalar<FloatDPBounds>(*)(const ScalarMultivariateFunction<P>&,const Vector<FloatDPBounds>&)) &evaluate);
    }

    export_scalar_function_evaluation(scalar_function_class);
}


template<class P> Void export_vector_function(pybind11::module& module)
{
    pybind11::class_<VectorMultivariateFunction<P>> vector_function_class(module,(class_name<P>()+"VectorMultivariateFunction").c_str());
    vector_function_class.def(pybind11::init<VectorMultivariateFunction<P>>());
    vector_function_class.def(pybind11::init<Nat,Nat>());
    if constexpr (IsSame<P,ValidatedTag>::value) {
        vector_function_class.def(pybind11::init<VectorMultivariateFunction<EffectiveTag>>());
//        pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ValidatedTag>>();
    }
    // NOTE: This must go *after* the conversion constructor
    vector_function_class.def(pybind11::init([](std::vector<ScalarMultivariateFunction<P>> const& lst){return VectorMultivariateFunction<P>(lst);}));

    vector_function_class.def("result_size", &VectorMultivariateFunction<P>::result_size);
    vector_function_class.def("argument_size", &VectorMultivariateFunction<P>::argument_size);
    vector_function_class.def("__getitem__", &VectorMultivariateFunction<P>::get);
    vector_function_class.def("__setitem__", &VectorMultivariateFunction<P>::set);

    // TODO: Put these in C++ API
    // define_vector_algebra_arithmetic<VectorMultivariateFunction<P>,ScalarMultivariateFunction<P>,Number<P>>(module,vector_function_class);

    // FIXME: Define vector function operations for Validated and Approximate
    if constexpr (IsSame<P,EffectiveTag>::value) {
        define_vector_algebra_arithmetic<VectorMultivariateFunction<P>,ScalarMultivariateFunction<P>>(module,vector_function_class);
    }
    export_vector_function_evaluation(vector_function_class);

    vector_function_class.def("jacobian", (Matrix<FloatDPApproximation>(VectorMultivariateFunction<P>::*)(const Vector<FloatDPApproximation>&)const) &VectorMultivariateFunction<P>::jacobian);
    if constexpr (not IsSame<P,ApproximateTag>::value) {
        vector_function_class.def("jacobian", (Matrix<FloatDPBounds>(VectorMultivariateFunction<P>::*)(const Vector<FloatDPBounds>&)const) &VectorMultivariateFunction<P>::jacobian);
    }

    vector_function_class.def("__str__", &__cstr__<VectorMultivariateFunction<P>>);
    vector_function_class.def("__repr__", &__crepr__<VectorMultivariateFunction<P>>);

    vector_function_class.def_static("identity", (VectorMultivariateFunction<P>(*)(SizeType)) &VectorMultivariateFunction<P>::identity);

    module.def("evaluate", (Vector<FloatDPApproximation>(*)(const VectorMultivariateFunction<P>&,const Vector<FloatDPApproximation>&)) &evaluate);
    if constexpr (not IsSame<P,ApproximateTag>::value) {
        module.def("evaluate", (Vector<FloatDPBounds>(*)(const VectorMultivariateFunction<P>&,const Vector<FloatDPBounds>&)) &evaluate);
    }

    module.def("join", (VectorMultivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&, const ScalarMultivariateFunction<P>&)) &join);
    module.def("join", (VectorMultivariateFunction<P>(*)(const VectorMultivariateFunction<P>&, const ScalarMultivariateFunction<P>&)) &join);
    module.def("join", (VectorMultivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&, const VectorMultivariateFunction<P>&)) &join);
    module.def("join", (VectorMultivariateFunction<P>(*)(const VectorMultivariateFunction<P>&, const VectorMultivariateFunction<P>&)) &join);

    module.def("compose", (ScalarMultivariateFunction<P>(*)(const ScalarMultivariateFunction<P>&,const VectorMultivariateFunction<P>&)) &compose);
    module.def("compose", (VectorMultivariateFunction<P>(*)(const VectorMultivariateFunction<P>&,const VectorMultivariateFunction<P>&)) &compose);

    export_vector_function_evaluation(vector_function_class);
}

template<class Y, class X> Void export_procedure(pybind11::module& module) {
    typedef Paradigm<Y> P;
    pybind11::class_<Procedure<Y>> procedure_class(module,(class_name<P>()+"Procedure").c_str());
    procedure_class.def("__str__", &__cstr__<Procedure<Y>>);
    module.def("make_procedure", (Procedure<Y>(*)(ScalarMultivariateFunction<P> const&)) &make_procedure);
    module.def("evaluate", (X(*)(Procedure<Y> const&, Vector<X> const&)) &evaluate);
    module.def("gradient", (Covector<X>(*)(Procedure<Y> const&, Vector<X> const&)) &gradient);
    module.def("hessian", (X(*)(Procedure<Y> const&, Vector<X> const&, Vector<X> const&)) &hessian);
}

Void export_scalar_functions(pybind11::module& module) {
    export_scalar_function<ApproximateTag>(module);
    export_scalar_function<ValidatedTag>(module);
    export_scalar_function<EffectiveTag>(module);
    pybind11::implicitly_convertible<ScalarMultivariateFunction<EffectiveTag>,ScalarMultivariateFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<ScalarMultivariateFunction<EffectiveTag>,ScalarMultivariateFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<ScalarMultivariateFunction<ValidatedTag>,ScalarMultivariateFunction<ApproximateTag>>();
    module.def("lie_derivative", (ScalarMultivariateFunction<EffectiveTag>(*)(const ScalarMultivariateFunction<EffectiveTag>&,const VectorMultivariateFunction<EffectiveTag>&)) &lie_derivative);
}

Void export_vector_functions(pybind11::module& module) {
    export_vector_function<ApproximateTag>(module);
    export_vector_function<ValidatedTag>(module);
    export_vector_function<EffectiveTag>(module);
    pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ValidatedTag>>();
    pybind11::implicitly_convertible<VectorMultivariateFunction<EffectiveTag>,VectorMultivariateFunction<ApproximateTag>>();
    pybind11::implicitly_convertible<VectorMultivariateFunction<ValidatedTag>,VectorMultivariateFunction<ApproximateTag>>();
}




Void function_submodule(pybind11::module& module) {

    export_multi_index(module);

    export_polynomials(module);

    export_univariate_function(module);
    export_scalar_functions(module);
    export_vector_functions(module);

    export_procedure<ApproximateNumber, FloatDPApproximation>(module);
    export_procedure<ValidatedNumber, FloatDPBounds>(module);
}

