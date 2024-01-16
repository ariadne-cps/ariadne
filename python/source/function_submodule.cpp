/***************************************************************************
 *            function_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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
#include "numeric_submodule.hpp"

#include <iostream>
#include <iomanip>

#include "utility/array.hpp"
#include "utility/container.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/symmetric_matrix.hpp"
#include "algebra/expansion.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/differential.hpp"
#include "function/formula.hpp"
#include "function/polynomial.hpp"
#include "function/affine.hpp"
#include "function/taylor_model.hpp"
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

template<> struct PythonTemplateName<MultivariatePolynomial> { static std::string get() { return "MultivariatePolynomial"; } };
template<> struct PythonTemplateName<Procedure> { static std::string get() { return "Procedure"; } };
template<> struct PythonTemplateName<Function> { static std::string get() { return "Function"; } };

template<class X> struct PythonClassName<MultivariatePolynomial<X>> { static std::string get() { return python_template_class_name<X>("MultivariatePolynomial"); } };

template struct PythonClassName<MultivariatePolynomial<FloatDPApproximation>>;
template struct PythonClassName<MultivariatePolynomial<FloatMPApproximation>>;
template struct PythonClassName<MultivariatePolynomial<FloatDPBounds>>;
template struct PythonClassName<MultivariatePolynomial<FloatMPBounds>>;


template<class K, class V> pybind11::dict to_python_dict(Map<K,V> const& map) {
    pybind11::dict pydict;
    for (auto item : map) {
        pydict[item.first]=item.second;
    }
    return pydict;
}

template<class K, class V> Map<K,V> map_from_python(pybind11::dict pydict) {
    Map<K,V> r;
    for (auto item : pydict) {
        auto k = pybind11::cast<K>(item.first);
        auto v = pybind11::cast<V>(item.second);
        r.insert(std::make_pair(k,v));
    }
    return r;
}

template Map<int,double> map_from_python(pybind11::dict);

MultiIndex multi_index_from_python_tuple(pybind11::tuple pytup) {
    MultiIndex res(static_cast<SizeType>(len(pytup)));
    for(SizeType i=0; i!=res.size(); ++i) { res.set(i,pybind11::cast<Nat>(pytup[i])); }
    return res;
}

MultiIndex multi_index_from_python_list(pybind11::list pylst) {
    MultiIndex res(static_cast<SizeType>(len(pylst)));
    for(SizeType i=0; i!=res.size(); ++i) { res.set(i,pybind11::cast<Nat>(pylst[i])); }
    return res;
}

pybind11::tuple multi_index_to_python_tuple(MultiIndex const& a) {
    pybind11::list pylst;
    for(SizeType i=0; i!=a.size(); ++i) { pylst.append(a[i]); }
    return pybind11::tuple(pylst);
}

template<class I, class X> pybind11::dict expansion_to_python_dict(Expansion<I,X> const& e) {
    pybind11::dict pydict;
    for (auto t : e) {
        I a=t.index();
        X c=t.coefficient();
        if constexpr (Same<I,MultiIndex>) {
            pydict[multi_index_to_python_tuple(a)]=c;
        } else {
            pybind11::handle k=pybind11::cast(a);
            pydict[k]=c;
        }
    }
    return pydict;
}

template<class I, class X> Expansion<I,X> expansion_from_python(pybind11::dict dct) {
    auto item0 = *dct.begin();
    auto a0 = pybind11::cast<I>(item0.first);
    auto c0 = pybind11::cast<X>(item0.second);
    Expansion<I,X> e(a0.size(),nul(c0));
    for (auto item : dct) {
        auto a = pybind11::cast<I>(item.first);
        auto c = pybind11::cast<X>(item.second);
        e.append(a,c);
    }
    return e;
//    return Expansion<I,X>(map_from_python<I,X>(dct))
}

template<class I, class X, class... PRS> Expansion<I,X> expansion_from_python(pybind11::dict dct, PRS... prs) {
    using Y=typename X::GenericType;
    static_assert(Constructible<X,Y,PRS...>);
    auto item0 = *dct.begin();
    auto a0 = pybind11::cast<I>(item0.first);
    auto z0 = X(prs...);
    Expansion<I,X> e(a0.size(),z0);
    for (auto item : dct) {
        auto a = pybind11::cast<I>(item.first);
        auto c = pybind11::cast<Y>(item.second);
        e.append(a,X(c,prs...));
    }
    return e;
}

template<class P, class ARG> pybind11::object function_from_python(typename DomainTraits<ARG>::EntireDomainType dom, pybind11::object pyf) {
    Function<P,ARG(ARG)> id = Function<P,ARG(ARG)>::identity(dom);
    pybind11::object pyobjf = pyf(id);
    try {
        return pybind11::cast(pybind11::cast<VectorFunction<P,ARG>>(pyobjf));
    } catch (...) { }
    try {
        pybind11::list pylstf = pybind11::cast<pybind11::list>(pyobjf);
        SizeType rs = len(pylstf);
        List<ScalarFunction<P,ARG>> lsf;
        for (SizeType i=0; i!=rs; ++i) {
            lsf.append(pybind11::cast<ScalarFunction<P,ARG>>(pylstf[i]));
        }
        if (lsf.size()==0) {
            return pybind11::cast(VectorFunction<P,ARG>(0,dom));
        } else {
            return pybind11::cast(VectorFunction<P,ARG>(lsf));
        }
    } catch (...) { }
    return pybind11::cast(pybind11::cast<ScalarFunction<P,ARG>>(pyobjf));
}

pybind11::object univariate_function_from_python(pybind11::object pyf) {
    return function_from_python<EffectiveTag,RealScalar>(RealDomain(),pyf);
}
pybind11::object multivariate_function_from_python(SizeType as, pybind11::object pyf) {
    return function_from_python<EffectiveTag,RealVector>(EuclideanDomain(as),pyf);
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

template<class P, class SIG> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Function<P,SIG> >& frepr) {
    Function<P,SIG> const& f = frepr.reference();
    os << "Function(";
    if constexpr (Same<decltype(f.argument_size()),SizeOne>) { }
    else { os << frepr.reference().argument_size() << ","; }
    return os << " lambda x: " << frepr.reference() << " )";
}

}// namespace Ariadne

using namespace Ariadne;

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


template<class AT, class F> void define_call(pybind11::class_<F>& f) {
    typedef typename F::template Argument<AT> A; typedef ResultOf<F(A)> R; f.def("__call__", (R(F::*)(A const&)const) &F::operator()); }

template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<ScalarUnivariateFunction<P>>& function_class) {
    typedef ScalarUnivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(slope(declval<F>(),declval<A>())) R;
    module.def("slope", (R(*)(F const&,A const&)) &slope);
    function_class.def("slope", (R(F::*)(A const&)const) &F::slope);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<ScalarMultivariateFunction<P>>& function_class) {
    typedef ScalarMultivariateFunction<P> F; typedef typename F::template Argument<AT> A;
    typedef decltype(gradient(declval<F>(),declval<A>())) R;
    module.def("gradient", (R(*)(F const&,A const&)) &gradient);
    function_class.def("gradient", (R(F::*)(A const&)const) &F::gradient);
    typedef decltype(hessian(declval<F>(),declval<A>())) H;
    module.def("hessian", (H(*)(F const&,A const&)) &hessian);
    //function_class.def("hessian", (H(F::*)(A const&)const) &F::hessian);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<VectorUnivariateFunction<P>>& function_class) {
    typedef VectorUnivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(tangent(declval<F>(),declval<A>())) R;
    module.def("tangent", (R(*)(F const&,A const&)) &tangent);
    function_class.def("tangent", (R(F::*)(A const&)const) &F::tangent);
}
template<class AT, class P> void define_named_derivatives(pybind11::module& module, pybind11::class_<VectorMultivariateFunction<P>>& function_class) {
    typedef VectorMultivariateFunction<P> F; typedef typename F::template Argument<AT> A; typedef decltype(jacobian(declval<F>(),declval<A>())) R;
    module.def("jacobian", (R(*)(F const&,A const&)) &jacobian);
    function_class.def("jacobian", (R(F::*)(A const&)const) &F::jacobian);
}


template<class F, class T> using ArgumentType = typename F::template Argument<T>;

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, ApproximateTag) {
    define_call<FloatDPApproximation>(function_class);
    define_call<FloatMPApproximation>(function_class);
    define_call<Differential<FloatDPApproximation>>(function_class);
    define_call<Differential<FloatMPApproximation>>(function_class);

    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatDPApproximation>>);
    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatMPApproximation>>);
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, ValidatedTag) {
    define_call<FloatDPBounds>(function_class);
    define_call<FloatMPBounds>(function_class);
    define_call<Differential<FloatDPBounds>>(function_class);
    define_call<Differential<FloatMPBounds>>(function_class);

    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatDPBounds>>);
    module.def("evaluate", _evaluate_<F,ArgumentType<F,FloatMPBounds>>);

    export_function_evaluation(module,function_class,ApproximateTag());
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class, EffectiveTag) {
    export_function_evaluation(module,function_class,ValidatedTag());
}

template<class F> Void export_function_evaluation(pybind11::module module, pybind11::class_<F>& function_class)
{
    export_function_evaluation(module, function_class, Paradigm<F>());
}


template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class, ApproximateTag) {
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatDPApproximation>,DegreeType>);
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatMPApproximation>,DegreeType>);

    module.def("differential", _differential_<F,ArgumentType<F,FloatDPApproximation>,DegreeType>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatMPApproximation>,DegreeType>);

    define_named_derivatives<FloatDPApproximation>(module,function_class);
    define_named_derivatives<FloatMPApproximation>(module,function_class);
}

template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class, ValidatedTag) {
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatDPBounds>,DegreeType>);
    function_class.def("differential", _differential_<F,ArgumentType<F,FloatMPBounds>,DegreeType>);

    module.def("differential", _differential_<F,ArgumentType<F,FloatDPBounds>,DegreeType>);
    module.def("differential", _differential_<F,ArgumentType<F,FloatMPBounds>,DegreeType>);

    define_named_derivatives<FloatDPBounds>(module,function_class);
    define_named_derivatives<FloatMPBounds>(module,function_class);

    export_function_derivatives(module,function_class,ApproximateTag());
}

template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class, EffectiveTag) {
    export_function_derivatives(module,function_class,ValidatedTag());
}

template<class F> Void export_function_derivatives(pybind11::module module, pybind11::class_<F>& function_class)
{
    export_function_derivatives(module, function_class, Paradigm<F>());
}


Void export_multi_index(pybind11::module& module)
{
    pybind11::class_< MultiIndex > multi_index_class(module,"MultiIndex");
    multi_index_class.def(pybind11::init<Nat>());
    multi_index_class.def(pybind11::init(&multi_index_from_python_tuple));
    multi_index_class.def(pybind11::init(&multi_index_from_python_list));
    multi_index_class.def(pybind11::init<MultiIndex>());
    multi_index_class.def("__len__",&MultiIndex::size);
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def("__str__", &__cstr__<MultiIndex>);
    multi_index_class.def("__repr__", &__repr__<MultiIndex>);

    pybind11::implicitly_convertible<pybind11::tuple,MultiIndex>();
    //multi_index_class.def("tuple",&multi_index_to_python_tuple);
}


template<class I, class X> struct PythonExpansionIterator : public Expansion<I,X>::ConstIterator {
    PythonExpansionIterator(typename Expansion<I,X>::ConstIterator iter)
        : Expansion<I,X>::ConstIterator(iter) { }
    Monomial<I,X> operator*() const {
        auto ac=this->Expansion<I,X>::ConstIterator::operator*();
        return Monomial<I,X>(ac.index(),ac.coefficient());
    }
};
template<class I, class X> decltype(auto) pybegin(Expansion<I,X> const& e) { return PythonExpansionIterator<I,X>(e.begin()); }
template<class I, class X> decltype(auto) pyend(Expansion<I,X> const& e) { return PythonExpansionIterator<I,X>(e.end()); }

#warning Use generic pow
namespace Ariadne {
template<class I, class X> Polynomial<I,X> pypow(Polynomial<I,X> p, Nat m) {
    if (m==0) { return nul(p)+1; }
    Polynomial<I,X> r=p;
    for (Nat i=1; i!=m; ++i) {
        r=r*p;
    }
    return r;
}
} // namespace Ariadne

template<class X>
pybind11::class_< MultivariatePolynomial<X> > export_polynomial(pybind11::module& module)
{
    pybind11::class_< MultivariateMonomial<X> > monomial_class(module,python_template_class_name<X>("MultivariateMonomial").c_str());
    monomial_class.def(pybind11::init<MultiIndex,X>());
    monomial_class.def("index", (MultiIndex const&(MultivariateMonomial<X>::*)()const) &MultivariateMonomial<X>::index);
    monomial_class.def("coefficient", (X const&(MultivariateMonomial<X>::*)()const) &MultivariateMonomial<X>::coefficient);
    monomial_class.def("__str__", &__cstr__<MultivariateMonomial<X>>);

    pybind11::class_< MultivariatePolynomial<X> > polynomial_class(module,python_template_class_name<X>("MultivariatePolynomial").c_str());
    polynomial_class.def(pybind11::init< MultivariatePolynomial<X> >());
    polynomial_class.def(pybind11::init([](pybind11::dict const& dct){return MultivariatePolynomial<X>(expansion_from_python<MultiIndex,X>(dct));}));
    polynomial_class.def( pybind11::init<Expansion<MultiIndex,X>>());
    polynomial_class.def_static("constant", (MultivariatePolynomial<X>(*)(SizeType,X const&)) &MultivariatePolynomial<X>::constant);
    if constexpr (HasPrecisionType<X>) {
        typedef typename X::PrecisionType PR;
        polynomial_class.def(pybind11::init<Nat,PR>());
        polynomial_class.def(pybind11::init([](pybind11::dict const& dct, PR pr){return MultivariatePolynomial<X>(expansion_from_python<MultiIndex,X>(dct,pr));}));
        //polynomial_class.def_static("variable", (MultivariatePolynomial<X>(*)(SizeType,SizeType,PR)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("coordinate", (MultivariatePolynomial<X>(*)(SizeType,SizeType,PR)) &MultivariatePolynomial<X>::coordinate);
        polynomial_class.def_static("coordinates", (Vector<MultivariatePolynomial<X>>(*)(SizeType,PR)) &MultivariatePolynomial<X>::coordinates);
    } else {
        polynomial_class.def(pybind11::init<Nat>());
        //polynomial_class.def_static("variable", (MultivariatePolynomial<X>(*)(SizeType,SizeType)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("coordinate", (MultivariatePolynomial<X>(*)(SizeType,SizeType)) &MultivariatePolynomial<X>::variable);
        polynomial_class.def_static("coordinates", [](Nat as){return MultivariatePolynomial<X>::coordinates(as).array();});
    }

    polynomial_class.def("argument_size", &MultivariatePolynomial<X>::argument_size);
    polynomial_class.def("insert", &MultivariatePolynomial<X>::insert);
    //polynomial_class.def("expansion", (Expansion<MultiIndex,X>const&(MultivariatePolynomial<X>::*)()const) &MultivariatePolynomial<X>::expansion, pybind11::return_value_policy::reference_internal);
    //polynomial_class.def("expansion", [](MultivariatePolynomial<X>const& p){return expansion_to_python_dict(p.expansion());});
    polynomial_class.def("coefficient", (X const&(MultivariatePolynomial<X>::*)(MultiIndex const&)const) &MultivariatePolynomial<X>::operator[]);
    polynomial_class.def("__getitem__", (X const&(MultivariatePolynomial<X>::*)(MultiIndex const&)const) &MultivariatePolynomial<X>::operator[]);
    polynomial_class.def("__iter__", [](MultivariatePolynomial<X> const& p){return pybind11::make_iterator(pybegin(p.expansion()),pyend(p.expansion()));});

    polynomial_class.def("__call__", (X(MultivariatePolynomial<X>::*)(Vector<X>const&)const) &MultivariatePolynomial<X>::operator());
    module.def("evaluate", (X(*)(MultivariatePolynomial<X>const&,Vector<X>const&)) &evaluate);

    define_algebra(module,polynomial_class);
    polynomial_class.def("__pow__", [](MultivariatePolynomial<X> const& p, Nat m){return pypow(p,m);});
    if constexpr (HasGenericType<X>) {
        typedef typename X::GenericType Y;
        define_mixed_arithmetic(module,polynomial_class,Tag<Y>());
    }

    module.def("compose", (MultivariatePolynomial<X>(*)(MultivariatePolynomial<X>const&,Vector<MultivariatePolynomial<X>>const&)) &compose);
    module.def("derivative", (MultivariatePolynomial<X>(*)(MultivariatePolynomial<X>,SizeType j)) &derivative);
    module.def("antiderivative", (MultivariatePolynomial<X>(*)(MultivariatePolynomial<X>,SizeType j)) &antiderivative);

    polynomial_class.def("__str__",&__cstr__<MultivariatePolynomial<X>>);

    export_vector<MultivariatePolynomial<X>>(module, (python_template_class_name<X>("MultivariatePolynomialVector")).c_str());

    return polynomial_class;
}

Void export_polynomials(pybind11::module& module)
{
    export_polynomial<Rational>(module);
    export_polynomial<FloatDPBounds>(module);
    export_polynomial<FloatDPApproximation>(module);
    export_polynomial<FloatMPBounds>(module);
    export_polynomial<FloatMPApproximation>(module);

    template_<MultivariatePolynomial> multivariate_polynomial_template(module);
    multivariate_polynomial_template.instantiate<Rational>();
    multivariate_polynomial_template.instantiate<FloatDPBounds>();
    multivariate_polynomial_template.instantiate<FloatDPApproximation>();
    multivariate_polynomial_template.instantiate<FloatMPBounds>();
    multivariate_polynomial_template.instantiate<FloatMPApproximation>();
}

Void export_domains(pybind11::module& module)
{
    pybind11::class_<RealDomain> real_domain_class(module,"RealDomain");
    real_domain_class.def(pybind11::init<>());
    //    real_domain_class.def(pybind11::init<SizeOne>());
    real_domain_class.def("__repr__", &__cstr__<RealDomain>);

    pybind11::class_<EuclideanDomain> euclidean_domain_class(module,"EuclideanDomain");
    euclidean_domain_class.def(pybind11::init<SizeType>());
    euclidean_domain_class.def("__repr__", &__cstr__<EuclideanDomain>);
}

template<class SIG> String signature_name() {
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    return String(Same<RES,RealScalar> ? "Scalar" : "Vector") + (Same<ARG,RealScalar> ? "Univariate" : "Multivariate");
}

template<class P, class SIG> Void export_function(pybind11::module& module) {
    using F=Function<P,SIG>;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using DomainType = typename SignatureTraits<SIG>::DomainType;
    using CodomainType = typename SignatureTraits<SIG>::CodomainType;
    using ResultSizeType = ElementSizeType<CodomainType>;
    using ArgumentSizeType = ElementSizeType<DomainType>;

    pybind11::class_<F> function_class(module,(class_name<P>()+signature_name<SIG>()+"Function").c_str());

    function_class.def(pybind11::init<F>());
    if constexpr (Same<P,ApproximateTag>) {
        function_class.def(pybind11::init<Function<ValidatedTag,SIG>>());
        function_class.def(pybind11::init<Function<EffectiveTag,SIG>>());
    }
    if constexpr (Same<P,ValidatedTag>) {
        function_class.def(pybind11::init<Function<EffectiveTag,SIG>>());
    }

    function_class.def(pybind11::init<ResultSizeType,ArgumentSizeType>());
    function_class.def(pybind11::init<ResultSizeType,DomainType>());
    if constexpr (Same<RES,RealScalar>) {
        function_class.def(pybind11::init<ArgumentSizeType>());
        function_class.def(pybind11::init<DomainType>());
    }

    function_class.def(pybind11::init<typename SignatureTraits<SIG>::ArgumentSpaceType, typename SignatureTraits<SIG>::template Result<RealExpression>>());

    if constexpr (Same<RES,RealVector>) {
        // NOTE: This must go *after* the conversion constructor
        function_class.def(pybind11::init([](std::vector<Function<P,RealScalar(ARG)>> const& lst){return F(lst);}));
    }

    function_class.def("result_size", &F::result_size);
    function_class.def("argument_size", &F::argument_size);

    if constexpr (Same<RES,RealVector>) {
        function_class.def("__getitem__", &F::get);
        function_class.def("__setitem__", &F::set);
    }

    if constexpr (Same<ARG,RealScalar>) {
        function_class.def("derivative", (F(F::*)(IndexZero)const) &F::derivative);
        function_class.def("derivative", (F(*)(const F&)) &derivative);
        module.def("derivative", (F(F::*)(IndexZero)const) &F::derivative);
        module.def("derivative", (F(*)(const F&)) &derivative);
    } else if constexpr (Same<ARG,RealVector>) {
        function_class.def("derivative", (F(F::*)(SizeType)const) &F::derivative);
        module.def("derivative", (F(F::*)(SizeType)const) &F::derivative);
    }


    if constexpr (Same<RES,RealScalar>) {
        define_elementary_algebra<F,Number<P>>(module,function_class);
        module.def("pow", [](F const& f, Int n){return pow(f,n);});
        function_class.def("__pow__", [](F const& f, Int n){return pow(f,n);});
    }

    // TODO: Put these in C++ API
    // define_vector_algebra_arithmetic<VectorMultivariateFunction<P>,ScalarMultivariateFunction<P>,Number<P>>(module,vector_function_class);
    // FIXME: Define vector function operations for Validated and Approximate
    if constexpr (Same<RES,RealVector>) {
        if constexpr (Same<ARG,RealVector> and Same<P,EffectiveTag>) {
            define_vector_algebra_arithmetic<VectorFunction<P,ARG>,ScalarFunction<P,ARG>>(module,function_class);
        }
    }

    export_function_evaluation(module,function_class);
    export_function_derivatives(module,function_class);

    if constexpr (Same<SIG,RealScalar(RealScalar)>) {
        if constexpr (not Same<P,ApproximateTag>) {
            module.def("derivative", [](const F& f, const FloatDPBounds& x){return static_cast<FloatDPBounds>(differential(f,x,1u).gradient()[0]);} );
            module.def("second_derivative", [](const F& f, const FloatDPBounds& x){return static_cast<FloatDPBounds>(differential(f,x,2u).hessian()[0][0]);} );
        }
        module.def("derivative", [](const F& f, const FloatDPApproximation& x){return static_cast<FloatDPApproximation>(differential(f,x,1u).gradient()[0]);} );
        module.def("second_derivative", [](const F& f, const FloatDPApproximation& x){return static_cast<FloatDPApproximation>(differential(f,x,2u).hessian()[0][0]);} );
    }
    if constexpr (Same<SIG,RealScalar(RealVector)>) {
        function_class.def("gradient", (Covector<FloatDPApproximation>(ScalarMultivariateFunction<P>::*)(const Vector<FloatDPApproximation>&)const) &ScalarMultivariateFunction<P>::gradient);
        if constexpr (not Same<P,ApproximateTag>) {
            function_class.def("gradient", (Covector<FloatDPBounds>(ScalarMultivariateFunction<P>::*)(const Vector<FloatDPBounds>&)const) &ScalarMultivariateFunction<P>::gradient);
        }
    }

    module.def("evaluate", (Scalar<FloatDPApproximation>(*)(const ScalarMultivariateFunction<P>&,const Vector<FloatDPApproximation>&)) &evaluate);
        if constexpr (not Same<P,ApproximateTag>) {
            module.def("evaluate", (Scalar<FloatDPBounds>(*)(const ScalarMultivariateFunction<P>&,const Vector<FloatDPBounds>&)) &evaluate);
        }

    if constexpr (Same<SIG,RealScalar(RealScalar)>) {
        function_class.def_static("constant", (F(*)(Number<P>)) &F::constant);
        function_class.def_static("constant", (F(*)(RealDomain,Number<P>)) &F::constant);
        function_class.def_static("coordinate", (F(*)()) &F::coordinate);
        function_class.def_static("coordinate", (F(*)(RealDomain)) &F::coordinate);
        function_class.def_static("identity", (F(*)()) &F::identity);
        function_class.def_static("identity", (F(*)(RealDomain)) &F::identity);
    }
    if constexpr (Same<SIG,RealVector(RealScalar)>) {
        function_class.def_static("zeros", (F(*)(SizeType)) &F::zeros);
        function_class.def_static("zeros", (F(*)(SizeType,RealDomain)) &F::zeros);
    }
    if constexpr (Same<SIG,RealScalar(RealVector)>) {
        function_class.def_static("constant", (F(*)(SizeType,Number<P>)) &F::constant);
        function_class.def_static("constant", (F(*)(EuclideanDomain,Number<P>)) &F::constant);
        function_class.def_static("coordinate", (F(*)(SizeType,SizeType)) &F::coordinate);
        function_class.def_static("coordinate", (F(*)(EuclideanDomain,SizeType)) &F::coordinate);
    }
    if constexpr (Same<SIG,RealVector(RealVector)>) {
        function_class.def_static("zeros", (F(*)(SizeType,SizeType)) &F::zeros);
        function_class.def_static("zeros", (F(*)(SizeType,EuclideanDomain)) &F::zeros);
        function_class.def_static("identity", (F(*)(SizeType)) &F::identity);
        function_class.def_static("identity", (F(*)(EuclideanDomain)) &F::identity);
    }

    function_class.def("__str__", &__cstr__<F>);
    function_class.def("__repr__", &__repr__<F>);
    function_class.def("__crepr__", &__crepr__<F>);

    if constexpr (Same<RES,RealVector>) if constexpr (Same<ARG,RealVector>) {
        module.def("join", (VectorFunction<P,ARG>(*)(const ScalarFunction<P,ARG>&, const ScalarFunction<P,ARG>&)) &join);
        module.def("join", (VectorFunction<P,ARG>(*)(const VectorFunction<P,ARG>&, const ScalarFunction<P,ARG>&)) &join);
        module.def("join", (VectorFunction<P,ARG>(*)(const ScalarFunction<P,ARG>&, const VectorFunction<P,ARG>&)) &join);
        module.def("join", (VectorFunction<P,ARG>(*)(const VectorFunction<P,ARG>&, const VectorFunction<P,ARG>&)) &join);
    }

    module.def("compose", (Function<P,RES(ARG)>(*)(const Function<P,RES(RealScalar)>&,const Function<P,RealScalar(ARG)>&)) &compose);
    module.def("compose", (Function<P,RES(ARG)>(*)(const Function<P,RES(RealVector)>&,const Function<P,RealVector(ARG)>&)) &compose);

    if constexpr (Same<P,EffectiveTag> and Same<SIG,RealScalar(RealVector)>) {
        function_class.def("__eq__", &__eq__<EffectiveScalarMultivariateFunction,EffectiveNumber>);
        function_class.def("__le__", &__le__<EffectiveScalarMultivariateFunction,EffectiveNumber>);
        function_class.def("__ge__", &__ge__<EffectiveScalarMultivariateFunction,EffectiveNumber>);
    }

}


template<class SIG> Void export_functions(pybind11::module& module)
{
    export_function<EffectiveTag,SIG>(module);
    export_function<ValidatedTag,SIG>(module);
    export_function<ApproximateTag,SIG>(module);

    pybind11::implicitly_convertible<Function<EffectiveTag,SIG>,Function<ValidatedTag,SIG>>();
    pybind11::implicitly_convertible<Function<EffectiveTag,SIG>,Function<ApproximateTag,SIG>>();
    pybind11::implicitly_convertible<Function<ValidatedTag,SIG>,Function<ApproximateTag,SIG>>();

    if constexpr (Same<SIG,RealScalar(RealVector)>) {
        module.def("derivatives", (VectorMultivariateFunction<EffectiveTag>(*)(const ScalarMultivariateFunction<EffectiveTag>&)) &derivatives);
        module.def("lie_derivative", (ScalarMultivariateFunction<EffectiveTag>(*)(const ScalarMultivariateFunction<EffectiveTag>&,const VectorMultivariateFunction<EffectiveTag>&)) &lie_derivative);
    }
}


template<class Y, class X> Void define_procedure_evaluation(pybind11::module& module) {
    module.def("evaluate", (X(*)(Procedure<Y> const&, Vector<X> const&)) &evaluate);
    module.def("gradient", (Covector<X>(*)(Procedure<Y> const&, Vector<X> const&)) &gradient);
    module.def("hessian", (X(*)(Procedure<Y> const&, Vector<X> const&, Vector<X> const&)) &hessian);
}

template<class Y> pybind11::class_<Procedure<Y>> export_procedure(pybind11::module& module) {
    typedef Paradigm<Y> P;
    pybind11::class_<Procedure<Y>> procedure_class(module,(class_name<P>()+"Procedure").c_str());
    procedure_class.def("__str__", &__cstr__<Procedure<Y>>);
    module.def("make_procedure", (Procedure<Y>(*)(ScalarMultivariateFunction<P> const&)) &make_procedure);
    if constexpr (Same<Y,ApproximateNumber>) {
        define_procedure_evaluation<Y,FloatDPApproximation>(module);
        //define_procedure_evaluation<Y,FloatMPApproximation>(module);
    }
    if constexpr (Same<Y,ValidatedNumber>) {
        define_procedure_evaluation<Y,FloatDPBounds>(module);
        //define_procedure_evaluation<Y,FloatMPBounds>(module);
    }
    return procedure_class;
}

Void export_procedures(pybind11::module& module) {
    export_procedure<ApproximateNumber>(module);
    export_procedure<ValidatedNumber>(module);
}



template<class P, class SIG> Void export_function_patch(pybind11::module& module) {
    using FP=FunctionPatch<P,SIG>;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using DomainType = typename FunctionPatch<P,SIG>::DomainType;
    using ArgumentIndexType = ElementIndexType<DomainType>;

    pybind11::class_<FP> function_patch_class(module,(class_name<P>()+signature_name<SIG>()+"FunctionPatch").c_str());
    function_patch_class.def(pybind11::init<FunctionPatch<P,SIG>>());

    function_patch_class.def("result_size", &FP::result_size);
    function_patch_class.def("argument_size", &FP::argument_size);
    function_patch_class.def("domain", &FP::domain);
    function_patch_class.def("codomain", &FP::codomain);
    function_patch_class.def("range", &FP::range);

    function_patch_class.def("__str__", &__cstr__<FP>);

    if constexpr (Same<RES,RealScalar>) {
        define_elementary_algebra<FP,Number<P>>(module,function_patch_class);
    }
    if constexpr (Same<RES,RealVector>) {
        define_vector_algebra_arithmetic<FunctionPatch<P,RealVector(ARG)>,FunctionPatch<P,RealScalar(ARG)>,Number<P>>(module,function_patch_class);
    }

    // TODO: Export all function operators
    // export_function_evaluation(module,function_patch_class);
    // export_function_derivatives(module,function_patch_class);


    if constexpr (IsStronger<P,ApproximateTag>::value) {
        define_call<FloatDPApproximation>(function_patch_class);
        define_call<FloatMPApproximation>(function_patch_class);
    }
    if constexpr (IsStronger<P,ValidatedTag>::value) {
        define_call<FloatDPBounds>(function_patch_class);
        define_call<FloatMPBounds>(function_patch_class);
    }

    if  constexpr (Same<RES,RealVector>) {
        module.def("antiderivative", [](FunctionPatch<P,SIG>const& fp, ArgumentIndexType j){return _antiderivative_(fp,j);});
        module.def("antiderivative", [](FunctionPatch<P,SIG>const& fp, ArgumentIndexType j, Number<P> const& c){return _antiderivative_(fp,j,c);});
    }

    module.def("compose", [](const Function<P,RealScalar(RES)>& f, const FunctionPatch<P,SIG>& gp){return compose(f,gp);});
    module.def("compose", [](const Function<P,RealVector(RES)>& f, const FunctionPatch<P,SIG>& gp){return compose(f,gp);});
    module.def("compose", [](const FunctionPatch<P,RealScalar(RES)>& fp, const FunctionPatch<P,SIG>& gp){return compose(fp,gp);});
    module.def("compose", [](const FunctionPatch<P,RealVector(RES)>& fp, const FunctionPatch<P,SIG>& gp){return compose(fp,gp);});
    module.def("unchecked_compose", [](const FunctionPatch<P,RealScalar(RES)>& fp, const FunctionPatch<P,SIG>& gp){return unchecked_compose(fp,gp);});
    module.def("unchecked_compose", [](const FunctionPatch<P,RealVector(RES)>& fp, const FunctionPatch<P,SIG>& gp){return unchecked_compose(fp,gp);});
}

template<class SIG> Void export_function_patches(pybind11::module& module) {
    export_function_patch<ValidatedTag,SIG>(module);
    export_function_patch<ApproximateTag,SIG>(module);
}

Void export_function_patches(pybind11::module& module) {
    // TODO: Correct export of univariate function patches.
    export_function_patches<RealScalar(RealScalar)>(module);
    // export_function_patches<RealVector(RealScalar)>(module);
    export_function_patches<RealScalar(RealVector)>(module);
    export_function_patches<RealVector(RealVector)>(module);
}


Void function_submodule(pybind11::module& module) {

    export_multi_index(module);

    export_polynomials(module);

    export_domains(module);

    export_functions<RealScalar(RealScalar)>(module);
    export_functions<RealVector(RealScalar)>(module);
    export_functions<RealScalar(RealVector)>(module);
    export_functions<RealVector(RealVector)>(module);

    export_function_patches(module);
    export_procedures(module);

    template_<Function> function_template(module,python_template_name<Function>().c_str());
    function_template.def_new([](RealVariable var, Scalar<RealExpression> e){return EffectiveScalarUnivariateFunction(var,e);});
    function_template.def_new([](RealVariable var, Vector<RealExpression> e){return EffectiveVectorUnivariateFunction(var,e);});
    function_template.def_new([](RealSpace spc, Scalar<RealExpression> e){return EffectiveScalarMultivariateFunction(spc,e);});
    function_template.def_new([](RealSpace spc, Vector<RealExpression> e){return EffectiveVectorMultivariateFunction(spc,e);});

    function_template.def_new([](pybind11::object pyf){return univariate_function_from_python(pyf);});
    function_template.def_new([](SizeType as, pybind11::object pyf){return multivariate_function_from_python(as,pyf);});

}
