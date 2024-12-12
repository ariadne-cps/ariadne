/***************************************************************************
 *            symbolic/expression_patch.hpp
 *
 *  Copyright  2024  Pieter Collins
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

/*! \file symbolic/expression_patch.hpp
 *  \brief Symbolic expressions in named variables over bounded domains
 */

#ifndef ARIADNE_EXPRESSION_PATCH_HPP
#define ARIADNE_EXPRESSION_PATCH_HPP

#include "symbolic/variable.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_set.hpp"

#include "algebra/sweeper.hpp"
#include "function/taylor_function.hpp"

namespace Ariadne {

template<class P, class R> class ExpressionPatch;

template<class X> class ScalarOrVector
    : public Variant<Scalar<X>,Vector<X>>
{
    template<class S> requires Convertible<S,Scalar<X>> ScalarOrVector(S const& s)
        : Variant<Scalar<X>,Vector<X>>(Scalar<X>(s)) { }
    template<class V> requires Convertible<V,Vector<X>> ScalarOrVector(V const& v)
        : Variant<Scalar<X>,Vector<X>>(Scalar<X>(v)) { }
};

Vector<RealExpression> call(EffectiveVectorMultivariateFunction const& f, Vector<RealExpression> const & ev) {
    Vector<ElementaryAlgebra<Real>> av=ev;
    Vector<ElementaryAlgebra<Real>> rav=f(av);
    return Vector<RealExpression>(rav.size(),[&rav](SizeType i){return rav[i].extract<RealExpression>();});
}


template<template<class>class SV> FunctionPatch(RealSpacePatch, ExpressionPatch<ValidatedTag,SV<Real>>) -> FunctionPatch<ValidatedTag,SV<Real>(Vector<Real>)>;


template<class P, class T> class ExpressionPatch;
using ValidatedScalarExpressionPatch = ExpressionPatch<ValidatedTag,RealScalar>;
using ValidatedVectorExpressionPatch = ExpressionPatch<ValidatedTag,RealVector>;

typedef VariableInterval<typename IntervalDomainType::UpperBoundType> VariableIntervalDomainType;
typedef VariablesBox<IntervalDomainType> VariablesBoxDomainType;
typedef VectorVariableBox<IntervalDomainType> VectorVariableBoxDomainType;



template<class X, class F> decltype(auto) vapply(F const& f, Vector<X> const& v) {
    using Y = decltype(f(declval<X>()));
    return Vector<Y>(v.size(),[&f,&v](SizeType i){return f(v[i]);});
}

template<class T> class ConstantOrVariable;

template<> class ConstantOrVariable<Real> : public Variant<Constant<Real>,Variable<Real>> {
    using Variant<Constant<Real>,Variable<Real>>::Variant;
};

using RealVectorConstant = Constant<Vector<Real>>;


template<class T> concept IsConstantOrVariable = OneOf<T,Dyadic,RealConstant,RealVariable,RealVectorConstant,RealVectorVariable>;

template<class... FS> struct AllConstantOrVariable;
template<> struct AllConstantOrVariable<> : True { };
template<class F0, class... FS> struct AllConstantOrVariable<F0,FS...> {
    static const bool value = IsConstantOrVariable<F0> && AllConstantOrVariable<FS...>::value;
};

template<class... FS> concept AllAreConstantOrVariable = AllConstantOrVariable<FS...>::value;



template<> class ConstantOrVariable<Vector<Real>> : public Vector<ConstantOrVariable<Real>> {
  public:
    using Vector<ConstantOrVariable<Real>>::Vector;
    template<class... AS> requires AllAreConstantOrVariable<AS...>
        ConstantOrVariable(AS... as) : ConstantOrVariable(_make_list(as...)) { }
  private:
    ConstantOrVariable(Constant<Vector<Real>>);
    ConstantOrVariable(Variable<Vector<Real>>);

    template<class... AS> static List<ConstantOrVariable<Real>> _make_list(AS... as) {
        List<ConstantOrVariable<Real>> lst; _append_list(lst,as...); return lst; }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Dyadic a) {
        lst.append(Constant<Real>(a)); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Constant<Real> a) {
        lst.append(a); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Variable<Real> a) {
        lst.append(a); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, ConstantOrVariable<Real> a) {
        lst.append(a); }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Variable<Vector<Real>> a) {
        for (SizeType i=0; i!=a.size(); ++i) { lst.append(a[i]); } }
    static Void _append(List<ConstantOrVariable<Real>>& lst, Constant<Vector<Real>> a) {
        for (SizeType i=0; i!=a.size(); ++i) { lst.append(Constant<Real>(a[i])); } }
    static Void _append(List<ConstantOrVariable<Real>>& lst, ConstantOrVariable<Vector<Real>> a) {
        for (SizeType i=0; i!=a.size(); ++i) { lst.append(a[i]); } }
    static Void _append_list(List<ConstantOrVariable<Real>>& lst) { }
    template<class A0, class... AS> static Void _append_list(List<ConstantOrVariable<Real>>& lst, A0 a0, AS... as) {
        _append(lst,a0); _append_list(lst,as...); }

    template<class... AS> friend List<ConstantOrVariable<Real>> make_constant_or_variable_list(ConstantOrVariable<AS> const& ... cvs);
};

template<class... AS> List<ConstantOrVariable<Real>> make_constant_or_variable_list(ConstantOrVariable<AS> const& ... cvs) {
    return ConstantOrVariable<Vector<Real>>::_make_list(cvs...);
}




template<class T> class ConstantOrVariablePatch;

template<> class ConstantOrVariablePatch<Real> : public Variant<RealConstant,VariableIntervalDomainType> { };

template<> class ConstantOrVariablePatch<Vector<Real>>
    : public Variant<RealVectorConstant,VectorVariableBoxDomainType> { };



class AnyConstantOrVariablePatch
    : public Variant<RealConstant,RealVectorConstant,VariableIntervalDomainType,VectorVariableBoxDomainType>
{
    typedef Variant<RealConstant,RealVectorConstant,VariableIntervalDomainType,VectorVariableBoxDomainType> Base;
  public:
    AnyConstantOrVariablePatch(Dyadic h) : Base(RealConstant(h)) { }
    AnyConstantOrVariablePatch(RealConstant c) : Base(c) { }
    AnyConstantOrVariablePatch(VariableIntervalDomainType vivl) : Base(vivl) { }
    AnyConstantOrVariablePatch(VariablesBoxDomainType vsbx);
    AnyConstantOrVariablePatch(VectorVariableBoxDomainType vvbx) : Base(vvbx) { }
};

template<class T>
Void _append(List<ConstantOrVariablePatch<Real>>& lst, T const& t) {
    if constexpr (OneOf<T,Constant<RealVector>,Variable<RealVector>,Vector<RealVariable>>) {
        for (SizeType i=0; i!=t.size(); ++i) { lst.append(t[i]); }
    } else {
        lst.append(t);
    }
}

template<class... AS> List<ConstantOrVariablePatch<Real>> make_constant_or_variable_patch_list(InitializerList<AnyConstantOrVariablePatch> const& cvplst) {
    List<ConstantOrVariablePatch<Real>> lst;
    for(auto cvp : cvplst) { std::visit(cvp,[&lst](AnyConstantOrVariablePatch const& cvp){_append(lst,cvp);}); }
    return lst;
}




template<class P> VectorMultivariateFunction<P> coordinates(SizeType n, Array<SizeType> is) {
    if (is.size()==0) { return VectorMultivariateFunction<P>(0,n); }
    return Vector<ScalarMultivariateFunction<P>>(is.size(), [n,&is](SizeType j){return ScalarMultivariateFunction<P>::coordinate(n,is[j]);});
}
template<class P> VectorMultivariateFunction<P> coordinates(SizeType n, Range rng) {
    if (rng.size()==0) { return VectorMultivariateFunction<P>(0,n); }
    return Vector<ScalarMultivariateFunction<P>>(rng.size(), [n,rng](SizeType j){return ScalarMultivariateFunction<P>::coordinate(n,rng[j]);});
}


template<> class SpacePatch<Real> {
    using T=Real;
    RealSpace _spc;
    BoxDomainType _dom;
  public:
    SpacePatch(InitializerList<SpacePatch<T>> lst) : SpacePatch(*lst.begin()) {
        assert (lst.size()!=0); bool first=true;
        for (auto sp : lst) {
            if (first) { first=false; } else { _spc=join(_spc,sp._spc); _dom=product(_dom,sp._dom); } } }
    SpacePatch(VariableIntervalDomainType vivl) : _spc({vivl.variable()}), _dom({vivl.interval()}) { }
    SpacePatch(VectorVariableBoxDomainType vvbx) : _spc(vvbx.variable()), _dom(vvbx.box()) { }
    SpacePatch(List<RealVariable> vars, BoxDomainType dom)
        : _spc(vars), _dom(dom) { }
    SpacePatch(List<RealVariable> vars, Map<RealVariable,IntervalDomainType> vivls)
        : _spc(vars), _dom(vars.size(),[&](SizeType i){return vivls[vars[i]];}) { }

    DimensionType dimension() const { return this->_spc.dimension(); }
    Space<T> const& space() const { return this->_spc; }
    List<RealVariable> variables() const { return this->_spc.variables(); }
    explicit operator Space<T> const& () const { return this->_spc; }
    BoxDomainType domain() const { return this->_dom; }
    friend OutputStream& operator<<(OutputStream& os, SpacePatch<T> const& spcdom) {
        return os <<"SpacePatch(spc="<<spcdom._spc<<", dom="<<spcdom._dom<<")"; }
};




template<class P> VectorMultivariateFunction<P> combine(VectorMultivariateFunction<P> f1, VectorMultivariateFunction<P> f2) {
//    return make_combined_function(f1,f2);
    auto n1=f1.argument_size(); auto n2=f2.argument_size(); auto n=n1+n2;
    auto p2=coordinates<P>(n,range(n1,n1+n2));
    auto p1=coordinates<P>(n,range(0,n1));
    return join(compose(f1,p1),compose(f2,p2));
}


template<> class ExpressionPatch<ValidatedTag, RealScalar> {
    using P=ValidatedTag;
    using T=RealScalar;
    using SIG=T(RealVector);
  private: public:
    Vector<RealVariable> _vars;
    FunctionPatch<P,RealScalar(RealVector)> _fp;
  private:
  public:
    IntervalRangeType range();
    SizeOne size() const { return _fp.result_size(); }

    ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, RealExpression expr);
    ExpressionPatch(RealSpacePatch spcptch, RealExpression expr);
    ExpressionPatch(Vector<RealVariable> vars, BoxDomainType dom, Function<P,SIG> f);
    ExpressionPatch(RealSpacePatch spcptch, Function<P,SIG> f);
    ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp);
    ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp);

    friend auto operator+(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> { return pos(ep); }
    friend auto operator-(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> { return neg(ep); }

    friend auto operator+(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return add(ep1,ep2); }
    friend auto operator-(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return sub(ep1,ep2); }
    friend auto operator*(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return mul(ep1,ep2); }
    friend auto operator/(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return div(ep1,ep2); }
    friend auto operator+(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> { return add(ep1,c2); }
    friend auto operator-(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> { return sub(ep1,c2); }
    friend auto operator*(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> { return mul(ep1,c2); }
    friend auto operator/(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> { return div(ep1,c2); }
    friend auto operator+(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return add(c1,ep2); }
    friend auto operator-(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return sub(c1,ep2); }
    friend auto operator*(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return mul(c1,ep2); }
    friend auto operator/(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> { return div(c1,ep2); }

    friend auto operator*(ExpressionPatch<P,Scalar<Real>> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,Vector<Real>>;
    friend auto operator*(Vector<Number<P>> c1, ExpressionPatch<P,Scalar<Real>> ep2) -> ExpressionPatch<P,Vector<Real>>;
    friend auto operator/(Vector<Number<P>> c1, ExpressionPatch<P,Scalar<Real>> ep2) -> ExpressionPatch<P,Vector<Real>>;

    friend auto add(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Add(),ep1,ep2); }
    friend auto sub(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Sub(),ep1,ep2); }
    friend auto mul(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Mul(),ep1,ep2); }
    friend auto div(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Div(),ep1,ep2); }
    friend auto add(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply([&c2](auto x1){return add(x1,c2);},ep1); }
    friend auto sub(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply([&c2](auto x1){return sub(x1,c2);},ep1); }
    friend auto mul(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply([&c2](auto x1){return mul(x1,c2);},ep1); }
    friend auto div(ExpressionPatch<P,T> ep1, Number<P> c2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply([&c2](auto x1){return div(x1,c2);},ep1); }
    friend auto add(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(std::bind(Add(),c1,std::placeholders::_1),ep2); }
    friend auto sub(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(std::bind(Sub(),c1,std::placeholders::_1),ep2); }
    friend auto mul(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(std::bind(Mul(),c1,std::placeholders::_1),ep2); }
    friend auto div(Number<P> c1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(std::bind(Div(),c1,std::placeholders::_1),ep2); }

    friend auto nul(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Nul(),ep); }
    friend auto pos(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Pos(),ep); }
    friend auto neg(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Neg(),ep); }
    friend auto sqr(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Sqr(),ep); }
    friend auto hlf(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Hlf(),ep); }
    friend auto rec(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Rec(),ep); }
    friend auto pow(ExpressionPatch<P,T> ep, Int n) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(std::bind(Pow(),std::placeholders::_1,n),ep); }
    friend auto sqrt(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Sqrt(),ep); }
    friend auto exp(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Exp(),ep); }
    friend auto log(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Log(),ep); }
    friend auto sin(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Sin(),ep); }
    friend auto cos(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Cos(),ep); }
    friend auto tan(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Tan(),ep); }
    friend auto atan(ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Atan(),ep); }

    friend auto mul(ExpressionPatch<P,RealScalar> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,RealVector>;
    friend auto mul(Vector<Number<P>> c1, ExpressionPatch<P,RealScalar> ep2) -> ExpressionPatch<P,RealVector>;
    friend auto div(Vector<Number<P>> c1, ExpressionPatch<P,RealScalar> ep2) -> ExpressionPatch<P,RealVector>;

    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<P,T> const& ep) {
        return ep._write(os); }
  private:
    OutputStream& _write(OutputStream& os) const;
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,T> ep);
};

template<> class ExpressionPatch<ValidatedTag, RealVector> {
    using P=ValidatedTag;
    using T=RealVector;
    using SIG=T(RealVector);
    using VT=RealVector;
    using ST=RealScalar;
  private: public:
    Vector<RealVariable> _vars;
    FunctionPatch<P,RealVector(RealVector)> _fp;
  private:
  public:
    Map<RealVariable,IntervalDomainType> domains() const;

    BoxRangeType range();
    SizeType size() const { return _fp.result_size(); }

    template<class VD> requires Same<VD,Map<RealVariable,IntervalDomainType>>
        ExpressionPatch(VD dom, Vector<RealExpression> expr)
            : ExpressionPatch(RealSpacePatch(List<RealVariable>(dom.keys()),dom),expr) { }
#warning Constructing from a domain map and a vector of expressions should be default
    //ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, Vector<RealExpression> expr);
        ExpressionPatch(Vector<Number<P>> c);
    ExpressionPatch(RealSpacePatch dom, Vector<RealExpression> expr);
    ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f);
    ExpressionPatch(RealSpacePatch dom, Function<P,SIG> f);
    ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> f);
    ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp);

    typedef Variant<ExpressionPatch<P,RealVector>,ExpressionPatch<P,RealScalar>,VariableIntervalDomainType,RealConstant,Real> PreExpressionPatch;

    ExpressionPatch(InitializerList<PreExpressionPatch> lst);

    friend auto call(FunctionPatch<P,RealVector(RealVector)>, ExpressionPatch<P,RealVector>) -> ExpressionPatch<P,RealVector>;

    friend auto nul(ExpressionPatch<P,VT> ep) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Nul(),ep); }
    friend auto pos(ExpressionPatch<P,VT> ep) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Pos(),ep); }
    friend auto neg(ExpressionPatch<P,VT> ep) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Neg(),ep); }

    friend auto add(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Add(),ep1,ep2); }
    friend auto sub(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Sub(),ep1,ep2); }
    friend auto mul(ExpressionPatch<P,ST> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),ep1,ep2); }
    friend auto mul(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),ep1,ep2); }
    friend auto div(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_apply(Div(),ep1,ep2); }

    friend auto add(ExpressionPatch<P,VT> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Add(),ep1,c2); }
    friend auto sub(ExpressionPatch<P,VT> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Sub(),ep1,c2); }
    friend auto mul(ExpressionPatch<P,ST> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),ep1,c2); }
    friend auto mul(ExpressionPatch<P,VT> ep1, Scalar<Number<P>> c2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),ep1,c2); }
    friend auto div(ExpressionPatch<P,VT> ep1, Scalar<Number<P>> c2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Div(),ep1,c2); }
    friend auto add(Vector<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Add(),c1,ep2); }
    friend auto sub(Vector<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Sub(),c1,ep2); }
    friend auto mul(Scalar<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),c1,ep2); }
    friend auto mul(Vector<Number<P>> c1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),c1,ep2); }
    friend auto div(Vector<Number<P>> c1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Div(),c1,ep2); }

    friend auto add(ExpressionPatch<P,VT> ep1, Expression<VT> e2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Add(),ep1,e2); }
    friend auto sub(ExpressionPatch<P,VT> ep1, Expression<VT> e2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Sub(),ep1,e2); }
    friend auto mul(ExpressionPatch<P,ST> ep1, Expression<VT> e2) -> ExpressionPatch<P,VT>;
    friend auto mul(ExpressionPatch<P,VT> ep1, Expression<ST> e2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),ep1,e2); }
    friend auto div(ExpressionPatch<P,VT> ep1, Expression<ST> e2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Div(),ep1,e2); }
    friend auto add(Expression<VT> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Add(),e1,ep2); }
    friend auto sub(Expression<VT> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Sub(),e1,ep2); }
    friend auto mul(Expression<ST> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
        return ExpressionPatch<P,T>::_apply(Mul(),e1,ep2); }
    friend auto mul(Expression<VT> e1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT>;
    friend auto div(Expression<VT> e1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT>;

    friend auto operator+(ExpressionPatch<P,VT> ep) -> ExpressionPatch<P,T> { return pos(ep); }
    friend auto operator-(ExpressionPatch<P,VT> ep) -> ExpressionPatch<P,T> { return neg(ep); }
    friend auto operator+(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return add(ep1,ep2); }
    friend auto operator-(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return sub(ep1,ep2); }
    friend auto operator*(ExpressionPatch<P,ST> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return mul(ep1,ep2); }
    friend auto operator*(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T> { return mul(ep1,ep2); }
    friend auto operator/(ExpressionPatch<P,VT> ep1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T> { return div(ep1,ep2); }

    friend auto operator+(ExpressionPatch<P,VT> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,T> { return add(ep1,c2); }
    friend auto operator-(ExpressionPatch<P,VT> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,T> { return sub(ep1,c2); }
    friend auto operator*(ExpressionPatch<P,ST> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,T> { return _apply(Mul(),ep1,c2); }
    friend auto operator*(ExpressionPatch<P,VT> ep1, Scalar<Number<P>> c2) -> ExpressionPatch<P,T> { return _apply(Mul(),ep1,c2); }
    friend auto operator/(ExpressionPatch<P,VT> ep1, Scalar<Number<P>> c2) -> ExpressionPatch<P,T> { return _apply(Div(),ep1,c2); }
    friend auto operator+(Vector<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return add(c1,ep2); }
    friend auto operator-(Vector<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return sub(c1,ep2); }
    friend auto operator*(Scalar<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return _apply(Mul(),c1,ep2); }
    friend auto operator*(Vector<Number<P>> c1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T> { return _apply(Mul(),c1,ep2); }
    friend auto operator/(Vector<Number<P>> c1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T> { return _apply(Div(),c1,ep2); }

    friend auto operator+(ExpressionPatch<P,VT> ep1, Expression<VT> e2) -> ExpressionPatch<P,T> { return add(ep1,e2); }
    friend auto operator-(ExpressionPatch<P,VT> ep1, Expression<VT> e2) -> ExpressionPatch<P,T> { return sub(ep1,e2); }
    friend auto operator*(ExpressionPatch<P,ST> ep1, Expression<VT> e2) -> ExpressionPatch<P,T>;
    friend auto operator*(ExpressionPatch<P,VT> ep1, Expression<ST> e2) -> ExpressionPatch<P,T> { return mul(ep1,e2); }
    friend auto operator/(ExpressionPatch<P,VT> ep1, Expression<ST> e2) -> ExpressionPatch<P,T> { return div(ep1,e2); }
    friend auto operator+(Expression<VT> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return add(e1,ep2); }
    friend auto operator-(Expression<VT> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return sub(e1,ep2); }
    friend auto operator*(Expression<ST> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,T> { return mul(e1,ep2); }
    friend auto operator*(Expression<VT> e1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T>;
    friend auto operator/(Expression<VT> e1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,T>;

    friend auto join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
        return ExpressionPatch<P,T>::_join(ep1,ep2); }

    friend OutputStream& operator<<(OutputStream& os, ExpressionPatch<P,T> const& ep) {
        return ep._write(os); }

    friend ExpressionPatch<ValidatedTag,RealVector> FunctionPatch<ValidatedTag,RealVector(RealVector)>::operator() (ExpressionPatch<ValidatedTag,RealVector> const& ep) const;

  private:
    static ExpressionPatch<P,T> _join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2);
    template<class OP> static ExpressionPatch<P,VT> _apply(OP op, ExpressionPatch<P,VT> ep1, ExpressionPatch<P,VT> ep2);
    template<class OP> static ExpressionPatch<P,VT> _apply(OP op, ExpressionPatch<P,ST> ep1, ExpressionPatch<P,VT> ep2);
    template<class OP> static ExpressionPatch<P,VT> _apply(OP op, ExpressionPatch<P,VT> ep1, ExpressionPatch<P,ST> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,VT> ep1, Vector<Number<P>> c2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,ST> ep1, Vector<Number<P>> c2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,VT> ep1, Scalar<Number<P>> c2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, Vector<Number<P>> c1, ExpressionPatch<P,VT> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, Scalar<Number<P>> c1, ExpressionPatch<P,VT> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, Vector<Number<P>> c1, ExpressionPatch<P,ST> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,VT> ep1, Expression<VT> e2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,VT> ep1, Expression<ST> e2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, Expression<VT> e1, ExpressionPatch<P,VT> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, Expression<ST> e1, ExpressionPatch<P,VT> ep2);
    template<class OP> static ExpressionPatch<P,T> _apply(OP op, ExpressionPatch<P,T> ep);

    OutputStream& _write(OutputStream& os) const;
};

template<class P, class T1, class T2>
Tuple<RealSpace,MultivariateFunctionPatch<P,T1>,MultivariateFunctionPatch<P,T2>>
make_common_variables(ExpressionPatch<P,T1> ep1, ExpressionPatch<P,T2> ep2);


ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp)
    : _vars(vars), _fp(fp)
{
    assert(fp.argument_size()==vars.size());
}

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(RealSpacePatch spcdom, RealExpression expr)
    : ExpressionPatch(spcdom,Function<P,SIG>(spcdom.variables(),expr)) {
}

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(RealSpacePatch spcdom, Function<P,SIG> f)
    : _vars(spcdom.variables()), _fp(spcdom.domain(),f)
{
    ARIADNE_PRECONDITION(spcdom.dimension()==f.argument_size());
}

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp)
    : _vars(spc.variables()), _fp(fp)
{
    ARIADNE_PRECONDITION(spc.dimension()==fp.argument_size());
}

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, RealExpression expr)
    : ExpressionPatch(RealSpacePatch(List<RealVariable>(dom.keys()),dom),expr) {
}

ExpressionPatch<ValidatedTag,RealScalar>::ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f)
    : _vars(vars), _fp(fdom,f)
{
    assert(fdom.dimension()==f.argument_size());
    assert(f.argument_size()==vars.size());
}


template<class OP> auto ExpressionPatch<ValidatedTag,RealScalar>::_apply(OP op, ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
    auto sff = make_common_variables(ep1,ep2);
    return ExpressionPatch<P,T>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealScalar>::_apply(OP op, ExpressionPatch<P,T> ep) -> ExpressionPatch<P,T> {
    return ExpressionPatch<P,T>(ep._vars,op(ep._fp));
}





//ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(Map<RealVariable,IntervalDomainType> dom, Vector<RealExpression> expr)
//    : ExpressionPatch(RealSpacePatch(List<RealVariable>(dom.keys()),dom),expr) {
//}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(RealSpacePatch spcdom, Vector<RealExpression> expr)
    : ExpressionPatch(spcdom,Function<P,SIG>(spcdom.variables(),expr)) {
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(RealSpacePatch spcdom, Function<P,SIG> f)
    : _vars(spcdom.variables()), _fp(spcdom.domain(),f)
{
    ARIADNE_PRECONDITION(spcdom.dimension()==f.argument_size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(RealSpace spc, FunctionPatch<P,SIG> fp)
    : _vars(spc.variables()), _fp(fp)
{
    ARIADNE_PRECONDITION(spc.dimension()==fp.argument_size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(Vector<RealVariable> vars, BoxDomainType fdom, Function<P,SIG> f)
    : _vars(vars), _fp(fdom,f)
{
    assert(fdom.dimension()==f.argument_size());
    assert(f.argument_size()==vars.size());
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp)
    : _vars(vars), _fp(fp)
{
    assert(fp.argument_size()==vars.size());
}

auto ExpressionPatch<ValidatedTag,RealVector>::domains() const -> Map<RealVariable,IntervalDomainType> {
    Map<RealVariable,IntervalDomainType> doms;
    for (SizeType i=0; i!=_vars.size(); ++i) {
        RealVariable v=_vars[i];
        IntervalDomainType dom=_fp.domain()[i];
        if (doms.has_key(v)) { doms[v]=intersection(doms[v],dom); }
        else { doms[v]=dom; }
    }
    return doms;
}

template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
    auto sff = make_common_variables(ep1,ep2);
    return ExpressionPatch<P,VT>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT> {
    auto sff = make_common_variables(ep1,ep2);
    return ExpressionPatch<P,VT>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,ST> ep1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
    auto sff = make_common_variables(ep1,ep2);
    return ExpressionPatch<P,VT>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep1._vars, op(ep1._fp,c2));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,ST> ep1, Vector<Number<P>> c2) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep1._vars, op(ep1._fp,c2));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep1, Scalar<Number<P>> c2) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep1._vars, op(ep1._fp,c2));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, Vector<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep2._vars, op(c1,ep2._fp));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, Scalar<Number<P>> c1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep2._vars, op(c1,ep2._fp));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, Vector<Number<P>> c1, ExpressionPatch<P,ST> ep2) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep2._vars, op(c1,ep2._fp));
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep1, Expression<VT> e2) -> ExpressionPatch<P,VT> {
    ExpressionPatch<ValidatedTag,RealVector> ep2(ep1.domains(),Vector<Expression<Real>>(e2));
    return _apply(op,ep1,ep2);
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep1, Expression<ST> e2) -> ExpressionPatch<P,VT> {
    ExpressionPatch<ValidatedTag,RealScalar> ep2(ep1.domains(),Scalar<Expression<Real>>(e2));
    return _apply(op,ep1,ep2);
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, Expression<VT> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
    ExpressionPatch<ValidatedTag,RealVector> ep1(ep2.domains(),Vector<Expression<Real>>(e1));
    return _apply(op,ep1,ep2);
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, Expression<ST> e1, ExpressionPatch<P,VT> ep2) -> ExpressionPatch<P,VT> {
    ExpressionPatch<ValidatedTag,RealScalar> ep1(ep2.domains(),Scalar<Expression<Real>>(e1));
    return _apply(op,ep1,ep2);
}
template<class OP> auto ExpressionPatch<ValidatedTag,RealVector>::_apply(OP op, ExpressionPatch<P,VT> ep) -> ExpressionPatch<P,VT> {
    return ExpressionPatch<P,VT>(ep._vars, op(ep._fp));
}



#warning Overloads provided ad-hoc to allow certain constructions for systems
ValidatedScalarExpressionPatch operator-(VariableIntervalDomainType xr, Real c) {
    RealSpacePatch spc(xr);
    ValidatedScalarExpressionPatch ep1(spc,xr.variable());
    ValidatedScalarExpressionPatch ep2(spc,c);
    return ep1-ep2;
}

ValidatedVectorExpressionPatch operator-(ValidatedVectorExpressionPatch ep1, VectorVariableBoxDomainType vp2) {
    RealSpacePatch spc(vp2);
    ValidatedVectorExpressionPatch ep2(spc,Vector<Expression<Real>>(Expression<RealVector>(vp2.variable())));
    return ep1-ep2;
}



template<class CVE> Void
append_scalars(List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarExpressionPatch>>& lst, CVE const& cve) {
    if constexpr (Same<CVE,VariablesBoxDomainType>) {
        for (SizeType i=0; i!=cve.size(); ++i) { lst.append(cve[i]); }
    } else if constexpr (Same<CVE,ValidatedVectorExpressionPatch>) {
        for (SizeType i=0; i!=cve.size(); ++i) { lst.append(ValidatedScalarExpressionPatch(cve._vars,cve._fp[i])); }
    } else if constexpr (Same<CVE,Real>) {
        lst.append(RealConstant(cve));
    } else {
        lst.append(cve);
    }
}

List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarExpressionPatch>> make_scalar_list(InitializerList<ValidatedVectorExpressionPatch::PreExpressionPatch> lst) {
    List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarExpressionPatch>> res;
    for (auto cve : lst) {
        std::visit([&res](auto cve){append_scalars(res,cve);},cve);
    }
    return res;
}

ExpressionPatch<ValidatedTag,RealVector>
make_expression_patch(List<Variant<RealConstant,VariableIntervalDomainType,ValidatedScalarExpressionPatch>> const& lst) {
    List<RealVariable> vars;
    Map<RealVariable,SizeType> inds;
    Map<RealVariable,IntervalDomainType> doms;
    for (auto cve : lst) {
        if (std::holds_alternative<VariableIntervalDomainType>(cve)) {
            VariableIntervalDomainType const& vi = std::get<VariableIntervalDomainType>(cve);
            RealVariable const& v=vi.variable();
            IntervalDomainType dom=vi.interval();
            if (doms.has_key(v)) {
                doms[v]=intersection(doms[v],dom);
            } else {
                inds[v]=vars.size();
                vars.append(v);
                doms[v]=dom;
            }
        } else if (std::holds_alternative<ValidatedScalarExpressionPatch>(cve)) {
            ValidatedScalarExpressionPatch const& ep = std::get<ValidatedScalarExpressionPatch>(cve);
            Vector<RealVariable> const& argvars=ep._vars;
            for (SizeType i=0; i!=argvars.size(); ++i) {
                RealVariable var = argvars[i];
                IntervalDomainType dom = ep._fp.domain()[i];
                if (doms.has_key(var)) {
                    doms[var]=intersection(doms[var],dom);
                } else {
                    inds[var]=vars.size();
                    vars.append(var);
                    doms[var]=dom;
                }
            }
        }
    }
    BoxDomainType bxdom(vars.size(), [&](SizeType i){return doms[vars[i]]; });

    ValidatedScalarMultivariateTaylorFunctionModelDP tfm(bxdom,ThresholdSweeper<FloatDP>(dp,1e-8));
    auto fctry = factory(tfm);

    List<ValidatedScalarMultivariateFunctionPatch> fps;
    for (auto cve : lst) {
        if (std::holds_alternative<RealConstant>(cve)) {
            RealConstant const& rc = std::get<RealConstant>(cve);
            fps.append(fctry.create_constant(rc));
        } else if (std::holds_alternative<VariableIntervalDomainType>(cve)) {
            VariableIntervalDomainType const& vi = std::get<VariableIntervalDomainType>(cve);
            RealVariable const& v=vi.variable();
            SizeType ind=inds[v];
            fps.append(fctry.create_coordinate(ind));
        } else if (std::holds_alternative<ValidatedScalarExpressionPatch>(cve)) {
            ValidatedScalarExpressionPatch const& ep = std::get<ValidatedScalarExpressionPatch>(cve);
            ValidatedVectorMultivariateFunctionPatch projection=fctry.create_zeros(ep._vars.size());
            Vector<RealVariable> const& argvars=ep._vars;
            for (SizeType k=0; k!=argvars.size(); ++k) {
                projection[k]=static_cast<ValidatedScalarMultivariateFunctionPatch>(fctry.create_coordinate(inds[argvars[k]]));
            }
            fps.append(compose(ep._fp,projection));
        }
    }
    return ValidatedVectorExpressionPatch(vars,fps);
}

ExpressionPatch<ValidatedTag,RealVector>::ExpressionPatch(InitializerList<PreExpressionPatch> lst)
    : ExpressionPatch(make_expression_patch(make_scalar_list(lst)))
{
}




auto call(FunctionPatch<ValidatedTag,RealVector(RealVector)>, ExpressionPatch<ValidatedTag,RealVector>)
    -> ExpressionPatch<ValidatedTag,RealVector>;


template<class P, class... ARGS> FunctionPatch<P,RealScalar(ARGS...)>::
FunctionPatch(RealSpacePatch const& spc, ExpressionPatch<P,RES> const& ep) {
    BoxDomainType dom=spc.domain();

    Array<SizeType> coords(ep._vars.size());
    for(SizeType i=0; i!=coords.size(); ++i) {
        coords[i]=spc.space()[ep._vars[i]];
    }
    VectorMultivariateFunction<P> pr=coordinates<P>(spc.dimension(),coords);

    ValidatedVectorMultivariateFunctionPatch proj=factory(ep._fp).create(dom,pr);

    *this = compose(ep._fp,proj);
}

template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::
FunctionPatch(RealSpacePatch const& spc, ExpressionPatch<P,RES> const& ep) {
    BoxDomainType dom=spc.domain();

    Array<SizeType> coords(ep._vars.size());
    for(SizeType i=0; i!=coords.size(); ++i) {
        coords[i]=spc.space()[ep._vars[i]];
    }
    VectorMultivariateFunction<P> pr=coordinates<P>(spc.dimension(),coords);

    ValidatedVectorMultivariateFunctionPatch proj=factory(ep._fp).create(dom,pr);

    *this = compose(ep._fp,proj);
}

template<class P, class RES, class... ARGS> auto
make_expression_patch(FunctionPatch<P,RES(ARGS...)> const& fp, List<ConstantOrVariable<Real>> const& cvl) {

    ARIADNE_ASSERT(cvl.size()==fp.argument_size());

    List<RealVariable> variable_list;
    Map<RealVariable,SizeType> indices;
    Map<RealVariable,IntervalDomainType> domains;
    SizeType k=0;
    for (SizeType i=0; i!=cvl.size(); ++i) {
        if (std::holds_alternative<Variable<Real>>(cvl[i])) {
            RealVariable v=std::get<RealVariable>(cvl[i]);
            if (not indices.has_key(v)) {
                variable_list.append(v);
                indices[v]=k;
                domains[v]=fp.domain()[i];
                ++k;
            } else {
                domains[v]=intersection(domains[v],fp.domain()[i]);
            }
        }
    }

    Vector<RealVariable> variables(variable_list);
    BoxDomainType domain(variables.size(),[&](SizeType i){return domains[variables[i]];});

    FunctionPatch<P,RealVector(ARGS...)> projection(domain, Function<P,RealVector(ARGS...)>::zeros(fp.argument_size(),domain.dimension()));
    for (SizeType i=0; i!=fp.argument_size(); ++i) {
        if (std::holds_alternative<Variable<Real>>(cvl[i])) {
            RealVariable v=std::get<RealVariable>(cvl[i]);
            projection[i]=Function<P,RealScalar(ARGS...)>::coordinate(domain.dimension(), indices[v]);
        } else {
            RealConstant c=std::get<RealConstant>(cvl[i]);
            projection[i]=Function<P,RealScalar(ARGS...)>::constant(domain.dimension(),c);
        }
    }

    return ExpressionPatch<P,RES>(variables,domain,compose(fp,projection));
}

template<class P,class T1, class T2>
Tuple<RealSpace,MultivariateFunctionPatch<P,T1>,MultivariateFunctionPatch<P,T2>>
make_common_variables(ExpressionPatch<P,T1> ep1, ExpressionPatch<P,T2> ep2) {

    // For each variable used in ep1 or ep2, give it an index,
    // and find its domain, which is the intersection of the two original domains for a common variable.
    List<RealVariable> variables;
    Map<RealVariable,SizeType> indices;
    Map<RealVariable,IntervalDomainType> domains;
    SizeType as1=ep1._vars.size();
    SizeType as2=ep2._vars.size();
    for (SizeType i=0; i!=as1; ++i) {
        auto v=ep1._vars[i];
        auto d=ep1._fp.domain()[i];
        variables.append(v);
        indices[v]=i;
        domains.insert(v,d);
    }
    SizeType asr=as1;
    for (SizeType i=0; i!=as2; ++i) {
        auto v=ep2._vars[i];
        auto d=ep2._fp.domain()[i];
        if (domains.has_key(v)) {
            domains[v]=intersection(domains[v],d);
        } else {
            variables.append(v);
            indices[v]=asr; ++asr;
            domains.insert(v,d);
        }
    }

    // Construct the ordered domain
    BoxDomainType domain(domains.size(),[&](SizeType i){return domains[variables[i]];});

    // Make projection functions from the common variables onto those used by each individual expression
    VectorMultivariateFunctionPatch<P> p1=factory(ep1._fp).create_zeros(as1,domain);
    VectorMultivariateFunctionPatch<P> p2=factory(ep2._fp).create_zeros(as2,domain);
    for (SizeType i=0; i!=as1; ++i) {
        auto v=ep1._vars[i]; p1[i]=ScalarMultivariateFunction<P>::coordinate(asr,indices[v]);
    }
    for (SizeType i=0; i!=as2; ++i) {
        auto v=ep2._vars[i]; p2[i]=ScalarMultivariateFunction<P>::coordinate(asr,indices[v]);
    }

    // Make projection functions from the common variables onto those used by each individual expression
    RealSpace rspc(variables);
    auto pf1=compose(ep1._fp,p1);
    auto pf2=compose(ep2._fp,p2);

    return std::make_tuple(rspc,pf1,pf2);
}

template<> inline String class_name<Vector<Real>>() { return "RealVector"; }

auto ExpressionPatch<ValidatedTag,RealVector>::_join(ExpressionPatch<P,T> ep1, ExpressionPatch<P,T> ep2) -> ExpressionPatch<P,T> {
    auto sff = make_common_variables(ep1,ep2); using std::get;
    auto const& dom=std::get<0>(sff);
    auto const& fp1=get<1>(sff);
    auto const& fp2=get<2>(sff);
    return ExpressionPatch<P,T>(dom,join(fp1,fp2));
}




auto ExpressionPatch<ValidatedTag,RealScalar>::_write(OutputStream& os) const -> OutputStream& {
    ExpressionPatch<P,T> const& ep=*this;
    os << "ExpressionPatch({";
    for (SizeType i=0; i!=ep._vars.size(); ++i) {
        if (i!=0) { os << ","; } os << ep._vars[i]<<"|"<<ep._fp.domain()[i];
    } os << std::flush;
    return os <<"}, f=" << static_cast<Function<P,SIG>>(ep._fp) << ")";
}

auto ExpressionPatch<ValidatedTag,RealVector>::_write(OutputStream& os) const -> OutputStream& {
    ExpressionPatch<P,T> const& ep=*this;
    os << "ExpressionPatch(" << this->_vars << ", " << this->_fp << "\n";
    os << "ExpressionPatch({";
    for (SizeType i=0; i!=ep._vars.size(); ++i) {
        if (i!=0) { os << ","; } os << ep._vars[i]<<"|"<<ep._fp.domain()[i];
    }
//    return os <<"}, f=" << static_cast<Function<P,SIG>>(ep._fp) << ")";
    Function<P,SIG> f=cast_unrestricted(ep._fp);
    return os <<"}," << f << ")";
}


template<class P, class... ARGS> auto FunctionPatch<P,RealScalar(ARGS...)>::operator() (InitializerList<AnyConstantOrVariablePatch> const& cvs) const -> ExpressionPatch<P,RealScalar> {
    List<ConstantOrVariablePatch<Real>> cvl=make_constant_or_variable_patch_list(cvs);
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealScalar(ARGS...)>::
operator() (ConstantOrVariable<ARGS> const& ... cvs) const -> ExpressionPatch<P,RealScalar> {
    List<ConstantOrVariable<Real>> cvl=make_constant_or_variable_list(cvs...);
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::operator() (InitializerList<AnyConstantOrVariablePatch> const& cvs) const -> ExpressionPatch<P,RealVector> {
    List<ConstantOrVariablePatch<Real>> cvl=make_constant_or_variable_patch_list(cvs);
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::
operator() (ConstantOrVariable<ARGS> const& ... cvs) const -> ExpressionPatch<P,RealVector> {
    List<ConstantOrVariable<Real>> cvl=make_constant_or_variable_list(cvs...);
    return make_expression_patch(*this,cvl);
}

template<class P, class... ARGS> auto FunctionPatch<P,RealScalar(ARGS...)>::
operator() (ExpressionPatch<P,ARGS...> const& ep) const -> ExpressionPatch<P,RealScalar> {
    static_assert (Same<Tuple<ARGS...>,Tuple<RealVector>>);
    return ExpressionPatch<P,RealScalar>(ep._vars,compose(*this,ep._fp));
}

template<class P, class... ARGS> auto FunctionPatch<P,RealVector(ARGS...)>::
operator() (ExpressionPatch<P,ARGS...> const& ep) const -> ExpressionPatch<P,RealVector> {
    static_assert (Same<Tuple<ARGS...>,Tuple<RealVector>>);
    return ExpressionPatch<P,RealVector>(ep._vars,compose(*this,ep._fp));
}

} // namespace Ariadne




#include "function/restricted_function.hpp"

namespace Ariadne {

#warning Move to FunctionPatch

template<class P, class... ARGS> FunctionPatch<P,RealScalar(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(RestrictedFunction<P,SIG>(f,dom)) { }
template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(RestrictedFunction<P,SIG>(f,dom)) { }

} // namespace Ariadne

/*
#include "function/taylor_function.hpp"
namespace Ariadne {
template<class P, class... ARGS> FunctionPatch<P,RealScalar(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(ValidatedScalarMultivariateTaylorFunctionModelDP(dom,f,ThresholdSweeper<FloatDP>(dp,1e-8))) { }
template<class P, class... ARGS> FunctionPatch<P,RealVector(ARGS...)>::FunctionPatch(const DomainType& dom, const Function<P,SIG>& f)
    : FunctionPatch(ValidatedVectorMultivariateTaylorFunctionModelDP(dom,f,ThresholdSweeper<FloatDP>(dp,1e-8))) { }
} // namespace Ariadne
*/

#endif /* ARIADNE_EXPRESSION_PATCH_HPP */
