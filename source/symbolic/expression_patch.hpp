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

#include "utility/typedefs.hpp"

#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/function.hpp"
//#include "function/restricted_function.hpp"
#include "function/function_patch.hpp"

#include "symbolic/constant.hpp"
#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_set.hpp"

namespace Ariadne {

template<class RES> FunctionPatch(RealSpacePatch, Expression<RES>) -> FunctionPatch<ValidatedTag,RES(RealVector)>;
template<class P, class RES> FunctionPatch(RealSpace, ExpressionPatch<P,RES>) -> FunctionPatch<P,RES(RealVector)>;

template<class P, class R> class RestrictedExpression;
template<class R> using ValidatedRestrictedExpression = RestrictedExpression<ValidatedTag,R>;
using ValidatedScalarRestrictedExpression = RestrictedExpression<ValidatedTag,RealScalar>;
using ValidatedVectorRestrictedExpression = RestrictedExpression<ValidatedTag,RealVector>;

template<class P, class R> class ExpressionPatch;
template<class R> using ValidatedExpressionPatch = ExpressionPatch<ValidatedTag,R>;
using ValidatedScalarExpressionPatch = ExpressionPatch<ValidatedTag,RealScalar>;
using ValidatedVectorExpressionPatch = ExpressionPatch<ValidatedTag,RealVector>;

template<class T> class SpacePatch;
using RealSpacePatch = SpacePatch<Real>;

template<template<class>class SV> FunctionPatch(SpacePatch<Real>, RestrictedExpression<ValidatedTag,SV<Real>>) -> FunctionPatch<ValidatedTag,SV<Real>(Vector<Real>)>;


typedef VariableInterval<typename IntervalDomainType::UpperBoundType> VariableIntervalDomainType;
typedef VariablesBox<IntervalDomainType> VariablesBoxDomainType;
typedef VectorVariableBox<IntervalDomainType> VectorVariableBoxDomainType;

static_assert(not Convertible<VariablesBoxDomainType,VectorVariableBoxDomainType>);
static_assert(Convertible<VectorVariableBoxDomainType,VariablesBoxDomainType>);

template<class R> struct VariablePatchTypedef;
template<class R> using VariablePatchType = typename VariablePatchTypedef<R>::Type;
template<class R> using VariablePatch = typename VariablePatchTypedef<R>::Type;
template<> struct VariablePatchTypedef<RealScalar> { typedef VariableIntervalDomainType Type; };
template<> struct VariablePatchTypedef<RealVector> { typedef VariablesBoxDomainType Type; };
using RealVariablePatchType = VariablePatchType<Real>;
using RealVariablePatch = VariablePatchType<Real>;



template<class X> class ScalarOrVector
    : public Variant<Scalar<X>,Vector<X>>
{
  public:
    template<class S> requires Convertible<S,Scalar<X>> ScalarOrVector(S const& s)
        : Variant<Scalar<X>,Vector<X>>(Scalar<X>(s)) { }
    template<class V> requires Convertible<V,Vector<X>> ScalarOrVector(V const& v)
        : Variant<Scalar<X>,Vector<X>>(Scalar<X>(v)) { }
};


template<template<class>class T, class R> class OfScalarOrVector
    : public Variant<T<Scalar<R>>,T<Vector<R>>>
{
    using Base = Variant<T<Scalar<R>>,T<Vector<R>>>;
  public:
    template<class S> requires Convertible<S,T<Scalar<R>>> OfScalarOrVector(S const& s)
        : Base(T<Scalar<R>>(s)) { }
    template<class V> requires Convertible<V,T<Vector<R>>> OfScalarOrVector(V const& v)
        : Base(T<Vector<R>>(v)) { }
};

template<class R> Vector<R> join(List<ScalarOrVector<R>> const& lst);
template<template<class>class T, class R> T<Vector<R>> join(List<OfScalarOrVector<T,R>> const&);

template<class R> class ConstantOrVariable
    : public Variant<Constant<R>,Variable<R>>
{
    using Base = Variant<Constant<R>,Variable<R>>;
  public:
    template<class C> requires Convertible<C,R> ConstantOrVariable(C const& c)
        : Base(Constant<R>(R(c))) { }
    template<class C> requires Convertible<C,Constant<R>> and (not Convertible<C,R>) ConstantOrVariable(C const& c)
        : Base(Constant<R>(c)) { }
    template<class V> requires Convertible<V,Variable<R>> ConstantOrVariable(V const& v)
        : Base(Variable<R>(v)) { }
};

template<class R> class ConstantOrVariablePatch
    : public Variant<Constant<R>,VariablePatch<R>>
{
    using Base = Variant<Constant<R>,VariablePatch<R>>;
  public:
    template<class C> requires Convertible<C,R> ConstantOrVariablePatch(C const& c)
        : Base(Constant<R>(R(c))) { }
    template<class C> requires Convertible<C,Constant<R>> and (not Convertible<C,R>) ConstantOrVariablePatch(C const& c)
        : Base(Constant<R>(c)) { }
    template<class V> requires Convertible<V,VariablePatch<R>> ConstantOrVariablePatch(V const& v)
        : Base(Variable<R>(v)) { }
};

template<class R> class ScalarOrVectorConstantOrVariable
    : public Variant<Constant<Scalar<R>>,Constant<Vector<R>>,Variable<Scalar<R>>,Variable<Vector<R>>>
{
    using Base = Variant<Constant<Scalar<R>>,Constant<Vector<R>>,Variable<Scalar<R>>,Variable<Vector<R>>>;
  public:
    template<class S> requires Convertible<S,Scalar<R>> ScalarOrVectorConstantOrVariable(S const& s)
        : Base(Constant<Scalar<R>>(Scalar<R>(s))) { }
    template<class CS> requires Convertible<CS,Constant<Scalar<R>>> and (not Convertible<CS,Scalar<R>>)
    ScalarOrVectorConstantOrVariable(CS const& cs)
        : Base(Constant<Scalar<R>>(cs)) { }
    template<class CV> requires Convertible<CV,Constant<Vector<R>>> ScalarOrVectorConstantOrVariable(CV const& cv)
        : Base(Constant<Vector<R>>(cv)) { }
    template<class VS> requires Convertible<VS,Variable<Scalar<R>>> ScalarOrVectorConstantOrVariable(VS const& vs)
        : Base(Variable<Scalar<R>>(vs)) { }
    template<class VV> requires Convertible<VV,Variable<Vector<R>>> ScalarOrVectorConstantOrVariable(VV const& vv)
        : Base(Variable<Vector<R>>(vv)) { }
};

template<class P, class R> class ScalarOrVectorExpressionPatch
    : public Variant<ExpressionPatch<P,Scalar<R>>,ExpressionPatch<P,Vector<R>>>
{
    using Base = Variant<ExpressionPatch<P,Scalar<R>>,ExpressionPatch<P,Vector<R>>>;
  public:
    template<class S> requires Convertible<S,Scalar<R>> ScalarOrVectorExpressionPatch(S const& s)
        : Base(ExpressionPatch<P,Scalar<R>>(Constant<Scalar<R>>(Scalar<R>(s)))) { }
    template<class CS> requires Convertible<CS,Constant<Scalar<R>>> and (not Convertible<CS,Scalar<R>>)
    ScalarOrVectorExpressionPatch(CS const& cs)
        : Base(ExpressionPatch<P,Scalar<R>>(Constant<Scalar<R>>(cs))) { }
    template<class CV> requires Convertible<CV,Constant<Vector<R>>> ScalarOrVectorExpressionPatch(CV const& cv)
        : Base(Constant<Vector<R>>(cv)) { }
    template<class VS> requires Convertible<VS,VariablePatch<Scalar<R>>> ScalarOrVectorExpressionPatch(VS const& vs)
        : Base(ExpressionPatch<P,Scalar<R>>(VariablePatch<Scalar<R>>(vs))) { }
    template<class VV> requires Convertible<VV,VariablePatch<Vector<R>>> ScalarOrVectorExpressionPatch(VV const& vv)
        : Base(ExpressionPatch<P,Vector<R>>(VariablePatch<Vector<R>>(vv))) { }
    ScalarOrVectorExpressionPatch(ExpressionPatch<P,Scalar<R>> const& es)
        : Base(es) { }
    ScalarOrVectorExpressionPatch(ExpressionPatch<P,Vector<R>> const& ev)
        : Base(ev) { }
};
using ValidatedRealScalarOrVectorExpressionPatch = ScalarOrVectorExpressionPatch<ValidatedTag,Real>;


template<class R> class ConstantOrVariable;
template<class R> class ConstantOrVariablePatch;
template<class R> class ConstantOrVariablePatchOrRestrictedExpression;

template<class R> class ScalarOrVectorConstantOrVariablePatch;


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

template<class R> class VariableDomainMap;
using RealVariableDomainMap = VariableDomainMap<Real>;
template<> class VariableDomainMap<Real>
    : public Map<RealVariable,IntervalDomainType>
{
    using R=Real;
  public:
    VariableDomainMap();
    VariableDomainMap(InitializerList<VariableDomainMap<R>> lst);
    VariableDomainMap(VariableIntervalDomainType vivl);
    VariableDomainMap(VectorVariableBoxDomainType vvbx);
    VariableDomainMap(Map<RealVariable,IntervalDomainType>);
};

template<> class RestrictedExpression<ValidatedTag, RealScalar> {
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

    template<class C> requires Convertible<C,Real> RestrictedExpression(C const& c)
        : RestrictedExpression(Constant<Real>(Real(c))) { }
    RestrictedExpression(Constant<Real> const&);
    RestrictedExpression(VariablePatchType<Real> const&);

    explicit RestrictedExpression(RealVariableDomainMap doms, RealExpression expr);
    explicit RestrictedExpression(Vector<RealVariablePatch> ptchs, Function<P,SIG> f);
    explicit RestrictedExpression(Vector<RealVariable> vars, BoxDomainType dom, Function<P,SIG> f);
    explicit RestrictedExpression(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp);
//    explicit RestrictedExpression(RealSpacePatch spcptch, Function<P,SIG> f);
//    explicit RestrictedExpression(RealSpace spc, FunctionPatch<P,SIG> fp);
  public:
    friend auto operator+(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> { return pos(ep); }
    friend auto operator-(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> { return neg(ep); }

    friend auto operator+(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return add(ep1,ep2); }
    friend auto operator-(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return sub(ep1,ep2); }
    friend auto operator*(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return mul(ep1,ep2); }
    friend auto operator/(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return div(ep1,ep2); }
    friend auto operator+(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> { return add(ep1,c2); }
    friend auto operator-(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> { return sub(ep1,c2); }
    friend auto operator*(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> { return mul(ep1,c2); }
    friend auto operator/(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> { return div(ep1,c2); }
    friend auto operator+(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return add(c1,ep2); }
    friend auto operator-(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return sub(c1,ep2); }
    friend auto operator*(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return mul(c1,ep2); }
    friend auto operator/(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> { return div(c1,ep2); }

    friend auto operator*(RestrictedExpression<P,Scalar<Real>> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,Vector<Real>>;
    friend auto operator*(Vector<Number<P>> c1, RestrictedExpression<P,Scalar<Real>> ep2) -> RestrictedExpression<P,Vector<Real>>;
    friend auto operator/(Vector<Number<P>> c1, RestrictedExpression<P,Scalar<Real>> ep2) -> RestrictedExpression<P,Vector<Real>>;

    friend auto add(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Add(),ep1,ep2); }
    friend auto sub(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Sub(),ep1,ep2); }
    friend auto mul(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Mul(),ep1,ep2); }
    friend auto div(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Div(),ep1,ep2); }
    friend auto add(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Add(),std::placeholders::_1,c2),ep1); }
    friend auto sub(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Sub(),std::placeholders::_1,c2),ep1); }
    friend auto mul(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Mul(),std::placeholders::_1,c2),ep1); }
    friend auto div(RestrictedExpression<P,T> ep1, Number<P> c2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Div(),std::placeholders::_1,c2),ep1); }
    friend auto add(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Add(),c1,std::placeholders::_1),ep2); }
    friend auto sub(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Sub(),c1,std::placeholders::_1),ep2); }
    friend auto mul(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Mul(),c1,std::placeholders::_1),ep2); }
    friend auto div(Number<P> c1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Div(),c1,std::placeholders::_1),ep2); }

    friend auto nul(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Nul(),ep); }
    friend auto pos(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Pos(),ep); }
    friend auto neg(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Neg(),ep); }
    friend auto sqr(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Sqr(),ep); }
    friend auto hlf(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Hlf(),ep); }
    friend auto rec(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Rec(),ep); }
    friend auto pow(RestrictedExpression<P,T> ep, Int n) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(std::bind(Pow(),std::placeholders::_1,n),ep); }
    friend auto sqrt(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Sqrt(),ep); }
    friend auto exp(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Exp(),ep); }
    friend auto log(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Log(),ep); }
    friend auto sin(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Sin(),ep); }
    friend auto cos(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Cos(),ep); }
    friend auto tan(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Tan(),ep); }
    friend auto atan(RestrictedExpression<P,T> ep) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Atan(),ep); }

    friend auto mul(RestrictedExpression<P,RealScalar> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,RealVector>;
    friend auto mul(Vector<Number<P>> c1, RestrictedExpression<P,RealScalar> ep2) -> RestrictedExpression<P,RealVector>;
    friend auto div(Vector<Number<P>> c1, RestrictedExpression<P,RealScalar> ep2) -> RestrictedExpression<P,RealVector>;
    friend OutputStream& operator<<(OutputStream& os, RestrictedExpression<P,T> const& ep) {
        return ep._write(os); }
  private:
    OutputStream& _write(OutputStream& os) const;
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,T> const& ep1, RestrictedExpression<P,T> const& ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,T> const& ep);

};


template<> class RestrictedExpression<ValidatedTag, RealVector> {
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

    template<class C> requires Convertible<C,Constant<RealVector>> RestrictedExpression(C const& c)
        : RestrictedExpression(Constant<Vector<Real>>(c)) { }
    RestrictedExpression(Constant<RealVector> const& v);
    RestrictedExpression(VariablePatchType<RealVector> const& v);

    RestrictedExpression(InitializerList<ValidatedRealScalarOrVectorExpressionPatch> lst);


    explicit RestrictedExpression(Vector<RealVariable> vars, FunctionPatch<P,SIG> fp);
    explicit RestrictedExpression(Vector<RealVariable> vars, BoxDomainType dom, Function<P,SIG> f);
    explicit RestrictedExpression(Vector<RealVariablePatch> ptchs, Function<P,SIG> f);
    explicit RestrictedExpression(RealVariableDomainMap doms, Vector<RealExpression> expr);
//    explicit RestrictedExpression(RealSpacePatch spcptch, Function<P,SIG> f);
//    explicit RestrictedExpression(RealSpace spc, FunctionPatch<P,SIG> fp);

    friend auto call(FunctionPatch<P,RealVector(RealVector)>, RestrictedExpression<P,RealVector>) -> RestrictedExpression<P,RealVector>;

    friend auto nul(RestrictedExpression<P,VT> ep) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Nul(),ep); }
    friend auto pos(RestrictedExpression<P,VT> ep) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Pos(),ep); }
    friend auto neg(RestrictedExpression<P,VT> ep) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Neg(),ep); }

    friend auto add(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Add(),ep1,ep2); }
    friend auto sub(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Sub(),ep1,ep2); }
    friend auto mul(RestrictedExpression<P,ST> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),ep1,ep2); }
    friend auto mul(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),ep1,ep2); }
    friend auto div(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_apply(Div(),ep1,ep2); }

    friend auto add(RestrictedExpression<P,VT> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Add(),ep1,c2); }
    friend auto sub(RestrictedExpression<P,VT> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Sub(),ep1,c2); }
    friend auto mul(RestrictedExpression<P,ST> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),ep1,c2); }
    friend auto mul(RestrictedExpression<P,VT> ep1, Scalar<Number<P>> c2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),ep1,c2); }
    friend auto div(RestrictedExpression<P,VT> ep1, Scalar<Number<P>> c2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Div(),ep1,c2); }
    friend auto add(Vector<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Add(),c1,ep2); }
    friend auto sub(Vector<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Sub(),c1,ep2); }
    friend auto mul(Scalar<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),c1,ep2); }
    friend auto mul(Vector<Number<P>> c1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),c1,ep2); }
    friend auto div(Vector<Number<P>> c1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Div(),c1,ep2); }

    friend auto add(RestrictedExpression<P,VT> ep1, Expression<VT> e2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Add(),ep1,e2); }
    friend auto sub(RestrictedExpression<P,VT> ep1, Expression<VT> e2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Sub(),ep1,e2); }
    friend auto mul(RestrictedExpression<P,ST> ep1, Expression<VT> e2) -> RestrictedExpression<P,VT>;
    friend auto mul(RestrictedExpression<P,VT> ep1, Expression<ST> e2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),ep1,e2); }
    friend auto div(RestrictedExpression<P,VT> ep1, Expression<ST> e2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Div(),ep1,e2); }
    friend auto add(Expression<VT> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Add(),e1,ep2); }
    friend auto sub(Expression<VT> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Sub(),e1,ep2); }
    friend auto mul(Expression<ST> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
        return RestrictedExpression<P,T>::_apply(Mul(),e1,ep2); }
    friend auto mul(Expression<VT> e1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT>;
    friend auto div(Expression<VT> e1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT>;

    friend auto operator+(RestrictedExpression<P,VT> ep) -> RestrictedExpression<P,T> { return pos(ep); }
    friend auto operator-(RestrictedExpression<P,VT> ep) -> RestrictedExpression<P,T> { return neg(ep); }
    friend auto operator+(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return add(ep1,ep2); }
    friend auto operator-(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return sub(ep1,ep2); }
    friend auto operator*(RestrictedExpression<P,ST> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return mul(ep1,ep2); }
    friend auto operator*(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T> { return mul(ep1,ep2); }
    friend auto operator/(RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T> { return div(ep1,ep2); }

    friend auto operator+(RestrictedExpression<P,VT> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,T> { return add(ep1,c2); }
    friend auto operator-(RestrictedExpression<P,VT> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,T> { return sub(ep1,c2); }
    friend auto operator*(RestrictedExpression<P,ST> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,T> { return _apply(Mul(),ep1,c2); }
    friend auto operator*(RestrictedExpression<P,VT> ep1, Scalar<Number<P>> c2) -> RestrictedExpression<P,T> { return _apply(Mul(),ep1,c2); }
    friend auto operator/(RestrictedExpression<P,VT> ep1, Scalar<Number<P>> c2) -> RestrictedExpression<P,T> { return _apply(Div(),ep1,c2); }
    friend auto operator+(Vector<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return add(c1,ep2); }
    friend auto operator-(Vector<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return sub(c1,ep2); }
    friend auto operator*(Scalar<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return _apply(Mul(),c1,ep2); }
    friend auto operator*(Vector<Number<P>> c1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T> { return _apply(Mul(),c1,ep2); }
    friend auto operator/(Vector<Number<P>> c1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T> { return _apply(Div(),c1,ep2); }

    friend auto operator+(RestrictedExpression<P,VT> ep1, Expression<VT> e2) -> RestrictedExpression<P,T> { return add(ep1,e2); }
    friend auto operator-(RestrictedExpression<P,VT> ep1, Expression<VT> e2) -> RestrictedExpression<P,T> { return sub(ep1,e2); }
    friend auto operator*(RestrictedExpression<P,ST> ep1, Expression<VT> e2) -> RestrictedExpression<P,T>;
    friend auto operator*(RestrictedExpression<P,VT> ep1, Expression<ST> e2) -> RestrictedExpression<P,T> { return mul(ep1,e2); }
    friend auto operator/(RestrictedExpression<P,VT> ep1, Expression<ST> e2) -> RestrictedExpression<P,T> { return div(ep1,e2); }
    friend auto operator+(Expression<VT> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return add(e1,ep2); }
    friend auto operator-(Expression<VT> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return sub(e1,ep2); }
    friend auto operator*(Expression<ST> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,T> { return mul(e1,ep2); }
    friend auto operator*(Expression<VT> e1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T>;
    friend auto operator/(Expression<VT> e1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,T>;

    friend auto join(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2) -> RestrictedExpression<P,T> {
        return RestrictedExpression<P,T>::_join(ep1,ep2); }

    friend OutputStream& operator<<(OutputStream& os, RestrictedExpression<P,T> const& ep) {
        return ep._write(os); }

#warning
//    friend RestrictedExpression<ValidatedTag,RealVector> FunctionPatch<ValidatedTag,RealVector(RealVector)>::operator() (RestrictedExpression<ValidatedTag,RealVector> const& ep) const;

  private:
    static RestrictedExpression<P,T> _join(RestrictedExpression<P,T> ep1, RestrictedExpression<P,T> ep2);
    template<class OP> static RestrictedExpression<P,VT> _apply(OP op, RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2);
    template<class OP> static RestrictedExpression<P,VT> _apply(OP op, RestrictedExpression<P,ST> ep1, RestrictedExpression<P,VT> ep2);
    template<class OP> static RestrictedExpression<P,VT> _apply(OP op, RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,VT> ep1, Vector<Number<P>> c2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,ST> ep1, Vector<Number<P>> c2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,VT> ep1, Scalar<Number<P>> c2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, Vector<Number<P>> c1, RestrictedExpression<P,VT> ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, Scalar<Number<P>> c1, RestrictedExpression<P,VT> ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, Vector<Number<P>> c1, RestrictedExpression<P,ST> ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,VT> ep1, Expression<VT> e2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,VT> ep1, Expression<ST> e2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, Expression<VT> e1, RestrictedExpression<P,VT> ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, Expression<ST> e1, RestrictedExpression<P,VT> ep2);
    template<class OP> static RestrictedExpression<P,T> _apply(OP op, RestrictedExpression<P,T> ep);

    OutputStream& _write(OutputStream& os) const;
};

template<class P> class ExpressionPatch<P,RealScalar>
    : public RestrictedExpression<P,RealScalar>
{
    using R=RealScalar;
  public:
    using RestrictedExpression<P,R>::RestrictedExpression;
};


template<class P> class ExpressionPatch<P,RealVector>
    : public RestrictedExpression<P,RealVector>
{
    using R=RealVector;
  public:
    using RestrictedExpression<P,R>::RestrictedExpression;



};

#warning
/*
template<class... TS> class EasyConversionVariant;

template<class T1, class T2> class EasyConversionVariant<T1,T2>
    : public Variant<T1,T2>
{
  public:
    template<class A1> requires Convertible<A1,T1> EasyConversionVariant(A1 const& a1) : Variant<T1,T2>(T1(a1)) { }
    template<class A2> requires Convertible<A2,T2> EasyConversionVariant(A2 const& a2) : Variant<T2,T2>(T2(a2)) { }
};

template<class T1, class T2, class T3> class EasyConversionVariant<T1,T2,T3>
    : public Variant<T1,T2,T3>
{
  public:
    template<class A1> requires Convertible<A1,T1> EasyConversionVariant(A1 const& a1) : Variant<T1,T2,T3>(T1(a1)) { }
    template<class A2> requires Convertible<A2,T2> EasyConversionVariant(A2 const& a2) : Variant<T1,T2,T3>(T2(a2)) { }
    template<class A3> requires Convertible<A3,T3> EasyConversionVariant(A3 const& a3) : Variant<T1,T1,T3>(T3(a3)) { }
};


template<class R> class ConstantOrVariable : public EasyConversionVariant<Constant<R>,Variable<R>> {
    using EasyConversionVariant<Constant<R>,Variable<R>>::EasyConversionVariant;
};

template<class R> class ConstantOrVariablePatch : public EasyConversionVariant<Constant<R>,VariablePatchType<R>> {
    using EasyConversionVariant<Constant<R>,VariablePatchType<R>>::EasyConversionVariant;
};
*/

template<class T> concept IsConstantOrVariable = OneOf<T,Dyadic,RealConstant,RealVariable,RealVectorConstant,RealVectorVariable>;

template<class... FS> struct AllConstantOrVariable;
template<> struct AllConstantOrVariable<> : True { };
template<class F0, class... FS> struct AllConstantOrVariable<F0,FS...> {
    static const bool value = IsConstantOrVariable<F0> && AllConstantOrVariable<FS...>::value;
};

template<class... FS> concept AllAreConstantOrVariable = AllConstantOrVariable<FS...>::value;




/*
template<> class ConstantOrVariable<Vector<Real>> : public Vector<ConstantOrVariable<Real>> {
  public:
    using Vector<ConstantOrVariable<Real>>::Vector;
    template<class... AS> requires AllAreConstantOrVariable<AS...>
        ConstantOrVariable(AS... as) : ConstantOrVariable(_make_list(as...)) { }
  private:
    ConstantOrVariable(Constant<Vector<Real>>);
    ConstantOrVariable(Variable<Vector<Real>>);
};

template<class... AS> List<ConstantOrVariable<Real>> make_constant_or_variable_list(ConstantOrVariable<AS> const& ... cvs) {
    return ConstantOrVariable<Vector<Real>>::_make_list(cvs...);
}
*/


} // namespace Ariadne

#include "expression_patch.inl.hpp"

#endif /* ARIADNE_EXPRESSION_PATCH_HPP */
