/***************************************************************************
 *            float-user.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file float.h
 *  \brief Inclusion header for floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_USER_H
#define ARIADNE_FLOAT_USER_H

#include "utility/macros.h"

#include "number.decl.h"
#include "float.decl.h"
#include "float64.h"
#include "floatmp.h"
#include "float-raw.h"

namespace Ariadne {



template<class X, class R=X> struct DeclareFloatOperations
    : public DeclareRealOperations<X,R>
    , public DeclareComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , public DeclareMixedFieldOperators<X,GenericTrait<X>,R>
{
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class X, class NX=OppositeTrait<X>> struct DeclareDirectedFloatOperations
    : DeclareDirectedNumericOperations<X,NX>
    , DeclareMixedDirectedGroupOperators<X,NX,GenericTrait<X>,GenericTrait<NX>>
{
    friend OutputStream& operator<<(OutputStream&, X const&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class PX, class QPX=OppositeTrait<PX>> struct DeclarePositiveDirectedFloatOperations
    : DeclarePositiveDirectedNumericOperations<PX,QPX>
    , DeclareMixedDirectedSemifieldOperators<PX,QPX,GenericTrait<PX>,GenericTrait<QPX>>
{
    friend OutputStream& operator<<(OutputStream&, PX const&);
    friend InputStream& operator>>(InputStream&, PX&);
};

template<class PX> struct DeclarePositiveFloatOperations
    : DeclarePositiveDirectedFloatOperations<PX,PX>
{
};

template<class X, class R=X> struct DispatchFloatOperations
    : DispatchNumericOperations<X,R>
    , DispatchComparisonOperations<X,LessTrait<X>,EqualsTrait<X>>
    , DefineConcreteGenericOperators<X>
//    , ProvideConcreteGenericFieldOperations<X,GenericTrait<X>,R>
//    , ProvideConcreteGenericComparisonOperations<X,GenericTrait<X>,LessTrait<X>,EqualsTrait<X>>
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class X> struct DispatchDirectedFloatOperations
    : DispatchDirectedNumericOperations<X,OppositeTrait<X>>
    , DispatchDirectedNumericOperations<OppositeTrait<X>,X>
    , DispatchDirectedComparisonOperations<X,OppositeTrait<X>,LessTrait<X>,EqualsTrait<X>>
    , DispatchDirectedComparisonOperations<OppositeTrait<X>,X,LessTrait<OppositeTrait<X>>,EqualsTrait<OppositeTrait<X>>>
    , DefineConcreteGenericOperators<X>
//    , ProvideConcreteGenericDirectedGroupOperations<X,OppositeTrait<X>,GenericTrait<X>,GenericTrait<OppositeTrait<X>>>
//    , ProvideConcreteGenericDirectedGroupOperations<OppositeTrait<X>,X,GenericTrait<OppositeTrait<X>>,GenericTrait<X>>
//    , ProvideConcreteGenericDirectedComparisonOperations<X,GenericTrait<OppositeTrait<X>>,LessTrait<X>,EqualsTrait<X>>
//    , ProvideConcreteGenericDirectedComparisonOperations<OppositeTrait<X>,GenericTrait<X>,LessTrait<OppositeTrait<X>>,EqualsTrait<OppositeTrait<X>>>
{
    friend OutputStream& operator<<(OutputStream& os, X const& x) { return Operations<X>::_write(os,x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is,x); }
};

template<class PX, class QPX=OppositeTrait<PX>> struct DispatchPositiveDirectedFloatOperations
    : public DispatchPositiveDirectedNumericOperations<PX,QPX>
    , public DispatchPositiveDirectedNumericOperations<QPX,PX>
    , DefineConcreteGenericOperators<PX>
//    , public ProvideConcreteGenericDirectedSemiFieldOperations<PX,QPX,Nat,Nat>
//    , public ProvideConcreteGenericDirectedSemiFieldOperations<QPX,PX,Nat,Nat>
{
};

template<class PX> struct DispatchPositiveFloatOperations
    : public DispatchPositiveDirectedNumericOperations<PX,PX>
    , DefineConcreteGenericOperators<PX>
//    , public ProvideConcreteGenericDirectedSemiFieldOperations<PX,PX,Nat,Nat>
{
};


template<class PR, class P1, class P2> using FloatWeakerType = Float<Weaker<P1,P2>,PR>;

template<class PR, class P> using NegatedFloatType = Float<Negated<P>,PR>;
template<class PR, class P> using FloatNegateType = Float<Negated<P>,PR>;

template<class PR, class P1, class P2> using FloatSumType = Float<Widen<Weaker<P1,P2>>,PR>;
template<class PR, class P1, class P2> using FloatDifferenceType = Float<Widen<Weaker<P1,Negated<P2>>>,PR>;
template<class PR, class P1, class P2> using FloatProductType = Float<Widen<Weaker<P1,P2>>,PR>;
template<class PR, class P1, class P2> using FloatQuotientType = Float<Widen<Weaker<P1,Inverted<P2>>>,PR>;

template<class PR, class P1, class P2> using FloatEqualsType = Logical<Equality<Weaker<P1,Negated<P2>>>>;
template<class PR, class P1, class P2> using FloatLessType = Logical<Generic<Weaker<P1,Negated<P2>>>>;

} // namespace Ariadne

#include "float_approximation.h"
#include "float_lower_bound.h"
#include "float_upper_bound.h"
#include "float_bounds.h"
#include "float_ball.h"
#include "float_value.h"
#include "float_error.h"

#include "float_factory.h"

namespace Ariadne {

template<class PR> inline FloatBounds<PR>::FloatBounds(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper)
    : _l(lower.raw()), _u(upper.raw()) { }

template<class PR> inline FloatLowerBound<PR> const FloatBounds<PR>::lower() const {
    return FloatLowerBound<PR>(lower_raw()); }

template<class PR> inline FloatUpperBound<PR> const FloatBounds<PR>::upper() const {
    return FloatUpperBound<PR>(upper_raw()); }

template<class PR> inline const FloatValue<PR> FloatBounds<PR>::value() const {
    return FloatValue<PR>(med_near(this->_l,this->_u)); }

template<class PR> inline const FloatError<PR> FloatBounds<PR>::error() const {
    RawFloat<PR> _v=med_near(this->_l,this->_u); return FloatError<PR>(max(sub_up(this->_u,_v),sub_up(_v,this->_l))); }

template<class PR> inline FloatValue<PR> value(FloatBounds<PR> const& x) {
    return x.value(); }

template<class PR> inline FloatError<PR> error(FloatBounds<PR> const& x) {
    return x.error(); }

template<class PR> inline FloatBounds<PR> make_bounds(FloatError<PR> const& e) {
    return FloatBounds<PR>(-e.raw(),+e.raw()); }


template<class PR> inline FloatBall<PR>::FloatBall(FloatValue<PR> const& value, FloatError<PR> const& error)
    : _v(value.raw()), _e(error.raw()) { }

template<class PR> inline FloatLowerBound<PR> const FloatBall<PR>::lower() const {
    return FloatLowerBound<PR>(lower_raw()); }

template<class PR> inline FloatUpperBound<PR> const FloatBall<PR>::upper() const {
    return FloatUpperBound<PR>(upper_raw()); }

template<class PR> inline const FloatValue<PR> FloatBall<PR>::value() const {
    return FloatValue<PR>(this->_v); }

template<class PR> inline const FloatError<PR> FloatBall<PR>::error() const {
    return FloatError<PR>(this->_e); }

extern const Float64Value infty;

// Literals operation
Float64Value operator"" _exact(long double lx);
Float64Error operator"" _error(long double lx);
Float64Ball operator"" _near(long double lx);
Float64UpperBound operator"" _upper(long double lx);
Float64LowerBound operator"" _lower(long double lx);
Float64Approximation operator"" _approx(long double lx);



inline Float64Value const& cast_exact(RawFloat64 const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Approximation const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Value const& x) { return reinterpret_cast<Float64Value const&>(x); }
inline Float64Value const& cast_exact(Float64Error const& x) { return reinterpret_cast<Float64Value const&>(x); }

template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }
template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<Float64Approximation>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }
template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<Float64Value>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }
template<template<class>class T> inline const T<Float64Value>& cast_exact(const T<Float64Error>& t) {
    return reinterpret_cast<const T<Float64Value>&>(t); }

inline RawFloat64 const& cast_raw(RawFloat64 const& x) { return reinterpret_cast<RawFloat64 const&>(x); }
inline RawFloat64 const& cast_raw(Float64Approximation const& x) { return reinterpret_cast<RawFloat64 const&>(x); }
inline RawFloat64 const& cast_raw(Float64Value const& x) { return reinterpret_cast<RawFloat64 const&>(x); }

template<template<class>class T> inline const T<RawFloat64>& cast_raw(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }
template<template<class>class T> inline const T<RawFloat64>& cast_raw(const T<Float64Approximation>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }
template<template<class>class T> inline const T<RawFloat64>& cast_raw(const T<Float64Value>& t) {
    return reinterpret_cast<const T<RawFloat64>&>(t); }

inline Float64Approximation const& cast_approximate(RawFloat64 const& x) { return reinterpret_cast<Float64Approximation const&>(x); }
inline Float64Approximation const& cast_approximate(Float64Approximation const& x) { return reinterpret_cast<Float64Approximation const&>(x); }
inline Float64Approximation const& cast_approximate(Float64Value const& x) { return reinterpret_cast<Float64Approximation const&>(x); }

template<template<class>class T> inline const T<Float64Approximation>& cast_approximate(const T<RawFloat64>& t) {
    return reinterpret_cast<const T<Float64Approximation>&>(t); }
template<template<class>class T> inline const T<Float64Approximation>& cast_approximate(const T<Float64Approximation>& t) {
    return reinterpret_cast<const T<Float64Approximation>&>(t); }
template<template<class>class T> inline const T<Float64Approximation>& cast_approximate(const T<Float64Value>& t) {
    return reinterpret_cast<const T<Float64Approximation>&>(t); }

inline FloatMPValue const& cast_exact(RawFloatMP const& x) { return reinterpret_cast<FloatMPValue const&>(x); }
inline FloatMPValue const& cast_exact(FloatMPApproximation const& x) { return reinterpret_cast<FloatMPValue const&>(x); }


} // namespace Ariadne

#endif
