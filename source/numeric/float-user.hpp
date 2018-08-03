/***************************************************************************
 *            float-user.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float.hpp
 *  \brief Inclusion header for floating-point numbers.
 */

#ifndef ARIADNE_FLOAT_USER_HPP
#define ARIADNE_FLOAT_USER_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float-raw.hpp"

#include "float_operations.hpp"

#include "float_approximation.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "float_bounds.hpp"
#include "float_ball.hpp"
#include "float_value.hpp"
#include "float_error.hpp"

#include "float_factory.hpp"

namespace Ariadne {

template<class F> inline Bounds<F>::Bounds(LowerBound<F> const& lower, UpperBound<F> const& upper)
    : _l(lower.raw()), _u(upper.raw()) { }

template<class F> inline LowerBound<F> const Bounds<F>::lower() const {
    return LowerBound<F>(lower_raw()); }

template<class F> inline UpperBound<F> const Bounds<F>::upper() const {
    return UpperBound<F>(upper_raw()); }

template<class F> inline const Value<F> Bounds<F>::value() const {
    return Value<F>(med(near,this->_l,this->_u)); }

template<class F> inline const Error<F> Bounds<F>::error() const {
    RawFloat<PR> _v=med(near,this->_l,this->_u); return Error<F>(max(sub(up,this->_u,_v),sub(up,_v,this->_l))); }

template<class F> inline Value<F> value(Bounds<F> const& x) {
    return x.value(); }

template<class F> inline Error<F> error(Bounds<F> const& x) {
    return x.error(); }

template<class F> inline Bounds<F> make_bounds(Error<F> const& e) {
    return Bounds<F>(-e.raw(),+e.raw()); }


template<class F, class FE> inline Ball<F,FE>::Ball(Value<F> const& value, Error<FE> const& error)
    : _v(value.raw()), _e(error.raw()) { }

template<class F, class FE> inline LowerBound<F> const Ball<F,FE>::lower() const {
    return LowerBound<F>(lower_raw()); }

template<class F, class FE> inline UpperBound<F> const Ball<F,FE>::upper() const {
    return UpperBound<F>(upper_raw()); }

template<class F, class FE> inline const Value<F> Ball<F,FE>::value() const {
    return Value<F>(this->_v); }

template<class F, class FE> inline const Error<FE> Ball<F,FE>::error() const {
    return Error<FE>(this->_e); }


template<class F> template<class FE> inline FloatBall<PrecisionType<F>,Ariadne::PrecisionType<FE>> Value<F>::pm(Error<FE> e) const {
    return FloatBall<Ariadne::PrecisionType<F>,Ariadne::PrecisionType<FE>>(*this,e);
}

template<class F> template<class FE> Approximation<F>::Approximation(Ball<F,FE> const& x) : _a(x.value_raw()) {
}

template<class F> template<class FE> Bounds<F>::Bounds(Ball<F,FE> const& x) : _l(x.lower_raw()), _u(x.upper_raw()) {
}


extern const FloatDPValue infty;

// Literals operation
FloatDPValue operator"" _exact(long double lx);
FloatDPError operator"" _error(long double lx);
FloatDPBall operator"" _near(long double lx);
FloatDPUpperBound operator"" _upper(long double lx);
FloatDPLowerBound operator"" _lower(long double lx);
FloatDPApproximation operator"" _approx(long double lx);



template<class F> inline Value<F> const& cast_exact(F const& x) { return reinterpret_cast<Value<F> const&>(x); }
template<class F> inline Value<F> const& cast_exact(Approximation<F> const& x) { return reinterpret_cast<Value<F> const&>(x); }
template<class F> inline Value<F> const& cast_exact(LowerBound<F> const& x) { return reinterpret_cast<Value<F> const&>(x); }
template<class F> inline Value<F> const& cast_exact(UpperBound<F> const& x) { return reinterpret_cast<Value<F> const&>(x); }
template<class F> inline Value<F> const cast_exact(Bounds<F> const& x) { return cast_exact(Approximation<F>(x)); }
template<class F, class FE> inline Value<F> const& cast_exact(Ball<F,FE> const& x) { return reinterpret_cast<Value<F> const&>(x); }
template<class F> inline Value<F> const& cast_exact(Value<F> const& x) { return reinterpret_cast<Value<F> const&>(x); }
template<class F> inline Value<F> const& cast_exact(Error<F> const& x) { return reinterpret_cast<Value<F> const&>(x); }


template<template<class>class T, class F> inline const T<Value<F>>& cast_exact(const T<F>& t) {
    return reinterpret_cast<const T<Value<F>>&>(t); }
template<template<class>class T, class F> inline const T<Value<F>>& cast_exact(const T<Approximation<F>>& t) {
    return reinterpret_cast<const T<Value<F>>&>(t); }
template<template<class>class T, class F> inline const T<Value<F>>& cast_exact(const T<LowerBound<F>>& t) {
    return reinterpret_cast<const T<Value<F>>&>(t); }
template<template<class>class T, class F> inline const T<Value<F>>& cast_exact(const T<UpperBound<F>>& t) {
    return reinterpret_cast<const T<Value<F>>&>(t); }
template<template<class>class T, class F> inline const T<Value<F>>& cast_exact(const T<Value<F>>& t) {
    return reinterpret_cast<const T<Value<F>>&>(t); }
template<template<class>class T, class F> inline const T<Value<F>>& cast_exact(const T<Error<F>>& t) {
    return reinterpret_cast<const T<Value<F>>&>(t); }


inline RawFloatDP const& cast_raw(RawFloatDP const& x) { return reinterpret_cast<RawFloatDP const&>(x); }
inline RawFloatDP const& cast_raw(FloatDPApproximation const& x) { return reinterpret_cast<RawFloatDP const&>(x); }
inline RawFloatDP const& cast_raw(FloatDPValue const& x) { return reinterpret_cast<RawFloatDP const&>(x); }

template<template<class>class T> inline const T<RawFloatDP>& cast_raw(const T<RawFloatDP>& t) {
    return reinterpret_cast<const T<RawFloatDP>&>(t); }
template<template<class>class T> inline const T<RawFloatDP>& cast_raw(const T<FloatDPApproximation>& t) {
    return reinterpret_cast<const T<RawFloatDP>&>(t); }
template<template<class>class T> inline const T<RawFloatDP>& cast_raw(const T<FloatDPValue>& t) {
    return reinterpret_cast<const T<RawFloatDP>&>(t); }

inline FloatDPApproximation const& cast_approximate(RawFloatDP const& x) { return reinterpret_cast<FloatDPApproximation const&>(x); }
inline FloatDPApproximation const& cast_approximate(FloatDPApproximation const& x) { return reinterpret_cast<FloatDPApproximation const&>(x); }
inline FloatDPApproximation const& cast_approximate(FloatDPValue const& x) { return reinterpret_cast<FloatDPApproximation const&>(x); }

template<template<class>class T> inline const T<FloatDPApproximation>& cast_approximate(const T<RawFloatDP>& t) {
    return reinterpret_cast<const T<FloatDPApproximation>&>(t); }
template<template<class>class T> inline const T<FloatDPApproximation>& cast_approximate(const T<FloatDPApproximation>& t) {
    return reinterpret_cast<const T<FloatDPApproximation>&>(t); }
template<template<class>class T> inline const T<FloatDPApproximation>& cast_approximate(const T<FloatDPValue>& t) {
    return reinterpret_cast<const T<FloatDPApproximation>&>(t); }

inline FloatMPValue const& cast_exact(RawFloatMP const& x) { return reinterpret_cast<FloatMPValue const&>(x); }
inline FloatMPValue const& cast_exact(FloatMPApproximation const& x) { return reinterpret_cast<FloatMPValue const&>(x); }


} // namespace Ariadne

#endif
