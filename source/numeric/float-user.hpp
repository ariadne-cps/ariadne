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

#include "utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"
#include "float64.hpp"
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

template<class PR> inline FloatBounds<PR>::FloatBounds(FloatLowerBound<PR> const& lower, FloatUpperBound<PR> const& upper)
    : _l(lower.raw()), _u(upper.raw()) { }

template<class PR> inline FloatLowerBound<PR> const FloatBounds<PR>::lower() const {
    return FloatLowerBound<PR>(lower_raw()); }

template<class PR> inline FloatUpperBound<PR> const FloatBounds<PR>::upper() const {
    return FloatUpperBound<PR>(upper_raw()); }

template<class PR> inline const FloatValue<PR> FloatBounds<PR>::value() const {
    return FloatValue<PR>(med(near,this->_l,this->_u)); }

template<class PR> inline const FloatError<PR> FloatBounds<PR>::error() const {
    RawFloat<PR> _v=med(near,this->_l,this->_u); return FloatError<PR>(max(sub(up,this->_u,_v),sub(up,_v,this->_l))); }

template<class PR> inline FloatValue<PR> value(FloatBounds<PR> const& x) {
    return x.value(); }

template<class PR> inline FloatError<PR> error(FloatBounds<PR> const& x) {
    return x.error(); }

template<class PR> inline FloatBounds<PR> make_bounds(FloatError<PR> const& e) {
    return FloatBounds<PR>(-e.raw(),+e.raw()); }


template<class PR, class PRE> inline FloatBall<PR,PRE>::FloatBall(FloatValue<PR> const& value, FloatError<PRE> const& error)
    : _v(value.raw()), _e(error.raw()) { }

template<class PR, class PRE> inline FloatLowerBound<PR> const FloatBall<PR,PRE>::lower() const {
    return FloatLowerBound<PR>(lower_raw()); }

template<class PR, class PRE> inline FloatUpperBound<PR> const FloatBall<PR,PRE>::upper() const {
    return FloatUpperBound<PR>(upper_raw()); }

template<class PR, class PRE> inline const FloatValue<PR> FloatBall<PR,PRE>::value() const {
    return FloatValue<PR>(this->_v); }

template<class PR, class PRE> inline const FloatError<PRE> FloatBall<PR,PRE>::error() const {
    return FloatError<PRE>(this->_e); }


template<class PR> template<class PRE> inline FloatBall<PR,PRE> FloatValue<PR>::pm(FloatError<PRE> e) const {
    return FloatBall<PR,PRE>(*this,e);
}

template<class PR> template<class PRE> FloatApproximation<PR>::FloatApproximation(FloatBall<PR,PRE> const& x) : _a(x.value_raw()) {
}

template<class PR> template<class PRE> FloatBounds<PR>::FloatBounds(FloatBall<PR,PRE> const& x) : _l(x.lower_raw()), _u(x.upper_raw()) {
}


extern const FloatDPValue infty;

// Literals operation
FloatDPValue operator"" _exact(long double lx);
FloatDPError operator"" _error(long double lx);
FloatDPBall operator"" _near(long double lx);
FloatDPUpperBound operator"" _upper(long double lx);
FloatDPLowerBound operator"" _lower(long double lx);
FloatDPApproximation operator"" _approx(long double lx);



inline FloatDPValue const& cast_exact(RawFloatDP const& x) { return reinterpret_cast<FloatDPValue const&>(x); }
inline FloatDPValue const& cast_exact(FloatDPApproximation const& x) { return reinterpret_cast<FloatDPValue const&>(x); }
inline FloatDPValue const& cast_exact(FloatDPValue const& x) { return reinterpret_cast<FloatDPValue const&>(x); }
inline FloatDPValue const& cast_exact(FloatDPError const& x) { return reinterpret_cast<FloatDPValue const&>(x); }

template<template<class>class T> inline const T<FloatDPValue>& cast_exact(const T<RawFloatDP>& t) {
    return reinterpret_cast<const T<FloatDPValue>&>(t); }
template<template<class>class T> inline const T<FloatDPValue>& cast_exact(const T<FloatDPApproximation>& t) {
    return reinterpret_cast<const T<FloatDPValue>&>(t); }
template<template<class>class T> inline const T<FloatDPValue>& cast_exact(const T<FloatDPValue>& t) {
    return reinterpret_cast<const T<FloatDPValue>&>(t); }
template<template<class>class T> inline const T<FloatDPValue>& cast_exact(const T<FloatDPError>& t) {
    return reinterpret_cast<const T<FloatDPValue>&>(t); }

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
