/***************************************************************************
 *            numeric/casts.hpp
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

/*! \file numeric/casts.hpp
 *  \brief Inclusion header for casting between different kinds of information.
 */

#ifndef ARIADNE_CASTS_HPP
#define ARIADNE_CASTS_HPP

#include "../utility/macros.hpp"

#include "number.decl.hpp"
#include "float.decl.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"

#include "float_approximation.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "float_bounds.hpp"
#include "float_ball.hpp"
#include "float_value.hpp"
#include "float_error.hpp"

namespace Ariadne {

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

template<class F> inline const Positive<Value<F>> cast_exact(const Positive<Bounds<F>>& t) {
    return Positive<Value<F>>(cast_exact(static_cast<Bounds<F>const&>(t))); }

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
