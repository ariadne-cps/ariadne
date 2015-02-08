/***************************************************************************
 *            numeric.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file numeric.h
 *  \brief Number classes. File suitable for use as a pre-compiled header.
 */

#ifndef ARIADNE_NUMERIC_H
#define ARIADNE_NUMERIC_H

#include "utility/standard.h"

#include "config.h"

#include "utility/declarations.h"

#include "numeric/logical.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "numeric/float.h"
#include "numeric/real.h"

namespace Ariadne {

template<class X> inline X generic_pow(X p, Nat m) {
    X r=p.create_constant(1); while(m!=0) { if(m%2==1) { r=r*p; } p=sqr(p); m=m/2; } return r; }

template<class X> inline X generic_pow(X p, Int n) {
    return n>=0 ? generic_pow(p,Nat(n)) : rec(generic_pow(p,Nat(-n))); }

//! \ingroup NumericModule
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline Int numeric_cast(const Float64& a) { return Int(a.get_d()); }
template<> inline double numeric_cast(const Float64& a) { return a.get_d(); }
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline double numeric_cast(const ExactFloat64& a) { return a.get_d(); }
template<> inline double numeric_cast(const ValidatedFloat64& a) { return a.get_d(); }
template<> inline double numeric_cast(const ApproximateFloat64& a) { return a.get_d(); }
template<> inline float numeric_cast(const double& a) { return a; }
template<> inline float numeric_cast(const Float64& a) { return a.get_d(); }
template<> inline float numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Float64 numeric_cast(const ExactFloat64& a) { return a.raw(); }

template<> inline Real numeric_cast(const Float64& a) { return Real(ExactFloat64(a)); }
template<> inline Real numeric_cast(const ExactFloat64& a) { return Real(a); }
template<> inline Real numeric_cast(const ValidatedFloat64& a) { return Real(make_exact(ApproximateFloat64(a))); }

template<class R, class A> inline R approx_cast(const A& a);
template<> inline double approx_cast(const Float64& a) { return a.get_d(); }
template<> inline double approx_cast(const ExactFloat64& a) { return a.get_d(); }

//! \ingroup NumericModule
//! \related Float64
//! \brief Converts \a e to an object of type \a X, which may either be an
//! \c Float64 or \c ExactInterval, with the semantics that \a e denotes and error bound.
//! Returns the float 0.0 (since floating-point computations do not keep track of errors)
//! and the interval [-e,+e].
template<class X> inline X convert_error(const Float64& e);
template<> inline Float64 convert_error<Float64>(const Float64& e) { return Float64(0.0); }

template<class X> inline X convert_error(const ApproximateFloat64& e);
template<> inline ApproximateFloat64 convert_error<ApproximateFloat64>(const ApproximateFloat64& e) { return ApproximateFloat64(0.0); }
template<> inline ValidatedFloat64 convert_error<ValidatedFloat64>(const ApproximateFloat64& e) { return ValidatedFloat64(-e.raw(),+e.raw()); }

// Use 'enable_if' style template to restrict allowable instances. See the Boost documentation
// for enable_if to see how this works.
template<class X> struct IsNumeric { static const Bool value = false; };
template<> struct IsNumeric<unsigned int> { static const Bool value = true; };
template<> struct IsNumeric<int> { static const Bool value = true; };
template<> struct IsNumeric<double> { static const Bool value = true; };
template<> struct IsNumeric<Float64> { static const Bool value = true; };
template<> struct IsNumeric<ExactFloat64> { static const Bool value = true; };
template<> struct IsNumeric<ValidatedFloat64> { static const Bool value = true; };
template<> struct IsNumeric<ApproximateFloat64> { static const Bool value = true; };
template<> struct IsNumeric<Real> { static const Bool value = true; };

template<> struct IsNumeric<Integer> { static const Bool value = true; };
template<> struct IsNumeric<Rational> { static const Bool value = true; };

template<class X, class T> using EnableIfNumeric = EnableIf<IsNumeric<X>,T>;

template<class F, class T> struct IsSafelyConvertible : False { };
template<> struct IsSafelyConvertible<ExactFloat64,ExactFloat64> : True { };
template<> struct IsSafelyConvertible<ExactFloat64,Real> : True { };
template<> struct IsSafelyConvertible<ExactFloat64,ValidatedFloat64> : True { };
template<> struct IsSafelyConvertible<ExactFloat64,ApproximateFloat64> : True { };
template<> struct IsSafelyConvertible<Real,Real> : True { };
template<> struct IsSafelyConvertible<Real,ValidatedFloat64> : True { };
template<> struct IsSafelyConvertible<Real,ApproximateFloat64> : True { };
template<> struct IsSafelyConvertible<ValidatedFloat64,ValidatedFloat64> : True { };
template<> struct IsSafelyConvertible<ValidatedFloat64,ApproximateFloat64> : True { };
template<> struct IsSafelyConvertible<ApproximateFloat64,ApproximateFloat64> : True { };

template<> struct IsSafelyConvertible<Real,UpperInterval> : True { };
template<> struct IsSafelyConvertible<Real,Float64> : True { };
template<> struct IsSafelyConvertible<ExactFloat64,UpperInterval> : True { };
template<> struct IsSafelyConvertible<UpperInterval,UpperInterval> : True { };
template<> struct IsSafelyConvertible<UpperInterval,Float64> : True { };
template<> struct IsSafelyConvertible<Float64,Float64> : True { };

template<> struct IsSafelyConvertible<Rational,Rational> : True { };
template<> struct IsSafelyConvertible<Rational,Real> : True { };
template<> struct IsSafelyConvertible<Rational,UpperInterval> : True { };
template<> struct IsSafelyConvertible<Rational,Float64> : True { };
template<> struct IsSafelyConvertible<Nat,Rational> : True { };
template<> struct IsSafelyConvertible<Int,Rational> : True { };

template<> struct IsSafelyConvertible<Nat,Real> : True { };
template<> struct IsSafelyConvertible<Nat,ExactFloat64> : True { };
template<> struct IsSafelyConvertible<Nat,UpperInterval> : True { };
template<> struct IsSafelyConvertible<Nat,Float64> : True { };
template<> struct IsSafelyConvertible<Int,Real> : True { };
template<> struct IsSafelyConvertible<Int,ExactFloat64> : True { };
template<> struct IsSafelyConvertible<Int,UpperInterval> : True { };
template<> struct IsSafelyConvertible<Int,Float64> : True { };
template<> struct IsSafelyConvertible<double,ExactFloat64> : True { };
template<> struct IsSafelyConvertible<double,UpperInterval> : True { };
template<> struct IsSafelyConvertible<double,Float64> : True { };

template<class FROM, class TO> struct IsNumericCastable : IsSafelyConvertible<FROM,TO> { };
template<> struct IsNumericCastable<Float64,ExactFloat64> : True { };
template<> struct IsNumericCastable<Float64,UpperInterval> : True { };

// Type deduction for numerical arithmetic
template<class X1, class X2> struct Arithmetic { };
template<> struct Arithmetic<Real,Real> { typedef Real ResultType; };
template<> struct Arithmetic<ExactFloat64,ExactFloat64> { typedef ValidatedFloat64 ResultType; };
template<> struct Arithmetic<ExactFloat64,ValidatedFloat64> { typedef ValidatedFloat64 ResultType; };
template<> struct Arithmetic<ExactFloat64,ApproximateFloat64> { typedef ApproximateFloat64 ResultType; };
template<> struct Arithmetic<ValidatedFloat64,ExactFloat64> { typedef ValidatedFloat64 ResultType; };
template<> struct Arithmetic<ValidatedFloat64,ValidatedFloat64> { typedef ValidatedFloat64 ResultType; };
template<> struct Arithmetic<ValidatedFloat64,ApproximateFloat64> { typedef ApproximateFloat64 ResultType; };
template<> struct Arithmetic<ApproximateFloat64,ExactFloat64> { typedef ApproximateFloat64 ResultType; };
template<> struct Arithmetic<ApproximateFloat64,ValidatedFloat64> { typedef ApproximateFloat64 ResultType; };
template<> struct Arithmetic<ApproximateFloat64,ApproximateFloat64> { typedef ApproximateFloat64 ResultType; };

template<> struct Arithmetic<UpperInterval,ExactFloat64> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<ExactFloat64,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<UpperInterval,Real> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Real,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<UpperInterval,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Float64,ExactFloat64> { typedef Float64 ResultType; };
template<> struct Arithmetic<ExactFloat64,Float64> { typedef Float64 ResultType; };
template<> struct Arithmetic<Float64,UpperInterval> { typedef Float64 ResultType; };
template<> struct Arithmetic<UpperInterval,Float64> { typedef Float64 ResultType; };
template<> struct Arithmetic<Float64,Float64> { typedef Float64 ResultType; };
template<> struct Arithmetic<Float64,Real> { typedef Float64 ResultType; };
template<> struct Arithmetic<Real,Float64> { typedef Float64 ResultType; };
template<> struct Arithmetic<double,Real> { typedef Real ResultType; };
template<> struct Arithmetic<double,Float64> { typedef Float64 ResultType; };
template<> struct Arithmetic<double,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Real,double> { typedef Real ResultType; };
template<> struct Arithmetic<Float64,double> { typedef Float64 ResultType; };
template<> struct Arithmetic<UpperInterval,double> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Rational,Rational> { typedef Rational ResultType; };
template<> struct Arithmetic<Rational,double> { typedef Rational ResultType; };
template<> struct Arithmetic<double,Rational> { typedef Rational ResultType; };
template<> struct Arithmetic<Rational,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<UpperInterval,Rational> { typedef UpperInterval ResultType; };


} // namespace Ariadne

#endif
