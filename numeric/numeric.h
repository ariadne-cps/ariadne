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

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "numeric/float.h"
#include "numeric/float-exact.h"
#include "numeric/float-validated.h"
#include "numeric/float-approximate.h"
#include "geometry/interval.h"
#include "numeric/real.h"

namespace Ariadne {


//! \ingroup NumericModule
//! \brief Cast one %Ariadne numerical type or builtin numerical type to another.
template<class R, class A> inline R numeric_cast(const A& a) { return R(a); }
template<> inline int numeric_cast(const Float& a) { return int(a.get_d()); }
template<> inline double numeric_cast(const Float& a) { return a.get_d(); }
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline double numeric_cast(const ExactFloat& a) { return a.get_d(); }
template<> inline double numeric_cast(const ValidatedFloat& a) { return a.get_d(); }
template<> inline double numeric_cast(const ApproximateFloat& a) { return a.get_d(); }
template<> inline float numeric_cast(const double& a) { return a; }
template<> inline float numeric_cast(const Float& a) { return a.get_d(); }
template<> inline float numeric_cast(const Real& a) { return a.get_d(); }

template<> inline Real numeric_cast(const Float& a) { return Real(ExactFloat(a)); }
template<> inline Real numeric_cast(const ExactFloat& a) { return Real(a); }
template<> inline Real numeric_cast(const ValidatedFloat& a) { return Real(make_exact(ApproximateFloat(a))); }

template<> inline UpperInterval numeric_cast(const Float& a) { return UpperInterval(a); }
template<> inline Float numeric_cast(const UpperInterval& a) { return midpoint(a).raw(); }

template<class R, class A> inline R approx_cast(const A& a);
template<> inline double approx_cast(const Float& a) { return a.get_d(); }
template<> inline double approx_cast(const ExactFloat& a) { return a.get_d(); }

//! \ingroup NumericModule
//! \related Float
//! \brief Converts \a e to an object of type \a X, which may either be an
//! \c Float or \c ExactInterval, with the semantics that \a e denotes and error bound.
//! Returns the float 0.0 (since floating-point computations do not keep track of errors)
//! and the interval [-e,+e].
template<class X> inline X convert_error(const Float& e);
template<> inline Float convert_error<Float>(const Float& e) { return 0.0; }

template<class X> inline X convert_error(const ApproximateFloat& e);
template<> inline ApproximateFloat convert_error<ApproximateFloat>(const ApproximateFloat& e) { return 0.0; }
template<> inline ValidatedFloat convert_error<ValidatedFloat>(const ApproximateFloat& e) { return ValidatedFloat(-e.value(),+e.value()); }

// Use 'enable_if' style template to restrict allowable instances. See the Boost documentation
// for enable_if to see how this works.
template<class X> struct IsNumeric { static const bool value = false; };
template<> struct IsNumeric<unsigned int> { static const bool value = true; };
template<> struct IsNumeric<int> { static const bool value = true; };
template<> struct IsNumeric<double> { static const bool value = true; };
template<> struct IsNumeric<Float> { static const bool value = true; };
template<> struct IsNumeric<ExactFloat> { static const bool value = true; };
template<> struct IsNumeric<ValidatedFloat> { static const bool value = true; };
template<> struct IsNumeric<ApproximateFloat> { static const bool value = true; };
template<> struct IsNumeric<Real> { static const bool value = true; };

#ifdef HAVE_GMPXX_H
template<> struct IsNumeric<Integer> { static const bool value = true; };
template<> struct IsNumeric<Rational> { static const bool value = true; };
#endif // HAVE_GMPXX_H

template<class X, class T> struct EnableIfNumeric : EnableIf<IsNumeric<X>,T> { };

template<class F, class T> struct IsSafelyConvertible : False { };
template<> struct IsSafelyConvertible<ExactFloat,ExactFloat> : True { };
template<> struct IsSafelyConvertible<ExactFloat,Real> : True { };
template<> struct IsSafelyConvertible<ExactFloat,ValidatedFloat> : True { };
template<> struct IsSafelyConvertible<ExactFloat,ApproximateFloat> : True { };
template<> struct IsSafelyConvertible<Real,Real> : True { };
template<> struct IsSafelyConvertible<Real,ValidatedFloat> : True { };
template<> struct IsSafelyConvertible<Real,ApproximateFloat> : True { };
template<> struct IsSafelyConvertible<ValidatedFloat,ValidatedFloat> : True { };
template<> struct IsSafelyConvertible<ValidatedFloat,ApproximateFloat> : True { };
template<> struct IsSafelyConvertible<ApproximateFloat,ApproximateFloat> : True { };

template<> struct IsSafelyConvertible<Real,UpperInterval> : True { };
template<> struct IsSafelyConvertible<Real,Float> : True { };
template<> struct IsSafelyConvertible<ExactFloat,UpperInterval> : True { };
template<> struct IsSafelyConvertible<UpperInterval,UpperInterval> : True { };
template<> struct IsSafelyConvertible<UpperInterval,Float> : True { };
template<> struct IsSafelyConvertible<Float,Float> : True { };
#ifdef HAVE_GMPXX_H
template<> struct IsSafelyConvertible<Rational,Rational> : True { };
template<> struct IsSafelyConvertible<Rational,Real> : True { };
template<> struct IsSafelyConvertible<Rational,UpperInterval> : True { };
template<> struct IsSafelyConvertible<Rational,Float> : True { };
template<> struct IsSafelyConvertible<uint,Rational> : True { };
template<> struct IsSafelyConvertible<int,Rational> : True { };
#endif // HAVE_GMPXX_H
template<> struct IsSafelyConvertible<uint,Real> : True { };
template<> struct IsSafelyConvertible<uint,ExactFloat> : True { };
template<> struct IsSafelyConvertible<uint,UpperInterval> : True { };
template<> struct IsSafelyConvertible<uint,Float> : True { };
template<> struct IsSafelyConvertible<int,Real> : True { };
template<> struct IsSafelyConvertible<int,ExactFloat> : True { };
template<> struct IsSafelyConvertible<int,UpperInterval> : True { };
template<> struct IsSafelyConvertible<int,Float> : True { };
template<> struct IsSafelyConvertible<double,ExactFloat> : True { };
template<> struct IsSafelyConvertible<double,UpperInterval> : True { };
template<> struct IsSafelyConvertible<double,Float> : True { };

template<class FROM, class TO> struct IsNumericCastable : IsSafelyConvertible<FROM,TO> { };
template<> struct IsNumericCastable<Float,ExactFloat> : True { };
template<> struct IsNumericCastable<Float,UpperInterval> : True { };

// Type deduction for numerical arithmetic
template<class X1, class X2> struct Arithmetic { };
template<> struct Arithmetic<Real,Real> { typedef Real ResultType; };
template<> struct Arithmetic<ExactFloat,ExactFloat> { typedef ValidatedFloat ResultType; };
template<> struct Arithmetic<ExactFloat,ValidatedFloat> { typedef ValidatedFloat ResultType; };
template<> struct Arithmetic<ExactFloat,ApproximateFloat> { typedef ApproximateFloat ResultType; };
template<> struct Arithmetic<ValidatedFloat,ExactFloat> { typedef ValidatedFloat ResultType; };
template<> struct Arithmetic<ValidatedFloat,ValidatedFloat> { typedef ValidatedFloat ResultType; };
template<> struct Arithmetic<ValidatedFloat,ApproximateFloat> { typedef ApproximateFloat ResultType; };
template<> struct Arithmetic<ApproximateFloat,ExactFloat> { typedef ApproximateFloat ResultType; };
template<> struct Arithmetic<ApproximateFloat,ValidatedFloat> { typedef ApproximateFloat ResultType; };
template<> struct Arithmetic<ApproximateFloat,ApproximateFloat> { typedef ApproximateFloat ResultType; };

template<> struct Arithmetic<UpperInterval,ExactFloat> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<ExactFloat,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<UpperInterval,Real> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Real,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<UpperInterval,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Float,ExactFloat> { typedef Float ResultType; };
template<> struct Arithmetic<ExactFloat,Float> { typedef Float ResultType; };
template<> struct Arithmetic<Float,UpperInterval> { typedef Float ResultType; };
template<> struct Arithmetic<UpperInterval,Float> { typedef Float ResultType; };
template<> struct Arithmetic<Float,Float> { typedef Float ResultType; };
template<> struct Arithmetic<Float,Real> { typedef Float ResultType; };
template<> struct Arithmetic<Real,Float> { typedef Float ResultType; };
template<> struct Arithmetic<double,Real> { typedef Real ResultType; };
template<> struct Arithmetic<double,Float> { typedef Float ResultType; };
template<> struct Arithmetic<double,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<Real,double> { typedef Real ResultType; };
template<> struct Arithmetic<Float,double> { typedef Float ResultType; };
template<> struct Arithmetic<UpperInterval,double> { typedef UpperInterval ResultType; };
#ifdef HAVE_GMPXX_H
template<> struct Arithmetic<Rational,Rational> { typedef Rational ResultType; };
template<> struct Arithmetic<Rational,double> { typedef Rational ResultType; };
template<> struct Arithmetic<double,Rational> { typedef Rational ResultType; };
template<> struct Arithmetic<Rational,UpperInterval> { typedef UpperInterval ResultType; };
template<> struct Arithmetic<UpperInterval,Rational> { typedef UpperInterval ResultType; };
#endif // HAVE_GMPXX_H


} // namespace Ariadne

#endif
