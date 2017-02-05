/***************************************************************************
 *            numeric/number.decl.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/number.decl.hpp
 *  \brief
 */

#ifndef ARIADNE_NUMBER_DECL_HPP
#define ARIADNE_NUMBER_DECL_HPP

#include <stdexcept>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"
#include "numeric/paradigm.hpp"

namespace Ariadne {

/************ Number *********************************************************/

class DivideByZeroError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};


template<class X> struct IsNumericType : False { };
template<class X, class T> using EnableIfNumericType = EnableIf<IsNumericType<X>,T>;

template<class X> using NumericType = typename X::NumericType;

template<class X> struct NumericTraits;
template<class X> using GenericTrait = typename NumericTraits<X>::GenericType;
template<class X> using PositiveTrait = typename NumericTraits<X>::PositiveType;
template<class X> using OppositeTrait = typename NumericTraits<X>::OppositeType;
template<class X> using LessTrait = typename NumericTraits<X>::LessType;
template<class X> using EqualsTrait = typename NumericTraits<X>::EqualsType;

template<class X> using PropertiesType = typename X::PropertiesType;
template<class X> using GenericType = typename X::GenericType;


typedef uint Nat;
typedef int Int;
typedef double Dbl;

class Nat32;
class Nat64;
class Int32;
class Int64;

class ExactDouble;

class Dyadic;
class Decimal;

class Integer;
class Rational;
class Real;

template<class R, class A> R integer_cast(const A& _a);

template<class X> class Positive;

template<class P=Void> class Number;

template<> struct IsNumericType<Nat>;
template<> struct IsNumericType<Int>;
template<> struct IsNumericType<Dbl>;

template<> struct IsNumericType<Integer>;
template<> struct IsNumericType<Rational>;
template<> struct IsNumericType<Real>;

template<class X> struct IsGenericNumericType : IsConvertible<X,Rational> { };
template<> struct IsGenericNumericType<ExactDouble> : True { };
template<> struct IsGenericNumericType<Integer> : True { };
template<> struct IsGenericNumericType<Dyadic> : True { };
template<> struct IsGenericNumericType<Rational> : True { };
template<> struct IsGenericNumericType<Real> : True { };
template<> struct IsGenericNumericType<Dbl> : True { };
template<class P> struct IsGenericNumericType<Number<P>> : True { };


using ExactNumber=Number<ExactTag>;
using EffectiveNumber=Number<EffectiveTag>;
using EffectiveUpperNumber=Number<EffectiveUpperTag>;
using EffectiveLowerNumber=Number<EffectiveLowerTag>;
using ValidatedNumber=Number<ValidatedTag>;
using ValidatedUpperNumber=Number<ValidatedUpperTag>;
using ValidatedLowerNumber=Number<ValidatedLowerTag>;
using ApproximateNumber=Number<ApproximateTag>;

} // namespace Ariadne

#include "float.decl.hpp"

namespace Ariadne {

using ExactNumericType=Float64Value;
using EffectiveNumericType=Real;
using ValidatedNumericType=Float64Bounds;
using UpperNumericType=Float64UpperBound;
using LowerNumericType=Float64LowerBound;
using ApproximateNumericType=Float64Approximation;
using ErrorNumericType=Float64Error;

} // namespace Ariadne

#endif
