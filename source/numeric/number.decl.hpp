/***************************************************************************
 *            numeric/number.decl.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/number.decl.hpp
 *  \brief
 */

#ifndef ARIADNE_NUMBER_DECL_HPP
#define ARIADNE_NUMBER_DECL_HPP

#include <stdexcept>

#include "../utility/metaprogramming.hpp"
#include "../utility/typedefs.hpp"
#include "../numeric/paradigm.hpp"

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

template<class X> struct IsConcrete : Has<GenericType,X> { };

template<class X> struct IsScalar;
template<class X> struct IsConcreteScalar : And<Has<GenericType,X>,IsScalar<X>> { };
template<class X> struct IsGenericScalar : And<Not<Has<GenericType,X>>,IsScalar<X>> { };


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

template<class F> class Bounds;
template<> class Bounds<Dyadic>;
using DyadicBounds = Bounds<Dyadic>;

template<class R, class A> R integer_cast(const A& _a);

//! \ingroup NumericModule
//! \brief A modifier declaring that a number is positive.
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


using ExactNumber=Number<ExactTag>; //!< Alias for generic exact numbers. \ingroup NumericModule
using EffectiveNumber=Number<EffectiveTag>; //!< Alias for generic effective numbers. \ingroup NumericModule
using EffectiveUpperNumber=Number<EffectiveUpperTag>; //!< Alias for generic effective upper numbers. \ingroup NumericModule
using EffectiveLowerNumber=Number<EffectiveLowerTag>; //!< Alias for generic effective lower numbers. \ingroup NumericModule
using ValidatedNumber=Number<ValidatedTag>; //!< Alias for generic validated numbers. \ingroup NumericModule
using ValidatedUpperNumber=Number<ValidatedUpperTag>; //!< Alias for generic validated upper numbers. \ingroup NumericModule
using ValidatedLowerNumber=Number<ValidatedLowerTag>; //!< Alias for generic validated lower numbers. \ingroup NumericModule
using ApproximateNumber=Number<ApproximateTag>; //!< Alias for generic approximate numbers. \ingroup NumericModule

template<class P> using PositiveNumber = Positive<Number<P>>;

using PositiveExactNumber=PositiveNumber<ExactTag>; //!< Alias for generic positive validated upper numbers. \ingroup NumericModule
using PositiveValidatedUpperNumber=PositiveNumber<ValidatedUpperTag>; //!< Alias for generic positive validated upper numbers. \ingroup NumericModule
using PositiveValidatedLowerNumber=PositiveNumber<ValidatedLowerTag>; //!< Alias for generic positive validated lower numbers. \ingroup NumericModule

using ValidatedErrorNumber = PositiveValidatedUpperNumber;

} // namespace Ariadne

#include "float.decl.hpp"

namespace Ariadne {

using ExactNumericType=FloatDPValue;
using EffectiveNumericType=Real;
using ValidatedNumericType=FloatDPBounds;
using UpperNumericType=FloatDPUpperBound;
using LowerNumericType=FloatDPLowerBound;
using ApproximateNumericType=FloatDPApproximation;
using ErrorNumericType=FloatDPError;

} // namespace Ariadne

#endif
