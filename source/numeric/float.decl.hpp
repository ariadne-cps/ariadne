/***************************************************************************
 *            numeric/float.decl.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/float.decl.hpp
 *  \brief
 */

#include "../utility/typedefs.hpp"
#include "../numeric/paradigm.hpp"

#ifndef ARIADNE_FLOAT_DECL_HPP
#define ARIADNE_FLOAT_DECL_HPP

namespace Ariadne {

class Real;
class TwoExp;

class DoublePrecision;
class MultiplePrecision;

using DP = DoublePrecision;
using MP = MultiplePrecision;

template<class F> using PrecisionType = typename F::PrecisionType;

struct DecimalPlaces { Nat _places; explicit DecimalPlaces(Nat plc) : _places(plc) { } operator uint() const { return _places; } };
struct DecimalPrecision { Nat _figures; explicit DecimalPrecision(Nat figs) : _figures(figs) { } operator uint() const { return _figures; } };

template<class X> class Positive;

class FloatDP;
class FloatMP;

using DoublePrecisionFloat = FloatDP;
using MultiplePrecisionFloat = FloatMP;

using RawFloatDP = FloatDP;
using RawFloatMP = FloatMP;

template<class PR> struct RawFloatTypedef;
template<> struct RawFloatTypedef<DoublePrecision> { typedef FloatDP Type; };
template<> struct RawFloatTypedef<MultiplePrecision> { typedef FloatMP Type; };
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));
template<class PR> using RawFloatType = typename RawFloatTypedef<PR>::Type;

RawFloatDP cast_raw_float(DoublePrecision);
RawFloatMP cast_raw_float(MultiplePrecision);
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));


/*
template<class PR> class FloatApproximation;
template<class PR> class FloatLowerBound;
template<class PR> class FloatUpperBound;
template<class PR> class FloatBounds;
template<class PR, class PRE=PR> class FloatBall;
template<class PR> class FloatValue;
template<class PR> class FloatError;
*/

template<class F> class Approximation;
template<class F> class LowerBound;
template<class F> class UpperBound;
template<class F> class Bounds;
template<class F, class FE=F> class Ball;
template<class F> class Value;

template<class F> class Error;

//@{
//! \relates Positive
//! \name Type shorthands
template<class F> using PositiveApproximation = Positive<Approximation<F>>; //!< .
template<class F> using PositiveLowerBound = Positive<LowerBound<F>>; //!< .
template<class F> using PositiveUpperBound = Positive<UpperBound<F>>; //!< .
template<class F> using PositiveBounds = Positive<Bounds<F>>; //!< .
template<class F, class FE=F> using PositiveBall = Positive<Ball<F>>; //!< .
template<class F> using PositiveValue = Positive<Value<F>>; //!< .
//@}

//@{
//! \ingroup NumericModule
//! \name Type shorthands
template<class PR> using FloatApproximation = Approximation<RawFloatType<PR>>; //!< .
template<class PR> using FloatLowerBound = LowerBound<RawFloatType<PR>>; //!< .
template<class PR> using FloatUpperBound = UpperBound<RawFloatType<PR>>; //!< .
template<class PR> using FloatBounds = Bounds<RawFloatType<PR>>; //!< .
template<class PR, class PRE=PR> using FloatBall = Ball<RawFloatType<PR>,RawFloatType<PRE>>; //!< .
template<class PR> using FloatValue = Value<RawFloatType<PR>>; //!< .

template<class PR> using PositiveFloatApproximation = Positive<FloatApproximation<PR>>; //!< .
template<class PR> using PositiveFloatLowerBound = Positive<FloatLowerBound<PR>>; //!< .
template<class PR> using PositiveFloatUpperBound = Positive<FloatUpperBound<PR>>; //!< .
template<class PR> using PositiveFloatBounds = Positive<FloatBounds<PR>>; //!< .
template<class PR, class PRE=PR> using PositiveFloatBall = Positive<FloatBall<PR,PRE>>; //!< .
template<class PR> using PositiveFloatValue = Positive<FloatValue<PR>>; //!< .

template<class PR> using FloatError = Error<RawFloatType<PR>>;
//@}

template<class P, class PR, class PRE=PR> struct FloatTypedef;
template<class PR> struct FloatTypedef<ApproximateTag,PR> { typedef FloatApproximation<PR> Type; };
template<class PR> struct FloatTypedef<LowerTag,PR> { typedef FloatLowerBound<PR> Type; };
template<class PR> struct FloatTypedef<UpperTag,PR> { typedef FloatUpperBound<PR> Type; };
template<class PR> struct FloatTypedef<BoundedTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR, class PRE> struct FloatTypedef<MetricTag,PR,PRE> { typedef FloatBall<PR,PRE> Type; };
template<class PR> struct FloatTypedef<ExactTag,PR> { typedef FloatValue<PR> Type; };
template<class PR> struct FloatTypedef<ErrorTag,PR> { typedef FloatError<PR> Type; };

template<class PR> struct FloatTypedef<ValidatedTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR, class PRE> struct FloatTypedef<EffectiveTag,PR,PRE> { typedef FloatBall<PR,PRE> Type; };

template<class P, class PR, class PRE=PR> using FloatType = typename FloatTypedef<P,PR,PRE>::Type;


//template<class P, class PR, class PRE=PR> using Float = typename FloatTypedef<P,PR,PRE>::Type;
//template<class P> using FloatDP=Float<P,DoublePrecision>;
//template<class P> using FloatMP=Float<P,MultiplePrecision>;

//@{
//! \ingroup NumericModule
//! \name Type synonyms for FloatDP based numbers.
using FloatDPApproximation = FloatApproximation<DoublePrecision>; //!< An approximation to a number, represented in double precision.
using FloatDPLowerBound = FloatLowerBound<DoublePrecision>; //!< A lower bound for a number, represented in double precision.
using FloatDPUpperBound = FloatUpperBound<DoublePrecision>; //!< An upper bound for a number, represented in double precision.
using FloatDPBounds = FloatBounds<DoublePrecision>; //!< Lower and upper bounds for a number, represented in double precision.
using FloatDPBall = FloatBall<DoublePrecision>; //!< A ball around a number, with the approximating value and error bound represented in double precision.
using FloatDPValue = FloatValue<DoublePrecision>; //!< A double-precision floating-point object representing a number exactly.
using FloatDPError = FloatError<DoublePrecision>; //!< An over-approximation to the error of a computation, represented in double precision.
using PositiveFloatDPApproximation = PositiveFloatApproximation<DoublePrecision>; //!< .
using PositiveFloatDPLowerBound = PositiveFloatLowerBound<DoublePrecision>; //!< .
using PositiveFloatDPUpperBound = PositiveFloatUpperBound<DoublePrecision>; //!< .
using PositiveFloatDPBounds = PositiveFloatBounds<DoublePrecision>; //!< .
using PositiveFloatDPBall = PositiveFloatBall<DoublePrecision>; //!< .
using PositiveFloatDPValue = PositiveFloatValue<DoublePrecision>; //!< .

using FloatMPDPBall = FloatBall<MultiplePrecision,DoublePrecision>; //!< A ball around a number, with the approximating value represented in multiple precision, and the error in double precision.
//@}

//@{
//! \ingroup NumericModule
//! \ingroup NumericModule
//! \name Type synonyms for FloatMP based numbers.
using FloatMPApproximation = FloatApproximation<MultiplePrecision>; //!< .
using FloatMPLowerBound = FloatLowerBound<MultiplePrecision>; //!< .
using FloatMPUpperBound = FloatUpperBound<MultiplePrecision>; //!< .
using FloatMPBounds = FloatBounds<MultiplePrecision>; //!< .
using FloatMPBall = FloatBall<MultiplePrecision>; //!< .
using FloatMPValue = FloatValue<MultiplePrecision>; //!< .
using FloatMPError = FloatError<MultiplePrecision>; //!< .
using PositiveFloatMPApproximation = PositiveFloatApproximation<MultiplePrecision>; //!< .
using PositiveFloatMPLowerBound = PositiveFloatLowerBound<MultiplePrecision>; //!< .
using PositiveFloatMPUpperBound = PositiveFloatUpperBound<MultiplePrecision>; //!< .
using PositiveFloatMPBounds = PositiveFloatBounds<MultiplePrecision>; //!< .
using PositiveFloatMPBall = PositiveFloatBall<MultiplePrecision>; //!< .
using PositiveFloatMPValue = PositiveFloatValue<MultiplePrecision>; //!< .

using FloatMDPBall = FloatBall<MultiplePrecision,DoublePrecision>; //!< .
//@}


template<class X> struct IsFloat : False { };
template<> struct IsFloat<FloatDP> : True { };
template<> struct IsFloat<FloatMP> : True { };
template<class F> struct IsFloat<Approximation<F>> : IsFloat<F> { };
template<class F> struct IsFloat<LowerBound<F>> : IsFloat<F> { };
template<class F> struct IsFloat<UpperBound<F>> : IsFloat<F> { };
template<class F> struct IsFloat<Bounds<F>> : IsFloat<F> { };
template<class F, class FE> struct IsFloat<Ball<F,FE>> : IsFloat<F> { };
template<class F> struct IsFloat<Value<F>> : IsFloat<F> { };
template<class F> struct IsFloat<Error<F>> : IsFloat<F> { };
template<> struct IsFloat<Dbl> : False { };


template<class T> struct IsNumericType;
template<> struct IsNumericType<Dbl> : True { };
template<> struct IsNumericType<FloatDP> : True { };
template<> struct IsNumericType<FloatMP> : True { };
template<class F> struct IsNumericType<Approximation<F>> : IsNumericType<F> { };
template<class F> struct IsNumericType<Bounds<F>> : IsNumericType<F> { };
template<class F, class FE> struct IsNumericType<Ball<F,FE>> : IsNumericType<F> { };
template<class F> struct IsNumericType<Value<F>> : IsNumericType<F> { };

template<class T> struct NumericTraits;

template<class T> struct NumericTraits<Positive<T>> : public NumericTraits<T> {
    typedef Positive<typename NumericTraits<T>::GenericType> GenericType;
    typedef Positive<typename NumericTraits<T>::OppositeType> OppositeType;
};

} // namespace Ariadne

#endif
