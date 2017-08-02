/***************************************************************************
 *            numeric/float.decl.hpp
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

/*! \file numeric/float.decl.hpp
 *  \brief
 */

#include "utility/typedefs.hpp"
#include "numeric/paradigm.hpp"

#ifndef ARIADNE_FLOAT_DECL_HPP
#define ARIADNE_FLOAT_DECL_HPP

namespace Ariadne {

class Real;
class TwoExp;

class DoublePrecision;
class MultiplePrecision;

using DoublePrecision = DoublePrecision;
using MultiplePrecision = MultiplePrecision;
using DP = DoublePrecision;
using MP = MultiplePrecision;

template<class F> using PrecisionType = typename F::PrecisionType;

struct DecimalPlaces { int _places; DecimalPlaces(int plc) : _places(plc) { } operator int() const { return _places; } };
struct DecimalPrecision { uint _figures; operator uint() const { return _figures; } };

template<class X> class Positive;

class FloatDP;
class FloatMP;

using RawFloatDP = FloatDP;
using RawFloatMP = FloatMP;

RawFloatDP cast_raw_float(DoublePrecision);
RawFloatMP cast_raw_float(MultiplePrecision);
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));

template<class PR> class FloatApproximation;
template<class PR> class FloatLowerBound;
template<class PR> class FloatUpperBound;
template<class PR> class FloatBounds;
template<class PR, class PRE=PR> class FloatBall;
template<class PR> class FloatValue;
template<class PR> class FloatError;

template<class PR> using PositiveFloatApproximation = Positive<FloatApproximation<PR>>;
template<class PR> using PositiveFloatLowerBound = Positive<FloatLowerBound<PR>>;
template<class PR> using PositiveFloatUpperBound = Positive<FloatUpperBound<PR>>;
template<class PR> using PositiveFloatBounds = Positive<FloatBounds<PR>>;
template<class PR, class PRE=PR> using PositiveFloatBall = Positive<FloatBall<PR,PRE>>;
template<class PR> using PositiveFloatValue = Positive<FloatValue<PR>>;

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

template<class P, class PR, class PRE=PR> using Float = typename FloatTypedef<P,PR,PRE>::Type;
//template<class P> using FloatDP=Float<P,DoublePrecision>;
//template<class P> using FloatMP=Float<P,MultiplePrecision>;

using FloatDPApproximation = Float<ApproximateTag,DoublePrecision>;
using FloatDPLowerBound = Float<LowerTag,DoublePrecision>;
using FloatDPUpperBound = Float<UpperTag,DoublePrecision>;
using FloatDPBounds = Float<BoundedTag,DoublePrecision>;
using FloatDPBall = Float<MetricTag,DoublePrecision>;
using FloatDPValue = Float<ExactTag,DoublePrecision>;
using FloatDPError = FloatError<DoublePrecision>;
using PositiveFloatDPApproximation = PositiveFloatApproximation<DoublePrecision>;
using PositiveFloatDPLowerBound = PositiveFloatLowerBound<DoublePrecision>;
using PositiveFloatDPUpperBound = PositiveFloatUpperBound<DoublePrecision>;
using PositiveFloatDPBounds = PositiveFloatBounds<DoublePrecision>;
using PositiveFloatDPBall = PositiveFloatBall<DoublePrecision>;
using PositiveFloatDPValue = PositiveFloatValue<DoublePrecision>;

using FloatMPApproximation = Float<ApproximateTag,MultiplePrecision>;
using FloatMPLowerBound = Float<LowerTag,MultiplePrecision>;
using FloatMPUpperBound = Float<UpperTag,MultiplePrecision>;
using FloatMPBounds= Float<BoundedTag,MultiplePrecision>;
using FloatMPBall = Float<MetricTag,MultiplePrecision>;
using FloatMPValue = Float<ExactTag,MultiplePrecision>;
using FloatMPError = FloatError<MultiplePrecision>;
using PositiveFloatMPApproximation = PositiveFloatApproximation<MultiplePrecision>;
using PositiveFloatMPLowerBound = PositiveFloatLowerBound<MultiplePrecision>;
using PositiveFloatMPUpperBound = PositiveFloatUpperBound<MultiplePrecision>;
using PositiveFloatMPBounds = PositiveFloatBounds<MultiplePrecision>;
using PositiveFloatMPBall = PositiveFloatBall<MultiplePrecision>;
using PositiveFloatMPValue = PositiveFloatValue<MultiplePrecision>;

using FloatMDPBall = Float<MetricTag,MultiplePrecision,DoublePrecision>;

template<class F> using Approximation = FloatApproximation<typename F::PrecisionType>;
template<class F> using LowerBound = FloatLowerBound<typename F::PrecisionType>;
template<class F> using UpperBound = FloatUpperBound<typename F::PrecisionType>;
template<class F> using Bounds = FloatBounds<typename F::PrecisionType>;
template<class F, class FE=F> using Ball = FloatBall<typename F::PrecisionType, typename FE::PrecisionType>;
template<class F> using Value = FloatValue<typename F::PrecisionType>;
template<class F> using Error = FloatError<typename F::PrecisionType>;


template<class X> struct IsFloat : False { };
template<class PR> struct IsFloat<FloatApproximation<PR>> : True { };
template<class PR> struct IsFloat<FloatLowerBound<PR>> : True { };
template<class PR> struct IsFloat<FloatUpperBound<PR>> : True { };
template<class PR> struct IsFloat<FloatBounds<PR>> : True { };
template<class PR, class PRE> struct IsFloat<FloatBall<PR,PRE>> : True { };
template<class PR> struct IsFloat<FloatValue<PR>> : True { };
template<class PR> struct IsFloat<FloatError<PR>> : True { };
template<> struct IsFloat<Dbl> : False { };


template<class T> struct IsNumericType;
template<> struct IsNumericType<Dbl> : True { };
template<class PR> struct IsNumericType<FloatApproximation<PR>> : True { };
template<class PR> struct IsNumericType<FloatBounds<PR>> : True { };
template<class PR, class PRE> struct IsNumericType<FloatBall<PR,PRE>> : True { };
template<class PR> struct IsNumericType<FloatValue<PR>> : True { };

} // namespace Ariadne

#endif
