/***************************************************************************
 *            numeric/float.decl.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/float.decl.h
 *  \brief
 */

#include "utility/typedefs.h"
#include "numeric/paradigm.h"

#ifndef ARIADNE_FLOAT_DECL_H
#define ARIADNE_FLOAT_DECL_H

namespace Ariadne {

class Real;
class TwoExp;

class Precision64;
class PrecisionMP;

template<class F> using PrecisionType = typename F::PrecisionType;

struct DecimalPlaces { int _places; DecimalPlaces(int plc) : _places(plc) { } operator int() const { return _places; } };
struct DecimalPrecision { uint _figures; operator uint() const { return _figures; } };

class Float64;
class FloatMP;

using RawFloat64 = Float64;
using RawFloatMP = FloatMP;

RawFloat64 cast_raw_float(Precision64);
RawFloatMP cast_raw_float(PrecisionMP);
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));

template<class PR> class FloatApproximation;
template<class PR> class FloatLowerBound;
template<class PR> class FloatUpperBound;
template<class PR> class FloatBounds;
template<class PR> class FloatBall;
template<class PR> class FloatValue;
template<class PR> class FloatError;

template<class PR> class PositiveFloatApproximation;
template<class PR> class PositiveFloatLowerBound;
template<class PR> class PositiveFloatUpperBound;
template<class PR> class PositiveFloatBounds;
template<class PR> class PositiveFloatBall;
template<class PR> class PositiveFloatValue;

template<class P, class PR> class FloatTypedef;
template<class PR> struct FloatTypedef<ApproximateTag,PR> { typedef FloatApproximation<PR> Type; };
template<class PR> struct FloatTypedef<LowerTag,PR> { typedef FloatLowerBound<PR> Type; };
template<class PR> struct FloatTypedef<UpperTag,PR> { typedef FloatUpperBound<PR> Type; };
template<class PR> struct FloatTypedef<BoundedTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR> struct FloatTypedef<MetricTag,PR> { typedef FloatBall<PR> Type; };
template<class PR> struct FloatTypedef<ExactTag,PR> { typedef FloatValue<PR> Type; };
template<class PR> struct FloatTypedef<ErrorTag,PR> { typedef FloatError<PR> Type; };

template<class PR> struct FloatTypedef<ValidatedTag,PR> { typedef FloatBounds<PR> Type; };
template<class PR> struct FloatTypedef<EffectiveTag,PR> { typedef FloatBall<PR> Type; };

template<class P, class PR> using Float = typename FloatTypedef<P,PR>::Type;
//template<class P> using Float64=Float<P,Precision64>;
//template<class P> using FloatMP=Float<P,PrecisionMP>;

using Float64Approximation = Float<ApproximateTag,Precision64>;
using Float64LowerBound = Float<LowerTag,Precision64>;
using Float64UpperBound = Float<UpperTag,Precision64>;
using Float64Bounds = Float<BoundedTag,Precision64>;
using Float64Ball = Float<MetricTag,Precision64>;
using Float64Value = Float<ExactTag,Precision64>;
using Float64Error = FloatError<Precision64>;
using PositiveFloat64Approximation = PositiveFloatApproximation<Precision64>;
using PositiveFloat64LowerBound = PositiveFloatLowerBound<Precision64>;
using PositiveFloat64UpperBound = PositiveFloatUpperBound<Precision64>;
using PositiveFloat64Bounds = PositiveFloatBounds<Precision64>;
using PositiveFloat64Ball = PositiveFloatBall<Precision64>;
using PositiveFloat64Value = PositiveFloatValue<Precision64>;

using FloatMPApproximation = Float<ApproximateTag,PrecisionMP>;
using FloatMPLowerBound = Float<LowerTag,PrecisionMP>;
using FloatMPUpperBound = Float<UpperTag,PrecisionMP>;
using FloatMPBounds= Float<BoundedTag,PrecisionMP>;
using FloatMPBall = Float<MetricTag,PrecisionMP>;
using FloatMPValue = Float<ExactTag,PrecisionMP>;
using FloatMPError = FloatError<PrecisionMP>;
using PositiveFloatMPApproximation = PositiveFloatApproximation<PrecisionMP>;
using PositiveFloatMPLowerBound = PositiveFloatLowerBound<PrecisionMP>;
using PositiveFloatMPUpperBound = PositiveFloatUpperBound<PrecisionMP>;
using PositiveFloatMPBounds = PositiveFloatBounds<PrecisionMP>;
using PositiveFloatMPBall = PositiveFloatBall<PrecisionMP>;
using PositiveFloatMPValue = PositiveFloatValue<PrecisionMP>;

template<class F> using Approximation = FloatApproximation<typename F::PrecisionType>;
template<class F> using LowerBound = FloatLowerBound<typename F::PrecisionType>;
template<class F> using UpperBound = FloatUpperBound<typename F::PrecisionType>;
template<class F> using Bounds = FloatBounds<typename F::PrecisionType>;
template<class F> using Ball = FloatBall<typename F::PrecisionType>;
template<class F> using Value = FloatValue<typename F::PrecisionType>;
template<class F> using Error = FloatError<typename F::PrecisionType>;


template<class X> struct IsFloat : False { };
template<class PR> struct IsFloat<FloatApproximation<PR>> : True { };
template<class PR> struct IsFloat<FloatLowerBound<PR>> : True { };
template<class PR> struct IsFloat<FloatUpperBound<PR>> : True { };
template<class PR> struct IsFloat<FloatBounds<PR>> : True { };
template<class PR> struct IsFloat<FloatBall<PR>> : True { };
template<class PR> struct IsFloat<FloatValue<PR>> : True { };
template<class PR> struct IsFloat<FloatError<PR>> : True { };
template<> struct IsFloat<Dbl> : False { };


template<class T> struct IsNumericType;
template<> struct IsNumericType<Dbl> : True { };
template<class PR> struct IsNumericType<FloatApproximation<PR>> : True { };
template<class PR> struct IsNumericType<FloatBounds<PR>> : True { };
template<class PR> struct IsNumericType<FloatBall<PR>> : True { };
template<class PR> struct IsNumericType<FloatValue<PR>> : True { };

} // namespace Ariadne

#endif
