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

struct DecimalPlaces { int _places; DecimalPlaces(int plc) : _places(plc) { } operator int() const { return _places; } };
struct DecimalPrecision { uint _figures; operator uint() const { return _figures; } };

class Float64;
class FloatMP;

using RawFloat64 = Float64;
using RawFloatMP = FloatMP;

//template<class PR> struct FloatTypedef;
//template<> struct FloatTypedef<Precision64> { typedef Float64 Type; };
//template<> struct FloatTypedef<PrecisionMP> { typedef FloatMP Type; };
//template<class PR> using RawFloat = typename FloatTypedef<PR>::Type;

RawFloat64 cast_raw_float(Precision64);
RawFloatMP cast_raw_float(PrecisionMP);
template<class PR> using RawFloat = decltype(cast_raw_float(declval<PR>()));


template<class P, class PR> class Float;
//template<class P> using Float64=Float<P,Precision64>;
//template<class P> using FloatMP=Float<P,PrecisionMP>;
template<class PR> using FloatApproximation=Float<ApproximateTag,PR>;
template<class PR> using FloatLowerBound=Float<LowerTag,PR>;
template<class PR> using FloatUpperBound=Float<UpperTag,PR>;
template<class PR> using FloatBounds=Float<BoundedTag,PR>;
template<class PR> using FloatBall=Float<MetricTag,PR>;
template<class PR> using FloatValue=Float<ExactTag,PR>;
template<class PR> using FloatError=Float<ErrorTag,PR>;

using Float64Approximation = Float<ApproximateTag,Precision64>;
using Float64LowerBound = Float<LowerTag,Precision64>;
using Float64UpperBound = Float<UpperTag,Precision64>;
using Float64Bounds = Float<BoundedTag,Precision64>;
using Float64Ball = Float<MetricTag,Precision64>;
using Float64Value = Float<ExactTag,Precision64>;
using PositiveFloat64Approximation = Float<PositiveApproximateTag,Precision64>;
using PositiveFloat64LowerBound = Float<PositiveLowerTag,Precision64>;
using PositiveFloat64UpperBound = Float<PositiveUpperTag,Precision64>;
using PositiveFloat64Bounds = Float<PositiveBoundedTag,Precision64>;
using PositiveFloat64Ball = Float<PositiveMetricTag,Precision64>;
using PositiveFloat64Value = Float<PositiveExactTag,Precision64>;
using Float64Error = Float<ErrorTag,Precision64>;

using FloatMPApproximation = Float<ApproximateTag,PrecisionMP>;
using FloatMPLowerBound = Float<LowerTag,PrecisionMP>;
using FloatMPUpperBound = Float<UpperTag,PrecisionMP>;
using FloatMPBounds= Float<BoundedTag,PrecisionMP>;
using FloatMPBall = Float<MetricTag,PrecisionMP>;
using FloatMPValue = Float<ExactTag,PrecisionMP>;
using PositiveFloatMPApproximation = Float<PositiveApproximateTag,PrecisionMP>;
using PositiveFloatMPLowerBound = Float<PositiveLowerTag,PrecisionMP>;
using PositiveFloatMPUpperBound = Float<PositiveUpperTag,PrecisionMP>;
using PositiveFloatMPBounds = Float<PositiveBoundedTag,PrecisionMP>;
using PositiveFloatMPBall = Float<PositiveMetricTag,PrecisionMP>;
using PositiveFloatMPValue = Float<PositiveExactTag,PrecisionMP>;
using FloatMPError = Float<ErrorTag,PrecisionMP>;


template<class X> struct IsFloat : False { };
template<class P, class PR> struct IsFloat<Float<P,PR>> : True { };
template<> struct IsFloat<Dbl> : False { };


template<class T> struct IsNumericType;
template<> struct IsNumericType<Dbl> : True { };
template<class P, class PR> struct IsNumericType<Float<P,PR>> : True { };

} // namespace Ariadne

#endif
