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
template<class PR> using ApproximateFloat=Float<ApproximateTag,PR>;
template<class PR> using LowerFloat=Float<LowerTag,PR>;
template<class PR> using UpperFloat=Float<UpperTag,PR>;
template<class PR> using BoundedFloat=Float<BoundedTag,PR>;
template<class PR> using MetricFloat=Float<MetricTag,PR>;
template<class PR> using ExactFloat=Float<ExactTag,PR>;
template<class PR> using ErrorFloat=Float<ErrorTag,PR>;

using ApproximateFloat64 = Float<ApproximateTag,Precision64>;
using LowerFloat64 = Float<LowerTag,Precision64>;
using UpperFloat64 = Float<UpperTag,Precision64>;
using BoundedFloat64 = Float<BoundedTag,Precision64>;
using MetricFloat64 = Float<MetricTag,Precision64>;
using ExactFloat64 = Float<ExactTag,Precision64>;
using PositiveApproximateFloat64 = Float<PositiveApproximateTag,Precision64>;
using PositiveLowerFloat64 = Float<PositiveLowerTag,Precision64>;
using PositiveUpperFloat64 = Float<PositiveUpperTag,Precision64>;
using PositiveBoundedFloat64 = Float<PositiveBoundedTag,Precision64>;
using PositiveMetricFloat64 = Float<PositiveMetricTag,Precision64>;
using PositiveExactFloat64 = Float<PositiveExactTag,Precision64>;
using ErrorFloat64 = Float<ErrorTag,Precision64>;

using ApproximateFloatMP = Float<ApproximateTag,PrecisionMP>;
using LowerFloatMP = Float<LowerTag,PrecisionMP>;
using UpperFloatMP = Float<UpperTag,PrecisionMP>;
using BoundedFloatMP= Float<BoundedTag,PrecisionMP>;
using MetricFloatMP = Float<MetricTag,PrecisionMP>;
using ExactFloatMP = Float<ExactTag,PrecisionMP>;
using PositiveApproximateFloatMP = Float<PositiveApproximateTag,PrecisionMP>;
using PositiveLowerFloatMP = Float<PositiveLowerTag,PrecisionMP>;
using PositiveUpperFloatMP = Float<PositiveUpperTag,PrecisionMP>;
using PositiveBoundedFloatMP = Float<PositiveBoundedTag,PrecisionMP>;
using PositiveMetricFloatMP = Float<PositiveMetricTag,PrecisionMP>;
using PositiveExactFloatMP = Float<PositiveExactTag,PrecisionMP>;
using ErrorFloatMP = Float<ErrorTag,PrecisionMP>;


template<class X> struct IsFloat : False { };
template<class P, class PR> struct IsFloat<Float<P,PR>> : True { };
template<> struct IsFloat<Dbl> : False { };


template<class T> struct IsNumericType;
template<> struct IsNumericType<Dbl> : True { };
template<class P, class PR> struct IsNumericType<Float<P,PR>> : True { };

} // namespace Ariadne

#endif
