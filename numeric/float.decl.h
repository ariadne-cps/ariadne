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
template<class PR> using ApproximateFloat=Float<Approximate,PR>;
template<class PR> using LowerFloat=Float<Lower,PR>;
template<class PR> using UpperFloat=Float<Upper,PR>;
template<class PR> using BoundedFloat=Float<Bounded,PR>;
template<class PR> using MetricFloat=Float<Metric,PR>;
template<class PR> using ExactFloat=Float<Exact,PR>;
template<class PR> using ErrorFloat=Float<Error,PR>;

using ApproximateFloat64 = Float<Approximate,Precision64>;
using LowerFloat64 = Float<Lower,Precision64>;
using UpperFloat64 = Float<Upper,Precision64>;
using BoundedFloat64 = Float<Bounded,Precision64>;
using MetricFloat64 = Float<Metric,Precision64>;
using ExactFloat64 = Float<Exact,Precision64>;
using PositiveApproximateFloat64 = Float<PositiveApproximate,Precision64>;
using PositiveLowerFloat64 = Float<PositiveLower,Precision64>;
using PositiveUpperFloat64 = Float<PositiveUpper,Precision64>;
using PositiveBoundedFloat64 = Float<PositiveBounded,Precision64>;
using PositiveMetricFloat64 = Float<PositiveMetric,Precision64>;
using PositiveExactFloat64 = Float<PositiveExact,Precision64>;
using ErrorFloat64 = Float<Error,Precision64>;

using ApproximateFloatMP = Float<Approximate,PrecisionMP>;
using LowerFloatMP = Float<Lower,PrecisionMP>;
using UpperFloatMP = Float<Upper,PrecisionMP>;
using BoundedFloatMP= Float<Bounded,PrecisionMP>;
using MetricFloatMP = Float<Metric,PrecisionMP>;
using ExactFloatMP = Float<Exact,PrecisionMP>;
using PositiveApproximateFloatMP = Float<PositiveApproximate,PrecisionMP>;
using PositiveLowerFloatMP = Float<PositiveLower,PrecisionMP>;
using PositiveUpperFloatMP = Float<PositiveUpper,PrecisionMP>;
using PositiveBoundedFloatMP = Float<PositiveBounded,PrecisionMP>;
using PositiveMetricFloatMP = Float<PositiveMetric,PrecisionMP>;
using PositiveExactFloatMP = Float<PositiveExact,PrecisionMP>;
using ErrorFloatMP = Float<Error,PrecisionMP>;


template<class X> struct IsFloat : False { };
template<class P, class PR> struct IsFloat<Float<P,PR>> : True { };
template<> struct IsFloat<Dbl> : False { };


template<class T> struct IsNumericType;
template<> struct IsNumericType<Dbl> : True { };
template<class P, class PR> struct IsNumericType<Float<P,PR>> : True { };

} // namespace Ariadne

#endif
