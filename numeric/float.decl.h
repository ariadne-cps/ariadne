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

RawFloat64 make_raw_float(Precision64);
RawFloatMP make_raw_float(PrecisionMP);
template<class PR> using RawFloat = decltype(make_raw_float(declval<PR>()));


template<class P, class PR> class Float;
template<class P, class PR> using FloatTemplate = Float<P,PR>;
template<class P> using Float64Template = Float<P,Precision64>;
template<class P> using FloatMPTemplate = Float<P,PrecisionMP>;

using ApproximateFloat64 = Float<Approximate,Precision64>;
using LowerFloat64 = Float<Lower,Precision64>;
using UpperFloat64 = Float<Upper,Precision64>;
using BoundedFloat64 = Float<Bounded,Precision64>;
using MetricFloat64 = Float<Metric,Precision64>;
using ExactFloat64 = Float<Exact,Precision64>;
using PositiveApproximateFloat64 = Float<Approximate,Precision64>;
using PositiveLowerFloat64 = Float<PositiveLower,Precision64>;
using PositiveUpperFloat64 = Float<PositiveUpper,Precision64>;
using PositiveExactFloat64 = Float<PositiveExact,Precision64>;
using ErrorFloat64 = PositiveUpperFloat64;
using ValidatedFloat64 = BoundedFloat64;
using BoundFloat64 = BoundedFloat64;
using MetrcFloat64 = MetricFloat64;
using ApprxFloat64 = ApproximateFloat64;

using ApproximateFloatMP = FloatMPTemplate<Approximate>;
using LowerFloatMP = FloatMPTemplate<Lower>;
using UpperFloatMP = FloatMPTemplate<Upper>;
using BoundedFloatMP= FloatMPTemplate<Bounded>;
using MetricFloatMP = FloatMPTemplate<Metric>;
using ExactFloatMP = FloatMPTemplate<Exact>;
using PositiveApproximateFloatMP = FloatMPTemplate<Approximate>;
using PositiveLowerFloatMP = FloatMPTemplate<PositiveLower>;
using PositiveUpperFloatMP = FloatMPTemplate<PositiveUpper>;
using PositiveExactFloatMP = FloatMPTemplate<PositiveExact>;
using ErrorFloatMP = PositiveUpperFloatMP;
using ValidatedFloatMP = BoundedFloatMP;
using BoundFloatMP = BoundedFloatMP;
using MetrcFloatMP = MetricFloatMP;
using ApprxFloatMP = ApproximateFloatMP;

/*
using ApproximateFloatMP = FloatMPTemplate<Approximate>;
using LowerFloatMP = FloatMPTemplate<Lower>;
using UpperFloatMP = FloatMPTemplate<Upper>;
using BoundedFloatMP = FloatMPTemplate<Bounded>;
using MetricFloatMP = FloatMPTemplate<Metric>;
using ExactFloatMP = FloatMPTemplate<Exact>;
using ErrorFloatMP = FloatMPTemplate<Error>;
using ValidatedFloatMP = BoundedFloatMP;
using PositiveUpperFloatMP = ErrorFloatMP;
using BoundFloatMP = BoundedFloatMP;
using MetrcFloatMP = MetricFloatMP;
using ApprxFloatMP = ApproximateFloatMP;
*/

template<class X> struct IsFloat : False { };
template<class P, class PR> struct IsFloat<Float<P,PR>> : True { };

template<class T> struct IsNumber;
template<> struct IsNumber<double> : True { };
template<class P, class PR> struct IsNumber<Float<P,PR>> : True { };

} // namespace Ariadne

#endif
