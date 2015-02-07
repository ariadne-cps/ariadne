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

template<class X> class Positive;

class Float64;
class FloatMP;
class TwoExp;

class Precision64;
class PrecisionMP;

typedef Float64 Float;
typedef Precision64 Precision;

template<class PR> struct PrecisionTraits;
template<> struct PrecisionTraits<Precision64> { typedef Float64 Type; };
template<> struct PrecisionTraits<PrecisionMP> { typedef FloatMP Type; };
template<class PR> using FloatType = typename PrecisionTraits<PR>::Type;
template<class PR> using RawFloatType = typename PrecisionTraits<PR>::Type;

template<class P, class PR> class FloatTemplate;
template<class P> using Float64Template = FloatTemplate<P,Precision64>;
template<class P> using FloatMPTemplate = FloatTemplate<P,PrecisionMP>;

using RawFloat = Float;
using ExactFloat = FloatTemplate<Exact,Precision>;
using ValidatedFloat = FloatTemplate<Validated,Precision>;
using UpperFloat = FloatTemplate<Upper,Precision>;
using LowerFloat = FloatTemplate<Lower,Precision>;
using ApproximateFloat = FloatTemplate<Approximate,Precision>;
using PositiveUpperFloat = FloatTemplate<PositiveUpper,Precision>;

//using PositiveExactFloat = Positive<ExactFloat>;
//using PositiveUpperFloat = Positive<UpperFloat>;
//using PositiveLowerFloat = Positive<LowerFloat>;
//using PositiveApproximateFloat = Positive<ApproximateFloat>;
class PositiveExactFloat;
class PositiveLowerFloat;
class PositiveApproximateFloat;

using ErrorFloat = PositiveUpperFloat;

using MetricFloat = ValidatedFloat;
using BoundedFloat = ValidatedFloat;
using BoundFloat = BoundedFloat;
using MetrcFloat = MetricFloat;
using ApprxFloat = ApproximateFloat;


using RawFloat64 = Float64;
using ExactFloat64 = FloatTemplate<Exact,Precision64>;
using ValidatedFloat64 = FloatTemplate<Validated,Precision64>;
using UpperFloat64 = FloatTemplate<Upper,Precision64>;
using LowerFloat64 = FloatTemplate<Lower,Precision64>;
using ApproximateFloat64 = FloatTemplate<Approximate,Precision64>;
using PositiveUpperFloat64 = FloatTemplate<PositiveUpper,Precision64>;
using MetricFloat64 = ValidatedFloat64;
using BoundedFloat64 = ValidatedFloat64;
using ErrorFloat64 = PositiveUpperFloat64;
using BoundFloat64 = BoundedFloat64;
using MetrcFloat64 = MetricFloat64;
using ApprxFloat64 = ApproximateFloat64;

using RawFloatMP = FloatMP;
using ExactFloatMP = FloatTemplate<Exact,PrecisionMP>;
using ValidatedFloatMP = FloatTemplate<Validated,PrecisionMP>;
using UpperFloatMP = FloatTemplate<Upper,PrecisionMP>;
using LowerFloatMP = FloatTemplate<Lower,PrecisionMP>;
using ApproximateFloatMP = FloatTemplate<Approximate,PrecisionMP>;
using PositiveUpperFloatMP = FloatTemplate<PositiveUpper,PrecisionMP>;
using MetricFloatMP = ValidatedFloatMP;
using BoundedFloatMP = ValidatedFloatMP;
using ErrorFloatMP = PositiveUpperFloatMP;
using BoundFloatMP = BoundedFloatMP;
using MetrcFloatMP = MetricFloatMP;
using ApprxFloatMP = ApproximateFloatMP;


template<> struct IsNumber<double> : True { };

template<class X> struct IsFloat : False { };
template<> struct IsFloat<double> : True { };
template<> struct IsFloat<ExactFloat> : True { };
template<> struct IsFloat<ValidatedFloat> : True { };
template<> struct IsFloat<UpperFloat> : True { };
template<> struct IsFloat<LowerFloat> : True { };
template<> struct IsFloat<ApproximateFloat> : True { };

template<class T> struct IsNumber;
template<> struct IsNumber<ExactFloat> : True { };
template<> struct IsNumber<ValidatedFloat> : True { };
template<> struct IsNumber<UpperFloat> : True { };
template<> struct IsNumber<LowerFloat> : True { };
template<> struct IsNumber<ApproximateFloat> : True { };

} // namespace Ariadne

#endif
