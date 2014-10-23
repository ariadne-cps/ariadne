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

class Flt64;
class FltMP;
class TwoExp;

class Precision64;
class PrecisionMP;

template<class P> class Float64Template;
using Float64 = Flt64;
using RawFloat64 = Float64;
using ExactFloat64 = Float64Template<Exact>;
using ErrorFloat64 = Float64Template<Error>;
using MetricFloat64 = Float64Template<Metric>;
using BoundedFloat64 = Float64Template<Bounded>;
using UpperFloat64 = Float64Template<Upper>;
using LowerFloat64 = Float64Template<Lower>;
using ApproximateFloat64 = Float64Template<Approximate>;

using ValidatedFloat64 = BoundedFloat64;

using MetrcFloat64 = MetricFloat64;
using BoundFloat64 = BoundedFloat64;
using ValidFloat64 = ValidatedFloat64;
using ApprxFloat64 = ApproximateFloat64;

template<class P> class FloatMPTemplate;
using FloatMP = FltMP;
using RawFloatMP = FltMP;
using ExactFloatMP = FloatMPTemplate<Exact>;
using ErrorFloatMP = FloatMPTemplate<Error>;
using MetrcFloatMP = FloatMPTemplate<Metrc>;
using BoundFloatMP = FloatMPTemplate<Bound>;
using UpperFloatMP = FloatMPTemplate<Upper>;
using LowerFloatMP = FloatMPTemplate<Lower>;
using ApprxFloatMP = FloatMPTemplate<Apprx>;

using ValidFloatMP = FloatMPTemplate<Bounded>;

using MetricFloatMP = FloatMPTemplate<Metric>;
using BoundedFloatMP = FloatMPTemplate<Bounded>;
using ValidatedFloatMP = FloatMPTemplate<Bounded>;
using ApproximateFloatMP = FloatMPTemplate<Approximate>;

class Float;
using RawFloat = Float;
class ExactFloat;
class ValidatedFloat;
class UpperFloat;
class LowerFloat;
class ApproximateFloat;
using MetricFloat = ValidatedFloat;
using BoundedFloat = ValidatedFloat;
using PositiveUpperFloat = UpperFloat;

typedef PositiveUpperFloat ErrorFloat;

/*
typedef RawFloat64 RawFloat;
typedef ExactFloat64 ExactFloat;
typedef ErrorFloat64 ErrorFloat;
typedef MetrcFloat64 MetricFloat;
typedef BoundFloat64 BoundedFloat;
typedef UpperFloat64 UpperFloat;
typedef LowerFloat64 LowerFloat;
typedef ApprxFloat64 ApproximateFloat;
*/

typedef MetricFloat MetrcFloat;
typedef BoundedFloat BoundFloat;
typedef ValidatedFloat ValidFloat;
typedef ApproximateFloat ApprxFloat;

template<class X> struct IsFloat : False { };
template<> struct IsFloat<ExactFloat> : True { };
template<> struct IsFloat<ValidatedFloat> : True { };
template<> struct IsFloat<UpperFloat> : True { };
template<> struct IsFloat<LowerFloat> : True { };
template<> struct IsFloat<ApproximateFloat> : True { };
template<class P> struct IsFloat<Float64Template<P>> : True { };
template<class P> struct IsFloat<FloatMPTemplate<P>> : True { };
template<class T> struct IsNumber;
template<> struct IsNumber<ExactFloat> : True { };
template<> struct IsNumber<ValidatedFloat> : True { };
template<> struct IsNumber<UpperFloat> : True { };
template<> struct IsNumber<LowerFloat> : True { };
template<> struct IsNumber<ApproximateFloat> : True { };
template<class P> struct IsNumber<FloatMPTemplate<P>> : True { };
template<class P> struct IsNumber<Float64Template<P>> : True { };

template<class P> struct Float64Typedef;
template<> struct Float64Typedef<Exact> { typedef ExactFloat64 Type; };
template<> struct Float64Typedef<Effective> { typedef ValidatedFloat64 Type; };
template<> struct Float64Typedef<Validated> { typedef ValidatedFloat64 Type; };
template<> struct Float64Typedef<Approximate> { typedef ApproximateFloat64 Type; };
template<class P> using Float64Type = typename Float64Typedef<P>::Type;

template<class P> struct FloatMPTypedef { typedef FloatMPTemplate<P> Type; };
template<> struct FloatMPTypedef<Effective> { typedef FloatMPTemplate<Validated> Type; };
template<class P> using FloatMPType = typename FloatMPTypedef<P>::Type;

template<class P, class A> struct FloatTypedef;
template<class P> struct FloatTypedef<P,Precision64> { typedef Float64Type<P> Type; };
template<class P> struct FloatTypedef<P,PrecisionMP> { typedef FloatMPType<P> Type; };
template<class P, class A=Precision64> using FloatType = typename FloatTypedef<P,A>::Type;

} // namespace Ariadne

#endif
