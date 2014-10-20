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

template<class P> class Float64;
using RawFloat64 = Flt64;
using ExactFloat64 = Float64<Exact>;
using ErrorFloat64 = Float64<Error>;
using MetrcFloat64 = Float64<Metrc>;
using BoundFloat64 = Float64<Bound>;
using UpperFloat64 = Float64<Upper>;
using LowerFloat64 = Float64<Lower>;
using ApprxFloat64 = Float64<Apprx>;

using ValidFloat64 = Float64<Bounded>;

using MetricFloat64 = Float64<Metric>;
using BoundedFloat64 = Float64<Bounded>;
using ValidatedFloat64 = Float64<Bounded>;
using ApproximateFloat64 = Float64<Approximate>;

template<class P> class FloatMP;
using RawFloatMP = FltMP;
using ExactFloatMP = FloatMP<Exact>;
using ErrorFloatMP = FloatMP<Error>;
using MetrcFloatMP = FloatMP<Metrc>;
using BoundFloatMP = FloatMP<Bound>;
using UpperFloatMP = FloatMP<Upper>;
using LowerFloatMP = FloatMP<Lower>;
using ApprxFloatMP = FloatMP<Apprx>;

using ValidFloatMP = FloatMP<Bounded>;

using MetricFloatMP = FloatMP<Metric>;
using BoundedFloatMP = FloatMP<Bounded>;
using ValidatedFloatMP = FloatMP<Bounded>;
using ApproximateFloatMP = FloatMP<Approximate>;

template<class P> using Float = Float64<P>;
typedef Flt64 Flt;
typedef RawFloat64 RawFloat;
typedef ExactFloat64 ExactFloat;
typedef ErrorFloat64 ErrorFloat;
typedef MetrcFloat64 MetrcFloat;
typedef BoundFloat64 BoundFloat;
typedef UpperFloat64 UpperFloat;
typedef LowerFloat64 LowerFloat;
typedef ApprxFloat64 ApprxFloat;

using ValidFloat = Float<Bounded>;

typedef MetricFloat64 MetricFloat;
typedef BoundedFloat64 BoundedFloat;
typedef ValidatedFloat64 ValidatedFloat;
typedef ApproximateFloat64 ApproximateFloat;

template<class X> struct IsFloat : False { };
template<class P> struct IsFloat<Float64<P>> : True { };
template<class P> struct IsFloat<FloatMP<P>> : True { };
template<class T> struct IsNumber;
template<class P> struct IsNumber<Float64<P>> : True { };
template<class P> struct IsNumber<FloatMP<P>> : True { };

template<class P> struct Float64Typedef { typedef Float64<P> Type; };
template<> struct Float64Typedef<Effective> { typedef Float64<Metric> Type; };
template<> struct Float64Typedef<Validated> { typedef Float64<Metric> Type; };
template<class P> using Float64Type = typename Float64Typedef<P>::Type;

template<class P> struct FloatMPTypedef { typedef FloatMP<P> Type; };
template<> struct FloatMPTypedef<Effective> { typedef FloatMP<Validated> Type; };
template<class P> using FloatMPType = typename FloatMPTypedef<P>::Type;

template<class P, class A> struct FloatTypedef;
template<class P> struct FloatTypedef<P,Precision64> { typedef Float64Type<P> Type; };
template<class P> struct FloatTypedef<P,PrecisionMP> { typedef FloatMPType<P> Type; };
template<class P, class A=Precision64> using FloatType = typename FloatTypedef<P,A>::Type;

} // namespace Ariadne

#endif
