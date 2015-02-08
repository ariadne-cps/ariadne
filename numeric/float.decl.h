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

class Float64;
class FloatMP;
class TwoExp;

class Precision64;
class PrecisionMP;

typedef Float64 Float64;
typedef Precision64 Precision64;

using RawFloat64 = Float64;
class ExactFloat64;
class ValidatedFloat64;
class UpperFloat64;
class LowerFloat64;
class ApproximateFloat64;
class PositiveUpperFloat64;
using MetricFloat64 = ValidatedFloat64;
using BoundedFloat64 = ValidatedFloat64;
using ErrorFloat64 = PositiveUpperFloat64;
using BoundFloat64 = BoundedFloat64;
using MetrcFloat64 = MetricFloat64;
using ApprxFloat64 = ApproximateFloat64;

template<class P> class FloatTemplate;

template<class P> class Float64Template;

template<class P> class FloatMPTemplate;
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

template<> struct IsNumber<double> : True { };

template<class X> struct IsFloat : False { };
template<> struct IsFloat<double> : True { };
template<> struct IsFloat<ExactFloat64> : True { };
template<> struct IsFloat<ValidatedFloat64> : True { };
template<> struct IsFloat<UpperFloat64> : True { };
template<> struct IsFloat<LowerFloat64> : True { };
template<> struct IsFloat<ApproximateFloat64> : True { };

template<class T> struct IsNumber;
template<> struct IsNumber<ExactFloat64> : True { };
template<> struct IsNumber<ValidatedFloat64> : True { };
template<> struct IsNumber<UpperFloat64> : True { };
template<> struct IsNumber<LowerFloat64> : True { };
template<> struct IsNumber<ApproximateFloat64> : True { };

} // namespace Ariadne

#endif
