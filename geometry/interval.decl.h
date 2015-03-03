/***************************************************************************
 *            interval.decl.h
 *
 *  Copyright 2008-13  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_INTERVAL_DECL_H
#define ARIADNE_INTERVAL_DECL_H

#include "numeric/float.decl.h"
#include "numeric/number.decl.h"

namespace Ariadne {

template<class U> class Interval;
typedef Interval<Rational> RationalInterval;
typedef Interval<Real> RealInterval;

typedef Interval<ExactFloat64> ExactFloat64Interval;
typedef Interval<MetricFloat64> MetricFloat64Interval;
typedef Interval<BoundedFloat64> BoundedFloat64Interval;
typedef Interval<UpperFloat64> UpperFloat64Interval;
typedef Interval<LowerFloat64> LowerFloat64Interval;
typedef Interval<ApproximateFloat64> ApproximateFloat64Interval;

typedef Interval<ExactNumericType> ExactIntervalType;
typedef Interval<EffectiveNumericType> EffectiveIntervalType;
typedef Interval<ValidatedNumericType> ValidatedIntervalType;
typedef Interval<UpperNumericType> UpperIntervalType;
typedef Interval<LowerNumericType> LowerIntervalType;
typedef Interval<ApproximateNumericType> ApproximateIntervalType;


} // namespace Ariadne


#endif
