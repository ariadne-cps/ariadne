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
template<class U> class VariableInterval;
typedef Interval<Rational> RationalInterval;
typedef Interval<Real> RealInterval;

template<class PR> using FloatExactInterval = Interval<FloatValue<PR>>;
template<class PR> using FloatBallInterval = Interval<FloatBall<PR>> ;
template<class PR> using FloatBoundsInterval = Interval<FloatBounds<PR>>;
template<class PR> using FloatUpperInterval = Interval<FloatUpperBound<PR>>;
template<class PR> using FloatLowerInterval = Interval<FloatLowerBound<PR>>;
template<class PR> using FloatApproximateInterval = Interval<FloatApproximation<PR>>;

using Float64ExactInterval = FloatExactInterval<Precision64>;
using Float64BallInterval = FloatBallInterval<Precision64> ;
using Float64BoundsInterval = FloatBoundsInterval<Precision64>;
using Float64UpperInterval = FloatUpperInterval<Precision64>;
using Float64LowerInterval = FloatLowerInterval<Precision64>;
using Float64ApproximateInterval = FloatApproximateInterval<Precision64>;

using FloatMPExactInterval = FloatExactInterval<PrecisionMP>;
using FloatMPBallInterval = FloatBallInterval<PrecisionMP> ;
using FloatMPBoundsInterval = FloatBoundsInterval<PrecisionMP>;
using FloatMPUpperInterval = FloatUpperInterval<PrecisionMP>;
using FloatMPLowerInterval = FloatLowerInterval<PrecisionMP>;
using FloatMPApproximateInterval = FloatApproximateInterval<PrecisionMP>;


typedef Interval<ExactNumericType> ExactIntervalType;
typedef Interval<EffectiveNumericType> EffectiveIntervalType;
typedef Interval<ValidatedNumericType> ValidatedIntervalType;
typedef Interval<UpperNumericType> UpperIntervalType;
typedef Interval<LowerNumericType> LowerIntervalType;
typedef Interval<ApproximateNumericType> ApproximateIntervalType;


} // namespace Ariadne


#endif
