/***************************************************************************
 *            interval.decl.hpp
 *
 *  Copyright 2008-17  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ARIADNE_INTERVAL_DECL_HPP
#define ARIADNE_INTERVAL_DECL_HPP

#include "../numeric/float.decl.hpp"
#include "../numeric/number.decl.hpp"

namespace Ariadne {

template<class U> class Interval;
template<class U> class VariableInterval;
typedef Interval<Dyadic> DyadicInterval;
typedef Interval<Rational> RationalInterval;
typedef Interval<Real> RealInterval;

template<class F> using ExactInterval = Interval<Value<F>>;
template<class F> using BallInterval = Interval<Ball<F>> ;
template<class F> using BoundsInterval = Interval<Bounds<F>>;
template<class F> using UpperInterval = Interval<UpperBound<F>>;
template<class F> using LowerInterval = Interval<LowerBound<F>>;
template<class F> using ApproximateInterval = Interval<Approximation<F>>;

template<class PR> using FloatExactInterval = Interval<FloatValue<PR>>;
template<class PR> using FloatBallInterval = Interval<FloatBall<PR>> ;
template<class PR> using FloatBoundsInterval = Interval<FloatBounds<PR>>;
template<class PR> using FloatUpperInterval = Interval<FloatUpperBound<PR>>;
template<class PR> using FloatLowerInterval = Interval<FloatLowerBound<PR>>;
template<class PR> using FloatApproximateInterval = Interval<FloatApproximation<PR>>;

using FloatDPExactInterval = FloatExactInterval<DoublePrecision>;
using FloatDPBallInterval = FloatBallInterval<DoublePrecision> ;
using FloatDPBoundsInterval = FloatBoundsInterval<DoublePrecision>;
using FloatDPUpperInterval = FloatUpperInterval<DoublePrecision>;
using FloatDPLowerInterval = FloatLowerInterval<DoublePrecision>;
using FloatDPApproximateInterval = FloatApproximateInterval<DoublePrecision>;

using FloatMPExactInterval = FloatExactInterval<MultiplePrecision>;
using FloatMPBallInterval = FloatBallInterval<MultiplePrecision> ;
using FloatMPBoundsInterval = FloatBoundsInterval<MultiplePrecision>;
using FloatMPUpperInterval = FloatUpperInterval<MultiplePrecision>;
using FloatMPLowerInterval = FloatLowerInterval<MultiplePrecision>;
using FloatMPApproximateInterval = FloatApproximateInterval<MultiplePrecision>;


typedef Interval<ExactNumericType> ExactIntervalType;
typedef Interval<EffectiveNumericType> EffectiveIntervalType;
typedef Interval<ValidatedNumericType> ValidatedIntervalType;
typedef Interval<UpperNumericType> UpperIntervalType;
typedef Interval<LowerNumericType> LowerIntervalType;
typedef Interval<ApproximateNumericType> ApproximateIntervalType;


} // namespace Ariadne


#endif
