/***************************************************************************
 *            interval.decl.hpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

template<class IVL> struct IsInterval : False { };
template<class UB> struct IsInterval<Interval<UB>> : True { };


//@{
//! \relates Interval
//! \name Type shorthands
template<class F> using ExactInterval = Interval<Value<F>>; //!< .
template<class F> using BallInterval = Interval<Ball<F>>; //!< .
template<class F> using BoundsInterval = Interval<Bounds<F>>; //!< .
template<class F> using UpperInterval = Interval<UpperBound<F>>; //!< .
template<class F> using LowerInterval = Interval<LowerBound<F>>; //!< .
template<class F> using ApproximateInterval = Interval<Approximation<F>>; //!< .

template<class PR> using FloatExactInterval = Interval<FloatValue<PR>>; //!< .
template<class PR> using FloatBallInterval = Interval<FloatBall<PR>>; //!< .
template<class PR> using FloatBoundsInterval = Interval<FloatBounds<PR>>; //!< .
template<class PR> using FloatUpperInterval = Interval<FloatUpperBound<PR>>; //!< .
template<class PR> using FloatLowerInterval = Interval<FloatLowerBound<PR>>; //!< .
template<class PR> using FloatApproximateInterval = Interval<FloatApproximation<PR>>; //!< .
//@}

//@{
//! \relates Interval
//! \name Type synonyms
using DyadicInterval = Interval<Dyadic>; //!< .
using RationalInterval = Interval<Rational>; //!< .
using RealInterval = Interval<Real>; //!< .

using FloatDPExactInterval = FloatExactInterval<DoublePrecision>; //!< .
using FloatDPBallInterval = FloatBallInterval<DoublePrecision>; //!< .
using FloatDPBoundsInterval = FloatBoundsInterval<DoublePrecision>; //!< .
using FloatDPUpperInterval = FloatUpperInterval<DoublePrecision>; //!< .
using FloatDPLowerInterval = FloatLowerInterval<DoublePrecision>; //!< .
using FloatDPApproximateInterval = FloatApproximateInterval<DoublePrecision>; //!< .

using FloatMPExactInterval = FloatExactInterval<MultiplePrecision>; //!< .
using FloatMPBallInterval = FloatBallInterval<MultiplePrecision>; //!< .
using FloatMPBoundsInterval = FloatBoundsInterval<MultiplePrecision>; //!< .
using FloatMPUpperInterval = FloatUpperInterval<MultiplePrecision>; //!< .
using FloatMPLowerInterval = FloatLowerInterval<MultiplePrecision>; //!< .
using FloatMPApproximateInterval = FloatApproximateInterval<MultiplePrecision>; //!< .
//@}


//@{
//! \ingroup GeometryModule
//! \name Type definitions

//! The type used for the domain of a univariate function. \ingroup GeometryModule
typedef FloatDPExactInterval IntervalDomainType;
//! The type used for an over-approximation to the range of a validated scalar function. \ingroup GeometryModule
typedef FloatDPUpperInterval IntervalValidatedRangeType; //!< \ingroup GeometryModule .
//! The type used for an approximation to the range of an approximate scalar function. \ingroup GeometryModule
typedef FloatDPApproximateInterval IntervalApproximateRangeType; //!< \ingroup GeometryModule .

//! \brief The type used for testing properties of one-dimensional sets. \ingroup GeometryModule
typedef FloatDPExactInterval ExactIntervalType;
//! \brief The type used for the bounding interval of validated one-dimensional sets. \ingroup GeometryModule
typedef FloatDPUpperInterval UpperIntervalType;
//! \brief The type used for testing boundedness of validated one-dimensional sets. \ingroup GeometryModule
typedef FloatDPLowerInterval LowerIntervalType;
//! \brief The type used for the bounding interval of approximate one-dimensional sets. \ingroup GeometryModule
typedef FloatDPApproximateInterval ApproximateIntervalType;
//@}


} // namespace Ariadne


#endif
