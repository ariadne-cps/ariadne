/***************************************************************************
 *            box.decl.hpp
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

#ifndef ARIADNE_BOX_DECL_HPP
#define ARIADNE_BOX_DECL_HPP

#include "../geometry/interval.decl.hpp"

namespace Ariadne {

template<class IVL> class Box;
template<class IVL> class VariablesBox;

//@{
//! \relates Box
//! \name Type synonyms
using DyadicBox = Box<DyadicInterval>; //!< .
using RationalBox = Box<RationalInterval>; //!< .
using RealBox = Box<RealInterval>; //!< .

using FloatDPExactBox = Box<FloatDPExactInterval>; //!< .
using FloatDPBallBox = Box<FloatDPBallInterval>; //!< .
using FloatDPBoundsBox = Box<FloatDPBoundsInterval>; //!< .
using FloatDPUpperBox = Box<FloatDPUpperInterval>; //!< .
using FloatDPLowerBox = Box<FloatDPLowerInterval>; //!< .
using FloatDPApproximateBox = Box<FloatDPApproximateInterval>; //!< .
//@}


//@{
//! \ingroup GeometryModule
//! \name Type definitions
//! The type used for the domain of a multivariate function. \ingroup GeometryModule
typedef FloatDPExactBox BoxDomainType;
//! The type used for an over-approximation to the range of a validated vector function. \ingroup GeometryModule
typedef FloatDPUpperBox BoxValidatedRangeType; //!< \ingroup GeometryModule .
//! The type used for an approximation to the range of an approximate vector function. \ingroup GeometryModule
typedef FloatDPApproximateBox BoxApproximateRangeType; //!< \ingroup GeometryModule .

//! \brief The type used for the bounding box of validated sets in Euclidean space. \ingroup GeometryModule
typedef FloatDPUpperBox BoundingBoxType;

//! \brief The type used for testing properties of sets in Euclidean space. \ingroup GeometryModule
typedef FloatDPExactBox ExactBoxType;
//! \brief The type used for the bounding box of validated sets in Euclidean space. \ingroup GeometryModule
typedef FloatDPUpperBox UpperBoxType;
//! \brief The type used for testing boundedness of sets in Euclidean space. \ingroup GeometryModule
typedef FloatDPLowerBox LowerBoxType;
//! \brief The type used for the bounding box of approximate sets in Euclidean space. \ingroup GeometryModule
typedef FloatDPApproximateBox ApproximateBoxType;
//@}

template<class IVL> class BoxSet;
//@{
//! \relates BoxSet
//! \name Type synonyms
typedef BoxSet<ExactIntervalType> ExactBoxSetType; //!< .
typedef BoxSet<UpperIntervalType> UpperBoxSetType; //!< .
typedef BoxSet<ApproximateIntervalType> ApproximateBoxSetType; //!< .
typedef BoxSet<RealInterval> RealBoxSet; //!< .
//@}

} // namespace Ariadne


#endif /* ARIADNE_BOX_HPP */
