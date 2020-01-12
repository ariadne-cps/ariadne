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
//! \relates Box
//! \name Standard type names (deprecated)
typedef Box<ExactIntervalType> ExactBoxType; //!< .
typedef Box<EffectiveIntervalType> EffectiveBoxType; //!< .
typedef Box<ValidatedIntervalType> ValidatedBoxType; //!< .
typedef Box<UpperIntervalType> UpperBoxType; //!< .
typedef Box<LowerIntervalType> LowerBoxType; //!< .
typedef Box<ApproximateIntervalType> ApproximateBoxType; //!< .
//@}


//@{
//! \ingroup GeometryModule
//! \name Type definitions
typedef FloatDPExactBox BoxDomainType; //!< . \ingroup GeometryModule
typedef FloatDPUpperBox BoxValidatedRangeType; //!< . \ingroup GeometryModule
typedef FloatDPApproximateBox BoxApproximateRangeType; //!< . \ingroup GeometryModule
//@}

template<class IVL> class BoxSet;
//@{
//! \relates BoxSet
//! \name Type synonyms
typedef BoxSet<ExactIntervalType> ExactBoxSetType; //!< .
typedef BoxSet<ApproximateIntervalType> ApproximateBoxSetType; //!< .
typedef BoxSet<RealInterval> RealBoxSet; //!< .
//@}

} // namespace Ariadne


#endif /* ARIADNE_BOX_HPP */
