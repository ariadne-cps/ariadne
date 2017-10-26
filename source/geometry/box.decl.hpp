/***************************************************************************
 *            box.decl.hpp
 *
 *  Copyright 2008-17  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_BOX_DECL_HPP
#define ARIADNE_BOX_DECL_HPP

#include "geometry/interval.decl.hpp"

namespace Ariadne {

template<class IVL> class Box;
template<class IVL> class VariablesBox;

typedef Box<RationalInterval> RationalBox;
typedef Box<RealInterval> RealBox;

typedef Box<FloatDPExactInterval> FloatDPExactBox;
typedef Box<FloatDPBallInterval> FloatDPBallBox;
typedef Box<FloatDPBoundsInterval> FloatDPBoundsBox;
typedef Box<FloatDPUpperInterval> FloatDPUpperBox;
typedef Box<FloatDPLowerInterval> FloatDPLowerBox;
typedef Box<FloatDPApproximateInterval> FloatDPApproximateBox;

typedef Box<ExactIntervalType> ExactBoxType;
typedef Box<EffectiveIntervalType> EffectiveBoxType;
typedef Box<ValidatedIntervalType> ValidatedBoxType;
typedef Box<UpperIntervalType> UpperBoxType;
typedef Box<LowerIntervalType> LowerBoxType;
typedef Box<ApproximateIntervalType> ApproximateBoxType;

template<class IVL> class BoxSet;
typedef BoxSet<ExactIntervalType> ExactBoxSet;
typedef BoxSet<ApproximateIntervalType> ApproximateBoxSet;
typedef BoxSet<RealInterval> RealBoxSet;

using BoxDomainType = FloatDPExactBox;

template<class X> class Vector;

typedef Vector<ExactIntervalType> ExactIntervalVectorType;
typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Vector<ApproximateIntervalType> ApproximateIntervalVectorType;

} // namespace Ariadne


#endif /* ARIADNE_BOX_HPP */
