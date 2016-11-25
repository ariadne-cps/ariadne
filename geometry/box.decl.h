/***************************************************************************
 *            box.decl.h
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

#ifndef ARIADNE_BOX_DECL_H
#define ARIADNE_BOX_DECL_H

#include "geometry/interval.decl.h"

namespace Ariadne {

template<class IVL> class Box;
typedef Box<RationalInterval> RationalBox;
typedef Box<RealInterval> RealBox;

typedef Box<Float64ExactInterval> Float64ExactBox;
typedef Box<Float64BallInterval> Float64BallBox;
typedef Box<Float64BoundsInterval> Float64BoundsBox;
typedef Box<Float64UpperInterval> Float64UpperBox;
typedef Box<Float64LowerInterval> Float64LowerBox;
typedef Box<Float64ApproximateInterval> Float64ApproximateBox;

typedef Box<ExactIntervalType> ExactBoxType;
typedef Box<EffectiveIntervalType> EffectiveBoxType;
typedef Box<ValidatedIntervalType> ValidatedBoxType;
typedef Box<UpperIntervalType> UpperBoxType;
typedef Box<LowerIntervalType> LowerBoxType;
typedef Box<ApproximateIntervalType> ApproximateBoxType;

using BoxDomain = Float64ExactBox;

template<class X> class Vector;

typedef Vector<ExactIntervalType> ExactIntervalVectorType;
typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Vector<ApproximateIntervalType> ApproximateIntervalVectorType;

} // namespace Ariadne


#endif /* ARIADNE_BOX_H */
