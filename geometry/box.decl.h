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

typedef Box<ExactFloatInterval> ExactFloatBox;
typedef Box<MetricFloatInterval> MetricFloatBox;
typedef Box<BoundedFloatInterval> BoundedFloatBox;
typedef Box<UpperFloatInterval> UpperFloatBox;
typedef Box<LowerFloatInterval> LowerFloatBox;
typedef Box<ApproximateFloatInterval> ApproximateFloatBox;

typedef Box<ExactInterval> ExactBox;
typedef Box<EffectiveInterval> EffectiveBox;
typedef Box<ValidatedInterval> ValidatedBox;
typedef Box<UpperInterval> UpperBox;
typedef Box<LowerInterval> LowerBox;
typedef Box<ApproximateInterval> ApproximateBox;

using BoxDomain = ExactFloatBox;

template<class X> class Vector;
/*
typedef Vector<ApproximateFloat> ApproximateFloatVector;
typedef Vector<ValidatedFloat> ValidatedFloatVector;
typedef Vector<ExactFloat> ExactFloatVector;
typedef Vector<UpperFloatInterval> UpperFloatIntervalVector;
*/
typedef Vector<ExactInterval> ExactIntervalVector;
typedef Vector<UpperInterval> UpperIntervalVector;
typedef Vector<ApproximateInterval> ApproximateIntervalVector;

} // namespace Ariadne


#endif /* ARIADNE_BOX_H */
