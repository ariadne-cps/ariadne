/***************************************************************************
 *            algebra.cc
 *
 *  Copyright 2011  Pieter Collins
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

#include "numeric.h"
#include "config.h"

#include "algebra_mixin.tcc"

namespace Ariadne {

template GradedAlgebra<ApproximateNumberType> compose(const Series<ApproximateNumberType>&, const GradedAlgebra<ApproximateNumberType>&);
template GradedAlgebra<ValidatedNumberType> compose(const Series<ValidatedNumberType>&, const GradedAlgebra<ValidatedNumberType>&);
template GradedAlgebra<EffectiveNumberType> compose(const Series<EffectiveNumberType>&, const GradedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> rec(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> rec(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> rec(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> sqrt(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> sqrt(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> sqrt(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> exp(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> exp(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> exp(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> log(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> log(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> log(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> sin(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> sin(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> sin(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> cos(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> cos(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> cos(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> tan(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> tan(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> tan(const NormedAlgebra<EffectiveNumberType>&);

/*
template NormedAlgebra<ApproximateNumberType> asin(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> asin(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> asin(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> acos(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> acos(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> acos(const NormedAlgebra<EffectiveNumberType>&);

template NormedAlgebra<ApproximateNumberType> atan(const NormedAlgebra<ApproximateNumberType>&);
template NormedAlgebra<ValidatedNumberType> atan(const NormedAlgebra<ValidatedNumberType>&);
template NormedAlgebra<EffectiveNumberType> atan(const NormedAlgebra<EffectiveNumberType>&);
*/

} // namespace Ariadne
