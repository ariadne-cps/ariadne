/***************************************************************************
 *            algebra.cc
 *
 *  Copyright 2011-15  Pieter Collins
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

#include "numeric/numeric.h"
#include "config.h"

#include "algebra.h"
#include "algebra_mixin.h"
#include "algebra_operations.tcc"

namespace Ariadne {

template GradedAlgebra<ApproximateNumericType> compose(const Series<ApproximateNumericType>&, const GradedAlgebra<ApproximateNumericType>&);
template GradedAlgebra<ValidatedNumericType> compose(const Series<ValidatedNumericType>&, const GradedAlgebra<ValidatedNumericType>&);
template GradedAlgebra<EffectiveNumericType> compose(const Series<EffectiveNumericType>&, const GradedAlgebra<EffectiveNumericType>&);

template NormedAlgebra<ApproximateNumericType> rec(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> rec(const NormedAlgebra<ValidatedNumericType>&);

template NormedAlgebra<ApproximateNumericType> sqrt(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> sqrt(const NormedAlgebra<ValidatedNumericType>&);

template NormedAlgebra<ApproximateNumericType> exp(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> exp(const NormedAlgebra<ValidatedNumericType>&);

template NormedAlgebra<ApproximateNumericType> log(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> log(const NormedAlgebra<ValidatedNumericType>&);

template NormedAlgebra<ApproximateNumericType> sin(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> sin(const NormedAlgebra<ValidatedNumericType>&);

template NormedAlgebra<ApproximateNumericType> cos(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> cos(const NormedAlgebra<ValidatedNumericType>&);

template NormedAlgebra<ApproximateNumericType> tan(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> tan(const NormedAlgebra<ValidatedNumericType>&);

/*
template NormedAlgebra<ApproximateNumericType> asin(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> asin(const NormedAlgebra<ValidatedNumericType>&);
template NormedAlgebra<EffectiveNumericType> asin(const NormedAlgebra<EffectiveNumericType>&);

template NormedAlgebra<ApproximateNumericType> acos(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> acos(const NormedAlgebra<ValidatedNumericType>&);
template NormedAlgebra<EffectiveNumericType> acos(const NormedAlgebra<EffectiveNumericType>&);

template NormedAlgebra<ApproximateNumericType> atan(const NormedAlgebra<ApproximateNumericType>&);
template NormedAlgebra<ValidatedNumericType> atan(const NormedAlgebra<ValidatedNumericType>&);
template NormedAlgebra<EffectiveNumericType> atan(const NormedAlgebra<EffectiveNumericType>&);
*/

} // namespace Ariadne
