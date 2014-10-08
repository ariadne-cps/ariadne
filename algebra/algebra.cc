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

#include "numeric/numeric.h"
#include "config.h"

#include "algebra_mixin.tcc"

namespace Ariadne {

template GradedAlgebra<ApproximateNumber> compose(const Series<ApproximateNumber>&, const GradedAlgebra<ApproximateNumber>&);
template GradedAlgebra<ValidatedNumber> compose(const Series<ValidatedNumber>&, const GradedAlgebra<ValidatedNumber>&);
template GradedAlgebra<EffectiveNumber> compose(const Series<EffectiveNumber>&, const GradedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> rec(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> rec(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> rec(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> sqrt(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> sqrt(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> sqrt(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> exp(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> exp(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> exp(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> log(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> log(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> log(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> sin(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> sin(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> sin(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> cos(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> cos(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> cos(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> tan(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> tan(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> tan(const NormedAlgebra<EffectiveNumber>&);

/*
template NormedAlgebra<ApproximateNumber> asin(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> asin(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> asin(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> acos(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> acos(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> acos(const NormedAlgebra<EffectiveNumber>&);

template NormedAlgebra<ApproximateNumber> atan(const NormedAlgebra<ApproximateNumber>&);
template NormedAlgebra<ValidatedNumber> atan(const NormedAlgebra<ValidatedNumber>&);
template NormedAlgebra<EffectiveNumber> atan(const NormedAlgebra<EffectiveNumber>&);
*/

} // namespace Ariadne
