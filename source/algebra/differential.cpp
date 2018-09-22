/***************************************************************************
 *            differential.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../numeric/numeric.hpp"
#include "../geometry/interval.hpp"

#include "../algebra/differential.hpp"
#include "../algebra/univariate_differential.hpp"
#include "../algebra/fixed_differential.hpp"
#include "../algebra/fixed_univariate_differential.hpp"

#include "operations.hpp"

#include "differential.tpl.hpp"
#include "univariate_differential.tpl.hpp"

namespace Ariadne {

template class UnivariateDifferential<FloatDP>;
template class UnivariateDifferential<FloatDPApproximation>;
template class UnivariateDifferential<FloatDPBounds>;
template class UnivariateDifferential<UpperIntervalType>;

template class UnivariateDifferential<FloatMPApproximation>;
template class UnivariateDifferential<FloatMPBounds>;


template class Differential<FloatDP>;
template class Differential<FloatDPBounds>;
template class Differential<FloatDPApproximation>;
template class Differential<UpperIntervalType>;

template struct AlgebraOperations<Differential<FloatDP>>;
template struct AlgebraOperations<Differential<FloatDPApproximation>>;
template struct AlgebraOperations<Differential<FloatDPBounds>>;
template class GradedAlgebraOperations<Differential<FloatDP>>;
template class GradedAlgebraOperations<Differential<FloatDPApproximation>>;
template class GradedAlgebraOperations<Differential<FloatDPBounds>>;

template class Vector<Differential<FloatDP>>;
template class Vector<Differential<FloatDPBounds>>;
template class Vector<Differential<FloatDPApproximation>>;
//template class Vector<Differential<UpperIntervalType>>;

template class Differential<FloatMPBounds>;
template class Differential<FloatMPApproximation>;
template struct AlgebraOperations<Differential<FloatMPApproximation>>;
template struct AlgebraOperations<Differential<FloatMPBounds>>;
template class GradedAlgebraOperations<Differential<FloatMPApproximation>>;
template class GradedAlgebraOperations<Differential<FloatMPBounds>>;

template class Vector<Differential<FloatMPBounds>>;
template class Vector<Differential<FloatMPApproximation>>;

}
