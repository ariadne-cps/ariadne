/***************************************************************************
 *            taylor_model.cpp
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
#include "../function/taylor_model.tpl.hpp"

namespace Ariadne {

template class Series<FloatDPBounds>;
template class Series<FloatMPBounds>;

template class SweeperBase<FloatDP>;
template class SweeperBase<FloatMP>;
template class RelativeSweeperBase<FloatDP>;
template class RelativeSweeperBase<FloatMP>;

template class TaylorModel<ValidatedTag,FloatDP>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatDP>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatDP>>;

template class TaylorModel<ApproximateTag,FloatDP>;
template struct AlgebraOperations<TaylorModel<ApproximateTag,FloatDP>>;
template class NormedAlgebraOperations<TaylorModel<ApproximateTag,FloatDP>>;


template class TaylorModel<ValidatedTag,FloatMP>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatMP>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatMP>>;

template class TaylorModel<ApproximateTag,FloatMP>;
template struct AlgebraOperations<TaylorModel<ApproximateTag,FloatMP>>;
template class NormedAlgebraOperations<TaylorModel<ApproximateTag,FloatMP>>;

} // namespace Ariadne
