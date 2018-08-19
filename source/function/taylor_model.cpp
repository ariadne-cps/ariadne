/***************************************************************************
 *            taylor_model.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "../numeric/numeric.hpp"
#include "../function/taylor_model.tpl.hpp"

namespace Ariadne {

template class Series<FloatDPBounds>;
template class Series<FloatMPBounds>;

template class SweeperBase<FloatDP>;
template class SweeperBase<FloatMP>;

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
