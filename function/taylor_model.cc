/***************************************************************************
 *            taylor_model.cc
 *
 *  Copyright 2008-15  Pieter Collins
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

#include "numeric/numeric.h"
#include "function/taylor_model.tcc"

namespace Ariadne {

template class Series<Float64Bounds>;
template class Series<FloatMPBounds>;

template class TaylorModel<ValidatedTag,Float64>;
template class AlgebraOperations<TaylorModel<ValidatedTag,Float64>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,Float64>>;

template class TaylorModel<ApproximateTag,Float64>;
template class AlgebraOperations<TaylorModel<ApproximateTag,Float64>>;
template class NormedAlgebraOperations<TaylorModel<ApproximateTag,Float64>>;


template class TaylorModel<ValidatedTag,FloatMP>;
template class AlgebraOperations<TaylorModel<ValidatedTag,FloatMP>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatMP>>;

template class TaylorModel<ApproximateTag,FloatMP>;
template class AlgebraOperations<TaylorModel<ApproximateTag,FloatMP>>;
template class NormedAlgebraOperations<TaylorModel<ApproximateTag,FloatMP>>;

} // namespace Ariadne
