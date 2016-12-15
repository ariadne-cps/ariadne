/***************************************************************************
 *            affine_model.cc
 *
 *  Copyright 2009-16  Pieter Collins
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
#include "config.h"

#include "algebra/vector.h"
#include "function/function.h"
#include "function/taylor_model.h"
#include "function/affine_model.h"

#include "function/affine.h"
#include "function/taylor_function.h"
#include "algebra/vector.h"

#include "function/affine_model.tcc"


namespace Ariadne {

template struct AffineModel<ApproximateTag,Float64>;
template struct AffineModel<ApproximateTag,FloatMP>;

template struct AffineModel<ValidatedTag,Float64>;
template struct AffineModel<ValidatedTag,FloatMP>;

template struct AlgebraOperations<AffineModel<ApproximateTag,Float64>,FloatBounds<Precision64>>;
template struct AlgebraOperations<AffineModel<ApproximateTag,FloatMP>,FloatBounds<PrecisionMP>>;

template struct AlgebraOperations<AffineModel<ValidatedTag,Float64>,FloatBounds<Precision64>>;
template struct AlgebraOperations<AffineModel<ValidatedTag,FloatMP>,FloatBounds<PrecisionMP>>;

template AffineModel<ValidatedTag,Float64> affine_model(const ExactBoxType&, const ScalarFunction<ValidatedTag>&, Precision64);

} //namespace Ariadne


