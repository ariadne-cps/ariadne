/***************************************************************************
 *            affine_model.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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
#include "../config.hpp"

#include "../algebra/algebra.hpp"
#include "../algebra/vector.hpp"
#include "../function/function.hpp"
#include "../function/formula.hpp"
#include "../function/taylor_model.hpp"
#include "../function/affine_model.hpp"

#include "../function/affine.hpp"
#include "../function/taylor_function.hpp"

#include "../algebra/expansion.inl.hpp"
#include "../function/affine_model.tpl.hpp"


namespace Ariadne {

template class AffineModel<ApproximateTag,FloatDP>;
template class AffineModel<ApproximateTag,FloatMP>;

template class AffineModel<ValidatedTag,FloatDP>;
template class AffineModel<ValidatedTag,FloatMP>;

template class AlgebraOperations<AffineModel<ApproximateTag,FloatDP>>;
template class AlgebraOperations<AffineModel<ApproximateTag,FloatMP>>;

template class AlgebraOperations<AffineModel<ValidatedTag,FloatDP>>;
template class AlgebraOperations<AffineModel<ValidatedTag,FloatMP>>;

template AffineModel<ValidatedTag,FloatDP> affine_model(const ExactBoxType&, const ScalarFunction<ValidatedTag>&, DoublePrecision);

} //namespace Ariadne


