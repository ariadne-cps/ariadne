/***************************************************************************
 *            taylor_function.cpp
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
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "function/functional.hpp"
#include "config.hpp"

#include <iostream>
#include <iomanip>

#include "utility/macros.hpp"
#include "utility/exceptions.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "algebra/multi_index.hpp"
#include "function/polynomial.hpp"
#include "algebra/differential.hpp"
#include "algebra/evaluate.hpp"
#include "function/taylor_model.hpp"

#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/scaled_function_patch.hpp"

#include "function/taylor_function.hpp"

#include "taylor_model.tpl.hpp"
#include "scaled_function_patch.tpl.hpp"
#include "function_mixin.tpl.hpp"

#define VOLATILE ;

namespace Ariadne {

static double TAYLOR_FUNCTION_WRITING_ACCURACY = 1e-8;

template class ScaledFunctionPatchFactory<ValidatedTaylorModelDP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelDP>>;

template class ScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,BoxDomainType,IntervalDomainType>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,BoxDomainType,IntervalDomainType>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,BoxDomainType,BoxDomainType>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,BoxDomainType,BoxDomainType>;

template class ScaledFunctionPatchFactory<ValidatedTaylorModelMP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelMP>>;
template class ScaledFunctionPatch<ValidatedTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelMP>;

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper<FloatDP>());
}

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold) {
    return new TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,sweep_threshold));
}


} // namespace Ariadne
