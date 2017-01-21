/***************************************************************************
 *            taylor_function.cc
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
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "function/functional.h"
#include "config.h"

#include <iostream>
#include <iomanip>

#include "utility/macros.h"
#include "utility/exceptions.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/algebra.h"
#include "algebra/multi_index.h"
#include "function/polynomial.h"
#include "algebra/differential.h"
#include "algebra/evaluate.h"
#include "function/taylor_model.h"

#include "function/function.h"
#include "function/function_mixin.h"
#include "function/function_patch.h"

#include "function/taylor_function.h"

#include "taylor_model.tpl.h"
#include "function_patch.tpl.h"
#include "function_mixin.tpl.h"

#define VOLATILE ;

namespace Ariadne {

static double TAYLOR_FUNCTION_WRITING_ACCURACY = 1e-8;

template class FunctionPatchFactory<ValidatedTaylorModel64>;
template class FunctionModelCreator<FunctionPatchFactory<ValidatedTaylorModel64>>;

template class FunctionPatch<ValidatedTaylorModel64>;
template class FunctionMixin<FunctionPatch<ValidatedTaylorModel64>,ApproximateTag,BoxDomain,IntervalDomain>;
template class FunctionMixin<FunctionPatch<ValidatedTaylorModel64>,ValidatedTag,BoxDomain,IntervalDomain>;
template class VectorFunctionPatch<ValidatedTaylorModel64>;
template class FunctionMixin<VectorFunctionPatch<ValidatedTaylorModel64>,ApproximateTag,BoxDomain,BoxDomain>;
template class FunctionMixin<VectorFunctionPatch<ValidatedTaylorModel64>,ValidatedTag,BoxDomain,BoxDomain>;

template class FunctionPatchFactory<ValidatedTaylorModelMP>;
template class FunctionModelCreator<FunctionPatchFactory<ValidatedTaylorModelMP>>;
template class FunctionPatch<ValidatedTaylorModelMP>;
template class VectorFunctionPatch<ValidatedTaylorModelMP>;

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper<Float64>());
}

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold) {
    return new TaylorFunctionFactory(ThresholdSweeper<Float64>(Precision64(),sweep_threshold));
}


} // namespace Ariadne
