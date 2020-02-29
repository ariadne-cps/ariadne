/***************************************************************************
 *            function/taylor_function.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include <iostream>
#include <iomanip>

#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/algebra.hpp"
#include "../algebra/multi_index.hpp"
#include "../function/polynomial.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/evaluate.hpp"
#include "../function/taylor_model.hpp"

#include "../function/function.hpp"
#include "../function/function_mixin.hpp"
#include "../function/scaled_function_patch.hpp"

#include "../function/taylor_function.hpp"

#include "taylor_model.tpl.hpp"
#include "scaled_function_patch.tpl.hpp"
#include "function_mixin.tpl.hpp"

#define VOLATILE ;

namespace Ariadne {

template class ScaledFunctionPatchFactory<ValidatedTaylorModelDP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelDP>,BoxDomainType>;

template class ScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,BoxDomainType,IntervalDomainType>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,BoxDomainType,IntervalDomainType>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,BoxDomainType,BoxDomainType>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,BoxDomainType,BoxDomainType>;

template class ScaledFunctionPatchFactory<ValidatedTaylorModelMP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelMP>,BoxDomainType>;

template class ScaledFunctionPatch<ValidatedTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelMP>;


template class ScaledFunctionPatch<ValidatedIntervalTaylorModelDP>;
template class VectorScaledFunctionPatch<ValidatedIntervalTaylorModelDP>;
template class ScaledFunctionPatch<ValidatedIntervalTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedIntervalTaylorModelMP>;

FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper<FloatDP>());
}

FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory(Sweeper<FloatDP> const& sweeper) {
    return new TaylorFunctionFactory(sweeper);
}

} // namespace Ariadne
