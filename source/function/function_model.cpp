/***************************************************************************
 *            function/function_model.cpp
 *
 *  Copyright  2011-20  Pieter Collins
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

#include "function/function_model.hpp"
#include "function/function_model.tpl.hpp"
#include "function/taylor_model.hpp"

#include "algebra/algebra.hpp"

#include "function/formula.hpp"

namespace Ariadne {

template class FunctionModelFactoryInterface<ValidatedTag,FloatDP,FloatDP>;
template class FunctionModelFactoryInterface<ValidatedTag,FloatMP,FloatMP>;

template class FunctionModel<ValidatedTag,RealScalar(RealScalar),FloatDP>;
template class FunctionModel<ValidatedTag,RealScalar(RealScalar),FloatMP>;
template class FunctionModel<ValidatedTag,RealVector(RealScalar),FloatDP>;
template class FunctionModel<ValidatedTag,RealVector(RealScalar),FloatMP>;

template class FunctionModel<ValidatedTag,RealScalar(RealVector),FloatDP>;
template class FunctionModel<ValidatedTag,RealScalar(RealVector),FloatMP>;
template class FunctionModel<ValidatedTag,RealVector(RealVector),FloatDP>;
template class FunctionModel<ValidatedTag,RealVector(RealVector),FloatMP>;

} // namespace Ariadne
