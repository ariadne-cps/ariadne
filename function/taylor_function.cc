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

template class FunctionPatch<ValidatedTaylorModel64>;
template class FunctionMixin<FunctionPatch<ValidatedTaylorModel64>,ApproximateTag,BoxDomain,IntervalDomain>;
template class FunctionMixin<FunctionPatch<ValidatedTaylorModel64>,ValidatedTag,BoxDomain,IntervalDomain>;
template class VectorFunctionPatch<ValidatedTaylorModel64>;
template class FunctionMixin<VectorFunctionPatch<ValidatedTaylorModel64>,ApproximateTag,BoxDomain,BoxDomain>;
template class FunctionMixin<VectorFunctionPatch<ValidatedTaylorModel64>,ValidatedTag,BoxDomain,BoxDomain>;

template class FunctionPatch<ValidatedTaylorModelMP>;
template class VectorFunctionPatch<ValidatedTaylorModelMP>;

CanonicalNumericType<ValidatedTag> TaylorFunctionFactory::_create(const Number<ValidatedTag>& number) const
{
    return CanonicalNumericType<ValidatedTag>(number,Precision64());
}

ScalarTaylorFunction TaylorFunctionFactory::create(const ExactBoxType& domain, const ValidatedScalarFunctionInterface& function) const
{
    return ScalarTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction* TaylorFunctionFactory::_create(const ExactBoxType& domain, const ValidatedScalarFunctionInterface& function) const
{
    return new ScalarTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction TaylorFunctionFactory::create(const ExactBoxType& domain, const ValidatedVectorFunctionInterface& function) const
{
    return VectorTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction* TaylorFunctionFactory::_create(const ExactBoxType& domain, const ValidatedVectorFunctionInterface& function) const
{
    return new VectorTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction TaylorFunctionFactory::create_zero(const ExactBoxType& domain) const
{
    return ScalarTaylorFunction::zero(domain,this->_sweeper);
}

ScalarTaylorFunction TaylorFunctionFactory::create_constant(const ExactBoxType& domain, ValidatedNumericType value) const
{
    return ScalarTaylorFunction::constant(domain,value,this->_sweeper);
}

ScalarTaylorFunction TaylorFunctionFactory::create_coordinate(const ExactBoxType& domain, SizeType k) const
{
    return ScalarTaylorFunction::coordinate(domain,k,this->_sweeper);
}


VectorTaylorFunction TaylorFunctionFactory::create_identity(const ExactBoxType& domain) const
{
    return VectorTaylorFunction::identity(domain,this->_sweeper);
}


FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper<Float64>());
}

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold) {
    return new TaylorFunctionFactory(ThresholdSweeper<Float64>(Precision64(),sweep_threshold));
}


} // namespace Ariadne
