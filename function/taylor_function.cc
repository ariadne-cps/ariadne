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
#include "algebra/multi_index.h"
#include "function/polynomial.h"
#include "algebra/differential.h"
#include "algebra/evaluate.h"
#include "function/taylor_model.h"

#include "function/function.h"
#include "function/function_patch.h"
#include "function/taylor_function.h"

#include "function_mixin.tcc"
#include "function_patch.tcc"

#define VOLATILE ;

namespace Ariadne {

static double TAYLOR_FUNCTION_WRITING_ACCURACY = 1e-8;

template class FunctionPatch<ValidatedTaylorModel>;
template class VectorFunctionPatch<ValidatedTaylorModel>;

ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ScalarTaylorFunction& f) {
    this->_ptr=clone_on_copy_ptr< ScalarFunctionModelInterface<ValidatedTag> >(new ScalarTaylorFunction(f)); return *this;
}



ScalarTaylorFunction TaylorFunctionFactory::create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const
{
    return ScalarTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction* TaylorFunctionFactory::_create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const
{
    return new ScalarTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction TaylorFunctionFactory::create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const
{
    return VectorTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction* TaylorFunctionFactory::_create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const
{
    return new VectorTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction TaylorFunctionFactory::create_zero(const ExactBox& domain) const
{
    return ScalarTaylorFunction::zero(domain,this->_sweeper);
}

ScalarTaylorFunction TaylorFunctionFactory::create_constant(const ExactBox& domain, ValidatedNumericType value) const
{
    return ScalarTaylorFunction::constant(domain,value,this->_sweeper);
}

ScalarTaylorFunction TaylorFunctionFactory::create_coordinate(const ExactBox& domain, SizeType k) const
{
    return ScalarTaylorFunction::coordinate(domain,k,this->_sweeper);
}


VectorTaylorFunction TaylorFunctionFactory::create_identity(const ExactBox& domain) const
{
    return VectorTaylorFunction::identity(domain,this->_sweeper);
}


FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper());
}

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold) {
    return new TaylorFunctionFactory(ThresholdSweeper(sweep_threshold));
}


} // namespace Ariadne
