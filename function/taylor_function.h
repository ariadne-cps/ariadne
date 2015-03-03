/***************************************************************************
 *            taylor_function.h
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

/*! \file taylor_function.h
 *  \brief Over-approximations of functions based on Taylor expansions.
 */

#ifndef ARIADNE_TAYLOR_FUNCTION_H
#define ARIADNE_TAYLOR_FUNCTION_H

#include <iosfwd>
#include "utility/container.h"
#include "utility/exceptions.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "function/taylor_model.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function_model.h"
#include "function/function_patch.h"

namespace Ariadne {

template<class P, class F> class TaylorModel;
class TaylorFunctionFactory;

typedef ScalarFunctionPatch<ValidatedTaylorModel> ScalarTaylorFunction;
typedef VectorFunctionPatch<ValidatedTaylorModel> VectorTaylorFunction;

/*
class ScalarTaylorFunction : public FunctionPatch<ValidatedTaylorModel> {
  public:
    using FunctionPatch<ValidatedTaylorModel>::FunctionPatch;
    ScalarTaylorFunction() : FunctionPatch<ValidatedTaylorModel>() { }
    ScalarTaylorFunction(FunctionPatch<ValidatedTaylorModel> const& f) : FunctionPatch<ValidatedTaylorModel>(f) { }
};

class VectorTaylorFunction : public VectorFunctionPatch<ValidatedTaylorModel> {
  public:
    using VectorFunctionPatch<ValidatedTaylorModel>::VectorFunctionPatch;
    VectorTaylorFunction() : VectorFunctionPatch<ValidatedTaylorModel>() { }
    VectorTaylorFunction(VectorFunctionPatch<ValidatedTaylorModel> const& f) : VectorFunctionPatch<ValidatedTaylorModel>(f) { }
};
*/

class TaylorFunctionFactory
    : public FunctionModelFactoryInterface<ValidatedTag>
{
    Sweeper _sweeper;
  public:
    explicit TaylorFunctionFactory(Sweeper sweeper) : _sweeper(sweeper) { }
    Sweeper sweeper() const { return this->_sweeper; }
    TaylorFunctionFactory* clone() const { return new TaylorFunctionFactory(this->_sweeper); }
    Void write(OutputStream& os) const { os << "TaylorFunctionFactory( sweeper=" << this->_sweeper << " )"; }
    ScalarTaylorFunction create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const;
    VectorTaylorFunction create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const;
    ScalarTaylorFunction create_zero(const ExactBox& domain) const;
    ScalarTaylorFunction create_constant(const ExactBox& domain, ValidatedNumericType c) const;
    ScalarTaylorFunction create_coordinate(const ExactBox& domain, SizeType k) const;
    VectorTaylorFunction create_zero(SizeType i, const ExactBox& domain) const;
    ScalarTaylorFunction create_identity(const ExactInterval& domain) const;
    VectorTaylorFunction create_identity(const ExactBox& domain) const;
  private:
    ScalarTaylorFunction* _create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const;
    VectorTaylorFunction* _create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const;
};



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
