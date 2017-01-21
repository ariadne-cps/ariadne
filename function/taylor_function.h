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

typedef ScalarFunctionPatch<ValidatedTaylorModel64> ScalarTaylorFunction;
typedef VectorFunctionPatch<ValidatedTaylorModel64> VectorTaylorFunction;

/*
class ScalarTaylorFunction : public FunctionPatch<ValidatedTaylorModel64> {
  public:
    using FunctionPatch<ValidatedTaylorModel64>::FunctionPatch;
    ScalarTaylorFunction() : FunctionPatch<ValidatedTaylorModel64>() { }
    ScalarTaylorFunction(FunctionPatch<ValidatedTaylorModel64> const& f) : FunctionPatch<ValidatedTaylorModel64>(f) { }
};

class VectorTaylorFunction : public VectorFunctionPatch<Vali_datedTaylorModel64> {
  public:
    using VectorFunctionPatch<ValidatedTaylorModel64>::VectorFunctionPatch;
    VectorTaylorFunction() : VectorFunctionPatch<ValidatedTaylorModel64>() { }
    VectorTaylorFunction(VectorFunctionPatch<ValidatedTaylorModel64> const& f) : VectorFunctionPatch<ValidatedTaylorModel64>(f) { }
};
*/

class TaylorFunctionFactory
    : public FunctionPatchFactory<TaylorModel<ValidatedTag,Float64>>
{
    typedef TaylorModel<ValidatedTag,Float64> M;
  public:
    typedef typename M::SweeperType SweeperType;
    using FunctionPatchFactory<M>::FunctionPatchFactory;
    explicit TaylorFunctionFactory(SweeperType sweeper) : FunctionPatchFactory<M>(sweeper) { }
    SweeperType sweeper() const { return this->properties(); }
    friend OutputStream& operator<<(OutputStream& os, TaylorFunctionFactory const& factory) { return os << "TaylorFunctionFactory( sweeper=" << factory.sweeper() << " )"; }
};



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
