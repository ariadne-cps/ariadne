/***************************************************************************
 *            taylor_function.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file taylor_function.hpp
 *  \brief Over-approximations of functions based on Taylor expansions.
 */

#ifndef ARIADNE_TAYLOR_FUNCTION_HPP
#define ARIADNE_TAYLOR_FUNCTION_HPP

#include <iosfwd>
#include "utility/container.hpp"
#include "utility/exceptions.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/taylor_model.hpp"

#include "function/function_interface.hpp"
#include "function/function_mixin.hpp"
#include "function/function_model.hpp"
#include "function/scaled_function_patch.hpp"

namespace Ariadne {

template<class P, class F> class TaylorModel;
class TaylorFunctionFactory;

typedef ScalarScaledFunctionPatch<ValidatedTaylorModel64> ValidatedScalarTaylorFunctionModel64;
typedef VectorScaledFunctionPatch<ValidatedTaylorModel64> ValidatedVectorTaylorFunctionModel64;

/*
class ValidatedScalarTaylorFunctionModel64 : public ScaledFunctionPatch<ValidatedTaylorModel64> {
  public:
    using ScaledFunctionPatch<ValidatedTaylorModel64>::ScaledFunctionPatch;
    ValidatedScalarTaylorFunctionModel64() : ScaledFunctionPatch<ValidatedTaylorModel64>() { }
    ValidatedScalarTaylorFunctionModel64(ScaledFunctionPatch<ValidatedTaylorModel64> const& f) : ScaledFunctionPatch<ValidatedTaylorModel64>(f) { }
};

class ValidatedVectorTaylorFunctionModel64 : public VectorScaledFunctionPatch<Vali_datedTaylorModel64> {
  public:
    using VectorScaledFunctionPatch<ValidatedTaylorModel64>::VectorScaledFunctionPatch;
    ValidatedVectorTaylorFunctionModel64() : VectorScaledFunctionPatch<ValidatedTaylorModel64>() { }
    ValidatedVectorTaylorFunctionModel64(VectorScaledFunctionPatch<ValidatedTaylorModel64> const& f) : VectorScaledFunctionPatch<ValidatedTaylorModel64>(f) { }
};
*/

class TaylorFunctionFactory
    : public ScaledFunctionPatchFactory<TaylorModel<ValidatedTag,Float64>>
{
    typedef TaylorModel<ValidatedTag,Float64> M;
  public:
    typedef typename M::SweeperType SweeperType;
    using ScaledFunctionPatchFactory<M>::ScaledFunctionPatchFactory;
    explicit TaylorFunctionFactory(SweeperType sweeper) : ScaledFunctionPatchFactory<M>(sweeper) { }
    SweeperType sweeper() const { return this->properties(); }
    friend OutputStream& operator<<(OutputStream& os, TaylorFunctionFactory const& factory) { return os << "TaylorFunctionFactory( sweeper=" << factory.sweeper() << " )"; }
};



} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_HPP
