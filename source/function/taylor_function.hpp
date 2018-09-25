/***************************************************************************
 *            taylor_function.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file taylor_function.hpp
 *  \brief Over-approximations of functions based on Taylor expansions.
 */

#ifndef ARIADNE_TAYLOR_FUNCTION_HPP
#define ARIADNE_TAYLOR_FUNCTION_HPP

#include <iosfwd>
#include "../utility/container.hpp"
#include "../utility/exceptions.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../function/taylor_model.hpp"

#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function_model.hpp"
#include "../function/scaled_function_patch.hpp"

namespace Ariadne {

template<class P, class F> class TaylorModel;
class TaylorFunctionFactory;

template<class P, class F> using ScalarMultivariateTaylorFunctionModel = ScalarScaledFunctionPatch<TaylorModel<P,F>>;
template<class P, class F> using VectorMultivariateTaylorFunctionModel = VectorScaledFunctionPatch<TaylorModel<P,F>>;

template<class F> using ValidatedScalarMultivariateTaylorFunctionModel = ScalarScaledFunctionPatch<ValidatedTaylorModel<F>>;
template<class F> using ValidatedVectorMultivariateTaylorFunctionModel = VectorScaledFunctionPatch<ValidatedTaylorModel<F>>;

using ValidatedScalarMultivariateTaylorFunctionModelDP = ScalarScaledFunctionPatch<ValidatedTaylorModelDP>;
using ValidatedVectorMultivariateTaylorFunctionModelDP = VectorScaledFunctionPatch<ValidatedTaylorModelDP>;


class TaylorFunctionFactory
    : public ScaledFunctionPatchFactory<TaylorModel<ValidatedTag,FloatDP>>
{
    typedef TaylorModel<ValidatedTag,FloatDP> M;
    typedef BoxDomainType D;
  public:
    typedef typename M::SweeperType SweeperType;
    using ScaledFunctionPatchFactory<M>::ScaledFunctionPatchFactory;
    explicit TaylorFunctionFactory(SweeperType sweeper) : ScaledFunctionPatchFactory<M>(sweeper) { }
    SweeperType sweeper() const { return this->properties(); }
    friend OutputStream& operator<<(OutputStream& os, TaylorFunctionFactory const& factory) { return os << "TaylorFunctionFactory( sweeper=" << factory.sweeper() << " )"; }
};

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory();
FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold_);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_HPP
