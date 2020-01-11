/***************************************************************************
 *            first_order_pde.hpp
 *
 *  Copyright  2018-20  Pieter Collins, Svetlana Selivanova
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

/*! \file first_order_pde.hpp
 *  \brief Rigorous solutions of first-order hyperbolic partial differential equations.
 */

#ifndef ARIADNE_FIRST_ORDER_PDE_HPP
#define ARIADNE_FIRST_ORDER_PDE_HPP

#include "numeric/float.decl.hpp"
#include "algebra/tensor.hpp"
#include "function/function.decl.hpp"

namespace Ariadne {

//! \brief Compute the value of the multiaffine interpolation of the data \a us at the point \a x.
template<class X> Vector<Bounds<X>> multiaffine_interpolate(Tensor<2,Vector<Bounds<X>>> const& us, Vector<Value<X>> const& x);

//! \brief Solve the first-order paritial differential equation
//! \f[ A u_{,t} + \sum_{j=1}^{n} B_j u_{,x_j} = f(u); \quad u(0,x) = \phi_0(x); \quad u:\R^n \to \R^m . \f]
template<class PR> Tuple< FloatValue<PR>, FloatValue<PR>, Tensor<3,Vector<FloatBounds<PR>>>, FloatUpperBound<PR> >
first_order_pde(Matrix<Real> rA, Array<Matrix<Real>> rBs, Array<DiagonalMatrix<Real>> const& rDs, Array<Matrix<Real>> const& rTs, EffectiveVectorMultivariateFunction const& f, EffectiveVectorMultivariateFunction const& phi0, PR pr);

} // namespace Ariadne

#endif // ARIADNE_FIRST_ORDER_PDE_HPP
