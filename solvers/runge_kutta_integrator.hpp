/***************************************************************************
 *            runge_kutta_integrator.hpp
 *
 *  Copyright  2010  Pieter Collins
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

#ifndef ARIADNE_RUNGE_KUTTA_INTEGRATOR_HPP
#define ARIADNE_RUNGE_KUTTA_INTEGRATOR_HPP

#include <iostream>

#include "utility/container.hpp"
#include "utility/declarations.hpp"
#include "solvers/integrator.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"

namespace Ariadne {

template<class T1, class T2> using Pair = std::pair<T1,T2>;
template<class T> class List;


class RungeKutta4Integrator
{
  public:
    RungeKutta4Integrator(double step_size);

    FloatApproximationVector
    step(const ApproximateVectorFunction& f, const FloatApproximationVector& x, const Float64Approximation& h) const;

    List< Pair<Float64Approximation,FloatApproximationVector> >
    evolve(const ApproximateVectorFunction& f, const FloatApproximationVector& x0, const Float64Approximation& tmax) const;
  private:
    double _step_size;
};

} // namespace Ariadne

#endif /* ARIADNE_RUNGE_KUTTA_INTEGRATOR_HPP */
