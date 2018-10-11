/***************************************************************************
 *            runge_kutta_integrator.cpp
 *
 *  Copyright  2010  Pieter Collins
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

#include "../utility/standard.hpp"
#include "../config.hpp"

#include <iostream>

#include "../solvers/runge_kutta_integrator.hpp"

#include "../utility/container.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/algebra.hpp"
#include "../function/function.hpp"
#include "../function/taylor_model.hpp"
#include "../function/formula.hpp"

namespace Ariadne {

inline auto operator*(double s, FloatApproximationVector v) -> decltype(FloatDPApproximation(s)*v) { return FloatDPApproximation(s)*v; }

RungeKutta4Integrator::RungeKutta4Integrator(double step_size)
    : _step_size(step_size)
{
}

FloatApproximationVector
RungeKutta4Integrator::step(const ApproximateVectorMultivariateFunction& f, const FloatApproximationVector& x, const FloatDPApproximation& h) const
{
    FloatApproximationVector k1=f(x);
    FloatApproximationVector k2=f(FloatApproximationVector(x+(h/2)*k1));
    FloatApproximationVector k3=f(FloatApproximationVector(x+(h/2)*k2));
    FloatApproximationVector k4=f(FloatApproximationVector(x+h*k3));
    //std::cerr<<"k1,2,3,4="<<k1<<k2<<k3<<k4<<"\n";
    return x+(h/6)*(k1+2.0*k3+2.0*k4+k2);
}

List< Pair<FloatDPApproximation,FloatApproximationVector> >
RungeKutta4Integrator::evolve(const ApproximateVectorMultivariateFunction& f, const FloatApproximationVector& x0, const FloatDPApproximation& tmax) const
{
    static const FloatDPApproximation h(this->_step_size,dp);

    List< Pair<FloatDPApproximation,FloatApproximationVector> > res;
    FloatDPApproximation t(0.0,dp);
    FloatApproximationVector x=x0;

    res.push_back(make_pair(t,x));
    while(decide(t<tmax)) {
        x=this->step(f,x,h);
        t+=h;
        res.push_back(make_pair(t,x));
    }
    return res;
}

} // namespace Ariadne
