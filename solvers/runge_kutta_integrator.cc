/***************************************************************************
 *            runge_kutta_integrator.cc
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

#include "utility/standard.h"
#include "config.h"

#include <iostream>

#include "solvers/runge_kutta_integrator.h"

#include "utility/container.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/algebra.h"
#include "function/function.h"
#include "function/taylor_model.h"
#include "expression/formula.h"

namespace Ariadne {

inline auto operator*(double s, FloatApproximationVector v) -> decltype(Float64Approximation(s)*v) { return Float64Approximation(s)*v; }

RungeKutta4Integrator::RungeKutta4Integrator(double step_size)
    : _step_size(step_size)
{
}

FloatApproximationVector
RungeKutta4Integrator::step(const ApproximateVectorFunction& f, const FloatApproximationVector& x, const Float64Approximation& h) const
{
    FloatApproximationVector k1=f(x);
    FloatApproximationVector k2=f(FloatApproximationVector(x+(h/2)*k1));
    FloatApproximationVector k3=f(FloatApproximationVector(x+(h/2)*k2));
    FloatApproximationVector k4=f(FloatApproximationVector(x+h*k3));
    //std::cerr<<"k1,2,3,4="<<k1<<k2<<k3<<k4<<"\n";
    return x+(h/6)*(k1+2.0*k3+2.0*k4+k2);
}

List< Pair<Float64Approximation,FloatApproximationVector> >
RungeKutta4Integrator::evolve(const ApproximateVectorFunction& f, const FloatApproximationVector& x0, const Float64Approximation& tmax) const
{
    static const Float64Approximation h(this->_step_size,Precision64());

    List< Pair<Float64Approximation,FloatApproximationVector> > res;
    Float64Approximation t(0.0,Precision64());
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
