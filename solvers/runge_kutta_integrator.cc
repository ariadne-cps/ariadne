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
#include "function/function.h"

namespace Ariadne {

inline auto operator*(double s, ApproximateFloatVector v) -> decltype(ApproximateFloat(s)*v) { return ApproximateFloat(s)*v; }

RungeKutta4Integrator::RungeKutta4Integrator(double step_size)
    : _step_size(step_size)
{
}

ApproximateFloatVector
RungeKutta4Integrator::step(const ApproximateVectorFunctionInterface& f, const ApproximateFloatVector& x, const ApproximateFloat& h) const
{
    ApproximateFloatVector k1=f.evaluate(x);
    ApproximateFloatVector k2=f.evaluate(ApproximateFloatVector(x+(h/2)*k1));
    ApproximateFloatVector k3=f.evaluate(ApproximateFloatVector(x+(h/2)*k2));
    ApproximateFloatVector k4=f.evaluate(ApproximateFloatVector(x+h*k3));
    //std::cerr<<"k1,2,3,4="<<k1<<k2<<k3<<k4<<"\n";
    return x+(h/6)*(k1+2.0*k3+2.0*k4+k2);
}

List< Pair<ApproximateFloat,ApproximateFloatVector> >
RungeKutta4Integrator::evolve(const ApproximateVectorFunctionInterface& f, const ApproximateFloatVector& x0, const ApproximateFloat& tmax) const
{
    static const ApproximateFloat h=this->_step_size;

    List< Pair<ApproximateFloat,ApproximateFloatVector> > res;
    ApproximateFloat t=0.0;
    ApproximateFloatVector x=x0;

    res.push_back(make_pair(t,x));
    while(t<tmax) {
        x=this->step(f,x,h);
        t+=h;
        res.push_back(make_pair(t,x));
    }
    return res;
}

} // namespace Ariadne
