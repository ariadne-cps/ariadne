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

#include "standard.h"
#include "config.h"

#include <iostream>

#include "runge_kutta_integrator.h"

#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "function.h"

namespace Ariadne {

RungeKutta4Integrator::RungeKutta4Integrator(double step_size)
    : _step_size(step_size)
{
}

FloatVector
RungeKutta4Integrator::step(const ApproximateVectorFunctionInterface& f, const FloatVector& x, const Float& h) const
{
    FloatVector k1=f.evaluate(x);
    FloatVector k2=f.evaluate(FloatVector(x+(h/2)*k1));
    FloatVector k3=f.evaluate(FloatVector(x+(h/2)*k2));
    FloatVector k4=f.evaluate(FloatVector(x+h*k3));
    //std::cerr<<"k1,2,3,4="<<k1<<k2<<k3<<k4<<"\n";
    return x+(h/6)*(k1+2.0*k3+2.0*k4+k2);
}

List< Pair<Float,FloatVector> >
RungeKutta4Integrator::evolve(const ApproximateVectorFunctionInterface& f, const FloatVector& x0, const Float& tmax) const
{
    static const Float h=this->_step_size;

    List< Pair<Float,FloatVector> > res;
    Float t=0.0;
    FloatVector x=x0;

    res.push_back(make_pair(t,x));
    while(t<tmax) {
        x=this->step(f,x,h);
        t+=h;
        res.push_back(make_pair(t,x));
    }
    return res;
}

} // namespace Ariadne
