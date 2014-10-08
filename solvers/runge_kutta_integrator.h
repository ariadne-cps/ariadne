/***************************************************************************
 *            runge_kutta_integrator.h
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

#ifndef ARIADNE_RUNGE_KUTTA_INTEGRATOR_H
#define ARIADNE_RUNGE_KUTTA_INTEGRATOR_H

#include <iostream>

#include "solvers/integrator.h"
#include "utility/container.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "function/function.h"

namespace Ariadne {

template<class T1, class T2> class Pair;
template<class T> class List;

class ApproximateFloat;
template<class X> class Vector;
typedef Vector<ApproximateFloat> ApproximateFloatVector;
template<class P> class VectorFunctionInterface;
typedef VectorFunctionInterface<ApproximateTag> ApproximateVectorFunctionInterface;

class RungeKutta4Integrator
{
  public:
    RungeKutta4Integrator(double step_size);

    ApproximateFloatVector
    step(const ApproximateVectorFunctionInterface& f, const ApproximateFloatVector& x, const ApproximateFloat& h) const;

    List< Pair<ApproximateFloat,ApproximateFloatVector> >
    evolve(const ApproximateVectorFunctionInterface& f, const ApproximateFloatVector& x0, const ApproximateFloat& tmax) const;
  private:
    double _step_size;
};

} // namespace Ariadne

#endif /* ARIADNE_RUNGE_KUTTA_INTEGRATOR_H */
