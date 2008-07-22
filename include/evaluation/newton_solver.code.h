/***************************************************************************
 *            newton.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "newton_solver.h"

#include "output/logging.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/function_interface.h"
#include "geometry/point.h"
#include "geometry/box.h"

#include "evaluation/exceptions.h"

namespace Ariadne {


template<class R>
Vector<typename IntervalNewtonSolver<R>::I>
IntervalNewtonSolver<R>::solve(const FunctionInterface<R>& f, 
                               const Vector<I>& ix)
{
  const R& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
  Point<I> x(ix);
  Box<R> r(x);
  while(n>0) {
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Point<R> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Point<I> im(m);
    Vector<I> w=f(im.position_vector());
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    Matrix<I> A=f.jacobian(x.position_vector());
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    Matrix<I> Ainv=inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    Vector<I> dx=Ainv * w;
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Point<I> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Box<R> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");

    if(subset(nr,r) && radius(nx) < e) {
      return nr.position_vectors();
    }
    if(disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=closed_intersection(nr,r);
    x=r;
    n=n-1;
  }
  throw EvaluationException("No result found -- maximum number of steps reached");
}

} // namespace Ariadne
