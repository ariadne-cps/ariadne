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
 
#include "newton.h"

#include "output/logging.h"
#include "geometry/point.h"

#include "system/vector_field.h"

#include "evaluation/exceptions.h"

namespace Ariadne {

namespace Evaluation { static int& verbosity = hybrid_evolver_verbosity; }

template<class R>
Geometry::Point<typename Evaluation::IntervalNewtonSolver<R>::I>
Evaluation::IntervalNewtonSolver<R>::solve(const Function::FunctionInterface<R>& f, 
                                           const Geometry::Point<I>& ix)
{
  const R& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
  Geometry::Point<I> x=ix;
  Geometry::Box<R> r(x);
  while(n>0) {
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<Geometry::radius(x)<<"  x="<<x<<"\n");
    Geometry::Point<R> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Geometry::Point<I> im(m);
    LinearAlgebra::Vector<I> w=f(im.position_vector());
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    LinearAlgebra::Matrix<I> A=f.jacobian(x.position_vector());
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    LinearAlgebra::Matrix<I> Ainv=LinearAlgebra::inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    LinearAlgebra::Vector<I> dx=Ainv * w;
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Geometry::Point<I> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Geometry::Box<R> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");

    if(Geometry::subset(nr,r) && Geometry::radius(nx) < e) {
      return nr;
    }
    if(Geometry::disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=Geometry::closed_intersection(nr,r);
    x=r;
    n=n-1;
  }
  throw EvaluationException("No result found -- maximum number of steps reached");
}

}
