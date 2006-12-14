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

#include "../debug.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"

#include "../system/vector_field.h"



template<class R>
Ariadne::Geometry::Rectangle<R>
Ariadne::Evaluation::IntervalNewtonSolver<R>::solve(const System::VectorField<R>& f, 
                                                    const Geometry::Rectangle<R>& x)
{
  const R& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  if(verbosity>1) { std::cerr << "verbosity=" << verbosity << "\n"; }
  Geometry::Rectangle<R> r=x;
  while(true) {
    if(verbosity>1) { std::cerr << "Testing for root in " << r << "\n"; }
    if(verbosity>1) { std::cerr << "  e=" << r.radius() << "  r=" << r << std::endl; }
    Geometry::Point<R> m=r.centre();
    if(verbosity>1) { std::cerr << "  m=" << m << std::endl; }
    Geometry::Rectangle<R> mr(m);
    if(verbosity>1) { std::cerr << "  mr=" << mr << std::endl; }
    LinearAlgebra::Vector< Interval<R> > w=f(mr);
    if(verbosity>1) { std::cerr << "  f(mr)=" << w << std::endl; }
    LinearAlgebra::Matrix< Interval<R> > A=f.jacobian(r);
    if(verbosity>1) { std::cerr << "  Df(r)=" << A << std::endl; }
    LinearAlgebra::Matrix< Interval<R> > Ainv=A.inverse();
    if(verbosity>1) { std::cerr << "  inverse(Df(r))=" << Ainv << std::endl; }
    LinearAlgebra::Vector< Interval<R> > dr=Ainv * w;
    if(verbosity>1) { std::cerr << "  dr=" << dr << std::endl; }
    Geometry::Rectangle<R> nr= mr - dr;
    if(verbosity>1) { std::cerr << "  nr=" << nr << std::endl; } 
    if(verbosity>1) {
      std::cerr << "  f(x)=" << f(r) << std::flush;
      std::cerr << "  f(m)=" << approximate_value(f(mr)) << std::flush;
      std::cerr << "  Df(x) =" << A << "  inv=" << inverse(A) << "  I=" << A*inverse(A) << std::flush;
      std::cerr << "  nx =" << nr << "\n\n" << std::flush;
      std::cerr << nr << " subset " << r << " ? " << Geometry::subset(nr,r) << "\n";
      std::cerr << nr.radius() << " < " << e << " ? " << (nr.radius() < e) << "\n";
    }
    if(Geometry::subset(nr,r) and (nr.radius() < e)) {   \
      return nr;
    }
    if(Geometry::disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=Geometry::closed_intersection(nr,r);
    n=n-1;
  }
  throw EvaluationException("No result found -- disjoint");
}

template<class R>
Ariadne::Geometry::Rectangle<R>
Ariadne::Evaluation::interval_newton(const System::VectorField<R>& f, 
                                     const Geometry::Rectangle<R>& x,
                                     const R& e, 
                                     uint s)
{
  return IntervalNewtonSolver<R>(e,s).solve(f,x);
}
