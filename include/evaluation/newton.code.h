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

#include "../logging.h"
#include "../geometry/point.h"

#include "../system/vector_field.h"

#include "../evaluation/exceptions.h"


template<class R>
//Ariadne::Geometry::Point<typename IntervalNewtonSolver<R>::I>
Ariadne::Geometry::Point<typename Ariadne::Evaluation::IntervalNewtonSolver<R>::I>
Ariadne::Evaluation::IntervalNewtonSolver<R>::solve(const System::VectorField<R>& f, 
                                                    const Geometry::Point<I>& ix)
{
  const R& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  if(verbosity>1) { std::cerr << "verbosity=" << verbosity << "\n"; }
  Geometry::Point<I> x=ix;
  Geometry::Rectangle<R> r(x);
  while(n>0) {
    if(verbosity>1) { std::cerr << "Testing for root in " << x << "\n"; }
    if(verbosity>1) { std::cerr << "  e=" << Geometry::error_bound(x) << "  x=" << x << std::endl; }
    Geometry::Point<R> m=approximate_value(x);
    if(verbosity>1) { std::cerr << "  m=" << m << std::endl; }
    Geometry::Point<I> im(m);
    LinearAlgebra::Vector<I> w=f(im);
    if(verbosity>1) { std::cerr << "  f(m)=" << w << std::endl; }
    LinearAlgebra::Matrix<I> A=f.jacobian(x);
    if(verbosity>1) { std::cerr << "  Df(r)=" << A << std::endl; }
    LinearAlgebra::Matrix<I> Ainv=A.inverse();
    if(verbosity>1) { std::cerr << "  inverse(Df(r))=" << Ainv << std::endl; }
    LinearAlgebra::Vector<I> dx=Ainv * w;
    if(verbosity>1) { std::cerr << "  dx=" << dx << std::endl; }
    Geometry::Point<I> nx= m - dx;
    if(verbosity>1) { std::cerr << "  nx=" << nx << std::endl; } 
    Geometry::Rectangle<R> nr(nx);
    if(verbosity>1) { std::cerr << "  nr=" << nr << std::endl; } 

    if(verbosity>1) {
      std::cerr << "  f(x)=" << f(x) << std::flush;
      std::cerr << "  f(m)=" << approximate_value(f(im)) << std::flush;
      std::cerr << "  Df(x) =" << A << "  inv=" << inverse(A) << "  I=" << A*inverse(A) << std::flush;
      std::cerr << "  nx =" << nx << "\n" << std::flush;
      std::cerr << "  nr =" << nr << "\n" << std::flush;
      std::cerr << "\n";
      std::cerr << error_bound(nx) << " < " << e << " ? " << (error_bound(nx) < e) << "\n\n";
    }
    
    if(Geometry::subset(nr,r) && Geometry::error_bound(nx) < e) {
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
