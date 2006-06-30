/***************************************************************************
 *            newton.tpl
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
 
#include "../declarations.h"

#include "../linear_algebra/interval_vector.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"

#include "../system/vector_field.h"

namespace Ariadne {
  static int debug_level=0;
}


namespace Ariadne {
  namespace Evaluation {

    template<typename Real>
    Geometry::Rectangle<Real>
    interval_newton(const System::VectorField<Real>& f, 
                    const Geometry::Rectangle<Real>& x, 
                    const Real& e,
                    uint max_steps)
    {
      uint n=max_steps;
      if(debug_level>0) { std::cerr << "debug_level=" << debug_level << "\n"; }
      Geometry::Rectangle<Real> r=x;
      while(true) {
        if(debug_level>0) { std::cerr << "Testing for root in " << r << "\n"; }
        Geometry::Point<Real> m=r.centre();
        Geometry::Rectangle<Real> mr=m;
        LinearAlgebra::IntervalVector<Real> w=f(mr);
        LinearAlgebra::IntervalMatrix<Real> A=f.derivative(r);
        LinearAlgebra::IntervalMatrix<Real> Ainv=A.inverse();
        LinearAlgebra::IntervalVector<Real> dr=Ainv * w;
        Geometry::Rectangle<Real> nr= mr - dr;
        if(debug_level>0) {
          std::cerr << "e=" << r.radius() << "  x=" << r << "  m=" << m << std::flush;
          std::cerr << "  f(x)=" << f(r) << std::flush;
          std::cerr << "  f(m)=" << f(mr).centre() << std::flush;
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
        r=Geometry::intersection(nr,r);
        n=n-1;
      }
      throw EvaluationException("No result found -- disjoint");
    }
    
      

  }
}
