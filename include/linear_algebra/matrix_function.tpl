/***************************************************************************
 *            matrix_function.tpl
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 
#include "matrix_function.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    Matrix<R>
    exp_approx(const Matrix<R>& A, 
               const R& e) 
    {
      Matrix<R> result=Matrix<R>::identity(A.size1());
      
      R norm_A=norm(A);
      Matrix<R> term=result;
      uint n=0;
      while(norm(term)*n >= e*(n-norm_A)) {
        ++n;
        term=(term*A)/n;
        result=result+term;
      }
      return result;
    }

    template<typename R>
    IntervalMatrix<R>
    exp(const IntervalMatrix<R>& A) 
    {
      assert(A.size1()==A.size2());
      R err=A.radius_norm()/65536;
      if(err==0) {
        err=A.upper_norm()/65536;
        err/=65536;
        err/=65536;
      }
            
      IntervalMatrix<R> result=IntervalMatrix<R>::identity(A.size1())+A;
      IntervalMatrix<R> term=A;
      unsigned int n=1;
      while(term.upper_norm()>err) {
        n=n+1;
        term=(term*A)/Interval<R>(R(n));
        result+=term;
      }
      term=Interval<R>(-1,1)*term;
      result+=term;
      
      return result;
    }
   
  }
}
