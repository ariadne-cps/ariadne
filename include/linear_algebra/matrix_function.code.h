/***************************************************************************
 *            matrix_function.code.h
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

#include "../numeric/interval.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    template<class R>
    Matrix< Numeric::Interval<R> >
    exp(const Matrix< Numeric::Interval<R> >& A) 
    {
      using Numeric::Interval;
      
      ARIADNE_CHECK_SQUARE(A,__PRETTY_FUNCTION__);
      R err=div_up(A.norm().upper(),R(65536));
      if(err==0) {
        err=div_up(A.norm().upper(),R(65536));
        err=div_up(err,R(65536));
        err=div_up(err,R(65536));
      }
            
      Matrix< Interval<R> > result=Matrix< Interval<R> >::identity(A.number_of_rows())+A;
      Matrix< Interval<R> > term=A;
      unsigned int n=1;
      while(term.norm().upper()>err) {
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
