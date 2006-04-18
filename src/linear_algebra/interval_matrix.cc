/***************************************************************************
 *            Matrix.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
#include "real_typedef.h"

#include "linear_algebra/interval_matrix.h"
#include "linear_algebra/interval_matrix.tpl"

namespace Ariadne {
  namespace LinearAlgebra {
    template class IntervalMatrix<Real>;
    template class IntervalMatrix<Field>;
      
    template std::ostream& operator<<(std::ostream&, const IntervalMatrix<Real>&);
    template std::ostream& operator<<(std::ostream&, const IntervalMatrix<Field>&);
           
    template IntervalVector<Real> prod(const Matrix<Real>&, const IntervalVector<Real>&);
    template IntervalVector<Field> prod(const Matrix<Field>&, const IntervalVector<Field>&);
    
    template IntervalVector<Real> prod(const IntervalMatrix<Real>&, const Vector<Real>&);
    template IntervalVector<Field> prod(const IntervalMatrix<Field>&, const Vector<Field>&);
    
    template IntervalVector<Real> prod(const IntervalMatrix<Real>&, const IntervalVector<Real>&);
    template IntervalVector<Field> prod(const IntervalMatrix<Field>&, const IntervalVector<Field>&);
    
    template IntervalMatrix<Real> prod(const Matrix<Real>&, const IntervalMatrix<Real>&);
    template IntervalMatrix<Field> prod(const Matrix<Field>&, const IntervalMatrix<Field>&);
    
    template IntervalMatrix<Real> prod(const IntervalMatrix<Real>&, const Matrix<Real>&);
    template IntervalMatrix<Field> prod(const IntervalMatrix<Field>&, const Matrix<Field>&);
    
    template IntervalMatrix<Real> prod(const IntervalMatrix<Real>&, const IntervalMatrix<Real>&);
    template IntervalMatrix<Field> prod(const IntervalMatrix<Field>&, const IntervalMatrix<Field>&);
    
    template IntervalMatrix<Real> fprod(const Matrix<Field>&, const IntervalMatrix<Real>&);
     
    template Matrix<Real> over_approximation(const IntervalMatrix<Real>&);
    template Matrix<Field> over_approximation(const IntervalMatrix<Field>&);

    template IntervalMatrix<Real> approximate(const Matrix<Field>& A,const Real& e); 

    template IntervalMatrix<Real> exp(const IntervalMatrix<Real>& A); 

  }
}
