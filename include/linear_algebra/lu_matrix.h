/***************************************************************************
 *            matrix.h
 *
 *  Mon May  3 12:31:15 2004
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
 
/*! \file LUmatrix.h
 *  \brief LU factorisation and factorised matrices.
 */

#ifndef _ARIADNE_LU_MATRIX_H
#define _ARIADNE_LU_MATRIX_H

#include "../linear_algebra/linear_algebra_declarations.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix stored in LU product form. */
    template<typename R>
    class LUMatrix {
      Matrix<R> L() const;
      Matrix<R> U() const;
      
      operator Matrix<R> () const;
      Matrix<R> inverse() const;
     private:
      Matrix<R> _elements;
      array<size_type> _row_permuation;
      array<size_type> _column_permuation;
    };
    
    template <typename R>
    void 
    lu_local_dec(Matrix<R>& A, 
                 const array<size_type>& row, 
                 const array<size_type>& col, 
                 const size_type& rows, 
                 const size_type& columns, 
                 const size_type& p);
   
    template <typename R>
    Matrix<R> 
    lu_decompose(const Matrix<R>& A, 
                 array<size_type>& p_col, 
                 array<size_type>& p_row);
                              
    /* PAY ATTENTION!!! 
     * I supose that matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename R>
    Matrix<R> 
    lu_decompose(const Matrix<R> &A, 
                 array<size_type>& p_array);

    
    /* PAY ATTENTION!!! 
     * I supose that boost::numeric::ublas::matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename R>
    Vector<R> 
    lu_solve(const Matrix<R>& A, 
             const array<size_type>& p_array, 
             const Vector<R>& b);
             

  }
}


#endif /* _ARIADNE_LU_MATRIX_H */
