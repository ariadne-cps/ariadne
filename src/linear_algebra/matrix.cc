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
 
#include <cstdlib>
#include <cstdio>

#include "real_typedef.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix.tpl"

namespace Ariadne {
  namespace LinearAlgebra {
    
    template class Matrix<Real>;
    template class Matrix<Field>;

    template std::ostream& operator<<(std::ostream&, const Matrix<Real>&);
    template std::ostream& operator<<(std::ostream&, const Matrix<Field>&);
    
    template std::istream& operator>>(std::istream&, Matrix<Real>&);
    template std::istream& operator>>(std::istream&, Matrix<Field>&);
    
    template Matrix<Real> zero_Matrix(size_type r, size_type c);
    template Matrix<Field> zero_Matrix(size_type r, size_type c);

    template Matrix<Real> concatenate_columns(const Matrix<Real>& A1,
                                                const Matrix<Real>& A2);
    template Matrix<Field> concatenate_columns(const Matrix<Field>& A1,
                                                  const Matrix<Field>& A2);

    template Matrix<Real> exp_approx(const Matrix<Real> &A, 
                                       const Real& e); 
    
    template Matrix<Field> exp_approx(const Matrix<Field> &A, 
                                         const Field& e); 

    
    template void lu_local_dec(Matrix<Field>& A, 
                               const array<size_type>& row, 
                               const array<size_type>& col, 
                               const size_type& rows, 
                               const size_type& columns, 
                               const size_type& p);
    
    template Matrix<Field> lu_decompose(const Matrix<Field>& A, 
                                           array<size_type>& p_col, 
                                           array<size_type>& p_row);
                              
    template Matrix<Field> lu_decompose(const Matrix<Field> &A, 
                                           array<size_type>& p_array);

    template Vector<Field> lu_solve(const Matrix<Field>& A, 
                                       const array<size_type>& p_array, 
                                       const Vector<Field>& b);
                                       
    template Matrix<Field> Householder_QR(const Matrix<Field> &A);

    template Matrix<Real> hermitian(const Matrix<Real>& m);
    
    template Integer common_denominator(const Matrix<Real>& A);
    
    template Vector<Integer> row_common_denominators(const Matrix<Real>& A);
    template Vector<Integer> row_common_denominators(const Matrix<Field>& A);

    template void transform_linear_inequalities(const Matrix<Real>& T, 
                                                Matrix<Real>& A, 
                                                Vector<Real>& b);
    
    template bool independent_rows(Matrix<Real> A);

    template bool have_same_dimensions(const Matrix<Real> &A,  const Matrix<Real> &B);
    
    template bool equivalent_columns(const Matrix<Real> &A, 
                                     const size_type &A_col, 
                                     const Matrix<Real> &B, 
                                     const size_type &B_col);
    
  
    template size_type find_first_not_null_in_col(const Matrix<Real> &A, 
                                                  const size_type &col);
    
    template Matrix<Real> remove_null_columns_but_one(const Matrix<Real> &A);
    
    template void remove_null_columns(const Matrix<Real>& A, 
                                      array<size_type>& row, 
                                      array<size_type>& col);

    template Matrix<Field> compute_space(const Matrix<Field>& SA, 
                                            array<size_type>& row,
                                            const array<size_type>& col);
    
  }
}
