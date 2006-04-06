/***************************************************************************
 *            matrix.cc
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
    
    template class matrix<Real>;
    template class matrix<Field>;

    template std::ostream& operator<<(std::ostream&, const matrix<Real>&);
    template std::ostream& operator<<(std::ostream&, const matrix<Field>&);
    
    template matrix<Real> zero_matrix(size_type r, size_type c);
    template matrix<Field> zero_matrix(size_type r, size_type c);

    template Real norm(const matrix<Real>& A);
    template Field norm(const matrix<Field>& A);

    template Real log_norm(const matrix<Real>& A);
    template Field log_norm(const matrix<Field>& A);

    template matrix<Real> concatenate_columns(const matrix<Real>& A1,
                                                const matrix<Real>& A2);
    template matrix<Field> concatenate_columns(const matrix<Field>& A1,
                                                  const matrix<Field>& A2);

    template matrix<Real> exp_approx(const matrix<Real> &A, 
                                       const Real& e); 
    
    template matrix<Field> exp_approx(const matrix<Field> &A, 
                                         const Field& e); 

    
    template void lu_local_dec(matrix<Field>& A, 
                               const array<size_type>& row, 
                               const array<size_type>& col, 
                               const size_type& rows, 
                               const size_type& columns, 
                               const size_type& p);
    
    template matrix<Field> lu_decompose(const matrix<Field>& A, 
                                           array<size_type>& p_col, 
                                           array<size_type>& p_row);
                              
    template matrix<Field> lu_decompose(const matrix<Field> &A, 
                                           array<size_type>& p_array);

    template vector<Field> lu_solve(const matrix<Field>& A, 
                                       const array<size_type>& p_array, 
                                       const vector<Field>& b);
                                       
    template matrix<Field> Householder_QR(const matrix<Field> &A);

    template matrix<Real> hermitian(const matrix<Real>& m);
    
    template matrix<Field> inverse(const matrix<Real> &A);
    template matrix<Field> inverse(const matrix<Field> &A);
    template matrix<Float64> inverse(const matrix<Float64> &A);
    
    template Integer common_denominator(const matrix<Real>& A);
    
    template vector<Integer> row_common_denominators(const matrix<Real>& A);
    template vector<Integer> row_common_denominators(const matrix<Field>& A);

    template void transform_linear_inequalities(const matrix<Real>& T, 
                                                matrix<Real>& A, 
                                                vector<Real>& b);
    
    template bool independent_rows(matrix<Real> A);

    template bool have_same_dimensions(const matrix<Real> &A,  const matrix<Real> &B);
    
    template bool equivalent_columns(const matrix<Real> &A, 
                                     const size_type &A_col, 
                                     const matrix<Real> &B, 
                                     const size_type &B_col);
    
  
    template size_type find_first_not_null_in_col(const matrix<Real> &A, 
                                                  const size_type &col);
    
    template matrix<Real> remove_null_columns_but_one(const matrix<Real> &A);
    
    template void remove_null_columns(const matrix<Real>& A, 
                                      array<size_type>& row, 
                                      array<size_type>& col);

    template matrix<Field> compute_space(const matrix<Field>& SA, 
                                            array<size_type>& row,
                                            const array<size_type>& col);
    
  }
}
