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
 
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix.tpl"

namespace boost {
  namespace numeric {
    namespace ublas {
      using namespace Ariadne;
      
      template class matrix<Dyadic>;
      template class matrix<Rational>;
      template class matrix<Float64>;

      template std::ostream& operator<<(std::ostream&, const matrix<Dyadic>&);
      template std::ostream& operator<<(std::ostream&, const matrix<Rational>&);
      template std::ostream& operator<<(std::ostream&, const matrix<Float64>&);
    }
  }
}

namespace Ariadne {
  namespace LinearAlgebra {
    
    template matrix<Dyadic> zero_matrix(size_type r, size_type c);
    template matrix<Rational> zero_matrix(size_type r, size_type c);

    template Dyadic norm(const matrix<Dyadic>& A);
    template Rational norm(const matrix<Rational>& A);

    template Dyadic log_norm(const matrix<Dyadic>& A);
    template Rational log_norm(const matrix<Rational>& A);

    template matrix<Dyadic> concatenate_columns(const matrix<Dyadic>& A1,
                                                const matrix<Dyadic>& A2);
    template matrix<Rational> concatenate_columns(const matrix<Rational>& A1,
                                                  const matrix<Rational>& A2);

    template matrix<Dyadic> exp_approx(const matrix<Dyadic> &A, 
                                       const Dyadic& e); 
    
    template matrix<Rational> exp_approx(const matrix<Rational> &A, 
                                         const Rational& e); 

    
    template void lu_local_dec(matrix<Dyadic>& A, 
                               const array<size_type>& row, 
                               const array<size_type>& col, 
                               const size_type& rows, 
                               const size_type& columns, 
                               const size_type& p);
    template void lu_local_dec(matrix<Rational>& A, 
                               const array<size_type>& row, 
                               const array<size_type>& col, 
                               const size_type& rows, 
                               const size_type& columns, 
                               const size_type& p);
    
    template matrix<Dyadic> lu_decompose(const matrix<Dyadic>& A, 
                                         array<size_type>& p_col, 
                                         array<size_type>& p_row);
    template matrix<Rational> lu_decompose(const matrix<Rational>& A, 
                                           array<size_type>& p_col, 
                                           array<size_type>& p_row);
                              
    template matrix<Dyadic> lu_decompose(const matrix<Dyadic> &A, 
                                         array<size_type>& p_array);
    template matrix<Rational> lu_decompose(const matrix<Rational> &A, 
                                           array<size_type>& p_array);

    template vector<Dyadic> lu_solve(const matrix<Dyadic>& A, 
                                     const array<size_type>& p_array, 
                                     const vector<Dyadic>& b);
             
    template vector<Rational> lu_solve(const matrix<Rational>& A, 
                                       const array<size_type>& p_array, 
                                       const vector<Rational>& b);
                                       
    template matrix<Dyadic> Householder_QR(const matrix<Dyadic> &A);

    template matrix<Dyadic> hermitian(const matrix<Dyadic>& m);
    
    template matrix<Dyadic> inverse(const matrix<Dyadic> &A);
    template matrix<Rational> inverse(const matrix<Rational> &A);
    template matrix<Float64> inverse(const matrix<Float64> &A);
    
    template Integer common_denominator(const matrix<Dyadic>& A);
    
    template vector<Integer> row_common_denominators(const matrix<Dyadic>& A);
    template vector<Integer> row_common_denominators(const matrix<Rational>& A);

    template void transform_linear_inequalities(const matrix<Dyadic>& T, 
                                                matrix<Dyadic>& A, 
                                                vector<Dyadic>& b);
    
    template bool independent_rows(matrix<Dyadic> A);

    template bool have_same_dimensions(const matrix<Dyadic> &A,  const matrix<Dyadic> &B);
    
    template bool equivalent_columns(const matrix<Dyadic> &A, 
                                     const size_type &A_col, 
                                     const matrix<Dyadic> &B, 
                                     const size_type &B_col);
    
  
    template size_type find_first_not_null_in_col(const matrix<Dyadic> &A, 
                                                  const size_type &col);
    
    template matrix<Dyadic> remove_null_columns_but_one(const matrix<Dyadic> &A);
    
    template void remove_null_columns(const matrix<Dyadic>& A, 
                                      array<size_type>& row, 
                                      array<size_type>& col);

    template matrix<Dyadic> compute_space(const matrix<Dyadic>& SA, 
                                          array<size_type>& row,
                                          const array<size_type>& col);
    
  }
}
