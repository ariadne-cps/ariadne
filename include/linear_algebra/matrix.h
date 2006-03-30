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
 
/*! \file matrix.h
 *  \brief Matrices and matrix operations.
 */

#ifndef _ARIADNE_MATRIX_H
#define _ARIADNE_MATRIX_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include "../base/basic_type.h"
#include "../base/numerical_type.h"
#include "../base/interval.h"
#include "../base/utility.h"

namespace boost {
  namespace numeric {
    namespace ublas {
            
      template<typename Real>
      inline
      vector<Real>
      operator*(const matrix<Real>& A, const vector<Real>& B) {
        return prod(A,B);
      }
           
      template<typename Real>
      inline
      matrix<Real>
      operator*(const matrix<Real>& A, const matrix<Real>& B) {
        return prod(A,B);
      }
      

      template <typename Real>
      std::ostream&
      operator<<(std::ostream& os, const matrix<Real>& A);
       
    }
  }
}
      

namespace Ariadne {
  namespace LinearAlgebra {

    using boost::numeric::ublas::identity_matrix;
    using boost::numeric::ublas::herm;
    using boost::numeric::ublas::vector;
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::matrix_column;
  
    template <typename Real>
    matrix<Real> zero_matrix(size_type r, size_type c);
  
    template <typename Real>
    inline
    matrix<Real> zero_matrix(size_type n)
    {
      return zero_matrix<Real>(n,n);
    }
    
  
    template<typename Real>
    Real
    norm(const matrix<Real>& A);
    
    template<typename Real>
    Real
    log_norm(const matrix<Real>& A);
    

    template <typename Real>
    matrix<Real> exp_Ah(const matrix<Real> &A, 
                        const Real& h, 
                        const Real& e); 

    template <typename Real> 
    vector<Real> exp_b_approx(const matrix<Real> &A, 
                                     const vector<Real> &b, 
                                     const Real h, const unsigned int n); 
   
    template <typename Real>
    void 
    lu_local_dec(matrix<Real>& A, 
                 const array<size_type>& row, const array<size_type>& col, 
                 const size_type& rows, const size_type& columns, 
                 const size_type& p);
   
    template <typename Real>
    matrix<Real> 
    lu_decompose(const matrix<Real>& A, 
                 array<size_type>& p_col, 
                 array<size_type>& p_row);
                              
    /* PAY ATTENTION!!! 
     * I supose that matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename Real>
    matrix<Real> 
    lu_decompose(const matrix<Real> &A, 
                 array<size_type>& p_array);

    
    /* PAY ATTENTION!!! 
     * I supose that boost::numeric::ublas::matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename Real>
    vector<Real> 
    lu_solve(const matrix<Real>& A, 
             const array<size_type>& p_array, 
             const vector<Real>& b);
             
    /* WARNING!!! The following function has some precision problems */ 
    template <typename Real>
    inline
    matrix<Real> 
    Householder_QR(const matrix<Real> &A);

    template <typename Real>
    matrix<Real>
    hermitian(const matrix<Real>& m);
    
    template <typename Real>
    matrix<Real> 
    inverse(const matrix<Real> &A);
    
    template <typename Real>
    Integer 
    common_denominator(const matrix<Real>& A);
    
    template <typename Real>
    vector<Integer> 
    row_common_denominators(const matrix<Real>& A);

    
    
    /* \brief Transforms the linear inequalities $Ax\leq b$ to $AT^{-1}y \leq b$. */
    template <typename Real>
    void 
    transform_linear_inequalities(const matrix<Real>& Real, 
                                  matrix<Real>& A, 
                                  vector<Real>& b);
    
    template <>
    void 
    transform_linear_inequalities<Dyadic>(const matrix<Dyadic>& Real, 
                                          matrix<Dyadic>& A, 
                                          vector<Dyadic>& b); 

    template <class Real>
    bool independent_rows(matrix<Real> A);

    template <typename Real>
    bool 
    have_same_dimensions(const matrix<Real> &A,  const matrix<Real> &B);
    
    template <typename Real>
    bool 
    equivalent_columns(const matrix<Real> &A, 
                       const size_type &A_col, 
                       const matrix<Real> &B, 
                       const size_type &B_col);
    
  
    template<typename Real>
    size_type 
    find_first_not_null_in_col(const matrix<Real> &A, 
                               const size_type &col);
    
    template <class Real>
    matrix<Real> 
    remove_null_columns_but_one(const matrix<Real> &A);
    
    template <typename Real>
    void 
    remove_null_columns(const matrix<Real>& A, 
                        array<size_type>& row, 
                        array<size_type>& col);

    template <typename Real>
    matrix<Real> 
    compute_space(const matrix<Real>& SA, 
                  array<size_type>& row,
                  const array<size_type>& col);
    
  }
}


#endif /* _ARIADNE_MATRIX_H */
