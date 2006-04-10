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

#include <iosfwd>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "../declarations.h"
#include "../numeric/integer.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix over \a R. */
    template<typename R>
    class matrix : public boost::numeric::ublas::matrix<R> 
    {
      typedef boost::numeric::ublas::matrix<R> Base;
      typedef typename numerical_traits<R>::field_extension_type F;
     public:
      matrix() : Base() { }
      matrix(const size_type& r, const size_type& c) : Base(r,c) { }
      template<typename E> matrix(const boost::numeric::ublas::matrix_expression<E>& A) : Base(A()) { }
      
      matrix(const std::string& s);

      R norm() const;
      R log_norm() const;
      
      matrix<F> inverse() const;
      vector<F> solve(const vector<R>& v) const;
    };
    
    using boost::numeric::ublas::identity_matrix;
    using boost::numeric::ublas::herm;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::matrix_column;
  
    template<typename R>
    inline
    vector<R>
    operator*(const matrix<R>& A, const vector<R>& B) {
      return boost::numeric::ublas::prod(A,B);
    }
           
    template<typename R>
    inline
    matrix<R>
    operator*(const matrix<R>& A, const matrix<R>& B) {
      return boost::numeric::ublas::prod(A,B);
    }
    

    template <typename R>
    matrix<R> zero_matrix(size_type r, size_type c);
  
    template <typename R>
    inline
    matrix<R> zero_matrix(size_type n)
    {
      return zero_matrix<R>(n,n);
    }
    
  
    template<typename R>
    inline
    R
    norm(const matrix<R>& A) 
    {
      return A.norm();
    }
    
    template<typename R>
    inline
    R
    log_norm(const matrix<R>& A)
    {
      return A.log_norm();
    }

    template <typename R>
    inline
    matrix<typename numerical_traits<R>::field_extension_type> 
    inverse(const matrix<R>& A)
    {
      return A.inverse();
    }
    
    template <typename R>
    inline
    vector<typename numerical_traits<R>::field_extension_type> 
    solve(const matrix<R>& A, const vector<R>& v)
    {
      return A.solve(v);
    }
    
    template<typename R>
    matrix<R>
    exp_approx(const matrix<R>& A, const R& e);
    
    template<typename R>
    matrix<R>
    concatenate_columns(const matrix<R>& A1, const matrix<R>& A2);

    template <typename R>
    void 
    lu_local_dec(matrix<R>& A, 
                 const array<size_type>& row, const array<size_type>& col, 
                 const size_type& rows, const size_type& columns, 
                 const size_type& p);
   
    template <typename R>
    matrix<R> 
    lu_decompose(const matrix<R>& A, 
                 array<size_type>& p_col, 
                 array<size_type>& p_row);
                              
    /* PAY ATTENTION!!! 
     * I supose that matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename R>
    matrix<R> 
    lu_decompose(const matrix<R> &A, 
                 array<size_type>& p_array);

    
    /* PAY ATTENTION!!! 
     * I supose that boost::numeric::ublas::matrix is row based i.e. 
     * A(i,j) is the element in the i-th row and in the j-th column 
     */
    template <typename R>
    vector<R> 
    lu_solve(const matrix<R>& A, 
             const array<size_type>& p_array, 
             const vector<R>& b);
             
    /* WARNING!!! The following function has some precision problems */ 
    template <typename R>
    inline
    matrix<R> 
    Householder_QR(const matrix<R> &A);

    template <typename R>
    matrix<R>
    hermitian(const matrix<R>& m);
    
    template <typename R>
    Integer 
    common_denominator(const matrix<R>& A);
    
    template <typename R>
    vector<Integer> 
    row_common_denominators(const matrix<R>& A);

    
    
    /* \brief Transforms the linear inequalities $Ax\leq b$ to $AT^{-1}y \leq b$. */
    template <typename R>
    void 
    transform_linear_inequalities(const matrix<R>& R, 
                                  matrix<R>& A, 
                                  vector<R>& b);
    
    template <class R>
    bool independent_rows(matrix<R> A);

    template <typename R>
    bool 
    have_same_dimensions(const matrix<R> &A,  const matrix<R> &B);
    
    template <typename R>
    bool 
    equivalent_columns(const matrix<R> &A, 
                       const size_type &A_col, 
                       const matrix<R> &B, 
                       const size_type &B_col);
    
  
    template<typename R>
    size_type 
    find_first_not_null_in_col(const matrix<R> &A, 
                               const size_type &col);
    
    template <class R>
    matrix<R> 
    remove_null_columns_but_one(const matrix<R> &A);
    
    template <typename R>
    void 
    remove_null_columns(const matrix<R>& A, 
                        array<size_type>& row, 
                        array<size_type>& col);

    template <typename R>
    matrix<R> 
    compute_space(const matrix<R>& SA, 
                  array<size_type>& row,
                  const array<size_type>& col);
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const matrix<R>& A);

    template <typename R>
    std::istream&
    operator>>(std::istream& is, matrix<R>& A);

  }
}


#endif /* _ARIADNE_MATRIX_H */
