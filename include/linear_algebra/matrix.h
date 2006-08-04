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
 *  \brief Matrices and Matrix operations.
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
    class Matrix : public boost::numeric::ublas::matrix<R> 
    {
      typedef boost::numeric::ublas::matrix<R> Base;
      typedef typename numerical_traits<R>::field_extension_type F;
     public:
      Matrix() : Base() { }
      Matrix(const size_type& r, const size_type& c) : Base(r,c) { }
      Matrix(const size_type& r, const size_type& c, 
             const R* ptr, const size_type& ld) : Base(r,c) 
      { 
        for(size_type i=0; i!=r; ++i) { 
          for(size_type j=0; j!=c; ++j) { 
            Base::operator()(i,j)=ptr[i*ld+j];
          } 
        } 
      }

      template<typename E> Matrix(const boost::numeric::ublas::matrix_expression<E>& A) : Base(A()) { }
            
      Matrix(const std::string& s);

      bool operator==(const Matrix<R>& other) const { 
        const Matrix<R>& self=*this;
        if(self.number_of_rows() != other.number_of_rows() ||
           self.number_of_columns() != other.number_of_columns()) 
        {
          return false; 
        }
        for(size_type i=0; i!=this->number_of_rows(); ++i) {
          for(size_type j=0; j!=this->number_of_columns(); ++j) {
            if(self(i,j)!=other(i,j)) {
              return false;
            }
          }
        }
        return true;
      }

      size_type number_of_rows() const { return Base::size1(); }
      size_type number_of_columns() const { return Base::size2(); }

      R norm() const;
      R log_norm() const;
      
      Matrix<R> transpose() const;

      bool singular() const;
      R determinant() const;
      
      Matrix<F> inverse() const;
      Vector<F> solve(const Vector<R>& v) const;

      R* begin() { return &(*this)(0,0); }
      const R* begin() const { return const_cast< Matrix<R>* >(this)->begin(); }
    };
    
    using boost::numeric::ublas::identity_matrix;
    using boost::numeric::ublas::herm;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::matrix_column;
  
    template<typename R>
    inline
    Vector<R>
    operator*(const Matrix<R>& A, const Vector<R>& B) {
      return boost::numeric::ublas::prod(A,B);
    }
           
    template<typename R>
    inline
    Matrix<R>
    operator*(const Matrix<R>& A, const Matrix<R>& B) {
      return boost::numeric::ublas::prod(A,B);
    }
    

    template <typename R>
    Matrix<R> zero_Matrix(size_type r, size_type c);
  
    template <typename R>
    inline
    Matrix<R> zero_Matrix(size_type n)
    {
      return zero_Matrix<R>(n,n);
    }
    
  
    template<typename R>
    inline
    R
    norm(const Matrix<R>& A) 
    {
      return A.norm();
    }
    
    template<typename R>
    inline
    R
    log_norm(const Matrix<R>& A)
    {
      return A.log_norm();
    }

    template <typename R>
    inline
    Matrix<typename numerical_traits<R>::field_extension_type> 
    inverse(const Matrix<R>& A)
    {
      return A.inverse();
    }
    
    template <typename R>
    inline
    Vector<typename numerical_traits<R>::field_extension_type> 
    solve(const Matrix<R>& A, const Vector<R>& v)
    {
      return A.solve(v);
    }
    
    template<typename R>
    Matrix<R>
    exp_approx(const Matrix<R>& A, const R& e);
    
    template<typename R>
    Matrix<R>
    concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2);

    template <typename R>
    void 
    lu_local_dec(Matrix<R>& A, 
                 const array<size_type>& row, const array<size_type>& col, 
                 const size_type& rows, const size_type& columns, 
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
             
    /* WARNING!!! The following function has some precision problems */ 
    template <typename R>
    inline
    Matrix<R> 
    Householder_QR(const Matrix<R> &A);

    template <typename R>
    Matrix<R>
    hermitian(const Matrix<R>& m);
    
    template <typename R>
    Integer 
    common_denominator(const Matrix<R>& A);
    
    template <typename R>
    Vector<Integer> 
    row_common_denominators(const Matrix<R>& A);

    
    
    /* \brief Transforms the linear inequalities $Ax\leq b$ to $AT^{-1}y \leq b$. */
    template <typename R>
    void 
    transform_linear_inequalities(const Matrix<R>& R, 
                                  Matrix<R>& A, 
                                  Vector<R>& b);
    
    template <class R>
    bool independent_rows(Matrix<R> A);

    template <typename R>
    bool 
    have_same_dimensions(const Matrix<R> &A,  const Matrix<R> &B);
    
    template <typename R>
    bool 
    equivalent_columns(const Matrix<R> &A, 
                       const size_type &A_col, 
                       const Matrix<R> &B, 
                       const size_type &B_col);
    
  
    template<typename R>
    size_type 
    find_first_not_null_in_col(const Matrix<R> &A, 
                               const size_type &col);
    
    template <class R>
    Matrix<R> 
    remove_null_columns_but_one(const Matrix<R> &A);
    
    template <typename R>
    void 
    remove_null_columns(const Matrix<R>& A, 
                        array<size_type>& row, 
                        array<size_type>& col);

    template <typename R>
    Matrix<R> 
    compute_space(const Matrix<R>& SA, 
                  array<size_type>& row,
                  const array<size_type>& col);
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Matrix<R>& A);

    template <typename R>
    std::istream&
    operator>>(std::istream& is, Matrix<R>& A);

  }
}


#endif /* _ARIADNE_MATRIX_H */
