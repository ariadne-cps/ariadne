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

    using boost::numeric::ublas::identity_matrix;
    using boost::numeric::ublas::herm;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::matrix_column;
  
    /*! \ingroup LinearAlgebra
     *  \brief A matrix over \a R. 
     */
    template<typename R>
    class Matrix : public boost::numeric::ublas::matrix<R> 
    {
      typedef boost::numeric::ublas::matrix<R> _Base;
      typedef boost::numeric::ublas::matrix<R> _boost_matrix;
      typedef typename numerical_traits<R>::field_extension_type F;
     public:
      /*! \brief Construct a 0 by 0 matrix. */
      Matrix() : _Base() { }
      /*! \brief Construct an \a r by \a c matrix, all of whose entries are zero. */
      Matrix(const size_type& r, const size_type& c) : _Base(r,c) { }
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr, 
       *  incrementing the input row elements by \a ri and input columns by \a ci. */
      Matrix(const size_type& r, const size_type& c, 
             const R* ptr, const size_type& ri, const size_type& ci=1) : _Base(r,c) 
      { 
        for(size_type i=0; i!=r; ++i) { 
          for(size_type j=0; j!=c; ++j) { 
            _Base::operator()(i,j)=ptr[i*r+j*ci];
          } 
        } 
      }

      template<typename E> Matrix(const boost::numeric::ublas::matrix_expression<E>& A) : _Base(A) { }
            
      /*! \brief Construct from a string literal of the form "[a11,a12,...,a1n; a21,a22,...,a2n;...;am1,am2,...amn]". */
      explicit Matrix(const std::string& s);

#ifdef DOXYGEN
      /*! \brief Copy constructor. */
      Matrix(const Matrix<R>& A);
      /*! \brief Copy assignment operator. */
      Matrix<R>& operator=(const Matrix<R>& A);
#endif
      
      /*! \brief The equality operator. */
      bool operator==(const Matrix<R>& A) const;
      /*! \brief The inequality operator. */
      bool operator!=(const Matrix<R>& A) const;

    
     
      /*! \brief The number of elements in the \a d th dimension. */
      size_type size(const size_type& d) const {
        if(d==0) { return number_of_rows(); } 
        else if(d==1) { return number_of_columns(); }
        else { assert(false); } 
      }
      /*! \brief The number of rows of the matrix. */
      size_type number_of_rows() const { return _Base::size1(); }
      /*! \brief The number of columns of the matrix. */
      size_type number_of_columns() const { return _Base::size2(); }
      
      /*! \brief A constant reference to \a i,\a j th element. */
      const R& operator() (const size_type& i, const size_type& j) const { 
        return this->_Base::operator()(i,j); }
      /*! \brief A reference to \a i,\a j th element. */
      R& operator() (const size_type& i, const size_type& j) { 
        return this->_Base::operator()(i,j); }
      
      /*! \brief An \a r by \a c matrix, all of whose entries are zero. */
      static Matrix<R> zero(const size_type r, const size_type c) {
        return Matrix<R>(r,c); }
      /*! \brief The \a n by \a n identity matrix. */
      static Matrix<R> identity(const size_type n) {
        Matrix<R> result(n,n); 
        for(size_type i=0; i!=n; ++i) { result(i,i)=R(1); } 
        return result;
      }
      
      /*! \brief Concatenate the columns of two matrices. */
      static Matrix<R> concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2);
      
      /*! \brief The operator norm with respect to the supremum norm on vectors. 
       * Equal to the supremum over all rows of the sum of absolute values.
       */
      R norm() const;

      /*! \brief The logarithmic norm. */
      R log_norm() const;
      
      /*! \brief True if the matrix is singular. */
      bool singular() const;
      /*! \brief The determinant of the matrix. */
      R determinant() const;
      
      /*! \brief The transposed matrix. */
      Matrix<R> transpose() const;

      /*! \brief The inverse of the matrix. */
      Matrix<R> inverse() const;
      /*! \brief The solution of the linear equation \f$ Ax=b\f$. */
      Vector<R> solve(const Vector<R>& b) const;

      /*! \brief A pointer to the first element of the array of values. */
      R* begin() { return &(*this)(0,0); }
      /*! \brief A constant pointer to the first element of the array of values. */
      const R* begin() const { return const_cast< Matrix<R>* >(this)->begin(); }

      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);

#ifdef DOXYGEN
      /*! \brief The additive inverse of the matrix \a A. */
      friend Matrix<R> operator-<>(const Matrix<R>& A);
      /*! \brief The sum of \a A1 and \a A2. */
      friend Matrix<R> operator+<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The difference of \a A1 and \a A2. */
      friend Matrix<R> operator-<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<R> operator*<>(const R& s, const Matrix<R>& A);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<R> operator*<>(const Matrix<R>& A, const R& s);
      /*! \brief The scalar product of \a A by the reciprocal of \a s. */
      friend Matrix<R> operator/<>(const Matrix<R>& A, const R& s);
      /*! \brief The product of matrix \a A1 with matrix \a A2. */
      friend Matrix<R> operator*<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The product of matrix \a A with vector \a v. */
      friend Vector<R> operator*<>(const Matrix<R>& A, const Vector<R>& v);

      /*! \brief The inverse of \a A. */
      friend Matrix<R> inverse<>(const Matrix<R>& A);
      /*! \brief Solve the linear system \f$Ax=b\f$. */
      friend Matrix<R> solve<>(const Matrix<R>& A, const Vector<R>& b);

      /*! \brief Catenate the columns of \a A1 with those of \a A2. */
      friend Matrix<R> concatenate_columns<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief Checks if the matrices have the same dimensions. */
      friend bool have_same_dimensions<>(const Matrix<R>& A1, const Matrix<R>& A2);
#endif 

    };
  
    template<typename R>
    inline
    bool
    have_same_dimensions(const Matrix<R>& A1, const Matrix<R>& A2) 
    {
      return A1.number_of_rows()==A2.number_of_rows() 
          && A1.number_of_columns()==A2.number_of_columns();
    }
    
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
    Matrix<R> 
    inverse(const Matrix<R>& A)
    {
      return A.inverse();
    }
    
    template <typename R>
    inline
    Vector<R> 
    solve(const Matrix<R>& A, const Vector<R>& v)
    {
      return A.solve(v);
    }
    
    template<typename R>
    inline
    Matrix<R>
    concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2) {
      return Matrix<R>::concatenate_columns(A1,A2);
    }
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Matrix<R>& A) {
      return A.write(os);
    }

    template <typename R>
    std::istream&
    operator>>(std::istream& is, Matrix<R>& A) {
      return A.read(is); 
    }

        
             
    
/*        
        
    
    
    template <class R>
    bool independent_rows(Matrix<R> A);

    
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
*/

  }
}


#endif /* _ARIADNE_MATRIX_H */
