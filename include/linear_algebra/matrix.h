/***************************************************************************
 *            matrix.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_MATRIX_H
#define ARIADNE_MATRIX_H

#include <iosfwd>

#include "base/types.h"
#include "base/array.h"
#include "numeric/declarations.h"
#include "numeric/traits.h"
#include "numeric/integer.h"
#include "numeric/interval.h"

#include "linear_algebra/exceptions.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix_expression.h"

namespace Ariadne {
  

    class Slice;
    template<class R> class Matrix;
    template<class R> class MatrixSlice;

    /*!\ingroup LinearAlgebra
     * \brief A matrix over \a R. 
     *
     * \internal
     * Agreed on number_of_rows() and number_of_columns() for size functions as 
     * this is clearest and best English.
     */
    template<class R>
    class Matrix : public MatrixExpression< Matrix<R> > 
    {
     private:
      size_type _nr;
      size_type _nc;
      array<R> _array;
     private:
      typedef typename traits<R>::arithmetic_type F;
     public:
      /*! \brief The type of real number stored by the matrix. */
      typedef R value_type;
       
      /*! \brief Construct a 0 by 0 matrix. */
      explicit Matrix();
      /*! \brief Construct an \a r by \a c matrix, all of whose entries are zero. */
      explicit Matrix(const size_type& r, const size_type& c);
      /*! \brief Construct an \a r by \a c matrix from the one-dimensional array (in row-major format) beginning at \a ptr. */
      template<class RR> explicit Matrix(const size_type& nr, const size_type& nc,const RR* ptr);
      /*! \brief Construct an \a r by \a c matrix from the two-dimensional C-style array \a ary. */
      template<int NC, class RR> explicit Matrix(const size_type& nr, const size_type& nc,const RR ary[][NC]);
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr, 
       *  incrementing the input row elements by \a ri and input columns by \a ci. */
      template<class RR> explicit Matrix(const size_type& nr, const size_type& nc, 
                                         const RR* ptr, const size_type& ri, const size_type& ci=1u);

      /*! \brief Convert from a matrix with a different numerical type. */
      template<class RR> Matrix(const Matrix<RR>& A);
      /*! \brief Convert from a matrix expression. */
      template<class E> Matrix(const MatrixExpression<E>& A);
        
      /*! \brief Construct from a string literal of the form "[a11,a12,...,a1n; a21,a22,...,a2n;...;am1,am2,...amn]". */
      explicit Matrix(const std::string& s);

      /*! \brief Copy constructor. */
      Matrix(const Matrix<R>& A);
      /*! \brief Copy assignment operator. */
      Matrix<R>& operator=(const Matrix<R>& A);
      /*! \brief Assignment operator. */
      template<class E> Matrix<R>& operator=(const MatrixExpression<E>& A);
      
      /*! \brief The equality operator. */
      bool operator==(const Matrix<R>& A) const;
      /*! \brief The inequality operator. */
      bool operator!=(const Matrix<R>& A) const;

      /*! A reference to the array holding the element data. */
      array<R>& data();
      /*! A constant reference to the array holding the element data. */
      const array<R>& data() const;
      
      /*! \brief Resize the matrix. */
      void resize(const size_type& nr, const size_type nc);
      
      /*! \brief The number of rows of the matrix. (Deprecated) \deprecated */
      size_type number_of_rows() const;
      /*! \brief The number of columns of the matrix. (Deprecated) \deprecated */
      size_type number_of_columns() const;
      /*! \brief The number of rows of the matrix. */
      size_type row_size() const;
      /*! \brief The number of columns of the matrix. */
      size_type column_size() const;
      /*! \brief A two-element array giving the number of rows and number of columns of the matrix. */
      array<size_type,2u> size() const;

      /*! \brief A constant reference to \a i,\a j th element. */
      const R& operator() (const size_type& i, const size_type& j) const;
      /*! \brief A reference to \a i,\a j th element. */
      R& operator() (const size_type& i, const size_type& j);
      
      /*! \brief A constant slice through the matrix. */
      MatrixSlice<const R> operator() (const Slice& i, const Slice& j) const;
      /*! \brief A slice through the matrix. */
      MatrixSlice<R> operator() (const Slice& i, const Slice& j);
      
      /*! \brief Subscripting operator; A[i][j] returns a reference to the (i,j)th element. */
      MatrixRow< Matrix<R> > operator[](const size_type& i);
        
      /*! \brief Constant subscripting operator; A[i][j] returns the (i,j)th element. */
      MatrixRow< const Matrix<R> > operator[](const size_type& i) const;
        
      /*! \brief An \a r by \a c matrix, all of whose entries are zero. */
      static Matrix<R> zero(const size_type r, const size_type c);
        
      /*! \brief An \a r by \a c matrix, all of whose entries are one. */
      static Matrix<R> one(const size_type r, const size_type c);
        
      /*! \brief The \a n by \a n identity matrix. */
      static Matrix<R> identity(const size_type n);
      
      /*! \brief The \a i th row. */
      VectorSlice<const R> row(const size_type& i) const;
      /*! \brief The \a j th column. */
      VectorSlice<const R> column(const size_type& j) const;
        
      /*! \brief The \a i th row. */
      VectorSlice<R> row(const size_type& i);
      /*! \brief The \a j th column. */
      VectorSlice<R> column(const size_type& j);
        
      /*! \brief A pointer to the first element of the array of values. */
      R* begin();
      /*! \brief A constant pointer to the first element of the array of values. */
      const R* begin() const;

      /*! \brief The increment to the array element needed to move one row down. */
      size_type row_increment() const;
      /*! \brief The increment to the array element needed to move one column right. */
      size_type column_increment() const;
      
      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);

#ifdef DOXYGEN
      typedef traits<R>::arithmetic_type AT;
      /*! \brief The additive inverse of the matrix \a A. */
      friend Matrix<R> operator-(const Matrix<R>& A);
      /*! \brief The sum of \a A1 and \a A2. */
      friend Matrix<AT> operator+(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The difference of \a A1 and \a A2. */
      friend Matrix<AT> operator-(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<AT> operator*(const R& s, const Matrix<R>& A);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<AT> operator*(const Matrix<R>& A, const R& s);
      /*! \brief The scalar product of \a A by the reciprocal of \a s. */
      friend Matrix<RAT> operator/(const Matrix<R>& A, const R& s);
      /*! \brief The product of matrix \a A1 with matrix \a A2. */
      friend Matrix<AT> operator*(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The product of matrix \a A with vector \a v. */
      friend Vector<AT> operator*(const Matrix<R>& A, const Vector<R>& v);
      /*! \brief The product of vector \a v with matrix \a A. */
      friend Vector<AT> operator*(const Vector<R>& v, const Matrix<R>& A);

      /*! \brief The inverse of \a A. */
      friend Matrix<AT> inverse(const Matrix<R>& A);
      /*! \brief Solve the linear system \f$Ax=b\f$. */
      friend Matrix<AT> solve(const Matrix<R>& A, const Vector<R>& b);

      /*! \brief Compute an approximate solution to the linear system \f$Ax=b\f$. */
      friend Matrix<R> solve_approx(const Matrix<R>& A, const Vector<R>& b);
      /*! \brief Compute an approximate inverse to the matrix \f$A\f$. */
      friend Matrix<R> inverse_approx(const Matrix<R>& A);
      /*! \brief Compute an approximate QR factorisation of the matrix \f$A\f$. */
      friend std::pair< Matrix<R>,Matrix<R> > qr_approx(const Matrix<R>& A);

      /*! \brief Catenate the rows of \a A with the row vector \a v. */
      friend Matrix<R> concatenate_rows(const Matrix<R>& A, const Vector<R>& v);
      /*! \brief Catenate the rows of \a A1 with those of \a A2. */
      friend Matrix<R> concatenate_rows(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief Catenate the rows of \a A with the row vector \a v. */
      friend Matrix<R> concatenate_columns(const Matrix<R>& A, const Vector<R>& v);
      /*! \brief Catenate the columns of \a A1 with those of \a A2. */
      friend Matrix<R> concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2);
#endif 
     private:
      static void instantiate();
    };
  


    template<class R> Matrix<R> approximation(const Matrix<R>& v);
    template<class R> Matrix<R> approximation(const Matrix< Interval<R> >& v);

    template<class R> Matrix<R> midpoint(const Matrix< Interval<R> >& iA); 
    template<class R> bool encloses(const Matrix< Interval<R> >& iA, const Matrix<R>& A); 
    template<class R> bool refines(const Matrix< Interval<R> >& iA1, const Matrix< Interval<R> >& iA2); 



    template<class R> Vector< Interval<R> > radius_row_sum(const Matrix< Interval<R> >& im); 
    
    template<class R> bool have_same_dimensions(const Matrix<R>& A1, const Matrix<R>& A2); 
    
    
    template<class R1, class R2>  Matrix<R1>& operator+=(Matrix<R1>& A1, const Matrix<R2>& A2); 
    template<class R1, class R2>  Matrix<R1>& operator-=(Matrix<R1>& A1, const Matrix<R2>& A2); 

    template<class R> Matrix<R> operator-(const Matrix<R>& A);
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator+(const Matrix<R1>& A1, const Matrix<R2>& A2); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator-(const Matrix<R1>& A1, const Matrix<R2>& A2); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const R1& s, const Matrix<R2>& A); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const Matrix<R1>& A, const R2& s); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator/(const Matrix<R1>& A, const R2& s); 
    template<class R1, class R2> Vector<typename traits<R1,R2>::arithmetic_type> operator*(const Matrix<R1>& A, const Vector<R2>& v); 
    template<class R1, class R2> Vector<typename traits<R1,R2>::arithmetic_type> operator*(const Vector<R1>& v, const Matrix<R2>& A); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const Matrix<R1>& A1, const Matrix<R2>& A2); 


    template<class R> Matrix<R> operator-(const MatrixSlice<R>& A);
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const MatrixSlice<R1>& A1, const Matrix<R2>& A2); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const MatrixSlice<R1>& A1, const MatrixSlice<R2>& A2); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const R1& s, const MatrixSlice<R2>& A); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const Matrix<R1>& A1, const MatrixSlice<R2>& A2); 
    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> operator*(const MatrixSlice<R1>& A, const R2& s); 
    template<class R1, class R2> Vector<typename traits<R1,R2>::arithmetic_type> operator*(const MatrixSlice<R1>& A, const Vector<R2>& v); 

    template<class R1, class R2> Matrix<typename traits<R1,R2>::arithmetic_type> outer_product(const Vector<R1>& v1, const Vector<R2>& v2); 
    
    template<class R> Vector<typename traits<R>::arithmetic_type> row_norms(const Matrix<R>& A); 
    
    template<class R> typename traits<R>::arithmetic_type norm(const Matrix<R>& A); 
    template<class R> typename traits<R>::arithmetic_type norm(const MatrixSlice<R>& A); 
    template<class R> typename traits<R>::arithmetic_type sup_norm(const Matrix<R>& A); 
    template<class R> typename traits<R>::arithmetic_type sup_norm(const MatrixSlice<R>& A); 
    template<class R> typename traits<R>::arithmetic_type log_norm(const Matrix<R>& A);

    template<class R> bool  singular(const Matrix<R>& A);
    template<class R> typename traits<R>::arithmetic_type  determinant(const Matrix<R>& A);
    template<class R> Matrix<R>  transpose(const Matrix<R>& A);
    template<class R> Matrix<typename traits<R>::arithmetic_type> inverse(const Matrix<R>& A);
    template<class R1, class R2> Vector<typename traits<R1,R2>::arithmetic_type> solve(const Matrix<R1>& A, const Vector<R2>& v);
    
    template<class T> Vector< Float<T> > mul_approx(const Matrix< Float<T> >& A, const Vector< Float<T> >& v);
    template<class T> Matrix< Float<T> > mul_approx(const Matrix< Float<T> >& A, const Matrix< Float<T> >& B);
    template<class T> Vector< Float<T> > solve_approx(const Matrix< Float<T> >& A, const Vector< Float<T> >& v);
    template<class T> Matrix< Float<T> > inverse_approx(const Matrix< Float<T> >& A);
    template<class T> std::pair<Matrix<Float<T> >,Matrix<Float<T> > > qr_approx(const Matrix< Float<T> >& A);

    template<class R> Matrix<R> direct_sum(const Matrix<R>& A1, const Matrix<R>& A2);

    template<class R> MatrixSlice<R> project(Matrix<R>& A, const Slice& rs, const Slice& cs);
    template<class R> MatrixSlice<const R> project(const Matrix<R>& A, const Slice& rs, const Slice& cs);
    template<class R> Matrix<R> concatenate(const Matrix<R>& A1, const Matrix<R>& A2);
    template<class R> Matrix<R> concatenate_rows(const Matrix<R>& A, const Vector<R>& v);
    template<class R> Matrix<R> concatenate_rows(const Matrix<R>& A1, const Matrix<R>& A2);
    template<class R> Matrix<R> concatenate_columns(const Matrix<R>& A, const Vector<R>& v);
    template<class R> Matrix<R> concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2);
    template<class R> Matrix<R> concatenate_columns(const MatrixSlice<R>& A1, const MatrixSlice<R>& A2);
    
    template<class R> Matrix<R> over_approximation(const Matrix< Interval<R> >& A);

    template<class R> std::ostream& operator<<(std::ostream& os, const Matrix<R>& A);
    template<class R> std::istream& operator>>(std::istream& is, Matrix<R>& A);
   
    template<class R> Matrix<typename traits<R>::arithmetic_type> schulz_inverse(const Matrix<R>& A);


} // namespace Ariadne

#include "matrix_slice.h"
#include "matrix.inline.h"
#include "matrix.template.h"

#endif /* ARIADNE_MATRIX_H */
