 /***************************************************************************
 *            diagonal_matrix.h
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
 
/*! \file diagonal_matrix.h
 *  \brief Diagonal matrices.
 */

#ifndef ARIADNE_DIAGONAL_MATRIX_H
#define ARIADNE_DIAGONAL_MATRIX_H

#include <iosfwd>

#include "declarations.h"
#include "base/stlio.h"
#include "numeric/integer.h"
#include "numeric/interval.h"
#include "linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*!\ingroup LinearAlgebra
     * \brief A matrix over \a R. 
     *
     * \internal
     * Agreed on number_of_rows() and number_of_columns() for size functions as 
     * this is clearest and best English.
     */
    template<class R>
    class DiagonalMatrix
      : public MatrixExpression< DiagonalMatrix<R> >
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     private:
      array<R> _data;
     public:
      /*! \brief Construct a 0 by 0 matrix. */
      DiagonalMatrix() { }
      /*! \brief Construct an \a r by \a c matrix, all of whose entries are zero. */
      DiagonalMatrix(const size_type& n) : _data(n) { }
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr, 
       *  incrementing the input row elements by \a ri and input columns by \a ci. */
      DiagonalMatrix(const size_type& n, const R* ptr, const size_type& i=1) : _data(n) {
        for(size_type k=0; k!=n; ++k) { _data[k]=ptr[k*i]; } }
      /*! \brief Construct an \a r by \a c matrix from a vector expression. */
      template<class E> DiagonalMatrix(const VectorExpression<E>& ve);

      /*! \brief The number of rows of the matrix. */
      size_type number_of_rows() const { return this->_data.size(); }
      /*! \brief The number of columns of the matrix. */
      size_type number_of_columns() const { return this->_data.size(); };
      
      /*! \brief The data array of the matrix. */
      array<R>& data() { return this->_data; }
      /*! \brief A constant reference to the data array of the matrix. */
      const array<R>& data() const { return this->_data; }
      /*! \brief A constant reference to \a i,\a j th element. */
      const R& operator() (const size_type& i, const size_type& j) const;
      /*! \brief A reference to \a i,\a j th element. */
      R& operator() (const size_type& i, const size_type& j);
      
      // /*! \brief A slice through the \a i th row. */
      //VectorSlice<R> operator[] (const size_type& i) {
      //  return VectorSlice<R>(this->number_of_columns(),this->begin()+i*this->column_increment(),this->row_increment()); }
        
      /*! \brief An \a r by \a c matrix, all of whose entries are zero. */
      static DiagonalMatrix<R> zero(const size_type n);
      /*! \brief The \a n by \a n identity matrix. */
      static DiagonalMatrix<R> identity(const size_type n);
      
      /*! \brief The operator norm with respect to the supremum norm on vectors. 
       * Equal to the supremum over all rows of the sum of absolute values.
       */
      F norm() const;

      /*! \brief The logarithmic norm. */
      F log_norm() const;
      
      /*! \brief True if the matrix is singular. */
      bool singular() const;
      /*! \brief The determinant of the matrix. */
      F determinant() const;
      
      /*! \brief The transposed matrix. */
      Matrix<R> transpose() const;

      /*! \brief The inverse of the matrix. */
      Matrix<F> inverse() const;
      /*! \brief The solution of the linear equation \f$ Ax=b\f$. */
      Vector<F> solve(const Vector<R>& b) const;

      /*! \brief A pointer to the first element of the array of values. */
      //R* begin() { return &(*this)(0,0); }
      R* begin() { return this->data().begin(); }
      /*! \brief A constant pointer to the first element of the array of values. */
      //const R* begin() const { return const_cast< Matrix<R>* >(this)->begin(); }
      const R* begin() const { return this->data().begin(); }

      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);
    };

    template<class R> template<class E> inline
    DiagonalMatrix<R>::DiagonalMatrix(const VectorExpression<E>& ve) 
      : _data(ve().size())
    {
      for(size_type i=0; i!=ve().size(); ++i) {
        this->_data[i]=ve()[i];
      }
    }

    template<class X1, class X2>
    Matrix<typename Numeric::traits<X1,X2>::arithmetic_type>
    operator*(const Matrix<X1>& A, const DiagonalMatrix<X2>& D)
    {
      typedef typename Numeric::traits<X1,X2>::arithmetic_type X0;
      assert(A.number_of_columns()==D.number_of_rows());
      Matrix<X0> R(A.number_of_rows(),D.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          R(i,j)=A(i,j)*D.data()[j];
        }
      }
      return R;
    }

    template<class T>
    Matrix< Numeric::Float<T> >
    mul_approx(const Matrix< Numeric::Float<T> >& A, const DiagonalMatrix< Numeric::Float<T> >& D)
    {
      assert(A.number_of_columns()==D.number_of_rows());
      Matrix< Numeric::Float<T> > R(A.number_of_rows(),D.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          R(i,j)=mul_approx(A(i,j),D.data()[j]);
        }
      }
      return R;
    }

    template<class R>
    std::ostream& operator<<(std::ostream& os, const DiagonalMatrix<R>& D) {
      return os <<"DiagonalMatrix("<<D.data()<<")";
    }

  }
}


#endif /* ARIADNE_DIAGONAL_MATRIX_H */
