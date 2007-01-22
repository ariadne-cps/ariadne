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

#ifndef _ARIADNE_DIAGONAL_MATRIX_H
#define _ARIADNE_DIAGONAL_MATRIX_H

#include <iosfwd>

#include "../declarations.h"
#include "../numeric/integer.h"
#include "../numeric/interval.h"
#include "../linear_algebra/matrix.h"

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
    class DiagonalMatrix : public MatrixExpression< DiagonalMatrix<R> >
    {
      public:
      /*! \brief Construct a 0 by 0 matrix. */
      DiagonalMatrix() { }
      /*! \brief Construct an \a r by \a c matrix, all of whose entries are zero. */
      DiagonalMatrix(const size_type& n);
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr, 
       *  incrementing the input row elements by \a ri and input columns by \a ci. */
      DiagonalMatrix(const size_type& n, const R* ptr, const size_type& i);
      /*! \brief Construct an \a r by \a c matrix from a vector expression. */
      template<class E> DiagonalMatrix(const VectorExpression<E>& ve);

      /*! \brief The number of rows of the matrix. */
      size_type number_of_rows() const;
      /*! \brief The number of columns of the matrix. */
      size_type number_of_columns() const;
      
      /*! \brief A constant reference to \a i,\a j th element. */
      const R& operator() (const size_type& i, const size_type& j);
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

  }
}


#endif /* _ARIADNE_DIAGONAL_MATRIX_H */
