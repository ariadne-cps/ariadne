/***************************************************************************
 *            matrix_slice.h
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
 
/*! \file matrix_slice.h
 *  \brief Matrix slices
 */

#ifndef ARIADNE_MATRIX_SLICE_H
#define ARIADNE_MATRIX_SLICE_H

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
    template<class R> class MatrixSlice<const R>;

    /*!\ingroup LinearAlgebra
     * \brief A slice through a matrix with equally spaced row and column increments. 
     */
    template<class R>
    class MatrixSlice
      : public MatrixExpression< MatrixSlice<R> >
    {
     public:
      typedef R value_type;
      MatrixSlice(const size_type& nr, const size_type& nc, R* ptr, const size_type& rinc, const size_type& cinc=1u) 
        : _rs(nr), _cs(nc), _ptr(ptr), _ri(rinc), _ci(cinc) { }
      MatrixSlice(const size_type& nr, const size_type& nc, R* ptr)
        : _rs(nr), _cs(nc), _ptr(ptr), _ri(nc), _ci(1u) { }
      MatrixSlice(Matrix<R>& m)
        : _rs(m.number_of_rows()), _cs(m.number_of_columns()), _ptr(m.data().begin()),
          _ri(m.row_increment()), _ci(m.column_increment()) { }

      size_type number_of_rows() const { return this->_rs; }
      size_type number_of_columns() const { return this->_cs; };
      array<size_type,2u> size() const { array<size_type,2u> r; r[0]=this->_rs; r[1]=this->_cs; return r; }
      size_type row_increment() const { return this->_ri; }
      size_type column_increment() const { return this->_ci; };
      R* begin() { return this->_ptr; }
      const R* begin() const { return this->_ptr; }
      R& operator() (const size_type& i, const size_type& j) { return this->_ptr[this->_ri*i+this->_ci*j]; }
      const R& operator() (const size_type& i, const size_type& j) const { return this->_ptr[this->_ri*i+this->_ci*j]; }
      VectorSlice<R> operator[] (const size_type& i) { return VectorSlice<R>(this->_nc,this->_ptr+this->_ri*i,this->_ci); }
      VectorSlice<const R> operator[] (const size_type& i) const { return VectorSlice<const R>(this->_nc,this->_ptr+this->_ri*i,this->_ci); }
      VectorSlice<R> row(const size_type& i) { return VectorSlice<R>(this->_nc,this->_ptr+this->_ri*i,this->_ci); }
      VectorSlice<R> column(const size_type& j) { return VectorSlice<R>(this->_nr,this->_ptr+this->_ci*j,this->_ri); }
      VectorSlice<const R> row(const size_type& i) const { return VectorSlice<const R>(this->_nc,this->_ptr+this->_ri*i,this->_ci); }
      VectorSlice<const R> column(const size_type& j) const { return VectorSlice<const R>(this->_nr,this->_ptr+this->_ci*j,this->_ri); }

      MatrixSlice<R>& operator=(const R& m);
      template<class E> MatrixSlice<R>& operator=(const MatrixExpression< E >& m);
     private:
      size_type _rs;
      size_type _cs;
      R* _ptr;
      size_type _ri;
      size_type _ci;
    };
    

    template<class R>
    class MatrixSlice<const R>
      : public MatrixExpression< MatrixSlice<const R> >
    {
     public:
      typedef R value_type;
      MatrixSlice(const size_type& nr, const size_type& nc, const R* ptr, const size_type& rinc, const size_type& cinc=1u) 
        : _rs(nr), _cs(nc), _ptr(ptr), _ri(rinc), _ci(cinc) { }
      MatrixSlice(const size_type& nr, const size_type& nc, const R* ptr)
        : _rs(nr), _cs(nc), _ptr(ptr), _ri(nc), _ci(1u) { }
      MatrixSlice(const Matrix<R>& m)
        : _rs(m.number_of_rows()), _cs(m.number_of_columns()), _ptr(m.data().begin()),
          _ri(m.row_increment()), _ci(m.column_increment()) { }

      size_type number_of_rows() const { return this->_rs; }
      size_type number_of_columns() const { return this->_cs; };
      array<size_type,2u> size() const { array<size_type,2u> r; r[0]=this->_rs; r[1]=this->_cs; return r; }
      size_type row_increment() const { return this->_ri; }
      size_type column_increment() const { return this->_ci; };
      const R* begin() const { return this->_ptr; }
      const R& operator() (const size_type& i, const size_type& j) const { return this->_ptr[this->_ri*i+this->_ci*j]; }
      VectorSlice<const R> operator[] (const size_type& i) const { return VectorSlice<const R>(this->_nc,this->_ptr+this->_ri*i,this->_ci); }
      VectorSlice<const R> row(const size_type& i) const { return VectorSlice<const R>(this->_nc,this->_ptr+this->_ri*i,this->_ci); }
      VectorSlice<const R> column(const size_type& j) const { return VectorSlice<const R>(this->_nr,this->_ptr+this->_ci*j,this->_ri); }
     private:
      size_type _rs;
      size_type _cs;
      const R* _ptr;
      size_type _ri;
      size_type _ci;
    };
    
    template<class R> template<class E> inline 
    MatrixSlice<R>& MatrixSlice<R>::operator=(const MatrixExpression<E>& e) {
      const E& A=e();
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          (*this)(i,j)=A(i,j);
        }
      }
      return *this;
    }

      
    template<class R> inline std::ostream& operator<<(std::ostream& os, const MatrixSlice<R>& A) {
      return os << Matrix<R>(A); }
    template<class R> inline std::ostream& operator<<(std::ostream& os, const MatrixSlice<const R>& A) {
      return os << Matrix<R>(A); }


} // namespace Ariadne

#endif /* ARIADNE_MATRIX_SLICE_H */
