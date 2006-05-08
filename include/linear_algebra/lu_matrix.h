/***************************************************************************
 *            lu_matrix.h
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
 
/*! \file lu_matrix.h
 *  \brief LU factorisation and factorised matrices.
 */

#ifndef _ARIADNE_LU_MATRIX_H
#define _ARIADNE_LU_MATRIX_H

#include <blas/geset.hpp>
#include <blas/trcpy.hpp>
#include <blas/gemm.hpp>
#include <blas/trmm.hpp>
#include <lapack/laswp.hpp>
#include <lapack/getrf.hpp>
#include <lapack/getrs.hpp>

#include "../declarations.h"
#include "../base/array.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix stored in LU product form. */
    template<typename R>
    class LUMatrix {
     public:
      LUMatrix(const Matrix<R>& A);
      
      size_type number_of_rows() const { return _elements.number_of_rows(); }
      size_type number_of_columns() const { return _elements.number_of_columns(); }
      R operator() (const size_type& i, const size_type& j) const { 
        R result=0; for(int k=0; k<=BLAS::min(i,j); ++k) { result(i,j)+= (i==k ? R(1) : _elements(i,k)) * _elements(k,j); } return result; }

      Matrix<R> P() const;
      Matrix<R> L() const;
      Matrix<R> U() const;
      
      R determinant() const;
      bool singular() const;

      operator Matrix<R> () const;
      Matrix<R> inverse() const;
     private:
      Matrix<R> _elements;
      array<int> _pivots;
    };
    
    template <typename R>
    inline 
    LUMatrix<R>::LUMatrix(const Matrix<R>& A) 
      : _elements(A), _pivots(A.number_of_rows()) 
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      LAPACK::getrf(BLAS::RowMajor,m,n,this->_elements.begin(),n,
                    this->_pivots.begin());
    }

    template <typename R>
    inline
    Matrix<R>
    LUMatrix<R>::P() const
    {
      int m=this->number_of_rows();
      Matrix<R> result(m,m);
      BLAS::geset(BLAS::RowMajor,m,m,R(0),R(1),result.data().begin(),m);
      LAPACK::laswp(BLAS::RowMajor,m,result.data().begin(),m,
                    0,m,this->_pivots.begin(),-1);

      return result;
    }

    template <typename R>
    inline
    Matrix<R>
    LUMatrix<R>::L() const
    {
      int m=this->number_of_rows();
      Matrix<R> result(m,m);
      BLAS::geset(BLAS::RowMajor,m,m,R(0),R(1),result.begin(),m);
      BLAS::trcpy(BLAS::RowMajor,BLAS::Lower,BLAS::Unit,
                  m,m,
                  this->_elements.begin(), m,
                  result.begin(), m);
      return result;
    }

    template <typename R>
    inline
    Matrix<R>
    LUMatrix<R>::U() const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      Matrix<R> result(m,n);
      BLAS::geset(BLAS::RowMajor,m,n,R(0),R(0),result.begin(),n);
      BLAS::trcpy(BLAS::RowMajor,BLAS::Upper,BLAS::NonUnit,
                  m,n,
                  this->_elements.begin(),n,
                  result.begin(),n);
      return result;
    }

    template <typename R>
    inline
    LUMatrix<R>::operator Matrix<R> () const
    {
      size_type m=this->number_of_rows();
      size_type n=this->number_of_columns();
      Matrix<R> result(m,n);
        
      /* No BLAS routine */
      for(size_type i=0; i!=m; ++i) {
        for(size_type j=0; j!=n; ++j) {
          result(i,j)=0;
          for(size_type k=0; k<=BLAS::min(i,j); ++k) {
            if(i==k) { result(i,j)+=_elements(k,j); }
            else { result(i,j)+=_elements(i,k)*_elements(k,j); }
          }
        }
      }
      /* Fixme: Apply permutation */
      LAPACK::laswp(BLAS::RowMajor,n,result.data().begin(),n,0,n,this->_pivots.begin(),-1);
      return result;
    }

    template <typename R>
    inline
    bool
    LUMatrix<R>::singular() const
    {
      return this->determinant()==0; 
    }

    template <typename R>
    inline
    R
    LUMatrix<R>::determinant() const
    {
      assert(this->number_of_rows()==this->number_of_columns());
      size_type n=this->number_of_rows();
      R result=0;
      for(size_type i=0; i!=n; ++i) { 
        result+=this->_elements(i,i);
      }
      return result;
    }

    template <typename R>
    inline
    Matrix<R>
    LUMatrix<R>::inverse() const
    {
      size_type n=this->number_of_rows();
      Matrix<R> result(n,n);
      BLAS::geset(BLAS::RowMajor,n,n,R(0),R(1),result.begin(),n);
      LAPACK::getrs(BLAS::RowMajor,BLAS::NoTrans,n,n,this->_elements.begin(),n,this->_pivots.begin(),result.begin(),n);
      return result;
    }

  }
}

#endif /* _ARIADNE_LU_MATRIX_H */
