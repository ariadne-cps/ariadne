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

#ifndef ARIADNE_LU_MATRIX_H
#define ARIADNE_LU_MATRIX_H

#include <algorithm>

#include <tblas/geset.hpp>
#include <tblas/trcpy.hpp>
#include <tblas/gemm.hpp>
#include <tblas/trmm.hpp>
#include <tlapack/laswp.hpp>
#include <tlapack/getrf.hpp>
#include <tlapack/getrs.hpp>

#include "../base/array.h"
#include "../linear_algebra/exceptions.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \ingroup LinearAlgebra
     *  \brief A matrix stored in LU product form. 
     */
    template<class Real>
    class LUMatrix {
     public:
      /*! \brief Construct from an ordinary matrix. */
      LUMatrix(const Matrix<Real>& A);
      
      /*! \brief The number of rows. */
      size_type number_of_rows() const { return _elements.number_of_rows(); }
      /*! \brief The number of columns. */
      size_type number_of_columns() const { return _elements.number_of_columns(); }
      /*! \brief The \a i,\a j th element. */
      Real operator() (const size_type& i, const size_type& j) const { 
        Real result=0; for(int k=0; k<=std::min(i,j); ++k) { result(i,j)+= (i==k ? Real(1) : _elements(i,k)) * _elements(k,j); } return result; }

      /*! \brief The pivoting permutation. */
      Matrix<Real> P() const;
      /*! \brief The lower-triangular factor. */
      Matrix<Real> L() const;
      /*! \brief The upper-triangular factor. */
      Matrix<Real> U() const;
      
      /*! \brief The determinant. */
      Real determinant() const;
      /*! \brief True if the matrix is singular (non-invertible). */
      bool singular() const;

      /*! \brief Convert to an ordinary matrix. */
      operator Matrix<Real> () const;
      /*! \brief The inverse. */
      Matrix<Real> inverse() const;
     private:
      Matrix<Real> _elements;
      array<int> _pivots;
    };
    
    template<class Real>
    inline 
    LUMatrix<Real>::LUMatrix(const Matrix<Real>& A) 
      : _elements(A), _pivots(A.number_of_rows()) 
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      TLAPACK::getrf(TBLAS::RowMajor,m,n,this->_elements.begin(),n,
                    this->_pivots.begin());      
    }

    template<class Real>
    inline
    Matrix<Real>
    LUMatrix<Real>::P() const
    {
      int m=this->number_of_rows();
      Matrix<Real> result(m,m);
      TBLAS::geset(TBLAS::RowMajor,m,m,Real(0),Real(1),result.data().begin(),m);
      TLAPACK::laswp(TBLAS::RowMajor,m,result.data().begin(),m,
                    0,m,this->_pivots.begin(),-1);

      return result;
    }

    template<class Real>
    inline
    Matrix<Real>
    LUMatrix<Real>::L() const
    {
      int m=this->number_of_rows();
      Matrix<Real> result(m,m);
      TBLAS::geset(TBLAS::RowMajor,m,m,Real(0),Real(1),result.begin(),m);
      TBLAS::trcpy(TBLAS::RowMajor,TBLAS::Lower,TBLAS::Unit,
                  m,m,
                  this->_elements.begin(), m,
                  result.begin(), m);
      return result;
    }

    template<class Real>
    inline
    Matrix<Real>
    LUMatrix<Real>::U() const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      Matrix<Real> result(m,n);
      TBLAS::geset(TBLAS::RowMajor,m,n,Real(0),Real(0),result.begin(),n);
      TBLAS::trcpy(TBLAS::RowMajor,TBLAS::Upper,TBLAS::NonUnit,
                  m,n,
                  this->_elements.begin(),n,
                  result.begin(),n);
      return result;
    }

    template<class Real>
    inline
    LUMatrix<Real>::operator Matrix<Real> () const
    {
      size_type m=this->number_of_rows();
      size_type n=this->number_of_columns();
      Matrix<Real> result(m,n);
        
      /* No TBLAS routine */
      for(size_type i=0; i!=m; ++i) {
        for(size_type j=0; j!=n; ++j) {
          result(i,j)=0;
          for(size_type k=0; k<=TBLAS::min(i,j); ++k) {
            if(i==k) { result(i,j)+=_elements(k,j); }
            else { result(i,j)+=_elements(i,k)*_elements(k,j); }
          }
        }
      }
      /* Fixme: Apply permutation */
      TLAPACK::laswp(TBLAS::RowMajor,n,result.data().begin(),n,0,n,this->_pivots.begin(),-1);
      return result;
    }

    template<class Real>
    inline
    bool
    LUMatrix<Real>::singular() const
    {
      return this->determinant()==Real(0); 
    }

    template<class Real>
    inline
    Real
    LUMatrix<Real>::determinant() const
    {
      check_square(*this,__PRETTY_FUNCTION__);
      
      size_type n=this->number_of_rows();
      Real result=1;
      for(size_type i=0; i!=n; ++i) { 
        result*=this->_elements(i,i);
      }
      return result;
    }

    template<class Real>
    inline
    Matrix<Real>
    LUMatrix<Real>::inverse() const
    {
      size_type n=this->number_of_rows();
      Matrix<Real> result(n,n);
      TBLAS::geset(TBLAS::RowMajor,n,n,Real(0),Real(1),result.begin(),n);
      TLAPACK::getrs(TBLAS::RowMajor,TBLAS::NoTrans,n,n,this->_elements.begin(),n,this->_pivots.begin(),result.begin(),n);
      return result;
    }

  }
}

#endif /* ARIADNE_LU_MATRIX_H */
