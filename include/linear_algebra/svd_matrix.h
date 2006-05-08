/***************************************************************************
 *            svd_matrix.h
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
 
/*! \file svd_matrix.h
 *  \brief Singular value decomposition.
 */

#ifndef _ARIADNE_SVD_MATRIX_H
#define _ARIADNE_SVD_MATRIX_H

#include <blas/copy.hpp>
#include <blas/geset.hpp>
#include <lapack/gesvd.hpp>

#include "../declarations.h"
#include "../base/array.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix stored in SVD product form. */
    template<typename Real>
    class SVDMatrix {
     public:
      SVDMatrix(const Matrix<Real>& A);
      
      size_type number_of_rows() const { return _U.number_of_rows(); }
      size_type number_of_columns() const { return _Vt.number_of_rows(); }
      Real operator() (const size_type& i, const size_type& j) const;

      const Vector<Real>& S() const { return _S; }
      const Matrix<Real>& U() const { return _U; }
      const Matrix<Real> V() const { return _Vt.transpose(); }
      const Matrix<Real>& Vt() const { return _Vt; }
      Matrix<Real> D() const;
      
      operator Matrix<Real> () const;
     private:
      Vector<Real> _S;
      Matrix<Real> _U;
      Matrix<Real> _Vt;
    };
    
    template <typename Real>
    inline 
    SVDMatrix<Real>::SVDMatrix(const Matrix<Real>& A) 
      : _S(BLAS::min(A.number_of_rows(),A.number_of_columns())),
        _U(A.number_of_rows(),A.number_of_rows()),
        _Vt(A.number_of_columns(),A.number_of_columns())
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      array<Real> work(A.begin(),A.begin()+m*n);
      LAPACK::gesvd(BLAS::RowMajor,m,n,work.begin(),n,_S.data().begin(),_U.data().begin(),m,_Vt.data().begin(),n);
    }

    template <typename Real>
    inline
    Matrix<Real>
    SVDMatrix<Real>::D() const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      int nsv=BLAS::min(m,n);
      Matrix<Real> result(m,n);
      BLAS::geset(BLAS::RowMajor,m,n,Real(0),Real(0),result.data().begin(),n);
      BLAS::copy(nsv,_S.data().begin(),1,result.data().begin(),n+1);
      return result;
    }

    template <typename Real>
    inline
    Real
    SVDMatrix<Real>::operator() (const size_type& i, const size_type& j) const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();

      Real result=0;
      for(int k=0; k!=BLAS::min(m,n); ++k) { 
        result += _U(i,k) * _S(k) * _Vt(j,k); 
      } 
      return result; 
    }


    template <typename Real>
    inline
    SVDMatrix<Real>::operator Matrix<Real> () const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();

      Matrix<Real> result(m,n);
      for(int i=0; i!=m; ++i) {
        for(int j=0; j!=n; ++j) {
          result(i,j)=0;
          for(int k=0; k!=BLAS::min(m,n); ++k) { 
            result(i,j) += _U(i,k) * _S(k) * _Vt(j,k); 
          } 
        }
      }
      return result; 
    }

  }
}

#endif /* _ARIADNE_SVD_MATRIX_H */
