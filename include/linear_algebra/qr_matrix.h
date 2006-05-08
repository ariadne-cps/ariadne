/***************************************************************************
 *            qr_matrix.h
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
 
/*! \file qr_matrix.h
 *  \brief QR factorisation and factorised matrices.
 */

#ifndef _ARIADNE_QR_MATRIX_H
#define _ARIADNE_QR_MATRIX_H

#include <blas/gecpy.hpp>
#include <blas/trset.hpp>
#include <blas/gemm.hpp>
#include <lapack/larf.hpp>
#include <lapack/larfg.hpp>
#include <lapack/geqrf.hpp>
#include <lapack/orgqr.hpp>

#include "../declarations.h"
#include "../base/array.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \brief A matrix stored in QR product form. */
    template<typename Real>
    class QRMatrix {
     public:
      QRMatrix(const Matrix<Real>& A);
      
      size_type number_of_rows() const { return _R.number_of_rows(); }
      size_type number_of_columns() const { return _R.number_of_columns(); }
      Real operator() (const size_type& i, const size_type& j) const { 
        Real result=0; for(int k=0; k!=j; ++k) { result(i,j)+= _Q(i,k) * _R(k,j); } return result; }

      const Matrix<Real>& Q() const { return _Q; }
      const Matrix<Real>& R() const { return _R; }
      
      operator Matrix<Real> () const;
     private:
      Matrix<Real> _Q;
      Matrix<Real> _R;
    };
    
    template <typename Real>
    inline 
    QRMatrix<Real>::QRMatrix(const Matrix<Real>& A) 
      : _Q(A.number_of_rows(),A.number_of_rows()), _R(A)
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      array<Real> tau(m);
      array<Real> work(BLAS::max(m,n));
      LAPACK::geqrf(BLAS::RowMajor,m,n,this->_R.data().begin(),n,tau.begin(),work.begin());
      BLAS::gecpy(BLAS::RowMajor,m,m,_R.begin(),n,_Q.begin(),m);
      //BLAS::copy(m*m,_R.data().begin(),1,_Q.data().begin(),1);
      LAPACK::orgqr(BLAS::RowMajor,m,m,m,this->_Q.data().begin(),m,tau.begin(),work.begin());
      BLAS::trset(BLAS::RowMajor,BLAS::Lower,BLAS::Unit,m,n,Real(0),_R.begin(),n);
    }

    template <typename Real>
    inline
    QRMatrix<Real>::operator Matrix<Real> () const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();

      Matrix<Real> result(m,n);
      BLAS::gemm(BLAS::RowMajor,BLAS::NoTrans,BLAS::NoTrans,m,n,m,Real(1),_Q.begin(),m,_R.begin(),n,Real(0),result.begin(),n);
      return result;
    }

  }
}

#endif /* _ARIADNE_QR_MATRIX_H */
