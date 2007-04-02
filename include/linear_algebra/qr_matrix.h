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

#ifndef ARIADNE_QR_MATRIX_H
#define ARIADNE_QR_MATRIX_H

#include <tblas/gecpy.hpp>
#include <tblas/trset.hpp>
#include <tblas/gemm.hpp>
#include <tlapack/larf.hpp>
#include <tlapack/larfg.hpp>
#include <tlapack/geqrf.hpp>
#include <tlapack/orgqr.hpp>

#include "../base/array.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \ingroup LinearAlgebra
     *  \brief A matrix stored in QR product form. */
    template<class Real>
    class QRMatrix {
     public:
      /*! \brief Construct from an ordinary matrix. */
      QRMatrix(const Matrix<Real>& A);
      
      /*! \brief The number of rows. */
      size_type number_of_rows() const { return _r.number_of_rows(); }
      /*! \brief The number of columns. */
      size_type number_of_columns() const { return _r.number_of_columns(); }
      /*! \brief The \a i,\a j th element. */
      Real operator() (const size_type& i, const size_type& j) const { 
        Real result=0; for(int k=0; k!=j; ++k) { result(i,j)+= _q(i,k) * _r(k,j); } return result; }

      /*! \brief The orthogonal factor. */
      const Matrix<Real>& Q() const { return _q; }
      /*! \brief The upper-triangular factor. */
      const Matrix<Real>& R() const { return _r; }
      
      /*! \brief Convert to an ordinary matrix. */
      operator Matrix<Real> () const;
     private:
      Matrix<Real> _q;
      Matrix<Real> _r;
    };
    
    template<class Real>
    inline 
    QRMatrix<Real>::QRMatrix(const Matrix<Real>& A) 
      : _q(A.number_of_rows(),A.number_of_rows()), _r(A)
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      array<Real> tau(m);
      array<Real> work(TBLAS::max(m,n));
      TLAPACK::geqrf(TBLAS::RowMajor,m,n,this->_r.data().begin(),n,tau.begin(),work.begin());
      TBLAS::gecpy(TBLAS::RowMajor,m,m,_r.begin(),n,_q.begin(),m);
      //TBLAS::copy(m*m,_r.data().begin(),1,_q.data().begin(),1);
      TLAPACK::orgqr(TBLAS::RowMajor,m,m,m,this->_q.data().begin(),m,tau.begin(),work.begin());
      TBLAS::trset(TBLAS::RowMajor,TBLAS::Lower,TBLAS::Unit,m,n,Real(0),_r.begin(),n);
    }

    template<class Real>
    inline
    QRMatrix<Real>::operator Matrix<Real> () const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();

      Matrix<Real> result(m,n);
      TBLAS::gemm(TBLAS::RowMajor,TBLAS::NoTrans,TBLAS::NoTrans,m,n,m,Real(1),_q.begin(),m,_r.begin(),n,Real(0),result.begin(),n);
      return result;
    }

  }
}

#endif /* ARIADNE_QR_MATRIX_H */
