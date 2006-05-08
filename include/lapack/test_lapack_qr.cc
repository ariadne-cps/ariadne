/*
 * Copyright (C) 2006 Pieter Collins <Pieter.Collins@cwi.nl>
 */

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <blas/blas_output.hpp>

#include "lapack.hpp"

#include "geqrf.hpp"
#include "orgqr.hpp"

#include "larf.hpp"
#include "larfg.hpp"

using std::cout;
using std::cerr;
using std::endl;

using BLAS::matrix;
using BLAS::vector;
using BLAS::RowMajor;
using BLAS::Trans;
using BLAS::NoTrans;
using BLAS::Left;
using BLAS::Right;
using BLAS::NonUnit;
using BLAS::Unit;
using BLAS::Lower;
using BLAS::Upper;

int main() {
  const double A[9]={2.0,-0.4,1.0,  2.0,1.0,0.0,  1.0,0.0,1.0};
  //const double A[9]={4.0,3.0,0.0,  3.0,1.0,0.0,  0.0,0.0,1.0};
  //const double A[9]={3.0,1.0,1.0,  0.0,2.0,1.0,  0.0,0.0,1.0};
  //const double A[9]={1.0,0.0,1.0,  1.0,2.0,0.0,  3.0,1.0,1.0};
  //const double A[9]={1.0,2.0,3.0,  4.0,5.0,6.0,  7.0,8.0,9.0};
  const double I[9]={1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0};

  double QxR[9];
  double QTxA[9];
  double R[9];
  double QR[9];
  double Q[9];
  double QT[9];
  double J[9];
  double Z[9];
  int piv[3];
  double E[9];

  int info=0;

  BLAS::copy(9,A,1,QR,1);
  BLAS::copy(9,I,1,J,1);
  
  double work[9];
  double tau[3];
  
  cout << matrix(3,3,QR) << endl;
  LAPACK::geqrf(RowMajor, 3,3,QR,3,tau,work);
  cout << "QR=\n" << matrix(3,3,QR) << endl;
  cout << vector(3,tau) << endl << endl;
  BLAS::copy(9,QR,1,Q,1);
  cout << "QR=\n" << matrix(3,3,QR) << endl;
  LAPACK::orgqr(RowMajor, 3,3,3,Q,3,tau,work);
  cout << "QR=\n" << matrix(3,3,QR) << endl;
  cout << matrix(3,3,Q) << endl;
  BLAS::copy(9,I,1,R,1);
  BLAS::trmm(RowMajor,Left,Upper,NoTrans,NonUnit,3,3,1.0,QR,3,R,3);
  cout << "R=\n" << matrix(3,3,R) << endl;
  BLAS::gemm(RowMajor,NoTrans,Trans,3,3,3,1.0,I,3,Q,3,0.0,QT,3);
  cout << matrix(3,3,QT) << endl;
  BLAS::gemm(RowMajor,NoTrans,NoTrans,3,3,3,1.0,Q,3,QT,3,0.0,Z,3);
  cout << matrix(3,3,Z) << endl;
  BLAS::copy(9,Z,1,E,1);
  BLAS::axpy(9,-1.0,I,1,E,1);
  cout << matrix(3,3,E) << endl;
  cout << BLAS::amax(9,E,1) << endl;
  BLAS::gemm(RowMajor,NoTrans,NoTrans,3,3,3,1.0,Q,3,R,3,0.0,QxR,3);
  cout << "QxR=\n" << matrix(3,3,QxR) << endl;
  BLAS::copy(9,A,1,E,1);
  BLAS::axpy(9,-1.0,QxR,1,E,1);
  cout << BLAS::amax(9,E,1) << endl;
  BLAS::gemm(RowMajor,Trans,NoTrans,3,3,3,1.0,Q,3,A,3,0.0,QTxA,3);
  cout << "Q'xA=\n" << matrix(3,3,QTxA) << endl;
  
  return 0;
}
