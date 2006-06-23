/***************************************************************************
 *            test_linear_algebra.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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

#define NO_CBLAS

#include <iostream>
#include <fstream>

#include <tblas/tblas.hpp>
#include <tblas/output.hpp>
#include <tlapack/tlapack.hpp>

#include <tblas/gemv.hpp>

#include "declarations.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/vector.tpl"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix.tpl"
#include "linear_algebra/lu_matrix.h"
#include "linear_algebra/qr_matrix.h"
#include "linear_algebra/svd_matrix.h"

#undef DEBUG

using namespace std;
using namespace Ariadne::LinearAlgebra;

int main() {
  cout << "test_matrix: " << flush;
  ofstream clog("test_matrix.log");
  
  int m=3;
  int n=3;
  double Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  
  Matrix<double> A(m,n,Aptr,n);

  clog << "Testing LU\n";
  {
    LUMatrix<double> LU(A);
    Matrix<double> P=LU.P();
    Matrix<double> L=LU.L();
    Matrix<double> U=LU.U();
    clog << "A=" << A << "\n";
    clog << "P=" << P << "\n";
    clog << "L=" << L << "\n";
    clog << "U=" << U << "\n";
    clog << "P*L*U=" << P*L*U << "\n";
    clog << "LU=" << Matrix<double>(LU) << "\n";
  }

  clog << "Testing QR\n";
  {
    QRMatrix<double> QR(A);
    Matrix<double> Q=QR.Q();
    Matrix<double> R=QR.R();
    clog << "A=" << A << "\n";
    clog << "Q=" << Q << "\n";
    clog << "R=" << R << "\n";
    clog << "Q*Q^T=" << Q*Q.transpose() << "\n";
    clog << "Q*R=" << Q*R << "\n";
    clog << "QR=" << Matrix<double>(QR) << "\n";
  }

  clog << "Testing SVD\n";
  {
    SVDMatrix<double> SVD(A);
    Vector<double> S=SVD.S();
    Matrix<double> U=SVD.U();
    Matrix<double> V=SVD.V();
    Matrix<double> D=SVD.D();
    clog << "A=" << A << "\n";
    clog << "S=" << S << "\n";
    clog << "D=" << D << "\n";
    clog << "U=" << U << "\n";
    clog << "V=" << V << "\n";
    clog << "U*U^T=" << U*U.transpose() << "\n";
    clog << "V*V^T=" << V*V.transpose() << "\n";
    clog << "U*D*V^T" << U*D*V.transpose() << "\n";
    clog << "SVD=" << Matrix<double>(SVD) << "\n";
  }

  clog.close();
  cout << "PASS\n";

  return 0;
}
