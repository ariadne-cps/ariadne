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

using namespace std;
using namespace Ariadne::LinearAlgebra;

int main() {
  int m=3;
  int n=3;
  double Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  
  Matrix<double> A(m,n,Aptr,n);

  std::cout << "Testing LU\n";
  {
    LUMatrix<double> LU(A);
    Matrix<double> P=LU.P();
    Matrix<double> L=LU.L();
    Matrix<double> U=LU.U();
    std::cout << "A=" << A << "\n";
    std::cout << "P=" << P << "\n";
    std::cout << "L=" << L << "\n";
    std::cout << "U=" << U << "\n";
    std::cout << "P*L*U=" << P*L*U << "\n";
    std::cout << "LU=" << Matrix<double>(LU) << "\n";
  }

  std::cout << "Testing QR\n";
  {
    QRMatrix<double> QR(A);
    Matrix<double> Q=QR.Q();
    Matrix<double> R=QR.R();
    std::cout << "A=" << A << "\n";
    std::cout << "Q=" << Q << "\n";
    std::cout << "R=" << R << "\n";
    std::cout << "Q*Q^T=" << Q*Q.transpose() << "\n";
    std::cout << "Q*R=" << Q*R << "\n";
    std::cout << "QR=" << Matrix<double>(QR) << "\n";
  }

  std::cout << "Testing SVD\n";
  {
    SVDMatrix<double> SVD(A);
    Vector<double> S=SVD.S();
    Matrix<double> U=SVD.U();
    Matrix<double> V=SVD.V();
    Matrix<double> D=SVD.D();
    std::cout << "A=" << A << "\n";
    std::cout << "S=" << S << "\n";
    std::cout << "D=" << D << "\n";
    std::cout << "U=" << U << "\n";
    std::cout << "V=" << V << "\n";
    std::cout << "U*U^T=" << U*U.transpose() << "\n";
    std::cout << "V*V^T=" << V*V.transpose() << "\n";
    std::cout << "U*D*V^T" << U*D*V.transpose() << "\n";
    std::cout << "SVD=" << Matrix<double>(SVD) << "\n";
  }

  std::cout << "Passed\n";

  return 0;
}
