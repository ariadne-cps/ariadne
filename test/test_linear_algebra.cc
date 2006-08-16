/***************************************************************************
 *            test_linear_algebra.cc
 *
 *  31 Jan 2006
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

#include <iostream>
#include <fstream>
#include <string>

#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/lu_matrix.h"
#include "linear_algebra/qr_matrix.h"
#include "linear_algebra/svd_matrix.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace std;

int main() {
  cout << "test_linear_algebra: " << flush;
  ofstream clog("test_linear_algebra.log");

  {
    
    Matrix< Rational > A(3,3), LU(3,3);
    Matrix< Rational > Aq(3,3), QT, R;
    Vector< Rational > b(3), x(3);
    Vector< dimension_type > pivot(3);
    
    A(0,0)=Rational(7,3);A(0,1)=4;A(0,2)=-1;
    A(1,0)=-13;A(1,1)=Rational(10,9);A(1,2)=6;
    A(2,0)=0;A(2,1)=Rational(1,13);A(2,2)=9;
   
    b(0)=Rational(16,3);
    b(1)=-Rational(53,9);
    b(2)=Rational(118,13);
   
    test_assert(1==13*A(2,1),"mutiplication");
   
    x=A.solve(b);
  
    test_assert(x(0)==1 && x(1)==1 && x(2)==1,"lu_solve");
  
  }
  
  Real Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  //Matrix<Real> A(3,3,Aptr,3);
  Matrix<Real> A(3,3,Aptr,3,1);

  clog << "Testing LU\n";
  {
    LUMatrix<Real> LU(A);
    Matrix<Real> P=LU.P();
    Matrix<Real> L=LU.L();
    Matrix<Real> U=LU.U();
    clog << "A=" << A << "\n";
    clog << "P=" << P << "\n";
    clog << "L=" << L << "\n";
    clog << "U=" << U << "\n";
    clog << "P*L*U=" << P*L*U << "\n";
    clog << "LU=" << Matrix<Real>(LU) << "\n";
  }

  clog << "Testing QR\n";
  {
    QRMatrix<Real> QR(A);
    Matrix<Real> Q=QR.Q();
    Matrix<Real> R=QR.R();
    clog << "A=" << A << "\n";
    clog << "Q=" << Q << "\n";
    clog << "R=" << R << "\n";
    clog << "Q*Q^T=" << Q*Q.transpose() << "\n";
    clog << "Q*R=" << Q*R << "\n";
    clog << "QR=" << Matrix<Real>(QR) << "\n";
  }

  clog << "Testing SVD\n";
  {
    SVDMatrix<Real> SVD(A);
    Vector<Real> S=SVD.S();
    Matrix<Real> U=SVD.U();
    Matrix<Real> V=SVD.V();
    Matrix<Real> D=SVD.D();
    clog << "A=" << A << "\n";
    clog << "S=" << S << "\n";
    clog << "D=" << D << "\n";
    clog << "U=" << U << "\n";
    clog << "V=" << V << "\n";
    clog << "U*U^T=" << U*U.transpose() << "\n";
    clog << "V*V^T=" << V*V.transpose() << "\n";
    clog << "U*D*V^T" << U*D*V.transpose() << "\n";
    clog << "SVD=" << Matrix<Real>(SVD) << "\n";
  }
    
  clog.close();
  cout << "INCOMPLETE\n";

    return 0;
}
