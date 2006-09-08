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
#include <iomanip>
#include <fstream>
#include <string>

#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/interval_matrix.h"
#include "linear_algebra/lu_matrix.h"
#include "linear_algebra/qr_matrix.h"
#include "linear_algebra/svd_matrix.h"

#include <boost/numeric/interval/interval.hpp>
#include <boost/numeric/interval/utility.hpp>

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace std;

bool contains(const Matrix< Interval<Real> >& A, Matrix<Real>& B) {
  return IntervalMatrix<Real>(A).contains(B);
}

int main() {

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
   
    cout << "A=" << A << endl;
    cout << "b=" << b << endl;
    
    //assert(1==13*A(2,1));
   
    x=A.solve(b);
  
    cout << "A.solve(b)= " << x << endl;
    cout << "Ax=" << A*b << endl;
    
    assert(x(0)==1 && x(1)==1 && x(2)==1);
    
  }
  
  Rational Aqptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  Real Arptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  Matrix<Rational> Aq(3,3,Aqptr,3,1);
  Matrix<Real> Ar(3,3,Arptr,3);

  cout << "Testing LU\n";
  {
    Matrix<Rational> A(Aq);
    LUMatrix<Rational> LU(A);
    Matrix<Rational> P=LU.P();
    Matrix<Rational> L=LU.L();
    Matrix<Rational> U=LU.U();
    cout << "A=" << A << "\n";
    cout << "P=" << P << "\n";
    cout << "L=" << L << "\n";
    cout << "U=" << U << "\n";
    cout << "P*L*U=" << P*L*U << "\n";
    cout << "LU=" << Matrix<Rational>(LU) << "\n";
    cout << flush;
    assert(Matrix<Rational>(LU)==A);
  }

  cout << "Testing Real QR\n";
  {
    Matrix<Real> A(Ar);
    QRMatrix<Real> QR(A);
    Matrix<Real> Q=QR.Q();
    Matrix<Real> R=QR.R();
    cout << "A=" << A << "\n";
    cout << "Q=" << Q << "\n";
    cout << "Q^T=" << Q.transpose() << "\n";
    cout << "R=" << R << "\n";
    cout << "Q*Q^T=" << Q*Q.transpose() << "\n";
    cout << setprecision(20);
    cout << "Q*R=" << Q*R << "\n";
    cout << "QR =" << Matrix<Real>(QR) << "\n";
    cout << "A  =" << A << "\n";
    cout << flush;
    //assert(Matrix<Real>(QR)==A);
  }

  cout << "Testing Interval<Real> QR\n";
  {
    Matrix<Real> I=Matrix<Real>::identity(3);
    Matrix<Real> A(Ar);
    QRMatrix< Interval<Real> > QR(A);
    Matrix< Interval<Real> > Q=QR.Q();
    Matrix< Interval<Real> > R=QR.R();
    Matrix< Interval<Real> > B=QR;
    cout << "A=" << A << "\n";
    cout << "Q=" << Q << "\n";
    cout << "Q^T=" << Q.transpose() << "\n";
    cout << "R=" << R << "\n";
    cout << "Q*Q^T=" << Q*Q.transpose() << "\n";
    cout << setprecision(20);
    cout << "Q*R=" << Q*R << "\n";
    cout << "QR =" << B << "\n";
    cout << "A  =" << A << "\n";
    cout << flush;
    assert(contains(Q*Q.transpose(),I));
    assert(contains(B,A));
  }

  /*
  cout << "Testing QR\n";
  {
    QRMatrix<Rational> QR(A);
    Matrix<Rational> Q=QR.Q();
    Matrix<Rational> R=QR.R();
    cout << "A=" << A << "\n";
    cout << "Q=" << Q << "\n";
    cout << "R=" << R << "\n";
    cout << "Q*Q^T=" << Q*Q.transpose() << "\n";
    cout << "Q*R=" << Q*R << "\n";
    cout << "QR=" << Matrix<Rational>(QR) << "\n";
    cout << flush;
    assert(Matrix<Rational>(QR)==A);
  }

  cout << "Testing SVD\n";
  {
    SVDMatrix<Rational> SVD(A);
    Vector<Rational> S=SVD.S();
    Matrix<Rational> U=SVD.U();
    Matrix<Rational> V=SVD.V();
    Matrix<Rational> D=SVD.D();
    cout << "A=" << A << "\n";
    cout << "S=" << S << "\n";
    cout << "D=" << D << "\n";
    cout << "U=" << U << "\n";
    cout << "V=" << V << "\n";
    cout << "U*U^T=" << U*U.transpose() << "\n";
    cout << "V*V^T=" << V*V.transpose() << "\n";
    cout << "U*D*V^T" << U*D*V.transpose() << "\n";
    cout << "SVD=" << Matrix<Rational>(SVD) << "\n";
    cout << flush;
    assert(Matrix<Rational>(SVD)==A);
  }
  */

  cerr << "INCOMPLETE ";

  return 0;
}
