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

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/lu_matrix.h"
#include "linear_algebra/qr_matrix.h"
#include "linear_algebra/svd_matrix.h"

#include <boost/numeric/interval/interval.hpp>
#include <boost/numeric/interval/utility.hpp>

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace std;

template<typename R> int test_linear_algebra();
template<> int test_linear_algebra<Rational>();

int main() {
  test_linear_algebra<Rational>();
  test_linear_algebra<Float64>();
  test_linear_algebra<MPFloat>();
}

template<> 
int 
test_linear_algebra<Rational>()
{
  Matrix< Rational > A(3,3);
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
  
  Rational Aqptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  Aq=Matrix<Rational>(3,3,Aqptr,3,1);

  cout << "Testing LU\n";
  A=Matrix<Rational>(Aq);
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

  return 0;
}

template<typename Rl>
int
test_linear_algebra()
{
  typedef typename Numeric::numerical_traits<Rl>::arithmetic_type F;
  Rl Arptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  Matrix<Rl> Ar(3,3,Arptr,3);

  cout << "Testing real QR\n";
  Matrix<Rl> A(Ar);
  QRMatrix<F> QR(A);
  Matrix<F> Q=QR.Q();
  Matrix<F> R=QR.R();
  cout << "A=" << A << "\n";
  cout << "Q=" << Q << "\n";
  cout << "Q^T=" << Q.transpose() << "\n";
  cout << "R=" << R << "\n";
  cout << "Q*Q^T=" << Q*Q.transpose() << "\n";
  cout << setprecision(20);
  cout << "Q*R=" << Q*R << "\n";
  cout << "QR =" << Matrix<F>(QR) << "\n";
  cout << "A  =" << A << "\n";
  cout << flush;
  //assert(Matrix<Rl>(QR)==A);

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
