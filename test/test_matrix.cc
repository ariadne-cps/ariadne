/***************************************************************************
 *            test_matrix.cc
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

#include "declarations.h"
#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/vector.tpl"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix.tpl"
#include "linear_algebra/interval_matrix.h"
#include "linear_algebra/interval_matrix.tpl"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;

int main() {
  cout << "test_matrix: " << flush;
  ofstream clog("test_matrix.log");
  
  Real x=2.25;
  Interval<Real> ix(1.5,2.25);
  Real Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  Interval<Real> iAptr[4]={-1.0,3.0, -1.0,1.0};
  
  Matrix<Real> A0;
  clog << "A0=" << A0 << endl;
  Matrix<Real> A1(3,2);
  clog << "A1=" << A1 << endl;
  Matrix<Real> A2(3,3,Aptr,3,1);
  clog << "A2=" << A2 << endl;
  Matrix<Real> A3("[-1.0,3.0,1.0; -1.0,1.0,2.0; 2.0,1.0,1.0]");
  clog << "A3=" << A3 << endl;

  clog << "A1(0,0)= " << A1(0,0) << endl;
  clog << "A2(0,0)= " << A2(0,0) << endl;
  assert(A2==A3);
  
  A0=Matrix<Real>::zero(2,3);
  clog << "A0= " << A0 << endl;
  A1=Matrix<Real>::identity(4);
  clog << "A1= " << A1 << endl;
  clog << endl;
  
  A1=Matrix<Real>("[2,1,1;1,1,-1;1,-1,3]");
  
  A0=-A1;
  clog << A0 << " = -" << A1 << endl;
  A0=A1+A2;
  clog << A0 << " = " << A1 << " + " << A2 << endl;
  A0=A1-A2;
  clog << A0 << " = " << A1 << " - " << A2 << endl;
  A0=x*A2;
  clog << A0 << " = " << x << " * " << A2 << endl;
  A0=A1*x;
  clog << A0 << " = " << A1 << " * " << x << endl;
  A0=A1/x;
  clog << A0 << " = " << A1 << " / " << x << endl;
  A0=A1*A2;
  clog << A0 << " = " << A1 << " * " << A2 << endl;
  clog << endl;

  
  IntervalMatrix<Real> iA0;
  clog << "iA1= " << iA0 << endl;
  IntervalMatrix<Real> iA1(3,2);
  clog << "iA2= " << iA1 << endl;
  IntervalMatrix<Real> iA2(2,2,iAptr,3,1);
  clog << "iA3= " << iA2 << endl;
  IntervalMatrix<Real> iA3(A1);
  clog << "iA4= " << iA3 << endl;
  IntervalMatrix<Real> iA4("[[1.875,2.125],[0.75,1.25];[-1.25,-0.875],[0.5,1.75]]");
  clog << "iA5= " << iA4 << endl;
  
  iA0=IntervalMatrix<Real>::zero(2,3);
  clog << "iA0= " << iA0 << endl;
  iA1=IntervalMatrix<Real>::identity(4);
  clog << "iA1= " << iA1 << endl;

  clog << "iA4.norm()=" << iA4.norm() << endl;
  clog << "iA4.upper_norm()=" << iA4.upper_norm() << endl;
  
  clog.close();
  cout << "PASS\n";

  return 0;
}
