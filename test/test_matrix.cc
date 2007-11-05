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

#include "numeric/rational.h"
#include "test/test_float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

// The following includes are not necessary,
// but allow testing without rebuilding the entire library.
#include "linear_algebra/vector.code.h"
#include "linear_algebra/matrix.code.h"

using namespace std;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;

template<class R> int test_matrix();

int main() {
  test_matrix<Flt>();
  //test_matrix<Rational>();
  return 0;
}

template<class R> 
int 
test_matrix()
{
  typedef typename Numeric::traits<R>::arithmetic_type F;

  R x=2.25;
  Interval<R> ix(1.5,2.25);
  R Aptr[9]={-1.0,3.0,1.0, -1.0,1.0,2.0, 2.0,1.0,1.0};
  Interval<R> iAptr[4]={-1.0,3.0, -1.0,1.0};
  
  Matrix<R> A0;
  cout << "A0=" << A0 << endl;
  Matrix<R> A1(3,2);
  cout << "A1=" << A1 << endl;
  Matrix<R> A2(3,3,Aptr,3,1);
  cout << "A2=" << A2 << endl;
  Matrix<R> A3("[-1.0,3.0,1.0; -1.0,1.0,2.0; 2.0,1.0,1.0]");
  cout << "A3=" << A3 << endl;

  cout << "A1(0,0)= " << A1(0,0) << endl;
  cout << "A2(0,0)= " << A2(0,0) << endl;
  assert(A2==A3);
  
  A1[0][0]=1.0;
  cout << "A1[0][0]= " << A1[0][0] << endl;
  assert(A1[0][0]==1.0);
  cout << "A2[0][0]= " << A2[0][0] << endl;
  assert(A2==A3);
  
  A1(0,0)=2.0;
  cout << "A1(0,0)= " << A1(0,0) << endl;
  assert(A1(0,0)==2.0);
  assert(A1[0][0]==A1(0,0));
  
  A0=Matrix<R>::zero(2,3);
  cout << "A0= " << A0 << endl;
  A1=Matrix<R>::identity(4);
  cout << "A1= " << A1 << endl;
  cout << endl;
  
  A1=Matrix<R>("[2,1,1;1,1,-1;1,-1,3]");
  Matrix<F> fA0;

  fA0=-A1;
  cout << A0 << " = -" << A1 << endl;
  fA0=A1+A2;
  cout << A0 << " = " << A1 << " + " << A2 << endl;
  fA0=A1-A2;
  cout << A0 << " = " << A1 << " - " << A2 << endl;
  fA0=x*A2;
  cout << A0 << " = " << x << " * " << A2 << endl;
  fA0=A1*x;
  cout << A0 << " = " << A1 << " * " << x << endl;
  fA0=Matrix<F>(A1)/x;
  cout << A0 << " = " << A1 << " / " << x << endl;
  fA0=Matrix<F>(A1)*A2;
  cout << A0 << " = " << A1 << " * " << A2 << endl;
  cout << endl;

  
  Matrix< Interval<R> > iA0;
  cout << "iA1= " << iA0 << endl;
  Matrix< Interval<R> > iA1(3,2);
  cout << "iA2= " << iA1 << endl;
  Matrix< Interval<R> > iA2(2,2,iAptr,3,1);
  cout << "iA3= " << iA2 << endl;
  Matrix< Interval<R> > iA3(A1);
  cout << "iA4= " << iA3 << endl;
  Matrix< Interval<R> > iA4("[[1.875,2.125],[0.75,1.25];[-1.25,-0.875],[0.5,1.75]]");
  cout << "iA5= " << iA4 << endl;
  
  iA0=Matrix< Interval<R> >::zero(2,3);
  cout << "iA0= " << iA0 << endl;
  iA1=Matrix< Interval<R> >::identity(4);
  cout << "iA1= " << iA1 << endl;

  cout << "norm(iA4)=" << norm(iA4) << endl;
  cout << "norm(iA4).upper()=" << norm(iA4).upper() << endl;
  
  return 0;
}
