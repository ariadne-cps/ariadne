/***************************************************************************
 *            test_vector.cc
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

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/vector.tpl"

#undef DEBUG

using namespace std;
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;

template<typename R> int test_vector();

int main() {
  test_vector<MPFloat>();
  test_vector<Float64>();
  //test_vector<Rational>();
  //test_vector< Interval<MPFloat> >();
  //test_vector< Interval<Rational> >();

  return 0;
}  

template<typename R>
int 
test_vector()
{
  std::cout << "\ntest_vector<" << name<R>() << ">()\n";
  
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  int n=3;
  R vptr[3]={-4.0,3.0,1.0};
  R x=1.5;

  Vector<R> v0;
  cout << "v0=" << v0 << endl;
  Vector<R> v1(n,vptr);
  cout << "v1=" << v1 << endl;
  Vector<R> v2("[2.375,4.25,-1.25]");
  cout << "v2=" << v2 << endl;
  cout << "v1.norm()=" << v1.norm() << "  v2.norm()=" << v2.norm() << endl;
  assert(v1.norm()==4);
  assert(v2.norm()==4.25);

  Vector<R> v3(1);
  cout << "v3=" << v3 << endl;
  Vector<R> v4=v2;
  cout << "v4=" << v4 << endl;
  cout << endl;

  Vector<F> vf0;
  v1=Vector<R>("[0.25,-1.5]");
  v2=Vector<R>("[-0.5,2.25]");
  vf0=-v1;
  cout << vf0 << " = -" << v1 << endl;
  vf0=Vector<F>(v1)+v2;
  cout << vf0 << " = " << v1 << " + " << v2 << endl;
  vf0=Vector<F>(v1)-v2;
  cout << vf0 << " = " << v1 << " - " << v2 << endl;
  vf0=x*Vector<F>(v2);
  cout << vf0 << " = " << x << " * " << v2 << endl;
  vf0=Vector<F>(v1)*x;
  cout << vf0 << " = " << v1 << " * " << x << endl;
  vf0=Vector<F>(v1)/x;
  cout << vf0 << " = " << v1 << " / " << x << endl;
  cout << endl;
  
  Vector< Interval<R> > iv1("[[0.99,1.01],[2.25,2.375],[4.0,4.375],[-0.02,0.01]]");
  cout << "iv1=" << iv1 << endl;
  cout << "iv1.norm()=" << iv1.norm() << endl;
  cout << "iv1.norm().upper()=" << iv1.norm().upper() << endl;

  Vector< Interval<R> > iv2("[[-1,1],[-1,1]]");
  cout << "iv2=" << iv2 << endl;
  Vector< Interval<R> > iv3(3);
  cout << "iv3=" << iv3 << endl;
  iv3=Vector< Interval<R> >("[[4.25,4.25],[2.375,2.375]]");
  cout << "iv3=" << iv3 << endl;
  Interval<R> ix=Interval<R>(-2,1);
 
  Vector< Interval<R> > iv0;
  cout << "iv0=" << iv0 << endl;
  iv1=iv0;
  cout << "iv1=" << iv1 << endl;
  iv1=iv2;
  cout << "iv1=" << iv1 << endl;
  cout << endl;

  Interval<R> ix2=iv2[0];
  Interval<R> ix3=iv3[0];
  Interval<R> ix1=ix2+ix3;
  ix1=ix2+ix3;
  
  cout << "iv2=" << iv2 << ", iv3=" << iv3 << endl;
  iv1=iv2+iv3;
  cout << iv1 << " = " << iv2 << " + " << iv3 << endl;
  iv1=iv2-iv3;
  cout << iv1 << " = " << iv2 << " - " << iv3 << endl;
  iv1=ix*iv3;
  cout << iv1 << " = " << ix << " * " << iv3 << endl;
  iv1=iv2*ix;
  cout << iv1 << " = " << iv2 << " * " << ix << endl;
  ix=Interval<R>(1,2);
  iv1=iv2/ix;
  cout << iv1 << " = " << iv2 << " / " << ix << endl;
  cout << endl;
   
  v1=Vector<R>("[-1.25,0.75]");
  iv0=iv1+v1;
  cout << iv0 << " = " << iv1 << " + " << v1 << endl;
  iv0=v1+iv1;
  cout << iv0 << " = " << v1 << " + " << iv1 << endl;
  iv0=iv1-v1;
  cout << iv0 << " = " << iv1 << " - " << v1 << endl;
  iv0=v1-iv1;
  cout << iv0 << " = " << v1 << " - " << iv1 << endl;
  iv0=x*iv1;
  cout << iv0 << " = " << x << " * " << iv1 << endl;
  iv0=ix*v1;
  cout << iv0 << " = " << ix << " * " << v1 << endl;
  iv0=iv1*x;
  cout << iv0 << " = " << iv1 << " * " << x << endl;
  iv0=v1*ix;
  cout << iv0 << " = " << v1 << " * " << ix << endl;
  iv0=iv1/x;
  cout << iv0 << " = " << iv1 << " / " << x << endl;
  iv0=v1/ix;
  cout << iv0 << " = " << v1 << " / " << ix << endl;

  return 0;
}
