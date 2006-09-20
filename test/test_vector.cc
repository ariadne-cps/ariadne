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

#include "real_typedef.h"
#include "numeric/numerical_types.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/vector.tpl"

#undef DEBUG

using namespace std;
using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;

int main() {

  
  
  int n=3;
  Real vptr[3]={-4.0,3.0,1.0};
  
  Real x=1.5;
  Vector<Real> v1(n,vptr);
  Vector<Real> v2("[2.375,4.25,-1.25]");
  cout << "v1=" << v1 << "  v2=" << v2 << endl;
  cout << "v1.norm()=" << v1.norm() << "  v2.norm()=" << v2.norm() << endl;
  assert(v1.norm()==4.0);
  assert(v2.norm()==4.25);

  Vector<Real> v3(v1+v2);
  cout << "v3=" << v3 << endl;
  
  Vector<Real> v4(v3);
  cout << "v4=" << v4 << endl;
  Vector<Real> v0;
  cout << "v0=" << v0 << endl;
  Vector<Real> v5(1);
  cout << "v5=" << v5 << endl;
  v5=v3;
  cout << "v5=" << v5 << endl;
  cout << endl;
  
  v2=Vector<Real>("[0.25,-1.5]");
  v3=Vector<Real>("[-0.5,2.25]");
  v1=-v2;
  cout << v1 << " = -" << v2 << endl;
  v1=v2+v3;
  cout << v1 << " = " << v2 << " + " << v3 << endl;
  v1=v2-v3;
  cout << v1 << " = " << v2 << " - " << v3 << endl;
  v1=x*v3;
  cout << v1 << " = " << x << " * " << v3 << endl;
  v1=v2*x;
  cout << v1 << " = " << v2 << " * " << x << endl;
  v1=v2/x;
  cout << v1 << " = " << v2 << " / " << x << endl;
  
  Vector< Interval<Real> > iv1("[[0.99,1.01],[2.25,2.375],[4.0,4.375],[-0.02,0.01]]");
  cout << "iv1=" << iv1 << endl;
  cout << "iv1.norm()=" << iv1.norm() << endl;
  cout << "iv1.upper_norm()=" << iv1.upper_norm() << endl;

  Vector< Interval<Real> > iv2("[[-1,1],[-1,1]]");
  cout << "iv2=" << iv2 << endl;
  Vector< Interval<Real> > iv3(3);
  cout << "iv3=" << iv3 << endl;
  iv3=Vector< Interval<Real> >("[[4.25,4.25],[2.375,2.375]]");
  cout << "iv3=" << iv3 << endl;
  Interval<Real> ix=Interval<Real>(-2,1);
 
  Vector< Interval<Real> > iv0;
  cout << "iv0=" << iv0 << endl;
  iv1=iv0;
  //cout << "iv1=" << iv1 << endl;
  iv1=iv2;
  cout << "iv1=" << iv1 << endl;
  cout << endl;
  
  //Vector< Interval<Real> > ivv=Interval<Real>(x)*iv1;
  //Vector<Real> vv=x*v1;
  
  boost::numeric::ublas::vector<Real> bv1(3);
  boost::numeric::ublas::vector<Real> bv2=x*bv1;
  
  boost::numeric::ublas::vector< Interval<Real> > biv1(3);
  boost::numeric::ublas::vector< Interval<Real> > biv2=Interval<Real>(x)*biv1;
  
  cout << "Tested ublas" << endl;
  Vector<Real> tv1(3);
  Vector<Real> tv2=(x*x)*tv1;
  cout << "Tested Vector" << endl;
  Vector< Interval<Real> > tiv1(3);
  tiv1=ix*tiv1;
  cout << "Tested IntervalVector" << endl;
  Vector< Interval<Real> > tiv2=ix*tiv1;
  cout << "Tested IntervalVector" << endl;
  Vector< Interval<Real> > tiv3=Interval<Real>(x)*tiv1;
  cout << "Tested IntervalVector" << endl;
 
  
  iv1=iv2+iv3;
  cout << iv1 << " = " << iv2 << " + " << iv3 << endl;
  iv1=iv2-iv3;
  cout << iv1 << " = " << iv2 << " - " << iv3 << endl;
  iv1=ix*iv3;
  cout << iv1 << " = " << ix << " * " << iv3 << endl;
  iv1=iv2*ix;
  cout << iv1 << " = " << iv2 << " * " << ix << endl;
  ix=Interval<Real>(1,2);
  iv1=iv2/ix;
  cout << iv1 << " = " << iv2 << " / " << ix << endl;
  cout << endl;
   
  v1=Vector<Real>("[-1.25,0.75]");
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
