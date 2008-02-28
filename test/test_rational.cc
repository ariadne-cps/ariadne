/***************************************************************************
 *            test_rational.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>

#include "numeric/rational.h"

#include "test.h"

using namespace std;
using namespace Ariadne::Numeric;

class TestRational
{
 public:
  void test();
 private:
  void test_concept();
  void test_constructors();
  void test_literal();
  void test_arithmetic();
  void test_comparison();
};


void
TestRational::test()
{
  void test_constructors();
  void test_literal();
  void test_arithmetic();
  void test_comparison();
}


void
TestRational::test_concept()
{
  int n; uint m; double d; Integer z; Rational q;
  
  // Constructors and assignment
  q=Rational(); q=Rational(n); q=Rational(m); q=Rational(d); q=Rational(z); q=Rational(q);
  q=m; q=n; q=d; q=z; q=q;
  
  // Arithmetic
  q=min(q,q); q=max(q,q); 
  q=abs(q); q=pos(q); q=neg(q);

  q=add(q,q); q=sub(q,q); q=mul(q,q); q=div(q,q); 
  q=pow(q,n); q=pow(q,m);
  
  // Mixed arithmetic
  q=div(z,z); q=div(z,n); q=div(z,m);
  
  // Comparisons
  q==q; q!=q; q<=q; q>=q; q<q; q>q;
  q==n; q!=n; q<=n; q>=n; q<n; q>n;
  n==q; n!=q; n<=q; n>=q; n<q; n>q;
  q==m; q!=m; q<=m; q>=m; q<m; q>m;
  m==q; m!=q; m<=q; m>=q; m<q; m>q;
  q==d; q!=d; q<=d; q>=d; q<d; q>d;
  d==q; d!=q; d<=q; d>=q; d<q; d>q;
  q==z; q!=z; q<=z; q>=z; q<z; q>z;
  z==q; z!=q; z<=q; z>=q; z<q; z>q;
}

void
TestRational::test_constructors()
{
  cout << __PRETTY_FUNCTION__ << endl;

  // Construct from an int
  Rational q1(2);
  ARIADNE_TEST_ASSERT(q1==2);
  // Construct from a double
  Rational q2(1.25);
  ARIADNE_TEST_ASSERT(q2==1.25);
  // Copy constructor
  Rational q3(q2);
  ARIADNE_TEST_ASSERT(q3==q2);
  
  // Assign from an int
  q1=3;
  ARIADNE_TEST_ASSERT(q1==3);
  // Assign from a double
  q2=2.25;
  ARIADNE_TEST_ASSERT(q2==2.25);
  // Copy assignment
  q3=q2;
  ARIADNE_TEST_ASSERT(q3==q2);
}

void
TestRational::test_literal()
{
  cout << "Testing stream input" << endl;

  stringstream ss("42 -23/5 8/2 1.25 -2.25");
  Rational q1,q2,q3,q4,q5;
  ss >> q1;
  cout << "q1 = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q1==42);
  ss >> q2;
  cout << "q2 = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q2==Rational(-23,5));
  ss >> q3;
  cout << "q3 = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(8,2));
  ARIADNE_TEST_ASSERT(Rational(12,3)==Rational(8,2));
  ARIADNE_TEST_ASSERT(q3==Rational(4));
  ss >> q4;
  cout << "q4 = " << q4 << endl;
  ARIADNE_TEST_ASSERT(q4==Rational(5,4));
  ss >> q5;
  cout << "q5 = " << q4 << endl;
  ARIADNE_TEST_ASSERT(q5==Rational(-9,4));
}
  

void
TestRational::test_arithmetic()
{
  cout << "Testing arithmetic" << endl;

  Rational q1(1.25);
  Rational q2(2.25);
  Rational q3;
    
  q3=add(q1,q2);
  cout << q1 << " + " << q2 << " = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(7,2));
  q3=sub(q1,q2);
  cout << q1 << " - " << q2 << " = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(-1,1));
  q3=mul(q1,q2);
  cout << q1 << " * " << q2 << " = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(45,16));
  q3=div(q1,q2);
  cout << q1 << " / " << q2 << " = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(5,9));
  
  q3=med(q1,q2);
  cout << "med( " << q1 << " , " << q2 << ") = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(7,4));
  
  q3=rad(q1,q2);
  cout << "rad( " << q1 << " , " << q2 << ") = " << q3 << endl;
  ARIADNE_TEST_ASSERT(q3==Rational(1,2));
  
  cout << endl;
}


void
TestRational::test_comparison()
{
  // Test comparison of two equal numbers
  Rational q1(1.25); Rational q2(-1.25); Rational q3(-2.25); Rational q4(1.25); // Rational q5(1,0);
  ARIADNE_TEST_ASSERT(!(q1==q2)); 
  ARIADNE_TEST_ASSERT(q1!=q2); 
  ARIADNE_TEST_ASSERT(!(q1<=q2)); 
  ARIADNE_TEST_ASSERT(q1> q2);
  ARIADNE_TEST_ASSERT(q1>=q2); 
  ARIADNE_TEST_ASSERT(!(q1< q2));
  
  // Test comparison of two different numbers
  ARIADNE_TEST_ASSERT(q1==q4); 
  ARIADNE_TEST_ASSERT(!(q1!=q4));
  ARIADNE_TEST_ASSERT(q1<=q4); 
  ARIADNE_TEST_ASSERT(!(q1> q4));
  ARIADNE_TEST_ASSERT(q1>=q4); 
  ARIADNE_TEST_ASSERT(!(q1< q4));
  
  // Test comparison with infinity
  std::cerr << "WARNING: comparison with infinity does not work\n";
  /*
    ARIADNE_TEST_ASSERT(q4<q5); 
    ARIADNE_TEST_ASSERT(q4!=q5);
    ARIADNE_TEST_ASSERT(!(q4>q5)); 
    ARIADNE_TEST_ASSERT(q5==q5);
  */
  
  // Test comparison with in integer
  int i2=1;
  ARIADNE_TEST_ASSERT(!(q1==i2)); 
  ARIADNE_TEST_ASSERT(q1!=i2); 
  ARIADNE_TEST_ASSERT(!(q1<=i2)); 
  ARIADNE_TEST_ASSERT(q1> i2);
  ARIADNE_TEST_ASSERT(q1>=i2); 
  ARIADNE_TEST_ASSERT(!(q1< i2));
  
  int i1=1;
  ARIADNE_TEST_ASSERT(!(i1==q2)); 
  ARIADNE_TEST_ASSERT(i1!=q2); 
  ARIADNE_TEST_ASSERT(!(i1<=q2)); 
  ARIADNE_TEST_ASSERT(i1> q2);
  ARIADNE_TEST_ASSERT(i1>=q2); 
  ARIADNE_TEST_ASSERT(!(i1< q2));
  
  // Test comparison with a double
  double x2=1.0;
  ARIADNE_TEST_ASSERT(!(q1==x2)); 
  ARIADNE_TEST_ASSERT(q1!=x2); 
  ARIADNE_TEST_ASSERT(!(q1<=x2)); 
  ARIADNE_TEST_ASSERT(q1> x2);
  ARIADNE_TEST_ASSERT(q1>=x2); 
  ARIADNE_TEST_ASSERT(!(q1< x2));
  
  double x1=1.0;
  ARIADNE_TEST_ASSERT(!(x1==q2)); 
  ARIADNE_TEST_ASSERT(x1!=q2); 
  ARIADNE_TEST_ASSERT(!(x1<=q2)); 
  ARIADNE_TEST_ASSERT(x1> q2);
  ARIADNE_TEST_ASSERT(x1>=q2); 
  ARIADNE_TEST_ASSERT(!(x1< q2));
}
  


int main() {
  cout << setprecision(20);
  TestRational().test();
  return ARIADNE_TEST_FAILURES;
}

