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

#include "test/test.h"

using namespace std;
using namespace Ariadne::Numeric;

class TestRational {
 public:
  TestRational() {
		cout << setprecision(20);
	}
	
	void test_constructors() {
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
	
	void test_stream() {
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

	void test_cmp()
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
      assert(q4<q5); 
      assert(q4!=q5);
      assert(!(q4>q5)); 
      assert(q5==q5);
    */

    // Test comparison with an integer
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

	void test_add()
  {
		ARIADNE_TEST_ASSERT(add(Rational(1.25),Rational(2.25))==Rational(7,2));
	}
	
	void test_sub()
	{
		ARIADNE_TEST_ASSERT(sub(Rational(1.25),Rational(2.25))==Rational(-1,1));
	}
	
	void test_mul() {
		ARIADNE_TEST_ASSERT(mul(Rational(1.25),Rational(2.25))==Rational(45,16));
	}
	
	void test_div() {
		ARIADNE_TEST_ASSERT(div(Rational(1.25),Rational(2.25))==Rational(5,9));
	}
	
	void test_med() {
		ARIADNE_TEST_ASSERT(med(Rational(1.25),Rational(2.25))==Rational(7,4));
	}
	
	void test_rad() {
		ARIADNE_TEST_ASSERT(rad(Rational(1.25),Rational(2.25))==Rational(1,2));
	}

	void test() {
		ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_stream());
    ARIADNE_TEST_CALL(test_cmp());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_med());
    ARIADNE_TEST_CALL(test_rad());
  }
};


int main() {
  TestRational().test(); 
  return ARIADNE_TEST_FAILURES;
}
