/***************************************************************************
 *            test_sparse_differential.cc
 *
 *  Copyright  2008 Pieter Collins
 *
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "differentiation/sparse_differential.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class X>
class TestSparseDifferential {
 private:
  X c1;
  SparseDifferential<X> x1,x2,x3;
 public:
  TestSparseDifferential() {
    double a1[15]={ 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a2[15]={ 3.0, 1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a3[15]={ 2.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    c1=3.0;
    x1=SparseDifferential<X>(2,4,a1);
    x2=SparseDifferential<X>(2,4,a2);
    x3=SparseDifferential<X>(1,4,a3);
    
    ARIADNE_TEST_CALL(test_degree());
    ARIADNE_TEST_CALL(test_neg());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_rec());
    ARIADNE_TEST_CALL(test_pow());
    ARIADNE_TEST_CALL(test_compose());
  }

  void test_degree() {
    ARIADNE_TEST_ASSERT(x1.degree()==4);
  }

  void test_neg() {
    cout << -x1 << " = " << -x1 << std::endl;
    //assert((x1+x2)==SparseDifferential<X>("[3,2,0,0]"));
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    ARIADNE_TEST_EVALUATE(x1+x2);
    ARIADNE_TEST_EVALUATE(x1+c1);
    ARIADNE_TEST_EVALUATE(c1+x1);
    //assert((x1+x2)==SparseDifferential<X>("[3,2,0,0]"));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    ARIADNE_TEST_EVALUATE(x1-x2);
    ARIADNE_TEST_EVALUATE(x1-c1);
    ARIADNE_TEST_EVALUATE(c1-x1);
    //assert((x1-x2)==SparseDifferential<X>("[-1,0,0,0]"));
  }

  void test_mul() {
    double a1[6]={ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    double a2[6]={ 2.0, 3.0, 5.0, 7.0, 11.0, 13.0 };
    double a1m2[6]={ 2.0, 7.0, 11.0, 21.0, 40.0, 40.0 };
    double acm2[6]={ 10.0, 15.0, 25.0, 35.0, 55.0, 65.0 };
    SparseDifferential<X> x1(2,2,a1);
    SparseDifferential<X> x2(2,2,a2);
    SparseDifferential<X> x1mx2(2,2,a1m2);
    SparseDifferential<X> cmx2(2,2,acm2);
    X c=5;
    ARIADNE_TEST_EQUAL(x1*x2,x1mx2);
    ARIADNE_TEST_EQUAL(c*x2,cmx2);
    ARIADNE_TEST_EQUAL(x2*c,cmx2);
  }

  void test_div() {
    ARIADNE_TEST_EVALUATE(x1/x2);
    ARIADNE_TEST_EVALUATE(x1/c1);
    ARIADNE_TEST_EVALUATE(c1/x1);
    /*
    SparseDifferential<X> x3("[2,3,4]");
    SparseDifferential<X> x4("[1,0,0]");
    cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
    assert((x3/x4)==x3);
    cout << x4 << "/" << x3 << " = " << x4/x3 << std::endl;
    assert((x4/x3)==SparseDifferential<X>("[0.5,-0.75,1.25]"));
    cout << 1 << "/" << x2 << " = " << 1/x2 << std::endl;
    assert((1/x2)==SparseDifferential<X>("[0.5,-0.25,0.25,-0.375]"));
    cout << x1 << "/" << x2 << " = " << x1/x2 << std::endl;
    assert((x1/x2)==SparseDifferential<X>("[0.5,0.25,-0.25,0.375]"));
    */
  }

  void test_rec() {
    double a1[6]={ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    ARIADNE_TEST_CONSTRUCT(SparseDifferential<X>,x1,(2,2,a1));
    ARIADNE_TEST_EQUAL(rec(rec(x1)),x1);
  }

  void test_pow() {
    cout << x2 << "^5 = " << pow(x2,5) << std::endl;
    //    assert(pow(x2,5)==SparseDifferential<X>("[32,80,160,240]"));
  }

  void test_compose() {
    //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ax[10] = { 3.0,  1.0, 2.0,  1.0, 0.5, 2.0,  0.0, 0.0, 0.0, 0.0 };
    double ay[4] = { 1.0, 2.0, -3.0, 5.0 };
    double ayx[10] = { 1.0,  2.0, 4.0,  -1.0, -11.0, -8.0,  -1.0, 15.0, 42.0, 16.0 };
    double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
    ARIADNE_TEST_CONSTRUCT(SparseDifferential<X>,x,(2,3,ax));
    ARIADNE_TEST_CONSTRUCT(SparseSeries<X>,y,(3,ay));
    ARIADNE_TEST_CONSTRUCT(SparseSeries<X>,id,(3,aid));
    ARIADNE_TEST_EQUAL(compose(y,x),SparseDifferential<X>(2,3,ayx));
    ARIADNE_TEST_EQUAL(compose(id,x),x);
  }


};


template<class X>
class TestSparseDifferentialVector {
  SparseDifferentialVector<X> x1,x2,x3;
 public:
  TestSparseDifferentialVector() {
    double a1[15]={ 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a2[15]={ 3.0, 1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a3[15]={ 2.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    x1=SparseDifferentialVector<X>(1,2,4,a1);
    x2=SparseDifferentialVector<X>(1,2,4,a2);
    x3=SparseDifferentialVector<X>(1,1,4,a3);
    
    ARIADNE_TEST_CALL(test_degree());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_inverse());
    ARIADNE_TEST_CALL(test_implicit());
  }

  void test_degree() {
    ARIADNE_TEST_ASSERT(x1.degree()==4);
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << SparseDifferentialVector<X>(x1+x2) << std::endl;
    //assert((x1+x2)==SparseDifferentialVector<X>("[3,2,0,0]"));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << SparseDifferentialVector<X>(x1-x2) << std::endl;
    //assert((x1-x2)==SparseDifferentialVector<X>("[-1,0,0,0]"));
  }

  void test_mul() {
    //X c=2;
    //cout << x1 << "*" << c << " = " << SparseDifferentialVector<X>(x1*c) << std::endl;
    //cout << c << "*" << x1 << " = " << SparseDifferentialVector<X>(c*x1) << std::endl;
    //assert((x1*x2)==SparseDifferentialVector<X>("[2,3,2,0]"));
  }

  void test_div() {
    //X c=2;
    //cout << x1 << "/" << c << " = " << SparseDifferentialVector<X>(x1/c) << std::endl;
  }

  void test_compose() {
    //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ax[10] = { 3.0, 1.0, 0.0, 0.0, 0.125, 0.25, 0.0, 0.0, 0.0, 0.0 };
    double ay[4] = { 1.0, -1.0, 0.5, -0.25 };
    double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
    SparseDifferentialVector<X> x(1,2,3,ax);
    SparseDifferential<X> y(1,3,ay);
    SparseDifferentialVector<X> id(1,1,3,aid);
    cout << "x=" << x << endl;
    cout << "y=" << y << endl;
    cout << "compose(y,x)=" << compose(y,x) << endl;
    cout << "compose(id,x)=" << compose(id,x) << endl;
    SparseDifferentialVector<X> z = compose(id,x); z-=x;
    cout << "compose(id,x)-x=" << z << endl;
    
  }

  void test_inverse() {
    double ax[12]={ 0.0, 2.0, 1.0, 3.0, 4.0, 5.0,   0.0, 1.0, 1.0, 2.0, 3.0, 4.0 };
    Vector<X> c(2);
    SparseDifferentialVector<X> id=SparseDifferentialVector<X>::variable(2,2,2,c);
    SparseDifferentialVector<X> x(2,2,2,ax);
    ARIADNE_TEST_PRINT(c);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_PRINT(inverse(x,c));
    ARIADNE_TEST_EQUAL(compose(x,inverse(x,c)),id);
    ARIADNE_TEST_EQUAL(compose(inverse(x,c),x),id);
    ARIADNE_TEST_EQUAL(inverse(inverse(x,c),c),x);
  }

  void test_implicit() {
    double ax[30]={ 0.0,  2.0,1.0,3.0,1.0, 4.0,5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                    0.0,  1.0,1.0,2.0,1.0, 3.0,4.0,0.0,6.0,0.0,7.0,0.0,0.0,0.0,0.0 };
    Vector<X> c(2);
    SparseDifferentialVector<X> id1=SparseDifferentialVector<X>::variable(1,1,2,Vector<X>(1));
    SparseDifferentialVector<X> id2=SparseDifferentialVector<X>::variable(2,2,2,Vector<X>(2));
    SparseDifferentialVector<X> id3=SparseDifferentialVector<X>::variable(3,3,2,Vector<X>(3));
    ARIADNE_TEST_PRINT(id3);
    SparseDifferentialVector<X> x(2,4,2,ax);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_PRINT(implicit(x));
    SparseDifferentialVector<X> y=implicit(x);
    SparseDifferentialVector<X> z=join(SparseDifferentialVector<X>::variable(2,2,2,Vector<X>(2)),y);
    ARIADNE_TEST_PRINT(z);
    ARIADNE_TEST_EQUAL(compose(x,z),SparseDifferentialVector<X>::constant(2,2,2,Vector<X>(2)));
    
  }
};


int main() {
  TestSparseDifferential<Rational> t1;
  TestSparseDifferentialVector<Rational> t2;
  cout << "INCOMPLETE " << flush;
  return ARIADNE_TEST_FAILURES;
}
