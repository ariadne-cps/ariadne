/***************************************************************************
 *            test_taylor_series.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "test_float.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "function/taylor_series.h"
#include "function/taylor_variable.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Function;
using namespace std;

template<class X>
bool operator==(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2) {
  return ts1.data()==ts2.data(); 
}

template<class X>
class TestTaylorSeries {
 public:
  TestTaylorSeries() {
    
    ARIADNE_TEST_CALL(test_degree());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_pow());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_inverse());
  }

  void test_degree() {
    TaylorSeries<X> x1=TaylorSeries<X>::variable(3,1.0);
    ARIADNE_TEST_ASSERT(x1.degree()==3);
  }

  void test_add() {
    TaylorSeries<X> x1=TaylorSeries<X>::variable(3,1.0);
    TaylorSeries<X> x2=TaylorSeries<X>::variable(3,2.0);
    double ax1px2[4]={3,2,0,0};
    ARIADNE_TEST_ASSERT((x1+x2)==TaylorSeries<X>(3,ax1px2));
  }

  void test_sub() {
    TaylorSeries<X> x1=TaylorSeries<X>::variable(3,1.0);
    TaylorSeries<X> x2=TaylorSeries<X>::variable(3,2.0);
    double ax1sx2[4]={-1,0,0,0};
    ARIADNE_TEST_ASSERT((x1-x2)==TaylorSeries<X>(3,ax1sx2));
  }

  void test_mul() {
    double ax1[5]={1,2,3,4,5};
    double ax2[5]={1,3,5,2,4};
    double ax1mx2[5]={1,5,14,25,40};
    TaylorSeries<X> x1(4,ax1);
    TaylorSeries<X> x2(4,ax2);
    TaylorSeries<X> x1mx2(4,ax1mx2);
    ARIADNE_TEST_EQUAL((x1*x2),x1mx2);
  }

  void test_div() {
    double ax1[4]={1,1,0,0};
    double ax2[4]={2,1,0,0};
    double ax3[3]={2,3,2};
    double ax4[3]={1,0,0};
    double ax4dx3[3]={0.5,-0.75,0.625};
    double aodx2[4]={0.5,-0.25,0.125,-0.0625};
    double ax1dx2[4]={0.5,0.25,-0.125,0.0625};
    TaylorSeries<X> x1(3,ax1);
    TaylorSeries<X> x2(3,ax2);
    TaylorSeries<X> x3(2,ax3);
    TaylorSeries<X> x4(2,ax4);
    ARIADNE_TEST_EQUAL((x3/x4),x3);
    ARIADNE_TEST_EQUAL((x4/x3),TaylorSeries<X>(2,ax4dx3));
    ARIADNE_TEST_EQUAL((x4/x3)*x3,x4);
    ARIADNE_TEST_EQUAL(((1/x2)*x2),TaylorSeries<X>::constant(3,1.0));
    ARIADNE_TEST_EQUAL((1/x2),TaylorSeries<X>(3,aodx2));
    ARIADNE_TEST_EQUAL(((x1/x2)*x2),x1);
    ARIADNE_TEST_EQUAL((x1/x2),TaylorSeries<X>(3,ax1dx2));
  }

  void test_pow() {
    TaylorSeries<X> x2=TaylorSeries<X>::variable(3,2.0);
    double ax2p5[4]={32,80,80,40};
    ARIADNE_TEST_EQUAL(pow(x2,5),TaylorSeries<X>(3,ax2p5));
  }

  void test_compose() {
    double ax[6]={7,2,3,4,5,6};
    double ay[6]={11,2,-3,5,-8,13};
    double az[6]={11,4,-6,12,-13,38};
    double aid[6]={ax[0],1};
    TaylorSeries<X> x(5,ax);
    TaylorSeries<X> y(5,ay);
    TaylorSeries<X> z(5,az);
    TaylorSeries<X> id(5,aid);
    cout << "x="<<x<<"\n";
    cout << "y="<<y<<"\n";
    cout << "z="<<z<<"\n";
    cout << "compose(y,x)="<<compose(y,x)<<"\n";
    cout << "compose(id,x)="<<compose(id,x)<<"\n";
    cout << "compose(x,id)="<<compose(x,id)<<"\n";
    ARIADNE_TEST_EQUAL(compose(y,x),z);
    ARIADNE_TEST_EQUAL(compose(id,x),x);
    ARIADNE_TEST_EQUAL(compose(x,id),x);
  }

  void test_inverse() {
    double a3[6]={0,2,3,4,5,6};
    TaylorSeries<X> x3(5,a3);
    cout << "x3="<<x3<<"\n";
    cout << "inverse(x3)="<<inverse(x3,X(0)) << "\n";
    cout << "inverse(inverse(x3))="<<inverse(inverse(x3,X(0)),X(0)) << "\n";
  }

};


int main() {
  TestTaylorSeries<Rational> t1;
  
  return ARIADNE_TEST_FAILURES;
}
