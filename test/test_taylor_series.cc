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
 private:
  TaylorSeries<X> x1,x2,x3;
 private:
 public:
  TestTaylorSeries() {
    x1=TaylorSeries<X>::variable(3,1.0);
    x2=TaylorSeries<X>::variable(3,2.0);
    double ax3[4]={2.0,1.0,0.1,0.0};
    x3=TaylorSeries<X>(3,ax3);
    
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
    ARIADNE_TEST_ASSERT(x1.degree()==3);
  }

  void test_add() {
    double ax1px2[4]={3,2,0,0};
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    ARIADNE_TEST_ASSERT((x1+x2)==TaylorSeries<X>(3,ax1px2));
  }

  void test_sub() {
    double ax1sx2[4]={-1,0,0,0};
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    ARIADNE_TEST_ASSERT((x1-x2)==TaylorSeries<X>(3,ax1sx2));
  }

  void test_mul() {
    double ax3tx4[4]={2,3,2,0};
    cout << x1 << "*" << x2 << " = " << x1*x2 << std::endl;
    ARIADNE_TEST_ASSERT((x1*x2)==TaylorSeries<X>(3,ax3tx4));
  }

  void test_div() {
    double ax3[3]={2,3,4};
    double ax4[3]={1,0,0};
    double ax4dx3[3]={0.5,-0.75,1.25};
    double aodx2[4]={0.5,-0.25,0.25,-0.375};
    double ax1dx2[4]={0.5,0.25,-0.25,0.375};
    TaylorSeries<X> x3(2,ax3);
    TaylorSeries<X> x4(2,ax4);
    cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
    ARIADNE_TEST_ASSERT((x3/x4)==x3);
    cout << x4 << "/" << x3 << " = " << x4/x3 << std::endl;
    ARIADNE_TEST_ASSERT((x4/x3)*x3==x4);
    ARIADNE_TEST_ASSERT((x4/x3)==TaylorSeries<X>(2,ax4dx3));
    cout << 1 << "/" << x2 << " = " << 1/x2 << std::endl;
    ARIADNE_TEST_ASSERT(((1/x2)*x2)==TaylorSeries<X>::constant(3,1.0));
    ARIADNE_TEST_ASSERT((1/x2)==TaylorSeries<X>(3,aodx2));
    cout << x1 << "/" << x2 << " = " << x1/x2 << std::endl;
    ARIADNE_TEST_ASSERT(((x1/x2)*x2)==x1);
    ARIADNE_TEST_ASSERT((x1/x2)==TaylorSeries<X>(3,ax1dx2));
  }

  void test_pow() {
    double ax2p5[4]={32,80,160,240};
    cout << x2 << "^5 = " << pow(x2,5) << std::endl;
    ARIADNE_TEST_ASSERT(pow(x2,5)==TaylorSeries<X>(3,ax2p5));
  }

  void test_compose() {
    double ax[6]={7,2,3,4,5,6};
    double ay[6]={11,2,-3,5,-8,13};
    double az[6]={11,4,-6,-6,65,148};
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
    ARIADNE_TEST_ASSERT(compose(y,x)==z);
    ARIADNE_TEST_ASSERT(compose(id,x)==x);
    ARIADNE_TEST_ASSERT(compose(x,id)==x);
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
