/***************************************************************************
 *            test_differential_vector.cc
 *
 *  Copyright  2007-8  Pieter Collins
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

#include "config.h"

#include "test.h"
#include "numeric.h"
#include "vector.h"
#include "series.h"
#include "sparse_differential.h"
#include "differential_vector.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class R, class A, class P>
void henon(R& r, const A& x, const P& p) 
{
  r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
  r[1]=x[0];
}

template<class DF>
DifferentialVector<DF> 
henon(const DifferentialVector<DF>& x, const Vector<typename DF::ScalarType>& p) 
{
  DifferentialVector<DF> r(2,2,x.degree()); henon(r,x,p); return r;
}



template<class DF>
class TestDifferentialVector {
  typedef typename DF::ScalarType X;
  typedef typename DF::ScalarType ScalarType;
  typedef Vector<X> VectorType;
  typedef Series<X> SeriesType;
  typedef DF DifferentialType;
  typedef DifferentialVector<DF> DifferentialVectorType;
 private:
  DifferentialVectorType x1,x2,x3;
 public:
  TestDifferentialVector() {
    double a1[15]={ 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a2[15]={ 3.0, 1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double a3[15]={ 2.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    x1=DifferentialVectorType(1,2,4,a1);
    x2=DifferentialVectorType(1,2,4,a2);
    x3=DifferentialVectorType(1,1,4,a3);
    
    ARIADNE_TEST_CALL(test_degree());
    ARIADNE_TEST_CALL(test_add());
    ARIADNE_TEST_CALL(test_sub());
    ARIADNE_TEST_CALL(test_mul());
    ARIADNE_TEST_CALL(test_div());
    ARIADNE_TEST_CALL(test_evaluate());
    ARIADNE_TEST_CALL(test_differentiate());
    ARIADNE_TEST_CALL(test_translate());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_mapping());
    ARIADNE_TEST_CALL(test_inverse());
    ARIADNE_TEST_CALL(test_implicit());
  }

  void test_degree() {
    ARIADNE_TEST_ASSERT(x1.degree()==4);
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    //assert((x1+x2)==DifferentialVectorType("[3,2,0,0]"));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    //assert((x1-x2)==DifferentialVectorType("[-1,0,0,0]"));
  }

  void test_mul() {
    X c=2;
    cout << x1 << "*" << c << " = " << x1*c << std::endl;
    cout << c << "*" << x1 << " = " << c*x1 << std::endl;
    //assert((x1*x2)==DifferentialVectorType("[2,3,2,0]"));
  }

  void test_div() {
    X c=2;
    cout << x1 << "/" << c << " = " << x1/c << std::endl;
  }

  void test_evaluate() {
    Float ac[2]={1,2}; Float adv[10]={1,2,3,4,5,6,7,8,9,10};
    Vector<X> c(2u,ac); 
    DifferentialVectorType dv(1u,2u,3u,adv);
    std::cout << "c=" << c << std::endl;
    std::cout << "dv="<< dv << std::endl;
    std::cout << "v(c)" << evaluate(dv,c) << std::endl;
    std::cout << std::endl;
  }

  void test_translate() {
    Float ac[2]={0,1}; Float adv[10]={1,2,3,4,5,6,7,8,9,10};
    Vector<X> c(2u,ac); 
    Vector<X> mc=-c;
    DifferentialVectorType dv(1u,2u,2u,adv);
    std::cout << "c=" << c << std::endl;
    std::cout << "dv="<< dv << std::endl;
    std::cout << "dv(x+c)" << translate(dv,c) << std::endl;
    std::cout << "dv((x+c)-c)" <<  translate(translate(dv,c),mc) << std::endl;
    std::cout << "dv(x-c)" << translate(dv,mc) << std::endl;
    std::cout << "dv((x+c)-c)" <<  translate(translate(dv,mc),c) << std::endl;
    ARIADNE_ASSERT_EQUAL(translate(translate(dv,c),mc),dv);
  }

  void test_differentiate() {  
    double a[]={ 1,2,3,4,5,6,7,8,9,10 };
    DifferentialType y(2,3,a);
    cout << "y=" << y << endl;
    ARIADNE_ASSERT_EQUAL(derivative(antiderivative(y,0),0),y);
  }

  void test_compose() {
    //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double ax[10] = { 3.0, 1.0, 0.0, 0.0, 0.125, 0.25, 0.0, 0.0, 0.0, 0.0 };
    double ay[4] = { 1.0, -1.0, 0.5, -0.25 };
    double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
    DifferentialVectorType x(1,2,3,ax);
    DifferentialType y(1,3,ay);
    DifferentialVectorType id(1,1,3,aid);
    cout << "x=" << x << endl;
    cout << "y=" << y << endl;
    cout << "compose(y,x)=" << compose(y,x) << endl;
    cout << "compose(id,x)=" << compose(id,x) << endl;
    ARIADNE_TEST_EQUAL(compose(id,x),x);
    cout << "compose(id,x)-x=" << compose(id,x)-x << endl;
  }

  void test_mapping() {
    DifferentialVectorType x(2,2,2); 
    x[0][MultiIndex::unit(2,0)]=1; x[1][MultiIndex::unit(2,1)]=1; 
    cout << "x=" << x << endl;
    Vector<Float> p(2); p[0]=1.5; p[1]=0.375;
    double ahxp[12]={ 1.5, 0.0, -0.375, -1.0, 0.0, 0.0,   0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
    DifferentialVectorType hxp(2,2,2,ahxp);
    ARIADNE_TEST_EQUAL(henon(x,p),hxp);
  }

  void test_inverse() {
    double ax[12]={ 0.0, 2.0, 1.0, 3.0, 4.0, 5.0,   0.0, 1.0, 1.0, 2.0, 3.0, 4.0 };
    VectorType c(2);
    DifferentialVectorType id=DifferentialVectorType::variable(2,2,2,c);
    DifferentialVectorType x(2,2,2,ax);
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
    VectorType c(2);
    DifferentialVectorType id1=DifferentialVectorType::variable(1,1,2,VectorType(1));
    DifferentialVectorType id2=DifferentialVectorType::variable(2,2,2,VectorType(2));
    DifferentialVectorType id3=DifferentialVectorType::variable(3,3,2,VectorType(3));
    ARIADNE_TEST_PRINT(id3);
    DifferentialVectorType x(2,4,2,ax);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_PRINT(implicit(x));
    DifferentialVectorType y=implicit(x);
    DifferentialVectorType z=join(DifferentialVectorType::variable(2,2,2,VectorType(2)),y);
    ARIADNE_TEST_PRINT(z);
    ARIADNE_TEST_EQUAL(compose(x,z),DifferentialVectorType::constant(2,2,2,VectorType(2)));
    
  }
};


int main() {
  TestDifferentialVector< SparseDifferential<Rational> > t1;
  cout << "INCOMPLETE " << flush;
  return ARIADNE_TEST_FAILURES;
}
