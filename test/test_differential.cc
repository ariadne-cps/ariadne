/***************************************************************************
 *            test_differential.cc
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

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/differential.h"

#include "test.h"


using namespace Ariadne;
using namespace std;

template<class R, class A, class P>
Void henon(R& r, const A& x, const P& p)
{
    r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
    r[1]=x[0];
}

template<class DF>
Vector<DF>
henon(const Vector<DF>& x, const Vector<typename DF::NumericType>& p)
{
    Vector<DF> r(2,2,x.degree()); henon(r,x,p); return r;
}

// FIXME: Operators needed to prevent bad overloads
template<class X> Differential<X> operator*(const typename Differential<X>::NumericType& c,const NonAssignableDifferential<X>& dx) {
    return c * static_cast<Differential<X>const&>(dx); }


template<class DF>
class TestDifferential {
    typedef typename DF::ValueType X;
    typedef X ScalarType;
    typedef typename DF::SeriesType SeriesType;
    typedef DF DifferentialType;
    typedef Vector<DF> DifferentialVectorType;
  private:
    X c1;
    DifferentialType x1,x2,x3;
  public:
    TestDifferential() {
        double a1[15]={ 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double a2[15]={ 3.0, 1.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double a3[15]={ 2.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        c1=3.0;
        x1=DifferentialType(2,4,a1);
        x2=DifferentialType(2,4,a2);
        x3=DifferentialType(1,4,a3);

        ARIADNE_TEST_PRINT(x1);
        ARIADNE_TEST_PRINT(x2);

        ARIADNE_TEST_CALL(test_degree());
        ARIADNE_TEST_CALL(test_neg());
        ARIADNE_TEST_CALL(test_add());
        ARIADNE_TEST_CALL(test_sub());
        ARIADNE_TEST_CALL(test_mul());
        ARIADNE_TEST_CALL(test_div());
        ARIADNE_TEST_CALL(test_rec());
        ARIADNE_TEST_CALL(test_pow());
        ARIADNE_TEST_CALL(test_compose());
        ARIADNE_TEST_CALL(test_gradient());
        ARIADNE_TEST_CALL(test_hessian());
    }

    Void test_degree() {
        ARIADNE_TEST_ASSERT(x1.degree()==4);
    }

    Void test_neg() {
        DifferentialType nx1(2,4, { {{0,0},-2.0}, {{1,0},-1.0}, {{2,0},-0.5} });
        ARIADNE_TEST_EQUALS(-x1,nx1);
    }

    Void test_add() {
        DifferentialType x1px2(2,4, { {{0,0},5.0}, {{1,0},2.0}, {{2,0},0.75} });
        ARIADNE_TEST_EQUALS(x1+x2,x1px2);
        ARIADNE_TEST_EVALUATE(x1+c1);
        ARIADNE_TEST_EVALUATE(c1+x1);
        //assert((x1+x2)==DifferentialType("[3,2,0,0]"));
    }

    Void test_sub() {
        DifferentialType x1mx2(2,4, { {{0,0},-1.0}, {{1,0},0.0}, {{2,0},0.25} });
        ARIADNE_TEST_EQUALS(x1-x2,x1mx2);
        ARIADNE_TEST_EVALUATE(x1-c1);
        ARIADNE_TEST_EVALUATE(c1-x1);
        //assert((x1-x2)==DifferentialType("[-1,0,0,0]"));
    }

    Void test_mul() {
        double a1[6]={ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        double a2[6]={ 2.0, 3.0, 5.0, 7.0, 11.0, 13.0 };
        double a1m2[6]={ 2.0, 7.0, 11.0, 21.0, 40.0, 40.0 };
        double acm2[6]={ 10.0, 15.0, 25.0, 35.0, 55.0, 65.0 };
        DifferentialType x1(2,2,a1);
        DifferentialType x2(2,2,a2);
        DifferentialType x1mx2(2,2,a1m2);
        DifferentialType cmx2(2,2,acm2);
        X c=5;
        ARIADNE_TEST_EQUAL(x1*x2,x1mx2);
        ARIADNE_TEST_EQUAL(c*x2,cmx2);
        ARIADNE_TEST_EQUAL(x2*c,cmx2);
    }

    Void test_div() {
        ARIADNE_TEST_EVALUATE(x1/x2);
        ARIADNE_TEST_EVALUATE(x1/c1);
        ARIADNE_TEST_EVALUATE(c1/x1);
        /*
          DifferentialType x3("[2,3,4]");
          DifferentialType x4("[1,0,0]");
          cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
          assert((x3/x4)==x3);
          cout << x4 << "/" << x3 << " = " << x4/x3 << std::endl;
          assert((x4/x3)==DifferentialType("[0.5,-0.75,1.25]"));
          cout << 1 << "/" << x2 << " = " << 1/x2 << std::endl;
          assert((1/x2)==DifferentialType("[0.5,-0.25,0.25,-0.375]"));
          cout << x1 << "/" << x2 << " = " << x1/x2 << std::endl;
          assert((x1/x2)==DifferentialType("[0.5,0.25,-0.25,0.375]"));
        */
    }

    Void test_rec() {
        double a1[6]={ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        ARIADNE_TEST_CONSTRUCT(DifferentialType,x1,(2,2,a1));
        ARIADNE_TEST_EQUAL(rec(rec(x1)),x1);
    }

    Void test_pow() {
        cout << x2 << "^5 = " << pow(x2,5) << std::endl;
        //    assert(pow(x2,5)==DifferentialType("[32,80,160,240]"));
    }

    Void test_compose() {
        //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double ax[10] = { 3.0,  1.0, 2.0,  1.0, 0.5, 2.0,  0.0, 0.0, 0.0, 0.0 };
        double ay[4] = { 1.0, 2.0, -3.0, 5.0 };
        double ayx[10] = { 1.0,  2.0, 4.0,  -1.0, -11.0, -8.0,  -1.0, 15.0, 42.0, 16.0 };
        double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
        ARIADNE_TEST_CONSTRUCT(DifferentialType,x,(2,3,ax));
        ARIADNE_TEST_CONSTRUCT(SeriesType,y,(3,ay));
        ARIADNE_TEST_CONSTRUCT(SeriesType,id,(3,aid));
        ARIADNE_TEST_EQUAL(compose(y,x),DifferentialType(2,3,ayx));
        ARIADNE_TEST_EQUAL(compose(id,x),x);
    }

    Void test_gradient() {
        // Regression test based on errors in Henon evaluation.
        Vector<ApproximateFloat> x={0.875,-0.125};
        Vector< Differential<ApproximateFloat> > dx=Differential<ApproximateFloat>::variables(1u,x);
        Differential<ApproximateFloat> dfx=1.5-dx[0]*dx[0]-0.25*dx[1];
        ARIADNE_TEST_PRINT(dfx);
        Covector<ApproximateFloat> g = dfx.gradient();
        ARIADNE_TEST_PRINT(g);
        ARIADNE_TEST_EQUALS(g[0],-1.75);
        ARIADNE_TEST_EQUALS(g[1],-0.25);
    }

    Void test_hessian() {
        // Test Hessian matrix of
        ApproximateFloat a00=1.5; ApproximateFloat a01=2.5; ApproximateFloat a11=3.5;
        double x0=0.875; double x1=-1.25;
        Vector<ApproximateFloat> x={x0,x1};
        Vector< Differential<ApproximateFloat> > dx=Differential<ApproximateFloat>::variables(2u,x);
        Differential<ApproximateFloat> dfx=a00*dx[0]*dx[0]+a01*dx[0]*dx[1]+a11*dx[1]*dx[1];
        ARIADNE_TEST_PRINT(dfx);
        Matrix<ApproximateFloat> H = dfx.hessian();
        ARIADNE_TEST_PRINT(H);
        ARIADNE_TEST_EQUAL(H[0][1],H[1][0]);
        ARIADNE_TEST_EQUALS(H[0][0],a00*2);
        ARIADNE_TEST_EQUALS(H[0][1],a01*2);
        ARIADNE_TEST_EQUALS(H[1][1],a11*2);
    }


};



template<class DF>
class TestDifferentialVector {
    typedef typename DF::ValueType X;
    typedef X ScalarType;
    typedef Vector<X> VectorType;
    typedef Series<X> SeriesType;
    typedef DF DifferentialType;
    typedef Vector<DF> DifferentialVectorType;
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
        ARIADNE_TEST_CALL(test_jacobian());
        ARIADNE_TEST_CALL(test_differentiate());
        ARIADNE_TEST_CALL(test_compose());
        ARIADNE_TEST_CALL(test_mapping());
    }

    Void test_degree() {
        ARIADNE_TEST_EQUAL(x1.degree(),4);

        // Regression test to check setting of degree for null vector
        DifferentialVectorType dx(0u,2u,3u);

        if(dx.argument_size()!=2u) {
            ARIADNE_TEST_WARN("Vector<Differential<X>>::argument_size() returns 0 if the vector of differentials has no elements.");
        }
        if(dx.degree()!=3u) {
            ARIADNE_TEST_WARN("Vector<Differential<X>>::degree() returns 0 if the vector of differentials has no elements.");
        }

    }

    Void test_add() {
        cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
        //assert((x1+x2)==DifferentialVectorType("[3,2,0,0]"));
    }

    Void test_sub() {
        cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
        //assert((x1-x2)==DifferentialVectorType("[-1,0,0,0]"));
    }

    Void test_mul() {
        X c=2;
        cout << x1 << "*" << c << " = " << x1*c << std::endl;
        cout << c << "*" << x1 << " = " << c*x1 << std::endl;
        //assert((x1*x2)==DifferentialVectorType("[2,3,2,0]"));
    }

    Void test_div() {
        X c=2;
        cout << x1 << "/" << c << " = " << x1/c << std::endl;
    }

    Void test_jacobian() {
        DifferentialVectorType dv(0u,2u,3u);
        Matrix<ApproximateFloat> J=dv.jacobian();
        ARIADNE_TEST_EQUALS(J.row_size(),dv.result_size());
        ARIADNE_TEST_EQUALS(J.column_size(),dv.argument_size());
    }

    Void test_differentiate() {
        double ax[]={ 2,12,7,33,24,13 };
        DifferentialType x(2,2,ax);
        double ay[]={ 1,2,3,6,7,9,11,12,13,14 };
        DifferentialType y(2,3,ay);
        double az[]={ 0, 1.0,0.0, 1.0,3.0,0.0, 2.0,3.5,9.0,0.0, 2.75,4.0,6.5,14.0,0.0 };
        DifferentialType z(2,4,az);
        ARIADNE_TEST_EQUAL(derivative(y,0),x);
        ARIADNE_TEST_EQUAL(antiderivative(y,0),z);
        ARIADNE_TEST_EQUAL(derivative(antiderivative(y,0),0),y);
    }

    Void test_compose() {
        //double ax[10] = { 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        double ax[10] = { 3.0, 1.0, 0.0, 0.0, 0.125, 0.25, 0.0, 0.0, 0.0, 0.0 };
        double ay[4] = { 1.0, -1.0, 0.5, -0.25 };
        double az[10] = { 1.0, -1.0,0.0, 0.5,-0.125,-0.25, -0.25,0.125,0.25,0.0 };
        double aid[4] = { ax[0], 1.0, 0.0, 0.0 };
        DifferentialVectorType x(1,2,3,ax);
        DifferentialType y(1,3,ay);
        DifferentialType z(2,3,az);
        DifferentialVectorType id(1,1,3,aid);
        ARIADNE_TEST_PRINT(x);
        ARIADNE_TEST_PRINT(y);
        ARIADNE_TEST_PRINT(z);
        ARIADNE_TEST_EQUAL(compose(y,x),z);
        ARIADNE_TEST_EQUAL(compose(id,x),x);
    }

    Void test_mapping() {
        DifferentialVectorType x(2,2,2);
        x[0][MultiIndex::unit(2,0)]=1; x[1][MultiIndex::unit(2,1)]=1;
        cout << "x=" << x << endl;
        Vector<ApproximateFloat> p(2); p[0]=1.5; p[1]=0.375;
        double ahxp[12]={ 1.5, 0.0, -0.375, -1.0, 0.0, 0.0,   0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
        DifferentialVectorType hxp(2,2,2,ahxp);
        ARIADNE_TEST_EQUAL(henon(x,p),hxp);
    }
};

Int main() {
    TestDifferential< Differential<ApproximateFloat> > tf;
    TestDifferentialVector< Differential<ApproximateFloat> > tfv;
    return ARIADNE_TEST_FAILURES;
}
