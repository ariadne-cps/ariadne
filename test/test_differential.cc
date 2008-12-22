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

#include "numeric.h"
#include "vector.h"
#include "sparse_differential.h"

#include "test.h"


using namespace Ariadne;
using namespace std;

template<class DF>
class TestDifferential {
    typedef typename DF::ScalarType X;
    typedef typename DF::ScalarType ScalarType;
    typedef Series<X> SeriesType;
    typedef DF DifferentialType;
    typedef DifferentialVector<DF> DifferentialVectorType;
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
        //assert((x1+x2)==DifferentialType("[3,2,0,0]"));
    }

    void test_add() {
        cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
        ARIADNE_TEST_EVALUATE(x1+x2);
        ARIADNE_TEST_EVALUATE(x1+c1);
        ARIADNE_TEST_EVALUATE(c1+x1);
        //assert((x1+x2)==DifferentialType("[3,2,0,0]"));
    }

    void test_sub() {
        cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
        ARIADNE_TEST_EVALUATE(x1-x2);
        ARIADNE_TEST_EVALUATE(x1-c1);
        ARIADNE_TEST_EVALUATE(c1-x1);
        //assert((x1-x2)==DifferentialType("[-1,0,0,0]"));
    }

    void test_mul() {
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

    void test_div() {
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

    void test_rec() {
        double a1[6]={ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
        ARIADNE_TEST_CONSTRUCT(DifferentialType,x1,(2,2,a1));
        ARIADNE_TEST_EQUAL(rec(rec(x1)),x1);
    }

    void test_pow() {
        cout << x2 << "^5 = " << pow(x2,5) << std::endl;
        //    assert(pow(x2,5)==DifferentialType("[32,80,160,240]"));
    }

    void test_compose() {
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


};


int main() {
#ifdef HAVE_GMPXX_H
    TestDifferential< SparseDifferential<Rational> > t1;
#endif
    cout << "INCOMPLETE " << flush;
    return ARIADNE_TEST_FAILURES;
}
