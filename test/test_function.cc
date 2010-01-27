/***************************************************************************
 *            test_function.cc
 *
 *  Copyright  2009  Pieter Collins
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
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.h"


#include "function.h"
#include "predicate.h"
#include "numeric.h"

#include "test.h"


using namespace std;
using namespace Ariadne;

class TestFunction 
{
  public:
    void test();
};

void TestFunction::test()
{

    ScalarFunction sf1(3);
    ScalarFunction sf2(3);
    ScalarFunction sf3(3);
    //Vector<ScalarFunction> ve=join(3,*e1,*e2,*e3);

    VectorFunction vf=join(sf1,sf2);

    ScalarFunction g(2);
    ScalarFunction h=compose(g,vf);

    VectorFunction jf=join(sf1,sf2);
    //ScalarFunction cf=combine(sf1,sf2);

    Polynomial<Real> p;
    ScalarFunction pf(p);
    ARIADNE_TEST_EQUAL(p,pf.polynomial());

    Vector<Float> b(1); Matrix<Float> A(1,1);
    b[0]=1.0/3.0; A[0][0]=1.0;
    VectorAffineFunction aff(A,b);
    ARIADNE_TEST_EQUAL(aff.ib(),b);
    ARIADNE_TEST_EQUAL(aff.iA(),A);
    
    Vector<Interval> ib(1); Matrix<Interval> iA(1,1);
    ib[0]=Interval(1.0)/Interval(3.0);
    iA[0][0]=1.0;
    VectorAffineFunction iaff(iA,ib);
    ARIADNE_TEST_EQUAL(iaff.ib(),ib);
    ARIADNE_TEST_EQUAL(iaff.iA(),iA);
    
#ifdef HAVE_RATIONAL
    Vector<Rational> rb(1); Matrix<Rational> rA(1,1);
    rb[0]=Rational(1.0,3.0);
    rA[0][0]=1.0;
    VectorAffineFunction raff(rA,rb);
    ARIADNE_TEST_EQUAL(raff.ib(),rb);
    ARIADNE_TEST_EQUAL(raff.iA(),rA);
#endif // HAVE_RATIONAL
    
}

int main() {
    TestFunction().test();

    return ARIADNE_TEST_FAILURES;
}

