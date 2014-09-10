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

#include "numeric.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestFunction
{
  public:
    void test();
  private:
    void test_concept();
    void test_scalar_function();
    void test_vector_function();
    void test_conversions();
    void test_differentiation();
};

void TestFunction::test()
{
    ARIADNE_TEST_CALL(test_scalar_function());
    ARIADNE_TEST_CALL(test_vector_function());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_differentiation());
}

void TestFunction::test_concept()
{

    EffectiveScalarFunction sf1(3);
    EffectiveScalarFunction sf2(3);
    EffectiveScalarFunction sf3(3);
    //Vector<EffectiveScalarFunction> ve=join(3,*e1,*e2,*e3);

    EffectiveVectorFunction vf=join(sf1,sf2);

    EffectiveScalarFunction g(2);
    EffectiveScalarFunction h=compose(g,vf);

    EffectiveVectorFunction jf=join(sf1,sf2);
    //EffectiveScalarFunction cf=combine(sf1,sf2);

    //Polynomial<Real> p;
    //EffectiveScalarFunction pf(p);

    //Vector<Float> b; Matrix<Float> A;
    //VectorAffineFunction aff(A,b);

}

void TestFunction::test_scalar_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,o,constant(3,1));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,x,coordinate(3,0));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,y,coordinate(3,1));

    ARIADNE_TEST_CONSTRUCT(EffectiveScalarFunction,f,(o+x*y));
    ARIADNE_TEST_CONSTRUCT(Vector<Float>,p,({2.0,3.0,5.0}));
    ARIADNE_TEST_EQUAL(f(p),7.0);

    ARIADNE_TEST_PRINT(cos(f));

    EffectiveScalarFunction df=f.derivative(1);
    ARIADNE_TEST_PRINT(df);
    ARIADNE_TEST_EQUAL(df(p),2.0);
}

void TestFunction::test_vector_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveVectorFunction,f,identity(3));

    // Regression tests for element
    ARIADNE_TEST_PRINT(f[0]);
    const EffectiveVectorFunction& fcr=f;
    ARIADNE_TEST_PRINT(fcr[0]);
    EffectiveVectorFunction& fr=f;
    ARIADNE_TEST_PRINT(fr[0]);

    ARIADNE_TEST_EQUAL(f[0](Vector<Float>{2.0,3.0,5.0}),2.0);
    ARIADNE_TEST_EXECUTE(f[0]=f[1]);
    ARIADNE_TEST_EQUAL(f[0](Vector<Float>{2.0,3.0,5.0}),3.0);
}





void TestFunction::test_conversions()
{

}

void TestFunction::test_differentiation()
{
    EffectiveScalarFunction z=EffectiveScalarFunction::constant(2,0);
    EffectiveScalarFunction o=EffectiveScalarFunction::constant(2,1);
    EffectiveScalarFunction x=EffectiveScalarFunction::coordinate(2,0);
    EffectiveScalarFunction y=EffectiveScalarFunction::coordinate(2,1);

    EffectiveScalarFunction af=3*x-2*y+1;
    EffectiveScalarFunction daf=af.derivative(1);
    ARIADNE_TEST_EQUAL(daf.evaluate(Vector<Float>{2.4,1.3}),-2.0);

    ARIADNE_TEST_EQUAL(x.derivative(0).evaluate(Vector<Float>{2.4,1.3}),1.0);
    ARIADNE_TEST_EQUAL(x.derivative(0).evaluate(Vector<Interval>{{2.4,2.4},{1.3,1.3}}),1.0);
    ARIADNE_TEST_EQUAL(x.derivative(1).evaluate(Vector<Float>{2.4, 1.3}),0.0);

}


int main() {
    TestFunction().test();

    return ARIADNE_TEST_FAILURES;
}

