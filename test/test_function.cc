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

#include "taylor_model.h"
#include "expression.h"
#include "assignment.h"
#include "function.h"
#include "predicate.h"

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
    void test_expression();
    void test_conversions();
    void test_differentiation();
};

void TestFunction::test()
{
    ARIADNE_TEST_CALL(test_scalar_function());
    ARIADNE_TEST_CALL(test_vector_function());
    ARIADNE_TEST_CALL(test_expression());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_differentiation());
}

void TestFunction::test_concept()
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

    //Vector<Float> b; Matrix<Float> A;
    //VectorAffineFunction aff(A,b);

}

void TestFunction::test_scalar_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(ScalarFunction,o,constant(3,1.0));
    ARIADNE_TEST_NAMED_CONSTRUCT(ScalarFunction,x,coordinate(3,0));
    ARIADNE_TEST_NAMED_CONSTRUCT(ScalarFunction,y,coordinate(3,1));

    ARIADNE_TEST_CONSTRUCT(ScalarFunction,f,(o+x*y));
    ARIADNE_TEST_CONSTRUCT(Vector<Float>,p,(3,2.0,3.0,5.0));
    ARIADNE_TEST_EQUAL(f(p),7.0);

    ARIADNE_TEST_PRINT(cos(f));

    ScalarFunction df=f.derivative(1);
    ARIADNE_TEST_PRINT(df);
    ARIADNE_TEST_EQUAL(df(p),2.0);
}

void TestFunction::test_vector_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(VectorFunction,f,identity(3));
    ARIADNE_TEST_EQUAL(f[0](Vector<Float>(3, 2.0,3.0,5.0)),2.0);
    ARIADNE_TEST_EXECUTE(f[0]=f[1]);
    ARIADNE_TEST_EQUAL(f[0](Vector<Float>(3, 2.0,3.0,5.0)),3.0);
}


void TestFunction::test_expression()
{
    // Test to ensure that constants are handled correctly.
    TaylorModel tc=TaylorModel::constant(3,5.0);
    TaylorModel tx=TaylorModel::variable(3,0);
    TaylorModel ty=TaylorModel::variable(3,1);
    TaylorModel tz=TaylorModel::variable(3,2);

    Vector<TaylorModel> tv=TaylorModel::variables(3);

    RealConstant c("5",5.0);
    RealVariable x("x");
    RealVariable y("y");
    RealVariable z("z");

    RealExpression e1=c;
    ScalarFunction f1(e1,(x,y,z));
    ARIADNE_TEST_PRINT(f1);
    ARIADNE_TEST_EQUAL(f1.evaluate(tv), tc);

    RealExpression e2=c+x;
    ScalarFunction f2(e2,(x,y,z));
    ARIADNE_TEST_PRINT(f2);
    ARIADNE_TEST_EQUAL(f2.evaluate(tv), tc+tx);

    RealExpression e3=c+x+c*y;
    ScalarFunction f3(e3,(x,y,z));
    ARIADNE_TEST_PRINT(f3);
    ARIADNE_TEST_EQUAL(f3.evaluate(tv), tc+tx+tc*ty);

    ScalarFunction df3=f3.derivative(1);
    ARIADNE_TEST_EQUAL(df3.evaluate(tv), tc);

    ARIADNE_TEST_EVALUATE(ScalarFunction(c+x+2*y*z*z,(x,y,z)).derivative(1));
    //ARIADNE_TEST_EQUAL(ScalarFunction(c+x+2*y*z*z,(x,y,z)).derivative(1),ScalarFunction(2*z*z,(x,y,z)));

    ARIADNE_TEST_EVALUATE(VectorFunction((x+y,y+z*z),(x,y,z))[0]);
    //ARIADNE_TEST_EQUAL(VectorFunction((x+y,y+z*z),(x,y,z))[0],ScalarFunction(x+y,(x,y,z)));

    ARIADNE_TEST_EVALUATE(VectorFunction((dot(x),dot(y)),(dot(x)=x+y,dot(y)=y+z*z),(x,y,z))[0]);
    //ARIADNE_TEST_EQUAL(VectorFunction((x+y,y+z*z),(x,y,z))[0],ScalarFunction(x+y,(x,y,z)));
}




void TestFunction::test_conversions()
{

}

void TestFunction::test_differentiation()
{
    ScalarFunction z=ScalarFunction::constant(2,0.0);
    ScalarFunction o=ScalarFunction::constant(2,1.0);
    ScalarFunction x=ScalarFunction::coordinate(2,0);
    ScalarFunction y=ScalarFunction::coordinate(2,1);

    ScalarAffineFunction af(Vector<Real>(2, 3.0, -2.0), 1.0);
    ScalarFunction daf=af.derivative(1);
    ARIADNE_TEST_EQUAL(daf.evaluate(Vector<Float>(2,2.4,1.3)),-2.0);

    ARIADNE_TEST_EQUAL(x.derivative(0).evaluate(Vector<Float>(2,2.4,1.3)),1.0);
    ARIADNE_TEST_EQUAL(x.derivative(0).evaluate(Vector<Interval>(2, 2.4,2.4, 1.3,1.3)),1.0);
    ARIADNE_TEST_EQUAL(x.derivative(1).evaluate(Vector<Float>(2, 2.4, 1.3)),0.0);

}


int main() {
    TestFunction().test();

    return ARIADNE_TEST_FAILURES;
}

