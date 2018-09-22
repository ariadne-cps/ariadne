/***************************************************************************
 *            test_function.cpp
 *
 *  Copyright  2009  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/symbolic_function.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"

#include "numeric/numeric.hpp"

#include "../test.hpp"

using namespace std;
using namespace Ariadne;

class TestFunction
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_scalar_univariate_function();
    Void test_vector_univariate_function();
    Void test_scalar_function();
    Void test_vector_function();
    Void test_conversions();
    Void test_differentiation();
};

Void TestFunction::test()
{
    ARIADNE_TEST_CALL(test_scalar_univariate_function());
    ARIADNE_TEST_CALL(test_vector_univariate_function());
    ARIADNE_TEST_CALL(test_scalar_function());
    ARIADNE_TEST_CALL(test_vector_function());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_differentiation());
}

Void TestFunction::test_concept()
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

    //Vector<FloatDPApproximation> b; Matrix<FloatDPApproximation> A;
    //VectorAffineFunction aff(A,b);

}

Void TestFunction::test_scalar_univariate_function()
{
    ARIADNE_TEST_CONSTRUCT(FloatDPApproximation,p,({3.0_approx}));

    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarUnivariateFunction,o,constant(1));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarUnivariateFunction,x,coordinate());

    ARIADNE_TEST_EQUAL(o(p),1.0_approx);
    ARIADNE_TEST_EQUAL(x(p),3.0_approx);

    ARIADNE_TEST_CONSTRUCT(EffectiveScalarUnivariateFunction,f,(o+x*x));
    ARIADNE_TEST_EQUAL(f(p),10.0_approx);

    ARIADNE_TEST_PRINT(cos(f));

    EffectiveScalarUnivariateFunction df=derivative(f);
    ARIADNE_TEST_PRINT(df);
    ARIADNE_TEST_EQUAL(df(p),6.0_approx);
}

Void TestFunction::test_vector_univariate_function()
{
    Vector<Real> c({7,11});

    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveVectorUnivariateFunction,f,constant(c));

    FloatDPApproximation x=2.0_approx;

    ARIADNE_TEST_EQUAL(f(x)[0],c[0]);
    ARIADNE_TEST_EQUAL(f[0](x),c[0]);
    ARIADNE_TEST_EQUAL(f(x)[1],c[1]);
    ARIADNE_TEST_EQUAL(f[1](x),c[1]);
}

Void TestFunction::test_scalar_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,o,constant(3,1));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,x,coordinate(3,0));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,y,coordinate(3,1));


    ARIADNE_TEST_CONSTRUCT(EffectiveScalarFunction,f,(o+x*y));
    ARIADNE_TEST_CONSTRUCT(Vector<FloatDPApproximation>,p,({2.0_approx,3.0_approx,5.0_approx}));
    ARIADNE_TEST_EQUAL(f(p),7.0_approx);

    ARIADNE_TEST_PRINT(cos(f));

    EffectiveScalarFunction df=f.derivative(1);
    ARIADNE_TEST_PRINT(df);
    ARIADNE_TEST_EQUAL(df(p),2.0_approx);
}

Void TestFunction::test_vector_function()
{
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveVectorFunction,id,identity(3));

    // Regression tests for element
    ARIADNE_TEST_PRINT(id[0]);
    const EffectiveVectorFunction& id_const_ref=id;
    ARIADNE_TEST_PRINT(id_const_ref[0]);
    EffectiveVectorFunction& id_ref=id;
    ARIADNE_TEST_PRINT(id_ref[0]);

    Vector<FloatDPApproximation> v={2.0_approx,3.0_approx,5.0_approx};

    ARIADNE_TEST_CONSTRUCT(EffectiveVectorFunction,f,(id));
    ARIADNE_TEST_EQUAL(f(v),v);
    ARIADNE_TEST_EQUAL(f[0](v),v[0]);
    ARIADNE_TEST_EXECUTE(f[0]=f[1]);
    ARIADNE_TEST_EQUAL(f[0](v),v[1]);
    ARIADNE_TEST_EQUAL(f(v)[0],v[1]);

    // Regression test for function - vector product
    EffectiveScalarFunction x1=EffectiveScalarFunction::coordinate(3,1);
//    EffectiveVector e0={1,0};
    Vector<EffectiveNumber> e0={1,0};
    ARIADNE_TEST_EQUAL((x1*e0)(v)[0],v[1]);
    ARIADNE_TEST_EQUAL((x1*e0)(v)[1],0);
}





Void TestFunction::test_conversions()
{

}

Void TestFunction::test_differentiation()
{
    EuclideanDomain dom(2);
    DegreeType deg(3);

    EffectiveScalarFunction z=EffectiveScalarFunction::constant(dom,0);
    EffectiveScalarFunction o=EffectiveScalarFunction::constant(dom,1);
    EffectiveScalarFunction x=EffectiveScalarFunction::coordinate(dom,0);
    EffectiveScalarFunction y=EffectiveScalarFunction::coordinate(dom,1);

    ARIADNE_TEST_CONSTRUCT(EffectiveScalarFunction,f,(3*x-2*y+1));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,df0,constant(dom,3));
    ARIADNE_TEST_NAMED_CONSTRUCT(EffectiveScalarFunction,df1,constant(dom,-2));

    ARIADNE_TEST_CONSTRUCT(Vector<FloatDPBounds>,v,({2.25_exact,1.375_exact}));
    ARIADNE_TEST_CONSTRUCT(Vector<FloatDPApproximation>,va,(v));

    ARIADNE_TEST_EVALUATE(f.differential(v,deg));
    ARIADNE_TEST_EVALUATE(f.differential(va,deg));

    ARIADNE_TEST_EQUAL(f.differential(v,deg).value(),f.evaluate(v));
    ARIADNE_TEST_EQUAL(f.differential(v,deg).gradient()[0],df0.evaluate(v));
    ARIADNE_TEST_EQUAL(f.differential(v,deg).gradient()[1],df1.evaluate(v));
    ARIADNE_TEST_EQUAL(differential(f,v,deg),f.differential(v,deg));
    ARIADNE_TEST_EQUAL(gradient(f,v),f.differential(v,deg).gradient());
    ARIADNE_TEST_EQUAL(f.gradient(v),gradient(f,v));

    ARIADNE_TEST_EQUAL(f.derivative(1).evaluate(v),df1.evaluate(v));
    ARIADNE_TEST_EQUAL(derivative(f,1).evaluate(v),df1.evaluate(v));

    ARIADNE_TEST_EQUAL(x.derivative(0).evaluate(v),1.0_exact);
    ARIADNE_TEST_EQUAL(x.derivative(1).evaluate(v),0.0_exact);

}


Int main() {
    TestFunction().test();

    return ARIADNE_TEST_FAILURES;
}

