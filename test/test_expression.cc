/***************************************************************************
 *            test_expression.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include <iostream>

#include "container.h"
#include "stlio.h"
#include "numeric.h"
#include "expression.h"
#include "assignment.h"
#include "space.h"
#include "function.h"

#include "test.h"

using namespace Ariadne;

class TestExpression {
    RealConstant o;
    RealVariable x,y,z;
  public:
    TestExpression()
        : o("1.0",1.0), x("x"), y("y"), z("z") {
    }

    void test_variables() {
        ARIADNE_TEST_CONSTRUCT(RealVariable,a,("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable(a));
        ARIADNE_TEST_ASSERT(a!=RealVariable("b"));
    }

    void test_evaluate() {
        ARIADNE_TEST_CONSTRUCT(RealExpression,g,(x+3*y*z*z));

        Map<RealVariable,Real> v;
        v[x]=2.0; v[y]=3.0; v[z]=5.0;

        ARIADNE_TEST_PRINT(v);
        //ARIADNE_TEST_EQUAL(evaluate(g,v),Real(227));
    }

    void test_derivative() {
        ARIADNE_TEST_CONSTRUCT(RealExpression,g,(x+3*y*z*z));
        ARIADNE_TEST_PRINT(derivative(g,y));
        ARIADNE_TEST_PRINT(derivative(g,z));
        //ARIADNE_TEST_EQUAL(derivative(g,y),3*z*z);
        //ARIADNE_TEST_EQUAL(derivative(g,z),6*y*z);
    }

    void test_function()
    {
        // Test to ensure that constants are handled correctly.
        Real tc=5.0;
        Real tx=1.125;
        Real ty=2.375;
        Real tz=3.750;

        Vector<Real> tv=Vector<Real>((tx,ty,tz));

        RealConstant c("5",5.0);
        RealVariable x("x");
        RealVariable y("y");
        RealVariable z("z");

        RealExpression e1=c;
        RealScalarFunction f1=make_function(e1,RealSpace((x,y,z)));
        ARIADNE_TEST_PRINT(f1);
        ARIADNE_TEST_EQUAL(f1.evaluate(tv), tc);

        RealExpression e2=c+x;
        RealScalarFunction f2=make_function(e2,RealSpace((x,y,z)));
        ARIADNE_TEST_PRINT(f2);
        ARIADNE_TEST_EQUAL(f2.evaluate(tv), tc+tx);

        RealExpression e3=c+x+c*y;
        RealScalarFunction f3=make_function(e3,(x,y,z));
        ARIADNE_TEST_PRINT(f3);
        ARIADNE_TEST_EQUAL(f3.evaluate(tv), tc+tx+tc*ty);

        //ARIADNE_TEST_EVALUATE(RealVectorFunction((x+y,y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(RealVectorFunction((x+y,y+z*z),(x,y,z))[0],RealScalarFunction(x+y,(x,y,z)));

        //ARIADNE_TEST_EVALUATE(RealVectorFunction((dot(x),dot(y)),(dot(x)=x+y,dot(y)=y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(RealVectorFunction((x+y,y+z*z),(x,y,z))[0],RealScalarFunction(x+y,(x,y,z)));
    }

    void test() {
        test_variables();
        test_function();
    }
};


int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
};
