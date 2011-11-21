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

#include "config.h"
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

    void test_parameters() {
    	RealVariable x("x");

    	RealExpression expr = x;//+u;

    	Map<Identifier,Real> valuation;
    	valuation[x.name()] = Real(-1.0,1.0);

    	ARIADNE_TEST_EQUALS(expr.kind(),VARIABLE);
    	ARIADNE_TEST_EQUALS(expr.var(),"x");
    	ARIADNE_TEST_EQUALS(valuation[x.name()],Real(-1.0,1.0));

    	Real result1 = evaluate(expr,valuation);

    	ARIADNE_TEST_EQUALS(result1,Real(-1.0,1.0));
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

        RealExpression e4=exp(c+x);
        RealScalarFunction f4=make_function(e4,(x,y,z));
        ARIADNE_TEST_PRINT(f4);
        ARIADNE_TEST_EQUAL(f4.evaluate(tv), exp(tc+tx));

        //ARIADNE_TEST_EVALUATE(RealVectorFunction((x+y,y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(RealVectorFunction((x+y,y+z*z),(x,y,z))[0],RealScalarFunction(x+y,(x,y,z)));

        //ARIADNE_TEST_EVALUATE(RealVectorFunction((dot(x),dot(y)),(dot(x)=x+y,dot(y)=y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(RealVectorFunction((x+y,y+z*z),(x,y,z))[0],RealScalarFunction(x+y,(x,y,z)));
    }

    void test() {
        test_variables();
        test_parameters();
        test_function();
    }
};


int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
};
