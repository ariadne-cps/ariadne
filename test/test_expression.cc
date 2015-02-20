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
#include "utility/container.h"
#include "utility/stlio.h"
#include "numeric/numeric.h"
#include "expression/expression.h"
#include "expression/assignment.h"
#include "expression/valuation.h"
#include "expression/space.h"
#include "function/function.h"

#include "test.h"

using namespace Ariadne;

class TestExpression {
    RealConstant o;
    RealVariable x,y,z;
  public:
    TestExpression()
        : o("1.0",1.0_q), x("x"), y("y"), z("z") {
    }

    Void test_variables() {
        ARIADNE_TEST_CONSTRUCT(RealVariable,a,("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable(a));
        ARIADNE_TEST_ASSERT(a!=RealVariable("b"));
    }

    Void test_assignment() {
        RealVariable x("x"), y("y");
        Real z(0), o(1);
        Assignment<Variable<Real>,Real> ac = (x=o);
        Assignment<Variable<Real>,Expression<Real> > a = (x=o);
        Valuation<Real> v(ac);

        List< Assignment<Variable<Real>,Real> > lac;
        lac = {x=z,y=o};
        List< Assignment<Variable<Real>,Expression<Real> > > la;
        la = {x=z,y=o};
        ARIADNE_TEST_PRINT( la );
        v = Valuation<Real>(lac);
        v = lac;
    }

    Void test_evaluate() {
        ARIADNE_TEST_CONSTRUCT(RealExpression,g,(x+3*y*z*z));

        Map<RealVariable,Real> v;
        v[x]=Real(2.0); v[y]=Real(3.0); v[z]=Real(5.0);

        ARIADNE_TEST_PRINT(v);
        //ARIADNE_TEST_EQUAL(evaluate(g,v),Real(227));
    }

    Void test_parameters() {
        RealVariable x("x");

        RealExpression expr = x;//+u;

        Map<Identifier,Real> valuation;
        Real uncertain_value = numeric_cast<Real>(ValidatedNumber(-1.0,1.0));
        valuation[x.name()] = uncertain_value;

        ARIADNE_TEST_EQUALS(expr.kind(),OperatorKind::VARIABLE);
        ARIADNE_TEST_EQUALS(expr.var(),"x");
        ARIADNE_TEST_EQUALS(valuation[x.name()],uncertain_value);

        Real result1 = evaluate(expr,valuation);

        ARIADNE_TEST_EQUALS(result1,uncertain_value);
    }

    Void test_function()
    {
        // Test to ensure that constants are handled correctly.
        Real tc=Dyadic(5.0);
        Real tx=Dyadic(1.125);
        Real ty=Dyadic(2.375);
        Real tz=Dyadic(3.750);

        Vector<Real> tv={tx,ty,tz};

        RealConstant c("5",tc);
        RealVariable x("x");
        RealVariable y("y");
        RealVariable z("z");

        RealExpression e1=c;
        EffectiveScalarFunction f1=make_function(e1,RealSpace({x,y,z}));
        ARIADNE_TEST_PRINT(f1);
        ARIADNE_TEST_EQUAL(f1.evaluate(tv), tc);

        RealExpression e2=c+x;
        EffectiveScalarFunction f2=make_function(e2,RealSpace({x,y,z}));
        ARIADNE_TEST_PRINT(f2);
        ARIADNE_TEST_EQUAL(f2.evaluate(tv), tc+tx);

        RealExpression e3=c+x+c*y;
        EffectiveScalarFunction f3=make_function(e3,{x,y,z});
        ARIADNE_TEST_PRINT(f3);
        ARIADNE_TEST_EQUAL(f3.evaluate(tv), tc+tx+tc*ty);

        RealExpression e4=exp(c+x);
        EffectiveScalarFunction f4=make_function(e4,{x,y,z});
        ARIADNE_TEST_PRINT(f4);
        ARIADNE_TEST_ASSERT(possibly((f4.evaluate(tv) == exp(tc+tx)).check(Effort(0))));

        //ARIADNE_TEST_EVALUATE(EffectiveVectorFunction((x+y,y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(EffectiveVectorFunction((x+y,y+z*z),(x,y,z))[0],EffectiveScalarFunction(x+y,(x,y,z)));

        //ARIADNE_TEST_EVALUATE(EffectiveVectorFunction((dot(x),dot(y)),(dot(x)=x+y,dot(y)=y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(EffectiveVectorFunction((x+y,y+z*z),(x,y,z))[0],EffectiveScalarFunction(x+y,(x,y,z)));
    }

    Void test() {
        test_variables();
        test_assignment();
        test_parameters();
        test_function();
    }
};


Int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
};
