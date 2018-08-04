/***************************************************************************
 *            test_expression.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "expression/expression.hpp"
#include "expression/assignment.hpp"
#include "expression/valuation.hpp"
#include "expression/space.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"

#include "test.hpp"

using namespace Ariadne;

typedef Algebra<EffectiveNumericType> EffectiveAlgebra;
typedef SymbolicAlgebra<EffectiveNumericType> EffectiveSymbolicAlgebra;

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
        Real zero(0), one(1);
        RealVariable x("x"), y("y"), z("x");
        RealExpression e(x*y+o);

        ARIADNE_TEST_ASSERT((not IsAssignable<Variable<Real>,Expression<Real>>::value));
        ARIADNE_TEST_ASSERT((not IsConstructible<Assignment<Variable<Real>,Expression<Real>>,Assignment<Variable<Real>,Real>>::value));

        typedef Assignment<Variable<Real>,Real> ConstantRealAssignment;
        ARIADNE_TEST_CONSTRUCT(ConstantRealAssignment,ac,(x=one));
        ARIADNE_TEST_CONSTRUCT(List<ConstantRealAssignment>,lac,({x=zero,y=one}));
        ARIADNE_TEST_CONSTRUCT(Valuation<Real>,va,(lac));
        //ARIADNE_TEST_CONSTRUCT(Valuation<Real>,va,({x=zero,y=one})); Fails due to ambiguous overload

        ARIADNE_TEST_CONSTRUCT(RealAssignment,a,(let(x)=one));
        ARIADNE_TEST_CONSTRUCT(PrimedRealAssignment,pa,(prime(x)=one));
        ARIADNE_TEST_CONSTRUCT(DottedRealAssignment,da,(dot(x)=one));

        ARIADNE_TEST_CONSTRUCT(List<RealAssignment>,la,(let({x,y,z})={zero,x,e}));
        ARIADNE_TEST_CONSTRUCT(List<PrimedRealAssignment>,lpa,(prime({x,y,z})={zero,x,e}));
        ARIADNE_TEST_CONSTRUCT(List<DottedRealAssignment>,lda,(dot({x,y,z})={zero,x,e}));
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
        Real value = Real(ExactNumericType(-0.0626));
        valuation[x.name()] = value;

        ARIADNE_TEST_EQUALS(expr.kind(),OperatorKind::VARIABLE);
        ARIADNE_TEST_EQUALS(expr.var(),"x");
        ARIADNE_TEST_EQUALS(valuation[x.name()],value);

        Real result1 = evaluate(expr,valuation);

        ARIADNE_TEST_EQUALS(result1,value);
    }

    Void test_simplify() {

        RealVariable x("x"), y("y"), u("u");
        RealExpression expr = -u*x*y+2*x;
        RealExpression simplification = simplify(derivative(expr,x));
        ARIADNE_TEST_ASSERT(identical(simplification,-u*y+2));
    }

    Void test_scalar_properties()
    {
        RealVariable x("x"), y("y");
        Real c(3);
        ARIADNE_TEST_ASSERT(is_constant_in(3*y,{x}));
        ARIADNE_TEST_ASSERT(not is_constant_in(3*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in(0*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in((sin(2*c)-2*sin(c)*cos(c))*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in((sin(2*x)-2*sin(x)*cos(x))*y,{y}));
        ARIADNE_TEST_ASSERT(is_constant_in(simplify(0*y),{y}));

        ARIADNE_TEST_ASSERT(is_affine_in(2+3*x-5*y-x,{x,y}));
        ARIADNE_TEST_ASSERT(is_affine_in(3*y,{x,y}));
        ARIADNE_TEST_ASSERT(is_affine_in(x*y,{x}));
        ARIADNE_TEST_ASSERT(is_affine_in(3*x/y,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*y,{x,y}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*x,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(0*x*x,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x/y,{y}));
    }

    Void test_vector_properties()
    {
        RealVariable x("x"), y("y"), u1("u1"), u2("u2");
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u1,y+u2}),{u1,u2}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u2,y+u1}),{u1,u2}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x+u1,y}),{u1}));
        ARIADNE_TEST_ASSERT(is_additive_in(Vector<RealExpression>({x,y+u1}),{u1}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x+u1,y+u1}),{u1}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x,y+2*u1}),{u1}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x+u1,y+2*u2}),{u1,u2}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x*u1,y+u2}),{u1,u2}));
        ARIADNE_TEST_ASSERT(not is_additive_in(Vector<RealExpression>({x+u1,y+sqr(u2)}),{u1,u2}));
    }

    Void test_function()
    {
        // Test to ensure that constants are handled correctly.
        Real tc=Dyadic(5.0);
        Real tx=Dyadic(1.125);
        Real ty=Dyadic(2.375);
        Real tz=Dyadic(3.750);

        Vector<Real> tv={tx,ty,tz};
        Valuation<Real> tw({{x,tx},{y,ty},{z,tz}});

        RealConstant c("5",tc);
        RealVariable x("x");
        RealVariable y("y");
        RealVariable z("z");

        RealExpression e1=c;
        EffectiveScalarFunction f1=make_function(e1,RealSpace(List<RealVariable>({x,y,z})));
        ARIADNE_TEST_PRINT(f1);
        ARIADNE_TEST_EQUAL(f1.evaluate(tv), tc);

        RealExpression e2=c+x;
        EffectiveScalarFunction f2=make_function(e2,RealSpace(List<RealVariable>({x,y,z})));
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

        EffectiveAlgebra ax=RealExpression(x);
        EffectiveAlgebra ay=RealExpression(y);
        EffectiveAlgebra az=RealExpression(z);
        Vector<EffectiveAlgebra> va={ax,ay,az};
        ARIADNE_TEST_PRINT(f3(va));
        ARIADNE_TEST_PRINT(f3(va).extract<RealExpression>());
        ARIADNE_TEST_EQUALS(evaluate(f3(va).extract<RealExpression>(),tw),evaluate(e3,tw));
    }

    Void test() {
        test_variables();
        test_assignment();
        test_parameters();
        test_simplify();
        test_scalar_properties();
        test_vector_properties();
        test_function();
    }
};


Int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
};
