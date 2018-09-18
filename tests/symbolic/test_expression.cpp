/***************************************************************************
 *            test_expression.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include <iostream>

#include "config.hpp"
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/valuation.hpp"
#include "symbolic/space.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"

#include "algebra/matrix.tpl.hpp"

#include "../test.hpp"

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

    Void test_expression() {
        // Regression test for constructing Expression from 0 without being an ambiguous nullptr;
        ARIADNE_TEST_CONSTRUCT(IntegerExpression,ze,(0));
        ARIADNE_TEST_CONSTRUCT(RealExpression,re,(0));
        RealExpression(0);
        RealExpression(nullptr);
    }

    Void test_assignment() {
        Real zero(0), one(1);
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

    Void test_derivative() {
        RealExpression expr = 2*x+y;
        ARIADNE_TEST_ASSERT(identical(simplify(derivative(expr,x)),RealExpression::constant(2)));
        RealExpression expr2 = pow(x,3);
        ARIADNE_TEST_ASSERT(identical(simplify(derivative(expr2,x)),3*pow(x,2)));
    }

    Void test_simplify() {

        RealVariable u("u");
        RealExpression expr = -u*x*y+2*x;
        RealExpression simplification = simplify(derivative(expr,x));
        ARIADNE_TEST_ASSERT(identical(simplification,-u*y+2));
    }

    Void test_substitute() {

        RealVariable u1("u1"), u2("u2");
        RealExpression expr = -u1*x*y+2*pow(x+u2,2);

        List<Assignment<RealVariable,RealExpression>> subs={{u1,u1+1},{u2,u1*x}};

        RealExpression substitution = substitute(expr,subs);

        ARIADNE_TEST_ASSERT(identical(substitution,-(u1+1)*x*y+2*pow(x+u1*x,2)));
    }

    Void test_scalar_properties()
    {
        RealVariable u("u");
        Real c(3);
        ARIADNE_TEST_ASSERT(is_constant_in(3*y,{x}));
        ARIADNE_TEST_ASSERT(is_constant_in(pow(x,2),{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in(pow(x,2),{x}));
        ARIADNE_TEST_ASSERT(not is_constant_in(3*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in(0*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in((sin(2*c)-2*sin(c)*cos(c))*y,{y}));
        ARIADNE_TEST_ASSERT(not is_constant_in((sin(2*x)-2*sin(x)*cos(x))*y,{y}));
        ARIADNE_TEST_ASSERT(is_constant_in(simplify(0*y),{y}));

        ARIADNE_TEST_ASSERT(is_affine_in(sqr(x),{y}));
        ARIADNE_TEST_ASSERT(is_affine_in(pow(x,3),{y}));
        ARIADNE_TEST_ASSERT(is_affine_in(pow(x,3)+y,{y}));
        ARIADNE_TEST_ASSERT(is_affine_in(2+3*x-5*y-x,{x,y}));
        ARIADNE_TEST_ASSERT(is_affine_in(3*y,{x,y}));
        ARIADNE_TEST_ASSERT(is_affine_in(x*y,{x}));
        ARIADNE_TEST_ASSERT(is_affine_in(3*x/y,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(pow(x,3),{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(sqr(x),{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*y,{x,y}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x*x,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(0*x*x,{x}));
        ARIADNE_TEST_ASSERT(not is_affine_in(x/y,{y}));
    }

    Void test_vector_properties()
    {
        RealVariable u1("u1"), u2("u2");
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

        RealExpression e1=c;
        EffectiveScalarMultivariateFunction f1=make_function(e1,RealSpace(List<RealVariable>({x,y,z})));
        ARIADNE_TEST_PRINT(f1);
        ARIADNE_TEST_EQUAL(f1.evaluate(tv), tc);

        RealExpression e2=c+x;
        EffectiveScalarMultivariateFunction f2=make_function(e2,RealSpace(List<RealVariable>({x,y,z})));
        ARIADNE_TEST_PRINT(f2);
        ARIADNE_TEST_EQUAL(f2.evaluate(tv), tc+tx);

        RealExpression e3=c+x+c*y;
        EffectiveScalarMultivariateFunction f3=make_function(e3,{x,y,z});
        ARIADNE_TEST_PRINT(f3);
        ARIADNE_TEST_EQUAL(f3.evaluate(tv), tc+tx+tc*ty);

        RealExpression e4=exp(c+x);
        EffectiveScalarMultivariateFunction f4=make_function(e4,{x,y,z});
        ARIADNE_TEST_PRINT(f4);
        ARIADNE_TEST_ASSERT(possibly((f4.evaluate(tv) == exp(tc+tx)).check(Effort(0))));

        //ARIADNE_TEST_EVALUATE(EffectiveVectorMultivariateFunction((x+y,y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(EffectiveVectorMultivariateFunction((x+y,y+z*z),(x,y,z))[0],EffectiveScalarMultivariateFunction(x+y,(x,y,z)));

        //ARIADNE_TEST_EVALUATE(EffectiveVectorMultivariateFunction((dot(x),dot(y)),(dot(x)=x+y,dot(y)=y+z*z),(x,y,z))[0]);
        //ARIADNE_TEST_EQUAL(EffectiveVectorMultivariateFunction((x+y,y+z*z),(x,y,z))[0],EffectiveScalarMultivariateFunction(x+y,(x,y,z)));

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
        test_expression();
        test_assignment();
        test_parameters();
        test_derivative();
        test_simplify();
        test_substitute();
        test_scalar_properties();
        test_vector_properties();
        test_function();
    }

};


Int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
}
