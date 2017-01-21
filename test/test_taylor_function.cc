/***************************************************************************
 *            test_taylor_function.cc
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
#include <iomanip>
#include "config.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/covector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "algebra/expansion.h"
#include "algebra/algebra.h"
#include "function/taylor_model.h"
#include "function/taylor_function.h"
#include "algebra/differential.h"
#include "function/polynomial.h"
#include "function/function.h"

#include "function/formula.h"
#include "function/symbolic_function.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

Vector<Real> e(Nat n, Nat i) { return Vector<Real>::unit(n,i); }
Polynomial<Float64> p(Nat n, Nat j) { return Polynomial<Float64>::variable(n,j); }
ScalarTaylorFunction t(ExactBoxType d, Nat j,Sweeper<Float64> swp) { return ScalarTaylorFunction::coordinate(d,j,swp); }

template<class X> Vector< Expansion<X> > operator*(const Expansion<X>& e, const Vector<Float64> v) {
    Vector< Expansion<X> > r(v.size(),Expansion<X>(e.argument_size()));
    for(Nat i=0; i!=r.size(); ++i) { ARIADNE_ASSERT(v[i]==0.0 || v[i]==1.0); if(v[i]==1.0) { r[i]=e; } }
    return r;
}

class TestScalarTaylorFunction
{
    Precision64 pr;
    Sweeper<Float64> swp;
  public:
    TestScalarTaylorFunction(Sweeper<Float64> sweeper);
    Void test();
  private:
    Void test_concept();
    Void test_create();
    Void test_constructors();
    Void test_predicates();
    Void test_approximation();
    Void test_evaluate();
    Void test_gradient();
    Void test_arithmetic();
    Void test_functions();
    Void test_compose();
    Void test_antiderivative();
    Void test_conversion();
    Void test_generic();
  private:
    ExactBoxType d(SizeType n) { return Vector<ExactIntervalType>(n,ExactIntervalType(-1,+1)); }
    typedef Expansion<RawFloat64> e;
    typedef TaylorModel<ValidatedTag,Float64> TM;
};


TestScalarTaylorFunction::TestScalarTaylorFunction(Sweeper<Float64> sweeper)
    : swp(sweeper)
{
}


Void TestScalarTaylorFunction::test()
{
    Float64Value::set_output_places(18);
    Float64Error::set_output_places(18);
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_create());
    ARIADNE_TEST_CALL(test_predicates());
    ARIADNE_TEST_CALL(test_approximation());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_functions());
    ARIADNE_TEST_CALL(test_evaluate());
    ARIADNE_TEST_CALL(test_gradient());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_conversion());
    ARIADNE_TEST_CALL(test_generic());
}


Void TestScalarTaylorFunction::test_concept()
{
    const Float64Value f={0,pr};
    const Float64Bounds i;
    const Vector<Float64Value> vf;
    const Vector<Float64Bounds> vi;
    const ScalarTaylorFunction  t;
    ScalarTaylorFunction tr;

    tr=t+f; tr=t-f; tr=t*f; tr=t/f;
    tr=f+t; tr=f-t; tr=f*t; tr=f/t;
    tr=t+i; tr=t-i; tr=t*i; tr=t/i;
    tr=i+t; tr=i-t; tr=i*t; tr=i/t;
    tr=t+t; tr=t-t; tr=t*t; tr=t/t;

    tr+=f; tr-=f; tr*=f; tr/=f;
    tr+=i; tr-=i; tr*=i; tr/=i;
    tr+=t; tr-=t;

    tr=pos(tr); tr=neg(tr); tr=sqr(tr);
    tr=rec(tr); tr=pow(tr,1u); tr=pow(tr,-1);
    tr=exp(t); tr=log(t); tr=sqrt(t);
    tr=sin(t); tr=cos(t); tr=tan(t);
    //tr=asin(t); tr=acos(t); tr=atan(t);
    tr=max(tr,tr); tr=min(tr,tr); tr=abs(tr);

    tr.simplify(); tr.clobber();

    t(vi); evaluate(t,vi);
    t.domain(); t.range(); t.expansion(); t.error();

}

Void TestScalarTaylorFunction::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv1,({{-1,+1},{-1,+1}},{{{0,0},1.},{{1,0},2.},{{0,1},3.},{{2,0},4.},{{1,1},5.},{{0,2},6.},{{3,0},7.},{{2,1},8.},{{1,2},9.},{{0,3},10.}},0.25,swp));

    ARIADNE_ASSERT_EQUAL(tv1.domain(),Box<ExactIntervalType>({{-1,+1},{-1,+1}}));
    ARIADNE_ASSERT_EQUAL(tv1.argument_size(),2);
    ARIADNE_ASSERT_EQUAL(tv1.number_of_nonzeros(),10);
    ARIADNE_ASSERT_EQUAL(tv1.value().raw(),1.0);
    ARIADNE_ASSERT_EQUAL(tv1.error().raw(),0.25);
}

Void TestScalarTaylorFunction::test_create()
{
    ExactBoxType D={{-0.5,0.5},{-1.0,2.0}};
    Float64Value c1={1,pr};
    Float64Value c2={1,pr};
    ScalarTaylorFunction pf1=ScalarTaylorFunction::constant(D,c1,swp);
    ScalarTaylorFunction pf2=ScalarTaylorFunction::constant(D,c2,swp);
    ValidatedScalarFunctionModel fm1=pf1;
    ValidatedScalarFunctionModel fm2=pf2;
    ValidatedScalarFunction f2=fm2;
    ARIADNE_TEST_EQUAL(factory(pf1).create(pf2).range(),c2);
    ARIADNE_TEST_EQUAL(factory(pf1).create(fm2).range(),c2);
    ARIADNE_TEST_EQUAL(factory(pf1).create(f2).range(),c2);
    ARIADNE_TEST_EQUAL(factory(fm1).create(fm2).range(),c2);
    ARIADNE_TEST_EQUAL(factory(fm1).create(f2).range(),c2);
}

Void TestScalarTaylorFunction::test_predicates()
{
    ScalarTaylorFunction x=ScalarTaylorFunction::coordinate(d(1),0,swp);
    ScalarTaylorFunction e=ScalarTaylorFunction::unit_ball(d(1),swp);

    ScalarTaylorFunction tv1(d(1),{{{0},1.00},{{1},2.00},{{2},3.00}}, 0.75, swp);
    ScalarTaylorFunction tv2(d(1),{{{0},1.00},{{1},1.75},{{2},3.25}}, 0.25, swp);
    ScalarTaylorFunction tv3(d(1),{{{0},1.125},{{1},1.75},{{2},3.25}}, 0.25, swp);
    ScalarTaylorFunction tv4(d(1),{{{0},1.00},{{1},2.25},{{2},3.00},{{3},-0.25}}, 0.25, swp);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

Void TestScalarTaylorFunction::test_approximation()
{
//    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv1,(d(2),e(2,3,{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0}),0.25,swp));
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv2,(d(1),{{{0},1.0},{{1},2.0},{{2},3.0}},0.25,swp));
}

Void TestScalarTaylorFunction::test_evaluate()
{
    Vector<Float64Bounds> iv({{0.25_exact,0.5_exact},{-0.75_exact,-0.5_exact}});
    ScalarTaylorFunction tv(d(2),{{{0,0},1.0},{{1,0},2.0},{{0,1},3.0},{{2,0},4.0},{{1,1},5.0},{{0,2},6.0}},0.25,swp);
    ARIADNE_TEST_SAME(evaluate(tv,iv),Float64Bounds(-0.375000,3.43750));
}

Void TestScalarTaylorFunction::test_gradient()
{
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    Real a(1.5); Real b(0.25);
    EffectiveScalarFunction quadratic = a-x[0]*x[0]+b*x[1];

    ExactBoxType domain1={{-1.0,+1.0},{-1.0,+1.0}};
    ExactBoxType domain2={{-0.5,+0.5},{-0.25,+0.25}};
    ExactBoxType domain3={{-0.25,+0.75},{0.0,+0.50}};

    Vector<Float64Bounds> point1=Vector<Float64Bounds>{0.0_exact,0.0_exact};
    Vector<Float64Bounds> point2=Vector<Float64Bounds>{0.0_exact,0.25_exact};

    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(domain1,quadratic,swp).gradient(point1),quadratic.gradient(point1));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(domain1,quadratic,swp).gradient(point2),quadratic.gradient(point2));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(domain2,quadratic,swp).gradient(point1),quadratic.gradient(point1));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(domain2,quadratic,swp).gradient(point2),quadratic.gradient(point2));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(domain3,quadratic,swp).gradient(point1),quadratic.gradient(point1));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(domain3,quadratic,swp).gradient(point2),quadratic.gradient(point2));
}


Void TestScalarTaylorFunction::test_arithmetic()
{
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)+(-3), ScalarTaylorFunction(d(1),{{{0},-2.0},{{1},-2.0},{{2},3.0}}, 0.75, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)-(-3), ScalarTaylorFunction(d(1),{{{0},4.0},{{1},-2.0},{{2},3.0}}, 0.75, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)*(-3), ScalarTaylorFunction(d(1),{{{0},-3.0},{{1},6.0},{{2},-9.0}}, 2.25, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)/(-4), ScalarTaylorFunction(d(1),{{{0},-0.25},{{1},0.5},{{2},-0.75}}, 0.1875, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)+Float64Bounds(-1,2), ScalarTaylorFunction(d(1),{{{0},1.5},{{1},-2.0},{{2},3.0}}, 2.25, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)-Float64Bounds(-1,2), ScalarTaylorFunction(d(1),{{{0},0.5},{{1},-2.0},{{2},3.0}}, 2.25, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)*Float64Bounds(-1,2), ScalarTaylorFunction(d(1),{{{0},0.5},{{1},-1.0},{{2},1.5}}, 10.5, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)/Float64Bounds(0.25,2.0), ScalarTaylorFunction(d(1),{{{0},2.25},{{1},-4.5},{{2},6.75}}, 13.5, swp));
    ARIADNE_TEST_SAME(+ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp), ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp));
    ARIADNE_TEST_SAME(-ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp), ScalarTaylorFunction(d(1),{{{0},-1.0},{{1},2.0},{{2},-3.0}}, 0.75, swp));

    // Regression test to check subtraction yielding zero coefficients
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)+ScalarTaylorFunction(d(1),{{{0},3.0},{{1},2.0},{{2},-4.0}}, 0.5, swp), ScalarTaylorFunction(d(1),{{{0},4.0},{{1},0.0},{{2},-1.0}}, 1.25, swp));

    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)-ScalarTaylorFunction(d(1),{{{0},3.0},{{1},2.0},{{2},-4.0}}, 0.5, swp), ScalarTaylorFunction(d(1),{{{0},-2.0},{{1},-4.0},{{2},7.0}}, 1.25, swp));
    ARIADNE_TEST_SAME(ScalarTaylorFunction(d(1),{{{0},1.0},{{1},-2.0},{{2},3.0}}, 0.75, swp)*ScalarTaylorFunction(d(1),{{{0},3.0},{{1},2.0},{{2},-4.0}}, 0.5, swp),
                       ScalarTaylorFunction(d(1),{{{0},3.0},{{1},-4.0},{{2},1.0},{{3},14.0},{{4},-12.0}}, 10.125, swp));

}

Void TestScalarTaylorFunction::test_functions()
{
    ScalarTaylorFunction xz(d(1),{{{0},0.0},{{1}, 0.5}}, 0.0, swp);
    ScalarTaylorFunction xo(d(1),{{{0},1.0},{{1}, 0.5}}, 0.0, swp);

    //Functions based on their natural defining points
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),ScalarTaylorFunction(d(1),{{{0},1.00000},{{1},0.50000},{{2},0.12500},{{3},0.02083},{{4},0.00260},{{5},0.00026},{{6},0.00002}}, 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),ScalarTaylorFunction(d(1),{{{0},0.00000},{{1},0.50000},{{2},0.0000},{{3},-0.02083},{{4},0.00000},{{5},0.00026},{{6},0.00000}}, 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),ScalarTaylorFunction(d(1),{{{0},1.00000},{{1},0.0000},{{2},-0.12500},{{3},0.00000},{{4},0.00260},{{5},0.0000},{{6},-0.00002}}, 0.00003, swp));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),ScalarTaylorFunction(d(1),{{{0},1.000000},{{1},-0.500000},{{2},0.250000},{{3},-0.125000},{{4}, 0.062500},{{5},-0.031250},{{6},0.015625}}, 0.018, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),ScalarTaylorFunction(d(1),{{{0},1.000000},{{1}, 0.250000},{{2},-0.031250},{{3}, 0.007813},{{4},-0.002441},{{5}, 0.000854},{{6},-0.000320}}, 0.0003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),ScalarTaylorFunction(d(1),{{{0},0.000000},{{1},0.500000},{{2},-0.125000},{{3}, 0.041667},{{4},-0.015625},{{5}, 0.006250},{{6},-0.002604}}, 0.003, swp));

}


Void TestScalarTaylorFunction::test_compose()
{
}


Void TestScalarTaylorFunction::test_antiderivative()
{
    ScalarTaylorFunction tm=ScalarTaylorFunction::constant(d(2),1.0_exact,swp);
    ScalarTaylorFunction atm=antiderivative(tm,1u);

    ARIADNE_TEST_CONSTRUCT(ExactBoxType,dom,({{1.0, 4.0}, {0.0, 2.0}}));
    VectorTaylorFunction x = VectorTaylorFunction::identity(dom,swp);

    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction, f, (1 + 2*x[0]+3*x[1]+5*x[0]*x[0]+4*x[0]*x[1]));
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction, g, (antiderivative(f,0,2.0_exact)) );
    ARIADNE_TEST_PRINT( derivative(g,0) );
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction, dg, (derivative(g,0)) );
    ARIADNE_TEST_LESS(norm(dg-f),1e-8);

    // We should have f(c,s)=0 for all x1
    ScalarTaylorFunction s = ScalarTaylorFunction::coordinate({dom[1]},0u,swp);
    ScalarTaylorFunction c = ScalarTaylorFunction::constant(s.domain(),2.0_exact,swp);
    ScalarTaylorFunction h=compose(g,VectorTaylorFunction({c,s}));

    Vector<ExactIntervalType> hdom=h.domain();
    Vector<Float64Bounds> domv(reinterpret_cast<Vector<Float64Bounds>const&>(hdom));
    ARIADNE_ASSERT(definitely(mag(h(domv))<Float64Value(1e-8)));

}

Void TestScalarTaylorFunction::test_conversion() {
    {
        SizeType as=2;
        ValidatedTaylorModel64 ty=ValidatedTaylorModel64::coordinate(as,1,swp);
        ValidatedTaylorModel64 tc=ValidatedTaylorModel64::constant(as, 2, swp);
        ARIADNE_TEST_PRINT(ty);
        ARIADNE_TEST_PRINT(tc);
        ARIADNE_TEST_PRINT(rec(tc));
        return;
        ARIADNE_TEST_PRINT(ty/tc);
        ARIADNE_TEST_PRINT(rec(tc+ty));

    }


    // Test conversion between ordinary functions and Taylor functions.
    ExactBoxType D={{-0.5,0.5},{-1.0,2.0}};
    Vector<Float64Value> pt={-0.25_exact,0.25_exact};
    Vector<Float64Bounds> ipt(pt);
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,tx,(D,x,swp));

    ARIADNE_TEST_CONSTRUCT(EffectiveScalarFunction,f,(1-x[0]*x[0]-x[1]/2));
    ARIADNE_TEST_PRINT(compose(f,tx));
    ARIADNE_TEST_CONSTRUCT(EffectiveScalarFunction,y,(x[1]/2));
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,ty,(D,x[1]/2,swp));
    ARIADNE_TEST_PRINT(ScalarTaylorFunction(D,x[1],swp)/2);
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tf,(D,f,swp));

    // Conversion to TaylorFunction should be exact in second component
    ARIADNE_TEST_BINARY_PREDICATE(refines,f(ipt),tf(ipt));
    ARIADNE_TEST_BINARY_PREDICATE(refines,tf(ipt),f(ipt)+Float64Bounds(-1e-15,1e-15));
}

Void TestScalarTaylorFunction::test_generic() {
    ExactBoxType D={{-0.5,0.5},{-1.0,2.0}};
    Float64Value c0={0,pr};
    Float64Value c1={0,pr};
    Float64Value c2={0,pr};
    ScalarTaylorFunction pf0=ScalarTaylorFunction::constant(D,c0,swp);
    ScalarTaylorFunction pf1=ScalarTaylorFunction::constant(D,c1,swp);
    ScalarTaylorFunction pf2=ScalarTaylorFunction::constant(D,c2,swp);
    ValidatedScalarFunctionModel fm1=pf1;
    ValidatedScalarFunctionModel fm2=pf2;
    ValidatedScalarFunctionModel fm4=pf0;
    ValidatedScalarFunction f1=fm1;
    ValidatedScalarFunction f2=fm2;
    ARIADNE_TEST_EQUAL((pf1+pf2).range(),c1+c2);
//    ARIADNE_TEST_EQUAL((pf1+fm2).range(),c1+c2);
    ARIADNE_TEST_EQUAL((pf1+ f2).range(),c1+c2);
//    ARIADNE_TEST_EQUAL((fm1+pf2).range(),c1+c2);
    ARIADNE_TEST_EQUAL((fm1+fm2).range(),c1+c2);
    ARIADNE_TEST_EQUAL((fm1+ f2).range(),c1+c2);
    ARIADNE_TEST_EQUAL(( f1+pf2).range(),c1+c2);
    ARIADNE_TEST_EQUAL(( f1+fm2).range(),c1+c2);
    ARIADNE_TEST_EQUAL(( f1+ f2)(D.centre()),c1+c2);
}

/*
VectorTaylorFunction henon(const VectorTaylorFunction& x, const Vector<Float64Value>& p)
{
    VectorTaylorFunction r(2,2,x.degree()); henon(r,x,p); return r;
}
*/

class TestVectorTaylorFunction
{
    Precision64 pr;
    Sweeper<Float64> swp;
  public:
    TestVectorTaylorFunction(Sweeper<Float64> sweeper);
    Void test();
  private:
    Void test_constructors();
    Void test_restrict();
    Void test_jacobian();
    Void test_compose();
    Void test_antiderivative();
    Void test_join();
    Void test_combine();
    Void test_conversion();
    Void test_domain();
};


TestVectorTaylorFunction::TestVectorTaylorFunction(Sweeper<Float64> sweeper)
    : swp(sweeper)
{
  std::cout<<std::setprecision(17);
  std::cerr<<std::setprecision(17);
}


Void
TestVectorTaylorFunction::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_jacobian());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_join());
    ARIADNE_TEST_CALL(test_combine());
    ARIADNE_TEST_CALL(test_conversion());
    ARIADNE_TEST_CALL(test_domain());
}

Bool operator==(Expansion<Float64Value> const& e1, Expansion<RawFloat64> const& e2) {
    return reinterpret_cast<Expansion<RawFloat64>const&>(e1)==e2;
}

Void TestVectorTaylorFunction::test_constructors()
{
    Vector< Expansion<RawFloat64> > expansion(2, Expansion<RawFloat64>(2));
    expansion[0]=Expansion<RawFloat64>({ {{0,0},1.125}, {{1,0},-0.75}, {{0,1},0.0625}, {{2,0},-0.25} });
    expansion[1]=Expansion<RawFloat64>({ {{0,0},0.750}, {{1,0},0.50} });
    expansion[0].reverse_lexicographic_sort(); expansion[1].reverse_lexicographic_sort();
    Vector< RawFloat64 > errors(2);

    ExactBoxType domain={{0.25,1.25},{0.5,1.0}};
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    Real a(1.5); Real b(0.25);
    EffectiveVectorFunction henon_function={a-x[0]*x[0]+b*x[1], x[0]*1};
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,henon_model,(domain,henon_function,swp));
    ARIADNE_TEST_EQUAL(henon_model.models()[0].expansion(),expansion[0])
    ARIADNE_TEST_EQUAL(henon_model.models()[1].expansion(),expansion[1])

    Vector<Float64Bounds> e0(e(2,0),pr); Vector<Float64Bounds> e1(e(2,1),pr);

    VectorTaylorFunction t=VectorTaylorFunction::identity(domain,swp);
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,variables_model,((1.5_exact-t[0]*t[0]+0.25_exact*t[1])*e0+t[0]*e1));
    ARIADNE_TEST_EXECUTE(variables_model.simplify());
    ARIADNE_TEST_SAME(variables_model,VectorTaylorFunction(domain,expansion,errors,swp));

}

Void TestVectorTaylorFunction::test_restrict()
{
    Vector<RawFloat64> unit0={1};

    ExactBoxType domain1={{-1.0,+1.0},{-1.0,+1.0}};
    Expansion<RawFloat64> expansion1({{{0,0},1.0}, {{1,0},2.0},{{0,1},3.0}, {{2,0},4.0},{{1,1},5.0},{{0,2},6.0}, {{3,0},7.0},{{2,1},8.0},{{1,2},9.0},{{0,3},10.0}});
    Vector<ExactIntervalType> subdomain1={{-0.25,0.75},{-0.5,0.0}};
    Expansion<RawFloat64> subexpansion1({{{0,0},1.031250},{{1,0},1.812500},{{0,1},0.625000}, {{2,0},1.812500},{{1,1},0.562500},{{0,2},0.0468750},
                                        {{3,0},0.875000},{{2,1},0.500000},{{1,2},0.281250},{{0,3},0.156250}});
    Vector<RawFloat64> error1={0.0};
    VectorTaylorFunction function1(domain1,expansion1*unit0,error1,swp);
    VectorTaylorFunction restricted_function1(subdomain1,subexpansion1*unit0,error1,swp);
    ARIADNE_TEST_SAME(restriction(function1,subdomain1),restricted_function1);

    ExactBoxType domain2={{-1.0,+1.0}};
    Expansion<RawFloat64> expansion2={{{0},0.0},{{1},1.0} };
    Vector<RawFloat64> error2={0.125};
    Vector<ExactIntervalType> subdomain2={{3e-16,1.0}};
    Expansion<RawFloat64> subexpansion2={{{0},0.50000000000000022},{{1},0.49999999999999989}};
    Vector<RawFloat64> suberror2={0.125000000000000028};
    VectorTaylorFunction function2(domain2,expansion2*unit0,error2,swp);
    VectorTaylorFunction restricted_function2(subdomain2,subexpansion2*unit0,suberror2,swp);
    Float64Value::set_output_places(18);
    Float64Error::set_output_places(18);
    ARIADNE_TEST_SAME(restriction(function2,subdomain2),restricted_function2);
    ARIADNE_TEST_SAME(restriction(function2,subdomain2).domain(),restricted_function2.domain());
    ARIADNE_TEST_SAME(restriction(function2,subdomain2).error(),restricted_function2.error());
}

Void TestVectorTaylorFunction::test_jacobian()
{
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    Real a(1.5); Real b(0.25);
    EffectiveVectorFunction henon={a-x[0]*x[0]+b*x[1], x[0]*1};
    ExactBoxType domain1={{-1.0,+1.0},{-1.0,+1.0}};
    ExactBoxType domain2={{-0.5,+0.5},{-0.25,+0.25}};
    ExactBoxType domain3={{-0.25,+0.75},{0.0,+0.50}};
    //Vector<Float64Approximation> point1={0.0,0.0};
    //Vector<Float64Approximation> point2={0.5,0.25};
    Vector<Float64Bounds> point1=Vector<Float64Bounds>{0.0_exact,0.0_exact};
    Vector<Float64Bounds> point2=Vector<Float64Bounds>{0.0_exact,0.25_exact};
    //Vector<Float64Value> point1=Vector<Float64Value>{0.0,0.0};
    //Vector<Float64Value> point2=Vector<Float64Value>{0.5,0.25};
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain1,henon,swp).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain1,henon,swp).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain2,henon,swp).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain2,henon,swp).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain3,henon,swp).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain3,henon,swp).jacobian(point2),henon.jacobian(point2));
}

Void TestVectorTaylorFunction::test_compose()
{
    Real a(1.5); Real b(0.25);
    EffectiveScalarFunction x=EffectiveScalarFunction::coordinate(2,0);
    EffectiveScalarFunction y=EffectiveScalarFunction::coordinate(2,1);
    EffectiveVectorFunction henon_polynomial=(a-x*x+b*y)*e(2,0)+x*e(2,1);
    EffectiveVectorFunction henon_square_polynomial=
        {a*(1-a)+b*x-2*a*b*y+2*a*x*x-b*b*y*y+2*b*x*x*y-x*x*x*x, a-x*x+b*y};
    //    compose(henon_polynomial,henon_polynomial);
    ExactBoxType domain1={{0.25,1.25},{0.5,1.0}};
    ExactBoxType domain2={{-1.5,2.5},{0.25,1.25}};
    ARIADNE_TEST_PRINT((a-x*x+b*y));
    ARIADNE_TEST_PRINT(e(2,0));
    ARIADNE_TEST_PRINT((a-x*x+b*y)*e(2,0));
    ARIADNE_TEST_PRINT(henon_polynomial);
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,function1,(domain1,henon_polynomial,swp));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,function2,(domain2,henon_polynomial,swp));

    VectorTaylorFunction composition1(domain1,henon_square_polynomial,swp);
    ARIADNE_TEST_SAME(compose(function2,function1),composition1);
}


Void TestVectorTaylorFunction::test_antiderivative()
{
    SizeType index0=0;
    SizeType index1=1;

    Vector<RawFloat64> unit0={1};
    ExactBoxType domain1={{-1,+1},{-1,+1}};
    Expansion<RawFloat64> expansion1={{{0,0},3.0}};
    VectorTaylorFunction function1(domain1,expansion1*unit0,swp);
    Expansion<RawFloat64> aexpansion1={{{0,1},3.0}};
    VectorTaylorFunction antiderivative1(domain1,aexpansion1*unit0,swp);
    ARIADNE_TEST_SAME(antiderivative(function1,index1),antiderivative1);

    ExactBoxType domain2={{-0.25,0.75},{0.0,0.5}};
    Expansion<RawFloat64> expansion2={{{0,0},3.0}};
    VectorTaylorFunction function2(domain2,expansion2*unit0,swp);
    Expansion<RawFloat64> aexpansion2={{{0,1},0.75}};
    VectorTaylorFunction antiderivative2(domain2,aexpansion2*unit0,swp);
    ARIADNE_TEST_SAME(antiderivative(function2,index1),antiderivative2);

    ExactBoxType domain3={{-0.25,0.75},{0.0,0.5}};
    Expansion<RawFloat64> expansion3={{{0,0},1.0},{{1,0},2.0},{{0,1},3.0},{{2,0},4.0},{{1,1},5.0},{{0,2},6.0}};
    VectorTaylorFunction function3(domain3,expansion3*unit0,swp);
    Expansion<RawFloat64> aexpansion30={{{1,0},0.5},{{2,0},0.5},{{1,1},1.5},{{3,0},0.66666666666666663},{{2,1},1.25},{{1,2},3.0}};
    Vector<Float64> aerror30={5.5511151231257827e-17};
    VectorTaylorFunction antiderivative30(domain3,aexpansion30*unit0,aerror30,swp);
    ARIADNE_TEST_SAME(antiderivative(function3,index0),antiderivative30);
    Expansion<RawFloat64> aexpansion31={{{0,1},0.25},{{1,1},0.5},{{0,2},0.375},{{2,1},1.0},{{1,2},0.625},{{0,3},0.5}};
    VectorTaylorFunction antiderivative31(domain3,aexpansion31*unit0,swp);
    ARIADNE_TEST_SAME(antiderivative(function3,index1),antiderivative31);

}

Void TestVectorTaylorFunction::test_join()
{
    ExactBoxType domain={{-0.25,+0.25},{-0.5,+0.5}};
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    EffectiveVectorFunction function1 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(1,0);
    EffectiveVectorFunction function2 = (4*x[0]*x[0]+5*x[0]*x[1]+6*x[1]*x[1])*e(2,1);
    EffectiveVectorFunction function3 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(3,0)
        + (4*x[0]*x[0]+5*x[0]*x[1]+6*x[1]*x[1])*e(3,2);

    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,taylorfunction1,(domain,function1,swp));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,taylorfunction2,(domain,function2,swp));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,taylorfunction3,(domain,function3,swp));
    ARIADNE_TEST_SAME(join(taylorfunction1,taylorfunction2),taylorfunction3);

}

Void TestVectorTaylorFunction::test_combine()
{
    // This test contains a regression test to check correct behaviour for a zero component.
    ExactBoxType domain1={{-0.25,+0.25},{-0.5,+0.5}};
    ExactBoxType domain2={{-0.75,+0.75},{-1.0,+1.0},{-1.25,+1.25}};
    ExactBoxType domain3={{-0.25,+0.25},{-0.5,+0.5},{-0.75,+0.75},{-1.0,+1.0},{-1.25,+1.25}};
    EffectiveVectorFunction x;
    x=EffectiveVectorFunction::identity(2);
    EffectiveVectorFunction function1 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(1,0);
    x=EffectiveVectorFunction::identity(3);
    EffectiveVectorFunction function2 = (4*x[0]*x[0]+5*x[0]*x[1]+6*x[1]*x[2])*e(2,1);
    x=EffectiveVectorFunction::identity(5);
    EffectiveVectorFunction function3 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(3,0)
        + (4*x[2]*x[2]+5*x[2]*x[3]+6*x[3]*x[4])*e(3,2);
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,taylorfunction1,(domain1,function1,swp));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,taylorfunction2,(domain2,function2,swp));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,taylorfunction3,(domain3,function3,swp));
    ARIADNE_TEST_SAME(combine(taylorfunction1,taylorfunction2),taylorfunction3);

}

Void TestVectorTaylorFunction::test_conversion()
{
    // Test conversion between ordinary functions and Taylor functions.
    ExactBoxType D={{-0.5,0.5},{-1.0,2.0}};
    Vector<RawFloat64> pt={-0.25,0.25};
    Vector<Float64Bounds> ipt(pt);
    Vector<Float64Approximation> apt(pt);
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);

    EffectiveVectorFunction h={1-x[0]*x[0]-x[1]/2,x[0]+Real(0)};
    VectorTaylorFunction th(D,h,swp);

    ARIADNE_TEST_PRINT(h);
    ARIADNE_TEST_PRINT(th);

    // Conversion to TaylorFunction should be exact in second component
    ARIADNE_TEST_SAME(th(apt)[1],h(apt)[1]);
    ARIADNE_TEST_SAME(th(ipt)[1],h(ipt)[1]);
    ARIADNE_TEST_BINARY_PREDICATE(refines,h[0](ipt),th[0](ipt));


}

// Regression test for domain with empty interior
Void TestVectorTaylorFunction::test_domain()
{
    EffectiveScalarFunction z=EffectiveScalarFunction::constant(2,0);
    EffectiveScalarFunction o=EffectiveScalarFunction::constant(2,1);
    EffectiveScalarFunction x0=EffectiveScalarFunction::coordinate(2,0);
    EffectiveScalarFunction x1=EffectiveScalarFunction::coordinate(2,1);

    ExactBoxType D1={{-1.0,1.0},{-1.0,1.0}};
    VectorTaylorFunction t1(D1, {o,x0+x1}, swp);
    ARIADNE_TEST_PRINT(t1);
    ARIADNE_TEST_PRINT(t1.codomain());
    ExactBoxType D2={{1.0,1.0},{-2.0,2.0}};
    ScalarTaylorFunction t2(D2,2*x0+x1*x1,swp);
    ARIADNE_TEST_PRINT(t2.domain());
    ARIADNE_TEST_PRINT(t2.model());
    ARIADNE_TEST_PRINT(t2.codomain());

    ARIADNE_TEST_PRINT(t2);
    ARIADNE_TEST_PRINT(compose(t2,t1));
    ScalarTaylorFunction t3(D1,2+(x0+x1)*(x0+x1),swp);
    ARIADNE_TEST_SAME(compose(t2,t1),t3);

    Vector<Float64Bounds> x={Float64Bounds{1.0,1.0},Float64Bounds{0.5,1.5}};
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_SAME(evaluate(t2,x),Float64Bounds(2.25,4.25));

    // Ensure evaluation and composition throw errors when expected
    Vector<Float64Bounds> xe={Float64Bounds{0.875,1.125},Float64Bounds{0.5,1.5}};
    ARIADNE_TEST_THROWS(t2(xe),DomainException);
    ARIADNE_TEST_THROWS(evaluate(t2,xe),DomainException);

    // Ensure evaluation and composition throw errors when expected
    VectorTaylorFunction vt2={t2};
    ARIADNE_TEST_THROWS(vt2(xe),DomainException);
    ARIADNE_TEST_THROWS(evaluate(vt2,xe),DomainException);

    VectorTaylorFunction te1=t1; te1[0]=te1[0]+Float64Bounds(-0.125,+0.125);
    ARIADNE_TEST_THROWS(compose(t2,te1),DomainException);
    ARIADNE_TEST_THROWS(compose(vt2,te1),DomainException);

    ARIADNE_TEST_SAME(unchecked_evaluate(t2,xe),Float64Bounds(2.25,4.25));

    // Regression test for printing functions with trivial domain component
    ExactBoxType D4={{1.0,1.0},{0.0,2.0}};
    ScalarTaylorFunction st40(D4, x0, swp);
    ScalarTaylorFunction st41(D4, x1, swp);
    ARIADNE_TEST_PRINT(st40);
    ARIADNE_TEST_PRINT(st41);
    VectorTaylorFunction t4(D4, {x0,x1}, swp);
    ARIADNE_TEST_PRINT(t4);
}


class TestTaylorFunctionFactory
{
  public:
    TestTaylorFunctionFactory();
    Void test();
  private:
    Void test_create();
};

TestTaylorFunctionFactory::TestTaylorFunctionFactory()
{
}

Void TestTaylorFunctionFactory::test()
{
    ARIADNE_TEST_CALL(test_create());
}

Void TestTaylorFunctionFactory::test_create()
{
    Sweeper<Float64> sweeper(new ThresholdSweeper<Float64>(Precision64(),1e-4));
    TaylorFunctionFactory factory(sweeper);

    Vector<ExactIntervalType> dom={{-1,+1},{0.5,3.5}};
    Vector<Float64Bounds> args=reinterpret_cast<Vector<Float64Bounds>const&>(dom);

    ScalarTaylorFunction stf=factory.create(dom, EffectiveScalarFunction::zero(dom.size()) );
    ARIADNE_TEST_PRINT(stf);
    ARIADNE_TEST_EQUALS(&stf.properties(),&sweeper);
    ARIADNE_TEST_EQUALS(stf(args),Float64Bounds(0.0));
    ARIADNE_TEST_EQUALS(evaluate(stf,args),Float64Bounds(0.0));

    VectorTaylorFunction vtf=factory.create(dom, EffectiveVectorFunction::identity(dom.size()) );
    ARIADNE_TEST_PRINT(vtf);

    // Test evaluation gives a superset with small additional error
    Vector<Float64Bounds> errs(2,Float64Bounds(-1e-15,+1e-15));
    ARIADNE_TEST_BINARY_PREDICATE(refines,args,vtf(args));
    ARIADNE_TEST_BINARY_PREDICATE(refines,vtf(args),Vector<Float64Bounds>(args+errs));
    Vector<Float64Bounds> pt(2); pt[0]=Float64Bounds(0.2); pt[1]=Float64Bounds(1.25);
    ARIADNE_TEST_BINARY_PREDICATE(refines,pt,vtf(pt));
}



Int main() {
    ThresholdSweeper<Float64> sweeper(Precision64(),std::numeric_limits<float>::epsilon());
    TestScalarTaylorFunction(sweeper).test();
    TestVectorTaylorFunction(sweeper).test();
    TestTaylorFunctionFactory().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}

