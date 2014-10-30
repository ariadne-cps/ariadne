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
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "algebra/expansion.h"
#include "function/taylor_model.h"
#include "function/taylor_function.h"
#include "algebra/differential.h"
#include "function/polynomial.h"
#include "function/function.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

Vector<ExactFloat> e(uint n, uint i) { return Vector<ExactFloat>::unit(n,i); }
Expansion<RawFloat> v(uint n, uint j) { return Expansion<RawFloat>::variable(n,j); }
Polynomial<Float> p(uint n, uint j) { return Polynomial<Float>::variable(n,j); }
ScalarTaylorFunction t(ExactBox d, uint j,Sweeper swp) { return ScalarTaylorFunction::variable(d,j,swp); }

namespace Ariadne {
std::pair<Float, Vector<ExactInterval> > flow_bounds(EffectiveVectorFunction const&,Vector<ExactInterval> const&,Float const&);
}



class TestScalarTaylorFunction
{
    Sweeper swp;
  public:
    TestScalarTaylorFunction(Sweeper sweeper);
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_predicates();
    void test_approximation();
    void test_evaluate();
    void test_arithmetic();
    void test_functions();
    void test_compose();
    void test_antiderivative();
    void test_substitute();
    void test_conversion();
  private:
    ExactBox d(unsigned int n) { return Vector<ExactInterval>(n,ExactInterval(-1,+1)); }
    typedef Expansion<RawFloat> e;
    typedef TaylorModel<ValidatedTag> TM;
};


TestScalarTaylorFunction::TestScalarTaylorFunction(Sweeper sweeper)
    : swp(sweeper)
{
}


void TestScalarTaylorFunction::test()
{
    std::clog<<std::setprecision(17);
    std::cerr<<std::setprecision(17);
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_predicates());
    ARIADNE_TEST_CALL(test_approximation());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_functions());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_substitute());
    ARIADNE_TEST_CALL(test_conversion());
}


void TestScalarTaylorFunction::test_concept()
{
    const ExactFloat f=0;
    const ValidatedFloat i;
    const Vector<ExactFloat> vf;
    const Vector<ValidatedFloat> vi;
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

    tr=exp(t); tr=log(t); tr=sqrt(t);
    tr=sin(t); tr=cos(t); tr=tan(t);
    //tr=asin(t); tr=acos(t); tr=atan(t);

    tr.sweep(); tr.clobber();

    t.evaluate(vi); evaluate(t,vi);
    t.domain(); t.range(); t.expansion(); t.error();

}

void TestScalarTaylorFunction::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv1,({{-1,+1},{-1,+1}},{{{0,0},1.},{{1,0},2.},{{0,1},3.},{{2,0},4.},{{1,1},5.},{{0,2},6.},{{3,0},7.},{{2,1},8.},{{1,2},9.},{{0,3},10.}},0.25,swp));

    ARIADNE_ASSERT_EQUAL(tv1.domain(),Vector<ExactInterval>({{-1,+1},{-1,+1}}));
    ARIADNE_ASSERT_EQUAL(tv1.argument_size(),2);
    ARIADNE_ASSERT_EQUAL(tv1.number_of_nonzeros(),10);
    ARIADNE_ASSERT_EQUAL(tv1.value(),1.0);
    ARIADNE_ASSERT_EQUAL(tv1.error(),0.25);
}

void TestScalarTaylorFunction::test_predicates()
{
    ScalarTaylorFunction tv1(d(1),e(1,2, {1.00,2.00,3.00}), 0.75, swp);
    ScalarTaylorFunction tv2(d(1),e(1,2, {1.00,1.75,3.25}), 0.25, swp);
    ScalarTaylorFunction tv3(d(1),e(1,2, {1.125,1.75,3.25}), 0.25, swp);
    ScalarTaylorFunction tv4(d(1),e(1,3, {1.00,2.25,3.00,-0.25}), 0.25, swp);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

void TestScalarTaylorFunction::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv1,(d(2),e(2,3,{1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0}),0.25,swp));
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv2,(d(1),e(1,2,{1.0,2.0,3.0}),0.25,swp));
}

void TestScalarTaylorFunction::test_evaluate()
{
    Vector<ValidatedFloat> iv({{0.25,0.5},{-0.75,-0.5}});
    ScalarTaylorFunction tv(d(2),e(2,2,{1.0,2.0,3.0,4.0,5.0,6.0}),0.25,swp);
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),ValidatedFloat(-1,1));
}

void TestScalarTaylorFunction::test_arithmetic()
{
    ARIADNE_TEST_EQUAL(d(1),d(1));
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)+(-3), ScalarTaylorFunction(d(1),e(1,2, {-2.0,-2.0,3.0}), 0.75, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)-(-3), ScalarTaylorFunction(d(1),e(1,2, {4.0,-2.0,3.0}), 0.75, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)*(-3), ScalarTaylorFunction(d(1),e(1,2, {-3.0,6.0,-9.0}), 2.25, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)/(-4), ScalarTaylorFunction(d(1),e(1,2, {-0.25,0.5,-0.75}), 0.1875, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)+ValidatedFloat(-1,2), ScalarTaylorFunction(d(1),e(1,2, {1.5,-2.0,3.0}), 2.25, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)-ValidatedFloat(-1,2), ScalarTaylorFunction(d(1),e(1,2, {0.5,-2.0,3.0}), 2.25, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)*ValidatedFloat(-1,2), ScalarTaylorFunction(d(1),e(1,2, {0.5,-1.0,1.5}), 10.5, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)/ValidatedFloat(0.25,2.0), ScalarTaylorFunction(d(1),e(1,2, {2.25,-4.5,6.75}), 13.5, swp));
    ARIADNE_TEST_EQUAL(+ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp), ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp));
    ARIADNE_TEST_EQUAL(-ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp), ScalarTaylorFunction(d(1),e(1,2, {-1.0,2.0,-3.0}), 0.75, swp));

    // Regression test to check subtraction yielding zero coefficients
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)+ScalarTaylorFunction(d(1),e(1,2, {3.0,2.0,-4.0}), 0.5, swp), ScalarTaylorFunction(d(1),e(1,2, {4.0,0.0,-1.0}), 1.25, swp));

    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)-ScalarTaylorFunction(d(1),e(1,2, {3.0,2.0,-4.0}), 0.5, swp), ScalarTaylorFunction(d(1),e(1,2, {-2.0,-4.0,7.0}), 1.25, swp));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, {1.0,-2.0,3.0}), 0.75, swp)*ScalarTaylorFunction(d(1),e(1,2, {3.0,2.0,-4.0}), 0.5, swp), ScalarTaylorFunction(d(1),e(1,4, {3.0,-4.0,1.0,14.0,-12.0}), 10.125, swp));

}

void TestScalarTaylorFunction::test_functions()
{
    ScalarTaylorFunction xz(d(1),e(1,1, {0.0, 0.5}), 0.0, swp);
    ScalarTaylorFunction xo(d(1),e(1,1, {1.0, 0.5}), 0.0, swp);

    //Functions based on their natural defining points
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),ScalarTaylorFunction(d(1),e(1,6, {1.00000,0.50000,0.12500,0.02083,0.00260,0.00026,0.00002}), 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),ScalarTaylorFunction(d(1),e(1,6, {0.00000,0.50000,0.0000,-0.02083,0.00000,0.00026,0.00000}), 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),ScalarTaylorFunction(d(1),e(1,6, {1.00000,0.0000,-0.12500,0.00000,0.00260,0.0000,-0.00002}), 0.00003, swp));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),ScalarTaylorFunction(d(1),e(1,6,  {1.000000,-0.500000, 0.250000,-0.125000, 0.062500,-0.031250, 0.015625}), 0.018, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),ScalarTaylorFunction(d(1),e(1,6, {1.000000, 0.250000,-0.031250, 0.007813,-0.002441, 0.000854,-0.000320}), 0.0003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),ScalarTaylorFunction(d(1),e(1,6,  {0.000000, 0.500000,-0.125000, 0.041667,-0.015625, 0.006250,-0.002604}), 0.003, swp));

}


void TestScalarTaylorFunction::test_compose()
{
}


void TestScalarTaylorFunction::test_antiderivative()
{
    ScalarTaylorFunction tm=ScalarTaylorFunction::constant(d(2),1.0,swp);
    ScalarTaylorFunction atm=antiderivative(tm,1u);

    ExactIntervalVector dom({{1.0, 4.0}, {0.0, 2.0}});
    VectorTaylorFunction x = VectorTaylorFunction::identity(dom,swp);
    ScalarTaylorFunction f = 1.0 + 2.0*x[0]+3.0*x[1]+5.0*x[0]*x[0]+4.0*x[0]*x[1];

    ARIADNE_TEST_PRINT(f);
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction, g, (antiderivative(f,0,2.0_exact)) );
    ARIADNE_ASSERT(norm(derivative(g,0)-f)<1e-8);

    // We should have f(c,s)=0 for all x1
    ScalarTaylorFunction s = ScalarTaylorFunction::identity(dom[1],swp);
    ScalarTaylorFunction c = ScalarTaylorFunction::constant(s.domain(),2.0_exact,swp);
    ScalarTaylorFunction h=compose(g,VectorTaylorFunction((c,s)));

    Vector<ExactInterval> hdom=h.domain();
    Vector<ValidatedFloat> domv(reinterpret_cast<Vector<ValidatedFloat>const&>(hdom));
    ARIADNE_ASSERT(mag(h(domv))<1e-8);

}

void TestScalarTaylorFunction::test_substitute()
{
    ExactBox d1={{-0.75,+0.75}};
    ScalarTaylorFunction tu=ScalarTaylorFunction::coordinate(d1,0,swp);

    ExactBox d2={{-0.75,+0.75},{-0.25,+0.75}};
    ScalarTaylorFunction tx=ScalarTaylorFunction::coordinate(d2,0,swp);
    ScalarTaylorFunction ty=ScalarTaylorFunction::coordinate(d2,1,swp);

    ScalarTaylorFunction tf=1+3*tx*ty;
    ScalarTaylorFunction ts=0.25+tu*tu;

    ScalarTaylorFunction tg=substitute(tf,1,ts);
    ScalarTaylorFunction tr=1+0.75*tu+3*tu*tu*tu;
    ARIADNE_TEST_EQUAL(tg,tr);
}

void TestScalarTaylorFunction::test_conversion() {
    // Test conversion between ordinary functions and Taylor functions.
    ExactBox D={{-0.5,0.5},{-1.0,2.0}};
    Vector<ExactFloat> pt={-0.25_exact,0.25_exact};
    Vector<ValidatedFloat> ipt(pt);
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);

    EffectiveScalarFunction f=(1-x[0]*x[0]-x[1]/2);
    ScalarTaylorFunction tf(D,f,swp);

    ARIADNE_TEST_PRINT(f);
    ARIADNE_TEST_PRINT(tf);

    // Conversion to TaylorFunction should be exact in second component
    ARIADNE_TEST_BINARY_PREDICATE(refines,f(ipt),tf(ipt));
    ARIADNE_TEST_BINARY_PREDICATE(refines,tf(ipt),f(ipt)+ValidatedFloat(-1e-15,1e-15));
}


/*
VectorTaylorFunction henon(const VectorTaylorFunction& x, const Vector<ExactFloat>& p)
{
    VectorTaylorFunction r(2,2,x.degree()); henon(r,x,p); return r;
}
*/

class TestVectorTaylorFunction
{
    Sweeper swp;
  public:
    TestVectorTaylorFunction(Sweeper sweeper);
    void test();
  private:
    void test_constructors();
    void test_restrict();
    void test_jacobian();
    void test_compose();
    void test_antiderivative();
    void test_join();
    void test_combine();
    void test_conversion();
    void test_domain();
};


TestVectorTaylorFunction::TestVectorTaylorFunction(Sweeper sweeper)
    : swp(sweeper)
{
  std::cout<<std::setprecision(17);
  std::cerr<<std::setprecision(17);
}


void
TestVectorTaylorFunction::test()
{
    ARIADNE_TEST_CALL(test_combine());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_jacobian());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_join());
    ARIADNE_TEST_CALL(test_conversion());
    ARIADNE_TEST_CALL(test_domain());
}

bool operator==(Expansion<ExactFloat> const& e1, Expansion<RawFloat> const& e2) {
    return reinterpret_cast<Expansion<RawFloat>const&>(e1)==e2;
}

void TestVectorTaylorFunction::test_constructors()
{
    Vector< Expansion<RawFloat> > expansion(2);
    expansion[0]=Expansion<RawFloat>({ {{0,0},1.125}, {{1,0},-0.75}, {{0,1},0.0625}, {{2,0},-0.25} });
    expansion[1]=Expansion<RawFloat>({ {{0,0},0.750}, {{1,0},0.50} });
    expansion[0].reverse_lexicographic_sort(); expansion[1].reverse_lexicographic_sort();
    Vector< RawFloat > errors(2);

    ExactBox domain={{0.25,1.25},{0.5,1.0}};
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    ExactFloat a(1.5); ExactFloat b(0.25);
    EffectiveVectorFunction henon_function( (a-x[0]*x[0]+b*x[1], x[0]*1) );
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,henon_model,(domain,henon_function,swp));
    ARIADNE_TEST_EQUAL(henon_model.models()[0].expansion(),expansion[0])
    ARIADNE_TEST_EQUAL(henon_model.models()[1].expansion(),expansion[1])

    Vector<ValidatedFloat> e0=e(2,0); Vector<ValidatedFloat> e1=e(2,1);

    VectorTaylorFunction t=VectorTaylorFunction::identity(domain,swp);
    //VectorTaylorFunction variables_model((1.5-t[0]*t[0]+0.25*t[1])*e0+t[0]*e1);
    VectorTaylorFunction variables_model(ScalarTaylorFunction(1.5-t[0]*t[0]+0.25*t[1])*e0+ScalarTaylorFunction(t[0])*e1);
    variables_model.sweep();
    ARIADNE_TEST_EQUAL(variables_model,VectorTaylorFunction(domain,expansion,errors,swp));

}

void TestVectorTaylorFunction::test_restrict()
{
    Vector<RawFloat> unit0={1};

    ExactBox domain1={{-1.0,+1.0},{-1.0,+1.0}};
    Expansion<RawFloat> expansion1(2,3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});
    Vector<ExactInterval> subdomain1={{-0.25,0.75},{-0.5,0.0}};
    Expansion<RawFloat> subexpansion1(2,3, {1.031250, 1.812500,0.625000, 1.812500,0.562500,0.0468750,
                                             0.875000,0.500000,0.281250,0.156250});
    Vector<RawFloat> error1={0.0};
    VectorTaylorFunction function1(domain1,expansion1*unit0,error1,swp);
    VectorTaylorFunction restricted_function1(subdomain1,subexpansion1*unit0,error1,swp);
    ARIADNE_TEST_EQUAL(restrict(function1,subdomain1),restricted_function1);

    ExactBox domain2={{-1.0,+1.0}};
    Expansion<RawFloat> expansion2={{{0},0.0},{{1},1.0} };
    Vector<RawFloat> error2={0.125};
    Vector<ExactInterval> subdomain2={{3e-16,1.0}};
    Expansion<RawFloat> subexpansion2={{{0},0.50000000000000022},{{1},0.49999999999999989}};
    Vector<RawFloat> suberror2={0.12500000000000008};
    VectorTaylorFunction function2(domain2,expansion2*unit0,error2,swp);
    VectorTaylorFunction restricted_function2(subdomain2,subexpansion2*unit0,suberror2,swp);
    ARIADNE_TEST_EQUAL(restrict(function2,subdomain2),restricted_function2);
}

void TestVectorTaylorFunction::test_jacobian()
{
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    ExactFloat a(1.5); ExactFloat b(0.25);
    EffectiveVectorFunction henon( (a-x[0]*x[0]+b*x[1], x[0]*1) );
    ExactBox domain1={{-1.0,+1.0},{-1.0,+1.0}};
    ExactBox domain2={{-0.5,+0.5},{-0.25,+0.25}};
    ExactBox domain3={{-0.25,+0.75},{0.0,+0.50}};
    //Vector<ApproximateFloat> point1={0.0,0.0};
    //Vector<ApproximateFloat> point2={0.5,0.25};
    Vector<ValidatedFloat> point1=Vector<ValidatedFloat>{0,0};
    Vector<ValidatedFloat> point2=Vector<ValidatedFloat>{0,0.25_exact};
    //Vector<ExactFloat> point1=Vector<ExactFloat>{0.0,0.0};
    //Vector<ExactFloat> point2=Vector<ExactFloat>{0.5,0.25};
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain1,henon,swp).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain1,henon,swp).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain2,henon,swp).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain2,henon,swp).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain3,henon,swp).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain3,henon,swp).jacobian(point2),henon.jacobian(point2));
}

void TestVectorTaylorFunction::test_compose()
{
    Real a=ExactFloat(1.5); Real b=ExactFloat(0.25);
    EffectiveScalarFunction x=EffectiveScalarFunction::coordinate(2,0);
    EffectiveScalarFunction y=EffectiveScalarFunction::coordinate(2,1);
    EffectiveVectorFunction henon_polynomial=(a-x*x+b*y)*e(2,0)+x*e(2,1);
    EffectiveVectorFunction henon_square_polynomial=
        (a*(1-a)+b*x-2*a*b*y+2*a*x*x-b*b*y*y+2*b*x*x*y-x*x*x*x)*e(2,0)
            + (a-x*x+b*y)*e(2,1);
    //    compose(henon_polynomial,henon_polynomial);
    ExactBox domain1={{0.25,1.25},{0.5,1.0}};
    ExactBox domain2={{-1.5,2.5},{0.25,1.25}};
    ARIADNE_TEST_PRINT((a-x*x+b*y));
    ARIADNE_TEST_PRINT(e(2,0));
    ARIADNE_TEST_PRINT((a-x*x+b*y)*e(2,0));
    ARIADNE_TEST_PRINT(henon_polynomial);
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,function1,(domain1,henon_polynomial,swp));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,function2,(domain2,henon_polynomial,swp));

    VectorTaylorFunction composition1(domain1,henon_square_polynomial,swp);
    ARIADNE_TEST_EQUAL(compose(function2,function1),composition1);
}


void TestVectorTaylorFunction::test_antiderivative()
{
    unsigned int index0=0;
    unsigned int index1=1;

    Vector<RawFloat> unit0={1};
    ExactBox domain1={{-1,+1},{-1,+1}};
    Expansion<RawFloat> expansion1={{{0,0},3.0}};
    VectorTaylorFunction function1(domain1,expansion1*unit0,swp);
    Expansion<RawFloat> aexpansion1={{{0,1},3.0}};
    VectorTaylorFunction antiderivative1(domain1,aexpansion1*unit0,swp);
    ARIADNE_TEST_EQUAL(antiderivative(function1,index1),antiderivative1);

    ExactBox domain2={{-0.25,0.75},{0.0,0.5}};
    Expansion<RawFloat> expansion2={{{0,0},3.0}};
    VectorTaylorFunction function2(domain2,expansion2*unit0,swp);
    Expansion<RawFloat> aexpansion2={{{0,1},0.75}};
    VectorTaylorFunction antiderivative2(domain2,aexpansion2*unit0,swp);
    ARIADNE_TEST_EQUAL(antiderivative(function2,index1),antiderivative2);

    ExactBox domain3={{-0.25,0.75},{0.0,0.5}};
    Expansion<RawFloat> expansion3={{{0,0},1.0},{{1,0},2.0},{{0,1},3.0},{{2,0},4.0},{{1,1},5.0},{{0,2},6.0}};
    VectorTaylorFunction function3(domain3,expansion3*unit0,swp);
    Expansion<RawFloat> aexpansion30={{{1,0},0.5},{{2,0},0.5},{{1,1},1.5},{{3,0},0.66666666666666663},{{2,1},1.25},{{1,2},3.0}};
    Vector<Float> aerror30={5.5511151231257827e-17};
    VectorTaylorFunction antiderivative30(domain3,aexpansion30*unit0,aerror30,swp);
    ARIADNE_TEST_EQUAL(antiderivative(function3,index0),antiderivative30);
    Expansion<RawFloat> aexpansion31={{{0,1},0.25},{{1,1},0.5},{{0,2},0.375},{{2,1},1.0},{{1,2},0.625},{{0,3},0.5}};
    VectorTaylorFunction antiderivative31(domain3,aexpansion31*unit0,swp);
    ARIADNE_TEST_EQUAL(antiderivative(function3,index1),antiderivative31);

}

void TestVectorTaylorFunction::test_join()
{
    ExactBox domain={{-0.25,+0.25},{-0.5,+0.5}};
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);
    EffectiveVectorFunction function1 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(1,0);
    EffectiveVectorFunction function2 = (4*x[0]*x[0]+5*x[0]*x[1]+6*x[1]*x[1])*e(2,1);
    EffectiveVectorFunction function3 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(3,0)
        + (4*x[0]*x[0]+5*x[0]*x[1]+6*x[1]*x[1])*e(3,2);
    VectorTaylorFunction taylorfunction1(domain,function1,swp);
    VectorTaylorFunction taylorfunction2(domain,function2,swp);
    VectorTaylorFunction taylorfunction3(domain,function3,swp);
    ARIADNE_TEST_EQUAL(join(taylorfunction1,taylorfunction2),taylorfunction3);

}

void TestVectorTaylorFunction::test_combine()
{
    // This test contains a regression test to check correct behaviour for a zero component.
    ExactBox domain1={{-0.25,+0.25},{-0.5,+0.5}};
    ExactBox domain2={{-0.75,+0.75},{-1.0,+1.0},{-1.25,+1.25}};
    ExactBox domain3={{-0.25,+0.25},{-0.5,+0.5},{-0.75,+0.75},{-1.0,+1.0},{-1.25,+1.25}};
    EffectiveVectorFunction x;
    x=EffectiveVectorFunction::identity(2);
    EffectiveVectorFunction function1 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(1,0);
    x=EffectiveVectorFunction::identity(3);
    EffectiveVectorFunction function2 = (4*x[0]*x[0]+5*x[0]*x[1]+6*x[1]*x[2])*e(2,1);
    x=EffectiveVectorFunction::identity(5);
    EffectiveVectorFunction function3 = (x[0]*x[0]+2*x[0]*x[1]+3*x[1]*x[1])*e(3,0)
        + (4*x[2]*x[2]+5*x[2]*x[3]+6*x[3]*x[4])*e(3,2);
    VectorTaylorFunction taylorfunction1(domain1,function1,swp);
    VectorTaylorFunction taylorfunction2(domain2,function2,swp);
    VectorTaylorFunction taylorfunction3(domain3,function3,swp);
    ARIADNE_TEST_EQUAL(combine(taylorfunction1,taylorfunction2),taylorfunction3);

}

void TestVectorTaylorFunction::test_conversion()
{
    // Test conversion between ordinary functions and Taylor functions.
    ExactBox D={{-0.5,0.5},{-1.0,2.0}};
    Vector<RawFloat> pt={-0.25,0.25};
    Vector<ValidatedFloat> ipt(pt);
    Vector<ApproximateFloat> apt(pt);
    EffectiveVectorFunction x=EffectiveVectorFunction::identity(2);

    EffectiveVectorFunction h=EffectiveVectorFunction((1-x[0]*x[0]-x[1]/2,x[0]+Real(0)));
    VectorTaylorFunction th(D,h,swp);

    ARIADNE_TEST_PRINT(h);
    ARIADNE_TEST_PRINT(th);

    // Conversion to TaylorFunction should be exact in second component
    ARIADNE_TEST_EQUAL(th(apt)[1],h(apt)[1]);
    ARIADNE_TEST_EQUAL(th(ipt)[1],h(ipt)[1]);
    ARIADNE_TEST_BINARY_PREDICATE(refines,h[0](ipt),th[0](ipt));


}

// Regression test for domain with empty interior
void TestVectorTaylorFunction::test_domain()
{
    EffectiveScalarFunction z=EffectiveScalarFunction::constant(2,0);
    EffectiveScalarFunction o=EffectiveScalarFunction::constant(2,1);
    EffectiveScalarFunction x0=EffectiveScalarFunction::coordinate(2,0);
    EffectiveScalarFunction x1=EffectiveScalarFunction::coordinate(2,1);

    ExactBox D1={{-1.0,1.0},{-1.0,1.0}};
    VectorTaylorFunction t1(D1, (o,x0+x1), swp);
    ARIADNE_TEST_PRINT(t1);
    ARIADNE_TEST_PRINT(t1.codomain());
    ExactBox D2={{1.0,1.0},{-2.0,2.0}};
    ScalarTaylorFunction t2(D2,2*x0+x1*x1,swp);
    ARIADNE_TEST_PRINT(t2.domain());
    ARIADNE_TEST_PRINT(t2.model());
    ARIADNE_TEST_PRINT(t2.codomain());

    ARIADNE_TEST_PRINT(t2);
    ARIADNE_TEST_PRINT(compose(t2,t1));
    ScalarTaylorFunction t3(D1,2+(x0+x1)*(x0+x1),swp);
    ARIADNE_TEST_EQUAL(compose(t2,t1),t3);

    Vector<ValidatedFloat> x={{1.0,1.0},{0.5,1.5}};
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_EQUAL(evaluate(t2,x),ValidatedFloat(2.25,4.25));

    // Ensure evaluation and composition throw errors when expected
    Vector<ValidatedFloat> xe={{0.875,1.125},{0.5,1.5}};
    ARIADNE_TEST_THROWS(evaluate(t2,xe),DomainException);

    VectorTaylorFunction te1=t1; te1[0]=te1[0]+ValidatedFloat(-0.125,+0.125);
    ARIADNE_TEST_THROWS(compose(t2,te1),DomainException);

    ARIADNE_TEST_EQUAL(unchecked_evaluate(t2,xe),ValidatedFloat(2.25,4.25));

    // Regression test for printing functions with trivial domain component
    ExactBox D4={{1.0,1.0},{0.0,2.0}};
    ScalarTaylorFunction st40(D4, x0, swp);
    ScalarTaylorFunction st41(D4, x1, swp);
    ARIADNE_TEST_PRINT(st40);
    ARIADNE_TEST_PRINT(st41);
    VectorTaylorFunction t4(D4, (x0,x1), swp);
    ARIADNE_TEST_PRINT(t4);
}


class TestTaylorFunctionFactory
{
  public:
    TestTaylorFunctionFactory();
    void test();
  private:
    void test_create();
};

TestTaylorFunctionFactory::TestTaylorFunctionFactory()
{
}

void TestTaylorFunctionFactory::test()
{
    ARIADNE_TEST_CALL(test_create());
}

void TestTaylorFunctionFactory::test_create()
{
    Sweeper sweeper(new ThresholdSweeper(1e-4));
    TaylorFunctionFactory factory(sweeper);

    Vector<ExactInterval> dom={{-1,+1},{0.5,3.5}};
    Vector<ValidatedFloat> args=reinterpret_cast<Vector<ValidatedFloat>const&>(dom);

    ScalarTaylorFunction stf=factory.create(dom, EffectiveScalarFunction::zero(dom.size()) );
    ARIADNE_TEST_PRINT(stf);
    ARIADNE_TEST_EQUALS(&stf.sweeper(),&sweeper);
    ARIADNE_TEST_EQUALS(stf.evaluate(args),ValidatedFloat(0.0));

    VectorTaylorFunction vtf=factory.create(dom, EffectiveVectorFunction::identity(dom.size()) );
    ARIADNE_TEST_PRINT(vtf);

    // Test evaluation gives a superset with small additional error
    Vector<ValidatedFloat> errs(2,ValidatedFloat(-1e-15,+1e-15));
    ARIADNE_TEST_BINARY_PREDICATE(refines,args,vtf.evaluate(args));
    ARIADNE_TEST_BINARY_PREDICATE(refines,vtf.evaluate(args),Vector<ValidatedFloat>(args+errs));
    Vector<ValidatedFloat> pt(2); pt[0]=ValidatedFloat(0.2); pt[1]=ValidatedFloat(1.25);
    ARIADNE_TEST_BINARY_PREDICATE(refines,pt,vtf.evaluate(pt));
}



int main() {
    ThresholdSweeper sweeper(std::numeric_limits<float>::epsilon());
    TestScalarTaylorFunction(sweeper).test();
    TestVectorTaylorFunction(sweeper).test();
    TestTaylorFunctionFactory().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
