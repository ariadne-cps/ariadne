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
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "expansion.h"
#include "taylor_model.h"
#include "taylor_function.h"
#include "differential.h"
#include "polynomial.h"
#include "function.h"
#include "models.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

Vector<Float> e(uint n, uint i) { return Vector<Float>::unit(n,i); }
Expansion<Float> v(uint n, uint j) { return Expansion<Float>::variable(n,j); }
Polynomial<Float> p(uint n, uint j) { return Polynomial<Float>::variable(n,j); }
ScalarTaylorFunction t(Vector<Interval> d, uint j) { return ScalarTaylorFunction::variable(d,j); }

namespace Ariadne {
std::pair<Float, Vector<Interval> > flow_bounds(RealVectorFunction const&,Vector<Interval> const&,Float const&);
typedef VectorUserFunction<Henon> HenonFunction;
}



class TestScalarTaylorFunction
{
  public:
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
    Vector<Interval> d(unsigned int n) { return Vector<Interval>(n,Interval(-1,+1)); }
    typedef Expansion<Float> e;
};


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
    const Float f=0.0;
    const Interval i;
    const Vector<Float> vf;
    const Vector<Interval> vi;
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

    tr.sweep(); tr.truncate(); tr.clean();

    t.evaluate(vi); evaluate(t,vi);
    t.domain(); t.range(); t.expansion(); t.error();

}

void TestScalarTaylorFunction::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv1,(d(2),e(2,3, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0), 0.25));

    ARIADNE_ASSERT_EQUAL(tv1.domain(),Vector<Interval>(2,Interval(-1,+1)));
    ARIADNE_ASSERT_EQUAL(tv1.argument_size(),2);
    ARIADNE_ASSERT_EQUAL(tv1.number_of_nonzeros(),10);
    ARIADNE_ASSERT_EQUAL(tv1.value(),1.0);
    ARIADNE_ASSERT_EQUAL(tv1.error(),0.25);
}

void TestScalarTaylorFunction::test_predicates()
{
    ScalarTaylorFunction tv1(d(1),e(1,2, 1.00,2.00,3.00), 0.75);
    ScalarTaylorFunction tv2(d(1),e(1,2, 1.00,1.75,3.25), 0.25);
    ScalarTaylorFunction tv3(d(1),e(1,2, 1.125,1.75,3.25), 0.25);
    ScalarTaylorFunction tv4(d(1),e(1,3, 1.00,2.25,3.00,-0.25), 0.25);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

void TestScalarTaylorFunction::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv1,(d(2),e(2,3,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0),0.25));
    ARIADNE_TEST_CONSTRUCT(ScalarTaylorFunction,tv2,(d(1),e(1,2,1.0,2.0,3.0),0.25));
}

void TestScalarTaylorFunction::test_evaluate()
{
    Vector<Interval> iv(2, 0.25,0.5, -0.75,-0.5);
    ScalarTaylorFunction tv(d(2),e(2,1.0,2.0,3.0,4.0,5.0,6.0),0.25);
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),Interval(-1,1));
}

void TestScalarTaylorFunction::test_arithmetic()
{
    ARIADNE_TEST_EQUAL(d(1),d(1));
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)+(-3), ScalarTaylorFunction(d(1),e(1,2, -2.0,-2.0,3.0), 0.75));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)-(-3), ScalarTaylorFunction(d(1),e(1,2, 4.0,-2.0,3.0), 0.75));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)*(-3), ScalarTaylorFunction(d(1),e(1,2, -3.0,6.0,-9.0), 2.25));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)/(-4), ScalarTaylorFunction(d(1),e(1,2, -0.25,0.5,-0.75), 0.1875));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)+Interval(-1,2), ScalarTaylorFunction(d(1),e(1,2, 1.5,-2.0,3.0), 2.25));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)-Interval(-1,2), ScalarTaylorFunction(d(1),e(1,2, 0.5,-2.0,3.0), 2.25));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)*Interval(-1,2), ScalarTaylorFunction(d(1),e(1,2, 0.5,-1.0,1.5), 10.5));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)/Interval(0.25,2.0), ScalarTaylorFunction(d(1),e(1,2, 2.25,-4.5,6.75), 13.5));
    ARIADNE_TEST_EQUAL(+ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75), ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75));
    ARIADNE_TEST_EQUAL(-ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75), ScalarTaylorFunction(d(1),e(1,2, -1.0,2.0,-3.0), 0.75));

    // Regression test to check subtraction yielding zero coefficients
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)+ScalarTaylorFunction(d(1),e(1,2, 3.0,2.0,-4.0), 0.5), ScalarTaylorFunction(d(1),e(1,2, 4.0,0.0,-1.0), 1.25));

    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)-ScalarTaylorFunction(d(1),e(1,2, 3.0,2.0,-4.0), 0.5), ScalarTaylorFunction(d(1),e(1,2, -2.0,-4.0,7.0), 1.25));
    ARIADNE_TEST_EQUAL(ScalarTaylorFunction(d(1),e(1,2, 1.0,-2.0,3.0), 0.75)*ScalarTaylorFunction(d(1),e(1,2, 3.0,2.0,-4.0), 0.5), ScalarTaylorFunction(d(1),e(1,4, 3.0,-4.0,1.0,14.0,-12.0), 10.125));

}

void TestScalarTaylorFunction::test_functions()
{
    ScalarTaylorFunction xz(d(1),e(1,1, 0.0, 0.5), 0.0);
    ScalarTaylorFunction xo(d(1),e(1,1, 1.0, 0.5), 0.0);

    //Functions based on their natural defining points
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),ScalarTaylorFunction(d(1),e(1,6, 1.00000,0.50000,0.12500,0.02083,0.00260,0.00026,0.00002), 0.00003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),ScalarTaylorFunction(d(1),e(1,6, 0.00000,0.50000,0.0000,-0.02083,0.00000,0.00026,0.00000), 0.00003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),ScalarTaylorFunction(d(1),e(1,6, 1.00000,0.0000,-0.12500,0.00000,0.00260,0.0000,-0.00002), 0.00003));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),ScalarTaylorFunction(d(1),e(1,6,  1.000000,-0.500000, 0.250000,-0.125000, 0.062500,-0.031250, 0.015625), 0.018));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),ScalarTaylorFunction(d(1),e(1,6, 1.000000, 0.250000,-0.031250, 0.007813,-0.002441, 0.000854,-0.000320), 0.0003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),ScalarTaylorFunction(d(1),e(1,6,  0.000000, 0.500000,-0.125000, 0.041667,-0.015625, 0.006250,-0.002604), 0.003));

}


void TestScalarTaylorFunction::test_compose()
{
}


void TestScalarTaylorFunction::test_antiderivative()
{
    ScalarTaylorFunction tm=ScalarTaylorFunction::constant(d(2),1.0);
    ScalarTaylorFunction atm=antiderivative(tm,1u);
}

void TestScalarTaylorFunction::test_substitute()
{
    Vector<Interval> d1(1, -0.75,+0.75);
    ScalarTaylorFunction tu=ScalarTaylorFunction::coordinate(d1,0);

    Vector<Interval> d2(2, -0.75,+0.75, -0.25,+0.75);
    ScalarTaylorFunction tx=ScalarTaylorFunction::coordinate(d2,0);
    ScalarTaylorFunction ty=ScalarTaylorFunction::coordinate(d2,1);

    ScalarTaylorFunction tf=1+3*tx*ty;
    ScalarTaylorFunction ts=0.25+tu*tu;

    ScalarTaylorFunction tg=substitute(tf,1,ts);
    ScalarTaylorFunction tr=1+0.75*tu+3*tu*tu*tu;
    ARIADNE_TEST_EQUAL(tg,tr);
}

void TestScalarTaylorFunction::test_conversion() {
    // Test conversion between ordinary functions and Taylor functions.
    Vector<Interval> D(2, -0.5,0.5, -1.0,2.0);
    Vector<Float> pt(2, -0.25,0.25);
    Vector<Interval> ipt(pt);
    RealVectorFunction x=RealVectorFunction::identity(2);

    RealScalarFunction f=(1-x[0]*x[0]-0.5*x[1]);
    ScalarTaylorFunction tf(D,f);

    ARIADNE_TEST_PRINT(f);
    ARIADNE_TEST_PRINT(tf);

    // Conversion to TaylorFunction should be exact in second component
    ARIADNE_TEST_BINARY_PREDICATE(subset,f(ipt),tf(ipt));
    ARIADNE_TEST_BINARY_PREDICATE(subset,tf(ipt),f(ipt)+Interval(-1e-15,1e-15));
}


/*
VectorTaylorFunction henon(const VectorTaylorFunction& x, const Vector<Float>& p)
{
    VectorTaylorFunction r(2,2,x.degree()); henon(r,x,p); return r;
}
*/

class TestVectorTaylorFunction
{
  public:
    TestVectorTaylorFunction();
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


TestVectorTaylorFunction::TestVectorTaylorFunction()
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


void TestVectorTaylorFunction::test_constructors()
{
    Vector< Expansion<Float> > expansion(2);
    expansion[0]=Expansion<Float>(2,4, 0,0,1.125, 1,0,-0.75, 0,1,0.0625, 2,0,-0.25);
    expansion[1]=Expansion<Float>(2,2, 0,0,0.750, 1,0,0.50);
    expansion[0].reverse_lexicographic_sort(); expansion[1].reverse_lexicographic_sort();

    Vector<Interval> domain(2,0.25,1.25,0.5,1.0);
    HenonFunction henon_function(Vector<Float>(2,1.5,0.25));
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,henon_model,(domain,henon_function));
    ARIADNE_TEST_EQUAL(henon_model.models()[0].expansion(),expansion[0])
    ARIADNE_TEST_EQUAL(henon_model.models()[1].expansion(),expansion[1])

    Vector<Float> e0=e(2,0); Vector<Float> e1=e(2,1);
    Polynomial<Float> x=p(2,0); Polynomial<Float> y=p(2,1);
    Vector< Polynomial<Float> > polynomial=(1.5-x*x+0.25*y)*e0+x*e1;
    ARIADNE_TEST_CONSTRUCT(VectorTaylorFunction,polynomial_model,(domain,polynomial));
    ARIADNE_TEST_EQUAL(polynomial_model,VectorTaylorFunction(domain,expansion))

    VectorTaylorFunction t=VectorTaylorFunction::identity(domain);
    //VectorTaylorFunction variables_model((1.5-t[0]*t[0]+0.25*t[1])*e0+t[0]*e1);
    VectorTaylorFunction variables_model(ScalarTaylorFunction(1.5-t[0]*t[0]+0.25*t[1])*e0+ScalarTaylorFunction(t[0])*e1);
    variables_model.sweep();
    ARIADNE_TEST_EQUAL(variables_model,VectorTaylorFunction(domain,expansion));

}

void TestVectorTaylorFunction::test_restrict()
{
    Vector<Interval> domain1(2, -1.0,+1.0, -1.0,+1.0);
    Expansion<Float> expansion1(2,3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
    Vector<Interval> subdomain1(2, -0.25,0.75, -0.5,0.0);
    Expansion<Float> subexpansion1(2,3, 1.031250, 1.812500,0.625000, 1.812500,0.562500,0.0468750,
                                             0.875000,0.500000,0.281250,0.156250);
    VectorTaylorFunction function1(domain1,expansion1*e(1,0));
    VectorTaylorFunction restricted_function1(subdomain1,subexpansion1*e(1,0));
    ARIADNE_TEST_EQUAL(restrict(function1,subdomain1),restricted_function1);

    Vector<Interval> domain2(1, -1.0,+1.0);
    Expansion<Float> expansion2(1,2, 0,0.0, 1,1.0);
    Vector<Float> error2(1, 0.125);
    Vector<Interval> subdomain2(1, 3e-16, 1.0);
    Expansion<Float> subexpansion2(1,2, 0,0.50000000000000022, 1,0.49999999999999989);
    Vector<Float> suberror2(1, 0.12500000000000008);
    VectorTaylorFunction function2(domain2,expansion2*e(1,0),error2);
    VectorTaylorFunction restricted_function2(subdomain2,subexpansion2*e(1,0),suberror2);
    ARIADNE_TEST_EQUAL(restrict(function2,subdomain2),restricted_function2);
}

void TestVectorTaylorFunction::test_jacobian()
{
    HenonFunction henon(Vector<Float>(2,1.5,-0.25));
    Vector<Interval> domain1(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> domain2(2, -0.5,+0.5, -0.25,+0.25);
    Vector<Interval> domain3(2, -0.25,+0.75, 0.0,+0.50);
    Vector<Float> point1(2, 0.0,0.0);
    Vector<Float> point2(2, 0.5,0.25);
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain1,henon).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain1,henon).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain2,henon).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain2,henon).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain3,henon).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(VectorTaylorFunction(domain3,henon).jacobian(point2),henon.jacobian(point2));
}

void TestVectorTaylorFunction::test_compose()
{
    Float a=1.5; Float b=0.25;
    Polynomial<Float> x=p(2,0);
    Polynomial<Float> y=p(2,1);
    Vector< Polynomial<Float> > henon_polynomial=(a-x*x+b*y)*e(2,0)+x*e(2,1);
    Vector< Polynomial<Float> > henon_square_polynomial=
        (a*(1-a)+b*x-2*a*b*y+2*a*x*x-b*b*y*y+2*b*x*x*y-x*x*x*x)*e(2,0)
            + (a-x*x+b*y)*e(2,1);
    //    compose(henon_polynomial,henon_polynomial);
    Vector<Interval> domain1(2,0.25,1.25,0.5,1.0);
    VectorTaylorFunction function1(domain1,henon_polynomial);
    Vector<Interval> domain2(2, -1.5,2.5, 0.25,1.25);
    VectorTaylorFunction function2(domain2,henon_polynomial);

    VectorTaylorFunction composition1(domain1,henon_square_polynomial);
    ARIADNE_TEST_EQUAL(compose(function2,function1),composition1);
}


void TestVectorTaylorFunction::test_antiderivative()
{
    unsigned int index0=0;
    unsigned int index1=1;

    Vector<Interval> domain1(2,Interval(-1,1));
    Expansion<Float> expansion1(2,0, 3.0);
    VectorTaylorFunction function1(domain1,expansion1*e(1,0));
    Expansion<Float> aexpansion1(2,1, 0.0, 0.0,3.0);
    VectorTaylorFunction antiderivative1(domain1,aexpansion1*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function1,index1),antiderivative1);

    Vector<Interval> domain2(2, -0.25,0.75, 0.0,0.5);
    Expansion<Float> expansion2(2,0, 3.0);
    VectorTaylorFunction function2(domain2,expansion2*e(1,0));
    Expansion<Float> aexpansion2(2,1, 0.0, 0.0,0.75);
    VectorTaylorFunction antiderivative2(domain2,aexpansion2*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function2,index1),antiderivative2);

    Vector<Interval> domain3(2, -0.25,0.75, 0.0,0.5);
    Expansion<Float> expansion3(2,2, 1.0,2.0,3.0,4.0,5.0,6.0);
    VectorTaylorFunction function3(domain3,expansion3*e(1,0));
    Expansion<Float> aexpansion30(2,3, 0.0, 0.5,0.0, 0.5,1.5,0.0,0.66666666666666663,1.25,3.0,0.0);
    Vector<Float> aerror30(1,5.5511151231257827e-17);
    VectorTaylorFunction antiderivative30(domain3,aexpansion30*e(1,0),aerror30);
    ARIADNE_TEST_EQUAL(antiderivative(function3,index0),antiderivative30);
    Expansion<Float> aexpansion31(2,3, 0.0, 0.0,0.25, 0.0,0.5,0.375, 0.0,1.0,0.625,0.5);
    VectorTaylorFunction antiderivative31(domain3,aexpansion31*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function3,index1),antiderivative31);

}

void TestVectorTaylorFunction::test_join()
{
    Vector<Interval> domain(2, -0.25,+0.25, -0.5,+0.5);
    Vector< Polynomial<Float> > polynomial1 = (p(2,0)*p(2,0)+2.0*p(2,0)*p(2,1)+3.0*p(2,1)*p(2,1))*e(1,0);
    Vector< Polynomial<Float> > polynomial2 = (4.0*p(2,0)*p(2,0)+5.0*p(2,0)*p(2,1)+6.0*p(2,1)*p(2,1))*e(2,1);
    Vector< Polynomial<Float> > polynomial3 = (p(2,0)*p(2,0)+2.0*p(2,0)*p(2,1)+3.0*p(2,1)*p(2,1))*e(3,0)
        + (4.0*p(2,0)*p(2,0)+5.0*p(2,0)*p(2,1)+6.0*p(2,1)*p(2,1))*e(3,2);
    VectorTaylorFunction function1(domain,polynomial1);
    VectorTaylorFunction function2(domain,polynomial2);
    VectorTaylorFunction function3(domain,polynomial3);
    ARIADNE_TEST_EQUAL(join(function1,function2),function3);

}

void TestVectorTaylorFunction::test_combine()
{
    // This test contains a regression test to check correct behaviour for a zero component.
    Vector<Interval> domain1(2, -0.25,+0.25, -0.5,+0.5);
    Vector<Interval> domain2(3, -0.75,+0.75, -1.0,+1.0, -1.25,+1.25);
    Vector<Interval> domain3(5, -0.25,+0.25, -0.5,+0.5, -0.75,+0.75, -1.0,+1.0, -1.25,+1.25);
    Vector< Polynomial<Float> > polynomial1 = (p(2,0)*p(2,0)+2.0*p(2,0)*p(2,1)+3.0*p(2,1)*p(2,1))*e(1,0);
    Vector< Polynomial<Float> > polynomial2 = (4.0*p(3,0)*p(3,0)+5.0*p(3,0)*p(3,1)+6.0*p(3,1)*p(3,2))*e(2,1);
    Vector< Polynomial<Float> > polynomial3 = (p(5,0)*p(5,0)+2.0*p(5,0)*p(5,1)+3.0*p(5,1)*p(5,1))*e(3,0)
        + (4.0*p(5,2)*p(5,2)+5.0*p(5,2)*p(5,3)+6.0*p(5,3)*p(5,4))*e(3,2);
    VectorTaylorFunction function1(domain1,polynomial1);
    VectorTaylorFunction function2(domain2,polynomial2);
    VectorTaylorFunction function3(domain3,polynomial3);
    ARIADNE_TEST_EQUAL(combine(function1,function2),function3);

}

void TestVectorTaylorFunction::test_conversion()
{
    // Test conversion between ordinary functions and Taylor functions.
    Vector<Interval> D(2, -0.5,0.5, -1.0,2.0);
    Vector<Float> pt(2, -0.25,0.25);
    Vector<Interval> ipt(pt);
    RealVectorFunction x=RealVectorFunction::identity(2);

    RealVectorFunction h=RealVectorFunction((1-x[0]*x[0]-0.5*x[1],x[0]+Real(0)));
    VectorTaylorFunction th(D,h);

    ARIADNE_TEST_PRINT(h);
    ARIADNE_TEST_PRINT(th);

    // Conversion to TaylorFunction should be exact in second component
    ARIADNE_TEST_EQUAL(th(pt)[1],h(pt)[1]);
    ARIADNE_TEST_EQUAL(th(ipt)[1],h(ipt)[1]);
    ARIADNE_TEST_BINARY_PREDICATE(subset,h[0](ipt),th[0](ipt));


}

// Regression test for domain with empty interior
void TestVectorTaylorFunction::test_domain()
{
    RealScalarFunction z=RealScalarFunction::constant(2,0.0);
    RealScalarFunction o=RealScalarFunction::constant(2,1.0);
    RealScalarFunction x0=RealScalarFunction::coordinate(2,0);
    RealScalarFunction x1=RealScalarFunction::coordinate(2,1);

    Vector<Interval> D1(2, -1.0,1.0, -1.0,1.0);
    VectorTaylorFunction t1(D1, (o,x0+x1));
    ARIADNE_TEST_PRINT(t1);
    ARIADNE_TEST_PRINT(t1.codomain());
    Vector<Interval> D2(2, 1.0,1.0, -2.0,2.0);
    ScalarTaylorFunction t2(D2,2*x0+x1*x1);
    ARIADNE_TEST_PRINT(t2.domain());
    ARIADNE_TEST_PRINT(t2.model());
    ARIADNE_TEST_PRINT(t2.codomain());

    ARIADNE_TEST_PRINT(t2);
    ARIADNE_TEST_PRINT(compose(t2,t1));
    ScalarTaylorFunction t3(D1,2+(x0+x1)*(x0+x1));
    ARIADNE_TEST_EQUAL(compose(t2,t1),t3);

    Vector<Interval> x(2, 1.0,1.0, 0.5,1.5);
    ARIADNE_TEST_PRINT(x);
    ARIADNE_TEST_EQUAL(evaluate(t2,x),Interval(2.25,4.25));

    // Ensure evaluation and composition throw errors when expected
    Vector<Interval> xe(2, 0.875,1.125, 0.5,1.5);
    ARIADNE_TEST_THROWS(evaluate(t2,xe),DomainException);

    VectorTaylorFunction te1=t1; te1[0]=te1[0]+Interval(-0.125,+0.125);
    ARIADNE_TEST_THROWS(compose(t2,te1),DomainException);

    ARIADNE_TEST_EQUAL(unchecked_evaluate(t2,xe),Interval(2.25,4.25));
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
    TaylorModelAccuracy accuracy(1e-4,3u);
    TaylorFunctionFactory factory(accuracy);

    IntervalVector dom(2, -1,+1, 0.5,3.5);

    ScalarTaylorFunction stf=factory.create_zero(dom);
    ARIADNE_TEST_PRINT(stf);
    ARIADNE_TEST_EQUALS(*stf.accuracy_ptr(),accuracy);
    ARIADNE_TEST_EQUALS(stf.evaluate(dom),Interval(0.0));

    VectorTaylorFunction vtf=factory.create_identity(dom);
    ARIADNE_TEST_PRINT(vtf);

    // Test evaluation gives a superset with small additional error
    ARIADNE_TEST_BINARY_PREDICATE(subset,dom,vtf.evaluate(dom));
    ARIADNE_TEST_BINARY_PREDICATE(subset,vtf.evaluate(dom),dom+IntervalVector(2,Interval(-1e-15,+1e-15)));
    IntervalVector pt(2); pt[0]=Interval(0.2); pt[1]=Interval(1.25);
    ARIADNE_TEST_BINARY_PREDICATE(subset,pt,vtf.evaluate(pt));
}



int main() {
    TestScalarTaylorFunction().test();
    TestVectorTaylorFunction().test();
    TestTaylorFunctionFactory().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
