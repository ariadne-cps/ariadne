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
#include "taylor_variable.h"
#include "taylor_function.h"
#include "function.h"
#include "polynomial.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

Vector<Float> e(uint n, uint i) { return Vector<Float>::unit(n,i); }
Expansion<Float> v(uint n, uint j) { return Expansion<Float>::variable(n,j); }
Polynomial<Float> p(uint n, uint j) { return Polynomial<Float>::variable(n,j); }
TaylorVariable t(Vector<Interval> d, uint j) { return TaylorVariable::variable(d,j); }

struct Henon {
    uint result_size() const { return 2; }
    uint argument_size() const { return 2; }
    uint parameter_size() const { return 2; }
    uint smoothness() const { return 255; }
    template<class R, class A, class P>
    void compute(R& r, const A& x, const P& p) const
    {
        r[0]=-(x[0]*x[0])+p[0]-p[1]*x[1];
        r[1]=x[0];
    }
};
typedef Function<Henon> HenonFunction;

/*
TaylorFunction henon(const TaylorFunction& x, const Vector<Float>& p)
{
    TaylorFunction r(2,2,x.degree()); henon(r,x,p); return r;
}
*/

class TestTaylorFunction
{
  public:
    TestTaylorFunction();
    void test();
  private:
    void test_constructors();
    void test_restrict();
    void test_jacobian();
    void test_compose();
    void test_antiderivative();
    void test_flow();
    void test_implicit();
    void test_join();
    void test_combine();
};


TestTaylorFunction::TestTaylorFunction()
{
  std::cout<<std::setprecision(17);
  std::cerr<<std::setprecision(17);
}


void
TestTaylorFunction::test()
{
    ARIADNE_TEST_CALL(test_combine());
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_jacobian());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_flow());
    ARIADNE_TEST_CALL(test_implicit());
    ARIADNE_TEST_CALL(test_join());
}


void TestTaylorFunction::test_constructors()
{
    HenonFunction henon_function(Vector<Float>(2,1.5,-0.25));
    Vector<Interval> domain(2,0.25,1.25,0.5,1.0);
    Vector< Expansion<Float> > expansion(2);
    expansion[0]=Expansion<Float>(1,2,2, 1.125, -0.75,0.0625, -0.25,0.00,0.00);
    expansion[1]=Expansion<Float>(1,2,2, 0.750,  0.50,0.0000,  0.00,0.00,0.00);

    ARIADNE_TEST_CONSTRUCT(TaylorFunction,henon_model,(domain,henon_function));
    ARIADNE_TEST_EQUAL(henon_model.models()[0].expansion(),expansion[0])
    ARIADNE_TEST_EQUAL(henon_model.models()[1].expansion(),expansion[1])

    Vector<Float> e0=e(2,0); Vector<Float> e1=e(2,1);
    Polynomial<Float> x=p(2,0); Polynomial<Float> y=p(2,1);
    Vector< Polynomial<Float> > polynomial=(1.5-x*x+0.25*y)*e0+x*e1;
    ARIADNE_TEST_CONSTRUCT(TaylorFunction,polynomial_model,(domain,polynomial));
    ARIADNE_TEST_EQUAL(polynomial_model,TaylorFunction(domain,expansion))

    Vector<TaylorExpression> t=TaylorExpression::variables(domain);
    TaylorFunction variables_model = (1.5-t[0]*t[0]+0.25*t[1])*e0+t[0]*e1;
    ARIADNE_TEST_EQUAL(variables_model,TaylorFunction(domain,expansion));

}

void TestTaylorFunction::test_restrict()
{
    Vector<Interval> domain1(2, -1.0,+1.0, -1.0,+1.0);
    Expansion<Float> expansion1(2,3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
    Vector<Interval> subdomain1(2, -0.25,0.75, -0.5,0.0);
    Expansion<Float> subexpansion1(2,3, 1.031250, 1.812500,0.625000, 1.812500,0.562500,0.0468750,
                                             0.875000,0.500000,0.281250,0.156250);
    TaylorFunction function1(domain1,expansion1*e(1,0));
    TaylorFunction restricted_function1(subdomain1,subexpansion1*e(1,0));
    ARIADNE_TEST_EQUAL(restrict(function1,subdomain1),restricted_function1);

    Vector<Interval> domain2(1, -1.0,+1.0);
    Expansion<Float> expansion2(1,1, 0.0, 1.0);
    Vector<Interval> subdomain2(1, 1e-16, 1.0);
    Expansion<Float> subexpansion2(1,1, 0.50000000000000000,0.49999999999999994);
    Vector<Float> error2(1, 1.6653345369377348e-16);
    TaylorFunction function2(domain2,expansion2*e(1,0),error2);
    TaylorFunction restricted_function2(subdomain2,subexpansion2*e(1,0),error2);
    ARIADNE_TEST_EQUAL(restrict(function2,subdomain2),restricted_function2);
}

void TestTaylorFunction::test_jacobian()
{
    HenonFunction henon(Vector<Float>(2,1.5,-0.25));
    Vector<Interval> domain1(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> domain2(2, -0.5,+0.5, -0.25,+0.25);
    Vector<Interval> domain3(2, -0.25,+0.75, 0.0,+0.50);
    Vector<Float> point1(2, 0.0,0.0);
    Vector<Float> point2(2, 0.5,0.25);
    ARIADNE_TEST_EQUAL(TaylorFunction(domain1,henon).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(TaylorFunction(domain1,henon).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(TaylorFunction(domain2,henon).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(TaylorFunction(domain2,henon).jacobian(point2),henon.jacobian(point2));
    ARIADNE_TEST_EQUAL(TaylorFunction(domain3,henon).jacobian(point1),henon.jacobian(point1));
    ARIADNE_TEST_EQUAL(TaylorFunction(domain3,henon).jacobian(point2),henon.jacobian(point2));
}

void TestTaylorFunction::test_compose()
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
    TaylorFunction function1(domain1,henon_polynomial);
    Vector<Interval> domain2(2, -1.5,2.5, 0.25,1.25);
    TaylorFunction function2(domain2,henon_polynomial);

    TaylorFunction composition1(domain1,henon_square_polynomial);
    ARIADNE_TEST_EQUAL(compose(function2,function1),composition1);
}


void TestTaylorFunction::test_antiderivative()
{
    unsigned int index0=0;
    unsigned int index1=1;

    Vector<Interval> domain1(2,Interval(-1,1));
    Expansion<Float> expansion1(2,0, 3.0);
    TaylorFunction function1(domain1,expansion1*e(1,0));
    Expansion<Float> aexpansion1(2,1, 0.0, 0.0,3.0);
    TaylorFunction antiderivative1(domain1,aexpansion1*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function1,index1),antiderivative1);

    Vector<Interval> domain2(2, -0.25,0.75, 0.0,0.5);
    Expansion<Float> expansion2(2,0, 3.0);
    TaylorFunction function2(domain2,expansion2*e(1,0));
    Expansion<Float> aexpansion2(2,1, 0.0, 0.0,0.75);
    TaylorFunction antiderivative2(domain2,aexpansion2*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function2,index1),antiderivative2);

    Vector<Interval> domain3(2, -0.25,0.75, 0.0,0.5);
    Expansion<Float> expansion3(2,2, 1.0,2.0,3.0,4.0,5.0,6.0);
    TaylorFunction function3(domain3,expansion3*e(1,0));
    Expansion<Float> aexpansion30(1,2,3, 0.0, 0.5,0.0, 0.5,1.5,0.0,0.66666666666666663,1.25,3.0,0.0);
    Vector<Float> aerror30(1,5.5511151231257827e-17);
    TaylorFunction antiderivative30(domain3,aexpansion30*e(1,0),aerror30);
    ARIADNE_TEST_EQUAL(antiderivative(function3,index0),antiderivative30);
    Expansion<Float> aexpansion31(2,3, 0.0, 0.0,0.25, 0.0,0.5,0.375, 0.0,1.0,0.625,0.5);
    TaylorFunction antiderivative31(domain3,aexpansion31*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function3,index1),antiderivative31);

}

void TestTaylorFunction::test_implicit()
{
    Vector<Interval> df1(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> dh1=project(df1,range(0,1));
    Vector< Polynomial<Float> > pf1=(0.125*p(2,0)-0.5*p(2,1))*e(1,0);
    Vector< Polynomial<Float> > ph1=0.25*p(1,0)*e(1,0);
    TaylorFunction f1(df1,pf1);
    TaylorFunction h1(dh1,ph1);
    ARIADNE_TEST_EQUAL(implicit(f1),h1);

    // Test computation of sqrt(4+x)-2 on [-1,+1] by solving 4+x-(y+2)*(y+2)=0
    Vector<Interval> df2(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> dh2(1, -1.0,+1.0);
    Vector< Polynomial<Float> > pf2=(4.0+p(2,0)-(p(2,1)+2.0)*(p(2,1)+2.0))*e(1,0);
    TaylorFunction f2(df2,pf2);
    TaylorFunction i2=TaylorFunction::identity(dh2);
    TaylorFunction h2=implicit(f2);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(compose(f2,join(i2,h2)).range()),1e-6);
}


void TestTaylorFunction::test_flow()
{
    Float a=-0.75; Float b=-0.5; Float c=-1.0;
    Vector<Float> ex=e(2,0); Vector<Float> ey=e(2,1);
    Polynomial<Float> x=p(2,0); Polynomial<Float> y=p(2,1);
    Vector< Polynomial<Float> > vfp=(c+a*x-b*y)*ex+(b*x+a*y)*ey;

    x=p(3,0); y=p(3,1); Polynomial<Float> t=p(3,2);
    Vector< Polynomial<Float> > flp=
        (1.0*x+(-1.0-0.75*x+0.5*y)*t+(0.375+0.15625*x-0.375*y)*t*t
            +(-0.052083+0.023437*x+0.119791*y)*t*t*t)*ex
        +(1.0*y+(-0.5*x-0.75*y)*t+(0.25+0.375*x+0.15625*y)*t*t
            +(-0.125-0.119791*x+0.023437*y)*t*t*t)*ey;
    Vector<Interval> fle(2,Interval(-1e-3,1e-3)); // extra flow error

    Vector<Interval> d(2, -1.0,+1.0, -1.0,+1.0);
    TaylorFunction vector_field(d,vfp);
    Vector<Interval> domain(2, -0.5,+0.5, -0.5,+0.5);
    Interval time(-0.25,0.25);
    uint order(6);
    TaylorFunction computed_flow;
    ARIADNE_TEST_TRY(
        computed_flow=flow(vector_field,domain,time,order);
        TaylorFunction expected_flow(join(domain,time),flp);
        expected_flow+=Vector<Interval>(2,Interval(-1e-3,1e-3));
        ARIADNE_TEST_BINARY_PREDICATE(refines,computed_flow,expected_flow);

        TaylorFunction flow_error=computed_flow-antiderivative(compose(vector_field,computed_flow),2);
        for(uint i=0; i!=2; ++i) { const_cast<TaylorModel&>(flow_error.models()[i]).sweep(1e-4); }
        ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow_error.range()),1e-4);
    );
}


void TestTaylorFunction::test_join()
{
    Vector<Interval> domain(2, -0.25,+0.25, -0.5,+0.5);
    Vector< Polynomial<Float> > polynomial1 = (p(2,0)*p(2,0)+2.0*p(2,0)*p(2,1)+3.0*p(2,1)*p(2,1))*e(1,0);
    Vector< Polynomial<Float> > polynomial2 = (4.0*p(2,0)*p(2,0)+5.0*p(2,0)*p(2,1)+6.0*p(2,1)*p(2,1))*e(2,1);
    Vector< Polynomial<Float> > polynomial3 = (p(2,0)*p(2,0)+2.0*p(2,0)*p(2,1)+3.0*p(2,1)*p(2,1))*e(3,0)
        + (4.0*p(2,0)*p(2,0)+5.0*p(2,0)*p(2,1)+6.0*p(2,1)*p(2,1))*e(3,2);
    TaylorFunction function1(domain,polynomial1);
    TaylorFunction function2(domain,polynomial2);
    TaylorFunction function3(domain,polynomial3);
    ARIADNE_TEST_EQUAL(join(function1,function2),function3);

}

void TestTaylorFunction::test_combine()
{
    // This test contains a regression test to check correct behaviour for a zero component.
    Vector<Interval> domain1(2, -0.25,+0.25, -0.5,+0.5);
    Vector<Interval> domain2(3, -0.75,+0.75, -1.0,+1.0, -1.25,+1.25);
    Vector<Interval> domain3(5, -0.25,+0.25, -0.5,+0.5, -0.75,+0.75, -1.0,+1.0, -1.25,+1.25);
    Vector< Polynomial<Float> > polynomial1 = (p(2,0)*p(2,0)+2.0*p(2,0)*p(2,1)+3.0*p(2,1)*p(2,1))*e(1,0);
    Vector< Polynomial<Float> > polynomial2 = (4.0*p(3,0)*p(3,0)+5.0*p(3,0)*p(3,1)+6.0*p(3,1)*p(3,2))*e(2,1);
    Vector< Polynomial<Float> > polynomial3 = (p(5,0)*p(5,0)+2.0*p(5,0)*p(5,1)+3.0*p(5,1)*p(5,1))*e(3,0)
        + (4.0*p(5,2)*p(5,2)+5.0*p(5,2)*p(5,3)+6.0*p(5,3)*p(5,4))*e(3,2);
    TaylorFunction function1(domain1,polynomial1);
    TaylorFunction function2(domain2,polynomial2);
    TaylorFunction function3(domain3,polynomial3);
    ARIADNE_TEST_EQUAL(combine(function1,function2),function3);

}


int main() {
    TestTaylorFunction().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
