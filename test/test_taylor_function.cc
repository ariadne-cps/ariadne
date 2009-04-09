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
#include "taylor_set.h"
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

namespace Ariadne {
std::pair<Float, Vector<Interval> > flow_bounds(FunctionInterface const&,Vector<Interval> const&,Float const&);
}
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
    void test_implicit();
    void test_join();
    void test_combine();
    void test_flow();
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
    ARIADNE_TEST_CALL(test_implicit());
    ARIADNE_TEST_CALL(test_join());
    ARIADNE_TEST_CALL(test_flow());
}


void TestTaylorFunction::test_constructors()
{
    Vector< Expansion<Float> > expansion(2);
    expansion[0]=Expansion<Float>(2,4, 0,0,1.125, 1,0,-0.75, 0,1,0.0625, 2,0,-0.25);
    expansion[1]=Expansion<Float>(2,2, 0,0,0.750, 1,0,0.50);

    Vector<Interval> domain(2,0.25,1.25,0.5,1.0);
    HenonFunction henon_function(Vector<Float>(2,1.5,-0.25));
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
    Expansion<Float> expansion2(1,2, 0,0.0, 1,1.0);
    Vector<Float> error2(1, 0.125);
    Vector<Interval> subdomain2(1, 3e-16, 1.0);
    Expansion<Float> subexpansion2(1,2, 0,0.50000000000000022, 1,0.49999999999999989);
    Vector<Float> suberror2(1, 0.12500000000000008);
    TaylorFunction function2(domain2,expansion2*e(1,0),error2);
    TaylorFunction restricted_function2(subdomain2,subexpansion2*e(1,0),suberror2);
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
    Expansion<Float> aexpansion30(2,3, 0.0, 0.5,0.0, 0.5,1.5,0.0,0.66666666666666663,1.25,3.0,0.0);
    Vector<Float> aerror30(1,5.5511151231257827e-17);
    TaylorFunction antiderivative30(domain3,aexpansion30*e(1,0),aerror30);
    ARIADNE_TEST_EQUAL(antiderivative(function3,index0),antiderivative30);
    Expansion<Float> aexpansion31(2,3, 0.0, 0.0,0.25, 0.0,0.5,0.375, 0.0,1.0,0.625,0.5);
    TaylorFunction antiderivative31(domain3,aexpansion31*e(1,0));
    ARIADNE_TEST_EQUAL(antiderivative(function3,index1),antiderivative31);

}

void TestTaylorFunction::test_implicit()
{
    // Test computation of solution of x^2/8-y/2=0 on [-1,+1]
    Vector<Interval> df1(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> dh1=project(df1,range(0,1));
    Vector< Polynomial<Float> > pf1=(0.125*p(2,0)-0.5*p(2,1))*e(1,0);
    Vector< Polynomial<Float> > ph1=0.25*p(1,0)*e(1,0);
    TaylorFunction f1(df1,pf1);
    TaylorFunction h1(dh1,ph1);
    ARIADNE_TEST_EQUAL(implicit(f1),h1);

    // Test computation of sqrt(4+x)-2 on [-1,+1] by solving 4+x-(y+2)*(y+2)=0
    Vector<Interval> df2(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> dh2=project(df2,range(0,1));
    Vector< Polynomial<Float> > pf2=(4.0+p(2,0)-(p(2,1)+2.0)*(p(2,1)+2.0))*e(1,0);
    TaylorFunction f2(df2,pf2);
    TaylorFunction i2=TaylorFunction::identity(dh2);
    TaylorFunction h2=implicit(f2);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(compose(f2,join(i2,h2)).range()),1e-6);

    // Test computation of sqrt(4+x)-2 on [0,3] by solving 4+x-(y+2)*(y+2)=0
    // Note that the solution lies in (-1,+1)
    Vector<Interval> df3(2, 0.0,+3.0, -1.0,+1.0);
    Vector<Interval> dh3=project(df3,range(0,1));
    Vector< Polynomial<Float> > pf3=(4.0+p(2,0)-(p(2,1)+2.0)*(p(2,1)+2.0))*e(1,0);
    TaylorFunction f3(df3,pf3);
    TaylorFunction i3=TaylorFunction::identity(dh3);
    TaylorFunction h3=implicit(f3);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(compose(f3,join(i3,h3)).range()),1e-6);
}


void TestTaylorFunction::test_flow()
{
    {
        Vector<Float> ex=e(1,0); 
        Polynomial<Float> x=p(1,0); 
        Polynomial<Float> x0=p(2,0); Polynomial<Float> t=p(2,1);

        // Test scaling of simple flow dx/dt=3 on the unit interval
        Vector<Interval> b(1u,Interval(-4,+4));
        Vector<Interval> d(1u,Interval(-2,0));
        Interval h(0,0.125);
        Vector<Polynomial<Float> > vfp=(0.0*x+0.5)*ex;
        Vector<Polynomial<Float> > flwp=(x0+0.5*t)*ex;
        TaylorFunction vf(b,vfp);
        TaylorFunction flw(join(d,h),flwp);
        //ARIADNE_TEST_PRINT(flwp);
        //ARIADNE_TEST_PRINT(flw.polynomial());
        //ARIADNE_TEST_PRINT(flw.polynomial()<<" "<<flwp);
        //ARIADNE_TEST_PRINT(flow(vf,d,h,4).polynomial());
        ARIADNE_TEST_EQUAL(flow(vf,d,h,4),flw);

        h=Interval(-0.125,0.125);
        flw=TaylorFunction(join(d,h),flwp);
        ARIADNE_TEST_EQUAL(flow(vf,d,h,4),flw);


    }

    {
        // Test solution of the affine flow
        //   dx/dt = Ax+b; x(0)=x0
        // with solution
        //   x=inv(A)(exp(At)-I)b+exp(At)x0

        // In scalar form, use specific flow
        //   dx/dt=ax+by+e; dy/dt=cx+dy+f; x(0)=x0; y(0)=y0;
        // The first terms in the expansion are
        //   x=(x0+e*t)+(a*(x0+t/2)+b*(y0+f/2))*t+((a*a+b*c)*x0+b*(a+d)*y0)*t*t
        //      +e*t+(a*e+b*f)*t*t/2+...
        Vector<Float> ex=e(2,0); Vector<Float> ey=e(2,1);
        Polynomial<Float> x=p(2,0); Polynomial<Float> y=p(2,1);
        Polynomial<Float> x0=p(3,0); Polynomial<Float> y0=p(3,1); Polynomial<Float> t=p(3,2);

        Float a=-0.75; Float b=-0.50; Float e=-1.0;
        Float c=+0.50; Float d=-0.75; Float f= 0.0;

        Vector<Interval> vf_domain(2, -1.0,+1.0, -1.0,+1.0);
        Vector< Polynomial<Float> > vf_poly=(e+a*x+b*y)*ex+(f+c*x+d*y)*ey;
        TaylorFunction vector_field(vf_domain,vf_poly);

        Vector<Interval> flow_domain(2, -0.5,+0.5, -0.5,+0.5);
        Vector< Polynomial<Float> > flow_poly=flow(vf_poly,6);
        Vector< Polynomial<Interval> > ivl_flow_poly=flow_poly;
        Vector<Interval> flow_error(2,Interval(-1e-1,1e-1)); // extra flow error


        Interval time(-0.25,0.25);
        uint order(6);
        TaylorFunction computed_flow;

        // Expect exception indicating that flow cannot be computed over this time step.
        ARIADNE_TEST_THROWS(flow(vector_field,flow_domain,time,order),FlowBoundsException);


        // Compute flow of vector field
        time/=2;
        computed_flow=flow(vector_field,flow_domain,time,order);
        TaylorFunction expected_flow(join(flow_domain,time),flow_poly);
        expected_flow+=flow_error;

        ARIADNE_TEST_PRINT(flow_poly);
        ARIADNE_TEST_PRINT(midpoint(truncate(computed_flow.polynomial(),6)));

        ARIADNE_TEST_PRINT(computed_flow);
        ARIADNE_TEST_PRINT(expected_flow);


        ARIADNE_TEST_PRINT(norm(Vector<TaylorModel>(expected_flow.models()-computed_flow.models())));
        for(uint i=0; i!=2; ++i) { const_cast<TaylorModel&>(computed_flow.models()[i]).sweep(1e-6); }
        ARIADNE_TEST_BINARY_PREDICATE(refines,computed_flow,expected_flow);
    }
/*
    // Test for contraction of single flow step
    // This is not strictly needed for a valid solution
    TaylorFunction contracted_flow=antiderivative(compose(vector_field,computed_flow),2);
    if(!refines(contracted_flow,computed_flow)) {
        ARIADNE_TEST_WARN("Perron-Frobenius operator is not a contraction on flow "<<
                            computed_flow<<"; operator yields "<<contracted_flow);
    }
*/

    {
        // Test a flow starting at a corner of the flow bounds with constant derivative 
        // This test is designed to catch problems due to the flow being computed on the
        // time interval [0,h], while intermediate steps use [-h,h]
        Vector< Polynomial< Float > > vector_field_poly = (1.5+p(1,0)*0.0)*e(1,0);
        Vector<Interval> bounding_box(1, Interval(-1,+1));
        Vector<Interval> initial_box(1, Interval(-1,-0.5));
        Interval time_interval(0,1.0);
        uint order=4;
        TaylorFunction vector_field(bounding_box,vector_field_poly);
        TaylorFunction flow_model=flow(vector_field,initial_box,time_interval,order);
        Vector< Polynomial<Float> > flow_poly = (p(2,0)+1.5*p(2,1))*e(1,0);
        ARIADNE_TEST_EQUAL(flow_model.polynomial(),flow_poly);
    }

    {
        // Test a flow over a small time interval and space domain
        // This test is designed to catch problems due to scaling of the domains

        // Vector field dtx=1.5
        Vector< Polynomial< Float > > vector_field_poly = (1.5+p(1,0)*0.0)*e(1,0);
        Vector<Interval> bounding_box(1, Interval(0.25,0.75));
        Vector<Interval> initial_box(1, Interval(0.375,0.50));
        Interval time_interval(0,0.03125);
        uint order=2;
        TaylorFunction vector_field(bounding_box,vector_field_poly);
        TaylorFunction flow_model=flow(vector_field,initial_box,time_interval,order);
        Vector< Polynomial<Float> > flow_poly = (p(2,0)+1.5*p(2,1))*e(1,0);
        ARIADNE_TEST_PRINT(vector_field_poly)
        ARIADNE_TEST_PRINT(flow_model.domain());
        ARIADNE_TEST_PRINT(flow_model.models());
        ARIADNE_TEST_EQUAL(flow_model.domain(),join(initial_box,time_interval));
        ARIADNE_TEST_EQUAL(flow_model.polynomial(),flow_poly);
    }


/*
    {
        // Test integration step starting from an initial set
        Vector<Float> ex=Vector<Float>::unit(2,0);
        Vector<Float> ey=Vector<Float>::unit(2,1);
        Polynomial<Float> x=Polynomial<Float>::variable(2,0);
        Polynomial<Float> y=Polynomial<Float>::variable(2,1);
        Polynomial<Float> x0=Polynomial<Float>::variable(3,0);
        Polynomial<Float> y0=Polynomial<Float>::variable(3,1);
        Polynomial<Float> t=p(3,2);
        TaylorSet initial_set_model(PolynomialFunction((0.264+0.0005*x)*ex+(0.046+0.0005*y)*ey), Vector<Interval>(2,Interval(-1,+1)));
        PolynomialFunction vector_field((0.3-0.02*x)*ex+(0.0*y)*ey);
        double step_size=0.125;
        TaylorModel integration_time_model=TaylorModel::constant(2,step_size);

        Vector<Interval> flow_domain=initial_set_model.range();
        Vector<Interval> flow_bound(2, 0.263538,0.264482, 0.0458876,0.0468876);
        TaylorFunction vector_field_model(flow_bound,vector_field);
        TaylorFunction flow_model=unchecked_flow(vector_field_model,flow_domain,Interval(0.0,step_size),6);
        TaylorSet final_set_model=apply(flow_model,TaylorSet(join(initial_set_model.models(),integration_time_model)));

        ARIADNE_TEST_PRINT(vector_field);
        ARIADNE_TEST_PRINT(initial_set_model);
        ARIADNE_TEST_PRINT(integration_time_model);
        ARIADNE_TEST_PRINT(flow_model);
        ARIADNE_TEST_PRINT(final_set_model);

        for(uint i=0; i!=100; ++i) {
            Vector<Interval> flow_domain=initial_set_model.range();
            Vector<Interval> flow_bound=flow_bounds(vector_field,flow_domain,step_size).second;
            TaylorFunction vector_field_model(flow_bound,vector_field);
            TaylorFunction flow_model=unchecked_flow(vector_field_model,flow_domain,Interval(0.0,step_size),6);
            final_set_model=apply(flow_model,TaylorSet(join(initial_set_model.models(),integration_time_model)));
            initial_set_model=final_set_model;
            std::cerr<<final_set_model<<"\n";
        }
        ARIADNE_TEST_PRINT(final_set_model);
    }
*/
/*
TaylorSet([({0,0;0:0.26401, 1,0;1:0.000472061},2.53714e-15),({0,0;0:0.0463876, 0,1;1:0.0005},3.48973e-16)])
flow_model=TaylorFunction( [[0.263538:0.264482],[0.0458876:0.0468876],[0:0.125]], [{ 0,0,0:[-6.25667e-16:-6.99588e-17], 1,0,0:[1:1], 0,0,1:[0.0175693:0.0175693], 1,0,1:[-0.02:-0.02], 0,0,2:[-0.000175693:-0.000175693], 1,0,2:[0.0002:0.0002], 0,0,3:[1.17128e-06:1.17128e-06], 1,0,3:[-1.33333e-06:-1.33333e-06], 0,0,4:[-5.85422e-09:-5.85422e-09], 1,0,4:[6.65833e-09:6.65833e-09], 0,0,5:[1.63649e-11:1.63649e-11] },{ 0,0,0:[-7.63278e-17:7.63278e-17], 0,1,0:[1:1], 0,0,1:[-0.0463876:-0.0463876] }] )
integration_time_model=({0,0;0:0.125},0)
final_set_model=TaylorSet([({0,0;0:0.265544, 1,0;1:0.000470882},2.71563e-15),({0,0;0:0.0405892, 0,1;1:0.0005},3.51577e-16)])
*/
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
