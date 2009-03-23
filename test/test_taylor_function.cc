/***************************************************************************
 *            test_polynomial.cc
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
#include "taylor_variable.h"
#include "taylor_function.h"
#include "function.h"
#include "polynomial.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

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
    void test_misc();
};


TestTaylorFunction::TestTaylorFunction()
{
  std::cout<<std::setprecision(17);
  std::cerr<<std::setprecision(17);
}


void
TestTaylorFunction::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_jacobian());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_flow());
    ARIADNE_TEST_CALL(test_implicit());
    ARIADNE_TEST_CALL(test_join());
    ARIADNE_TEST_CALL(test_combine());
    ARIADNE_TEST_CALL(test_misc());
}


void TestTaylorFunction::test_constructors()
{
    HenonFunction henon_function(Vector<Float>(2,1.5,-0.25));
    Vector<Interval> domain(2,0.25,1.25,0.5,1.0);
    Vector<TaylorVariable> expansion(2,2,2,
        1.125, -0.75,0.0625, -0.25,0.00,0.00,   0.0,
        0.750,  0.50,0.0000,  0.00,0.00,0.00,   0.0);
    ARIADNE_TEST_CONSTRUCT(TaylorFunction,henon_model,(domain,henon_function));
    ARIADNE_TEST_EQUAL(henon_model,TaylorFunction(domain,expansion))
}

void TestTaylorFunction::test_restrict()
{
    {
        Vector<Interval> domain(2, -1.0,+1.0, -1.0,+1.0);
        Vector<TaylorVariable> expansion(1,2,3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 0.0);
        Vector<Interval> subdomain(2, -0.25,0.75, -0.5,0.0);
        Vector<TaylorVariable> subexpansion(1,2,3, 1.031250, 1.812500,0.625000, 1.812500,0.562500,0.0468750,
                                                0.875000,0.500000,0.281250,0.156250,  0.0);
        TaylorFunction function(domain,expansion);
        TaylorFunction restricted_function(subdomain,subexpansion);
        ARIADNE_TEST_EQUAL(restrict(function,subdomain),restricted_function);
    }
    
    {
        Vector<Interval> domain(1, -1.0,+1.0);
        Vector<TaylorVariable> expansion(1,1,1, 0.0, 1.0, 0.0);
        Vector<Interval> subdomain(1, 1e-16, 1.0);
        Vector<TaylorVariable> subexpansion(1,1,1, 0.50000000000000000,0.49999999999999994, 1.6653345369377348e-16);
        TaylorFunction function(domain,expansion);
        TaylorFunction restricted_function(subdomain,subexpansion);
        ARIADNE_TEST_EQUAL(restrict(function,subdomain),restricted_function);
    }
       
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
    HenonFunction henon_function(Vector<Float>(2,1.5,0.25));
    Vector<Interval> domain1(2,0.0,1.25,0.5,1.0);
    TaylorFunction henon_model1(domain1,henon_function);
    
    Vector<Interval> domain2=henon_model1.range();
    TaylorFunction henon_model2(domain2,henon_function);

    TaylorFunction henon_square=compose(henon_model2,henon_model1);
}


void TestTaylorFunction::test_antiderivative() 
{
    Vector<Interval> domain(3,Interval(-1,1));
    TaylorFunction tm=TaylorFunction::constant(domain,Vector<Float>(2,1.0,2.0));
    TaylorFunction atm=antiderivative(tm,1);
}

void TestTaylorFunction::test_implicit() 
{
    Vector<Interval> df1(2, -1.0,+1.0, -1.0,+1.0);
    Vector<Interval> dh1=project(df1,range(0,1));
    Vector<TaylorVariable> ef1(1,2,2, 0.0, 0.125,-0.50, 0.0,0.0,0.0,  0.0);
    Vector<TaylorVariable> eh1(1,1,2, 0.0, 0.25, 0.0,  0.0);
    Vector<TaylorVariable> e2(1,2,2, 2.0, 0.25,-0.50, 0.0,0.0,-0.0625,  0.0);
    Vector<TaylorVariable> e3(1,2,2, 2.0, 0.25,-1.0, 0.0,0.0,-1.0,  0.0);
    TaylorFunction f1(df1,ef1);
    TaylorFunction h1(dh1,eh1);
    ARIADNE_ASSERT(implicit(f1)==h1);
}


void TestTaylorFunction::test_flow() 
{
    Vector< TaylorVariable > f(2,2,1, -1.0,-0.8,-0.4, 0.0, 0.0,-0.4,-0.8, 0.0);
    TaylorFunction vector_field(Vector<Interval>(2,Interval(-2,2)),f);
    Vector<Interval> domain(2,Interval(-1,1));
    Interval time(-0.25,0.25);
    cout << flow(vector_field,domain,time) << std::endl << std::endl;
}

void TestTaylorFunction::test_join() 
{
    
}

void TestTaylorFunction::test_combine() 
{
    
}

void TestTaylorFunction::test_misc()
{
    MultiIndex a(4);
    for(uint i=0; i!=100; ++i) { cout << a << "\n"; ++a; }

    a=2*MultiIndex::unit(4,2)+MultiIndex::unit(4,1); cout << a << "\n";

    TaylorVariable y(2);
    a=MultiIndex::zero(2);   y[a]=2.0;
    a=MultiIndex::unit(2,0); y[a]=1.0;
    a=MultiIndex::unit(2,1); y[a]=1.0;
    TaylorFunction x(Vector<Interval>::unit_box(2),Vector< TaylorVariable >(1,y));
    a=MultiIndex::zero(2);   y[a]=3.0;
    a=MultiIndex::unit(2,0); y[a]=0.0;
    a=MultiIndex::unit(2,1); y[a]=1.0;
    cout << x << endl;
    cout << x+x << endl;
}


int main() {
    TestTaylorFunction().test();
    std::cerr<<"INCOMPLETE "<<std::flush;
    return ARIADNE_TEST_FAILURES;
}
