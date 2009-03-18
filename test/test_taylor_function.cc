/***************************************************************************
 *            test_taylor_function.cc
 *
 *  Copyright 2008  Pieter Collins
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
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "taylor_variable.h"
#include "taylor_function.h"
#include "function.h"
#include "models.h"

#include "test.h"
using namespace std;
using namespace Ariadne;

template<class R, class A, class P>
void henon(R& r, const A& x, const P& p) 
{
    r[0]=-(x[0]*x[0])+p[0]-p[1]*x[1];
    r[1]=x[0];
}

template<class R, class A>
void spiral(R& r, const A& x) 
{
    r[0]=-0.8*x[0]+0.4*x[1]-1.0;
    r[1]=-0.4*x[0]-0.8*x[1];
}

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
    void test();
  private:
    void test_constructors();
    void test_compose();
    void test_antiderivative();
    void test_flow();
    void test_implicit();
    void test_join();
    void test_combine();
    void test_restrict();
    void test_misc();
};


void TestTaylorFunction::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_flow());
    ARIADNE_TEST_CALL(test_implicit());
    ARIADNE_TEST_CALL(test_join());
    ARIADNE_TEST_CALL(test_combine());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_misc());
}


void TestTaylorFunction::test_constructors()
{
    HenonFunction henon_function(Vector<Float>(2,1.5,0.25));
    Vector<Interval> domain(2,0.0,1.25,0.5,1.0);
    Vector<TaylorVariable> expansion(2,2,2, 0.0,0.0,0.0,0.0,0.0,0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0, 0.0);

    ARIADNE_TEST_CONSTRUCT(TaylorFunction,henon_model,(domain,henon_function));
    ARIADNE_TEST_EQUAL(henon_model,TaylorFunction(domain,expansion))
}


void TestTaylorFunction::test_compose()
{
    HenonFunction henon_function(Vector<Float>(2,1.5,0.25));
    Vector<Interval> domain1(2,0.0,1.25,0.5,1.0);
    TaylorFunction henon_model1(domain1,henon_function);
    std::cerr<<henon_model1<<std::endl;
    
    Vector<Interval> domain2=henon_model1.range();
    TaylorFunction henon_model2(domain2,henon_function);

    std::cerr<<henon_model2<<std::endl;
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
    Vector<Interval> d(2, -1.0,+1.0, -1.0,+1.0);
    Vector<TaylorVariable> e(1,2,2, 0.0,0.0,0.0,1.0,0.0,1.0, 0.0);
    TaylorFunction tf(d,e);
    cout << "implicit(x)=" << implicit(tf) << endl;
    cout << endl;
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

void TestTaylorFunction::test_restrict() 
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
    return ARIADNE_TEST_FAILURES;
}
