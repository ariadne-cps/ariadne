/***************************************************************************
 *            test_taylor_model.cc
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
#include <iomanip>
#include "config.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "taylor_model.h"
#include "differential.h"
#include "function.h"
#include "polynomial.h"

#include "test.h"
using std::cout; using std::cerr; using std::endl;
using namespace Ariadne;

Vector<Float> v(uint n, uint i) { return Vector<Float>::unit(n,i); }
IntervalTaylorModel ctm(uint m, double c, Sweeper swp) { return IntervalTaylorModel::constant(m,c,swp); }
IntervalTaylorModel ctm(uint m, Sweeper swp) { return IntervalTaylorModel::constant(m,1.0,swp); }
//IntervalTaylorModel tm(uint m, uint i, Sweeper swp) { return IntervalTaylorModel::variable(m,i,swp); }


class TestTaylorModel
{
    typedef MultiIndex MI;
    typedef Expansion<Float> E;
    typedef Polynomial<Float> P;
    typedef IntervalTaylorModel T;
  public:
    Sweeper swp;
  public:
    TestTaylorModel(Sweeper swp);
    void test();
  private:
    void test_concept();
    void test_constructors();
    void test_predicates();
    void test_approximation();
    void test_unscale();
    void test_evaluate();
    void test_arithmetic();
    void test_range();
    void test_functions();
    void test_rescale();
    void test_restrict();
    void test_intersection();
    void test_split();
    void test_antiderivative();
    void test_compose();
};


TestTaylorModel::TestTaylorModel(Sweeper sweeper)
    : swp(sweeper)
{
}

void TestTaylorModel::test()
{
    std::cerr<<std::setprecision(17);
    std::cout<<std::setprecision(17);
    std::clog<<std::setprecision(17);

    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_predicates());
    ARIADNE_TEST_CALL(test_approximation());
    ARIADNE_TEST_CALL(test_unscale());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_range());
    ARIADNE_TEST_CALL(test_functions());
    ARIADNE_TEST_CALL(test_rescale());
    ARIADNE_TEST_CALL(test_restrict());
    ARIADNE_TEST_CALL(test_intersection());
    ARIADNE_TEST_CALL(test_split());
    ARIADNE_TEST_CALL(test_antiderivative());
    ARIADNE_TEST_CALL(test_compose());
}


void TestTaylorModel::test_concept()
{
    const Float f=0.0;
    const Interval i;
    const Vector<Float> vf;
    const Vector<Interval> vi;
    const IntervalTaylorModel  t(0,swp);
    IntervalTaylorModel tr(0,swp);

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

    evaluate(t,vi);
    t.domain(); t.range(); t.expansion(); t.error();

}

void TestTaylorModel::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(IntervalTaylorModel,tv1,(E(2,3, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0), 0.25, swp));
    ARIADNE_TEST_CONSTRUCT(IntervalTaylorModel,tv2,(E(2,10, 0,0,1.0, 1,0,2.0, 0,1,3.0, 2,0,4.0, 1,1,5.0, 0,2,6.0, 3,0,7.0, 2,1,8.0, 1,2,9.0, 0,3,10.0), 0.25, swp));

    ARIADNE_ASSERT_EQUAL(tv1.value(),1.0);
    ARIADNE_ASSERT_EQUAL(tv1.error(),0.25);
    ARIADNE_ASSERT_EQUAL(tv1.norm(),55.25);

    ARIADNE_ASSERT_EQUAL(tv2,tv1);
}

void TestTaylorModel::test_predicates()
{
    IntervalTaylorModel tv1(E(1,2, 1.00,2.00,3.00), 0.75, swp);
    IntervalTaylorModel tv2(E(1,2, 1.00,1.75,3.25), 0.25, swp);
    IntervalTaylorModel tv3(E(1,2, 1.125,1.75,3.25), 0.25, swp);
    IntervalTaylorModel tv4(E(1,3, 1.00,2.25,3.00,-0.25), 0.25, swp);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

void TestTaylorModel::test_approximation()
{
    ARIADNE_TEST_CONSTRUCT(IntervalTaylorModel,tv2,(E(1,2,1.0,2.0,3.0),0.25, swp));
}

void TestTaylorModel::test_unscale()
{
    if(unscale(IntervalTaylorModel(E(1,0, 3.0),0.0,swp),Interval(1.0)).range()!=Interval(1.0)) {
        ARIADNE_TEST_WARN("Unscaling over singleton domain does not yield constant");
    }
}

void TestTaylorModel::test_evaluate()
{
    Vector<Interval> iv(2, 0.25,0.5, -0.75,-0.5);
    IntervalTaylorModel tv(E(2,2,1.0,2.0,3.0,4.0,5.0,6.0),0.25,swp);
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),Interval(-1,1));
}

void TestTaylorModel::test_arithmetic()
{
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)+(-3), IntervalTaylorModel(E(1,2, -2.0,-2.0,3.0), 0.75,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)-(-3), IntervalTaylorModel(E(1,2, 4.0,-2.0,3.0), 0.75,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)*(-3), IntervalTaylorModel(E(1,2, -3.0,6.0,-9.0), 2.25,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)/(-4), IntervalTaylorModel(E(1,2, -0.25,0.5,-0.75), 0.1875,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)+Interval(-1,2), IntervalTaylorModel(E(1,2, 1.5,-2.0,3.0), 2.25,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)-Interval(-1,2), IntervalTaylorModel(E(1,2, 0.5,-2.0,3.0), 2.25,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)*Interval(-1,2), IntervalTaylorModel(E(1,2, 0.5,-1.0,1.5), 10.5,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)/Interval(0.25,2.0), IntervalTaylorModel(E(1,2, 2.25,-4.5,6.75), 13.5,swp));
    ARIADNE_TEST_EQUAL(+IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp), IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp));
    ARIADNE_TEST_EQUAL(-IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp), IntervalTaylorModel(E(1,2, -1.0,2.0,-3.0), 0.75,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)+IntervalTaylorModel(E(1,2, 3.0,2.0,-4.0),0.5,swp), IntervalTaylorModel(E(1,2, 4.0,0.0,-1.0), 1.25,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)-IntervalTaylorModel(E(1,2, 3.0,2.0,-4.0),0.5,swp), IntervalTaylorModel(E(1,2, -2.0,-4.0,7.0), 1.25,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 0.0,0.0,3.0), 0.75,swp)*IntervalTaylorModel(E(1,2, 3.0,2.0,-4.0),0.5,swp), IntervalTaylorModel(E(1,4, 0.0,0.0,9.0,6.0,-12.0), 8.625,swp));
    ARIADNE_TEST_EQUAL(IntervalTaylorModel(E(1,2, 1.0,-2.0,3.0), 0.75,swp)*IntervalTaylorModel(E(1,2, 3.0,2.0,-4.0),0.5,swp), IntervalTaylorModel(E(1,4, 3.0,-4.0,1.0,14.0,-12.0), 10.125,swp));

    IntervalTaylorModel tm_inf(E(2),+infty,swp);
    if(isnan(numeric_cast<double>((tm_inf * 0.0).error()))) {
        ARIADNE_TEST_WARN("Multiplying 0+/-infty by 0 yields 0+/-NaN");
    } else if((tm_inf * 0.0).error()==+infty) {
        ARIADNE_TEST_WARN("Multiplying 0+/-inf by 0 yields 0+/-inf");
    } else if((tm_inf * 0.0).error()==0.0) {
        ARIADNE_TEST_PRINT("Multiplying 0+/-inf by 0 yields 0+/-0");
    }
}

void TestTaylorModel::test_range()
{
    IntervalTaylorModel x0 = IntervalTaylorModel::variable(2,0,swp);
    IntervalTaylorModel x1 = IntervalTaylorModel::variable(2,1,swp);

    // Test range of cubic, which should be exact
    IntervalTaylorModel t1 = x0*x0*x0+x0;
    ARIADNE_TEST_ASSERT(t1.range()==Interval(-2,+2));

    // Test range of quadratic, which could be exact, but need not be
    IntervalTaylorModel t2 = x0*x0+x0;
    ARIADNE_TEST_BINARY_PREDICATE(subset,t2.range(),Interval(-2,+2));
    ARIADNE_TEST_BINARY_PREDICATE(subset,Interval(-0.25,+2),t2.range());
    if(!subset(t2.range(),Interval(-0.2578125,2.0))) { ARIADNE_TEST_WARN("IntervalTaylorModel::range() not exact for quadratic functions."); }
}

void TestTaylorModel::test_functions()
{
    IntervalTaylorModel x(E(1,1, 0.0, 1.0), 0.0, swp);
    IntervalTaylorModel xz(E(1,1, 0.0, 0.5), 0.0, swp);
    IntervalTaylorModel xo(E(1,1, 1.0, 0.5), 0.0, swp);

    ARIADNE_TEST_PRINT(exp(x));
    ARIADNE_TEST_PRINT(sin(x));
    ARIADNE_TEST_PRINT(cos(x));

    //Functions based on their natural defining points with variable dependence 0.5
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(T(E(1,1,0.0,1.0),0.0,swp)),T(E(1,6, 1.0,1.00,0.500,0.1667,0.0417,0.0083,0.0014),0.0004,swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(x),T(E(1,6, 0.0,1.0000,0.0,-0.1667,0.0,0.0083,0.0),0.0003,swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(x),T(E(1,6, 1.0000,0.0,-0.5000,0.0,0.0417,0.0,-0.0014),0.0003,swp));

    //Functions based on their natural defining points with variable dependence 0.5
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),IntervalTaylorModel(E(1,6, 1.00000,0.50000,0.12500,0.02083,0.00260,0.00026,0.00002), 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),IntervalTaylorModel(E(1,6, 0.00000,0.50000,0.0000,-0.02083,0.00000,0.00026,0.00000), 0.00003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),IntervalTaylorModel(E(1,6, 1.00000,0.0000,-0.12500,0.00000,0.00260,0.0000,-0.00002), 0.00003, swp));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),IntervalTaylorModel(E(1,6,  1.000000,-0.500000, 0.250000,-0.125000, 0.062500,-0.031250, 0.015625), 0.018, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),IntervalTaylorModel(E(1,6, 1.000000, 0.250000,-0.031250, 0.007813,-0.002441, 0.000854,-0.000320), 0.0003, swp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),IntervalTaylorModel(E(1,6,  0.000000, 0.500000,-0.125000, 0.041667,-0.015625, 0.006250,-0.002604), 0.003, swp));

    // Test exponential based at log2
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(T(E(1,1,0.693147,0.5),0.0,swp)),
                                  T(E(1,6, 2.00000,1.00000,0.25000,0.04166,0.00520,0.00052,0.00004), 0.00006, swp));

}


void TestTaylorModel::test_rescale()
{
}

void TestTaylorModel::test_restrict()
{
}

void TestTaylorModel::test_intersection()
{
    IntervalTaylorModel x=IntervalTaylorModel::variable(2,0,swp);
    IntervalTaylorModel y=IntervalTaylorModel::variable(2,1,swp);
    IntervalTaylorModel e=IntervalTaylorModel::error(2,1.0,swp);

    // Test intersection with no roundoff errors
    ARIADNE_TEST_EQUAL(intersection(T(E(1,4, 1.0,-0.75,0.0,3.0,3.25),0.5,swp),T(E(1,4, 1.0,0.0,0.25,2.0,3.0),1.0,swp)),
        T(E(1,4, 1.0,-0.625,0.0,2.75,3.25),0.50,swp));

    // Test intersection with roundoff errors
    ARIADNE_TEST_EQUAL(intersection(T(E(1,0, 2./3),0.5,swp),T(E(1,0, 6./5),0.25,swp)),
        T(E(1,0, 1.0583333333333331),0.10833333333333339,swp));
}

void TestTaylorModel::test_split()
{
    IntervalTaylorModel x=IntervalTaylorModel::variable(2,0,swp);
    IntervalTaylorModel y=IntervalTaylorModel::variable(2,1,swp);
    IntervalTaylorModel z=IntervalTaylorModel::zero(2,swp);
    IntervalTaylorModel t=1+3*x+2*y-5*x*x-7*x*y+11*y*y;
    IntervalTaylorModel es1=-1.75+4*x+5.5*y-1.25*x*x-3.5*x*y+11*y*y;
    IntervalTaylorModel es2=1.25-1*x-1.5*y-1.25*x*x-3.5*x*y+11*y*y;
    IntervalTaylorModel es3=1+1.5*x+2*y-1.25*x*x-3.5*x*y+11*y*y;

    ARIADNE_TEST_PRINT(t);
    ARIADNE_TEST_EQUAL(split(t,0).first,es1);
    ARIADNE_TEST_EQUAL(split(t,0).second,es2);
    ARIADNE_TEST_EQUAL(split(t,0,false),es1);
    ARIADNE_TEST_EQUAL(split(t,0,true),es2);
    ARIADNE_TEST_EQUAL(split(t,0,indeterminate),es3);

    // Test to make sure split(tm,k) and split(tm,k,lr) use same code
    ARIADNE_TEST_EQUAL(split(t,0).first,split(t,0,false));
    ARIADNE_TEST_EQUAL(split(t,0).second,split(t,0,true));
}


void TestTaylorModel::test_antiderivative()
{
    Interval unit_interval(-1,+1);
    IntervalTaylorModel tm=IntervalTaylorModel::constant(2,1.0,swp);
    IntervalTaylorModel atm=antiderivative(tm,1);

    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 0,0,2.0),0.,swp),0),T(E(2,1, 1,0,2.0),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 0,0,2.0),0.,swp),1),T(E(2,1, 0,1,2.0),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 1,0,3.0),0.,swp),0),T(E(2,1, 2,0,1.5),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 2,0,7.5),0.,swp),0),T(E(2,1, 3,0,2.5),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 2,4,7.5),0.,swp),0),T(E(2,1, 3,4,2.5),0.,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(2,1, 2,4,7.5),0.,swp),1),T(E(2,1, 2,5,1.5),0.,swp));

    // Test error control
    ARIADNE_TEST_EQUAL(antiderivative(T(E(1,1, 2,2.0),0.,swp),0),T(E(1,1, 3,0.66666666666666663),5.5511151231257827021e-17,swp));
    ARIADNE_TEST_EQUAL(antiderivative(T(E(1,1, 2,2.0),1.,swp),0),T(E(1,1, 3,0.66666666666666663),1.0000000000000002,swp));

    // Regression test for
    T t1=T(E(2,6, 0,0,1., 1,0,2., 0,1,3., 2,0,4., 1,1,5., 0,2,6.), 0., swp);
    T at1=T(E(2,6, 1,0,1., 2,0,1., 1,1,3., 3,0,1.33333333333333333, 2,1,2.5, 1,2,6.), 1.1102230246251565404e-16, swp);
    ARIADNE_TEST_EQUAL(antiderivative(t1,0),at1);
}


void TestTaylorModel::test_compose()
{

}



namespace Ariadne {
Vector<Interval> range(const Vector<IntervalTaylorModel>& tm) {
    Vector<Interval> r(tm.size()); for(uint i=0; i!=tm.size(); ++i) { r[i]=tm[i].range(); } return r; }
}

int main() {
    ThresholdSweeper sweeper(1e-8);
    TestTaylorModel(sweeper).test();

    return ARIADNE_TEST_FAILURES;
}
