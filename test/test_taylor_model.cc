/***************************************************************************
 *            test_taylor_variable.cc
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
#include "differential.h"
#include "taylor_variable.h"
#include "function.h"
#include "models.h"

#include "test.h"
using namespace std;
using namespace Ariadne;


class TestTaylorVariable
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
};


void TestTaylorVariable::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_predicates());
    ARIADNE_TEST_CALL(test_approximation());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_functions());
    ARIADNE_TEST_CALL(test_compose());
    ARIADNE_TEST_CALL(test_antiderivative());
}


void TestTaylorVariable::test_concept()
{
    const Float f=0.0;
    const Interval i;
    const Vector<Float> vf;
    const Vector<Interval> vi;
    const TaylorVariable  t;
    TaylorVariable tr;

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
    t.domain(); t.range(); t.polynomial(); t.error();
    
}

void TestTaylorVariable::test_constructors()
{
    double a[]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
    ARIADNE_TEST_CONSTRUCT(TaylorVariable,tv1,(2,3,a,0.25));
    ARIADNE_TEST_CONSTRUCT(TaylorVariable,tv2,(2,3, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 0.25));

    ARIADNE_ASSERT_EQUAL(tv1,tv2);
    ARIADNE_ASSERT_EQUAL(tv1.value(),1.0);
}

void TestTaylorVariable::test_predicates()
{
    TaylorVariable tv1(1,2, 1.00,2.00,3.00, 0.75);
    TaylorVariable tv2(1,2, 1.00,1.75,3.25, 0.25);
    TaylorVariable tv3(1,2, 1.125,1.75,3.25, 0.25);
    TaylorVariable tv4(1,3, 1.00,2.25,3.00, -0.25, 0.25);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv4,tv1);
}

void TestTaylorVariable::test_approximation()
{
    double a[]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
    ARIADNE_TEST_CONSTRUCT(TaylorVariable,tv1,(2,3,a,0.25));

    ARIADNE_TEST_CONSTRUCT(TaylorVariable,tv2,(1,2,1.0,2.0,3.0,0.25));
}

void TestTaylorVariable::test_evaluate()
{
    Vector<Interval> iv(2, 0.25,0.5, -0.75,-0.5);
    TaylorVariable tv(2,2,1.0,2.0,3.0,4.0,5.0,6.0,0.25);
    ARIADNE_TEST_EQUAL(evaluate(tv,iv),Interval(-1,1));
}

void TestTaylorVariable::test_arithmetic()
{
    //Operations which can be performed exactly with floating-point arithmetic.
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)+(-3), TaylorVariable(1,2, -2.0,-2.0,3.0, 0.75));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)-(-3), TaylorVariable(1,2, 4.0,-2.0,3.0, 0.75));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)*(-3), TaylorVariable(1,2, -3.0,6.0,-9.0, 2.25));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)/(-4), TaylorVariable(1,2, -0.25,0.5,-0.75, 0.1875));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)+Interval(-1,2), TaylorVariable(1,2, 1.5,-2.0,3.0, 2.25));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)-Interval(-1,2), TaylorVariable(1,2, 0.5,-2.0,3.0, 2.25));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)*Interval(-1,2), TaylorVariable(1,2, 0.5,-1.0,1.5, 10.5));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)/Interval(0.25,2.0), TaylorVariable(1,2, 2.25,-4.5,6.75, 13.5));
    ARIADNE_TEST_EQUAL(+TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75), TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75));
    ARIADNE_TEST_EQUAL(-TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75), TaylorVariable(1,2, -1.0,2.0,-3.0, 0.75));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)+TaylorVariable(1,2, 3.0,2.0,-4.0, 0.5), TaylorVariable(1,2, 4.0,0.0,-1.0, 1.25));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)-TaylorVariable(1,2, 3.0,2.0,-4.0, 0.5), TaylorVariable(1,2, -2.0,-4.0,7.0, 1.25));
    ARIADNE_TEST_EQUAL(TaylorVariable(1,2, 1.0,-2.0,3.0, 0.75)*TaylorVariable(1,2, 3.0,2.0,-4.0, 0.5), TaylorVariable(1,4, 3.0,-4.0,1.0,14.0,-12.0, 10.125));
}

void TestTaylorVariable::test_functions()
{
    TaylorVariable xz(1,1, 0.0, 0.5, 0.0);
    TaylorVariable xo(1,1, 1.0, 0.5, 0.0);

    //Functions based on their natural defining points
    ARIADNE_TEST_BINARY_PREDICATE(refines,exp(xz),TaylorVariable(1,6, 1.00000,0.50000,0.12500,0.02083,0.00260,0.00026,0.00002, 0.00003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sin(xz),TaylorVariable(1,6, 0.00000,0.50000,0.0000,-0.02083,0.00000,0.00026,0.00000, 0.00003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,cos(xz),TaylorVariable(1,6, 1.00000,0.0000,-0.12500,0.00000,0.00260,0.0000,-0.00002, 0.00003));

    ARIADNE_TEST_BINARY_PREDICATE(refines,rec(xo),TaylorVariable(1,6,  1.000000,-0.500000, 0.250000,-0.125000, 0.062500,-0.031250, 0.015625, 0.018));
    ARIADNE_TEST_BINARY_PREDICATE(refines,sqrt(xo),TaylorVariable(1,6, 1.000000, 0.250000,-0.031250, 0.007813,-0.002441, 0.000854,-0.000320, 0.0003));
    ARIADNE_TEST_BINARY_PREDICATE(refines,log(xo),TaylorVariable(1,6,  0.000000, 0.500000,-0.125000, 0.041667,-0.015625, 0.006250,-0.002604, 0.003));

}


void TestTaylorVariable::test_compose()
{
}


void TestTaylorVariable::test_antiderivative() 
{
    Interval unit_interval(-1,+1);
    TaylorVariable tm=TaylorVariable::constant(2,1.0);
    TaylorVariable atm=antiderivative(tm,unit_interval,1);
}



int main() {
    TestTaylorVariable().test();
    std::cout << "INCOMPLETE " << std::flush;
    return ARIADNE_TEST_FAILURES;
}
