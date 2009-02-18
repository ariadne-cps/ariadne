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
#include "sparse_differential.h"
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
    void test_evaluate();
    void test_predicates();
    void test_arithmetic();
    void test_functions();
    void test_compose();
    void test_antiderivative();
};


void TestTaylorVariable::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_predicates());
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
}

void TestTaylorVariable::test_constructors()
{
    double a[]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
    ARIADNE_TEST_CONSTRUCT(TaylorVariable,tv1,(2,3,a,0.25));

    ARIADNE_TEST_CONSTRUCT(TaylorVariable,tv2,(1,2,1.0,2.0,3.0,0.25));
}

void TestTaylorVariable::test_predicates()
{
    TaylorVariable tv1(1,2, 1.00,2.00,3.00, 0.75);
    TaylorVariable tv2(1,2, 1.00,1.75,3.25, 0.25);
    TaylorVariable tv3(1,2, 1.125,1.75,3.25, 0.25);

    ARIADNE_TEST_BINARY_PREDICATE(refines,tv1,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(refines,tv2,tv1);
    ARIADNE_TEST_BINARY_PREDICATE(!refines,tv3,tv1);
}

void TestTaylorVariable::test_arithmetic()
{
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
