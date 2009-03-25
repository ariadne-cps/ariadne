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
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "polynomial.h"

#include "test.h"
using namespace std;
using namespace Ariadne;


class TestPolynomial
{
  public:
    void test();
  private:
    void test_concept();
    void test_iterator_concept();
    void test_constructors();
    void test_variables();
};


void TestPolynomial::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_variables());
}


void TestPolynomial::test_concept()
{
    Float x=0;
    Vector<Float> v(3);
    MultiIndex a(3);
    Polynomial<Float> p(3);
    const Polynomial<Float> cp(3);

    p=Polynomial<Float>();
    p=Polynomial<Float>(3);
    p=Polynomial<Float>(cp);

    p=Polynomial<Float>(3,1, 0.0, 0.0,0.0,0.0);
    p=Polynomial<Float>(3,1, 1,2,3,5.0);

    //p=Polynomial<Float>::variable(3u,0u);
    //p=Polynomial<Float>::variables(3u)[0u];

    p=x;

    p.reserve(2u);
    p.insert(a,x);
    p.append(a,x);
    p.append(a,a,x);

    x=cp[a];
    p[a]=1.0;

    cp.argument_size();
    cp.evaluate(v);

    p.erase(p.begin());

    p.check();
    p.cleanup();
    p.clear();
}

void TestPolynomial::test_iterator_concept()
{
    MultiIndex a(3);
    Polynomial<Float> p(3);
    const Polynomial<Float> cp(3);

    Polynomial<Float>::iterator iter=p.begin(); iter=p.end(); iter=p.find(a);
    Polynomial<Float>::const_iterator citer=p.begin(); citer=p.end(); citer=p.find(a);
    citer=p.begin(); citer=cp.end(); citer=cp.find(a);

    Polynomial<Float>::value_type val=*iter;
    Polynomial<Float>::reference ref=*iter;
    //Polynomial<Float>::pointer ptr=iter.operator->();
    Polynomial<Float>::const_reference ncref=*iter;

    // WARNING: Cannot convert non-constant pointer to constant pointer
    //Polynomial<Float>::const_pointer ncptr=iter.operator->();

    Polynomial<Float>::value_type cval=*citer;
    Polynomial<Float>::const_reference cref=*citer;
    //Polynomial<Float>::const_pointer cptr=citer.operator->();

    ++iter; --iter;
    ++citer; --citer;

    iter==iter; iter!=iter; citer==citer; citer!=citer;
    citer==iter; citer!=iter; iter==citer; iter!=citer;

}

void TestPolynomial::test_constructors()
{
    ARIADNE_TEST_CONSTRUCT(Polynomial<Float>,p1,(3));
    ARIADNE_TEST_CONSTRUCT(Polynomial<Float>,p2,(3,2, 0.,0.,0.,0., 1.,3.,0.,0.,1.,0.));
    ARIADNE_TEST_CONSTRUCT(Polynomial<Float>,p3,(2,5, 1,2,5.0, 0,0,2.0, 1,0,3.0, 3,0,7.0, 0,1,11.0));
    //TODO: Case for repeated index
}

void TestPolynomial::test_variables()
{
    Vector< Polynomial<Float> > x=Polynomial<Float>::variables(3);
    array< Vector<Float> > e=Vector<Float>::units(2);

    ARIADNE_TEST_EVALUATE(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]);
    ARIADNE_TEST_EQUAL((x[0]*(x[1]*3.0+x[0])+x[1]*x[2]), Polynomial<Float>(3,3, 1,1,0,3.0, 2,0,0,1.0, 0,1,1,1.0));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[1], Polynomial<Float>(3,3, 1,1,0,3.0, 2,0,0,1.0, 0,1,1,1.0));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[0], Polynomial<Float>(3));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[0], Polynomial<Float>(3,0, 0.0));

}



int main() {
    TestPolynomial().test();
    return ARIADNE_TEST_FAILURES;
}
