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
    void test_constructors();
    void test_variables();
};


void TestPolynomial::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_variables());
}


void TestPolynomial::test_constructors()
{
    Polynomial<Float> p(3);
}

void TestPolynomial::test_variables()
{
    Vector< Polynomial<Float> > x=variables<Float>(3);
    array< Vector<Float> > e=unit_vectors<Float>(2);

    ARIADNE_TEST_EQUAL((x[0]*(x[1]*3.0+x[0])+x[1]*x[2]),Polynomial<Float>(3,2, 0.,0.,0.,0., 1.,3.,0.,0.,1.,0.));
    ARIADNE_TEST_EQUAL((e[1]*(x[0]*(x[1]*3.0+x[0])+x[1]*x[2]))[1],Polynomial<Float>(3,2, 0.,0.,0.,0., 1.,3.,0.,0.,1.,0.));
    
}



int main() {
    TestPolynomial().test();
    return ARIADNE_TEST_FAILURES;
}
