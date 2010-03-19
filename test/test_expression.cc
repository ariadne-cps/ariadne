/***************************************************************************
 *            test_expression.cc
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
#include "formula.h"
#include "real.h"

#include "test.h"

using namespace Ariadne;

class TestExpression {
    RealConstant o;
    RealVariable x,y,z;
  public:
    TestExpression()
        : o("1.0",1.0), x("x"), y("y"), z("z") {
    }

    void test_variables() {
        ARIADNE_TEST_CONSTRUCT(RealVariable,a,("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable("a"));
        ARIADNE_TEST_ASSERT(a==RealVariable(a));
        ARIADNE_TEST_ASSERT(a!=RealVariable("b"));
    }

    void test_evaluate() {
        ARIADNE_TEST_CONSTRUCT(RealExpression,g,(x+3*y*z*z));

        Map<RealVariable,Real> v;
        v[x]=2.0; v[y]=3.0; v[z]=5.0;

        ARIADNE_TEST_PRINT(v);
        //ARIADNE_TEST_EQUAL(evaluate(g,v),Real(227));
    }

    void test_derivative() {
        ARIADNE_TEST_CONSTRUCT(RealExpression,g,(x+3*y*z*z));
        ARIADNE_TEST_PRINT(derivative(g,y));
        ARIADNE_TEST_PRINT(derivative(g,z));
        //ARIADNE_TEST_EQUAL(derivative(g,y),3*z*z);
        //ARIADNE_TEST_EQUAL(derivative(g,z),6*y*z);
    }

    void test() {
        test_variables();
        test_derivative();
    }
};


int main() {
    TestExpression().test();
    return ARIADNE_TEST_FAILURES;
};