/***************************************************************************
 *            test_procedure.cc
 *
 *  Copyright  2010  Pieter Collins
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>

#include "procedure.h"

#include "numeric.h"
#include "vector.h"
#include "formula.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestProcedure
{
    //static Formula<Float> o;;
    //static Formula<Float> x;
    //static Formula<Float> y;
  public:
    TestProcedure();
    void test();
  private:
    void test_formula();
    void test_construct_from_formula();
    void test_construct_from_expansion();
    void test_evaluate();
};

//Formula<Float> TestProcedure::o(Formula<Float>::constant(1.0));
//Formula<Float> TestProcedure::x(Formula<Float>::coordinate(0));
//Formula<Float> TestProcedure::y(Formula<Float>::coordinate(1));


TestProcedure::TestProcedure()
{
}

void TestProcedure::test()
{
    ARIADNE_TEST_CALL(test_formula());
    ARIADNE_TEST_CALL(test_construct_from_formula());
    ARIADNE_TEST_CALL(test_construct_from_expansion());
    ARIADNE_TEST_CALL(test_evaluate());
}

void TestProcedure::test_formula()
{
    Formula<Float> o(Formula<Float>::constant(1.0));
    Formula<Float> x(Formula<Float>::coordinate(0));
    Formula<Float> y(Formula<Float>::coordinate(1));

    Formula<Float> r(sqrt(pow(x,2)+pow(y,2)));
    ARIADNE_TEST_PRINT(r);

    //Vector< Formula<Float> > f((sqrt(pow(x,2)+pow(y,2)), atan(y/x)));
}

void TestProcedure::test_construct_from_formula()
{
    Formula<Float> o(Formula<Float>::constant(1.0));
    Formula<Float> x(Formula<Float>::coordinate(0));
    Formula<Float> y(Formula<Float>::coordinate(1));

    Formula<Float> xs=pow(x,2);
    Formula<Float> ys=pow(y,2);

    Vector< Formula<Float> > f((sqrt(xs+ys), atan(y/x), xs-ys));
    ARIADNE_TEST_PRINT(f);

    Procedure<Float> p0(f[0]);
    ARIADNE_TEST_PRINT(p0);
    Vector< Procedure<Float> > p(f);
    ARIADNE_TEST_PRINT(p);

    p0+=Float(5);
    ARIADNE_TEST_PRINT(p0);
}

void TestProcedure::test_construct_from_expansion()
{
    Expansion<Float> e(2, 4, 0,0,1.0, 1,0,2.0, 0,2,3.0, 1,4,4.0);
    ARIADNE_TEST_PRINT(e);
    e.reverse_lexicographic_sort();
    Procedure<Float> p(e);
    ARIADNE_TEST_PRINT(p);
    Vector<Float> x(2, 2.0,3.0);
    ARIADNE_TEST_EQUAL(evaluate(p,x),simple_evaluate(e,x));
}


void TestProcedure::test_evaluate()
{
    Procedure<Float> p;
    p.new_instruction(IND,0u);
    p.new_instruction(IND,1u);
    p.new_instruction(POW,0u,2);
    p.new_instruction(SQR,1u);
    p.new_instruction(ADD,2u,3u);
    p.new_instruction(SQRT,4u);
    ARIADNE_TEST_PRINT(p);

    Vector<Float> x(2, 3.0, 4.0);
    ARIADNE_TEST_PRINT(x);

    ARIADNE_TEST_EQUALS(evaluate(p,x),5.0);
}



int main() {
    TestProcedure().test();
    return ARIADNE_TEST_FAILURES;
}

