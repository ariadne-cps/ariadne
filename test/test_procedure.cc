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

#include "config.h"

#include "function/procedure.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "expression/formula.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestProcedure
{
    //static ApproximateFormula o;;
    //static ApproximateFormula x;
    //static ApproximateFormula y;
  public:
    TestProcedure();
    Void test();
  private:
    Void test_formula();
    Void test_construct_from_formula();
    Void test_construct_from_expansion();
    Void test_evaluate();
    Void test_propagate();
};

//ApproximateFormula TestProcedure::o(ApproximateFormula::constant(1.0));
//ApproximateFormula TestProcedure::x(ApproximateFormula::coordinate(0));
//ApproximateFormula TestProcedure::y(ApproximateFormula::coordinate(1));


TestProcedure::TestProcedure()
{
}

Void TestProcedure::test()
{
    ARIADNE_TEST_CALL(test_formula());
    ARIADNE_TEST_CALL(test_construct_from_formula());
    ARIADNE_TEST_CALL(test_construct_from_expansion());
    ARIADNE_TEST_CALL(test_evaluate());
    ARIADNE_TEST_CALL(test_propagate());
}

Void TestProcedure::test_formula()
{
    ApproximateFormula o(ApproximateFormula::constant(1.0));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ApproximateFormula r(sqrt(pow(x,2)+pow(y,2)));
    ARIADNE_TEST_PRINT(r);

    //Vector< ApproximateFormula > f((sqrt(pow(x,2)+pow(y,2)), atan(y/x)));
}

Void TestProcedure::test_construct_from_formula()
{
    ApproximateFormula o(ApproximateFormula::constant(1.0));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ApproximateFormula xs=pow(x,2);
    ApproximateFormula ys=pow(y,2);

    Vector< ApproximateFormula > f={sqrt(xs+ys), atan(y/x), xs-ys};
    ARIADNE_TEST_PRINT(f);

    ApproximateProcedure p0(f[0]);
    ARIADNE_TEST_PRINT(p0);
    Vector< ApproximateProcedure > p(f);
    ARIADNE_TEST_PRINT(p);

    p0+=ApproximateFloat64(5);
    ARIADNE_TEST_PRINT(p0);
}

Void TestProcedure::test_construct_from_expansion()
{
    {
        Expansion<ApproximateFloat64> e({ {{0,0},1.0}, {{1,0},2.0}, {{0,2},3.0}, {{1,4},4.0} });
        ARIADNE_TEST_PRINT(e);
        e.reverse_lexicographic_sort();
        ARIADNE_TEST_PRINT(e);
        ApproximateProcedure p(e);
        ARIADNE_TEST_PRINT(p);
        Vector<ApproximateFloat64> x={2.0,3.0};
        ARIADNE_TEST_EQUAL(evaluate(p,x),simple_evaluate(e,x));
    }

    {
        Expansion<ApproximateFloat64> e({ {{0,0},1.0}, {{1,0},2.0}, {{0,1},3.0}, {{2,0},4.0}, {{1,1},5.0}, {{0,2},6.0} });
        e.reverse_lexicographic_sort();
        ApproximateProcedure p(e);
        ARIADNE_TEST_PRINT(p);
    }
}


Void TestProcedure::test_evaluate()
{
    ApproximateProcedure p;
    p.new_unary_instruction(OperatorCode::IND,0u);
    p.new_unary_instruction(OperatorCode::IND,1u);
    p.new_power_instruction(OperatorCode::POW,0u,2);
    p.new_unary_instruction(OperatorCode::SQR,1u);
    p.new_binary_instruction(OperatorCode::ADD,2ul,3ul);
    p.new_unary_instruction(OperatorCode::SQRT,4u);
    ARIADNE_TEST_PRINT(p);

    Vector<ApproximateFloat64> x={3.0,4.0};
    ARIADNE_TEST_PRINT(x);

    ARIADNE_TEST_EQUALS(evaluate(p,x),5.0);
}

Void TestProcedure::test_propagate()
{
    {
        ValidatedProcedure p;
        p.new_unary_instruction(OperatorCode::IND,0u);
        p.new_unary_instruction(OperatorCode::IND,1u);
        p.new_unary_instruction(OperatorCode::SQR,0u);
        p.new_unary_instruction(OperatorCode::SQR,1u);
        p.new_binary_instruction(OperatorCode::ADD,2u,3u);
        p.new_unary_instruction(OperatorCode::SQRT,4u);
        ARIADNE_TEST_PRINT(p);

        UpperBoxType x=ExactBoxType{ {0.25,2.0}, {0.5,3.0} };
        ARIADNE_TEST_PRINT(x);

        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_PRINT(x);
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_PRINT(x);
    }

    ValidatedFormula x(ValidatedFormula::coordinate(0));
    ValidatedFormula y(ValidatedFormula::coordinate(1));

    Vector< ValidatedFormula > ff={sqrt(sqr(x)+sqr(y)),2*x-y};
    ARIADNE_TEST_PRINT(ff);
    Vector< ValidatedProcedure > pp(ff);
    ARIADNE_TEST_PRINT(pp);
    UpperBoxType xx=ExactBoxType{ {0.125,2.0}, {0.25,3.0} };
    ExactBoxType cc=ExactBoxType{ {1.0,1.0}, {1.0,1.0} };
    ARIADNE_TEST_PRINT(xx);
    ARIADNE_TEST_PRINT(evaluate(pp,xx));
    simple_hull_reduce(xx,pp,cc);
    ARIADNE_TEST_PRINT(xx);

}



Int main() {
    TestProcedure().test();
    return ARIADNE_TEST_FAILURES;
}

