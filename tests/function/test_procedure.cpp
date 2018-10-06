/***************************************************************************
 *            test_procedure.cpp
 *
 *  Copyright  2010  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>

#include "config.hpp"

#include "function/procedure.hpp"
#include "function/procedure.tpl.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/evaluate.hpp"
#include "function/formula.hpp"

#include "../test.hpp"

#include "algebra/covector.hpp"
#include "function/function.hpp"
#include "algebra/expansion.inl.hpp"

using namespace std;
using namespace Ariadne;

template<class X> decltype(auto) mag(Covector<X> const& u) { return norm(transpose(u)); }

class TestProcedure
{
    //static ApproximateFormula o;
    //static ApproximateFormula x;
    //static ApproximateFormula y;
    DoublePrecision pr;
  public:
    TestProcedure();
    Void test();
  private:
    Void test_formula();
    Void test_construct_from_formula();
    Void test_construct_from_expansion();
    Void test_evaluate();
    Void test_propagate();
    Void test_derivative();
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
    ARIADNE_TEST_CALL(test_derivative());


}


Void TestProcedure::test_formula()
{
    ApproximateFormula o(ApproximateFormula::constant(1.0));
    ApproximateFormula t(ApproximateFormula::constant(2.0));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ApproximateFormula r(sqrt(pow(x,2)+pow(y,2))/t);
    ARIADNE_TEST_PRINT(r);

    //Vector<ApproximateFormula> f((sqrt(pow(x,2)+pow(y,2)), atan(y/x)));
}

Void TestProcedure::test_construct_from_formula()
{
    ApproximateNumber c(2.0);
    ApproximateFormula o(ApproximateFormula::constant(1.0));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ApproximateFormula xs=pow(x,2);
    ApproximateFormula ys=pow(y,2);

    ARIADNE_TEST_PRINT(y*c);
    ARIADNE_TEST_PRINT(c*y);
    Vector<ApproximateFormula> f={sqrt(xs+c*ys), atan(y/x), xs*c-ys};
    ARIADNE_TEST_PRINT(f);

    ApproximateProcedure p0(f[0]);
    ARIADNE_TEST_PRINT(p0);
    Vector<ApproximateProcedure> p(f);
    ARIADNE_TEST_PRINT(p);
}


Void TestProcedure::test_construct_from_expansion()
{
    {
        Expansion<MultiIndex,FloatDPApproximation> e({ {{0,0},1.0}, {{1,0},2.0}, {{0,2},3.0}, {{1,4},4.0} },pr);
        ARIADNE_TEST_PRINT(e);
        e.reverse_lexicographic_sort();
        ARIADNE_TEST_PRINT(e);
        Procedure<ApproximateNumber> p(e);
        ARIADNE_TEST_PRINT(p);
        Vector<FloatDPApproximation> x({2.0,3.0},pr);
        ARIADNE_TEST_EQUAL(evaluate(p,x),simple_evaluate(e,x));
    }

    {
        Expansion<MultiIndex,FloatDPApproximation> e({ {{0,0},1.0}, {{1,0},2.0}, {{0,1},3.0}, {{2,0},4.0}, {{1,1},5.0}, {{0,2},6.0} },pr);
        e.reverse_lexicographic_sort();
        Procedure<ApproximateNumber> p(e);
        ARIADNE_TEST_PRINT(p);
        Vector<FloatDPApproximation> x({2.0,3.0},pr);
        ARIADNE_TEST_EQUAL(evaluate(p,x),simple_evaluate(e,x));
    }
}


Void TestProcedure::test_evaluate()
{
    ApproximateProcedure p;
    p.new_unary_instruction(OperatorCode::IND,0ul);
    p.new_unary_instruction(OperatorCode::IND,1ul);
    p.new_graded_instruction(OperatorCode::POW,0ul,2);
    p.new_unary_instruction(OperatorCode::SQR,1ul);
    p.new_binary_instruction(OperatorCode::ADD,2ul,3ul);
    p.new_constant(9.0);
    p.new_scalar_instruction(OperatorCode::SMUL,0ul,4ul);
    p.new_unary_instruction(OperatorCode::SQRT,5ul);
    ARIADNE_TEST_PRINT(p);

    Vector<FloatDPApproximation> x({3.0,4.0},pr);
    ARIADNE_TEST_PRINT(x);

    ARIADNE_TEST_EQUALS(evaluate(p,x),15.0);
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
        p.new_constant(1.125_decimal);
        p.new_scalar_instruction(OperatorCode::SMUL,0ul,4ul);
        p.new_unary_instruction(OperatorCode::SQRT,5u);
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

    Vector<ValidatedFormula> ff={sqrt(sqr(x)+sqr(y)),2*x-y};
    ARIADNE_TEST_PRINT(ff);
    Vector<ValidatedProcedure> pp(ff);
    ARIADNE_TEST_PRINT(pp);
    UpperBoxType xx=ExactBoxType{ {0.125,2.0}, {0.25,3.0} };
    ExactBoxType cc=ExactBoxType{ {1.0,1.0}, {1.0,1.0} };
    ARIADNE_TEST_PRINT(xx);
    ARIADNE_TEST_PRINT(evaluate(pp,xx));
    simple_hull_reduce(xx,pp,cc);
    ARIADNE_TEST_PRINT(xx);
}

Void TestProcedure::test_derivative()
{
    typedef ApproximateTag P;
    typedef Number<P> Y;
    typedef FloatDPApproximation X;

    Y c(2);
    Formula<Y> o(Formula<Y>::constant(1));
    Formula<Y> x(Formula<Y>::coordinate(0));
    Formula<Y> y(Formula<Y>::coordinate(1));

    auto xs=sqr(x); auto ys=sqr(y);
    Formula<Y> e=sqrt(xs+c*ys)+sin(y/x)*xs*c-ys;
    ARIADNE_TEST_PRINT(e);
    ScalarFunction<P> f(EuclideanDomain(2),e);
    ARIADNE_TEST_PRINT(f);
    Procedure<Y> p(e);
    ARIADNE_TEST_PRINT(p);

    X zero(0,dp);

    Vector<X> q({2,1},dp);
    Vector<X> s({-1,3},dp);
    Vector<Differential<X>> ds=Differential<X>::variable(1,2,zero,0)*s+q;
    ARIADNE_TEST_WITHIN(gradient(p,q),f.gradient(q),8e-16);
    ARIADNE_TEST_EQUALS(hessian(p,q,s),f(ds).hessian().get(0,0));
}




Int main() {
    TestProcedure().test();
    return ARIADNE_TEST_FAILURES;
}

