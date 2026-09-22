/***************************************************************************
 *            test_smt_theory.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 */

#include <sstream>

#include "numeric/numeric.hpp"
#include "solvers/smt_theory.hpp"
#include "symbolic/variable.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestSmtTheory {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_relations());
        ARIADNE_TEST_CALL(test_negations());
        ARIADNE_TEST_CALL(test_normalization());
        ARIADNE_TEST_CALL(test_weakening());
        ARIADNE_TEST_CALL(test_reject_non_atom());
    }

  private:
    Void check_relation(ContinuousPredicate const& predicate,
                        SmtTheoryRelation expected,
                        String const& label) {
        std::cout << "[smt-theory] relation " << label << std::endl;
        SmtTheoryLiteral literal=make_smt_theory_literal(predicate);
        ARIADNE_TEST_EQUAL(literal.relation(),expected);
    }

    Void test_relations() {
        RealVariable x("x"), y("y");
        RealExpression ex=x;
        RealExpression ey=y;

        check_relation(ex==ey,SmtTheoryRelation::EQ,"x==y");
        check_relation(ex!=ey,SmtTheoryRelation::NEQ,"x!=y");
        check_relation(ex<=ey,SmtTheoryRelation::LEQ,"x<=y");
        check_relation(ex>=ey,SmtTheoryRelation::GEQ,"x>=y");
        check_relation(ex<ey,SmtTheoryRelation::LT,"x<y");
        check_relation(ex>ey,SmtTheoryRelation::GT,"x>y");
        check_relation(sgn(ex),SmtTheoryRelation::GT,"sgn(x) -> x>0");
    }

    Void check_negation(ContinuousPredicate const& predicate,
                        SmtTheoryRelation expected,
                        String const& label) {
        std::cout << "[smt-theory] negation " << label << std::endl;
        SmtTheoryLiteral literal=make_smt_theory_literal(predicate).negated();
        ARIADNE_TEST_EQUAL(literal.relation(),expected);
    }

    Void test_negations() {
        RealVariable x("x"), y("y");
        RealExpression ex=x;
        RealExpression ey=y;

        check_negation(ex==ey,SmtTheoryRelation::NEQ,"!(x==y) -> x!=y");
        check_negation(ex!=ey,SmtTheoryRelation::EQ,"!(x!=y) -> x==y");
        check_negation(ex<=ey,SmtTheoryRelation::GT,"!(x<=y) -> x>y");
        check_negation(ex>=ey,SmtTheoryRelation::LT,"!(x>=y) -> x<y");
        check_negation(ex<ey,SmtTheoryRelation::GEQ,"!(x<y) -> x>=y");
        check_negation(ex>ey,SmtTheoryRelation::LEQ,"!(x>y) -> x<=y");
        check_negation(sgn(ex),SmtTheoryRelation::LEQ,"!sgn(x) -> x<=0");
    }

    Void test_normalization() {
        RealVariable x("x"), y("y");
        RealExpression ex=x;
        RealExpression ey=y;

        auto check_single = [&](ContinuousPredicate const& predicate,
                                SmtTheoryPrimitiveRelation expected_relation,
                                RealExpression const& expected_expression,
                                String const& label) {
            std::cout << "[smt-theory] normalize " << label << std::endl;
            auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0][0].relation(),expected_relation);
            ARIADNE_TEST_SAME(alternatives[0][0].expression(),expected_expression);
        };

        check_single(ex==ey,SmtTheoryPrimitiveRelation::EQ_ZERO,ex-ey,"x==y -> x-y=0");
        check_single(ex>=ey,SmtTheoryPrimitiveRelation::GEQ_ZERO,ex-ey,"x>=y -> x-y>=0");
        check_single(ex>ey,SmtTheoryPrimitiveRelation::GT_ZERO,ex-ey,"x>y -> x-y>0");
        check_single(ex<=ey,SmtTheoryPrimitiveRelation::GEQ_ZERO,ey-ex,"x<=y -> y-x>=0");
        check_single(ex<ey,SmtTheoryPrimitiveRelation::GT_ZERO,ey-ex,"x<y -> y-x>0");
        check_single(
            sgn(ex),
            SmtTheoryPrimitiveRelation::GT_ZERO,
            ex-RealExpression(0),
            "sgn(x) -> x>0");

        std::cout << "[smt-theory] normalize x!=y -> (x-y>0) or (y-x>0)" << std::endl;
        auto neq=normalize_smt_theory_literal(make_smt_theory_literal(ex!=ey));
        ARIADNE_TEST_EQUAL(neq.size(),2u);
        ARIADNE_TEST_EQUAL(neq[0].size(),1u);
        ARIADNE_TEST_EQUAL(neq[1].size(),1u);
        ARIADNE_TEST_EQUAL(neq[0][0].relation(),SmtTheoryPrimitiveRelation::GT_ZERO);
        ARIADNE_TEST_EQUAL(neq[1][0].relation(),SmtTheoryPrimitiveRelation::GT_ZERO);
        ARIADNE_TEST_SAME(neq[0][0].expression(),ex-ey);
        ARIADNE_TEST_SAME(neq[1][0].expression(),ey-ex);
    }

    Void test_weakening() {
        RealVariable x("x"), y("y");
        RealExpression ex=x;
        RealExpression ey=y;
        ExactDouble epsilon=0.125_x;

        auto check = [&](ContinuousPredicate const& predicate,
                         SmtTheoryWeakRelation expected,
                         String const& label) {
            std::cout << "[smt-theory] weaken " << label << " epsilon=" << epsilon << std::endl;
            auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            auto weakened=weaken_smt_theory_literal(alternatives[0][0],epsilon);
            ARIADNE_TEST_EQUAL(weakened.relation(),expected);
            ARIADNE_TEST_EQUAL(weakened.epsilon(),epsilon);
        };

        check(ex==ey,SmtTheoryWeakRelation::ABS_LEQ_EPSILON,"x==y -> |x-y|<=epsilon");
        check(ex>=ey,SmtTheoryWeakRelation::GEQ_MINUS_EPSILON,"x>=y -> x-y>=-epsilon");
        check(ex>ey,SmtTheoryWeakRelation::GT_MINUS_EPSILON,"x>y -> x-y>-epsilon");
        check(sgn(ex),SmtTheoryWeakRelation::GT_MINUS_EPSILON,"sgn(x) -> x>-epsilon");

        std::cout << "[smt-theory] weakening requires positive epsilon" << std::endl;
        auto primitive=normalize_smt_theory_literal(make_smt_theory_literal(ex==ey))[0][0];
        ARIADNE_TEST_THROWS(weaken_smt_theory_literal(primitive,0.0_x),std::runtime_error);
    }

    Void test_reject_non_atom() {
        RealVariable x("x"), y("y");
        RealExpression ex=x;
        RealExpression ey=y;
        ContinuousPredicate formula=(ex<=0)&&(ey>=1);

        std::cout << "[smt-theory] reject non-atomic conjunction" << std::endl;
        ARIADNE_TEST_THROWS(make_smt_theory_literal(formula),std::runtime_error);
    }
};

Int main() {
    TestSmtTheory().test();
    return ARIADNE_TEST_FAILURES;
}
