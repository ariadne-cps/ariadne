/***************************************************************************
 *            test_smt_theory.cpp
 *
 *  Copyright  2026  Ariadne developers
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

        check_relation(x==y,SmtTheoryRelation::EQ,"x==y");
        check_relation(x!=y,SmtTheoryRelation::NEQ,"x!=y");
        check_relation(x<=y,SmtTheoryRelation::LEQ,"x<=y");
        check_relation(x>=y,SmtTheoryRelation::GEQ,"x>=y");
        check_relation(x<y,SmtTheoryRelation::LT,"x<y");
        check_relation(x>y,SmtTheoryRelation::GT,"x>y");
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

        check_negation(x==y,SmtTheoryRelation::NEQ,"!(x==y) -> x!=y");
        check_negation(x!=y,SmtTheoryRelation::EQ,"!(x!=y) -> x==y");
        check_negation(x<=y,SmtTheoryRelation::GT,"!(x<=y) -> x>y");
        check_negation(x>=y,SmtTheoryRelation::LT,"!(x>=y) -> x<y");
        check_negation(x<y,SmtTheoryRelation::GEQ,"!(x<y) -> x>=y");
        check_negation(x>y,SmtTheoryRelation::LEQ,"!(x>y) -> x<=y");
    }

    Void test_reject_non_atom() {
        RealVariable x("x"), y("y");
        ContinuousPredicate formula=(x<=0)&&(y>=1);

        std::cout << "[smt-theory] reject non-atomic conjunction" << std::endl;
        ARIADNE_TEST_THROWS(make_smt_theory_literal(formula),std::runtime_error);
    }
};

Int main() {
    TestSmtTheory().test();
    return ARIADNE_TEST_FAILURES;
}
