/***************************************************************************
 *            test_smt_boolean.cpp
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
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "solvers/smt_boolean.hpp"
#include "symbolic/variable.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestSmtBoolean {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_atom());
        ARIADNE_TEST_CALL(test_negation());
        ARIADNE_TEST_CALL(test_conjunction());
        ARIADNE_TEST_CALL(test_disjunction());
        ARIADNE_TEST_CALL(test_shared_atom());
        ARIADNE_TEST_CALL(test_nested_formula());
    }

  private:
    Void test_atom() {
        RealVariable x("x");
        ContinuousPredicate predicate=(x<=0);

        std::cout << "[smt-boolean] atom x<=0" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(predicate);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses()[0].size(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses()[0][0],1);
        ARIADNE_TEST_EQUAL(encoding.atom_variable(0),1u);
        ARIADNE_TEST_EQUAL(encoding.atom(0).code(),OperatorCode::LEQ);
    }

    Void test_negation() {
        RealVariable x("x");
        ContinuousPredicate atom=(x<=0);

        std::cout << "[smt-boolean] negation !(x<=0)" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(!atom);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses()[0][0],-1);
    }

    Void test_conjunction() {
        RealVariable x("x"), y("y");
        ContinuousPredicate lhs=(x<=0);
        ContinuousPredicate rhs=(y>=1);

        std::cout << "[smt-boolean] conjunction (x<=0)&&(y>=1)" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(lhs&&rhs);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),2u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),3u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),4u);

        std::cout << "[smt-boolean] conjunction variables=" << encoding.variable_count()
                  << " clauses=" << encoding.clauses().size() << std::endl;
    }

    Void test_disjunction() {
        RealVariable x("x"), y("y");
        ContinuousPredicate lhs=(x<=0);
        ContinuousPredicate rhs=(y>=1);

        std::cout << "[smt-boolean] disjunction (x<=0)||(y>=1)" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(lhs||rhs);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),2u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),3u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),4u);
    }

    Void test_shared_atom() {
        RealVariable x("x");
        ContinuousPredicate atom=(x<=0);

        std::cout << "[smt-boolean] shared atom (x<=0)||!(x<=0)" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(atom||!atom);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),2u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),4u);
    }

    Void test_nested_formula() {
        RealVariable x("x"), y("y"), z("z");
        ContinuousPredicate a=(x<=0);
        ContinuousPredicate b=(y>=1);
        ContinuousPredicate c=(z==2);

        std::cout << "[smt-boolean] nested ((x<=0)||(y>=1))&&!(z==2)" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode((a||b)&&!c);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),3u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),5u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),7u);

        std::cout << "[smt-boolean] nested atoms=" << encoding.atom_count()
                  << " variables=" << encoding.variable_count()
                  << " clauses=" << encoding.clauses().size() << std::endl;
    }
};

Int main() {
    TestSmtBoolean().test();
    return ARIADNE_TEST_FAILURES;
}
