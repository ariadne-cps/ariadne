/***************************************************************************
 *            test_smt_boolean.cpp
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
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/numeric.hpp"
#include "solvers/smt_boolean.hpp"
#include "symbolic/variable.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestSmtBoolean {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_atom());
        ARIADNE_TEST_CALL(test_constants());
        ARIADNE_TEST_CALL(test_sign_atom());
        ARIADNE_TEST_CALL(test_negation());
        ARIADNE_TEST_CALL(test_conjunction());
        ARIADNE_TEST_CALL(test_disjunction());
        ARIADNE_TEST_CALL(test_shared_atom());
        ARIADNE_TEST_CALL(test_structurally_shared_atom());
        ARIADNE_TEST_CALL(test_semantically_shared_atom());
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

    Void test_constants() {
        {
            std::cout << "[smt-boolean] constant true" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(ContinuousPredicate(true));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),0u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.clauses().size(),2u);
            ARIADNE_TEST_EQUAL(encoding.clauses()[0].size(),1u);
            ARIADNE_TEST_EQUAL(encoding.clauses()[0][0],1);
            ARIADNE_TEST_EQUAL(encoding.clauses()[1][0],1);
        }

        {
            std::cout << "[smt-boolean] constant false" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(ContinuousPredicate(false));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),0u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.clauses().size(),2u);
            ARIADNE_TEST_EQUAL(encoding.clauses()[0].size(),1u);
            ARIADNE_TEST_EQUAL(encoding.clauses()[0][0],-1);
            ARIADNE_TEST_EQUAL(encoding.clauses()[1][0],1);
        }

        {
            RealVariable x("x");
            ContinuousPredicate atom=(x<=0);
            std::cout << "[smt-boolean] constants composed with atom" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(
                (ContinuousPredicate(true)&&atom)||ContinuousPredicate(false));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.clauses().size(),1u);
            ARIADNE_TEST_EQUAL(encoding.clauses()[0][0],1);
        }

        {
            RealVariable x("x");
            ContinuousPredicate atom=(x<=0);

            std::cout << "[smt-boolean] fold false&&atom" << std::endl;
            SmtBooleanEncoding and_false=SmtBooleanEncoder().encode(
                ContinuousPredicate(false)&&atom);
            ARIADNE_TEST_EQUAL(and_false.atom_count(),0u);
            ARIADNE_TEST_EQUAL(and_false.variable_count(),1u);
            ARIADNE_TEST_EQUAL(and_false.clauses().size(),2u);
            ARIADNE_TEST_EQUAL(and_false.clauses()[0][0],-1);
            ARIADNE_TEST_EQUAL(and_false.clauses()[1][0],1);

            std::cout << "[smt-boolean] fold true||atom" << std::endl;
            SmtBooleanEncoding or_true=SmtBooleanEncoder().encode(
                ContinuousPredicate(true)||atom);
            ARIADNE_TEST_EQUAL(or_true.atom_count(),0u);
            ARIADNE_TEST_EQUAL(or_true.variable_count(),1u);
            ARIADNE_TEST_EQUAL(or_true.clauses().size(),2u);
            ARIADNE_TEST_EQUAL(or_true.clauses()[0][0],1);
            ARIADNE_TEST_EQUAL(or_true.clauses()[1][0],1);

            std::cout << "[smt-boolean] fold negated constant" << std::endl;
            SmtBooleanEncoding negated=SmtBooleanEncoder().encode(
                !ContinuousPredicate(false));
            ARIADNE_TEST_EQUAL(negated.atom_count(),0u);
            ARIADNE_TEST_EQUAL(negated.variable_count(),1u);
            ARIADNE_TEST_EQUAL(negated.clauses().size(),2u);
            ARIADNE_TEST_EQUAL(negated.clauses()[0][0],1);
            ARIADNE_TEST_EQUAL(negated.clauses()[1][0],1);
        }
    }

    Void test_sign_atom() {
        RealVariable x("x");
        RealExpression ex=x;
        ContinuousPredicate predicate=sgn(ex);

        std::cout << "[smt-boolean] sign atom sgn(x)" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(predicate);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses()[0].size(),1u);
        ARIADNE_TEST_EQUAL(encoding.clauses()[0][0],1);
        ARIADNE_TEST_EQUAL(encoding.atom_variable(0),1u);
        ARIADNE_TEST_EQUAL(encoding.atom(0).code(),OperatorCode::SGN);
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

    Void test_structurally_shared_atom() {
        RealVariable x("x");
        RealExpression ex=x;
        ContinuousPredicate first=(ex<=0);
        ContinuousPredicate second=(ex<=0);

        ARIADNE_TEST_ASSERT(first.node_raw_ptr()!=second.node_raw_ptr());
        std::cout << "[smt-boolean] structurally shared atoms built as distinct nodes" << std::endl;
        SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(first&&!second);

        ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
        ARIADNE_TEST_EQUAL(encoding.variable_count(),2u);
        ARIADNE_TEST_EQUAL(encoding.clauses().size(),4u);
        ARIADNE_TEST_EQUAL(encoding.atom_variable(0),1u);
    }

    Void test_semantically_shared_atom() {
        RealVariable x("x"), y("y");
        RealExpression ex=x;
        RealExpression ey=y;

        {
            std::cout << "[smt-boolean] x<=y and y>=x share an atom" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode((ex<=ey)&&!(ey>=ex));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),2u);
        }

        {
            std::cout << "[smt-boolean] x<y and y>x share an atom" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode((ex<ey)&&!(ey>ex));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),2u);
        }

        {
            std::cout << "[smt-boolean] symmetric equality shares an atom" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode((ex==ey)&&!(ey==ex));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),2u);
        }

        {
            std::cout << "[smt-boolean] sgn(x) and x>0 share an atom" << std::endl;
            SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(sgn(ex)&&!(ex>0));
            ARIADNE_TEST_EQUAL(encoding.atom_count(),1u);
            ARIADNE_TEST_EQUAL(encoding.variable_count(),2u);
        }
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
