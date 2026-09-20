/***************************************************************************
 *            test_smt_solver.cpp
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

#include <sstream>

#include "solvers/smt_solver.hpp"

#include "../test.hpp"

using namespace Ariadne;

class TestSmtSolver {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_configuration());
        ARIADNE_TEST_CALL(test_result());
    }

  private:
    Void test_configuration() {
        SmtSolverConfiguration configuration(0.125_x);
        ARIADNE_TEST_EQUAL(configuration.epsilon(),0.125_x);

        ARIADNE_TEST_THROWS(SmtSolverConfiguration(0.0_x),std::runtime_error);
        ARIADNE_TEST_THROWS(SmtSolverConfiguration(-0.125_x),std::runtime_error);
    }

    Void test_result() {
        SmtResult unsat=SmtResult::unsat();
        ARIADNE_TEST_ASSERT(unsat.is_unsat());
        ARIADNE_TEST_ASSERT(not unsat.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(not unsat.has_witness());
        ARIADNE_TEST_THROWS(unsat.witness(),std::runtime_error);

        UpperBoxType witness=ExactBoxType({ExactIntervalType(-1,1),ExactIntervalType(2,3)});
        SmtResult epsilon_sat=SmtResult::epsilon_sat(witness);
        ARIADNE_TEST_ASSERT(not epsilon_sat.is_unsat());
        ARIADNE_TEST_ASSERT(epsilon_sat.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(epsilon_sat.has_witness());
        ARIADNE_TEST_EQUAL(epsilon_sat.witness(),witness);

        std::ostringstream oss;
        oss << unsat.status() << " " << epsilon_sat.status();
        ARIADNE_TEST_EQUAL(oss.str(),String("UNSAT EPSILON_SAT"));
    }
};

Int main() {
    TestSmtSolver().test();
    return ARIADNE_TEST_FAILURES;
}
