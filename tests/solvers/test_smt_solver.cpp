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
        ARIADNE_TEST_CALL(test_solve());
    }

  private:
    Void test_configuration() {
        std::cout << "[smt-config] positive epsilon=0.125" << std::endl;
        SmtSolverConfiguration configuration(0.125_x);
        ARIADNE_TEST_EQUAL(configuration.epsilon(),0.125_x);

        std::cout << "[smt-config] reject zero epsilon" << std::endl;
        ARIADNE_TEST_THROWS(SmtSolverConfiguration(0.0_x),std::runtime_error);

        std::cout << "[smt-config] reject negative epsilon=-0.125" << std::endl;
        ARIADNE_TEST_THROWS(SmtSolverConfiguration(-0.125_x),std::runtime_error);
    }

    Void test_result() {
        std::cout << "[smt-result] construct UNSAT result" << std::endl;
        SmtResult unsat=SmtResult::unsat();
        ARIADNE_TEST_ASSERT(unsat.is_unsat());
        ARIADNE_TEST_ASSERT(not unsat.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(not unsat.has_witness());

        std::cout << "[smt-result] reject witness access for UNSAT" << std::endl;
        ARIADNE_TEST_THROWS(unsat.witness(),std::runtime_error);

        UpperBoxType witness=ExactBoxType({ExactIntervalType(-1,1),ExactIntervalType(2,3)});
        std::cout << "[smt-result] construct EPSILON_SAT result witness=" << witness << std::endl;
        SmtResult epsilon_sat=SmtResult::epsilon_sat(witness);
        ARIADNE_TEST_ASSERT(not epsilon_sat.is_unsat());
        ARIADNE_TEST_ASSERT(epsilon_sat.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(epsilon_sat.has_witness());

        std::cout << "[smt-result] verify witness dimension=" << witness.dimension() << std::endl;
        ARIADNE_TEST_EQUAL(epsilon_sat.witness().dimension(),witness.dimension());
        for(SizeType i=0; i!=witness.dimension(); ++i) {
            std::cout << "[smt-result] verify witness coordinate=" << i
                      << " expected=" << witness[i]
                      << " actual=" << epsilon_sat.witness()[i] << std::endl;
            ARIADNE_TEST_EQUAL(epsilon_sat.witness()[i].lower_bound().raw(),witness[i].lower_bound().raw());
            ARIADNE_TEST_EQUAL(epsilon_sat.witness()[i].upper_bound().raw(),witness[i].upper_bound().raw());
        }

        std::cout << "[smt-result] stream statuses" << std::endl;
        std::ostringstream oss;
        oss << unsat.status() << " " << epsilon_sat.status();
        ARIADNE_TEST_EQUAL(oss.str(),String("UNSAT EPSILON_SAT"));
    }

    Void test_solve() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        SmtSolver solver(SmtSolverConfiguration(0.125_x));

        {
            std::cout << "[smt-solve] linear UNSAT: x in [0,1], x=2" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(2),x[0],ValidatedNumber(2))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-solve] linear EPSILON_SAT: x in [0,1], 2*x=1" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(1),2*x[0],ValidatedNumber(1))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] transcendental EPSILON_SAT: x in [3,4], sin(x)=0" << std::endl;
            ExactBoxType domain({ExactIntervalType(3,4)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] transcendental UNSAT: x in [3,4], sin(x)=2" << std::endl;
            ExactBoxType domain({ExactIntervalType(3,4)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(2),sin(x[0]),ValidatedNumber(2))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-solve] even power EPSILON_SAT: x in [-2,0], x^2=1" << std::endl;
            ExactBoxType domain({ExactIntervalType(-2,0)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(1),sqr(x[0]),ValidatedNumber(1))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }
    }
};

Int main() {
    TestSmtSolver().test();
    return ARIADNE_TEST_FAILURES;
}
