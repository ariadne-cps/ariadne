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

        std::cout << "[smt-result] default search statistics are zero" << std::endl;
        ARIADNE_TEST_EQUAL(unsat.statistics().boxes_processed,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().boxes_pruned,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().boxes_split,0u);
        ARIADNE_TEST_EQUAL(epsilon_sat.statistics().boxes_processed,0u);
        ARIADNE_TEST_EQUAL(epsilon_sat.statistics().boxes_pruned,0u);
        ARIADNE_TEST_EQUAL(epsilon_sat.statistics().boxes_split,0u);

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
            std::cout << "[smt-stats] linear UNSAT processed="
                      << solve_result.statistics().boxes_processed
                      << " pruned=" << solve_result.statistics().boxes_pruned
                      << " split=" << solve_result.statistics().boxes_split << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_pruned>=1u);
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

        {
            std::cout << "[smt-solve] conjunction EPSILON_SAT: x in [0,1], x>=0.25 and x<=0.75" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0.25_x),x[0],ValidatedNumber(+infty)),
                ValidatedConstraint(ValidatedNumber(-infty),x[0],ValidatedNumber(0.75_x))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] conjunction UNSAT: x in [0,1], x<=0.25 and x>=0.75" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(-infty),x[0],ValidatedNumber(0.25_x)),
                ValidatedConstraint(ValidatedNumber(0.75_x),x[0],ValidatedNumber(+infty))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-solve] inequality EPSILON_SAT: x in [0,1], x<=0.25" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(-infty),x[0],ValidatedNumber(0.25_x))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] inequality EPSILON_SAT: x in [0,1], x>=0.75" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0.75_x),x[0],ValidatedNumber(+infty))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] epsilon boundary EPSILON_SAT: x in [0,1], x=1.125 with epsilon=0.125" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(1.125_x),x[0],ValidatedNumber(1.125_x))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] two-dimensional EPSILON_SAT: x+y=1, x-y=0" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            ExactBoxType domain({ExactIntervalType(0,1),ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(1),xy[0]+xy[1],ValidatedNumber(1)),
                ValidatedConstraint(ValidatedNumber(0),xy[0]-xy[1],ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.witness().dimension(),2u);
        }

        {
            std::cout << "[smt-solve] two-dimensional UNSAT: x+y=3 on [0,1]^2" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            ExactBoxType domain({ExactIntervalType(0,1),ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(3),xy[0]+xy[1],ValidatedNumber(3))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-solve] multiple transcendental branches EPSILON_SAT: sin(x)=0 on [0,7]" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,7)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            std::cout << "[smt-stats] multiple branches processed="
                      << solve_result.statistics().boxes_processed
                      << " pruned=" << solve_result.statistics().boxes_pruned
                      << " split=" << solve_result.statistics().boxes_split << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed>=1u);
        }

        {
            std::cout << "[smt-solve] empty constraint list: bounded nonempty domain is EPSILON_SAT" << std::endl;
            ExactBoxType domain({ExactIntervalType(-1,1)});
            List<ValidatedConstraint> constraints;
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-solve] empty domain is UNSAT" << std::endl;
            ExactBoxType domain({ExactIntervalType::empty_interval()});
            List<ValidatedConstraint> constraints;
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-solve] reject constraint dimension mismatch" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),xy[0],ValidatedNumber(1))
            });
            ARIADNE_TEST_THROWS(solver.solve(domain,constraints),std::runtime_error);
        }

        {
            std::cout << "[smt-solve] reject unbounded domain" << std::endl;
            ExactBoxType domain({ExactIntervalType(-infty,+infty)});
            List<ValidatedConstraint> constraints;
            ARIADNE_TEST_THROWS(solver.solve(domain,constraints),std::runtime_error);
        }
    }
};

Int main() {
    TestSmtSolver().test();
    return ARIADNE_TEST_FAILURES;
}
