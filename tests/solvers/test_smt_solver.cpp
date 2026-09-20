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
#include "betterthreads/thread_manager.hpp"

#include "../test.hpp"

using namespace Ariadne;

class ConcurrencyGuard {
  public:
    explicit ConcurrencyGuard(BetterThreads::ThreadManager& thread_manager)
        : _thread_manager(thread_manager), _original(thread_manager.concurrency()) { }

    ~ConcurrencyGuard() { _thread_manager.set_concurrency(_original); }

  private:
    BetterThreads::ThreadManager& _thread_manager;
    SizeType _original;
};

class TestSmtSolver {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_configuration());
        ARIADNE_TEST_CALL(test_result());
        ARIADNE_TEST_CALL(test_solve());
        ARIADNE_TEST_CALL(test_theory_solve());
        ARIADNE_TEST_CALL(test_boolean_theory_solve());
        ARIADNE_TEST_CALL(test_parallel_solve());
    }

  private:
    Void test_configuration() {
        std::cout << "[smt-config] positive epsilon=0.125" << std::endl;
        SmtSolverConfiguration configuration(0.125_x);
        ARIADNE_TEST_EQUAL(configuration.epsilon(),0.125_x);
        ARIADNE_TEST_EQUAL(
            configuration.theory_minimization_budget(),
            std::numeric_limits<SizeType>::max());

        std::cout << "[smt-config] bounded theory minimization budget=1" << std::endl;
        SmtSolverConfiguration bounded_configuration(0.125_x,1u);
        ARIADNE_TEST_EQUAL(bounded_configuration.theory_minimization_budget(),1u);
        ARIADNE_TEST_EQUAL(
            bounded_configuration.learned_clause_limit(),
            std::numeric_limits<SizeType>::max());

        std::cout << "[smt-config] learned clause limit=1" << std::endl;
        SmtSolverConfiguration pruning_configuration(
            0.125_x,std::numeric_limits<SizeType>::max(),1u);
        ARIADNE_TEST_EQUAL(pruning_configuration.learned_clause_limit(),1u);

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
            std::cout << "[smt-solve] forced splitting UNSAT: sin(x)=0 and cos(x)=0 on [0,7]" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,7)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0)),
                ValidatedConstraint(ValidatedNumber(0),cos(x[0]),ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-stats] forced splitting processed="
                      << solve_result.statistics().boxes_processed
                      << " pruned=" << solve_result.statistics().boxes_pruned
                      << " split=" << solve_result.statistics().boxes_split << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_split>0u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed>1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_pruned>0u);
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

    Void test_theory_solve() {
        RealVariable x("x");
        RealExpression ex=x;
        RealSpace space({x});
        SmtSolver solver(SmtSolverConfiguration(0.125_x));

        auto primitive = [&](ContinuousPredicate const& predicate) {
            auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            return alternatives[0][0];
        };

        {
            std::cout << "[smt-theory-solve] EQ at epsilon boundary: x=0 weakened on x=0.125" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex==0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(0.125_x,0.125_x)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
        }

        {
            std::cout << "[smt-theory-solve] EQ outside epsilon: x=0 weakened on x=0.25" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex==0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(0.25_x,0.25_x)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-theory-solve] GEQ at epsilon boundary: x>=0 weakened on x=-0.125" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex>=0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(-0.125_x,-0.125_x)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
        }

        {
            std::cout << "[smt-theory-solve] GT rejects epsilon boundary: x>0 weakened on x=-0.125" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex>0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(-0.125_x,-0.125_x)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-theory-solve] GT accepts strict interior: x>0 weakened on x=-0.0625" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex>0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(-0.0625_x,-0.0625_x)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
        }

        {
            std::cout << "[smt-theory-solve] transcendental GT: sin(x)>0 on [3,4]" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(sin(ex)>0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(3,4)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-theory-solve] parallel/sequential agreement for strict primitive" << std::endl;
            auto& thread_manager=BetterThreads::ThreadManager::instance();
            ConcurrencyGuard concurrency_guard(thread_manager);
            SizeType parallel_concurrency=thread_manager.maximum_concurrency()>=2u ? 2u : thread_manager.maximum_concurrency();
            if(parallel_concurrency>0u) {
                List<SmtTheoryPrimitiveLiteral> literals({primitive(sin(ex)>0)});
                ExactBoxType domain({ExactIntervalType(3,4)});
                thread_manager.set_concurrency(0);
                SmtResult sequential=solver.solve(space,domain,literals);
                thread_manager.set_concurrency(parallel_concurrency);
                SmtResult parallel=solver.solve_parallel(space,domain,literals);
                ARIADNE_TEST_EQUAL(sequential.status(),parallel.status());
                ARIADNE_TEST_ASSERT(sequential.has_witness()==parallel.has_witness());
            }
        }
    }

    Void test_boolean_theory_solve() {
        RealVariable x("x");
        RealExpression ex=x;
        RealSpace space({x});
        SmtSolver solver(SmtSolverConfiguration(0.125_x));

        {
            std::cout << "[smt-dpll] unit propagation on single atom" << std::endl;
            ContinuousPredicate formula=(ex>=0);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            std::cout << "[smt-dpll-stats] decisions="
                      << solve_result.statistics().boolean_decisions
                      << " propagations=" << solve_result.statistics().boolean_propagations
                      << " reasoned=" << solve_result.statistics().boolean_reasoned_propagations
                      << " conflicts=" << solve_result.statistics().boolean_conflicts
                      << " theory_checks=" << solve_result.statistics().theory_checks
                      << " theory_conflicts=" << solve_result.statistics().theory_conflicts << std::endl;
            ARIADNE_TEST_EQUAL(solve_result.statistics().boolean_decisions,0u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_propagations>=1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().boolean_reasoned_propagations,
                solve_result.statistics().boolean_propagations);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_checks,1u);
        }

        {
            std::cout << "[smt-dpll] Boolean contradiction pruned before theory" << std::endl;
            ContinuousPredicate atom=(ex>=0);
            ContinuousPredicate formula=atom&&(!atom);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] contradiction decisions="
                      << solve_result.statistics().boolean_decisions
                      << " propagations=" << solve_result.statistics().boolean_propagations
                      << " conflicts=" << solve_result.statistics().boolean_conflicts
                      << " theory_checks=" << solve_result.statistics().theory_checks << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_conflicts>=1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_checks,0u);
        }

        {
            std::cout << "[smt-dpll] partial theory conflict before free Boolean branch" << std::endl;
            ContinuousPredicate forced_conflict=(ex>=1)&&(ex<=0);
            ContinuousPredicate free_branch=(ex>=-1)||(ex<=2);
            ContinuousPredicate formula=forced_conflict&&free_branch;
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] partial conflict decisions="
                      << solve_result.statistics().boolean_decisions
                      << " propagations=" << solve_result.statistics().boolean_propagations
                      << " theory_checks=" << solve_result.statistics().theory_checks
                      << " theory_conflicts=" << solve_result.statistics().theory_conflicts
                      << " theory_learned=" << solve_result.statistics().theory_learned_clauses
                      << " theory_learned_literals="
                      << solve_result.statistics().theory_learned_clause_literals
                      << " minimization_checks="
                      << solve_result.statistics().theory_minimization_checks
                      << " raw_nogood_literals="
                      << solve_result.statistics().theory_nogood_raw_literals
                      << " minimized_nogood_literals="
                      << solve_result.statistics().theory_nogood_minimized_literals << std::endl;
            ARIADNE_TEST_EQUAL(solve_result.statistics().boolean_decisions,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_checks,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_conflicts,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_learned_clauses,1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_learned_clause_literals>=1u);
        }

        {
            std::cout << "[smt-dpll] minimize theory conflict nogood" << std::endl;
            ContinuousPredicate conflict=(ex>=1)&&(ex<=0);
            ContinuousPredicate irrelevant=(ex>=-100);
            ContinuousPredicate free_branch=(ex>=-2)||(ex<=2);
            ContinuousPredicate formula=conflict&&irrelevant&&free_branch;
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] minimization_checks="
                      << solve_result.statistics().theory_minimization_checks
                      << " raw=" << solve_result.statistics().theory_nogood_raw_literals
                      << " minimized="
                      << solve_result.statistics().theory_nogood_minimized_literals
                      << " removed="
                      << solve_result.statistics().theory_nogood_literals_removed << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_minimization_checks>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().theory_nogood_raw_literals
                > solve_result.statistics().theory_nogood_minimized_literals);
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_nogood_literals_removed>=1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().theory_nogood_raw_literals
                    - solve_result.statistics().theory_nogood_minimized_literals,
                solve_result.statistics().theory_nogood_literals_removed);
        }

        {
            std::cout << "[smt-dpll] minimize theory conflict nogood with multiple irrelevant literals" << std::endl;
            ContinuousPredicate conflict=(ex>=1)&&(ex<=0);
            ContinuousPredicate irrelevant_lower=(ex>=-100);
            ContinuousPredicate irrelevant_upper=(ex<=100);
            ContinuousPredicate free_branch=(ex>=-2)||(ex<=2);
            ContinuousPredicate formula=
                conflict&&irrelevant_lower&&irrelevant_upper&&free_branch;
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] minimization_checks="
                      << solve_result.statistics().theory_minimization_checks
                      << " raw=" << solve_result.statistics().theory_nogood_raw_literals
                      << " minimized="
                      << solve_result.statistics().theory_nogood_minimized_literals
                      << " removed="
                      << solve_result.statistics().theory_nogood_literals_removed << std::endl;
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_nogood_raw_literals,4u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_nogood_minimized_literals,2u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_nogood_literals_removed,2u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_minimization_checks>=4u);
        }

        {
            std::cout << "[smt-dpll] bound theory nogood minimization work" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(0.125_x,1u));
            ContinuousPredicate conflict=(ex>=1)&&(ex<=0);
            ContinuousPredicate irrelevant_lower=(ex>=-100);
            ContinuousPredicate irrelevant_upper=(ex<=100);
            ContinuousPredicate free_branch=(ex>=-2)||(ex<=2);
            ContinuousPredicate formula=
                conflict&&irrelevant_lower&&irrelevant_upper&&free_branch;
            SmtResult solve_result=bounded_solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] budget=1 minimization_checks="
                      << solve_result.statistics().theory_minimization_checks
                      << " raw=" << solve_result.statistics().theory_nogood_raw_literals
                      << " minimized="
                      << solve_result.statistics().theory_nogood_minimized_literals
                      << " removed="
                      << solve_result.statistics().theory_nogood_literals_removed
                      << " budget_exhaustions="
                      << solve_result.statistics().theory_minimization_budget_exhaustions
                      << " first_candidate_trail_rank="
                      << solve_result.statistics().first_minimization_candidate_trail_rank
                      << std::endl;
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_minimization_checks,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().theory_minimization_budget_exhaustions,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_nogood_raw_literals,4u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().first_minimization_candidate_trail_rank>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().theory_nogood_minimized_literals>=2u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().theory_nogood_minimized_literals<=4u);
        }

        {
            std::cout << "[smt-dpll] analyze Boolean conflict with 1-UIP" << std::endl;
            ContinuousPredicate a=(ex>=0);
            ContinuousPredicate b=(ex<=0);
            ContinuousPredicate formula=
                (a||b)&&(a||(!b))&&((!a)||b)&&((!a)||(!b));
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,0)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] analyzed="
                      << solve_result.statistics().boolean_conflicts_analyzed
                      << " learned_clauses=" << solve_result.statistics().learned_clauses
                      << " learned_literals_total=" << solve_result.statistics().learned_clause_literals
                      << " last_learned_literals=" << solve_result.statistics().last_learned_clause_literals
                      << " last_current_level_literals="
                      << solve_result.statistics().last_learned_current_level_literals
                      << " learned_propagations="
                      << solve_result.statistics().learned_clause_propagations
                      << " backjump_level=" << solve_result.statistics().last_backjump_level
                      << " nonchronological_backjumps="
                      << solve_result.statistics().nonchronological_backjumps
                      << " max_level=" << solve_result.statistics().max_decision_level << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_conflicts_analyzed>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clauses>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().last_learned_clause_literals>=1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().last_learned_current_level_literals,1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clause_propagations>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().last_backjump_level
                < solve_result.statistics().max_decision_level);
        }

        {
            std::cout << "[smt-dpll] learned clause propagation after backjump" << std::endl;
            ContinuousPredicate a=(ex>=-2);
            ContinuousPredicate b=(ex>=-1);
            ContinuousPredicate c=(ex>=0);
            ContinuousPredicate formula=
                (a||b)&&(a||(!b))&&((!a)||c)&&((!a)||(!c));
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,0)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] learned="
                      << solve_result.statistics().learned_clauses
                      << " learned_propagations="
                      << solve_result.statistics().learned_clause_propagations
                      << " backtracks=" << solve_result.statistics().boolean_backtracks
                      << " max_level=" << solve_result.statistics().max_decision_level << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clauses>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clause_propagations>=1u);
        }

        {
            std::cout << "[smt-dpll] prune inactive low-activity learned clauses" << std::endl;
            SmtSolver pruning_solver(SmtSolverConfiguration(
                0.125_x,std::numeric_limits<SizeType>::max(),1u));
            ContinuousPredicate a=(ex>=0);
            ContinuousPredicate b=(ex<=0);
            ContinuousPredicate c=(ex==0);
            ContinuousPredicate formula=
                ( a|| b|| c)&&
                ( a|| b||(!c))&&
                ( a||(!b)|| c)&&
                ( a||(!b)||(!c))&&
                ((!a)|| b|| c)&&
                ((!a)|| b||(!c))&&
                ((!a)||(!b)|| c)&&
                ((!a)||(!b)||(!c));
            SmtResult solve_result=pruning_solver.solve(
                space,ExactBoxType({ExactIntervalType(-1,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] learned="
                      << solve_result.statistics().learned_clauses
                      << " activity_bumps="
                      << solve_result.statistics().learned_clause_activity_bumps
                      << " pruning_runs="
                      << solve_result.statistics().learned_clause_pruning_runs
                      << " pruned=" << solve_result.statistics().learned_clauses_pruned
                      << " peak_active="
                      << solve_result.statistics().peak_active_non_theory_learned_clauses
                      << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clauses>=2u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clause_activity_bumps>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clause_pruning_runs>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clauses_pruned>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().peak_active_non_theory_learned_clauses>=2u);
        }

        {
            std::cout << "[smt-dpll] conjunction UNSAT: x>=1 and x<=0 on [0,1]" << std::endl;
            ContinuousPredicate formula=(ex>=1)&&(ex<=0);
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-dpll] disjunction EPSILON_SAT: x<0 or x>1 on [0,0.0625]" << std::endl;
            ContinuousPredicate formula=(ex<0)||(ex>1);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0.0_x,0.0625_x)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-dpll] theory conflict clause reused by Boolean propagation" << std::endl;
            ContinuousPredicate left=(ex<0);
            ContinuousPredicate right=(ex>1);
            ContinuousPredicate formula=left||right;
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0.5_x,0.5_x)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] theory_learned="
                      << solve_result.statistics().theory_learned_clauses
                      << " theory_learned_literals="
                      << solve_result.statistics().theory_learned_clause_literals
                      << " theory_learned_propagations="
                      << solve_result.statistics().theory_learned_clause_propagations
                      << " learned_propagations="
                      << solve_result.statistics().learned_clause_propagations << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_learned_clauses>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().theory_learned_clause_propagations>=1u);
        }

        {
            std::cout << "[smt-dpll] trail backtracking across Boolean branches" << std::endl;
            ContinuousPredicate formula=(ex<0)||(ex>1);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0.5_x,0.5_x)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            std::cout << "[smt-dpll-stats] trail decisions="
                      << solve_result.statistics().boolean_decisions
                      << " propagations=" << solve_result.statistics().boolean_propagations
                      << " reasoned=" << solve_result.statistics().boolean_reasoned_propagations
                      << " backtracks=" << solve_result.statistics().boolean_backtracks
                      << " max_level=" << solve_result.statistics().max_decision_level
                      << " theory_checks=" << solve_result.statistics().theory_checks
                      << " theory_conflicts=" << solve_result.statistics().theory_conflicts << std::endl;
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_decisions>=1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().boolean_reasoned_propagations,
                solve_result.statistics().boolean_propagations);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_backtracks>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().max_decision_level>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_conflicts>=1u);
        }

        {
            std::cout << "[smt-dpll] disjunction UNSAT: x<0 or x>1 on [0.375,0.625]" << std::endl;
            ContinuousPredicate formula=(ex<0)||(ex>1);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0.375_x,0.625_x)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-dpll] negated atom preserves strict boundary: !(x<=0) at x=-epsilon" << std::endl;
            ContinuousPredicate formula=!(ex<=0);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(-0.125_x,-0.125_x)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        {
            std::cout << "[smt-dpll] disequality expands to theory alternatives: x!=0 at x=0" << std::endl;
            ContinuousPredicate formula=(ex!=0);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,0)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
        }

        {
            std::cout << "[smt-dpll] mixed transcendental Boolean formula" << std::endl;
            ContinuousPredicate formula=((sin(ex)==0)||(cos(ex)==0))&&(ex>=3);
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(3,4)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-dpll] sequential/parallel theory agreement" << std::endl;
            auto& thread_manager=BetterThreads::ThreadManager::instance();
            ConcurrencyGuard concurrency_guard(thread_manager);
            SizeType parallel_concurrency=thread_manager.maximum_concurrency()>=2u
                ? 2u : thread_manager.maximum_concurrency();
            if(parallel_concurrency>0u) {
                ContinuousPredicate formula=((sin(ex)==0)||(cos(ex)==0))&&(ex>=3);
                ExactBoxType domain({ExactIntervalType(3,4)});

                thread_manager.set_concurrency(0);
                SmtResult sequential=solver.solve(space,domain,formula);

                thread_manager.set_concurrency(parallel_concurrency);
                SmtResult parallel=solver.solve_parallel(space,domain,formula);

                ARIADNE_TEST_EQUAL(sequential.status(),parallel.status());
                ARIADNE_TEST_ASSERT(sequential.has_witness()==parallel.has_witness());
            }
        }
    }

    Void test_parallel_solve() {
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        SmtSolver solver(SmtSolverConfiguration(0.125_x));
        auto& thread_manager=BetterThreads::ThreadManager::instance();
        ConcurrencyGuard concurrency_guard(thread_manager);

        std::cout << "[smt-parallel] sequential BetterThreads mode concurrency=0" << std::endl;
        thread_manager.set_concurrency(0);
        {
            ExactBoxType domain({ExactIntervalType(3,4)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve_parallel(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            std::cout << "[smt-parallel] concurrency=0 processed="
                      << solve_result.statistics().boxes_processed
                      << " pruned=" << solve_result.statistics().boxes_pruned
                      << " split=" << solve_result.statistics().boxes_split << std::endl;
        }

        SizeType parallel_concurrency=thread_manager.maximum_concurrency()>=2u ? 2u : thread_manager.maximum_concurrency();
        if(parallel_concurrency>0u) {
            std::cout << "[smt-parallel] concurrent BetterThreads mode concurrency="
                      << parallel_concurrency << std::endl;
            thread_manager.set_concurrency(parallel_concurrency);

            {
                ExactBoxType domain({ExactIntervalType(3,4)});
                List<ValidatedConstraint> constraints({
                    ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0))
                });
                SmtResult solve_result=solver.solve_parallel(domain,constraints);
                ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
                ARIADNE_TEST_ASSERT(solve_result.has_witness());
                std::cout << "[smt-parallel] EPSILON_SAT processed="
                          << solve_result.statistics().boxes_processed
                          << " pruned=" << solve_result.statistics().boxes_pruned
                          << " split=" << solve_result.statistics().boxes_split << std::endl;
            }

            {
                std::cout << "[smt-parallel] concurrent forced splitting UNSAT: sin(x)=0 and cos(x)=0 on [0,7]" << std::endl;
                ExactBoxType domain({ExactIntervalType(0,7)});
                List<ValidatedConstraint> constraints({
                    ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0)),
                    ValidatedConstraint(ValidatedNumber(0),cos(x[0]),ValidatedNumber(0))
                });
                SmtResult solve_result=solver.solve_parallel(domain,constraints);
                ARIADNE_TEST_ASSERT(solve_result.is_unsat());
                std::cout << "[smt-parallel] forced splitting processed="
                          << solve_result.statistics().boxes_processed
                          << " pruned=" << solve_result.statistics().boxes_pruned
                          << " split=" << solve_result.statistics().boxes_split << std::endl;
                ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_split>0u);
                ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed>1u);
                ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_pruned>0u);
            }

            {
                std::cout << "[smt-parallel] concurrent UNSAT: sin(x)=2 on [3,4]" << std::endl;
                ExactBoxType domain({ExactIntervalType(3,4)});
                List<ValidatedConstraint> constraints({
                    ValidatedConstraint(ValidatedNumber(2),sin(x[0]),ValidatedNumber(2))
                });
                SmtResult solve_result=solver.solve_parallel(domain,constraints);
                ARIADNE_TEST_ASSERT(solve_result.is_unsat());
                ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_pruned>=1u);
            }
        } else {
            std::cout << "[smt-parallel] hardware reports no worker concurrency; concurrent case skipped" << std::endl;
        }

        std::cout << "[smt-parallel] compare sequential and parallel logical status" << std::endl;
        auto compare_status = [&](String const& label,
                                  ExactBoxType const& domain,
                                  List<ValidatedConstraint> const& constraints) {
            thread_manager.set_concurrency(0);
            SmtResult sequential_result=solver.solve(domain,constraints);

            if(parallel_concurrency>0u) {
                thread_manager.set_concurrency(parallel_concurrency);
                SmtResult parallel_result=solver.solve_parallel(domain,constraints);

                std::cout << "[smt-parallel-compare] " << label
                          << " sequential=" << sequential_result.status()
                          << " parallel=" << parallel_result.status() << std::endl;

                ARIADNE_TEST_EQUAL(sequential_result.status(),parallel_result.status());
                ARIADNE_TEST_ASSERT(sequential_result.has_witness()==parallel_result.has_witness());
            }
        };

        compare_status(
            "linear epsilon-sat",
            ExactBoxType({ExactIntervalType(0,1)}),
            List<ValidatedConstraint>({
                ValidatedConstraint(ValidatedNumber(1),2*x[0],ValidatedNumber(1))
            }));

        compare_status(
            "linear unsat",
            ExactBoxType({ExactIntervalType(0,1)}),
            List<ValidatedConstraint>({
                ValidatedConstraint(ValidatedNumber(2),x[0],ValidatedNumber(2))
            }));

        compare_status(
            "transcendental epsilon-sat",
            ExactBoxType({ExactIntervalType(3,4)}),
            List<ValidatedConstraint>({
                ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0))
            }));

        compare_status(
            "transcendental unsat",
            ExactBoxType({ExactIntervalType(3,4)}),
            List<ValidatedConstraint>({
                ValidatedConstraint(ValidatedNumber(2),sin(x[0]),ValidatedNumber(2))
            }));

        compare_status(
            "forced splitting transcendental unsat",
            ExactBoxType({ExactIntervalType(0,7)}),
            List<ValidatedConstraint>({
                ValidatedConstraint(ValidatedNumber(0),sin(x[0]),ValidatedNumber(0)),
                ValidatedConstraint(ValidatedNumber(0),cos(x[0]),ValidatedNumber(0))
            }));
    }
};

Int main() {
    TestSmtSolver().test();
    return ARIADNE_TEST_FAILURES;
}
