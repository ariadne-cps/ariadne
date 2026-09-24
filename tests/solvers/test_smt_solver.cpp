/***************************************************************************
 *            test_smt_solver.cpp
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
        ARIADNE_TEST_CALL(test_public_preconditions());
        ARIADNE_TEST_CALL(test_result());
        ARIADNE_TEST_CALL(test_learned_clause_pruning_policy());
        ARIADNE_TEST_CALL(test_statistics_aggregation());
        ARIADNE_TEST_CALL(test_parallel_state_transitions());
        ARIADNE_TEST_CALL(test_epsilon_witness_candidate_limits());
        ARIADNE_TEST_CALL(test_candidate_witness_outcome());
        ARIADNE_TEST_CALL(test_sensitivity_split_selection());
        ARIADNE_TEST_CALL(test_monotone_coordinate_gating());
        ARIADNE_TEST_CALL(test_search_outcome());
        ARIADNE_TEST_CALL(test_cdcl_helpers());
        ARIADNE_TEST_CALL(test_invalid_internal_relations());
        ARIADNE_TEST_CALL(test_epsilon_predicates());
        ARIADNE_TEST_CALL(test_box_processing_statistics());
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

        std::cout << "[smt-config] box processing limit=0" << std::endl;
        SmtSolverConfiguration bounded_search_configuration(
            0.125_x,
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            0u);
        ARIADNE_TEST_EQUAL(bounded_search_configuration.box_processing_limit(),0u);
        ARIADNE_TEST_ASSERT(bounded_search_configuration.candidate_search_enabled());

        std::cout << "[smt-config] disable candidate witness search" << std::endl;
        SmtSolverConfiguration no_candidate_configuration(
            0.125_x,
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            false);
        ARIADNE_TEST_ASSERT(not no_candidate_configuration.candidate_search_enabled());
        ARIADNE_TEST_ASSERT(not configuration.monotone_reduction_enabled());
        std::cout << "[smt-config] enable monotone reduction" << std::endl;
        SmtSolverConfiguration monotone_configuration(
            0.125_x,
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            std::numeric_limits<SizeType>::max(),
            true,
            true);
        ARIADNE_TEST_ASSERT(monotone_configuration.monotone_reduction_enabled());

        std::cout << "[smt-config] reject zero epsilon" << std::endl;
        ARIADNE_TEST_THROWS(SmtSolverConfiguration(0.0_x),std::runtime_error);

        std::cout << "[smt-config] reject negative epsilon=-0.125" << std::endl;
        ARIADNE_TEST_THROWS(SmtSolverConfiguration(-0.125_x),std::runtime_error);
    }

    Void test_public_preconditions() {
        std::cout << "[smt-preconditions] reject invalid public inputs" << std::endl;

        SmtSolver solver(SmtSolverConfiguration(0.125_x));
        ExactBoxType bounded({ExactIntervalType(0,1)});
        ExactBoxType unbounded({ExactIntervalType(-inf,inf)});

        auto x1=ValidatedScalarMultivariateFunction::coordinates(1);
        auto x2=ValidatedScalarMultivariateFunction::coordinates(2);
        List<ValidatedConstraint> one_dimensional({
            ValidatedConstraint(
                ValidatedNumber(0),x1[0],ValidatedNumber(1))
        });
        List<ValidatedConstraint> two_dimensional({
            ValidatedConstraint(
                ValidatedNumber(0),x2[0],ValidatedNumber(1))
        });

        ARIADNE_TEST_THROWS(
            solver.solve(unbounded,one_dimensional),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve(bounded,two_dimensional),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve_parallel(unbounded,one_dimensional),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve_parallel(bounded,two_dimensional),
            std::runtime_error);

        RealVariable x("precondition_x"), y("precondition_y");
        RealExpression ex=x;
        RealSpace one_space({x});
        RealSpace two_space({x,y});
        List<SmtTheoryPrimitiveLiteral> literals({
            normalize_smt_theory_literal(
                make_smt_theory_literal(ex>=0))[0][0]
        });

        ARIADNE_TEST_THROWS(
            solver.solve(one_space,unbounded,literals),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve(two_space,bounded,literals),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve_parallel(one_space,unbounded,literals),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve_parallel(two_space,bounded,literals),
            std::runtime_error);

        ContinuousPredicate predicate=(ex>=0);
        ARIADNE_TEST_THROWS(
            solver.solve(one_space,unbounded,predicate),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve(two_space,bounded,predicate),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve_parallel(one_space,unbounded,predicate),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            solver.solve_parallel(two_space,bounded,predicate),
            std::runtime_error);
    }

    Void test_result() {
        std::cout << "[smt-result] construct UNSAT result" << std::endl;
        SmtResult unsat=SmtResult::unsat();
        ARIADNE_TEST_ASSERT(unsat.is_unsat());
        ARIADNE_TEST_ASSERT(not unsat.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(not unsat.has_witness());
        ARIADNE_TEST_ASSERT(not unsat.is_unknown());

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
        ARIADNE_TEST_EQUAL(unsat.statistics().box_budget_exhaustions,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().non_splittable_uncertified_boxes,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().non_splittable_epsilon_overlap_boxes,0u);
        ARIADNE_TEST_EQUAL(epsilon_sat.statistics().boxes_processed,0u);
        ARIADNE_TEST_EQUAL(epsilon_sat.statistics().boxes_pruned,0u);
        ARIADNE_TEST_EQUAL(epsilon_sat.statistics().boxes_split,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().hull_reduction_rounds,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().hull_effective_reductions,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().shaving_reduction_rounds,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().shaving_effective_reductions,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().sensitivity_guided_splits,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().sensitivity_overrides_geometric_splits,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().epsilon_box_certifications,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().candidate_witness_searches,0u);
        ARIADNE_TEST_EQUAL(unsat.statistics().candidate_witness_successes,0u);

        std::cout << "[smt-result] construct UNKNOWN result" << std::endl;
        SmtResult unknown=SmtResult::unknown();
        ARIADNE_TEST_ASSERT(unknown.is_unknown());
        ARIADNE_TEST_ASSERT(not unknown.is_unsat());
        ARIADNE_TEST_ASSERT(not unknown.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(not unknown.has_witness());

        std::cout << "[smt-result] stream statuses" << std::endl;
        std::ostringstream oss;
        oss << unsat.status() << " " << epsilon_sat.status() << " " << unknown.status();
        ARIADNE_TEST_EQUAL(oss.str(),String("UNSAT EPSILON_SAT UNKNOWN"));

        std::cout << "[smt-result] reject invalid status" << std::endl;
        auto invalid_status=static_cast<SmtResultStatus>(999);
        ARIADNE_TEST_THROWS(
            static_cast<Void>(
                static_cast<OutputStream&>(std::cout) << invalid_status),
            std::runtime_error);
    }

    Void test_learned_clause_pruning_policy() {
        std::cout << "[smt-dpll] deterministic learned clause pruning policy" << std::endl;
        using Entry=SmtSolverTestSupport::LearnedClausePruningEntry;
        std::vector<Entry> entries({
            {false,false,false,false,false,false,false,1u,5u},
            {true,true,false,false,false,false,false,1u,5u},
            {true,false,true,false,false,false,false,1u,5u},
            {true,false,false,true,false,false,false,1u,2u},
            {true,false,false,false,true,false,false,2u,5u},
            {true,false,false,false,false,true,false,1u,5u},
            {true,false,false,false,false,false,true,1u,5u},
            {true,false,false,false,false,false,false,2u,7u},
            {true,false,false,false,false,false,false,1u,4u},
            {true,false,false,false,false,false,false,1u,6u}
        });
        std::vector<SizeType> candidates=
            SmtSolverTestSupport::learned_clause_pruning_candidates(entries);
        ARIADNE_TEST_EQUAL(candidates.size(),3u);
        ARIADNE_TEST_EQUAL(candidates[0],9u);
        ARIADNE_TEST_EQUAL(candidates[1],8u);
        ARIADNE_TEST_EQUAL(candidates[2],7u);

        SizeType pruned=SmtSolverTestSupport::apply_learned_clause_pruning(
            candidates,3u,1u);
        ARIADNE_TEST_EQUAL(pruned,0u);

        SizeType none=SmtSolverTestSupport::apply_learned_clause_pruning(
            candidates,1u,1u);
        ARIADNE_TEST_EQUAL(none,0u);
    }

    Void test_statistics_aggregation() {
        std::cout << "[smt-stats] deterministic aggregation semantics" << std::endl;

        SmtSearchStatistics target;
        target.last_learned_clause_literals=11u;
        target.last_learned_current_level_literals=7u;
        target.last_backjump_level=5u;

        SmtSearchStatistics no_conflict;
        no_conflict.boxes_processed=3u;
        no_conflict.monotone_reduction_rounds=4u;
        no_conflict.monotone_effective_reductions=2u;
        no_conflict.last_learned_clause_literals=99u;
        no_conflict.last_learned_current_level_literals=99u;
        no_conflict.last_backjump_level=99u;
        SmtSolverTestSupport::accumulate_statistics(target,no_conflict);
        ARIADNE_TEST_EQUAL(target.boxes_processed,3u);
        ARIADNE_TEST_EQUAL(target.monotone_reduction_rounds,4u);
        ARIADNE_TEST_EQUAL(target.monotone_effective_reductions,2u);
        ARIADNE_TEST_EQUAL(target.last_learned_clause_literals,11u);
        ARIADNE_TEST_EQUAL(target.last_learned_current_level_literals,7u);
        ARIADNE_TEST_EQUAL(target.last_backjump_level,5u);

        SmtSearchStatistics empty_target;
        SmtSearchStatistics empty_source;
        SmtSolverTestSupport::accumulate_statistics(empty_target,empty_source);
        ARIADNE_TEST_EQUAL(empty_target.first_minimization_candidate_trail_rank,0u);

        SmtSearchStatistics first_rank;
        SmtSolverTestSupport::record_first_minimization_candidate_trail_rank(
            first_rank,7u);
        SmtSolverTestSupport::record_first_minimization_candidate_trail_rank(
            first_rank,11u);
        ARIADNE_TEST_EQUAL(
            first_rank.first_minimization_candidate_trail_rank,7u);

        SmtSearchStatistics conflict;
        conflict.boolean_conflicts_analyzed=1u;
        conflict.last_learned_clause_literals=4u;
        conflict.last_learned_current_level_literals=1u;
        conflict.last_backjump_level=2u;
        conflict.first_minimization_candidate_trail_rank=9u;
        SmtSolverTestSupport::accumulate_statistics(target,conflict);
        ARIADNE_TEST_EQUAL(target.last_learned_clause_literals,4u);
        ARIADNE_TEST_EQUAL(target.last_learned_current_level_literals,1u);
        ARIADNE_TEST_EQUAL(target.last_backjump_level,2u);
        ARIADNE_TEST_EQUAL(target.first_minimization_candidate_trail_rank,9u);

        SmtSearchStatistics later;
        later.first_minimization_candidate_trail_rank=13u;
        SmtSolverTestSupport::accumulate_statistics(target,later);
        ARIADNE_TEST_EQUAL(target.first_minimization_candidate_trail_rank,9u);
    }

    Void test_parallel_state_transitions() {
        std::cout << "[smt-parallel] deterministic shared-state transitions" << std::endl;
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::parallel_stop_condition(false,false));
        ARIADNE_TEST_ASSERT(SmtSolverTestSupport::parallel_stop_condition(true,false));
        ARIADNE_TEST_ASSERT(SmtSolverTestSupport::parallel_stop_condition(false,true));
        ARIADNE_TEST_ASSERT(SmtSolverTestSupport::parallel_stop_condition(true,true));

        UpperBoxType left({UpperIntervalType(ExactIntervalType(-1,0))});
        UpperBoxType right({UpperIntervalType(ExactIntervalType(0,1))});
        Pair<UpperBoxType,UpperBoxType> children(left,right);
        auto append_children=
            SmtSolverTestSupport::parallel_children_to_append(false,children);
        ARIADNE_TEST_EQUAL(append_children.size(),2u);
        auto stopped_children=
            SmtSolverTestSupport::parallel_children_to_append(true,children);
        ARIADNE_TEST_ASSERT(stopped_children.empty());

        auto claims=SmtSolverTestSupport::parallel_witness_claim_sequence();
        ARIADNE_TEST_ASSERT(claims.first);
        ARIADNE_TEST_ASSERT(not claims.second);
    }

    Void test_epsilon_witness_candidate_limits() {
        std::cout << "[smt-candidate] deterministic corner candidate limits" << std::endl;

        UpperIntervalType unit(ExactIntervalType(0,1));
        UpperBoxType six_dimensions(6u,unit);
        UpperBoxType seven_dimensions(7u,unit);
        UpperBoxType sixty_four_dimensions(64u,unit);

        ARIADNE_TEST_EQUAL(
            SmtSolverTestSupport::epsilon_witness_candidate_count(six_dimensions),
            79u);
        ARIADNE_TEST_EQUAL(
            SmtSolverTestSupport::epsilon_witness_candidate_count(seven_dimensions),
            17u);
        ARIADNE_TEST_EQUAL(
            SmtSolverTestSupport::epsilon_witness_candidate_count(sixty_four_dimensions),
            131u);
    }

    Void test_candidate_witness_outcome() {
        std::cout << "[smt-candidate] deterministic candidate outcome" << std::endl;
        UpperBoxType witness({
            UpperIntervalType(ExactIntervalType(1,1))
        });

        auto disabled=SmtSolverTestSupport::candidate_witness_outcome(
            false,std::nullopt);
        ARIADNE_TEST_ASSERT(not disabled.attempted);
        ARIADNE_TEST_ASSERT(not disabled.certified);
        ARIADNE_TEST_ASSERT(not disabled.witness.has_value());

        auto failed=SmtSolverTestSupport::candidate_witness_outcome(
            true,std::nullopt);
        ARIADNE_TEST_ASSERT(failed.attempted);
        ARIADNE_TEST_ASSERT(not failed.certified);
        ARIADNE_TEST_ASSERT(not failed.witness.has_value());

        auto success=SmtSolverTestSupport::candidate_witness_outcome(
            true,witness);
        ARIADNE_TEST_ASSERT(success.attempted);
        ARIADNE_TEST_ASSERT(success.certified);
        ARIADNE_TEST_ASSERT(success.witness.has_value());
        ARIADNE_TEST_EQUAL(success.witness->dimension(),1u);
    }

    Void test_sensitivity_split_selection() {
        std::cout << "[smt-split] deterministic sensitivity selection" << std::endl;
        auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
        UpperBoxType domain({
            UpperIntervalType(ExactIntervalType(0,1)),
            UpperIntervalType(ExactIntervalType(0,4))
        });

        {
            std::vector<ValidatedScalarMultivariateFunction> functions({
                sqr(xy[0])
            });
            auto selection=SmtSolverTestSupport::sensitivity_split_selection(
                domain,functions);
            ARIADNE_TEST_ASSERT(selection.guided);
            ARIADNE_TEST_EQUAL(selection.coordinate,0u);
            ARIADNE_TEST_ASSERT(selection.overrode_geometric);
        }

        {
            std::vector<ValidatedScalarMultivariateFunction> functions({
                ValidatedScalarMultivariateFunction::zero(2u)
            });
            auto selection=SmtSolverTestSupport::sensitivity_split_selection(
                domain,functions);
            ARIADNE_TEST_ASSERT(not selection.guided);
            ARIADNE_TEST_EQUAL(selection.coordinate,1u);
            ARIADNE_TEST_ASSERT(not selection.overrode_geometric);
        }

        {
            std::vector<ValidatedScalarMultivariateFunction> functions({
                xy[0]+2*xy[1]
            });
            auto selection=SmtSolverTestSupport::sensitivity_split_selection(
                domain,functions);
            ARIADNE_TEST_ASSERT(selection.guided);
            ARIADNE_TEST_EQUAL(selection.coordinate,1u);
            ARIADNE_TEST_ASSERT(not selection.overrode_geometric);
        }

        {
            UpperBoxType wide_first({
                UpperIntervalType(ExactIntervalType(0,4)),
                UpperIntervalType(ExactIntervalType(0,1))
            });
            std::vector<ValidatedScalarMultivariateFunction> functions({
                2*xy[0]+xy[1]
            });
            auto selection=SmtSolverTestSupport::sensitivity_split_selection(
                wide_first,functions);
            ARIADNE_TEST_ASSERT(selection.guided);
            ARIADNE_TEST_EQUAL(selection.coordinate,0u);
            ARIADNE_TEST_ASSERT(not selection.overrode_geometric);
        }
    }

    Void test_monotone_coordinate_gating() {
        std::cout << "[smt-icp] validated monotone coordinate gating" << std::endl;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        UpperBoxType positive_domain({UpperIntervalType(ExactIntervalType(0,2))});
        UpperBoxType monotone_sine_domain({UpperIntervalType(ExactIntervalType(3,4))});
        UpperBoxType nonmonotone_sine_domain({UpperIntervalType(ExactIntervalType(0,4))});

        ARIADNE_TEST_ASSERT(
            SmtSolverTestSupport::monotone_coordinate_is_safe(
                exp(x[0])+x[0],positive_domain,0u));
        ARIADNE_TEST_ASSERT(
            SmtSolverTestSupport::monotone_coordinate_is_safe(
                -exp(x[0])-x[0],positive_domain,0u));
        ARIADNE_TEST_ASSERT(
            SmtSolverTestSupport::monotone_coordinate_is_safe(
                sin(x[0]),monotone_sine_domain,0u));
        ARIADNE_TEST_ASSERT(
            not SmtSolverTestSupport::monotone_coordinate_is_safe(
                sin(x[0]),nonmonotone_sine_domain,0u));
    }

    Void test_search_outcome() {
        std::cout << "[smt-dpll] deterministic search outcome classification" << std::endl;
        using Outcome=SmtSolverTestSupport::SearchOutcome;

        Outcome exhausted=Outcome::exhausted();
        ARIADNE_TEST_ASSERT(not exhausted.witness.has_value());
        ARIADNE_TEST_ASSERT(not exhausted.backjump_level.has_value());

        Outcome backjump=Outcome::backjump(3u);
        ARIADNE_TEST_ASSERT(not backjump.witness.has_value());
        ARIADNE_TEST_ASSERT(backjump.backjump_level.has_value());
        ARIADNE_TEST_EQUAL(*backjump.backjump_level,3u);

        UpperBoxType witness({
            UpperIntervalType(ExactIntervalType(0,0))
        });
        Outcome found=Outcome::found(witness);
        ARIADNE_TEST_ASSERT(found.witness.has_value());
        ARIADNE_TEST_ASSERT(not found.backjump_level.has_value());

        SmtSearchStatistics statistics;
        SmtResult sat=SmtSolverTestSupport::finalize_search_outcome(
            found,false,statistics);
        ARIADNE_TEST_ASSERT(sat.is_epsilon_sat());
        ARIADNE_TEST_ASSERT(sat.has_witness());

        SmtResult unknown=SmtSolverTestSupport::finalize_search_outcome(
            exhausted,true,statistics);
        ARIADNE_TEST_ASSERT(unknown.is_unknown());

        SmtResult unsat=SmtSolverTestSupport::finalize_search_outcome(
            exhausted,false,statistics);
        ARIADNE_TEST_ASSERT(unsat.is_unsat());
    }

    Void test_cdcl_helpers() {
        std::cout << "[smt-cdcl] deterministic clause resolution and nogood ordering" << std::endl;

        std::vector<Int> resolved=SmtSolverTestSupport::resolve_clause_on_variable(
            std::vector<Int>({1,2,3,3}),
            std::vector<Int>({-1,3,4,4}),
            1u);
        ARIADNE_TEST_EQUAL(resolved.size(),3u);
        ARIADNE_TEST_EQUAL(resolved[0],2);
        ARIADNE_TEST_EQUAL(resolved[1],3);
        ARIADNE_TEST_EQUAL(resolved[2],4);

        std::vector<Int> clause({1,-2,3,-4});
        std::vector<SizeType> decision_levels({0u,1u,2u,2u,1u});
        std::vector<SizeType> trail_rank({0u,1u,2u,4u,3u});
        SmtSolverTestSupport::order_theory_nogood(
            clause,decision_levels,trail_rank);
        ARIADNE_TEST_EQUAL(clause[0],3);
        ARIADNE_TEST_EQUAL(clause[1],-2);
        ARIADNE_TEST_EQUAL(clause[2],-4);
        ARIADNE_TEST_EQUAL(clause[3],1);

        {
            std::cout << "[smt-theory-propagation] validated atom classification" << std::endl;
            RealVariable tx("classify_x");
            RealExpression etx=tx;
            RealSpace tspace({tx});
            ExactBoxType zero({ExactIntervalType(0,0)});
            ExactBoxType positive({ExactIntervalType(1,2)});
            ExactBoxType negative({ExactIntervalType(-2,-1)});
            ExactBoxType crossing({ExactIntervalType(-1,1)});
            using Truth=SmtSolverTestSupport::TheoryAtomTruth;

            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,zero,(etx==0))==Truth::TRUE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,positive,(etx==0))==Truth::FALSE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,crossing,(etx==0))==Truth::UNKNOWN);

            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,positive,(etx!=0))==Truth::TRUE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,zero,(etx!=0))==Truth::FALSE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,crossing,(etx!=0))==Truth::UNKNOWN);

            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,positive,(etx>=0))==Truth::TRUE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,negative,(etx>=0))==Truth::FALSE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,crossing,(etx>=0))==Truth::UNKNOWN);

            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,negative,(etx<=0))==Truth::TRUE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,positive,(etx<=0))==Truth::FALSE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,crossing,(etx<=0))==Truth::UNKNOWN);

            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,positive,(etx>0))==Truth::TRUE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,zero,(etx>0))==Truth::FALSE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,crossing,(etx>0))==Truth::UNKNOWN);

            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,negative,(etx<0))==Truth::TRUE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,zero,(etx<0))==Truth::FALSE_VALUE);
            ARIADNE_TEST_ASSERT(
                SmtSolverTestSupport::classify_theory_atom(
                    tspace,crossing,(etx<0))==Truth::UNKNOWN);
        }

        {
            std::cout << "[smt-theory-propagation] epsilon-safe domain implications" << std::endl;
            RealVariable ix("implication_x");
            RealExpression eix=ix;
            RealSpace ispace({ix});
            ExactBoxType domain({ExactIntervalType(0,1)});
            SmtSolver implication_solver(SmtSolverConfiguration(0.125_x));

            auto false_implication=SmtSolverTestSupport::domain_theory_implication(
                implication_solver,ispace,domain,(eix>2));
            ARIADNE_TEST_ASSERT(not false_implication.force_true);
            ARIADNE_TEST_ASSERT(false_implication.force_false);

            auto true_implication=SmtSolverTestSupport::domain_theory_implication(
                implication_solver,ispace,domain,(eix>=-2));
            ARIADNE_TEST_ASSERT(true_implication.force_true);
            ARIADNE_TEST_ASSERT(not true_implication.force_false);

            auto relaxed_strict=SmtSolverTestSupport::domain_theory_implication(
                implication_solver,
                ispace,
                ExactBoxType({ExactIntervalType(0,0)}),
                (eix>0));
            ARIADNE_TEST_ASSERT(not relaxed_strict.force_true);
            ARIADNE_TEST_ASSERT(not relaxed_strict.force_false);

            auto undecided=SmtSolverTestSupport::domain_theory_implication(
                implication_solver,ispace,domain,(eix==0.5_x));
            ARIADNE_TEST_ASSERT(not undecided.force_true);
            ARIADNE_TEST_ASSERT(not undecided.force_false);
        }

        std::vector<Bool> theory_flags({false,true});
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::clause_is_learned(1u,2u));
        ARIADNE_TEST_ASSERT(SmtSolverTestSupport::clause_is_learned(2u,2u));
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::learned_clause_is_theory(
            1u,2u,theory_flags));
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::learned_clause_is_theory(
            2u,2u,theory_flags));
        ARIADNE_TEST_ASSERT(SmtSolverTestSupport::learned_clause_is_theory(
            3u,2u,theory_flags));
        ARIADNE_TEST_ASSERT(SmtSolverTestSupport::assignment_locks_clause(
            1,std::optional<SizeType>(7u),7u));
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::assignment_locks_clause(
            -1,std::optional<SizeType>(7u),7u));
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::assignment_locks_clause(
            1,std::nullopt,7u));
        ARIADNE_TEST_ASSERT(not SmtSolverTestSupport::assignment_locks_clause(
            1,std::optional<SizeType>(6u),7u));

        UpperBoxType witness({
            UpperIntervalType(ExactIntervalType(0,0))
        });
        auto sat=SmtSolverTestSupport::interpret_theory_result(
            SmtResult::epsilon_sat(witness));
        ARIADNE_TEST_ASSERT(sat.consistent);
        ARIADNE_TEST_ASSERT(not sat.unknown);
        ARIADNE_TEST_ASSERT(sat.witness.has_value());

        auto unknown=SmtSolverTestSupport::interpret_theory_result(
            SmtResult::unknown());
        ARIADNE_TEST_ASSERT(unknown.consistent);
        ARIADNE_TEST_ASSERT(unknown.unknown);
        ARIADNE_TEST_ASSERT(not unknown.witness.has_value());

        auto unsat=SmtSolverTestSupport::interpret_theory_result(
            SmtResult::unsat());
        ARIADNE_TEST_ASSERT(not unsat.consistent);
        ARIADNE_TEST_ASSERT(not unsat.unknown);
        ARIADNE_TEST_ASSERT(not unsat.witness.has_value());
    }

    Void test_invalid_internal_relations() {
        std::cout << "[smt-internal] reject invalid primitive relation" << std::endl;
        SmtSolver solver(SmtSolverConfiguration(0.125_x));
        auto invalid=static_cast<SmtTheoryPrimitiveRelation>(999);
        ARIADNE_TEST_THROWS(
            SmtSolverTestSupport::original_bounds(solver,invalid),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            SmtSolverTestSupport::epsilon_bounds(solver,invalid),
            std::runtime_error);
        ARIADNE_TEST_THROWS(
            SmtSolverTestSupport::epsilon_primitive_image_infeasible(
                invalid,
                UpperIntervalType(ExactIntervalType(0,0)),
                FloatDP(0.125_x,dp)),
            std::runtime_error);
        SmtTheoryRelation invalid_theory_relation=
            static_cast<SmtTheoryRelation>(999);
        ARIADNE_TEST_THROWS(
            SmtSolverTestSupport::classify_theory_relation(
                invalid_theory_relation,
                UpperIntervalType(ExactIntervalType(0,0))),
            std::runtime_error);
    }

    Void test_epsilon_predicates() {
        std::cout << "[smt-epsilon] deterministic satisfaction" << std::endl;
        SmtSolver solver(SmtSolverConfiguration(0.125_x));
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        UpperBoxType point({
            UpperIntervalType(ExactIntervalType(0,0))
        });

        List<ValidatedConstraint> outside({
            ValidatedConstraint(
                ValidatedNumber(1),
                x[0],
                ValidatedNumber(1))
        });
        ARIADNE_TEST_ASSERT(
            not SmtSolverTestSupport::epsilon_satisfied(
                solver,point,outside));

        List<ValidatedConstraint> inside({
            ValidatedConstraint(
                ValidatedNumber(0),
                x[0],
                ValidatedNumber(0))
        });
        ARIADNE_TEST_ASSERT(
            SmtSolverTestSupport::epsilon_satisfied(
                solver,point,inside));

    }

    Void test_box_processing_statistics() {
        std::cout << "[smt-stats] box processing status accounting" << std::endl;
        using Status=SmtSolverTestSupport::BoxProcessingStatus;
        using Input=SmtSolverTestSupport::BoxProcessingStatisticsInput;

        SmtSearchStatistics statistics;
        SmtSolverTestSupport::accumulate_box_processing_statistics(
            statistics,Input{Status::PRUNED,1u,1u,2u,1u,3u,1u,true,true,true,true,true});
        ARIADNE_TEST_EQUAL(statistics.boxes_pruned,1u);
        ARIADNE_TEST_EQUAL(statistics.hull_reduction_rounds,1u);
        ARIADNE_TEST_EQUAL(statistics.hull_effective_reductions,1u);
        ARIADNE_TEST_EQUAL(statistics.shaving_reduction_rounds,2u);
        ARIADNE_TEST_EQUAL(statistics.shaving_effective_reductions,1u);
        ARIADNE_TEST_EQUAL(statistics.monotone_reduction_rounds,3u);
        ARIADNE_TEST_EQUAL(statistics.monotone_effective_reductions,1u);
        ARIADNE_TEST_EQUAL(statistics.sensitivity_guided_splits,1u);
        ARIADNE_TEST_EQUAL(statistics.sensitivity_overrides_geometric_splits,1u);
        ARIADNE_TEST_EQUAL(statistics.epsilon_box_certifications,1u);
        ARIADNE_TEST_EQUAL(statistics.candidate_witness_searches,1u);
        ARIADNE_TEST_EQUAL(statistics.candidate_witness_successes,1u);

        SmtSolverTestSupport::accumulate_box_processing_statistics(
            statistics,Input{Status::SPLIT});
        ARIADNE_TEST_EQUAL(statistics.boxes_split,1u);

        SmtSolverTestSupport::accumulate_box_processing_statistics(
            statistics,Input{Status::UNKNOWN});
        ARIADNE_TEST_EQUAL(statistics.boxes_unknown,1u);
        ARIADNE_TEST_EQUAL(statistics.non_splittable_uncertified_boxes,1u);
        ARIADNE_TEST_EQUAL(statistics.non_splittable_epsilon_overlap_boxes,1u);

        SmtSolverTestSupport::accumulate_box_processing_statistics(
            statistics,Input{Status::EPSILON_SAT});

        ARIADNE_TEST_THROWS(
            SmtSolverTestSupport::accumulate_box_processing_statistics(
                statistics,Input{static_cast<Status>(999)}),
            std::runtime_error);
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
            std::cout << "[smt-solve] zero box budget returns UNKNOWN" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0.5_x),x[0],ValidatedNumber(0.5_x))
            });
            SmtResult solve_result=bounded_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_uncertified_boxes,0u);
        }


        {
            std::cout << "[smt-solve] classify non-splittable uncertified singleton" << std::endl;
            auto x=ValidatedScalarMultivariateFunction::coordinates(1);
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(
                1e-30_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            ExactBoxType domain({ExactIntervalType(1,1)});
            ValidatedScalarMultivariateFunction residual=sin(x[0])-sin(x[0]);
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),residual,ValidatedNumber(0))
            });
            SmtResult solve_result=tiny_epsilon_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_uncertified_boxes,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_epsilon_overlap_boxes,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
        }

        {
            std::cout << "[smt-solve] constant function UNSAT after direct range check" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            ValidatedScalarMultivariateFunction constant=
                ValidatedScalarMultivariateFunction::constant(
                    1u,ValidatedNumber(2));
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(0),
                    constant,
                    ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_pruned,1u);
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
            std::cout << "[smt-solve] endpoint witness avoids split for x^2=4 on [-2,2]" << std::endl;
            ExactBoxType domain({ExactIntervalType(-2,2)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(4),sqr(x[0]),ValidatedNumber(4))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_EQUAL(
                solve_result.witness()[0].lower_bound().raw(),
                UpperIntervalType(ExactIntervalType(-2,-2)).lower_bound().raw());
        }

        {
            std::cout << "[smt-solve] mixed-corner witness avoids split for x*y=-1" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            ExactBoxType domain({
                ExactIntervalType(-1,1),
                ExactIntervalType(-1,1)
            });
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(-1),
                    xy[0]*xy[1],
                    ValidatedNumber(-1))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_ASSERT(
                solve_result.witness()[0].lower_bound().raw()
                != solve_result.witness()[1].lower_bound().raw());
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
            std::cout << "[smt-solve] contractor fixpoint closes chained UNSAT at root" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            ExactBoxType domain({
                ExactIntervalType(0,1),
                ExactIntervalType(0,1)
            });
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),xy[0]-xy[1],ValidatedNumber(0)),
                ValidatedConstraint(ValidatedNumber(0),xy[1],ValidatedNumber(0)),
                ValidatedConstraint(ValidatedNumber(1),xy[0],ValidatedNumber(1))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_pruned,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
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
            std::cout << "[smt-solve] whole-box epsilon certification" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,0.0625_x)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),x[0],ValidatedNumber(0))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().epsilon_box_certifications,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
        }

        {
            std::cout << "[smt-solve] epsilon boundary EPSILON_SAT: x in [0,1], x=1.125 with epsilon=0.125" << std::endl;
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(1.125_x),x[0],ValidatedNumber(1.125_x))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
            if(solve_result.is_epsilon_sat()) {
                ARIADNE_TEST_ASSERT(solve_result.has_witness());
            }
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
            std::cout << "[smt-solve] shaving contractor handles repeated nonlinear occurrence" << std::endl;
            ExactBoxType domain({ExactIntervalType(-2,2)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(0.75_x),
                    sqr(x[0])+x[0],
                    ValidatedNumber(0.75_x))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_ASSERT(solve_result.statistics().hull_reduction_rounds>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().shaving_reduction_rounds>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().hull_effective_reductions
                <= solve_result.statistics().hull_reduction_rounds);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().shaving_effective_reductions
                <= solve_result.statistics().shaving_reduction_rounds);
            std::cout << "[smt-icp-stats] hull="
                      << solve_result.statistics().hull_effective_reductions
                      << "/" << solve_result.statistics().hull_reduction_rounds
                      << " shaving="
                      << solve_result.statistics().shaving_effective_reductions
                      << "/" << solve_result.statistics().shaving_reduction_rounds
                      << std::endl;
        }

        {
            std::cout << "[smt-solve] monotone contractor on validated smooth constraint" << std::endl;
            SmtSolver monotone_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                true,
                true));
            ExactBoxType domain({ExactIntervalType(0,2)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(3),
                    exp(x[0])+x[0],
                    ValidatedNumber(3))
            });
            SmtResult result=monotone_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(not result.is_unknown());
            ARIADNE_TEST_EQUAL(result.statistics().monotone_reduction_rounds,1u);
            ARIADNE_TEST_EQUAL(result.statistics().monotone_effective_reductions,1u);
        }

        {
            std::cout << "[smt-solve] monotone contractor skips nonmonotone validated coordinate" << std::endl;
            SmtSolver monotone_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false,
                true));
            ExactBoxType domain({ExactIntervalType(0,4)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(0),
                    sin(x[0]),
                    ValidatedNumber(0))
            });
            SmtResult result=monotone_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(not result.is_unknown());
            ARIADNE_TEST_EQUAL(result.statistics().monotone_reduction_rounds,1u);
            ARIADNE_TEST_EQUAL(result.statistics().monotone_effective_reductions,0u);
        }

        {
            std::cout << "[smt-solve] monotone contractor on smooth composed constraint" << std::endl;
            RealVariable mx("monotone_x");
            RealExpression emx=mx;
            RealSpace mspace({mx});
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal((exp(emx)+emx==3)));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> mliterals;
            mliterals.append(alternatives[0][0]);
            SmtSolver monotone_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                true,
                true));
            SmtResult monotone_result=monotone_solver.solve(
                mspace,ExactBoxType({ExactIntervalType(0,2)}),mliterals);
            ARIADNE_TEST_ASSERT(not monotone_result.is_unknown());
            ARIADNE_TEST_EQUAL(
                monotone_result.statistics().monotone_reduction_rounds,1u);
            ARIADNE_TEST_EQUAL(
                monotone_result.statistics().monotone_effective_reductions,1u);
        }

        {
            std::cout << "[smt-solve] monotone contractor skips nonmonotone theory coordinate" << std::endl;
            RealVariable mx("nonmonotone_theory_x");
            RealExpression emx=mx;
            RealSpace mspace({mx});
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(sin(emx)==0));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> mliterals;
            mliterals.append(alternatives[0][0]);
            SmtSolver monotone_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false,
                true));
            SmtResult result=monotone_solver.solve(
                mspace,ExactBoxType({ExactIntervalType(0,4)}),mliterals);
            ARIADNE_TEST_ASSERT(not result.is_unknown());
            ARIADNE_TEST_EQUAL(result.statistics().monotone_reduction_rounds,1u);
            ARIADNE_TEST_EQUAL(result.statistics().monotone_effective_reductions,0u);
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
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_pruned>0u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().hull_reduction_rounds>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().shaving_reduction_rounds>=1u);
        }

        {
            std::cout << "[smt-solve] split ignores wider inactive coordinate" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.01_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                1u,
                false));
            ExactBoxType domain({
                ExactIntervalType(-1,1),
                ExactIntervalType(-10000,10000)
            });
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(0.3_x),
                    sin(10*xy[0]),
                    ValidatedNumber(0.3_x))
            });
            SmtResult solve_result=bounded_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().sensitivity_guided_splits,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().sensitivity_overrides_geometric_splits,1u);
        }

        {
            std::cout << "[smt-solve] sensitivity-guided split is exercised" << std::endl;
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            SmtSolver split_solver(SmtSolverConfiguration(
                0.01_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                1u,
                false));
            ExactBoxType domain({
                ExactIntervalType(-1,1),
                ExactIntervalType(-10000,10000)
            });
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(0.3_x),
                    sin(10*xy[0])+0.000001_x*xy[1],
                    ValidatedNumber(0.3_x))
            });
            SmtResult solve_result=split_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().sensitivity_guided_splits,1u);
        }

        {
            std::cout << "[smt-solve] interior-point candidate search on 7D epsilon-SAT" << std::endl;
            auto coordinates=ValidatedScalarMultivariateFunction::coordinates(7);
            ValidatedScalarMultivariateFunction sum=coordinates[0];
            for(SizeType i=1u; i!=7u; ++i) {
                sum=sum+coordinates[i];
            }
            ExactBoxType domain({
                ExactIntervalType(0,1),ExactIntervalType(0,1),
                ExactIntervalType(0,1),ExactIntervalType(0,1),
                ExactIntervalType(0,1),ExactIntervalType(0,1),
                ExactIntervalType(0,1)
            });
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(1.3_x),sum,ValidatedNumber(1.3_x))
            });
            SmtResult solve_result=solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_successes,1u);
        }

        {
            std::cout << "[smt-solve] empty constraint list ignores zero box budget" << std::endl;
            SmtSolver zero_budget_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            ExactBoxType domain({ExactIntervalType(-1,1)});
            List<ValidatedConstraint> constraints;
            SmtResult solve_result=zero_budget_solver.solve(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
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
            std::cout << "[smt-solve] non-splittable uncertified singleton returns UNKNOWN" << std::endl;
            auto sx=ValidatedScalarMultivariateFunction::coordinates(1);
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(
                1e-30_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            ValidatedScalarMultivariateFunction residual=
                sqr(sin(sx[0]))+sqr(cos(sx[0]))-1;
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),residual,ValidatedNumber(0))
            });
            SmtResult solve_result=tiny_epsilon_solver.solve(
                ExactBoxType({ExactIntervalType(1,1)}),constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_unknown,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_uncertified_boxes,1u);
        }

        {
            std::cout << "[smt-solve] terminal box skips interior-point candidate search" << std::endl;
            auto sx=ValidatedScalarMultivariateFunction::coordinates(1);
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(1e-30_x));
            ValidatedScalarMultivariateFunction residual=
                sqr(sin(sx[0]))+sqr(cos(sx[0]))-1;
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0),residual,ValidatedNumber(0))
            });
            SmtResult solve_result=tiny_epsilon_solver.solve(
                ExactBoxType({ExactIntervalType(1,1)}),constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_successes,0u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_uncertified_boxes,1u);
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

        {
            std::cout << "[smt-theory-solve] empty domain is UNSAT" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals;
            SmtResult solve_result=solver.solve(
                space,
                ExactBoxType({ExactIntervalType::empty_interval()}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
        }

        auto primitive = [&](ContinuousPredicate const& predicate) {
            auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            return alternatives[0][0];
        };

        {
            std::cout << "[smt-theory-solve] empty primitive list ignores zero box budget" << std::endl;
            SmtSolver zero_budget_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            List<SmtTheoryPrimitiveLiteral> literals;
            SmtResult solve_result=zero_budget_solver.solve(
                space,
                ExactBoxType({ExactIntervalType(-1,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
        }

        {
            std::cout << "[smt-theory-solve] simplified zero primitive ignores zero box budget" << std::endl;
            RealVariable zero_x("zero_x");
            RealExpression zero_ex=zero_x;
            RealSpace zero_space({zero_x});
            SmtSolver zero_budget_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u,
                false));
            RealExpression residual=sin(zero_ex)-sin(zero_ex);
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(residual>0));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=zero_budget_solver.solve(
                zero_space,
                ExactBoxType({ExactIntervalType(1,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
        }

        {
            std::cout << "[smt-theory-solve] zero box budget propagates UNKNOWN" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex==0)});
            SmtResult solve_result=bounded_solver.solve(
                space,
                ExactBoxType({ExactIntervalType(0,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_uncertified_boxes,0u);
        }

        {
            std::cout << "[smt-theory-solve] sequential terminal uncertainty propagates UNKNOWN" << std::endl;
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(1e-30_x));
            RealExpression residual=sqr(sin(ex))+sqr(cos(ex))-1;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(residual==0)});
            SmtResult solve_result=tiny_epsilon_solver.solve(
                space,ExactBoxType({ExactIntervalType(1,1)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_unknown,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_uncertified_boxes,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().candidate_witness_searches,0u);
        }

        {
            std::cout << "[smt-theory-solve] shaving proves dependency-hidden strict UNSAT" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({
                primitive(ex*(1-ex)>0.32_x)
            });
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_pruned,1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().shaving_reduction_rounds>=1u);
        }

        {
            std::cout << "[smt-theory-solve] later hull reduction invalidates earlier strict literal" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({
                primitive(sin(ex)>0),
                primitive(ex==4)
            });
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,7)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_pruned,1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().hull_effective_reductions>=1u);
        }

        {
            std::cout << "[smt-theory-solve] strict epsilon predicate rejects negative point" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex>0)});
            ARIADNE_TEST_ASSERT(
                not SmtSolverTestSupport::epsilon_satisfied(
                    solver,
                    space,
                    UpperBoxType({UpperIntervalType(ExactIntervalType(-0.25_x,-0.25_x))}),
                    literals));
        }


        {
            std::cout << "[smt-theory-solve] simplify repeated expression before interval solving" << std::endl;
            RealVariable singleton_x("singleton_x");
            RealExpression singleton_ex=singleton_x;
            RealSpace singleton_space({singleton_x});
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(
                1e-30_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            RealExpression residual=sin(singleton_ex)-sin(singleton_ex);
            ARIADNE_TEST_ASSERT(identical(simplify(residual),RealExpression(0)));
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(residual==0));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=tiny_epsilon_solver.solve(
                singleton_space,
                ExactBoxType({ExactIntervalType(1,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_unknown,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
        }

        {
            std::cout << "[smt-theory-solve] EQ at epsilon boundary: x=0 weakened on x=0.125" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(ex==0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(0.125_x,0.125_x)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
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
            ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
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
            ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
        }

        {
            std::cout << "[smt-theory-solve] whole-box epsilon certification" << std::endl;
            RealVariable box_x("box_x");
            RealSpace box_space({box_x});
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(RealExpression(box_x)==0));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=solver.solve(
                box_space,
                ExactBoxType({ExactIntervalType(0,0.0625_x)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().epsilon_box_certifications,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
        }

        {
            std::cout << "[smt-theory-solve] contractor fixpoint closes chained UNSAT at root" << std::endl;
            RealVariable y("y");
            RealExpression ey=y;
            RealSpace xy_space({x,y});
            auto xy_primitive = [&](ContinuousPredicate const& predicate) {
                auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
                ARIADNE_TEST_EQUAL(alternatives.size(),1u);
                ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
                return alternatives[0][0];
            };
            List<SmtTheoryPrimitiveLiteral> literals({
                xy_primitive(ex-ey==0),
                xy_primitive(ey==0),
                xy_primitive(ex==1)
            });
            SmtResult solve_result=solver.solve(
                xy_space,
                ExactBoxType({ExactIntervalType(0,1),ExactIntervalType(0,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_pruned,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
        }

        {
            std::cout << "[smt-theory-solve] shaving contractor handles repeated nonlinear occurrence" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(sqr(ex)+ex==0.75_x)});
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(-2,2)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_ASSERT(solve_result.statistics().hull_reduction_rounds>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().shaving_reduction_rounds>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().hull_effective_reductions
                <= solve_result.statistics().hull_reduction_rounds);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().shaving_effective_reductions
                <= solve_result.statistics().shaving_reduction_rounds);
            std::cout << "[smt-icp-stats] hull="
                      << solve_result.statistics().hull_effective_reductions
                      << "/" << solve_result.statistics().hull_reduction_rounds
                      << " shaving="
                      << solve_result.statistics().shaving_effective_reductions
                      << "/" << solve_result.statistics().shaving_reduction_rounds
                      << std::endl;
        }

        {
            std::cout << "[smt-theory-solve] transcendental GT: sin(x)>0 on [3,4]" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(sin(ex)>0)});
            SmtResult solve_result=solver.solve(space,ExactBoxType({ExactIntervalType(3,4)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
        }

        {
            std::cout << "[smt-theory-solve] strict GT rejects zero-only image" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(sqr(ex)>0)});
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(0,0)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_pruned,1u);
        }

        {
            std::cout << "[smt-theory-solve] endpoint witness avoids split for x^2=4 on [-2,2]" << std::endl;
            List<SmtTheoryPrimitiveLiteral> literals({primitive(sqr(ex)==4)});
            SmtResult solve_result=solver.solve(
                space,ExactBoxType({ExactIntervalType(-2,2)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
        }

        {
            std::cout << "[smt-theory-solve] mixed-corner witness avoids split for x*y=-1" << std::endl;
            RealVariable y("y");
            RealExpression ey=y;
            RealSpace xy_space({x,y});
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(ex*ey==-1));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=solver.solve(
                xy_space,
                ExactBoxType({ExactIntervalType(-1,1),ExactIntervalType(-1,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
        }

        {
            std::cout << "[smt-theory-solve] cover supported RealExpression operations" << std::endl;

            auto expect_status = [&](String const& label,
                                     ExactIntervalType const& interval,
                                     ContinuousPredicate const& predicate,
                                     SmtResultStatus expected) {
                std::cout << "[smt-real-op] " << label << std::endl;
                List<SmtTheoryPrimitiveLiteral> literals({primitive(predicate)});
                SmtResult solve_result=solver.solve(space,ExactBoxType({interval}),literals);
                ARIADNE_TEST_EQUAL(solve_result.status(),expected);
            };

            expect_status("neg SAT",ExactIntervalType(1,1),neg(ex)==-1,SmtResultStatus::EPSILON_SAT);
            expect_status("neg UNSAT",ExactIntervalType(1,1),neg(ex)==1,SmtResultStatus::UNSAT);

            expect_status("sqr SAT",ExactIntervalType(2,2),sqr(ex)==4,SmtResultStatus::EPSILON_SAT);
            expect_status("sqr UNSAT",ExactIntervalType(2,2),sqr(ex)==5,SmtResultStatus::UNSAT);
            {
                std::cout << "[smt-real-op] sqr epsilon-overlap delta-decision" << std::endl;
                auto alternatives=normalize_smt_theory_literal(
                    make_smt_theory_literal(sqr(ex)==4.0625_x));
                ARIADNE_TEST_EQUAL(alternatives.size(),1u);
                ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
                List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
                SmtResult solve_result=solver.solve(
                    space,ExactBoxType({ExactIntervalType(2,2)}),literals);
                ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
                ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
            }

            expect_status("pow SAT",ExactIntervalType(2,2),pow(ex,3)==8,SmtResultStatus::EPSILON_SAT);
            expect_status("pow UNSAT",ExactIntervalType(2,2),pow(ex,3)==9,SmtResultStatus::UNSAT);

            expect_status("rec SAT",ExactIntervalType(2,2),rec(ex)==0.5_x,SmtResultStatus::EPSILON_SAT);
            expect_status("rec UNSAT",ExactIntervalType(2,2),rec(ex)==1,SmtResultStatus::UNSAT);

            expect_status("sqrt SAT",ExactIntervalType(4,4),sqrt(ex)==2,SmtResultStatus::EPSILON_SAT);
            expect_status("sqrt UNSAT",ExactIntervalType(4,4),sqrt(ex)==3,SmtResultStatus::UNSAT);

            expect_status("exp SAT",ExactIntervalType(0,0),exp(ex)==1,SmtResultStatus::EPSILON_SAT);
            expect_status("exp UNSAT",ExactIntervalType(0,0),exp(ex)==2,SmtResultStatus::UNSAT);

            expect_status("log SAT",ExactIntervalType(1,1),log(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("log UNSAT",ExactIntervalType(1,1),log(ex)==1,SmtResultStatus::UNSAT);

            expect_status("sin SAT",ExactIntervalType(0,0),sin(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("sin UNSAT",ExactIntervalType(0,0),sin(ex)==1,SmtResultStatus::UNSAT);

            expect_status("cos SAT",ExactIntervalType(0,0),cos(ex)==1,SmtResultStatus::EPSILON_SAT);
            expect_status("cos UNSAT",ExactIntervalType(0,0),cos(ex)==0,SmtResultStatus::UNSAT);

            expect_status("tan SAT",ExactIntervalType(0,0),tan(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("tan UNSAT",ExactIntervalType(0,0),tan(ex)==1,SmtResultStatus::UNSAT);

            expect_status("asin SAT",ExactIntervalType(0,0),asin(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("asin UNSAT",ExactIntervalType(0,0),asin(ex)==1,SmtResultStatus::UNSAT);

            expect_status("acos SAT",ExactIntervalType(1,1),acos(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("acos UNSAT",ExactIntervalType(1,1),acos(ex)==1,SmtResultStatus::UNSAT);

            expect_status("atan SAT",ExactIntervalType(0,0),atan(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("atan UNSAT",ExactIntervalType(0,0),atan(ex)==1,SmtResultStatus::UNSAT);

            expect_status("max SAT",ExactIntervalType(-1,-1),max(ex,RealExpression(0))==0,SmtResultStatus::EPSILON_SAT);
            expect_status("max UNSAT",ExactIntervalType(-1,-1),max(ex,RealExpression(0))==1,SmtResultStatus::UNSAT);

            expect_status("min SAT",ExactIntervalType(1,1),min(ex,RealExpression(0))==0,SmtResultStatus::EPSILON_SAT);
            expect_status("min UNSAT",ExactIntervalType(1,1),min(ex,RealExpression(0))==1,SmtResultStatus::UNSAT);

            expect_status("abs SAT",ExactIntervalType(-1,-1),abs(ex)==1,SmtResultStatus::EPSILON_SAT);
            expect_status("abs UNSAT",ExactIntervalType(-1,-1),abs(ex)==0,SmtResultStatus::UNSAT);

            expect_status("nul SAT",ExactIntervalType(2,2),nul(ex)==0,SmtResultStatus::EPSILON_SAT);
            expect_status("nul UNSAT",ExactIntervalType(2,2),nul(ex)==1,SmtResultStatus::UNSAT);

            expect_status("pos SAT",ExactIntervalType(2,2),pos(ex)==2,SmtResultStatus::EPSILON_SAT);
            expect_status("pos UNSAT",ExactIntervalType(2,2),pos(ex)==3,SmtResultStatus::UNSAT);

            expect_status("hlf SAT",ExactIntervalType(2,2),hlf(ex)==1,SmtResultStatus::EPSILON_SAT);
            expect_status("hlf UNSAT",ExactIntervalType(2,2),hlf(ex)==2,SmtResultStatus::UNSAT);

            expect_status("add SAT",ExactIntervalType(2,2),add(ex,RealExpression(3))==5,SmtResultStatus::EPSILON_SAT);
            expect_status("add UNSAT",ExactIntervalType(2,2),add(ex,RealExpression(3))==6,SmtResultStatus::UNSAT);

            expect_status("sub SAT",ExactIntervalType(2,2),sub(ex,RealExpression(3))==-1,SmtResultStatus::EPSILON_SAT);
            expect_status("sub UNSAT",ExactIntervalType(2,2),sub(ex,RealExpression(3))==0,SmtResultStatus::UNSAT);

            expect_status("mul SAT",ExactIntervalType(2,2),mul(ex,RealExpression(3))==6,SmtResultStatus::EPSILON_SAT);
            expect_status("mul UNSAT",ExactIntervalType(2,2),mul(ex,RealExpression(3))==7,SmtResultStatus::UNSAT);

            expect_status("div SAT",ExactIntervalType(6,6),div(ex,RealExpression(3))==2,SmtResultStatus::EPSILON_SAT);
            expect_status("div UNSAT",ExactIntervalType(6,6),div(ex,RealExpression(3))==3,SmtResultStatus::UNSAT);

            expect_status("operator add SAT",ExactIntervalType(2,2),(ex+RealExpression(3))==5,SmtResultStatus::EPSILON_SAT);
            expect_status("operator sub SAT",ExactIntervalType(2,2),(ex-RealExpression(3))==-1,SmtResultStatus::EPSILON_SAT);
            expect_status("operator mul SAT",ExactIntervalType(2,2),(ex*RealExpression(3))==6,SmtResultStatus::EPSILON_SAT);
            expect_status("operator div SAT",ExactIntervalType(6,6),(ex/RealExpression(3))==2,SmtResultStatus::EPSILON_SAT);
        }

        {
            std::cout << "[smt-theory-solve] cover binary RealExpression operations" << std::endl;
            RealVariable y("y");
            RealExpression ey=y;
            RealSpace xy_space({x,y});

            auto xy_primitive = [&](ContinuousPredicate const& predicate) {
                auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
                ARIADNE_TEST_EQUAL(alternatives.size(),1u);
                ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
                return alternatives[0][0];
            };

            auto expect_xy_status = [&](String const& label,
                                        ExactIntervalType const& x_interval,
                                        ExactIntervalType const& y_interval,
                                        ContinuousPredicate const& predicate,
                                        SmtResultStatus expected) {
                std::cout << "[smt-real-binary-op] " << label << std::endl;
                List<SmtTheoryPrimitiveLiteral> literals({xy_primitive(predicate)});
                SmtResult solve_result=solver.solve(
                    xy_space,ExactBoxType({x_interval,y_interval}),literals);
                ARIADNE_TEST_EQUAL(solve_result.status(),expected);
            };

            expect_xy_status("add SAT",ExactIntervalType(2,2),ExactIntervalType(3,3),ex+ey==5,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("add UNSAT",ExactIntervalType(2,2),ExactIntervalType(3,3),ex+ey==6,SmtResultStatus::UNSAT);

            expect_xy_status("sub SAT",ExactIntervalType(2,2),ExactIntervalType(3,3),ex-ey==-1,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("sub UNSAT",ExactIntervalType(2,2),ExactIntervalType(3,3),ex-ey==0,SmtResultStatus::UNSAT);

            expect_xy_status("mul SAT",ExactIntervalType(2,2),ExactIntervalType(3,3),ex*ey==6,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("mul UNSAT",ExactIntervalType(2,2),ExactIntervalType(3,3),ex*ey==7,SmtResultStatus::UNSAT);

            expect_xy_status("div SAT",ExactIntervalType(6,6),ExactIntervalType(3,3),ex/ey==2,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("div UNSAT",ExactIntervalType(6,6),ExactIntervalType(3,3),ex/ey==3,SmtResultStatus::UNSAT);

            expect_xy_status("max left SAT",ExactIntervalType(3,3),ExactIntervalType(2,2),max(ex,ey)==3,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("max right SAT",ExactIntervalType(2,2),ExactIntervalType(3,3),max(ex,ey)==3,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("max UNSAT",ExactIntervalType(2,2),ExactIntervalType(3,3),max(ex,ey)==2,SmtResultStatus::UNSAT);

            expect_xy_status("min left SAT",ExactIntervalType(2,2),ExactIntervalType(3,3),min(ex,ey)==2,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("min right SAT",ExactIntervalType(3,3),ExactIntervalType(2,2),min(ex,ey)==2,SmtResultStatus::EPSILON_SAT);
            expect_xy_status("min UNSAT",ExactIntervalType(2,2),ExactIntervalType(3,3),min(ex,ey)==3,SmtResultStatus::UNSAT);
        }

        {
            std::cout << "[smt-theory-solve] cover RealExpression domain edge cases" << std::endl;

            auto expect_edge_status = [&](String const& label,
                                          ExactIntervalType const& interval,
                                          ContinuousPredicate const& predicate,
                                          SmtResultStatus expected) {
                std::cout << "[smt-real-domain] " << label << std::endl;
                List<SmtTheoryPrimitiveLiteral> literals({primitive(predicate)});
                SmtResult solve_result=solver.solve(space,ExactBoxType({interval}),literals);
                ARIADNE_TEST_EQUAL(solve_result.status(),expected);
            };

            expect_edge_status("sqrt negative domain",ExactIntervalType(-1,-1),sqrt(ex)==0,SmtResultStatus::UNSAT);
            expect_edge_status("sqrt zero boundary",ExactIntervalType(0,0),sqrt(ex)==0,SmtResultStatus::EPSILON_SAT);

            expect_edge_status("log negative domain",ExactIntervalType(-1,-1),log(ex)==0,SmtResultStatus::UNSAT);
            expect_edge_status("log zero domain",ExactIntervalType(0,0),log(ex)==0,SmtResultStatus::UNSAT);

            expect_edge_status("asin below domain",ExactIntervalType(-2,-2),asin(ex)==0,SmtResultStatus::UNSAT);
            expect_edge_status("asin above domain",ExactIntervalType(2,2),asin(ex)==0,SmtResultStatus::UNSAT);
            expect_edge_status("acos below domain",ExactIntervalType(-2,-2),acos(ex)==0,SmtResultStatus::UNSAT);
            expect_edge_status("acos above domain",ExactIntervalType(2,2),acos(ex)==0,SmtResultStatus::UNSAT);

            expect_edge_status("rec positive near zero",ExactIntervalType(0.5_x,0.5_x),rec(ex)==2,SmtResultStatus::EPSILON_SAT);
            expect_edge_status("rec negative",ExactIntervalType(-1,-1),rec(ex)==-1,SmtResultStatus::EPSILON_SAT);
            expect_edge_status("rec singleton zero undefined",ExactIntervalType(0,0),rec(ex)==0,SmtResultStatus::UNSAT);

            expect_edge_status("tan regular branch",ExactIntervalType(0,0),tan(ex)==0,SmtResultStatus::EPSILON_SAT);

            RealVariable y("y");
            RealExpression ey=y;
            RealSpace xy_space({x,y});
            auto xy_primitive = [&](ContinuousPredicate const& predicate) {
                auto alternatives=normalize_smt_theory_literal(make_smt_theory_literal(predicate));
                ARIADNE_TEST_EQUAL(alternatives.size(),1u);
                ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
                return alternatives[0][0];
            };
            {
                std::cout << "[smt-real-domain] division singleton zero denominator" << std::endl;
                List<SmtTheoryPrimitiveLiteral> literals({xy_primitive(ex/ey==1)});
                SmtResult solve_result=solver.solve(
                    xy_space,
                    ExactBoxType({ExactIntervalType(1,1),ExactIntervalType(0,0)}),
                    literals);
                ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            }
        }

        {
            std::cout << "[smt-theory-solve] split ignores wider inactive coordinate" << std::endl;
            RealVariable y("y");
            RealSpace xy_space({x,y});
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.01_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                1u,
                false));
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(sin(10*ex)==0.3_x));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=bounded_solver.solve(
                xy_space,
                ExactBoxType({ExactIntervalType(-1,1),ExactIntervalType(-10000,10000)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().sensitivity_guided_splits,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().sensitivity_overrides_geometric_splits,1u);
        }

        {
            std::cout << "[smt-theory-solve] sensitivity-guided split is exercised" << std::endl;
            RealVariable y("y");
            RealExpression ey=y;
            RealSpace xy_space({x,y});
            SmtSolver split_solver(SmtSolverConfiguration(
                0.01_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                1u,
                false));
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(sin(10*ex)+0.000001_x*ey==0.3_x));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=split_solver.solve(
                xy_space,
                ExactBoxType({ExactIntervalType(-1,1),ExactIntervalType(-10000,10000)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().sensitivity_guided_splits,1u);
        }

        {
            std::cout << "[smt-theory-solve] interior-point candidate search on 7D epsilon-SAT" << std::endl;
            RealVariable x0("x0"); RealVariable x1("x1"); RealVariable x2("x2");
            RealVariable x3("x3"); RealVariable x4("x4"); RealVariable x5("x5");
            RealVariable x6("x6");
            RealSpace seven_space({x0,x1,x2,x3,x4,x5,x6});
            RealExpression sum=
                RealExpression(x0)+RealExpression(x1)+RealExpression(x2)
                +RealExpression(x3)+RealExpression(x4)+RealExpression(x5)
                +RealExpression(x6);
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(sum==1.3_x));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=solver.solve(
                seven_space,
                ExactBoxType({
                    ExactIntervalType(0,1),ExactIntervalType(0,1),
                    ExactIntervalType(0,1),ExactIntervalType(0,1),
                    ExactIntervalType(0,1),ExactIntervalType(0,1),
                    ExactIntervalType(0,1)
                }),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_successes,1u);
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
            std::cout << "[smt-dpll] symbolic simplification removes dependency before theory solving" << std::endl;
            RealExpression residual=sin(ex)-sin(ex);
            SmtResult solve_result=solver.solve(
                space,
                ExactBoxType({ExactIntervalType(1,1)}),
                residual==0);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_unknown,0u);
        }

        {
            std::cout << "[smt-dpll] simplified zero theory atom ignores zero box budget" << std::endl;
            SmtSolver zero_budget_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u,
                false));
            RealExpression residual=sin(ex)-sin(ex);
            SmtResult solve_result=zero_budget_solver.solve(
                space,
                ExactBoxType({ExactIntervalType(1,1)}),
                residual>0);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_checks,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
        }

        {
            std::cout << "[smt-dpll] Boolean theory solve preserves disabled candidate search" << std::endl;
            SmtSolver no_candidate_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            ContinuousPredicate formula=(sqr(ex)==0.75_x);
            SmtResult solve_result=no_candidate_solver.solve(
                space,ExactBoxType({ExactIntervalType(-2,2)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_successes,0u);
        }

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
            std::cout << "[smt-dpll] zero theory box budget still permits Boolean UNSAT" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            ContinuousPredicate atom=(ex>=0);
            SmtResult solve_result=bounded_solver.solve(
                space,
                ExactBoxType({ExactIntervalType(-1,1)}),
                atom&&!atom);
            ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_checks,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_conflicts>=1u);
        }

        {
            std::cout << "[smt-dpll] zero global box budget stops before theory" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            ContinuousPredicate formula=(ex>=0);
            SmtResult solve_result=bounded_solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_checks,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_conflicts,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().theory_learned_clauses,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_uncertified_boxes,0u);
        }

        {
            std::cout << "[smt-dpll] one-box global budget permits one theory check" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                1u));
            ContinuousPredicate formula=(ex>=0);
            SmtResult solve_result=bounded_solver.solve(
                space,ExactBoxType({ExactIntervalType(0,1)}),formula);
            ARIADNE_TEST_ASSERT(
                solve_result.is_epsilon_sat() || solve_result.is_unknown());
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_checks>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed<=1u);
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
            std::cout << "[smt-dpll] bounded learned database remains stable" << std::endl;
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
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clauses>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clause_activity_bumps>=1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().peak_active_non_theory_learned_clauses>=1u);
        }

        {
            std::cout << "[smt-dpll] conservative pruning preserves protected learned clauses" << std::endl;
            SmtSolver pruning_solver(SmtSolverConfiguration(
                0.125_x,std::numeric_limits<SizeType>::max(),1u));
            ContinuousPredicate a=(ex>=-0.75_x);
            ContinuousPredicate b=(ex>=-0.25_x);
            ContinuousPredicate c=(ex>=0.25_x);
            ContinuousPredicate d=(ex>=0.75_x);
            ContinuousPredicate formula=
                ( a|| b|| c|| d)&&
                ( a|| b|| c||(!d))&&
                ( a|| b||(!c)|| d)&&
                ( a|| b||(!c)||(!d))&&
                ( a||(!b)|| c|| d)&&
                ( a||(!b)|| c||(!d))&&
                ( a||(!b)||(!c)|| d)&&
                ( a||(!b)||(!c)||(!d))&&
                ((!a)|| b|| c|| d)&&
                ((!a)|| b|| c||(!d))&&
                ((!a)|| b||(!c)|| d)&&
                ((!a)|| b||(!c)||(!d))&&
                ((!a)||(!b)|| c|| d)&&
                ((!a)||(!b)|| c||(!d))&&
                ((!a)||(!b)||(!c)|| d)&&
                ((!a)||(!b)||(!c)||(!d));
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
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clauses>=4u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().learned_clause_pruning_runs>=1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().learned_clauses_pruned,0u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().peak_active_non_theory_learned_clauses>=2u);
        }

        {
            std::cout << "[smt-dpll] semantically equivalent atoms share one Boolean variable" << std::endl;
            RealVariable y("y");
            RealExpression ey=y;
            RealSpace xy_space({x,y});
            ExactBoxType xy_domain({ExactIntervalType(-1,1),ExactIntervalType(-1,1)});

            SmtResult reversed_inequality=solver.solve(
                xy_space,xy_domain,(ex<=ey)&&!(ey>=ex));
            ARIADNE_TEST_ASSERT(reversed_inequality.is_unsat());
            ARIADNE_TEST_EQUAL(reversed_inequality.statistics().theory_checks,0u);

            SmtResult sign_comparison=solver.solve(
                space,
                ExactBoxType({ExactIntervalType(-1,1)}),
                sgn(ex)&&!(ex>0));
            ARIADNE_TEST_ASSERT(sign_comparison.is_unsat());
            ARIADNE_TEST_EQUAL(sign_comparison.statistics().theory_checks,0u);
        }

        {
            std::cout << "[smt-dpll] structurally identical atoms share one Boolean variable" << std::endl;
            ContinuousPredicate first=(ex<=0);
            ContinuousPredicate second=(ex<=0);
            ARIADNE_TEST_ASSERT(first.node_raw_ptr()!=second.node_raw_ptr());

            SmtResult duplicate_atom_contradiction=solver.solve(
                space,
                ExactBoxType({ExactIntervalType(-1,1)}),
                first&&!second);
            ARIADNE_TEST_ASSERT(duplicate_atom_contradiction.is_unsat());
            ARIADNE_TEST_EQUAL(duplicate_atom_contradiction.statistics().theory_checks,0u);
        }

        {
            std::cout << "[smt-dpll] Boolean constant folding avoids theory work" << std::endl;
            ContinuousPredicate atom=(ex>=0);

            SmtResult folded_false=solver.solve(
                space,
                ExactBoxType({ExactIntervalType(-1,1)}),
                ContinuousPredicate(false)&&atom);
            ARIADNE_TEST_ASSERT(folded_false.is_unsat());
            ARIADNE_TEST_EQUAL(folded_false.statistics().theory_checks,0u);

            SmtResult folded_true=solver.solve(
                space,
                ExactBoxType({ExactIntervalType(-1,1)}),
                ContinuousPredicate(true)||atom);
            ARIADNE_TEST_ASSERT(folded_true.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(folded_true.has_witness());
            ARIADNE_TEST_EQUAL(folded_true.statistics().theory_checks,0u);

            SmtResult folded_true_empty=solver.solve(
                space,
                ExactBoxType({ExactIntervalType::empty_interval()}),
                ContinuousPredicate(true)||atom);
            ARIADNE_TEST_ASSERT(folded_true_empty.is_unsat());
            ARIADNE_TEST_EQUAL(folded_true_empty.statistics().theory_checks,0u);
        }

        {
            std::cout << "[smt-dpll] Boolean constants" << std::endl;
            ExactBoxType domain({ExactIntervalType(-1,1)});

            SmtResult true_result=solver.solve(space,domain,ContinuousPredicate(true));
            ARIADNE_TEST_ASSERT(true_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(true_result.has_witness());

            SmtResult false_result=solver.solve(space,domain,ContinuousPredicate(false));
            ARIADNE_TEST_ASSERT(false_result.is_unsat());

            ContinuousPredicate atom=(ex>=0);
            SmtResult mixed_true=solver.solve(
                space,domain,ContinuousPredicate(true)&&atom);
            ARIADNE_TEST_ASSERT(mixed_true.is_epsilon_sat());

            SmtResult mixed_false=solver.solve(
                space,domain,ContinuousPredicate(false)&&atom);
            ARIADNE_TEST_ASSERT(mixed_false.is_unsat());

            SmtResult neutral_or=solver.solve(
                space,ExactBoxType({ExactIntervalType(-1,-0.5_x)}),
                ContinuousPredicate(false)||atom);
            ARIADNE_TEST_ASSERT(neutral_or.is_unsat());
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
            ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
            if(solve_result.is_epsilon_sat()) {
                ARIADNE_TEST_ASSERT(solve_result.has_witness());
            }
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
            std::cout << "[smt-dpll] both Boolean branches exhaust on theory UNKNOWN" << std::endl;
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(
                1e-30_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            ContinuousPredicate uncertain=
                (sqr(sin(ex))+sqr(cos(ex))-1==0);
            SmtResult solve_result=tiny_epsilon_solver.solve(
                space,
                ExactBoxType({ExactIntervalType(1,1)}),
                uncertain||(!uncertain));
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_decisions>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_backtracks>=2u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_checks>=2u);
        }

        {
            std::cout << "[smt-dpll] global box budget stops further Boolean search" << std::endl;
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                2u));
            ContinuousPredicate formula=(ex<0)||(ex>1);
            SmtResult solve_result=bounded_solver.solve(
                space,ExactBoxType({ExactIntervalType(0.5_x,0.5_x)}),formula);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_ASSERT(solve_result.statistics().theory_checks>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed<=2u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_processed>=1u);
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
            ARIADNE_TEST_ASSERT(solve_result.statistics().boolean_conflicts_analyzed>=1u);
            ARIADNE_TEST_ASSERT(solve_result.statistics().last_learned_clause_literals>=1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().last_learned_current_level_literals,1u);
            ARIADNE_TEST_ASSERT(
                solve_result.statistics().last_backjump_level
                < solve_result.statistics().max_decision_level);
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
            ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
        }

        {
            std::cout << "[smt-dpll] audit relation and negation delta-boundaries" << std::endl;

            auto require_epsilon_sat = [&](String const& label,
                                           ContinuousPredicate const& formula,
                                           ExactDouble value) {
                std::cout << "[smt-relation] " << label << " value=" << value << std::endl;
                SmtResult solve_result=solver.solve(
                    space,ExactBoxType({ExactIntervalType(value,value)}),formula);
                ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
                ARIADNE_TEST_ASSERT(solve_result.has_witness());
            };

            auto require_unsat = [&](String const& label,
                                     ContinuousPredicate const& formula,
                                     ExactDouble value) {
                std::cout << "[smt-relation] " << label << " value=" << value << std::endl;
                SmtResult solve_result=solver.solve(
                    space,ExactBoxType({ExactIntervalType(value,value)}),formula);
                ARIADNE_TEST_ASSERT(solve_result.is_unsat());
            };

            auto require_delta_overlap = [&](String const& label,
                                             ContinuousPredicate const& formula,
                                             ExactDouble value) {
                std::cout << "[smt-relation] " << label << " value=" << value << std::endl;
                SmtResult solve_result=solver.solve(
                    space,ExactBoxType({ExactIntervalType(value,value)}),formula);
                ARIADNE_TEST_ASSERT(solve_result.is_unsat() || solve_result.is_epsilon_sat());
                ARIADNE_TEST_ASSERT(not solve_result.is_unknown());
            };

            // Equality: exact models must be delta-sat; outside the epsilon band must be UNSAT.
            require_epsilon_sat("x=0 original true",ex==0,0.0_x);
            require_delta_overlap("x=0 positive epsilon boundary",ex==0,0.125_x);
            require_delta_overlap("x=0 negative epsilon boundary",ex==0,-0.125_x);
            require_unsat("x=0 outside epsilon",ex==0,0.25_x);

            // Non-strict inequalities use closed epsilon boundaries.
            require_epsilon_sat("x>=0 original true",ex>=0,0.25_x);
            require_delta_overlap("x>=0 epsilon boundary",ex>=0,-0.125_x);
            require_unsat("x>=0 outside epsilon",ex>=0,-0.25_x);

            require_epsilon_sat("x<=0 original true",ex<=0,-0.25_x);
            require_delta_overlap("x<=0 epsilon boundary",ex<=0,0.125_x);
            require_unsat("x<=0 outside epsilon",ex<=0,0.25_x);

            // Strict inequalities retain strictness after weakening.
            require_epsilon_sat("x>0 original true",ex>0,0.25_x);
            require_delta_overlap("x>0 weakened interior",ex>0,-0.0625_x);
            require_unsat("x>0 epsilon boundary excluded",ex>0,-0.125_x);

            require_epsilon_sat("x<0 original true",ex<0,-0.25_x);
            require_delta_overlap("x<0 weakened interior",ex<0,0.0625_x);
            require_unsat("x<0 epsilon boundary excluded",ex<0,0.125_x);

            // Sign predicates are strict positivity tests and use the same weakening as x>0.
            require_epsilon_sat("sgn(x) original true",sgn(ex),0.25_x);
            require_delta_overlap("sgn(x) weakened interior",sgn(ex),-0.0625_x);
            require_unsat("sgn(x) epsilon boundary excluded",sgn(ex),-0.125_x);

            // Disequality normalizes to x>0 OR -x>0. Its delta weakening is
            // therefore satisfied at equality as well.
            require_epsilon_sat("x!=0 original true positive",ex!=0,0.25_x);
            require_epsilon_sat("x!=0 original true negative",ex!=0,-0.25_x);
            require_delta_overlap("x!=0 equality delta-overlap",ex!=0,0.0_x);

            // Audit Boolean negation against the corresponding complementary relation.
            require_epsilon_sat("!(x<=0) original true",!(ex<=0),0.25_x);
            require_delta_overlap("!(x<=0) weakened interior",!(ex<=0),-0.0625_x);
            require_unsat("!(x<=0) strict epsilon boundary",!(ex<=0),-0.125_x);

            require_epsilon_sat("!(x>=0) original true",!(ex>=0),-0.25_x);
            require_delta_overlap("!(x>=0) weakened interior",!(ex>=0),0.0625_x);
            require_unsat("!(x>=0) strict epsilon boundary",!(ex>=0),0.125_x);

            require_epsilon_sat("!(x<0) original true",!(ex<0),0.25_x);
            require_delta_overlap("!(x<0) closed epsilon boundary",!(ex<0),-0.125_x);
            require_unsat("!(x<0) outside epsilon",!(ex<0),-0.25_x);

            require_epsilon_sat("!(x>0) original true",!(ex>0),-0.25_x);
            require_delta_overlap("!(x>0) closed epsilon boundary",!(ex>0),0.125_x);
            require_unsat("!(x>0) outside epsilon",!(ex>0),0.25_x);

            require_epsilon_sat("!sgn(x) original true",!sgn(ex),-0.25_x);
            require_delta_overlap("!sgn(x) closed epsilon boundary",!sgn(ex),0.125_x);
            require_unsat("!sgn(x) outside epsilon",!sgn(ex),0.25_x);

            require_epsilon_sat("!(x!=0) equality",!(ex!=0),0.0_x);
            require_delta_overlap("!(x!=0) epsilon boundary",!(ex!=0),0.125_x);
            require_unsat("!(x!=0) outside epsilon",!(ex!=0),0.25_x);

            require_delta_overlap("!(x==0) equality delta-overlap",!(ex==0),0.0_x);
            require_epsilon_sat("!(x==0) original true",!(ex==0),0.25_x);
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
        {
            std::cout << "[smt-parallel] empty domains are UNSAT" << std::endl;
            ExactBoxType empty_domain({ExactIntervalType::empty_interval()});
            SmtSolver solver(SmtSolverConfiguration(0.125_x));

            List<ValidatedConstraint> constraints;
            SmtResult constraints_result=solver.solve_parallel(
                empty_domain,constraints);
            ARIADNE_TEST_ASSERT(constraints_result.is_unsat());

            RealVariable x("empty_parallel_x");
            RealSpace space({x});
            List<SmtTheoryPrimitiveLiteral> literals;
            SmtResult theory_result=solver.solve_parallel(
                space,empty_domain,literals);
            ARIADNE_TEST_ASSERT(theory_result.is_unsat());
        }

        {
            std::cout << "[smt-parallel] empty conjunctions ignore zero box budget" << std::endl;
            SmtSolver zero_budget_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            ExactBoxType domain({ExactIntervalType(-1,1)});

            List<ValidatedConstraint> constraints;
            SmtResult constraints_result=zero_budget_solver.solve_parallel(
                domain,constraints);
            ARIADNE_TEST_ASSERT(constraints_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(constraints_result.has_witness());
            ARIADNE_TEST_EQUAL(constraints_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(
                constraints_result.statistics().box_budget_exhaustions,0u);

            RealVariable x("x");
            RealSpace space({x});
            List<SmtTheoryPrimitiveLiteral> literals;
            SmtResult theory_result=zero_budget_solver.solve_parallel(
                space,domain,literals);
            ARIADNE_TEST_ASSERT(theory_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(theory_result.has_witness());
            ARIADNE_TEST_EQUAL(theory_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(
                theory_result.statistics().box_budget_exhaustions,0u);
        }

        {
            std::cout << "[smt-parallel] zero box budget returns UNKNOWN" << std::endl;
            auto x=ValidatedScalarMultivariateFunction::coordinates(1);
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            ExactBoxType domain({ExactIntervalType(0,1)});
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(ValidatedNumber(0.5_x),x[0],ValidatedNumber(0.5_x))
            });
            SmtResult solve_result=bounded_solver.solve_parallel(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_uncertified_boxes,0u);
        }

        {
            std::cout << "[smt-parallel] non-splittable epsilon-overlap counted once" << std::endl;
            auto x=ValidatedScalarMultivariateFunction::coordinates(1);
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(
                1e-30_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            ExactBoxType domain({ExactIntervalType(1,1)});
            ValidatedScalarMultivariateFunction residual=sin(x[0])-sin(x[0]);
            List<ValidatedConstraint> constraints({
                ValidatedConstraint(
                    ValidatedNumber(0),
                    residual,
                    ValidatedNumber(0))
            });
            SmtResult solve_result=tiny_epsilon_solver.solve_parallel(domain,constraints);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_uncertified_boxes,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_epsilon_overlap_boxes,1u);
        }

        {
            std::cout << "[smt-parallel] theory epsilon-overlap counted once" << std::endl;
            RealVariable x("x");
            RealExpression ex=x;
            RealSpace space({x});
            SmtSolver tiny_epsilon_solver(SmtSolverConfiguration(
                1e-30_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                false));
            RealExpression residual=sqr(sin(ex))+sqr(cos(ex))-1;
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(residual==0));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=tiny_epsilon_solver.solve_parallel(
                space,
                ExactBoxType({ExactIntervalType(1,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_uncertified_boxes,1u);
            ARIADNE_TEST_EQUAL(
                solve_result.statistics().non_splittable_epsilon_overlap_boxes,1u);
        }

        {
            std::cout << "[smt-parallel] simplified zero primitive ignores zero box budget" << std::endl;
            RealVariable x("zero_parallel_x");
            RealExpression ex=x;
            RealSpace space({x});
            SmtSolver zero_budget_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u,
                false));
            RealExpression residual=sin(ex)-sin(ex);
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(residual>=0));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=zero_budget_solver.solve_parallel(
                space,
                ExactBoxType({ExactIntervalType(1,1)}),
                literals);
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,0u);
        }

        {
            std::cout << "[smt-parallel] theory zero box budget returns UNKNOWN" << std::endl;
            RealVariable x("x");
            RealSpace space({x});
            SmtSolver bounded_solver(SmtSolverConfiguration(
                0.125_x,
                std::numeric_limits<SizeType>::max(),
                std::numeric_limits<SizeType>::max(),
                0u));
            auto alternatives=normalize_smt_theory_literal(
                make_smt_theory_literal(RealExpression(x)==0.5_x));
            ARIADNE_TEST_EQUAL(alternatives.size(),1u);
            ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
            List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
            SmtResult solve_result=bounded_solver.solve_parallel(
                space,ExactBoxType({ExactIntervalType(0,1)}),literals);
            ARIADNE_TEST_ASSERT(solve_result.is_unknown());
            ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,0u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
            ARIADNE_TEST_EQUAL(solve_result.statistics().non_splittable_uncertified_boxes,0u);
        }

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
            SmtSolverTestSupport::begin_parallel_execution_observation();
            SmtResult solve_result=solver.solve_parallel(domain,constraints);
            auto observation=SmtSolverTestSupport::end_parallel_execution_observation();
            ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
            ARIADNE_TEST_ASSERT(solve_result.has_witness());
            ARIADNE_TEST_ASSERT(observation.calling_thread_observed);
            ARIADNE_TEST_EQUAL(observation.worker_thread_count,0u);
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
                SmtSolverTestSupport::begin_parallel_execution_observation();
                SmtResult solve_result=solver.solve_parallel(domain,constraints);
                auto observation=SmtSolverTestSupport::end_parallel_execution_observation();
                ARIADNE_TEST_ASSERT(solve_result.is_epsilon_sat());
                ARIADNE_TEST_ASSERT(solve_result.has_witness());
                ARIADNE_TEST_ASSERT(not observation.calling_thread_observed);
                ARIADNE_TEST_ASSERT(observation.worker_thread_count>=1u);
                std::cout << "[smt-parallel] EPSILON_SAT processed="
                          << solve_result.statistics().boxes_processed
                          << " pruned=" << solve_result.statistics().boxes_pruned
                          << " split=" << solve_result.statistics().boxes_split << std::endl;
            }

            {
                std::cout << "[smt-parallel] concurrent split with one-box budget" << std::endl;
                auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
                SmtSolver split_solver(SmtSolverConfiguration(
                    0.01_x,
                    std::numeric_limits<SizeType>::max(),
                    std::numeric_limits<SizeType>::max(),
                    1u,
                    false));
                ExactBoxType domain({
                    ExactIntervalType(-1,1),
                    ExactIntervalType(-10000,10000)
                });
                List<ValidatedConstraint> constraints({
                    ValidatedConstraint(
                        ValidatedNumber(0.3_x),
                        sin(10*xy[0]),
                        ValidatedNumber(0.3_x))
                });
                SmtResult solve_result=split_solver.solve_parallel(domain,constraints);
                ARIADNE_TEST_ASSERT(solve_result.is_unknown());
                ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().sensitivity_guided_splits,1u);
                ARIADNE_TEST_EQUAL(
                    solve_result.statistics().sensitivity_overrides_geometric_splits,1u);
            }

            {
                std::cout << "[smt-parallel] concurrent theory split with one-box budget" << std::endl;
                RealVariable tx("parallel_theory_x"), ty("parallel_theory_y");
                RealExpression etx=tx;
                RealSpace theory_space({tx,ty});
                SmtSolver split_solver(SmtSolverConfiguration(
                    0.01_x,
                    std::numeric_limits<SizeType>::max(),
                    std::numeric_limits<SizeType>::max(),
                    1u,
                    false));
                auto alternatives=normalize_smt_theory_literal(
                    make_smt_theory_literal(sin(10*etx)==0.3_x));
                ARIADNE_TEST_EQUAL(alternatives.size(),1u);
                ARIADNE_TEST_EQUAL(alternatives[0].size(),1u);
                List<SmtTheoryPrimitiveLiteral> literals({alternatives[0][0]});
                SmtResult solve_result=split_solver.solve_parallel(
                    theory_space,
                    ExactBoxType({
                        ExactIntervalType(-1,1),
                        ExactIntervalType(-10000,10000)
                    }),
                    literals);
                ARIADNE_TEST_ASSERT(solve_result.is_unknown());
                ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_processed,1u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().boxes_split,1u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().box_budget_exhaustions,1u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().candidate_witness_searches,0u);
                ARIADNE_TEST_EQUAL(solve_result.statistics().sensitivity_guided_splits,1u);
                ARIADNE_TEST_EQUAL(
                    solve_result.statistics().sensitivity_overrides_geometric_splits,1u);
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
                ARIADNE_TEST_ASSERT(solve_result.statistics().boxes_pruned>0u);
                ARIADNE_TEST_ASSERT(solve_result.statistics().hull_reduction_rounds>=1u);
                ARIADNE_TEST_ASSERT(solve_result.statistics().shaving_reduction_rounds>=1u);
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
