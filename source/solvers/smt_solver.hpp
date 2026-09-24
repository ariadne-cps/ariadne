/***************************************************************************
 *            solvers/smt_solver.hpp
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

/*! \file solvers/smt_solver.hpp
 *  \brief Common types for solving bounded real SMT problems.
 */

#ifndef ARIADNE_SMT_SOLVER_HPP
#define ARIADNE_SMT_SOLVER_HPP

#include <optional>
#include <vector>
#include <mutex>
#include <limits>

#include "geometry/box.hpp"
#include "numeric/numeric.hpp"
#include "function/constraint.hpp"
#include "solvers/smt_theory.hpp"
#include "solvers/smt_boolean.hpp"
#include "symbolic/space.hpp"

namespace Ariadne {

//! \ingroup Solvers
//! \brief Logical status returned by a bounded real epsilon-SMT solver.
enum class SmtResultStatus {
    UNSAT,
    EPSILON_SAT,
    UNKNOWN
};

//! \ingroup Solvers
//! \brief Configuration shared by bounded real epsilon-SMT solvers.
class SmtSolverConfiguration {
  public:
    explicit SmtSolverConfiguration(
        ExactDouble epsilon,
        SizeType theory_minimization_budget=std::numeric_limits<SizeType>::max(),
        SizeType learned_clause_limit=std::numeric_limits<SizeType>::max(),
        SizeType box_processing_limit=std::numeric_limits<SizeType>::max(),
        Bool candidate_search_enabled=true);

    //! \brief The logical epsilon used for weakening constraints.
    ExactDouble epsilon() const { return _epsilon; }

    //! \brief Maximum theory checks used to minimize one learned theory nogood.
    SizeType theory_minimization_budget() const { return _theory_minimization_budget; }

    //! \brief Maximum number of active non-theory learned clauses.
    SizeType learned_clause_limit() const { return _learned_clause_limit; }

    //! \brief Maximum number of search boxes processed by one theory solve.
    //! \details Reaching the limit yields UNKNOWN unless an epsilon witness was found first.
    SizeType box_processing_limit() const { return _box_processing_limit; }

    //! \brief Whether heuristic interior-point witness candidate generation is enabled.
    Bool candidate_search_enabled() const { return _candidate_search_enabled; }

  private:
    ExactDouble _epsilon;
    SizeType _theory_minimization_budget;
    SizeType _learned_clause_limit;
    SizeType _box_processing_limit;
    Bool _candidate_search_enabled;
};

//! \ingroup Solvers
//! \brief Search statistics for a bounded real epsilon-SMT query.
struct SmtSearchStatistics {
    SizeType boxes_processed = 0u;
    SizeType boxes_pruned = 0u;
    SizeType boxes_split = 0u;
    SizeType boxes_unknown = 0u;
    SizeType box_budget_exhaustions = 0u;
    SizeType non_splittable_uncertified_boxes = 0u;
    SizeType non_splittable_epsilon_overlap_boxes = 0u;
    SizeType hull_reduction_rounds = 0u;
    SizeType hull_effective_reductions = 0u;
    SizeType shaving_reduction_rounds = 0u;
    SizeType shaving_effective_reductions = 0u;
    SizeType sensitivity_guided_splits = 0u;
    SizeType sensitivity_overrides_geometric_splits = 0u;
    SizeType epsilon_box_certifications = 0u;
    SizeType candidate_witness_searches = 0u;
    SizeType candidate_witness_successes = 0u;
    SizeType boolean_decisions = 0u;
    SizeType boolean_propagations = 0u;
    SizeType boolean_reasoned_propagations = 0u;
    SizeType boolean_conflicts = 0u;
    SizeType boolean_backtracks = 0u;
    SizeType max_decision_level = 0u;
    SizeType boolean_conflicts_analyzed = 0u;
    SizeType learned_clause_literals = 0u;
    SizeType last_learned_clause_literals = 0u;
    SizeType last_learned_current_level_literals = 0u;
    SizeType last_backjump_level = 0u;
    SizeType learned_clauses = 0u;
    SizeType learned_clause_propagations = 0u;
    SizeType nonchronological_backjumps = 0u;
    SizeType theory_checks = 0u;
    SizeType theory_conflicts = 0u;
    SizeType theory_learned_clauses = 0u;
    SizeType theory_learned_clause_literals = 0u;
    SizeType theory_learned_clause_propagations = 0u;
    SizeType theory_minimization_checks = 0u;
    SizeType theory_nogood_raw_literals = 0u;
    SizeType theory_nogood_minimized_literals = 0u;
    SizeType theory_nogood_literals_removed = 0u;
    SizeType theory_minimization_budget_exhaustions = 0u;
    SizeType first_minimization_candidate_trail_rank = 0u;
    SizeType learned_clause_activity_bumps = 0u;
    SizeType learned_clause_pruning_runs = 0u;
    SizeType learned_clauses_pruned = 0u;
    SizeType peak_active_non_theory_learned_clauses = 0u;
};

//! \ingroup Solvers
//! \brief Result of a bounded real epsilon-SMT query.
class SmtResult {
  public:
    //! \brief Construct an UNSAT result.
    static SmtResult unsat(SmtSearchStatistics statistics = {});

    //! \brief Construct an EPSILON_SAT result with a validated witness box.
    static SmtResult epsilon_sat(UpperBoxType const& witness,
                                 SmtSearchStatistics statistics = {});

    //! \brief Construct an inconclusive result.
    static SmtResult unknown(SmtSearchStatistics statistics = {});

    SmtResultStatus status() const { return _status; }
    Bool is_unsat() const { return _status==SmtResultStatus::UNSAT; }
    Bool is_epsilon_sat() const { return _status==SmtResultStatus::EPSILON_SAT; }
    Bool is_unknown() const { return _status==SmtResultStatus::UNKNOWN; }

    //! \brief True iff the result contains a witness box.
    Bool has_witness() const { return _witness.has_value(); }

    //! \brief Return the witness box.
    //! \pre The result is EPSILON_SAT.
    UpperBoxType const& witness() const;

    SmtSearchStatistics const& statistics() const { return _statistics; }

  private:
    SmtResult(SmtResultStatus status, SmtSearchStatistics statistics);
    SmtResult(SmtResultStatus status, UpperBoxType const& witness,
              SmtSearchStatistics statistics);

    SmtResultStatus _status;
    std::optional<UpperBoxType> _witness;
    SmtSearchStatistics _statistics;
};

OutputStream& operator<<(OutputStream& os, SmtResultStatus status);

class SmtSolver;

namespace SmtSolverTestSupport {

struct LearnedClausePruningEntry {
    Bool active = true;
    Bool theory = false;
    Bool recent = false;
    Bool short_clause = false;
    Bool useful = false;
    Bool protected_clause = false;
    Bool locked = false;
    SizeType activity = 1u;
    SizeType size = 0u;
};

std::vector<SizeType> learned_clause_pruning_candidates(
    std::vector<LearnedClausePruningEntry> const& entries);

SizeType apply_learned_clause_pruning(
    std::vector<SizeType> const& candidates,
    SizeType active_count,
    SizeType limit);

Void accumulate_statistics(
    SmtSearchStatistics& target,
    SmtSearchStatistics const& source);

Void record_first_minimization_candidate_trail_rank(
    SmtSearchStatistics& statistics,
    SizeType trail_rank);

struct ParallelExecutionObservation {
    SizeType observed_thread_count = 0u;
    SizeType worker_thread_count = 0u;
    Bool calling_thread_observed = false;
};

Void begin_parallel_execution_observation();
ParallelExecutionObservation end_parallel_execution_observation();
Void record_parallel_processing_thread();

Bool parallel_stop_condition(Bool found, Bool limit_reached);
std::vector<UpperBoxType> parallel_children_to_append(
    Bool found,
    Pair<UpperBoxType,UpperBoxType> const& children);
Pair<Bool,Bool> parallel_witness_claim_sequence();

struct CandidateWitnessOutcome {
    Bool attempted = false;
    Bool certified = false;
    std::optional<UpperBoxType> witness;
};

CandidateWitnessOutcome candidate_witness_outcome(
    Bool attempted,
    std::optional<UpperBoxType> const& certified_witness);

SizeType epsilon_witness_candidate_count(UpperBoxType const& domain);

struct SensitivitySplitSelection {
    SizeType coordinate = 0u;
    Bool guided = false;
    Bool overrode_geometric = false;
};

SensitivitySplitSelection sensitivity_split_selection(
    UpperBoxType const& domain,
    std::vector<ValidatedScalarMultivariateFunction> const& functions);

struct SearchOutcome {
    std::optional<UpperBoxType> witness;
    std::optional<SizeType> backjump_level;

    static SearchOutcome exhausted();
    static SearchOutcome found(UpperBoxType const& witness);
    static SearchOutcome backjump(SizeType level);
};

SmtResult finalize_search_outcome(
    SearchOutcome const& outcome,
    Bool theory_unknown_seen,
    SmtSearchStatistics const& statistics);

enum class BoxProcessingStatus {
    PRUNED,
    EPSILON_SAT,
    SPLIT,
    UNKNOWN
};

struct BoxProcessingStatisticsInput {
    BoxProcessingStatus status;
    SizeType hull_rounds = 0u;
    SizeType hull_effective = 0u;
    SizeType shaving_rounds = 0u;
    SizeType shaving_effective = 0u;
    Bool sensitivity_guided_split = false;
    Bool sensitivity_overrode_geometric_split = false;
    Bool epsilon_box_certification = false;
    Bool candidate_witness_search = false;
    Bool candidate_witness_success = false;
};

Void accumulate_box_processing_statistics(
    SmtSearchStatistics& statistics,
    BoxProcessingStatisticsInput const& input);

Void validate_primitive_relation(SmtTheoryPrimitiveRelation relation);

std::vector<Int> resolve_clause_on_variable(
    std::vector<Int> const& lhs,
    std::vector<Int> const& rhs,
    SizeType variable);

Void order_theory_nogood(
    std::vector<Int>& clause,
    std::vector<SizeType> const& decision_levels,
    std::vector<SizeType> const& trail_rank);

Bool clause_is_learned(SizeType index, SizeType original_clause_count);
Bool learned_clause_is_theory(
    SizeType index,
    SizeType original_clause_count,
    std::vector<Bool> const& theory_flags);
Bool assignment_locks_clause(
    int8_t assignment_value,
    std::optional<SizeType> const& reason_clause,
    SizeType clause_index);

enum class TheoryAtomTruth {
    FALSE_VALUE,
    UNKNOWN,
    TRUE_VALUE
};

TheoryAtomTruth classify_theory_relation(
    SmtTheoryRelation relation,
    UpperIntervalType const& image);

TheoryAtomTruth classify_theory_atom(
    RealSpace const& space,
    ExactBoxType const& domain,
    ContinuousPredicate const& atom);

struct TheoryResultInterpretation {
    Bool consistent = false;
    Bool unknown = false;
    std::optional<UpperBoxType> witness;
};

TheoryResultInterpretation interpret_theory_result(SmtResult const& result);

ExactIntervalType original_bounds(
    SmtSolver const& solver,
    SmtTheoryPrimitiveRelation relation);

ExactIntervalType epsilon_bounds(
    SmtSolver const& solver,
    SmtTheoryPrimitiveRelation relation);

Bool epsilon_satisfied(
    SmtSolver const& solver,
    UpperBoxType const& domain,
    List<ValidatedConstraint> const& constraints);

Bool epsilon_satisfied(
    SmtSolver const& solver,
    RealSpace const& space,
    UpperBoxType const& domain,
    List<SmtTheoryPrimitiveLiteral> const& literals);

} // namespace SmtSolverTestSupport

//! \ingroup Solvers
//! \brief Sequential reference epsilon-SMT solver for bounded conjunctions.
class SmtSolver {
  public:
    explicit SmtSolver(SmtSolverConfiguration configuration)
        : _configuration(configuration) { }

    //! \brief Solve a bounded conjunction of validated real constraints.
    SmtResult solve(ExactBoxType const& domain,
                    List<ValidatedConstraint> const& constraints) const;

    //! \brief Solve using BetterThreads dynamic workload processing.
    //! \details The actual concurrency is controlled by BetterThreads::ThreadManager.
    SmtResult solve_parallel(ExactBoxType const& domain,
                             List<ValidatedConstraint> const& constraints) const;

    //! \brief Solve a conjunction of normalized real theory primitives.
    SmtResult solve(RealSpace const& space,
                    ExactBoxType const& domain,
                    List<SmtTheoryPrimitiveLiteral> const& literals) const;

    //! \brief Solve normalized theory primitives using BetterThreads.
    SmtResult solve_parallel(RealSpace const& space,
                             ExactBoxType const& domain,
                             List<SmtTheoryPrimitiveLiteral> const& literals) const;

    //! \brief Solve a bounded Boolean combination of real predicates.
    SmtResult solve(RealSpace const& space,
                    ExactBoxType const& domain,
                    ContinuousPredicate const& predicate) const;

    //! \brief Solve a bounded Boolean combination using parallel theory search.
    SmtResult solve_parallel(RealSpace const& space,
                             ExactBoxType const& domain,
                             ContinuousPredicate const& predicate) const;

    SmtSolverConfiguration const& configuration() const { return _configuration; }

  private:
    friend struct SmtParallelTask;
    friend ExactIntervalType SmtSolverTestSupport::original_bounds(
        SmtSolver const&, SmtTheoryPrimitiveRelation);
    friend ExactIntervalType SmtSolverTestSupport::epsilon_bounds(
        SmtSolver const&, SmtTheoryPrimitiveRelation);
    friend Bool SmtSolverTestSupport::epsilon_satisfied(
        SmtSolver const&, UpperBoxType const&, List<ValidatedConstraint> const&);
    friend Bool SmtSolverTestSupport::epsilon_satisfied(
        SmtSolver const&, RealSpace const&, UpperBoxType const&,
        List<SmtTheoryPrimitiveLiteral> const&);
    using BoxProcessingStatus=SmtSolverTestSupport::BoxProcessingStatus;

    struct ReductionStatistics {
        SizeType hull_rounds = 0u;
        SizeType hull_effective = 0u;
        SizeType shaving_rounds = 0u;
        SizeType shaving_effective = 0u;
    };

    struct BoxProcessingResult {
        BoxProcessingStatus status;
        std::optional<UpperBoxType> witness;
        std::optional<Pair<UpperBoxType,UpperBoxType>> children;
        ReductionStatistics reductions;
        Bool sensitivity_guided_split = false;
        Bool sensitivity_overrode_geometric_split = false;
        Bool epsilon_box_certification = false;
        Bool candidate_witness_search = false;
        Bool candidate_witness_success = false;
        Bool non_splittable_epsilon_overlap = false;
    };

    struct CompiledTheoryLiteral {
        ValidatedScalarMultivariateFunction function;
        SmtTheoryPrimitiveRelation relation;
    };
    using CompiledTheoryLiterals = std::vector<CompiledTheoryLiteral>;

    struct ConjunctionReference {
        List<ValidatedConstraint> const* constraints;
        CompiledTheoryLiterals const* theory_literals;

        explicit ConjunctionReference(List<ValidatedConstraint> const& conjunction)
            : constraints(&conjunction), theory_literals(nullptr) { }

        explicit ConjunctionReference(CompiledTheoryLiterals const& conjunction)
            : constraints(nullptr), theory_literals(&conjunction) { }
    };

    template<class Conjunction>
    BoxProcessingResult _process_box(
        UpperBoxType domain,
        Conjunction const& conjunction) const;

    BoxProcessingResult _process_box(
        UpperBoxType domain,
        ConjunctionReference const& conjunction) const;

    Void _accumulate_box_processing_statistics(
        SmtSearchStatistics& statistics,
        BoxProcessingResult const& processing) const;

    SmtResult _solve_sequential_conjunction(
        ExactBoxType const& domain,
        ConjunctionReference const& conjunction) const;

    SmtResult _solve_parallel_conjunction(
        ExactBoxType const& domain,
        ConjunctionReference const& conjunction) const;

    ExactIntervalType _original_bounds(ValidatedConstraint const& constraint) const;
    ExactIntervalType _original_bounds(SmtTheoryPrimitiveRelation relation) const;
    ExactIntervalType _epsilon_bounds(ValidatedConstraint const& constraint) const;
    ExactIntervalType _epsilon_bounds(SmtTheoryPrimitiveRelation relation) const;
    ExactIntervalType _epsilon_bounds(CompiledTheoryLiteral const& literal) const;
    Bool _original_reduce(UpperBoxType& domain,
                          List<ValidatedConstraint> const& constraints,
                          ReductionStatistics& statistics) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            List<ValidatedConstraint> const& constraints) const;

    CompiledTheoryLiterals _compile_theory_literals(
        RealSpace const& space,
        List<SmtTheoryPrimitiveLiteral> const& literals) const;
    Bool _original_reduce(UpperBoxType& domain,
                          CompiledTheoryLiterals const& literals,
                          ReductionStatistics& statistics) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            CompiledTheoryLiterals const& literals) const;

    template<class Conjunction>
    std::optional<UpperBoxType> _epsilon_witness(
        UpperBoxType const& domain,
        Conjunction const& conjunction) const;
    template<class Conjunction>
    UpperBoxType _epsilon_candidate_witness(
        UpperBoxType const& domain,
        Conjunction const& conjunction) const;
    template<class Conjunction>
    Pair<Pair<UpperBoxType,UpperBoxType>,Pair<Bool,Bool>> _split_box(
        UpperBoxType const& domain,
        Conjunction const& conjunction) const;

    ValidatedScalarMultivariateFunction const& _function(
        ValidatedConstraint const& constraint) const;
    ValidatedScalarMultivariateFunction const& _function(
        CompiledTheoryLiteral const& literal) const;


    SmtSolverConfiguration _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_SMT_SOLVER_HPP
