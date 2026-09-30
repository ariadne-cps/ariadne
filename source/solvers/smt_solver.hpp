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
#include "solvers/constraint_solver.hpp"
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
//! \brief Reason why an SMT search returned UNKNOWN.
enum class SmtUnknownReason {
    NONE,
    RESOURCE_EXHAUSTED,
    DP_RESOLUTION_EXHAUSTED,
    MIXED
};

//! \ingroup Solvers
//! \brief Outcome of a cheap preclassification attempt before contractors.
enum class SmtPreclassificationOutcome {
    PRUNED,
    EPSILON_SAT,
    UNRESOLVED
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
        Bool candidate_search_enabled=true,
        Bool monotone_reduction_enabled=false,
        Bool sensitivity_split_enabled=true,
        Bool deterministic_witness_probing_enabled=true,
        Bool shaving_reduction_enabled=true,
        Bool hull_reduction_enabled=true,
        Bool interval_lookahead_split_enabled=false,
        Bool upper_child_first=false,
        Bool interval_newton_reduction_enabled=false,
        Bool preclassification_enabled=false,
        Bool adaptive_preclassification_enabled=false);

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

    //! \brief Whether validated monotone/Newton contraction is enabled.
    Bool monotone_reduction_enabled() const { return _monotone_reduction_enabled; }

    //! \brief Whether sensitivity-guided splitting is enabled.
    Bool sensitivity_split_enabled() const { return _sensitivity_split_enabled; }

    //! \brief Whether interval lookahead splitting is enabled.
    Bool interval_lookahead_split_enabled() const {
        return _interval_lookahead_split_enabled;
    }

    //! \brief Whether sequential DFS visits the upper split child first.
    Bool upper_child_first() const { return _upper_child_first; }

    //! \brief Whether deterministic midpoint/endpoint/corner witness probing is enabled.
    Bool deterministic_witness_probing_enabled() const {
        return _deterministic_witness_probing_enabled;
    }

    //! \brief Whether coordinate shaving is enabled during validated propagation.
    Bool shaving_reduction_enabled() const { return _shaving_reduction_enabled; }

    //! \brief Whether hull reduction is enabled during validated propagation.
    Bool hull_reduction_enabled() const { return _hull_reduction_enabled; }

    //! \brief Whether square EQ_ZERO conjunctions use Interval Newton contraction.
    Bool interval_newton_reduction_enabled() const {
        return _interval_newton_reduction_enabled;
    }

    //! \brief Whether compiled theory boxes are cheaply classified before contractors.
    Bool preclassification_enabled() const {
        return _preclassification_enabled;
    }

    //! \brief Whether sequential theory search adaptively suspends unproductive preclassification.
    Bool adaptive_preclassification_enabled() const {
        return _adaptive_preclassification_enabled;
    }

  private:
    ExactDouble _epsilon;
    SizeType _theory_minimization_budget;
    SizeType _learned_clause_limit;
    SizeType _box_processing_limit;
    Bool _candidate_search_enabled;
    Bool _monotone_reduction_enabled;
    Bool _sensitivity_split_enabled;
    Bool _interval_lookahead_split_enabled;
    Bool _upper_child_first;
    Bool _deterministic_witness_probing_enabled;
    Bool _shaving_reduction_enabled;
    Bool _hull_reduction_enabled;
    Bool _interval_newton_reduction_enabled;
    Bool _preclassification_enabled;
    Bool _adaptive_preclassification_enabled;
};

//! \ingroup Solvers
//! \brief Search statistics for a bounded real epsilon-SMT query.
struct SmtSearchStatistics {
    SizeType boxes_processed = 0u;
    SizeType boxes_pruned = 0u;
    SizeType boxes_split = 0u;
    SizeType boxes_unknown = 0u;
    SizeType box_budget_exhaustions = 0u;
    SizeType dp_resolution_exhaustions = 0u;
    SizeType non_splittable_uncertified_boxes = 0u;
    SizeType non_splittable_epsilon_overlap_boxes = 0u;
    SizeType hull_reduction_rounds = 0u;
    SizeType hull_effective_reductions = 0u;
    SizeType hull_procedure_builds = 0u;
    double hull_procedure_build_seconds = 0.0;
    double hull_contraction_seconds = 0.0;
    double hull_direct_rejection_seconds = 0.0;
    double hull_temporary_allocation_seconds = 0.0;
    double hull_forward_execution_seconds = 0.0;
    double hull_backward_propagation_seconds = 0.0;
    SizeType interval_newton_attempts = 0u;
    SizeType interval_newton_effective_reductions = 0u;
    SizeType interval_newton_infeasible = 0u;
    SizeType interval_newton_singular = 0u;
    double interval_newton_seconds = 0.0;
    SizeType shaving_reduction_rounds = 0u;
    SizeType shaving_effective_reductions = 0u;
    SizeType shaving_coordinate_attempts = 0u;
    SizeType shaving_coordinate_effective = 0u;
    SizeType shaving_dependency_skipped = 0u;
    SizeType shaving_adaptive_skipped = 0u;
    SizeType shaving_refresh_rounds = 0u;
    SizeType shaving_active_rounds = 0u;
    SizeType shaving_function_evaluations = 0u;
    double shaving_seconds = 0.0;
    SizeType monotone_reduction_rounds = 0u;
    SizeType monotone_effective_reductions = 0u;
    SizeType sensitivity_guided_splits = 0u;
    SizeType sensitivity_overrides_geometric_splits = 0u;
    SizeType sensitivity_derivatives_built = 0u;
    SizeType sensitivity_derivative_evaluations = 0u;
    SizeType interval_lookahead_guided_splits = 0u;
    SizeType interval_lookahead_overrides_geometric_splits = 0u;
    SizeType interval_lookahead_function_evaluations = 0u;
    SizeType epsilon_box_certifications = 0u;
    SizeType fused_direct_classification_boxes = 0u;
    SizeType fused_direct_literal_evaluations = 0u;
    SizeType preclassification_boxes = 0u;
    SizeType preclassification_pruned_boxes = 0u;
    SizeType preclassification_epsilon_boxes = 0u;
    SizeType preclassification_literal_evaluations = 0u;
    std::vector<SmtPreclassificationOutcome> preclassification_outcomes;
    SizeType preclassification_adaptive_skipped_boxes = 0u;
    SizeType preclassification_adaptive_suspensions = 0u;
    SizeType preclassification_adaptive_reactivations = 0u;
    SizeType candidate_witness_searches = 0u;
    SizeType candidate_witness_successes = 0u;
    double theory_compile_seconds = 0.0;
    double reduction_seconds = 0.0;
    double epsilon_check_seconds = 0.0;
    double witness_probe_seconds = 0.0;
    double split_seconds = 0.0;
    double sensitivity_derivative_build_seconds = 0.0;
    double sensitivity_derivative_evaluation_seconds = 0.0;
    double interval_lookahead_evaluation_seconds = 0.0;
    double preclassification_seconds = 0.0;
    double candidate_search_seconds = 0.0;
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

    //! \brief Construct an inconclusive result with an explicit cause.
    static SmtResult unknown(
        SmtUnknownReason reason,
        SmtSearchStatistics statistics = {});

    SmtResultStatus status() const { return _status; }
    Bool is_unsat() const { return _status==SmtResultStatus::UNSAT; }
    Bool is_epsilon_sat() const { return _status==SmtResultStatus::EPSILON_SAT; }
    Bool is_unknown() const { return _status==SmtResultStatus::UNKNOWN; }

    //! \brief Cause of UNKNOWN; NONE for conclusive results.
    SmtUnknownReason unknown_reason() const { return _unknown_reason; }

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
    SmtUnknownReason _unknown_reason = SmtUnknownReason::NONE;
    std::optional<UpperBoxType> _witness;
    SmtSearchStatistics _statistics;
};

OutputStream& operator<<(OutputStream& os, SmtResultStatus status);
OutputStream& operator<<(OutputStream& os, SmtUnknownReason reason);

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

struct PreclassificationShadowSummary {
    SizeType checks = 0u;
    SizeType skipped = 0u;
    SizeType observed_hits = 0u;
    SizeType skipped_hits = 0u;
    SizeType suspensions = 0u;
    SizeType reactivations = 0u;
};

PreclassificationShadowSummary preclassification_shadow_summary(
    std::vector<SmtPreclassificationOutcome> const& outcomes,
    SizeType initial_window=16u,
    SizeType active_window=8u,
    SizeType refresh_period=8u);

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
std::vector<UpperBoxType> sequential_children_to_push(
    Bool upper_child_first,
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

Void validate_differentiable_expression_kind(OperatorKind kind);

Bool expression_is_differentiable(RealExpression const& expression);

std::optional<ValidatedScalarMultivariateFunction> optional_derivative(
    RealExpression const& expression,
    ValidatedScalarMultivariateFunction const& function,
    SizeType variable);

//! \brief Number of cached coordinate derivatives produced while compiling
//! theory literals under the solver's current configuration.
SizeType compiled_theory_derivative_count(
    SmtSolver const& solver,
    RealSpace const& space,
    List<SmtTheoryPrimitiveLiteral> const& literals);

struct SearchOutcome {
    std::optional<UpperBoxType> witness;
    std::optional<SizeType> backjump_level;

    static SearchOutcome exhausted();
    static SearchOutcome found(UpperBoxType const& witness);
    static SearchOutcome backjump(SizeType level);
};

SmtUnknownReason combine_unknown_reasons(
    SmtUnknownReason first,
    SmtUnknownReason second);

SmtResult finalize_search_outcome(
    SearchOutcome const& outcome,
    SmtUnknownReason theory_unknown_reason,
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
    SizeType monotone_rounds = 0u;
    SizeType monotone_effective = 0u;
    Bool sensitivity_guided_split = false;
    Bool sensitivity_overrode_geometric_split = false;
    Bool epsilon_box_certification = false;
    Bool dp_resolution_exhausted = false;
    Bool candidate_witness_search = false;
    Bool candidate_witness_success = false;
};

Void accumulate_box_processing_statistics(
    SmtSearchStatistics& statistics,
    BoxProcessingStatisticsInput const& input);

struct DirectBoxProcessingObservation {
    BoxProcessingStatus status = BoxProcessingStatus::UNKNOWN;
    Bool fused_direct_classification = false;
    Bool epsilon_box_certification = false;
};

DirectBoxProcessingObservation process_compiled_box(
    SmtSolver const& solver,
    RealSpace const& space,
    UpperBoxType const& domain,
    List<SmtTheoryPrimitiveLiteral> const& literals);

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

struct TheoryAtomImplication {
    Bool force_true = false;
    Bool force_false = false;
};

Bool epsilon_primitive_image_infeasible(
    SmtTheoryPrimitiveRelation relation,
    UpperIntervalType const& image,
    FloatDP const& epsilon);

Bool epsilon_theory_literal_infeasible(
    SmtSolver const& solver,
    RealSpace const& space,
    ExactBoxType const& domain,
    SmtTheoryLiteral const& literal);

TheoryAtomImplication domain_theory_implication(
    SmtSolver const& solver,
    RealSpace const& space,
    ExactBoxType const& domain,
    ContinuousPredicate const& atom);

struct TheoryResultInterpretation {
    Bool consistent = false;
    SmtUnknownReason unknown_reason = SmtUnknownReason::NONE;
    std::optional<UpperBoxType> witness;
};

TheoryResultInterpretation interpret_theory_result(SmtResult const& result);

ExactIntervalType original_bounds(
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
    friend Bool SmtSolverTestSupport::epsilon_satisfied(
        SmtSolver const&, UpperBoxType const&, List<ValidatedConstraint> const&);
    friend Bool SmtSolverTestSupport::epsilon_satisfied(
        SmtSolver const&, RealSpace const&, UpperBoxType const&,
        List<SmtTheoryPrimitiveLiteral> const&);
    friend SizeType SmtSolverTestSupport::compiled_theory_derivative_count(
        SmtSolver const&, RealSpace const&,
        List<SmtTheoryPrimitiveLiteral> const&);
    friend SmtSolverTestSupport::DirectBoxProcessingObservation
        SmtSolverTestSupport::process_compiled_box(
            SmtSolver const&, RealSpace const&, UpperBoxType const&,
            List<SmtTheoryPrimitiveLiteral> const&);
    using BoxProcessingStatus=SmtSolverTestSupport::BoxProcessingStatus;

    struct ReductionStatistics {
        SizeType hull_rounds = 0u;
        SizeType hull_effective = 0u;
        SizeType hull_procedure_builds = 0u;
        double hull_procedure_build_seconds = 0.0;
        double hull_contraction_seconds = 0.0;
        double hull_direct_rejection_seconds = 0.0;
        double hull_temporary_allocation_seconds = 0.0;
        double hull_forward_execution_seconds = 0.0;
        double hull_backward_propagation_seconds = 0.0;
        SizeType interval_newton_attempts = 0u;
        SizeType interval_newton_effective = 0u;
        SizeType interval_newton_infeasible = 0u;
        SizeType interval_newton_singular = 0u;
        double interval_newton_seconds = 0.0;
        SizeType shaving_rounds = 0u;
        SizeType shaving_effective = 0u;
        SizeType shaving_coordinate_attempts = 0u;
        SizeType shaving_coordinate_effective = 0u;
        SizeType shaving_dependency_skipped = 0u;
        SizeType shaving_adaptive_skipped = 0u;
        SizeType shaving_refresh_rounds = 0u;
        SizeType shaving_active_rounds = 0u;
        SizeType shaving_function_evaluations = 0u;
        double shaving_seconds = 0.0;
        SizeType monotone_rounds = 0u;
        SizeType monotone_effective = 0u;
    };

    struct BoxProcessingResult {
        BoxProcessingStatus status;
        std::optional<UpperBoxType> witness;
        std::optional<Pair<UpperBoxType,UpperBoxType>> children;
        ReductionStatistics reductions;
        Bool sensitivity_guided_split = false;
        Bool sensitivity_overrode_geometric_split = false;
        SizeType sensitivity_derivatives_built = 0u;
        SizeType sensitivity_derivative_evaluations = 0u;
        Bool interval_lookahead_guided_split = false;
        Bool interval_lookahead_overrode_geometric_split = false;
        SizeType interval_lookahead_function_evaluations = 0u;
        Bool epsilon_box_certification = false;
        Bool fused_direct_classification = false;
        SizeType fused_direct_literal_evaluations = 0u;
        Bool preclassification = false;
        Bool preclassification_pruned = false;
        Bool preclassification_epsilon_satisfied = false;
        SizeType preclassification_literal_evaluations = 0u;
        double preclassification_seconds = 0.0;
        Bool dp_resolution_exhausted = false;
        Bool candidate_witness_search = false;
        Bool candidate_witness_success = false;
        Bool non_splittable_epsilon_overlap = false;
        double reduction_seconds = 0.0;
        double epsilon_check_seconds = 0.0;
        double witness_probe_seconds = 0.0;
        double split_seconds = 0.0;
        double sensitivity_derivative_build_seconds = 0.0;
        double sensitivity_derivative_evaluation_seconds = 0.0;
        double interval_lookahead_evaluation_seconds = 0.0;
        double candidate_search_seconds = 0.0;
    };

    using CompiledTheoryLiteral = ConstraintPropagationConstraint;
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
        Conjunction const& conjunction,
        Bool allow_preclassification=true) const;

    BoxProcessingResult _process_box(
        UpperBoxType domain,
        ConjunctionReference const& conjunction,
        Bool allow_preclassification=true) const;

    Bool _preclassification_eligible(
        UpperBoxType const& domain,
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

    ExactIntervalType _original_bounds(SmtTheoryPrimitiveRelation relation) const;
    ExactIntervalType _epsilon_bounds(ValidatedConstraint const& constraint) const;
    ExactIntervalType _epsilon_bounds(CompiledTheoryLiteral const& literal) const;
    Bool _original_reduce(UpperBoxType& domain,
                          List<ValidatedConstraint> const& constraints,
                          ReductionStatistics& statistics) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            ValidatedConstraint const& constraint) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            List<ValidatedConstraint> const& constraints) const;

    CompiledTheoryLiterals _compile_theory_literals(
        RealSpace const& space,
        List<SmtTheoryPrimitiveLiteral> const& literals) const;
    Bool _original_reduce(UpperBoxType& domain,
                          CompiledTheoryLiterals const& literals,
                          ReductionStatistics& statistics) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            CompiledTheoryLiteral const& literal) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            CompiledTheoryLiterals const& literals) const;

    struct DirectClassification {
        Bool used = false;
        Bool preclassification = false;
        Bool pruned = false;
        Bool epsilon_satisfied = false;
        SizeType literal_evaluations = 0u;
        double seconds = 0.0;
    };

    DirectClassification _direct_classification(
        UpperBoxType const& domain,
        List<ValidatedConstraint> const& constraints,
        ReductionStatistics& statistics,
        Bool allow_preclassification) const;
    DirectClassification _direct_classification(
        UpperBoxType const& domain,
        CompiledTheoryLiterals const& literals,
        ReductionStatistics& statistics,
        Bool allow_preclassification) const;

    template<class Conjunction>
    std::optional<UpperBoxType> _epsilon_witness(
        UpperBoxType const& domain,
        Conjunction const& conjunction) const;

    template<class Conjunction>
    UpperBoxType _epsilon_candidate_witness(
        UpperBoxType const& domain,
        Conjunction const& conjunction) const;
    struct SplitBoxResult {
        Pair<UpperBoxType,UpperBoxType> children;
        Bool sensitivity_guided = false;
        Bool sensitivity_overrode_geometric = false;
        Bool interval_lookahead_guided = false;
        Bool interval_lookahead_overrode_geometric = false;
        SizeType interval_lookahead_function_evaluations = 0u;
        double interval_lookahead_evaluation_seconds = 0.0;
    };

    template<class Conjunction>
    SplitBoxResult _split_box(
        UpperBoxType const& domain,
        Conjunction const& conjunction,
        SizeType& derivatives_built,
        SizeType& derivative_evaluations,
        double& derivative_build_seconds,
        double& derivative_evaluation_seconds) const;

    ValidatedScalarMultivariateFunction const& _function(
        ValidatedConstraint const& constraint) const;
    ValidatedScalarMultivariateFunction const& _function(
        CompiledTheoryLiteral const& literal) const;


    SmtSolverConfiguration _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_SMT_SOLVER_HPP
