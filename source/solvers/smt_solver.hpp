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
        SizeType box_processing_limit=std::numeric_limits<SizeType>::max());

    //! \brief The logical epsilon used for weakening constraints.
    ExactDouble epsilon() const { return _epsilon; }

    //! \brief Maximum theory checks used to minimize one learned theory nogood.
    SizeType theory_minimization_budget() const { return _theory_minimization_budget; }

    //! \brief Maximum number of active non-theory learned clauses.
    SizeType learned_clause_limit() const { return _learned_clause_limit; }

    //! \brief Maximum number of search boxes processed by one theory solve.
    //! \details Reaching the limit yields UNKNOWN unless an epsilon witness was found first.
    SizeType box_processing_limit() const { return _box_processing_limit; }

  private:
    ExactDouble _epsilon;
    SizeType _theory_minimization_budget;
    SizeType _learned_clause_limit;
    SizeType _box_processing_limit;
};

//! \ingroup Solvers
//! \brief Search statistics for a bounded real epsilon-SMT query.
struct SmtSearchStatistics {
    SizeType boxes_processed = 0u;
    SizeType boxes_pruned = 0u;
    SizeType boxes_split = 0u;
    SizeType boxes_unknown = 0u;
    SizeType hull_reduction_rounds = 0u;
    SizeType hull_effective_reductions = 0u;
    SizeType shaving_reduction_rounds = 0u;
    SizeType shaving_effective_reductions = 0u;
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
    enum class BoxProcessingStatus {
        PRUNED,
        EPSILON_SAT,
        SPLIT,
        UNKNOWN
    };

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
    };

    struct CompiledTheoryLiteral {
        ValidatedScalarMultivariateFunction function;
        SmtTheoryPrimitiveRelation relation;
    };
    using CompiledTheoryLiterals = std::vector<CompiledTheoryLiteral>;

    ExactIntervalType _epsilon_bounds(ValidatedConstraint const& constraint) const;
    ExactIntervalType _epsilon_bounds(SmtTheoryPrimitiveRelation relation) const;
    Bool _epsilon_reduce(UpperBoxType& domain,
                         List<ValidatedConstraint> const& constraints,
                         ReductionStatistics& statistics) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            List<ValidatedConstraint> const& constraints) const;
    std::optional<UpperBoxType> _epsilon_witness(
        UpperBoxType const& domain,
        List<ValidatedConstraint> const& constraints) const;
    BoxProcessingResult _process_box(UpperBoxType domain,
                                     List<ValidatedConstraint> const& constraints) const;
    Pair<UpperBoxType,UpperBoxType> _split_box(
        UpperBoxType const& domain,
        List<ValidatedConstraint> const& constraints) const;

    CompiledTheoryLiterals _compile_theory_literals(
        RealSpace const& space,
        List<SmtTheoryPrimitiveLiteral> const& literals) const;
    Bool _epsilon_reduce(UpperBoxType& domain,
                         CompiledTheoryLiterals const& literals,
                         ReductionStatistics& statistics) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            CompiledTheoryLiterals const& literals) const;
    std::optional<UpperBoxType> _epsilon_witness(
        UpperBoxType const& domain,
        CompiledTheoryLiterals const& literals) const;
    BoxProcessingResult _process_box(UpperBoxType domain,
                                     CompiledTheoryLiterals const& literals) const;
    Pair<UpperBoxType,UpperBoxType> _split_box(
        UpperBoxType const& domain,
        CompiledTheoryLiterals const& literals) const;

    SmtSolverConfiguration _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_SMT_SOLVER_HPP
