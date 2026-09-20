/***************************************************************************
 *            solvers/smt_solver.hpp
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

/*! \file solvers/smt_solver.hpp
 *  \brief Common types for solving bounded real SMT problems.
 */

#ifndef ARIADNE_SMT_SOLVER_HPP
#define ARIADNE_SMT_SOLVER_HPP

#include <optional>
#include <mutex>

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
    EPSILON_SAT
};

//! \ingroup Solvers
//! \brief Configuration shared by bounded real epsilon-SMT solvers.
class SmtSolverConfiguration {
  public:
    explicit SmtSolverConfiguration(ExactDouble epsilon);

    //! \brief The logical epsilon used for weakening constraints.
    ExactDouble epsilon() const { return _epsilon; }

  private:
    ExactDouble _epsilon;
};

//! \ingroup Solvers
//! \brief Search statistics for a bounded real epsilon-SMT query.
struct SmtSearchStatistics {
    SizeType boxes_processed = 0u;
    SizeType boxes_pruned = 0u;
    SizeType boxes_split = 0u;
    SizeType boolean_decisions = 0u;
    SizeType boolean_propagations = 0u;
    SizeType boolean_reasoned_propagations = 0u;
    SizeType boolean_conflicts = 0u;
    SizeType boolean_backtracks = 0u;
    SizeType max_decision_level = 0u;
    SizeType theory_checks = 0u;
    SizeType theory_conflicts = 0u;
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

    SmtResultStatus status() const { return _status; }
    Bool is_unsat() const { return _status==SmtResultStatus::UNSAT; }
    Bool is_epsilon_sat() const { return _status==SmtResultStatus::EPSILON_SAT; }

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
        SPLIT
    };

    struct BoxProcessingResult {
        BoxProcessingStatus status;
        std::optional<UpperBoxType> witness;
        std::optional<Pair<UpperBoxType,UpperBoxType>> children;
    };

    struct CompiledTheoryLiteral {
        ValidatedScalarMultivariateFunction function;
        SmtTheoryPrimitiveRelation relation;
    };
    using CompiledTheoryLiterals = std::vector<CompiledTheoryLiteral>;

    ExactIntervalType _epsilon_bounds(ValidatedConstraint const& constraint) const;
    ExactIntervalType _epsilon_bounds(SmtTheoryPrimitiveRelation relation) const;
    Bool _epsilon_reduce(UpperBoxType& domain,
                         List<ValidatedConstraint> const& constraints) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            List<ValidatedConstraint> const& constraints) const;
    BoxProcessingResult _process_box(UpperBoxType domain,
                                     List<ValidatedConstraint> const& constraints) const;

    CompiledTheoryLiterals _compile_theory_literals(
        RealSpace const& space,
        List<SmtTheoryPrimitiveLiteral> const& literals) const;
    Bool _epsilon_reduce(UpperBoxType& domain,
                         CompiledTheoryLiterals const& literals) const;
    Bool _epsilon_satisfied(UpperBoxType const& domain,
                            CompiledTheoryLiterals const& literals) const;
    BoxProcessingResult _process_box(UpperBoxType domain,
                                     CompiledTheoryLiterals const& literals) const;

    SmtSolverConfiguration _configuration;
};

} // namespace Ariadne

#endif // ARIADNE_SMT_SOLVER_HPP
