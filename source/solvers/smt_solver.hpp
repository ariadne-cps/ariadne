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

#include "geometry/box.hpp"
#include "numeric/numeric.hpp"

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
//! \brief Result of a bounded real epsilon-SMT query.
class SmtResult {
  public:
    //! \brief Construct an UNSAT result.
    static SmtResult unsat();

    //! \brief Construct an EPSILON_SAT result with a validated witness box.
    static SmtResult epsilon_sat(UpperBoxType const& witness);

    SmtResultStatus status() const { return _status; }
    Bool is_unsat() const { return _status==SmtResultStatus::UNSAT; }
    Bool is_epsilon_sat() const { return _status==SmtResultStatus::EPSILON_SAT; }

    //! \brief True iff the result contains a witness box.
    Bool has_witness() const { return _witness.has_value(); }

    //! \brief Return the witness box.
    //! \pre The result is EPSILON_SAT.
    UpperBoxType const& witness() const;

  private:
    explicit SmtResult(SmtResultStatus status);
    SmtResult(SmtResultStatus status, UpperBoxType const& witness);

    SmtResultStatus _status;
    std::optional<UpperBoxType> _witness;
};

OutputStream& operator<<(OutputStream& os, SmtResultStatus status);

} // namespace Ariadne

#endif // ARIADNE_SMT_SOLVER_HPP
