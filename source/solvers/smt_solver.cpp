/***************************************************************************
 *            solvers/smt_solver.cpp
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

#include "solvers/smt_solver.hpp"

#include "utility/exceptions.hpp"

namespace Ariadne {

SmtSolverConfiguration::SmtSolverConfiguration(ExactDouble epsilon)
    : _epsilon(epsilon)
{
    ARIADNE_PRECONDITION(epsilon>ExactDouble(0));
}

SmtResult::SmtResult(SmtResultStatus status)
    : _status(status), _witness()
{
}

SmtResult::SmtResult(SmtResultStatus status, UpperBoxType const& witness)
    : _status(status), _witness(witness)
{
}

SmtResult SmtResult::unsat()
{
    return SmtResult(SmtResultStatus::UNSAT);
}

SmtResult SmtResult::epsilon_sat(UpperBoxType const& witness)
{
    return SmtResult(SmtResultStatus::EPSILON_SAT,witness);
}

UpperBoxType const& SmtResult::witness() const
{
    ARIADNE_PRECONDITION(this->has_witness());
    return *_witness;
}

OutputStream& operator<<(OutputStream& os, SmtResultStatus status)
{
    switch(status) {
        case SmtResultStatus::UNSAT: return os << "UNSAT";
        case SmtResultStatus::EPSILON_SAT: return os << "EPSILON_SAT";
    }
    ARIADNE_FAIL_MSG("Unknown SmtResultStatus");
}

} // namespace Ariadne
