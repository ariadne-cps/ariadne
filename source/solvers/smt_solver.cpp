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

#include <vector>

#include "solvers/constraint_solver.hpp"
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
        default: ARIADNE_FAIL_MSG("Unknown SmtResultStatus");
    }
}

ExactIntervalType SmtSolver::_epsilon_bounds(ValidatedConstraint const& constraint) const
{
    FloatDP epsilon(_configuration.epsilon(),dp);
    ExactIntervalType bounds=constraint.bounds();
    return ExactIntervalType(
        sub(down,bounds.lower_bound(),epsilon),
        add(up,bounds.upper_bound(),epsilon));
}

Bool SmtSolver::_epsilon_reduce(UpperBoxType& domain,
                                List<ValidatedConstraint> const& constraints) const
{
    ConstraintSolver contractor;
    for(SizeType i=0; i!=constraints.size(); ++i) {
        if(contractor.hull_reduce(domain,constraints[i].function(),this->_epsilon_bounds(constraints[i]))) {
            return true;
        }
    }
    return definitely(domain.is_empty());
}

Bool SmtSolver::_epsilon_satisfied(UpperBoxType const& domain,
                                   List<ValidatedConstraint> const& constraints) const
{
    UpperBoxType point(domain.dimension(),[&](SizeType i) {
        auto m=domain[i].midpoint();
        return UpperIntervalType(m,m);
    });

    for(SizeType i=0; i!=constraints.size(); ++i) {
        auto const& constraint=constraints[i];
        UpperIntervalType image=apply(constraint.function(),point);
        if(not definitely(subset(image,this->_epsilon_bounds(constraint)))) {
            return false;
        }
    }
    return true;
}

SmtResult SmtSolver::solve(ExactBoxType const& domain,
                           List<ValidatedConstraint> const& constraints) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    for(SizeType i=0; i!=constraints.size(); ++i) {
        ARIADNE_PRECONDITION(constraints[i].argument_size()==domain.dimension());
    }

    if(domain.is_empty()) {
        return SmtResult::unsat();
    }

    std::vector<UpperBoxType> pending;
    pending.emplace_back(domain);

    while(not pending.empty()) {
        UpperBoxType current=std::move(pending.back());
        pending.pop_back();

        if(this->_epsilon_reduce(current,constraints)) {
            continue;
        }

        if(this->_epsilon_satisfied(current,constraints)) {
            UpperBoxType witness(current.dimension(),[&](SizeType i) {
                auto m=current[i].midpoint();
                return UpperIntervalType(m,m);
            });
            return SmtResult::epsilon_sat(witness);
        }

        Pair<UpperBoxType,UpperBoxType> children=current.split();

        Bool first_same=true;
        Bool second_same=true;
        for(SizeType i=0; i!=current.dimension(); ++i) {
            first_same = first_same
                and children.first[i].lower_bound().raw()==current[i].lower_bound().raw()
                and children.first[i].upper_bound().raw()==current[i].upper_bound().raw();
            second_same = second_same
                and children.second[i].lower_bound().raw()==current[i].lower_bound().raw()
                and children.second[i].upper_bound().raw()==current[i].upper_bound().raw();
        }
        ARIADNE_ASSERT_MSG(not (first_same and second_same),
                           "SMT search reached a non-splittable uncertified box: "<<current);

        pending.push_back(std::move(children.second));
        pending.push_back(std::move(children.first));
    }

    return SmtResult::unsat();
}

} // namespace Ariadne
