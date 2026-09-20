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
#include <utility>
#include <atomic>
#include <memory>
#include <mutex>

#include "betterthreads/workload.hpp"

#include "solvers/constraint_solver.hpp"
#include "utility/exceptions.hpp"

namespace Ariadne {

namespace {

class SequentialSmtWorkQueue {
  public:
    Void push(UpperBoxType box) { _boxes.push_back(std::move(box)); }

    UpperBoxType pop() {
        ARIADNE_PRECONDITION(not _boxes.empty());
        UpperBoxType box=std::move(_boxes.back());
        _boxes.pop_back();
        return box;
    }

    Bool empty() const { return _boxes.empty(); }

  private:
    std::vector<UpperBoxType> _boxes;
};

} // namespace

SmtSolverConfiguration::SmtSolverConfiguration(ExactDouble epsilon)
    : _epsilon(epsilon)
{
    ARIADNE_PRECONDITION(epsilon>ExactDouble(0));
}

SmtResult::SmtResult(SmtResultStatus status, SmtSearchStatistics statistics)
    : _status(status), _witness(), _statistics(statistics)
{
}

SmtResult::SmtResult(SmtResultStatus status, UpperBoxType const& witness,
                     SmtSearchStatistics statistics)
    : _status(status), _witness(witness), _statistics(statistics)
{
}

SmtResult SmtResult::unsat(SmtSearchStatistics statistics)
{
    return SmtResult(SmtResultStatus::UNSAT,statistics);
}

SmtResult SmtResult::epsilon_sat(UpperBoxType const& witness,
                                 SmtSearchStatistics statistics)
{
    return SmtResult(SmtResultStatus::EPSILON_SAT,witness,statistics);
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

SmtSolver::BoxProcessingResult
SmtSolver::_process_box(UpperBoxType domain,
                        List<ValidatedConstraint> const& constraints) const
{
    if(this->_epsilon_reduce(domain,constraints)) {
        return {BoxProcessingStatus::PRUNED,std::nullopt,std::nullopt};
    }

    if(this->_epsilon_satisfied(domain,constraints)) {
        UpperBoxType witness(domain.dimension(),[&](SizeType i) {
            auto m=domain[i].midpoint();
            return UpperIntervalType(m,m);
        });
        return {BoxProcessingStatus::EPSILON_SAT,witness,std::nullopt};
    }

    Pair<UpperBoxType,UpperBoxType> children=domain.split();

    Bool first_same=true;
    Bool second_same=true;
    for(SizeType i=0; i!=domain.dimension(); ++i) {
        first_same = first_same
            and children.first[i].lower_bound().raw()==domain[i].lower_bound().raw()
            and children.first[i].upper_bound().raw()==domain[i].upper_bound().raw();
        second_same = second_same
            and children.second[i].lower_bound().raw()==domain[i].lower_bound().raw()
            and children.second[i].upper_bound().raw()==domain[i].upper_bound().raw();
    }
    ARIADNE_ASSERT_MSG(not (first_same and second_same),
                       "SMT search reached a non-splittable uncertified box: "<<domain);

    return {BoxProcessingStatus::SPLIT,std::nullopt,children};
}

SmtResult SmtSolver::solve(ExactBoxType const& domain,
                           List<ValidatedConstraint> const& constraints) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    for(SizeType i=0; i!=constraints.size(); ++i) {
        ARIADNE_PRECONDITION(constraints[i].argument_size()==domain.dimension());
    }

    SmtSearchStatistics statistics;

    if(domain.is_empty()) {
        return SmtResult::unsat(statistics);
    }

    SequentialSmtWorkQueue pending;
    pending.push(UpperBoxType(domain));

    while(not pending.empty()) {
        UpperBoxType current=pending.pop();
        ++statistics.boxes_processed;

        BoxProcessingResult processing=this->_process_box(std::move(current),constraints);
        switch(processing.status) {
            case BoxProcessingStatus::PRUNED:
                ++statistics.boxes_pruned;
                break;

            case BoxProcessingStatus::EPSILON_SAT:
                ARIADNE_ASSERT(processing.witness.has_value());
                return SmtResult::epsilon_sat(*processing.witness,statistics);

            case BoxProcessingStatus::SPLIT:
                ARIADNE_ASSERT(processing.children.has_value());
                ++statistics.boxes_split;
                pending.push(std::move(processing.children->second));
                pending.push(std::move(processing.children->first));
                break;

            default:
                ARIADNE_FAIL_MSG("Unknown BoxProcessingStatus");
        }
    }

    return SmtResult::unsat(statistics);
}


namespace {

struct ParallelSmtSearchState {
    std::mutex mutex;
    SmtSearchStatistics statistics;
    std::optional<UpperBoxType> witness;
    std::atomic<bool> found{false};
};

using ParallelSmtWorkload = BetterThreads::DynamicWorkload<UpperBoxType>;

} // namespace

SmtResult SmtSolver::solve_parallel(ExactBoxType const& domain,
                                    List<ValidatedConstraint> const& constraints) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    for(SizeType i=0; i!=constraints.size(); ++i) {
        ARIADNE_PRECONDITION(constraints[i].argument_size()==domain.dimension());
    }

    auto state=std::make_shared<ParallelSmtSearchState>();
    if(domain.is_empty()) {
        return SmtResult::unsat(state->statistics);
    }

    ParallelSmtWorkload workload(
        [](UpperBoxType const&, std::shared_ptr<ConcLog::ProgressIndicator>) { },
        [this,&constraints,state](ParallelSmtWorkload::Access& access, UpperBoxType const& box) {
            if(state->found.load()) {
                return;
            }

            BoxProcessingResult processing=this->_process_box(box,constraints);
            {
                std::lock_guard<std::mutex> lock(state->mutex);
                ++state->statistics.boxes_processed;
                if(processing.status==BoxProcessingStatus::PRUNED) {
                    ++state->statistics.boxes_pruned;
                } else if(processing.status==BoxProcessingStatus::SPLIT) {
                    ++state->statistics.boxes_split;
                }
            }

            switch(processing.status) {
                case BoxProcessingStatus::PRUNED:
                    return;

                case BoxProcessingStatus::EPSILON_SAT: {
                    ARIADNE_ASSERT(processing.witness.has_value());
                    bool expected=false;
                    if(state->found.compare_exchange_strong(expected,true)) {
                        std::lock_guard<std::mutex> lock(state->mutex);
                        state->witness=*processing.witness;
                    }
                    return;
                }

                case BoxProcessingStatus::SPLIT:
                    ARIADNE_ASSERT(processing.children.has_value());
                    if(not state->found.load()) {
                        access.append(processing.children->first);
                        access.append(processing.children->second);
                    }
                    return;

                default:
                    ARIADNE_FAIL_MSG("Unknown BoxProcessingStatus");
            }
        });

    workload.append(UpperBoxType(domain));
    workload.process();

    std::lock_guard<std::mutex> lock(state->mutex);
    if(state->found.load()) {
        ARIADNE_ASSERT(state->witness.has_value());
        return SmtResult::epsilon_sat(*state->witness,state->statistics);
    }
    return SmtResult::unsat(state->statistics);
}

} // namespace Ariadne
