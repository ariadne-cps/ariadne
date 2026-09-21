/***************************************************************************
 *            solvers/smt_solver.cpp
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

#include "solvers/smt_solver.hpp"

#include <vector>
#include <utility>
#include <atomic>
#include <memory>
#include <mutex>
#include <cstdint>
#include <algorithm>

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

Bool same_box(UpperBoxType const& first, UpperBoxType const& second)
{
    ARIADNE_ASSERT(first.dimension()==second.dimension());
    for(SizeType i=0; i!=first.dimension(); ++i) {
        if(first[i].lower_bound().raw()!=second[i].lower_bound().raw()
           || first[i].upper_bound().raw()!=second[i].upper_bound().raw()) {
            return false;
        }
    }
    return true;
}

Pair<SizeType,Bool> sensitivity_split_coordinate(
    UpperBoxType const& domain,
    std::vector<ValidatedScalarMultivariateFunction> const& functions)
{
    auto widths=domain.widths();

    SizeType geometric=0u;
    for(SizeType variable=1u; variable!=domain.dimension(); ++variable) {
        if(definitely(widths[variable]>widths[geometric])) {
            geometric=variable;
        }
    }

    std::optional<SizeType> selected;
    std::optional<PositiveFloatDPUpperBound> selected_score;
    for(SizeType variable=0u; variable!=domain.dimension(); ++variable) {
        PositiveFloatDPUpperBound sensitivity(0u,dp);
        Bool active=false;
        for(auto const& function:functions) {
            UpperIntervalType derivative_image=apply(function.derivative(variable),domain);
            if(not definitely(derivative_image.lower_bound()==0)
               || not definitely(derivative_image.upper_bound()==0)) {
                active=true;
                PositiveFloatDPUpperBound candidate=
                    domain[variable].width()*mag(derivative_image);
                if(definitely(candidate>sensitivity)) {
                    sensitivity=candidate;
                }
            }
        }
        if(active && (not selected_score.has_value()
                      || definitely(sensitivity>*selected_score))) {
            selected=variable;
            selected_score=sensitivity;
        }
    }

    SizeType coordinate=selected.has_value() ? *selected : geometric;
    return {coordinate,coordinate!=geometric};
}

std::vector<UpperBoxType> epsilon_witness_candidates(UpperBoxType const& domain)
{
    constexpr SizeType max_corner_candidates=64u;

    std::vector<UpperBoxType> candidates;
    candidates.reserve(3u+2u*domain.dimension()+max_corner_candidates);

    auto midpoint_point=[&]() {
        return UpperBoxType(domain.dimension(),[&](SizeType i) {
            auto m=domain[i].midpoint();
            return UpperIntervalType(m,m);
        });
    };
    UpperBoxType midpoint=midpoint_point();
    candidates.push_back(midpoint);

    UpperBoxType lower=midpoint;
    UpperBoxType upper=midpoint;
    for(SizeType i=0u; i!=domain.dimension(); ++i) {
        auto l=domain[i].lower_bound().raw();
        auto u=domain[i].upper_bound().raw();
        lower[i]=UpperIntervalType(l,l);
        upper[i]=UpperIntervalType(u,u);
    }
    candidates.push_back(lower);
    candidates.push_back(upper);

    for(SizeType i=0u; i!=domain.dimension(); ++i) {
        UpperBoxType low_axis=midpoint;
        UpperBoxType high_axis=midpoint;
        auto l=domain[i].lower_bound().raw();
        auto u=domain[i].upper_bound().raw();
        low_axis[i]=UpperIntervalType(l,l);
        high_axis[i]=UpperIntervalType(u,u);
        candidates.push_back(std::move(low_axis));
        candidates.push_back(std::move(high_axis));
    }

    if(domain.dimension()<std::numeric_limits<std::uint64_t>::digits) {
        std::uint64_t corner_count=std::uint64_t(1) << domain.dimension();
        if(corner_count<=max_corner_candidates) {
            for(std::uint64_t mask=0u; mask<corner_count; ++mask) {
                UpperBoxType corner=midpoint;
                for(SizeType i=0u; i!=domain.dimension(); ++i) {
                    Bool use_upper=(mask & (std::uint64_t(1) << i))!=0u;
                    auto endpoint=use_upper
                        ? domain[i].upper_bound().raw()
                        : domain[i].lower_bound().raw();
                    corner[i]=UpperIntervalType(endpoint,endpoint);
                }
                candidates.push_back(std::move(corner));
            }
        }
    }

    return candidates;
}

} // namespace

SmtSolverConfiguration::SmtSolverConfiguration(
    ExactDouble epsilon, SizeType theory_minimization_budget,
    SizeType learned_clause_limit, SizeType box_processing_limit)
    : _epsilon(epsilon),
      _theory_minimization_budget(theory_minimization_budget),
      _learned_clause_limit(learned_clause_limit),
      _box_processing_limit(box_processing_limit)
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

SmtResult SmtResult::unknown(SmtSearchStatistics statistics)
{
    return SmtResult(SmtResultStatus::UNKNOWN,statistics);
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
        case SmtResultStatus::UNKNOWN: return os << "UNKNOWN";
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

ExactIntervalType SmtSolver::_epsilon_bounds(SmtTheoryPrimitiveRelation relation) const
{
    FloatDP epsilon(_configuration.epsilon(),dp);
    switch(relation) {
        case SmtTheoryPrimitiveRelation::EQ_ZERO:
            return ExactIntervalType(-epsilon,+epsilon);
        case SmtTheoryPrimitiveRelation::GEQ_ZERO:
        case SmtTheoryPrimitiveRelation::GT_ZERO:
            return ExactIntervalType(-epsilon,+infty);
        default:
            ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
    }
}

Bool SmtSolver::_epsilon_reduce(UpperBoxType& domain,
                                List<ValidatedConstraint> const& constraints,
                                ReductionStatistics& statistics) const
{
    ConstraintSolver contractor;
    while(true) {
        UpperBoxType previous=domain;
        ++statistics.hull_rounds;
        for(SizeType i=0; i!=constraints.size(); ++i) {
            if(contractor.hull_reduce(
                    domain,constraints[i].function(),this->_epsilon_bounds(constraints[i]))) {
                return true;
            }
        }
        if(definitely(domain.is_empty())) {
            return true;
        }

        if(not same_box(domain,previous)) {
            ++statistics.hull_effective;
        }

        if(same_box(domain,previous)) {
            UpperBoxType before_shaving=domain;
            ++statistics.shaving_rounds;
            for(SizeType i=0; i!=constraints.size(); ++i) {
                for(SizeType variable=0u; variable!=domain.dimension(); ++variable) {
                    if(contractor.box_reduce(
                            domain,
                            constraints[i].function(),
                            this->_epsilon_bounds(constraints[i]),
                            variable)) {
                        return true;
                    }
                }
            }
            if(definitely(domain.is_empty())) {
                return true;
            }
            if(not same_box(domain,before_shaving)) {
                ++statistics.shaving_effective;
            }
            if(same_box(domain,before_shaving)) {
                for(SizeType i=0; i!=constraints.size(); ++i) {
                    UpperIntervalType image=apply(constraints[i].function(),domain);
                    if(definitely(disjoint(image,this->_epsilon_bounds(constraints[i])))) {
                        return true;
                    }
                }
                return false;
            }
            continue;
        }

        for(SizeType i=0; i!=constraints.size(); ++i) {
            UpperIntervalType image=apply(constraints[i].function(),domain);
            if(definitely(disjoint(image,this->_epsilon_bounds(constraints[i])))) {
                return true;
            }
        }
    }
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

std::optional<UpperBoxType>
SmtSolver::_epsilon_witness(UpperBoxType const& domain,
                            List<ValidatedConstraint> const& constraints) const
{
    for(UpperBoxType const& candidate:epsilon_witness_candidates(domain)) {
        if(this->_epsilon_satisfied(candidate,constraints)) {
            return candidate;
        }
    }
    return std::nullopt;
}

Pair<Pair<UpperBoxType,UpperBoxType>,Bool>
SmtSolver::_split_box(UpperBoxType const& domain,
                      List<ValidatedConstraint> const& constraints) const
{
    std::vector<ValidatedScalarMultivariateFunction> functions;
    functions.reserve(constraints.size());
    for(SizeType i=0u; i!=constraints.size(); ++i) {
        functions.push_back(constraints[i].function());
    }
    auto selection=sensitivity_split_coordinate(domain,functions);
    return {domain.split(selection.first),selection.second};
}

SmtSolver::BoxProcessingResult
SmtSolver::_process_box(UpperBoxType domain,
                        List<ValidatedConstraint> const& constraints) const
{
    ReductionStatistics reductions;
    if(this->_epsilon_reduce(domain,constraints,reductions)) {
        return {BoxProcessingStatus::PRUNED,std::nullopt,std::nullopt,reductions};
    }

    if(auto witness=this->_epsilon_witness(domain,constraints); witness.has_value()) {
        return {BoxProcessingStatus::EPSILON_SAT,*witness,std::nullopt,reductions};
    }

    auto split_result=this->_split_box(domain,constraints);
    Pair<UpperBoxType,UpperBoxType> children=split_result.first;

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
    if(first_same and second_same) {
        return {BoxProcessingStatus::UNKNOWN,std::nullopt,std::nullopt,reductions};
    }

    return {BoxProcessingStatus::SPLIT,std::nullopt,children,reductions,split_result.second};
}

SmtSolver::CompiledTheoryLiterals
SmtSolver::_compile_theory_literals(RealSpace const& space,
                                    List<SmtTheoryPrimitiveLiteral> const& literals) const
{
    CompiledTheoryLiterals result;
    result.reserve(literals.size());
    for(SizeType i=0; i!=literals.size(); ++i) {
        result.push_back({
            ValidatedScalarMultivariateFunction(space,literals[i].expression()),
            literals[i].relation()
        });
    }
    return result;
}

Bool SmtSolver::_epsilon_reduce(UpperBoxType& domain,
                                CompiledTheoryLiterals const& literals,
                                ReductionStatistics& statistics) const
{
    ConstraintSolver contractor;
    FloatDP epsilon(_configuration.epsilon(),dp);
    while(true) {
        UpperBoxType previous=domain;
        ++statistics.hull_rounds;
        for(auto const& literal:literals) {
            if(contractor.hull_reduce(
                    domain,literal.function,this->_epsilon_bounds(literal.relation))) {
                return true;
            }
            if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
                UpperIntervalType image=apply(literal.function,domain);
                if(definitely(image.upper_bound()<=-epsilon)) {
                    return true;
                }
            }
        }
        if(definitely(domain.is_empty())) {
            return true;
        }

        if(not same_box(domain,previous)) {
            ++statistics.hull_effective;
        }

        if(same_box(domain,previous)) {
            UpperBoxType before_shaving=domain;
            ++statistics.shaving_rounds;
            for(auto const& literal:literals) {
                for(SizeType variable=0u; variable!=domain.dimension(); ++variable) {
                    if(contractor.box_reduce(
                            domain,
                            literal.function,
                            this->_epsilon_bounds(literal.relation),
                            variable)) {
                        return true;
                    }
                }
            }
            if(definitely(domain.is_empty())) {
                return true;
            }
            if(not same_box(domain,before_shaving)) {
                ++statistics.shaving_effective;
            }
            if(same_box(domain,before_shaving)) {
                for(auto const& literal:literals) {
                    UpperIntervalType image=apply(literal.function,domain);
                    switch(literal.relation) {
                        case SmtTheoryPrimitiveRelation::EQ_ZERO:
                        case SmtTheoryPrimitiveRelation::GEQ_ZERO:
                            if(definitely(disjoint(
                                    image,this->_epsilon_bounds(literal.relation)))) {
                                return true;
                            }
                            break;
                        case SmtTheoryPrimitiveRelation::GT_ZERO:
                            if(definitely(image.upper_bound()<=-epsilon)) {
                                return true;
                            }
                            break;
                        default:
                            ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
                    }
                }
                return false;
            }
            continue;
        }

        for(auto const& literal:literals) {
            UpperIntervalType image=apply(literal.function,domain);
            switch(literal.relation) {
                case SmtTheoryPrimitiveRelation::EQ_ZERO:
                case SmtTheoryPrimitiveRelation::GEQ_ZERO:
                    if(definitely(disjoint(image,this->_epsilon_bounds(literal.relation)))) {
                        return true;
                    }
                    break;
                case SmtTheoryPrimitiveRelation::GT_ZERO:
                    if(definitely(image.upper_bound()<=-epsilon)) {
                        return true;
                    }
                    break;
                default:
                    ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
            }
        }
    }
}

Bool SmtSolver::_epsilon_satisfied(UpperBoxType const& domain,
                                   CompiledTheoryLiterals const& literals) const
{
    UpperBoxType point(domain.dimension(),[&](SizeType i) {
        auto m=domain[i].midpoint();
        return UpperIntervalType(m,m);
    });

    FloatDP epsilon(_configuration.epsilon(),dp);
    for(auto const& literal:literals) {
        UpperIntervalType image=apply(literal.function,point);
        switch(literal.relation) {
            case SmtTheoryPrimitiveRelation::EQ_ZERO:
            case SmtTheoryPrimitiveRelation::GEQ_ZERO:
                if(not definitely(subset(image,this->_epsilon_bounds(literal.relation)))) {
                    return false;
                }
                break;
            case SmtTheoryPrimitiveRelation::GT_ZERO:
                if(not definitely(image.lower_bound()>-epsilon)) {
                    return false;
                }
                break;
            default:
                ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
        }
    }
    return true;
}

std::optional<UpperBoxType>
SmtSolver::_epsilon_witness(UpperBoxType const& domain,
                            CompiledTheoryLiterals const& literals) const
{
    for(UpperBoxType const& candidate:epsilon_witness_candidates(domain)) {
        if(this->_epsilon_satisfied(candidate,literals)) {
            return candidate;
        }
    }
    return std::nullopt;
}

Pair<Pair<UpperBoxType,UpperBoxType>,Bool>
SmtSolver::_split_box(UpperBoxType const& domain,
                      CompiledTheoryLiterals const& literals) const
{
    std::vector<ValidatedScalarMultivariateFunction> functions;
    functions.reserve(literals.size());
    for(auto const& literal:literals) {
        functions.push_back(literal.function);
    }
    auto selection=sensitivity_split_coordinate(domain,functions);
    return {domain.split(selection.first),selection.second};
}

SmtSolver::BoxProcessingResult
SmtSolver::_process_box(UpperBoxType domain,
                        CompiledTheoryLiterals const& literals) const
{
    ReductionStatistics reductions;
    if(this->_epsilon_reduce(domain,literals,reductions)) {
        return {BoxProcessingStatus::PRUNED,std::nullopt,std::nullopt,reductions};
    }

    if(auto witness=this->_epsilon_witness(domain,literals); witness.has_value()) {
        return {BoxProcessingStatus::EPSILON_SAT,*witness,std::nullopt,reductions};
    }

    auto split_result=this->_split_box(domain,literals);
    Pair<UpperBoxType,UpperBoxType> children=split_result.first;

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
    if(first_same and second_same) {
        return {BoxProcessingStatus::UNKNOWN,std::nullopt,std::nullopt,reductions};
    }

    return {BoxProcessingStatus::SPLIT,std::nullopt,children,reductions,split_result.second};
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
    Bool unknown_seen=false;

    while(not pending.empty()) {
        if(statistics.boxes_processed>=_configuration.box_processing_limit()) {
            unknown_seen=true;
            break;
        }
        UpperBoxType current=pending.pop();
        ++statistics.boxes_processed;

        BoxProcessingResult processing=this->_process_box(std::move(current),constraints);
        statistics.hull_reduction_rounds+=processing.reductions.hull_rounds;
        statistics.hull_effective_reductions+=processing.reductions.hull_effective;
        statistics.shaving_reduction_rounds+=processing.reductions.shaving_rounds;
        statistics.shaving_effective_reductions+=processing.reductions.shaving_effective;
        if(processing.sensitivity_guided_split) {
            ++statistics.sensitivity_guided_splits;
        }
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

            case BoxProcessingStatus::UNKNOWN:
                ++statistics.boxes_unknown;
                unknown_seen=true;
                break;

            default:
                ARIADNE_FAIL_MSG("Unknown BoxProcessingStatus");
        }
    }

    if(unknown_seen) {
        return SmtResult::unknown(statistics);
    }
    return SmtResult::unsat(statistics);
}


SmtResult SmtSolver::solve(RealSpace const& space,
                           ExactBoxType const& domain,
                           List<SmtTheoryPrimitiveLiteral> const& literals) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    ARIADNE_PRECONDITION(space.size()==domain.dimension());

    SmtSearchStatistics statistics;
    if(domain.is_empty()) {
        return SmtResult::unsat(statistics);
    }

    CompiledTheoryLiterals compiled=this->_compile_theory_literals(space,literals);
    SequentialSmtWorkQueue pending;
    pending.push(UpperBoxType(domain));
    Bool unknown_seen=false;

    while(not pending.empty()) {
        if(statistics.boxes_processed>=_configuration.box_processing_limit()) {
            unknown_seen=true;
            break;
        }
        UpperBoxType current=pending.pop();
        ++statistics.boxes_processed;

        BoxProcessingResult processing=this->_process_box(std::move(current),compiled);
        statistics.hull_reduction_rounds+=processing.reductions.hull_rounds;
        statistics.hull_effective_reductions+=processing.reductions.hull_effective;
        statistics.shaving_reduction_rounds+=processing.reductions.shaving_rounds;
        statistics.shaving_effective_reductions+=processing.reductions.shaving_effective;
        if(processing.sensitivity_guided_split) {
            ++statistics.sensitivity_guided_splits;
        }
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
            case BoxProcessingStatus::UNKNOWN:
                ++statistics.boxes_unknown;
                unknown_seen=true;
                break;
            default:
                ARIADNE_FAIL_MSG("Unknown BoxProcessingStatus");
        }
    }

    if(unknown_seen) {
        return SmtResult::unknown(statistics);
    }
    return SmtResult::unsat(statistics);
}


namespace {

struct ParallelSmtSearchState {
    std::mutex mutex;
    SmtSearchStatistics statistics;
    std::optional<UpperBoxType> witness;
    std::atomic<bool> found{false};
    std::atomic<bool> unknown{false};
    std::atomic<bool> limit_reached{false};
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
            if(state->found.load() || state->limit_reached.load()) {
                return;
            }

            {
                std::lock_guard<std::mutex> lock(state->mutex);
                if(state->statistics.boxes_processed>=_configuration.box_processing_limit()) {
                    state->unknown.store(true);
                    state->limit_reached.store(true);
                    return;
                }
                ++state->statistics.boxes_processed;
            }

            BoxProcessingResult processing=this->_process_box(box,constraints);
            {
                std::lock_guard<std::mutex> lock(state->mutex);
                state->statistics.hull_reduction_rounds+=processing.reductions.hull_rounds;
                state->statistics.hull_effective_reductions+=processing.reductions.hull_effective;
                state->statistics.shaving_reduction_rounds+=processing.reductions.shaving_rounds;
                state->statistics.shaving_effective_reductions+=processing.reductions.shaving_effective;
                if(processing.sensitivity_guided_split) {
                    ++state->statistics.sensitivity_guided_splits;
                }
                if(processing.sensitivity_guided_split) {
                    ++state->statistics.sensitivity_guided_splits;
                }
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

                case BoxProcessingStatus::UNKNOWN:
                    state->unknown.store(true);
                    {
                        std::lock_guard<std::mutex> lock(state->mutex);
                        ++state->statistics.boxes_unknown;
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
    if(state->unknown.load()) {
        return SmtResult::unknown(state->statistics);
    }
    return SmtResult::unsat(state->statistics);
}

SmtResult SmtSolver::solve_parallel(RealSpace const& space,
                                    ExactBoxType const& domain,
                                    List<SmtTheoryPrimitiveLiteral> const& literals) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    ARIADNE_PRECONDITION(space.size()==domain.dimension());

    auto state=std::make_shared<ParallelSmtSearchState>();
    if(domain.is_empty()) {
        return SmtResult::unsat(state->statistics);
    }

    CompiledTheoryLiterals compiled=this->_compile_theory_literals(space,literals);
    ParallelSmtWorkload workload(
        [](UpperBoxType const&, std::shared_ptr<ConcLog::ProgressIndicator>) { },
        [this,&compiled,state](ParallelSmtWorkload::Access& access, UpperBoxType const& box) {
            if(state->found.load() || state->limit_reached.load()) {
                return;
            }

            {
                std::lock_guard<std::mutex> lock(state->mutex);
                if(state->statistics.boxes_processed>=_configuration.box_processing_limit()) {
                    state->unknown.store(true);
                    state->limit_reached.store(true);
                    return;
                }
                ++state->statistics.boxes_processed;
            }

            BoxProcessingResult processing=this->_process_box(box,compiled);
            {
                std::lock_guard<std::mutex> lock(state->mutex);
                state->statistics.hull_reduction_rounds+=processing.reductions.hull_rounds;
                state->statistics.hull_effective_reductions+=processing.reductions.hull_effective;
                state->statistics.shaving_reduction_rounds+=processing.reductions.shaving_rounds;
                state->statistics.shaving_effective_reductions+=processing.reductions.shaving_effective;
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
                case BoxProcessingStatus::UNKNOWN:
                    state->unknown.store(true);
                    {
                        std::lock_guard<std::mutex> lock(state->mutex);
                        ++state->statistics.boxes_unknown;
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


namespace {

Void add_statistics(SmtSearchStatistics& target, SmtSearchStatistics const& source)
{
    target.boxes_processed+=source.boxes_processed;
    target.boxes_pruned+=source.boxes_pruned;
    target.boxes_split+=source.boxes_split;
    target.boxes_unknown+=source.boxes_unknown;
    target.hull_reduction_rounds+=source.hull_reduction_rounds;
    target.hull_effective_reductions+=source.hull_effective_reductions;
    target.shaving_reduction_rounds+=source.shaving_reduction_rounds;
    target.shaving_effective_reductions+=source.shaving_effective_reductions;
    target.sensitivity_guided_splits+=source.sensitivity_guided_splits;
    target.boolean_decisions+=source.boolean_decisions;
    target.boolean_propagations+=source.boolean_propagations;
    target.boolean_reasoned_propagations+=source.boolean_reasoned_propagations;
    target.boolean_conflicts+=source.boolean_conflicts;
    target.boolean_backtracks+=source.boolean_backtracks;
    target.max_decision_level=std::max(target.max_decision_level,source.max_decision_level);
    target.boolean_conflicts_analyzed+=source.boolean_conflicts_analyzed;
    target.learned_clause_literals+=source.learned_clause_literals;
    target.last_learned_clause_literals=source.last_learned_clause_literals;
    target.last_learned_current_level_literals=source.last_learned_current_level_literals;
    target.last_backjump_level=source.last_backjump_level;
    target.learned_clauses+=source.learned_clauses;
    target.learned_clause_propagations+=source.learned_clause_propagations;
    target.nonchronological_backjumps+=source.nonchronological_backjumps;
    target.theory_checks+=source.theory_checks;
    target.theory_conflicts+=source.theory_conflicts;
    target.theory_learned_clauses+=source.theory_learned_clauses;
    target.theory_learned_clause_literals+=source.theory_learned_clause_literals;
    target.theory_learned_clause_propagations+=source.theory_learned_clause_propagations;
    target.theory_minimization_checks+=source.theory_minimization_checks;
    target.theory_nogood_raw_literals+=source.theory_nogood_raw_literals;
    target.theory_nogood_minimized_literals+=source.theory_nogood_minimized_literals;
    target.theory_nogood_literals_removed+=source.theory_nogood_literals_removed;
    target.theory_minimization_budget_exhaustions+=
        source.theory_minimization_budget_exhaustions;
    if(target.first_minimization_candidate_trail_rank==0u
       && source.first_minimization_candidate_trail_rank!=0u) {
        target.first_minimization_candidate_trail_rank=
            source.first_minimization_candidate_trail_rank;
    }
    target.learned_clause_activity_bumps+=source.learned_clause_activity_bumps;
    target.learned_clause_pruning_runs+=source.learned_clause_pruning_runs;
    target.learned_clauses_pruned+=source.learned_clauses_pruned;
    target.peak_active_non_theory_learned_clauses=std::max(
        target.peak_active_non_theory_learned_clauses,
        source.peak_active_non_theory_learned_clauses);
}

class SmtDpllSearch {
  public:
    SmtDpllSearch(SmtSolver const& solver,
                  RealSpace const& space,
                  ExactBoxType const& domain,
                  SmtBooleanEncoding const& encoding,
                  Bool parallel)
        : _solver(solver),
          _space(space),
          _domain(domain),
          _encoding(encoding),
          _parallel(parallel),
          _assignment(encoding.variable_count()+1u)
    {
        _trail.reserve(encoding.variable_count());
        _decision_level_markers.reserve(encoding.variable_count()+1u);
        _decision_level_markers.push_back(0u);
    }

    SmtResult solve()
    {
        SearchOutcome outcome=this->_search_boolean();
        if(outcome.witness.has_value()) {
            return SmtResult::epsilon_sat(*outcome.witness,_statistics);
        }
        if(_theory_unknown_seen) {
            return SmtResult::unknown(_statistics);
        }
        return SmtResult::unsat(_statistics);
    }

  private:
    struct AssignmentInfo {
        int8_t value = -1;
        SizeType decision_level = 0u;
        std::optional<SizeType> reason_clause;
    };

    struct ConflictAnalysis {
        std::vector<Int> learned_clause;
        SizeType backjump_level = 0u;
    };

    struct SearchOutcome {
        std::optional<UpperBoxType> witness;
        std::optional<SizeType> backjump_level;

        static SearchOutcome exhausted() { return {}; }
        static SearchOutcome found(UpperBoxType const& witness) {
            SearchOutcome result;
            result.witness=witness;
            return result;
        }
        static SearchOutcome backjump(SizeType level) {
            SearchOutcome result;
            result.backjump_level=level;
            return result;
        }
    };

    Bool _literal_true(Int literal) const
    {
        SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
        int8_t value=_assignment[variable].value;
        ARIADNE_ASSERT(value>=0);
        return literal>0 ? value==1 : value==0;
    }

    Bool _assign_literal(Int literal, std::optional<SizeType> reason_clause = std::nullopt)
    {
        SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
        int8_t value=literal>0 ? 1 : 0;
        AssignmentInfo& assignment=_assignment[variable];
        if(assignment.value<0) {
            assignment.value=value;
            assignment.decision_level=this->_decision_level();
            assignment.reason_clause=reason_clause;
            _trail.push_back(variable);
            return true;
        }
        return assignment.value==value;
    }

    SizeType _original_clause_count() const
    {
        return _encoding.clauses().size();
    }

    SizeType _clause_count() const
    {
        return this->_original_clause_count()+_learned_clauses.size();
    }

    SmtBooleanEncoding::Clause const& _clause(SizeType index) const
    {
        if(index<this->_original_clause_count()) {
            return _encoding.clauses()[index];
        }
        return _learned_clauses[index-this->_original_clause_count()];
    }

    Bool _is_learned_clause(SizeType index) const
    {
        return index>=this->_original_clause_count();
    }

    SizeType _add_learned_clause(std::vector<Int> const& clause, Bool theory_clause = false)
    {
        _learned_clauses.emplace_back(clause.begin(),clause.end());
        _learned_clause_is_theory.push_back(theory_clause);
        _learned_clause_active.push_back(true);
        _learned_clause_activity.push_back(1u);
        _learned_clause_generation.push_back(_statistics.learned_clauses);
        ++_statistics.learned_clauses;
        if(theory_clause) {
            ++_statistics.theory_learned_clauses;
            _statistics.theory_learned_clause_literals+=clause.size();
        } else {
            _statistics.peak_active_non_theory_learned_clauses=std::max(
                _statistics.peak_active_non_theory_learned_clauses,
                this->_active_non_theory_learned_clause_count());
        }
        return this->_original_clause_count()+_learned_clauses.size()-1u;
    }

    Bool _is_theory_learned_clause(SizeType index) const
    {
        if(not this->_is_learned_clause(index)) {
            return false;
        }
        return _learned_clause_is_theory[index-this->_original_clause_count()];
    }

    Bool _is_active_clause(SizeType index) const
    {
        if(not this->_is_learned_clause(index)) {
            return true;
        }
        return _learned_clause_active[index-this->_original_clause_count()];
    }

    SizeType _active_non_theory_learned_clause_count() const
    {
        SizeType count=0u;
        for(SizeType i=0u; i<_learned_clauses.size(); ++i) {
            if(_learned_clause_active[i] && not _learned_clause_is_theory[i]) {
                ++count;
            }
        }
        return count;
    }

    Bool _learned_clause_locked(SizeType index) const
    {
        if(not this->_is_learned_clause(index)) {
            return true;
        }
        for(AssignmentInfo const& assignment:_assignment) {
            if(assignment.value>=0
               && assignment.reason_clause.has_value()
               && *assignment.reason_clause==index) {
                return true;
            }
        }
        return false;
    }

    Void _bump_learned_clause_activity(SizeType index)
    {
        if(not this->_is_learned_clause(index)) {
            return;
        }
        SizeType learned_index=index-this->_original_clause_count();
        if(not _learned_clause_active[learned_index]) {
            return;
        }
        ++_learned_clause_activity[learned_index];
        ++_statistics.learned_clause_activity_bumps;
    }

    Void _maybe_prune_learned_clauses(std::optional<SizeType> protected_clause=std::nullopt)
    {
        SizeType const limit=_solver.configuration().learned_clause_limit();
        SizeType active=this->_active_non_theory_learned_clause_count();
        if(active<=limit) {
            return;
        }

        ++_statistics.learned_clause_pruning_runs;
        std::vector<SizeType> candidates;
        for(SizeType i=0u; i<_learned_clauses.size(); ++i) {
            SizeType clause_index=this->_original_clause_count()+i;
            SizeType const current_generation=_statistics.learned_clauses;
            SizeType const clause_generation=_learned_clause_generation[i];
            Bool const recent=(current_generation<=clause_generation+2u);
            Bool const short_clause=(_learned_clauses[i].size()<=2u);
            Bool const useful=(_learned_clause_activity[i]>1u);
            if(not _learned_clause_active[i]
               || _learned_clause_is_theory[i]
               || recent
               || short_clause
               || useful
               || (protected_clause.has_value() && clause_index==*protected_clause)
               || this->_learned_clause_locked(clause_index)) {
                continue;
            }
            candidates.push_back(clause_index);
        }

        std::stable_sort(candidates.begin(),candidates.end(),[this](SizeType lhs, SizeType rhs) {
            SizeType li=lhs-this->_original_clause_count();
            SizeType ri=rhs-this->_original_clause_count();
            if(_learned_clause_activity[li]!=_learned_clause_activity[ri]) {
                return _learned_clause_activity[li]<_learned_clause_activity[ri];
            }
            return _learned_clauses[li].size()>_learned_clauses[ri].size();
        });

        for(SizeType clause_index:candidates) {
            if(active<=limit) {
                break;
            }
            SizeType learned_index=clause_index-this->_original_clause_count();
            _learned_clause_active[learned_index]=false;
            --active;
            ++_statistics.learned_clauses_pruned;
        }
    }

    Bool _unit_propagate()
    {
        _last_boolean_conflict_clause.reset();
        Bool changed=true;
        while(changed) {
            changed=false;
            for(SizeType clause_index=0u; clause_index<this->_clause_count(); ++clause_index) {
                if(not this->_is_active_clause(clause_index)) {
                    continue;
                }
                auto const& clause=this->_clause(clause_index);
                Bool satisfied=false;
                SizeType unassigned_count=0u;
                Int unit_literal=0;

                for(Int literal:clause) {
                    SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
                    int8_t value=_assignment[variable].value;
                    if(value<0) {
                        ++unassigned_count;
                        unit_literal=literal;
                    } else if(this->_literal_true(literal)) {
                        satisfied=true;
                        break;
                    }
                }

                if(satisfied) {
                    continue;
                }

                if(unassigned_count==0u) {
                    ++_statistics.boolean_conflicts;
                    this->_bump_learned_clause_activity(clause_index);
                    _last_boolean_conflict_clause=clause_index;
                    return false;
                }

                if(unassigned_count==1u) {
                    SizeType variable=static_cast<SizeType>(unit_literal>0 ? unit_literal : -unit_literal);
                    if(_assignment[variable].value<0) {
                        Bool assigned=this->_assign_literal(unit_literal,clause_index);
                        ARIADNE_ASSERT(assigned);
                        ARIADNE_ASSERT(_assignment[variable].reason_clause.has_value());
                        ARIADNE_ASSERT(*_assignment[variable].reason_clause==clause_index);
                        ++_statistics.boolean_propagations;
                        ++_statistics.boolean_reasoned_propagations;
                        if(this->_is_learned_clause(clause_index)) {
                            this->_bump_learned_clause_activity(clause_index);
                            ++_statistics.learned_clause_propagations;
                            if(this->_is_theory_learned_clause(clause_index)) {
                                ++_statistics.theory_learned_clause_propagations;
                            }
                        }
                        changed=true;
                    } else if(not this->_literal_true(unit_literal)) {
                        ++_statistics.boolean_conflicts;
                        this->_bump_learned_clause_activity(clause_index);
                        _last_boolean_conflict_clause=clause_index;
                        return false;
                    }
                }
            }
        }
        return true;
    }

    Bool _clause_contains_variable(std::vector<Int> const& clause, SizeType variable) const
    {
        for(Int literal:clause) {
            SizeType literal_variable=static_cast<SizeType>(literal>0 ? literal : -literal);
            if(literal_variable==variable) {
                return true;
            }
        }
        return false;
    }

    SizeType _current_level_literal_count(std::vector<Int> const& clause) const
    {
        SizeType count=0u;
        for(Int literal:clause) {
            SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
            if(_assignment[variable].decision_level==this->_decision_level()) {
                ++count;
            }
        }
        return count;
    }

    std::vector<Int> _resolve_on_variable(
        std::vector<Int> const& lhs,
        SmtBooleanEncoding::Clause const& rhs,
        SizeType variable) const
    {
        std::vector<Int> result;
        result.reserve(lhs.size()+rhs.size());

        auto append_unique=[&](Int literal) {
            SizeType literal_variable=static_cast<SizeType>(literal>0 ? literal : -literal);
            if(literal_variable==variable) {
                return;
            }
            for(Int existing:result) {
                if(existing==literal) {
                    return;
                }
            }
            result.push_back(literal);
        };

        for(Int literal:lhs) { append_unique(literal); }
        for(Int literal:rhs) { append_unique(literal); }
        return result;
    }

    ConflictAnalysis _analyze_boolean_conflict(SizeType conflict_clause_index)
    {
        ARIADNE_ASSERT(conflict_clause_index<this->_clause_count());
        this->_bump_learned_clause_activity(conflict_clause_index);
        ConflictAnalysis analysis;
        auto const& conflict_clause=this->_clause(conflict_clause_index);
        analysis.learned_clause.assign(conflict_clause.begin(),conflict_clause.end());

        while(this->_current_level_literal_count(analysis.learned_clause)>1u) {
            std::optional<SizeType> pivot;
            for(auto iter=_trail.rbegin(); iter!=_trail.rend(); ++iter) {
                SizeType variable=*iter;
                AssignmentInfo const& assignment=_assignment[variable];
                if(assignment.decision_level==this->_decision_level()
                   && assignment.reason_clause.has_value()
                   && this->_clause_contains_variable(analysis.learned_clause,variable)) {
                    pivot=variable;
                    break;
                }
            }

            ARIADNE_ASSERT(pivot.has_value());
            SizeType reason_index=*_assignment[*pivot].reason_clause;
            ARIADNE_ASSERT(reason_index<this->_clause_count());
            this->_bump_learned_clause_activity(reason_index);
            analysis.learned_clause=this->_resolve_on_variable(
                analysis.learned_clause,this->_clause(reason_index),*pivot);
        }

        SizeType current_level=this->_decision_level();
        SizeType backjump_level=0u;
        for(Int literal:analysis.learned_clause) {
            SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
            SizeType level=_assignment[variable].decision_level;
            if(level!=current_level) {
                backjump_level=std::max(backjump_level,level);
            }
        }
        analysis.backjump_level=backjump_level;
        ARIADNE_ASSERT(this->_current_level_literal_count(analysis.learned_clause)==1u);
        return analysis;
    }

    SizeType _next_unassigned_variable() const
    {
        for(SizeType variable=1u; variable<=_encoding.variable_count(); ++variable) {
            if(_assignment[variable].value<0) {
                return variable;
            }
        }
        return 0u;
    }

    SizeType _decision_level() const
    {
        return _decision_level_markers.size()-1u;
    }

    Void _push_decision_level()
    {
        _decision_level_markers.push_back(_trail.size());
        _statistics.max_decision_level=std::max(
            _statistics.max_decision_level,this->_decision_level());
    }

    Void _backtrack_to_level(SizeType level)
    {
        ARIADNE_ASSERT(level<this->_decision_level_markers.size());
        SizeType marker=_decision_level_markers[level];
        while(_trail.size()>marker) {
            SizeType variable=_trail.back();
            _trail.pop_back();
            _assignment[variable]=AssignmentInfo();
        }
        _decision_level_markers.resize(level+1u);
    }

    SearchOutcome _search_boolean()
    {
        SizeType const box_limit=_solver.configuration().box_processing_limit();
        if(_statistics.boxes_processed>=box_limit) {
            _theory_unknown_seen=true;
            return SearchOutcome::exhausted();
        }

        if(not this->_unit_propagate()) {
            if(not _last_boolean_conflict_clause.has_value()) {
                return SearchOutcome::exhausted();
            }

            if(this->_decision_level()==0u) {
                _last_boolean_conflict_clause.reset();
                return SearchOutcome::exhausted();
            }

            SizeType conflict_level=this->_decision_level();
            ConflictAnalysis analysis=this->_analyze_boolean_conflict(
                *_last_boolean_conflict_clause);
            ++_statistics.boolean_conflicts_analyzed;
            _statistics.learned_clause_literals+=analysis.learned_clause.size();
            _statistics.last_learned_clause_literals=analysis.learned_clause.size();
            _statistics.last_learned_current_level_literals=
                this->_current_level_literal_count(analysis.learned_clause);
            _statistics.last_backjump_level=analysis.backjump_level;

            SizeType learned_index=this->_add_learned_clause(analysis.learned_clause);
            _last_boolean_conflict_clause.reset();

            ARIADNE_ASSERT(analysis.backjump_level<conflict_level);
            if(analysis.backjump_level+1u<conflict_level) {
                ++_statistics.nonchronological_backjumps;
            }
            this->_backtrack_to_level(analysis.backjump_level);
            ++_statistics.boolean_backtracks;
            this->_maybe_prune_learned_clauses(learned_index);
            return SearchOutcome::backjump(analysis.backjump_level);
        }

        SizeType variable=this->_next_unassigned_variable();
        if(variable==0u) {
            std::optional<UpperBoxType> witness=this->_check_theory_assignment();
            if(witness.has_value()) {
                return SearchOutcome::found(*witness);
            }
            return SearchOutcome::exhausted();
        }

        if(not this->_check_partial_theory_consistency()) {
            if(not _last_theory_conflict_clause.has_value()) {
                return SearchOutcome::exhausted();
            }

            if(this->_decision_level()==0u) {
                _last_theory_conflict_clause.reset();
                return SearchOutcome::exhausted();
            }

            SizeType conflict_level=this->_decision_level();
            ConflictAnalysis analysis=this->_analyze_boolean_conflict(
                *_last_theory_conflict_clause);
            ++_statistics.boolean_conflicts_analyzed;
            _statistics.learned_clause_literals+=analysis.learned_clause.size();
            _statistics.last_learned_clause_literals=analysis.learned_clause.size();
            _statistics.last_learned_current_level_literals=
                this->_current_level_literal_count(analysis.learned_clause);
            _statistics.last_backjump_level=analysis.backjump_level;

            SizeType learned_index=this->_add_learned_clause(analysis.learned_clause);
            _last_theory_conflict_clause.reset();

            ARIADNE_ASSERT(analysis.backjump_level<conflict_level);
            if(analysis.backjump_level+1u<conflict_level) {
                ++_statistics.nonchronological_backjumps;
            }
            this->_backtrack_to_level(analysis.backjump_level);
            ++_statistics.boolean_backtracks;
            this->_maybe_prune_learned_clauses(learned_index);
            return SearchOutcome::backjump(analysis.backjump_level);
        }

        ++_statistics.boolean_decisions;
        SizeType parent_level=this->_decision_level();

        this->_push_decision_level();
        ARIADNE_ASSERT(this->_assign_literal(-static_cast<Int>(variable)));
        ARIADNE_ASSERT(_assignment[variable].decision_level==this->_decision_level());
        ARIADNE_ASSERT(not _assignment[variable].reason_clause.has_value());

        SearchOutcome first=this->_search_boolean();
        if(first.witness.has_value()) {
            return first;
        }
        if(first.backjump_level.has_value()) {
            if(*first.backjump_level<parent_level) {
                return first;
            }
            ARIADNE_ASSERT(*first.backjump_level==parent_level);
            return this->_search_boolean();
        }

        this->_backtrack_to_level(parent_level);
        ++_statistics.boolean_backtracks;

        this->_push_decision_level();
        ARIADNE_ASSERT(this->_assign_literal(static_cast<Int>(variable)));
        ARIADNE_ASSERT(_assignment[variable].decision_level==this->_decision_level());
        ARIADNE_ASSERT(not _assignment[variable].reason_clause.has_value());

        SearchOutcome second=this->_search_boolean();
        if(second.witness.has_value()) {
            return second;
        }
        if(second.backjump_level.has_value()) {
            if(*second.backjump_level<parent_level) {
                return second;
            }
            ARIADNE_ASSERT(*second.backjump_level==parent_level);
            return this->_search_boolean();
        }

        this->_backtrack_to_level(parent_level);
        ++_statistics.boolean_backtracks;
        return SearchOutcome::exhausted();
    }

    std::vector<Int> _current_theory_nogood() const
    {
        std::vector<Int> clause;
        clause.reserve(_encoding.atom_count());

        for(SizeType i=0u; i!=_encoding.atom_count(); ++i) {
            SizeType variable=_encoding.atom_variable(i);
            AssignmentInfo const& assignment=_assignment[variable];
            if(assignment.value<0) {
                continue;
            }

            Int literal=static_cast<Int>(variable);
            if(assignment.value==1) {
                literal=-literal;
            }

            Bool duplicate=false;
            for(Int existing:clause) {
                if(existing==literal) {
                    duplicate=true;
                    break;
                }
            }
            if(not duplicate) {
                clause.push_back(literal);
            }
        }

        return clause;
    }

    std::vector<SmtTheoryAlternatives> _theory_alternatives_for_nogood(
        std::vector<Int> const& clause) const
    {
        std::vector<SmtTheoryAlternatives> alternatives;
        alternatives.reserve(clause.size());

        for(Int nogood_literal:clause) {
            SizeType variable=static_cast<SizeType>(
                nogood_literal>0 ? nogood_literal : -nogood_literal);

            std::optional<SizeType> atom_index;
            for(SizeType i=0u; i!=_encoding.atom_count(); ++i) {
                if(_encoding.atom_variable(i)==variable) {
                    atom_index=i;
                    break;
                }
            }
            ARIADNE_ASSERT(atom_index.has_value());

            SmtTheoryLiteral literal=make_smt_theory_literal(_encoding.atom(*atom_index));
            Bool assignment_value=nogood_literal<0;
            if(not assignment_value) {
                literal=literal.negated();
            }
            alternatives.push_back(normalize_smt_theory_literal(literal));
        }

        return alternatives;
    }

    Bool _nogood_theory_consistent(std::vector<Int> const& clause)
    {
        if(clause.empty()) {
            return true;
        }

        ++_statistics.theory_minimization_checks;
        std::vector<SmtTheoryAlternatives> alternatives=
            this->_theory_alternatives_for_nogood(clause);
        List<SmtTheoryPrimitiveLiteral> literals;
        return this->_theory_alternatives_consistent(alternatives,0u,literals);
    }

    std::vector<Int> _minimize_theory_nogood(std::vector<Int> clause)
    {
        std::vector<SizeType> trail_rank(_assignment.size(),0u);
        for(SizeType rank=0u; rank<_trail.size(); ++rank) {
            trail_rank[_trail[rank]]=rank+1u;
        }

        std::stable_sort(clause.begin(),clause.end(),[this,&trail_rank](Int lhs, Int rhs) {
            SizeType lhs_variable=static_cast<SizeType>(lhs>0 ? lhs : -lhs);
            SizeType rhs_variable=static_cast<SizeType>(rhs>0 ? rhs : -rhs);
            SizeType lhs_level=_assignment[lhs_variable].decision_level;
            SizeType rhs_level=_assignment[rhs_variable].decision_level;
            if(lhs_level!=rhs_level) {
                return lhs_level>rhs_level;
            }
            return trail_rank[lhs_variable]>trail_rank[rhs_variable];
        });

        if(not clause.empty()
           && _statistics.first_minimization_candidate_trail_rank==0u) {
            SizeType first_variable=static_cast<SizeType>(
                clause.front()>0 ? clause.front() : -clause.front());
            _statistics.first_minimization_candidate_trail_rank=
                trail_rank[first_variable];
        }

        SizeType const budget=_solver.configuration().theory_minimization_budget();
        SizeType checks=0u;
        SizeType i=0u;
        while(i<clause.size()) {
            if(checks>=budget) {
                ++_statistics.theory_minimization_budget_exhaustions;
                break;
            }

            std::vector<Int> candidate=clause;
            candidate.erase(candidate.begin()+static_cast<std::ptrdiff_t>(i));
            if(candidate.empty()) {
                ++i;
                continue;
            }

            ++checks;
            if(not this->_nogood_theory_consistent(candidate)) {
                clause=std::move(candidate);
            } else {
                ++i;
            }
        }
        return clause;
    }

    SizeType _learn_current_theory_nogood()
    {
        std::vector<Int> clause=this->_current_theory_nogood();
        ARIADNE_ASSERT(not clause.empty());

        SizeType raw_size=clause.size();
        _statistics.theory_nogood_raw_literals+=raw_size;
        clause=this->_minimize_theory_nogood(std::move(clause));
        ARIADNE_ASSERT(not clause.empty());
        _statistics.theory_nogood_minimized_literals+=clause.size();
        _statistics.theory_nogood_literals_removed+=raw_size-clause.size();

        return this->_add_learned_clause(clause,true);
    }

    Bool _check_partial_theory_consistency()
    {
        _last_theory_conflict_clause.reset();
        std::vector<SmtTheoryAlternatives> alternatives;
        alternatives.reserve(_encoding.atom_count());

        for(SizeType i=0; i!=_encoding.atom_count(); ++i) {
            SizeType variable=_encoding.atom_variable(i);
            if(_assignment[variable].value<0) {
                continue;
            }

            SmtTheoryLiteral literal=make_smt_theory_literal(_encoding.atom(i));
            if(_assignment[variable].value==0) {
                literal=literal.negated();
            }
            alternatives.push_back(normalize_smt_theory_literal(literal));
        }

        if(alternatives.empty()) {
            return true;
        }

        ++_statistics.theory_checks;
        List<SmtTheoryPrimitiveLiteral> literals;
        Bool consistent=this->_theory_alternatives_consistent(alternatives,0u,literals);
        if(not consistent) {
            ++_statistics.theory_conflicts;
            _last_theory_conflict_clause=this->_learn_current_theory_nogood();
        }
        return consistent;
    }

    SmtResult _solve_theory_literals(List<SmtTheoryPrimitiveLiteral> const& literals)
    {
        SizeType const limit=_solver.configuration().box_processing_limit();
        SizeType const processed=_statistics.boxes_processed;
        SizeType const remaining=processed>=limit ? 0u : limit-processed;
        SmtSolver theory_solver(SmtSolverConfiguration(
            _solver.configuration().epsilon(),
            _solver.configuration().theory_minimization_budget(),
            _solver.configuration().learned_clause_limit(),
            remaining));
        return _parallel
            ? theory_solver.solve_parallel(_space,_domain,literals)
            : theory_solver.solve(_space,_domain,literals);
    }

    Bool _theory_alternatives_consistent(
        std::vector<SmtTheoryAlternatives> const& alternatives,
        SizeType atom,
        List<SmtTheoryPrimitiveLiteral>& literals)
    {
        if(atom==alternatives.size()) {
            SmtResult result=this->_solve_theory_literals(literals);
            add_statistics(_statistics,result.statistics());
            if(result.is_unknown()) {
                _theory_unknown_seen=true;
                return true;
            }
            return result.is_epsilon_sat();
        }

        for(auto const& alternative:alternatives[atom]) {
            SizeType old_size=literals.size();
            for(auto const& primitive:alternative) {
                literals.append(primitive);
            }

            Bool consistent=this->_theory_alternatives_consistent(
                alternatives,atom+1u,literals);
            literals.erase(
                literals.begin()+static_cast<std::ptrdiff_t>(old_size),
                literals.end());

            if(consistent) {
                return true;
            }
        }
        return false;
    }

    std::optional<UpperBoxType> _check_theory_assignment()
    {
        ++_statistics.theory_checks;

        std::vector<SmtTheoryAlternatives> alternatives;
        alternatives.reserve(_encoding.atom_count());

        for(SizeType i=0; i!=_encoding.atom_count(); ++i) {
            SizeType variable=_encoding.atom_variable(i);
            ARIADNE_ASSERT(_assignment[variable].value>=0);

            SmtTheoryLiteral literal=make_smt_theory_literal(_encoding.atom(i));
            if(_assignment[variable].value==0) {
                literal=literal.negated();
            }
            alternatives.push_back(normalize_smt_theory_literal(literal));
        }

        List<SmtTheoryPrimitiveLiteral> literals;
        return this->_search_theory_alternatives(alternatives,0u,literals);
    }

    std::optional<UpperBoxType> _search_theory_alternatives(
        std::vector<SmtTheoryAlternatives> const& alternatives,
        SizeType atom,
        List<SmtTheoryPrimitiveLiteral>& literals)
    {
        if(atom==alternatives.size()) {
            SmtResult result=this->_solve_theory_literals(literals);
            add_statistics(_statistics,result.statistics());
            if(result.is_epsilon_sat()) {
                return result.witness();
            }
            if(result.is_unknown()) {
                _theory_unknown_seen=true;
            }
            return std::nullopt;
        }

        for(auto const& alternative:alternatives[atom]) {
            SizeType old_size=literals.size();
            for(auto const& primitive:alternative) {
                literals.append(primitive);
            }

            if(auto witness=this->_search_theory_alternatives(
                    alternatives,atom+1u,literals); witness.has_value()) {
                literals.erase(literals.begin()+static_cast<std::ptrdiff_t>(old_size),literals.end());
                return witness;
            }
            literals.erase(literals.begin()+static_cast<std::ptrdiff_t>(old_size),literals.end());
        }
        return std::nullopt;
    }

    SmtSolver const& _solver;
    RealSpace const& _space;
    ExactBoxType const& _domain;
    SmtBooleanEncoding const& _encoding;
    Bool _parallel;
    std::vector<AssignmentInfo> _assignment;
    std::vector<SizeType> _trail;
    std::vector<SizeType> _decision_level_markers;
    std::vector<SmtBooleanEncoding::Clause> _learned_clauses;
    std::vector<Bool> _learned_clause_is_theory;
    std::vector<Bool> _learned_clause_active;
    std::vector<SizeType> _learned_clause_activity;
    std::vector<SizeType> _learned_clause_generation;
    std::optional<SizeType> _last_boolean_conflict_clause;
    std::optional<SizeType> _last_theory_conflict_clause;
    Bool _theory_unknown_seen=false;
    SmtSearchStatistics _statistics;
};
} // namespace

SmtResult SmtSolver::solve(RealSpace const& space,
                           ExactBoxType const& domain,
                           ContinuousPredicate const& predicate) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    ARIADNE_PRECONDITION(space.size()==domain.dimension());

    SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(predicate);
    return SmtDpllSearch(*this,space,domain,encoding,false).solve();
}

SmtResult SmtSolver::solve_parallel(RealSpace const& space,
                                    ExactBoxType const& domain,
                                    ContinuousPredicate const& predicate) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    ARIADNE_PRECONDITION(space.size()==domain.dimension());

    SmtBooleanEncoding encoding=SmtBooleanEncoder().encode(predicate);
    return SmtDpllSearch(*this,space,domain,encoding,true).solve();
}

} // namespace Ariadne
