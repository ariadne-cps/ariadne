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
#include <set>
#include <thread>

#include "betterthreads/workload.hpp"

#include "solvers/constraint_solver.hpp"
#include "solvers/nonlinear_programming.hpp"
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

Pair<SizeType,Pair<Bool,Bool>> sensitivity_split_coordinate(
    UpperBoxType const& domain,
    std::vector<ValidatedScalarMultivariateFunction> const& functions)
{
    auto widths=domain.widths();

    SizeType geometric=0u;
    for(SizeType variable=1u; variable!=domain.dimension(); ++variable) {
        if(widths[variable].raw()>widths[geometric].raw()) {
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
                sensitivity+=candidate;
            }
        }
        if(active && (not selected_score.has_value()
                      || sensitivity.raw()>selected_score->raw())) {
            selected=variable;
            selected_score=sensitivity;
        }
    }

    Bool guided=selected.has_value();
    SizeType coordinate=guided ? *selected : geometric;
    Bool overrode=guided && coordinate!=geometric;
    return {coordinate,{guided,overrode}};
}

UpperBoxType singleton_box(
    ConstraintSolverInterface::ExactPointType const& point)
{
    return UpperBoxType(point.dimension(),[&](SizeType i) {
        return UpperIntervalType(ExactIntervalType(point[i],point[i]));
    });
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
    SizeType learned_clause_limit, SizeType box_processing_limit,
    Bool candidate_search_enabled)
    : _epsilon(epsilon),
      _theory_minimization_budget(theory_minimization_budget),
      _learned_clause_limit(learned_clause_limit),
      _box_processing_limit(box_processing_limit),
      _candidate_search_enabled(candidate_search_enabled)
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

ExactIntervalType SmtSolver::_original_bounds(ValidatedConstraint const& constraint) const
{
    return constraint.bounds();
}

ExactIntervalType SmtSolver::_original_bounds(SmtTheoryPrimitiveRelation relation) const
{
    SmtSolverTestSupport::validate_primitive_relation(relation);
    if(relation==SmtTheoryPrimitiveRelation::EQ_ZERO) {
        return ExactIntervalType(0,0);
    }
    return ExactIntervalType(0,+infty);
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
    SmtSolverTestSupport::validate_primitive_relation(relation);
    FloatDP epsilon(_configuration.epsilon(),dp);
    if(relation==SmtTheoryPrimitiveRelation::EQ_ZERO) {
        return ExactIntervalType(-epsilon,+epsilon);
    }
    return ExactIntervalType(-epsilon,+infty);
}
ExactIntervalType
SmtSolver::_epsilon_bounds(CompiledTheoryLiteral const& literal) const
{
    return this->_epsilon_bounds(literal.relation);
}


Bool SmtSolver::_original_reduce(UpperBoxType& domain,
                                 List<ValidatedConstraint> const& constraints,
                                 ReductionStatistics& statistics) const
{
    ConstraintSolver contractor;
    while(true) {
        UpperBoxType previous=domain;
        ++statistics.hull_rounds;
        for(SizeType i=0; i!=constraints.size(); ++i) {
            if(contractor.hull_reduce(
                    domain,constraints[i].function(),this->_original_bounds(constraints[i]))) {
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
                            this->_original_bounds(constraints[i]),
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
                    if(definitely(disjoint(image,this->_original_bounds(constraints[i])))) {
                        return true;
                    }
                }
                return false;
            }
            continue;
        }

        for(SizeType i=0; i!=constraints.size(); ++i) {
            UpperIntervalType image=apply(constraints[i].function(),domain);
            if(definitely(disjoint(image,this->_original_bounds(constraints[i])))) {
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

Bool SmtSolver::_epsilon_overlaps(
    UpperBoxType const& domain,
    List<ValidatedConstraint> const& constraints) const
{
    UpperBoxType point(domain.dimension(),[&](SizeType i) {
        auto m=domain[i].midpoint();
        return UpperIntervalType(m,m);
    });
    for(SizeType i=0u; i!=constraints.size(); ++i) {
        UpperIntervalType image=apply(constraints[i].function(),point);
        UpperIntervalType bounds(this->_epsilon_bounds(constraints[i]));
        if(definitely(disjoint(image,bounds))) {
            return false;
        }
    }
    return true;
}

SmtSolver::CompiledTheoryLiterals
SmtSolver::_compile_theory_literals(RealSpace const& space,
                                    List<SmtTheoryPrimitiveLiteral> const& literals) const
{
    CompiledTheoryLiterals result;
    result.reserve(literals.size());
    for(SizeType i=0; i!=literals.size(); ++i) {
        RealExpression expression=simplify(literals[i].expression());
        if(is_constant(expression,Real(0))) {
            continue;
        }
        result.push_back({
            ValidatedScalarMultivariateFunction(space,expression),
            literals[i].relation()
        });
    }
    return result;
}

Bool SmtSolver::_original_reduce(UpperBoxType& domain,
                                 CompiledTheoryLiterals const& literals,
                                 ReductionStatistics& statistics) const
{
    ConstraintSolver contractor;
    while(true) {
        UpperBoxType previous=domain;
        ++statistics.hull_rounds;
        for(auto const& literal:literals) {
            if(contractor.hull_reduce(
                    domain,literal.function,this->_original_bounds(literal.relation))) {
                return true;
            }
            if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
                UpperIntervalType image=apply(literal.function,domain);
                if(definitely(image.upper_bound()<=0)) {
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
                            this->_original_bounds(literal.relation),
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
                    SmtSolverTestSupport::validate_primitive_relation(literal.relation);
                    if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
                        if(definitely(image.upper_bound()<=0)) {
                            return true;
                        }
                    } else if(definitely(disjoint(
                            image,this->_original_bounds(literal.relation)))) {
                        return true;
                    }
                }
                return false;
            }
            continue;
        }

        for(auto const& literal:literals) {
            UpperIntervalType image=apply(literal.function,domain);
            SmtSolverTestSupport::validate_primitive_relation(literal.relation);
            if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
                if(definitely(image.upper_bound()<=0)) {
                    return true;
                }
            } else if(definitely(disjoint(
                    image,this->_original_bounds(literal.relation)))) {
                return true;
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
        SmtSolverTestSupport::validate_primitive_relation(literal.relation);
        if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
            if(not definitely(image.lower_bound()>-epsilon)) {
                return false;
            }
        } else if(not definitely(
                subset(image,this->_epsilon_bounds(literal.relation)))) {
            return false;
        }
    }
    return true;
}

Bool SmtSolver::_epsilon_overlaps(
    UpperBoxType const& domain,
    CompiledTheoryLiterals const& literals) const
{
    UpperBoxType point(domain.dimension(),[&](SizeType i) {
        auto m=domain[i].midpoint();
        return UpperIntervalType(m,m);
    });
    FloatDP epsilon(_configuration.epsilon(),dp);
    for(auto const& literal:literals) {
        UpperIntervalType image=apply(literal.function,point);
        SmtSolverTestSupport::validate_primitive_relation(literal.relation);
        if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
            if(definitely(image.upper_bound()<=-epsilon)) {
                return false;
            }
        } else if(definitely(disjoint(
                image,UpperIntervalType(this->_epsilon_bounds(literal.relation))))) {
            return false;
        }
    }
    return true;
}

ValidatedScalarMultivariateFunction const&
SmtSolver::_function(ValidatedConstraint const& constraint) const
{
    return constraint.function();
}

ValidatedScalarMultivariateFunction const&
SmtSolver::_function(CompiledTheoryLiteral const& literal) const
{
    return literal.function;
}

template<class Conjunction>
std::optional<UpperBoxType>
SmtSolver::_epsilon_witness(
    UpperBoxType const& domain,
    Conjunction const& conjunction) const
{
    for(UpperBoxType const& candidate:epsilon_witness_candidates(domain)) {
        if(this->_epsilon_satisfied(candidate,conjunction)) {
            return candidate;
        }
    }
    return std::nullopt;
}

template<class Conjunction>
std::optional<UpperBoxType>
SmtSolver::_epsilon_candidate_witness(
    UpperBoxType const& domain,
    Conjunction const& conjunction) const
{
    if(conjunction.empty()) {
        return std::nullopt;
    }

    ValidatedVectorMultivariateFunction function(
        conjunction.size(),this->_function(conjunction[0]).domain());
    ExactBoxType codomain(conjunction.size());
    for(SizeType i=0u; i!=conjunction.size(); ++i) {
        function[i]=this->_function(conjunction[i]);
        codomain[i]=this->_epsilon_bounds(conjunction[i]);
    }

    NonlinearInfeasibleInteriorPointOptimiser candidate_solver;
    auto result=candidate_solver.feasible_candidate(
        cast_exact_box(domain),function,codomain);
    UpperBoxType witness=singleton_box(cast_exact(result.second));
    if(this->_epsilon_satisfied(witness,conjunction)) {
        return witness;
    }
    return std::nullopt;
}

template<class Conjunction>
Pair<Pair<UpperBoxType,UpperBoxType>,Pair<Bool,Bool>>
SmtSolver::_split_box(
    UpperBoxType const& domain,
    Conjunction const& conjunction) const
{
    std::vector<ValidatedScalarMultivariateFunction> functions;
    functions.reserve(conjunction.size());
    for(auto const& item:conjunction) {
        functions.push_back(this->_function(item));
    }
    auto selection=sensitivity_split_coordinate(domain,functions);
    return {domain.split(selection.first),selection.second};
}

template<class Conjunction>
SmtSolver::BoxProcessingResult
SmtSolver::_process_box(
    UpperBoxType domain,
    Conjunction const& conjunction) const
{
    ReductionStatistics reductions;
    if(this->_original_reduce(domain,conjunction,reductions)) {
        return {BoxProcessingStatus::PRUNED,std::nullopt,std::nullopt,reductions};
    }

    if(this->_epsilon_satisfied(domain,conjunction)) {
        UpperBoxType witness(domain.dimension(),[&](SizeType i) {
            auto m=domain[i].midpoint();
            return UpperIntervalType(m,m);
        });
        BoxProcessingResult result{
            BoxProcessingStatus::EPSILON_SAT,witness,std::nullopt,reductions};
        result.epsilon_box_certification=true;
        return result;
    }

    if(auto witness=this->_epsilon_witness(domain,conjunction); witness.has_value()) {
        return {BoxProcessingStatus::EPSILON_SAT,*witness,std::nullopt,reductions};
    }

    std::optional<UpperBoxType> candidate;
    Bool candidate_certified=false;
    if(_configuration.candidate_search_enabled()) {
        candidate=this->_epsilon_candidate_witness(domain,conjunction);
        candidate_certified=candidate.has_value();
    }
    auto candidate_outcome=SmtSolverTestSupport::candidate_witness_outcome(
        _configuration.candidate_search_enabled(),candidate,candidate_certified);
    if(candidate_outcome.certified) {
        ARIADNE_ASSERT(candidate_outcome.witness.has_value());
        BoxProcessingResult result{
            BoxProcessingStatus::EPSILON_SAT,
            *candidate_outcome.witness,
            std::nullopt,
            reductions};
        result.candidate_witness_search=true;
        result.candidate_witness_success=true;
        return result;
    }

    auto split_result=this->_split_box(domain,conjunction);
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
        BoxProcessingResult result{
            BoxProcessingStatus::UNKNOWN,std::nullopt,std::nullopt,reductions};
        result.candidate_witness_search=candidate_outcome.attempted;
        result.non_splittable_epsilon_overlap=
            this->_epsilon_overlaps(domain,conjunction);
        return result;
    }

    BoxProcessingResult result{
        BoxProcessingStatus::SPLIT,
        std::nullopt,
        children,
        reductions,
        split_result.second.first,
        split_result.second.second};
    result.candidate_witness_search=candidate_outcome.attempted;
    return result;
}

Void
SmtSolver::_accumulate_box_processing_statistics(
    SmtSearchStatistics& statistics,
    BoxProcessingResult const& processing) const
{
    SmtSolverTestSupport::accumulate_box_processing_statistics(
        statistics,{
            processing.status,
            processing.reductions.hull_rounds,
            processing.reductions.hull_effective,
            processing.reductions.shaving_rounds,
            processing.reductions.shaving_effective,
            processing.sensitivity_guided_split,
            processing.sensitivity_overrode_geometric_split,
            processing.epsilon_box_certification,
            processing.candidate_witness_search,
            processing.candidate_witness_success,
            processing.non_splittable_epsilon_overlap
        });
}

template<class Conjunction>
SmtResult
SmtSolver::_solve_sequential_conjunction(
    ExactBoxType const& domain,
    Conjunction const& conjunction) const
{
    SmtSearchStatistics statistics;
    SequentialSmtWorkQueue pending;
    pending.push(UpperBoxType(domain));
    Bool unknown_seen=false;

    while(not pending.empty()) {
        if(statistics.boxes_processed>=_configuration.box_processing_limit()) {
            ++statistics.box_budget_exhaustions;
            unknown_seen=true;
            break;
        }
        UpperBoxType current=pending.pop();
        ++statistics.boxes_processed;

        BoxProcessingResult processing=this->_process_box(
            std::move(current),conjunction);
        this->_accumulate_box_processing_statistics(statistics,processing);

        if(processing.status==BoxProcessingStatus::EPSILON_SAT) {
            ARIADNE_ASSERT(processing.witness.has_value());
            return SmtResult::epsilon_sat(*processing.witness,statistics);
        }
        if(processing.status==BoxProcessingStatus::SPLIT) {
            ARIADNE_ASSERT(processing.children.has_value());
            pending.push(std::move(processing.children->second));
            pending.push(std::move(processing.children->first));
        } else if(processing.status==BoxProcessingStatus::UNKNOWN) {
            unknown_seen=true;
        }
    }

    return unknown_seen
        ? SmtResult::unknown(statistics)
        : SmtResult::unsat(statistics);
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
    if(constraints.empty()) {
        return SmtResult::epsilon_sat(singleton_box(domain.midpoint()));
    }
    return this->_solve_sequential_conjunction(domain,constraints);
}

SmtResult SmtSolver::solve(RealSpace const& space,
                           ExactBoxType const& domain,
                           List<SmtTheoryPrimitiveLiteral> const& literals) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    ARIADNE_PRECONDITION(space.size()==domain.dimension());

    if(domain.is_empty()) {
        return SmtResult::unsat();
    }
    if(literals.empty()) {
        return SmtResult::epsilon_sat(singleton_box(domain.midpoint()));
    }

    CompiledTheoryLiterals compiled=this->_compile_theory_literals(space,literals);
    if(compiled.empty()) {
        return SmtResult::epsilon_sat(singleton_box(domain.midpoint()));
    }
    return this->_solve_sequential_conjunction(domain,compiled);
}


namespace {

std::mutex parallel_execution_observation_mutex;
Bool parallel_execution_observation_enabled=false;
std::thread::id parallel_execution_calling_thread;
std::set<std::thread::id> parallel_execution_threads;

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

template<class Conjunction>
SmtResult
SmtSolver::_solve_parallel_conjunction(
    ExactBoxType const& domain,
    Conjunction const& conjunction) const
{
    auto state=std::make_shared<ParallelSmtSearchState>();
    ParallelSmtWorkload workload(
        [](UpperBoxType const&, std::shared_ptr<ConcLog::ProgressIndicator>) { },
        [this,&conjunction,state](
                ParallelSmtWorkload::Access& access,
                UpperBoxType const& box) {
            SmtSolverTestSupport::record_parallel_processing_thread();
            if(state->found.load() || state->limit_reached.load()) {
                return;
            }

            {
                std::lock_guard<std::mutex> lock(state->mutex);
                if(state->statistics.boxes_processed>=
                        _configuration.box_processing_limit()) {
                    ++state->statistics.box_budget_exhaustions;
                    state->unknown.store(true);
                    state->limit_reached.store(true);
                    return;
                }
                ++state->statistics.boxes_processed;
            }

            BoxProcessingResult processing=this->_process_box(box,conjunction);
            {
                std::lock_guard<std::mutex> lock(state->mutex);
                this->_accumulate_box_processing_statistics(
                    state->statistics,processing);
            }

            if(processing.status==BoxProcessingStatus::PRUNED) {
                return;
            }
            if(processing.status==BoxProcessingStatus::EPSILON_SAT) {
                ARIADNE_ASSERT(processing.witness.has_value());
                bool expected=false;
                if(state->found.compare_exchange_strong(expected,true)) {
                    std::lock_guard<std::mutex> lock(state->mutex);
                    state->witness=*processing.witness;
                }
                return;
            }
            if(processing.status==BoxProcessingStatus::SPLIT) {
                ARIADNE_ASSERT(processing.children.has_value());
                if(not state->found.load()) {
                    access.append(processing.children->first);
                    access.append(processing.children->second);
                }
                return;
            }
            ARIADNE_ASSERT(processing.status==BoxProcessingStatus::UNKNOWN);
            state->unknown.store(true);
            return;
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

SmtResult SmtSolver::solve_parallel(
    ExactBoxType const& domain,
    List<ValidatedConstraint> const& constraints) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    for(SizeType i=0; i!=constraints.size(); ++i) {
        ARIADNE_PRECONDITION(constraints[i].argument_size()==domain.dimension());
    }

    if(domain.is_empty()) {
        return SmtResult::unsat();
    }
    if(constraints.empty()) {
        return SmtResult::epsilon_sat(singleton_box(domain.midpoint()));
    }
    return this->_solve_parallel_conjunction(domain,constraints);
}

SmtResult SmtSolver::solve_parallel(
    RealSpace const& space,
    ExactBoxType const& domain,
    List<SmtTheoryPrimitiveLiteral> const& literals) const
{
    ARIADNE_PRECONDITION(domain.is_bounded());
    ARIADNE_PRECONDITION(space.size()==domain.dimension());

    if(domain.is_empty()) {
        return SmtResult::unsat();
    }
    if(literals.empty()) {
        return SmtResult::epsilon_sat(singleton_box(domain.midpoint()));
    }

    CompiledTheoryLiterals compiled=this->_compile_theory_literals(space,literals);
    if(compiled.empty()) {
        return SmtResult::epsilon_sat(singleton_box(domain.midpoint()));
    }
    return this->_solve_parallel_conjunction(domain,compiled);
}


namespace SmtSolverTestSupport {

Void begin_parallel_execution_observation()
{
    std::lock_guard<std::mutex> lock(parallel_execution_observation_mutex);
    parallel_execution_threads.clear();
    parallel_execution_calling_thread=std::this_thread::get_id();
    parallel_execution_observation_enabled=true;
}

ParallelExecutionObservation end_parallel_execution_observation()
{
    std::lock_guard<std::mutex> lock(parallel_execution_observation_mutex);
    ParallelExecutionObservation result;
    result.observed_thread_count=parallel_execution_threads.size();
    result.calling_thread_observed=
        parallel_execution_threads.find(parallel_execution_calling_thread)
        != parallel_execution_threads.end();
    result.worker_thread_count=result.observed_thread_count
        - (result.calling_thread_observed ? 1u : 0u);
    parallel_execution_observation_enabled=false;
    parallel_execution_threads.clear();
    return result;
}

Void record_parallel_processing_thread()
{
    std::lock_guard<std::mutex> lock(parallel_execution_observation_mutex);
    if(parallel_execution_observation_enabled) {
        parallel_execution_threads.insert(std::this_thread::get_id());
    }
}

Void accumulate_statistics(SmtSearchStatistics& target, SmtSearchStatistics const& source)
{
    target.boxes_processed+=source.boxes_processed;
    target.boxes_pruned+=source.boxes_pruned;
    target.boxes_split+=source.boxes_split;
    target.boxes_unknown+=source.boxes_unknown;
    target.box_budget_exhaustions+=source.box_budget_exhaustions;
    target.non_splittable_uncertified_boxes+=source.non_splittable_uncertified_boxes;
    target.non_splittable_epsilon_overlap_boxes+=source.non_splittable_epsilon_overlap_boxes;
    target.hull_reduction_rounds+=source.hull_reduction_rounds;
    target.hull_effective_reductions+=source.hull_effective_reductions;
    target.shaving_reduction_rounds+=source.shaving_reduction_rounds;
    target.shaving_effective_reductions+=source.shaving_effective_reductions;
    target.sensitivity_guided_splits+=source.sensitivity_guided_splits;
    target.sensitivity_overrides_geometric_splits+=source.sensitivity_overrides_geometric_splits;
    target.epsilon_box_certifications+=source.epsilon_box_certifications;
    target.candidate_witness_searches+=source.candidate_witness_searches;
    target.candidate_witness_successes+=source.candidate_witness_successes;
    target.boolean_decisions+=source.boolean_decisions;
    target.boolean_propagations+=source.boolean_propagations;
    target.boolean_reasoned_propagations+=source.boolean_reasoned_propagations;
    target.boolean_conflicts+=source.boolean_conflicts;
    target.boolean_backtracks+=source.boolean_backtracks;
    target.max_decision_level=std::max(target.max_decision_level,source.max_decision_level);
    target.boolean_conflicts_analyzed+=source.boolean_conflicts_analyzed;
    target.learned_clause_literals+=source.learned_clause_literals;
    if(source.boolean_conflicts_analyzed!=0u) {
        target.last_learned_clause_literals=source.last_learned_clause_literals;
        target.last_learned_current_level_literals=source.last_learned_current_level_literals;
        target.last_backjump_level=source.last_backjump_level;
    }
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


std::vector<SizeType> learned_clause_pruning_candidates(
    std::vector<LearnedClausePruningEntry> const& entries)
{
    std::vector<SizeType> candidates;
    for(SizeType i=0u; i!=entries.size(); ++i) {
        auto const& entry=entries[i];
        if(not entry.active
           || entry.theory
           || entry.recent
           || entry.short_clause
           || entry.useful
           || entry.protected_clause
           || entry.locked) {
            continue;
        }
        candidates.push_back(i);
    }
    std::stable_sort(candidates.begin(),candidates.end(),
        [&entries](SizeType lhs, SizeType rhs) {
            if(entries[lhs].activity!=entries[rhs].activity) {
                return entries[lhs].activity<entries[rhs].activity;
            }
            return entries[lhs].size>entries[rhs].size;
        });
    return candidates;
}

SizeType apply_learned_clause_pruning(
    std::vector<Bool>& active,
    std::vector<SizeType> const& candidates,
    SizeType active_count,
    SizeType limit)
{
    SizeType pruned=0u;
    for(SizeType index:candidates) {
        if(active_count<=limit) {
            break;
        }
        active[index]=false;
        --active_count;
        ++pruned;
    }
    return pruned;
}


CandidateWitnessOutcome candidate_witness_outcome(
    Bool enabled,
    std::optional<UpperBoxType> const& candidate,
    Bool certified)
{
    if(not enabled) {
        return {};
    }
    if(candidate.has_value() && certified) {
        return {true,true,candidate};
    }
    return {true,false,std::nullopt};
}

SearchOutcome SearchOutcome::exhausted()
{
    return {};
}

SearchOutcome SearchOutcome::found(UpperBoxType const& witness)
{
    SearchOutcome result;
    result.witness=witness;
    return result;
}

SearchOutcome SearchOutcome::backjump(SizeType level)
{
    SearchOutcome result;
    result.backjump_level=level;
    return result;
}

SmtResult finalize_search_outcome(
    SearchOutcome const& outcome,
    Bool theory_unknown_seen,
    SmtSearchStatistics const& statistics)
{
    if(outcome.witness.has_value()) {
        return SmtResult::epsilon_sat(*outcome.witness,statistics);
    }
    if(theory_unknown_seen) {
        return SmtResult::unknown(statistics);
    }
    return SmtResult::unsat(statistics);
}

Void validate_primitive_relation(SmtTheoryPrimitiveRelation relation)
{
    switch(relation) {
        case SmtTheoryPrimitiveRelation::EQ_ZERO:
        case SmtTheoryPrimitiveRelation::GEQ_ZERO:
        case SmtTheoryPrimitiveRelation::GT_ZERO:
            return;
        default:
            ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
    }
}

std::vector<Int> resolve_clause_on_variable(
    std::vector<Int> const& lhs,
    std::vector<Int> const& rhs,
    SizeType variable)
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

Void order_theory_nogood(
    std::vector<Int>& clause,
    std::vector<SizeType> const& decision_levels,
    std::vector<SizeType> const& trail_rank)
{
    std::stable_sort(clause.begin(),clause.end(),
        [&decision_levels,&trail_rank](Int lhs, Int rhs) {
            SizeType lhs_variable=static_cast<SizeType>(lhs>0 ? lhs : -lhs);
            SizeType rhs_variable=static_cast<SizeType>(rhs>0 ? rhs : -rhs);
            SizeType lhs_level=decision_levels[lhs_variable];
            SizeType rhs_level=decision_levels[rhs_variable];
            if(lhs_level!=rhs_level) {
                return lhs_level>rhs_level;
            }
            return trail_rank[lhs_variable]>trail_rank[rhs_variable];
        });
}

AssignmentDecision assignment_decision(
    int8_t current_value,
    int8_t requested_value)
{
    if(current_value<0) {
        return {true,true};
    }
    return {current_value==requested_value,false};
}

Bool clause_is_learned(SizeType index, SizeType original_clause_count)
{
    return index>=original_clause_count;
}

Bool learned_clause_is_theory(
    SizeType index,
    SizeType original_clause_count,
    std::vector<Bool> const& theory_flags)
{
    if(not clause_is_learned(index,original_clause_count)) {
        return false;
    }
    return theory_flags[index-original_clause_count];
}

Bool clause_is_active(
    SizeType index,
    SizeType original_clause_count,
    std::vector<Bool> const& active_flags)
{
    if(not clause_is_learned(index,original_clause_count)) {
        return true;
    }
    return active_flags[index-original_clause_count];
}

Bool assignment_locks_clause(
    int8_t assignment_value,
    std::optional<SizeType> const& reason_clause,
    SizeType clause_index)
{
    return assignment_value>=0
        && reason_clause.has_value()
        && *reason_clause==clause_index;
}

Bool should_bump_learned_clause(
    SizeType index,
    SizeType original_clause_count,
    std::vector<Bool> const& active_flags)
{
    if(not clause_is_learned(index,original_clause_count)) {
        return false;
    }
    return active_flags[index-original_clause_count];
}

TheoryResultInterpretation interpret_theory_result(SmtResult const& result)
{
    if(result.is_epsilon_sat()) {
        return {true,false,result.witness()};
    }
    if(result.is_unknown()) {
        return {true,true,std::nullopt};
    }
    return {false,false,std::nullopt};
}

ChildSearchAction classify_child_search_outcome(
    SearchOutcome const& outcome,
    SizeType parent_level,
    Bool has_alternative_branch)
{
    if(outcome.witness.has_value()) {
        return ChildSearchAction::RETURN_OUTCOME;
    }
    if(outcome.backjump_level.has_value()) {
        if(*outcome.backjump_level<parent_level) {
            return ChildSearchAction::RETURN_OUTCOME;
        }
        ARIADNE_ASSERT(*outcome.backjump_level==parent_level);
        return ChildSearchAction::RESTART_AT_PARENT;
    }
    return has_alternative_branch
        ? ChildSearchAction::TRY_ALTERNATIVE
        : ChildSearchAction::EXHAUSTED;
}

ExactIntervalType original_bounds(
    SmtSolver const& solver,
    SmtTheoryPrimitiveRelation relation)
{
    return solver._original_bounds(relation);
}

ExactIntervalType epsilon_bounds(
    SmtSolver const& solver,
    SmtTheoryPrimitiveRelation relation)
{
    return solver._epsilon_bounds(relation);
}

Bool epsilon_satisfied(
    SmtSolver const& solver,
    UpperBoxType const& domain,
    List<ValidatedConstraint> const& constraints)
{
    return solver._epsilon_satisfied(domain,constraints);
}

Bool epsilon_overlaps(
    SmtSolver const& solver,
    UpperBoxType const& domain,
    List<ValidatedConstraint> const& constraints)
{
    return solver._epsilon_overlaps(domain,constraints);
}

std::optional<UpperBoxType> empty_candidate_witness(
    SmtSolver const& solver,
    UpperBoxType const& domain)
{
    List<ValidatedConstraint> constraints;
    return solver._epsilon_candidate_witness(domain,constraints);
}

Void accumulate_box_processing_statistics(
    SmtSearchStatistics& statistics,
    BoxProcessingStatisticsInput const& input)
{
    statistics.hull_reduction_rounds+=input.hull_rounds;
    statistics.hull_effective_reductions+=input.hull_effective;
    statistics.shaving_reduction_rounds+=input.shaving_rounds;
    statistics.shaving_effective_reductions+=input.shaving_effective;

    if(input.sensitivity_guided_split) {
        ++statistics.sensitivity_guided_splits;
    }
    if(input.sensitivity_overrode_geometric_split) {
        ++statistics.sensitivity_overrides_geometric_splits;
    }
    if(input.epsilon_box_certification) {
        ++statistics.epsilon_box_certifications;
    }
    if(input.candidate_witness_search) {
        ++statistics.candidate_witness_searches;
    }
    if(input.candidate_witness_success) {
        ++statistics.candidate_witness_successes;
    }

    switch(input.status) {
        case BoxProcessingStatus::PRUNED:
            ++statistics.boxes_pruned;
            break;
        case BoxProcessingStatus::SPLIT:
            ++statistics.boxes_split;
            break;
        case BoxProcessingStatus::UNKNOWN:
            ++statistics.boxes_unknown;
            ++statistics.non_splittable_uncertified_boxes;
            if(input.non_splittable_epsilon_overlap) {
                ++statistics.non_splittable_epsilon_overlap_boxes;
            }
            break;
        case BoxProcessingStatus::EPSILON_SAT:
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown BoxProcessingStatus");
    }
}

} // namespace SmtSolverTestSupport

namespace {

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
        SmtSolverTestSupport::SearchOutcome outcome=this->_search_boolean();
        return SmtSolverTestSupport::finalize_search_outcome(
            outcome,_theory_unknown_seen,_statistics);
    }

  private:
    using SearchOutcome=SmtSolverTestSupport::SearchOutcome;

    struct AssignmentInfo {
        int8_t value = -1;
        SizeType decision_level = 0u;
        std::optional<SizeType> reason_clause;
    };

    struct ConflictAnalysis {
        std::vector<Int> learned_clause;
        SizeType backjump_level = 0u;
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
        auto decision=SmtSolverTestSupport::assignment_decision(
            assignment.value,value);
        if(decision.newly_assigned) {
            assignment.value=value;
            assignment.decision_level=this->_decision_level();
            assignment.reason_clause=reason_clause;
            _trail.push_back(variable);
        }
        return decision.accepted;
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
        return SmtSolverTestSupport::clause_is_learned(
            index,this->_original_clause_count());
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
        return SmtSolverTestSupport::learned_clause_is_theory(
            index,this->_original_clause_count(),_learned_clause_is_theory);
    }

    Bool _is_active_clause(SizeType index) const
    {
        return SmtSolverTestSupport::clause_is_active(
            index,this->_original_clause_count(),_learned_clause_active);
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
            if(SmtSolverTestSupport::assignment_locks_clause(
                    assignment.value,assignment.reason_clause,index)) {
                return true;
            }
        }
        return false;
    }

    Void _bump_learned_clause_activity(SizeType index)
    {
        if(not SmtSolverTestSupport::should_bump_learned_clause(
                index,this->_original_clause_count(),_learned_clause_active)) {
            return;
        }
        SizeType learned_index=index-this->_original_clause_count();
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
        std::vector<SmtSolverTestSupport::LearnedClausePruningEntry> entries;
        entries.reserve(_learned_clauses.size());
        for(SizeType i=0u; i<_learned_clauses.size(); ++i) {
            SizeType clause_index=this->_original_clause_count()+i;
            SizeType const current_generation=_statistics.learned_clauses;
            SizeType const clause_generation=_learned_clause_generation[i];
            entries.push_back({
                _learned_clause_active[i],
                _learned_clause_is_theory[i],
                current_generation<=clause_generation+2u,
                _learned_clauses[i].size()<=2u,
                _learned_clause_activity[i]>1u,
                protected_clause.has_value() && clause_index==*protected_clause,
                this->_learned_clause_locked(clause_index),
                _learned_clause_activity[i],
                _learned_clauses[i].size()
            });
        }
        std::vector<SizeType> candidates=
            SmtSolverTestSupport::learned_clause_pruning_candidates(entries);

        _statistics.learned_clauses_pruned+=
            SmtSolverTestSupport::apply_learned_clause_pruning(
                _learned_clause_active,candidates,active,limit);
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
        return SmtSolverTestSupport::resolve_clause_on_variable(lhs,rhs,variable);
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
        switch(SmtSolverTestSupport::classify_child_search_outcome(
                   first,parent_level,true)) {
            case SmtSolverTestSupport::ChildSearchAction::RETURN_OUTCOME:
                return first;
            case SmtSolverTestSupport::ChildSearchAction::RESTART_AT_PARENT:
                return this->_search_boolean();
            case SmtSolverTestSupport::ChildSearchAction::TRY_ALTERNATIVE:
                break;
            case SmtSolverTestSupport::ChildSearchAction::EXHAUSTED:
            default:
                ARIADNE_FAIL_MSG("Invalid first child search action");
        }

        this->_backtrack_to_level(parent_level);
        ++_statistics.boolean_backtracks;

        this->_push_decision_level();
        ARIADNE_ASSERT(this->_assign_literal(static_cast<Int>(variable)));
        ARIADNE_ASSERT(_assignment[variable].decision_level==this->_decision_level());
        ARIADNE_ASSERT(not _assignment[variable].reason_clause.has_value());

        SearchOutcome second=this->_search_boolean();
        switch(SmtSolverTestSupport::classify_child_search_outcome(
                   second,parent_level,false)) {
            case SmtSolverTestSupport::ChildSearchAction::RETURN_OUTCOME:
                return second;
            case SmtSolverTestSupport::ChildSearchAction::RESTART_AT_PARENT:
                return this->_search_boolean();
            case SmtSolverTestSupport::ChildSearchAction::EXHAUSTED:
                break;
            case SmtSolverTestSupport::ChildSearchAction::TRY_ALTERNATIVE:
            default:
                ARIADNE_FAIL_MSG("Invalid second child search action");
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

            clause.push_back(literal);
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
        ARIADNE_PRECONDITION(not clause.empty());
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

        std::vector<SizeType> decision_levels(_assignment.size(),0u);
        for(SizeType variable=0u; variable!=_assignment.size(); ++variable) {
            decision_levels[variable]=_assignment[variable].decision_level;
        }
        SmtSolverTestSupport::order_theory_nogood(
            clause,decision_levels,trail_rank);

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
            remaining,
            _solver.configuration().candidate_search_enabled()));
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
            SmtSolverTestSupport::accumulate_statistics(_statistics,result.statistics());
            auto interpretation=SmtSolverTestSupport::interpret_theory_result(result);
            _theory_unknown_seen=_theory_unknown_seen || interpretation.unknown;
            return interpretation.consistent;
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
        if(_encoding.atom_count()==0u) {
            if(_domain.is_empty()) {
                return std::nullopt;
            }
            return singleton_box(_domain.midpoint());
        }

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
            SmtSolverTestSupport::accumulate_statistics(_statistics,result.statistics());
            auto interpretation=SmtSolverTestSupport::interpret_theory_result(result);
            _theory_unknown_seen=_theory_unknown_seen || interpretation.unknown;
            return interpretation.witness;
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
