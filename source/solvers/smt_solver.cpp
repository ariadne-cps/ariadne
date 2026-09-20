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
#include <cstdint>

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
                                CompiledTheoryLiterals const& literals) const
{
    ConstraintSolver contractor;
    FloatDP epsilon(_configuration.epsilon(),dp);
    for(auto const& literal:literals) {
        if(contractor.hull_reduce(domain,literal.function,this->_epsilon_bounds(literal.relation))) {
            return true;
        }
        if(literal.relation==SmtTheoryPrimitiveRelation::GT_ZERO) {
            UpperIntervalType image=apply(literal.function,domain);
            if(definitely(image.upper_bound()<=-epsilon)) {
                return true;
            }
        }
    }
    return definitely(domain.is_empty());
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

SmtSolver::BoxProcessingResult
SmtSolver::_process_box(UpperBoxType domain,
                        CompiledTheoryLiterals const& literals) const
{
    if(this->_epsilon_reduce(domain,literals)) {
        return {BoxProcessingStatus::PRUNED,std::nullopt,std::nullopt};
    }

    if(this->_epsilon_satisfied(domain,literals)) {
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
                       "SMT theory search reached a non-splittable uncertified box: "<<domain);

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

    while(not pending.empty()) {
        UpperBoxType current=pending.pop();
        ++statistics.boxes_processed;

        BoxProcessingResult processing=this->_process_box(std::move(current),compiled);
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
            if(state->found.load()) {
                return;
            }

            BoxProcessingResult processing=this->_process_box(box,compiled);
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


namespace {

Void add_statistics(SmtSearchStatistics& target, SmtSearchStatistics const& source)
{
    target.boxes_processed+=source.boxes_processed;
    target.boxes_pruned+=source.boxes_pruned;
    target.boxes_split+=source.boxes_split;
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
          _assignment(encoding.variable_count()+1u,-1)
    {
    }

    SmtResult solve()
    {
        std::optional<UpperBoxType> witness=this->_search_boolean(1u);
        if(witness.has_value()) {
            return SmtResult::epsilon_sat(*witness,_statistics);
        }
        return SmtResult::unsat(_statistics);
    }

  private:
    Bool _literal_true(Int literal) const
    {
        SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
        int8_t value=_assignment[variable];
        ARIADNE_ASSERT(value>=0);
        return literal>0 ? value==1 : value==0;
    }

    Bool _has_clause_conflict() const
    {
        for(auto const& clause:_encoding.clauses()) {
            Bool satisfied=false;
            Bool undecided=false;
            for(Int literal:clause) {
                SizeType variable=static_cast<SizeType>(literal>0 ? literal : -literal);
                int8_t value=_assignment[variable];
                if(value<0) {
                    undecided=true;
                } else if(this->_literal_true(literal)) {
                    satisfied=true;
                    break;
                }
            }
            if(not satisfied and not undecided) {
                return true;
            }
        }
        return false;
    }

    std::optional<UpperBoxType> _search_boolean(SizeType variable)
    {
        if(this->_has_clause_conflict()) {
            return std::nullopt;
        }

        if(variable>_encoding.variable_count()) {
            return this->_check_theory_assignment();
        }

        _assignment[variable]=0;
        if(auto witness=this->_search_boolean(variable+1u); witness.has_value()) {
            _assignment[variable]=-1;
            return witness;
        }

        _assignment[variable]=1;
        if(auto witness=this->_search_boolean(variable+1u); witness.has_value()) {
            _assignment[variable]=-1;
            return witness;
        }

        _assignment[variable]=-1;
        return std::nullopt;
    }

    std::optional<UpperBoxType> _check_theory_assignment()
    {
        std::vector<SmtTheoryAlternatives> alternatives;
        alternatives.reserve(_encoding.atom_count());

        for(SizeType i=0; i!=_encoding.atom_count(); ++i) {
            SizeType variable=_encoding.atom_variable(i);
            ARIADNE_ASSERT(_assignment[variable]>=0);

            SmtTheoryLiteral literal=make_smt_theory_literal(_encoding.atom(i));
            if(_assignment[variable]==0) {
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
            SmtResult result=_parallel
                ? _solver.solve_parallel(_space,_domain,literals)
                : _solver.solve(_space,_domain,literals);
            add_statistics(_statistics,result.statistics());
            if(result.is_epsilon_sat()) {
                return result.witness();
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
                literals.resize(old_size);
                return witness;
            }
            literals.resize(old_size);
        }
        return std::nullopt;
    }

    SmtSolver const& _solver;
    RealSpace const& _space;
    ExactBoxType const& _domain;
    SmtBooleanEncoding const& _encoding;
    Bool _parallel;
    std::vector<int8_t> _assignment;
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
