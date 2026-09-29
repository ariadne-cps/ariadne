/***************************************************************************
 *            solvers/constraint_solver.cpp
 *
 *  Copyright  2000-20  Pieter Collins
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

#include "function/functional.hpp"
#include "config.hpp"

#include <chrono>
#include <deque>

#include "utility/macros.hpp"
#include "utility/tuple.hpp"
#include "paradigm/logical.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "geometry/box.hpp"
#include "geometry/grid_paving.hpp"
#include "function/polynomial.hpp"
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/procedure.hpp"
#include "function/constraint.hpp"
#include "solvers/nonlinear_programming.hpp"
#include "function/function_mixin.hpp"
#include "function/taylor_function.hpp"

#include "solvers/constraint_solver.hpp"
#include "solvers/solver.hpp"

namespace Ariadne {

namespace {

inline double constraint_elapsed_seconds(
    std::chrono::steady_clock::time_point const& start)
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
}

} // namespace

template<class X> inline Approximation<X> affine(Approximation<X> l, Approximation<X> u, Nat i, Nat n) {
    return (l*(n-i)+u*i)/n;
}

typedef Vector<FloatDPApproximation> FloatApproximationVector;
typedef Vector<FloatDP> ExactFloatVector;

inline Sweeper<FloatDP> default_sweeper() { return Sweeper<FloatDP>(); }


auto ConstraintSolver::feasible(const ExactBoxType& domain,
                                const List<ValidatedConstraint>& constraints) const
    -> Pair<ValidatedKleenean,ExactPointType>
{
    if(constraints.empty()) { return make_pair(!domain.is_empty(),domain.midpoint()); }

    ValidatedVectorMultivariateFunction function(constraints.size(),constraints[0].function().domain());
    ExactBoxType bounds(constraints.size());

    for(SizeType i=0; i!=constraints.size(); ++i) {
        function[i]=constraints[i].function();
        bounds[i]=constraints[i].bounds();
    }
    return this->feasible(domain,function,bounds);
}


auto ConstraintSolver::feasible(const ExactBoxType& domain,
                                const ValidatedVectorMultivariateFunction& function,
                                const ExactBoxType& codomain) const
    -> Pair<ValidatedKleenean,ExactPointType>
{
    LOGGING_SCOPE_CREATE;

    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("function="<<function);
    LOGGING_PRINTLN("codomain="<<codomain);
    ARIADNE_ASSERT(codomain.dimension()>0);

    UpperBoxType image=apply(function,domain);
    LOGGING_PRINTLN_AT(1,"image="<<image);
    for(SizeType i=0; i!=image.size(); ++i) {
        if(definitely(disjoint(image[i],codomain[i]))) {
            LOGGING_PRINTLN("Proved disjointness using direct evaluation");
            return make_pair(false,ExactPointType(0u,dp));
        }
    }

    NonlinearInfeasibleInteriorPointOptimiser optimiser;
    Pair<ValidatedKleenean,FloatDPApproximationVector> candidate_result;
    try {
        candidate_result=optimiser.feasible_candidate(
            domain,function,codomain);
    } catch(const SingularMatrixException&) {
        CONCLOG_PRINTLN(
            "Interior-point candidate search encountered a singular system");
        return make_pair(indeterminate,ExactPointType());
    }

    if(definitely(candidate_result.first)) {
        return make_pair(true,cast_exact(candidate_result.second));
    }

    if(not possibly(candidate_result.first)) {
        return make_pair(false,ExactPointType());
    }

    return make_pair(indeterminate,ExactPointType());
}

Bool ConstraintSolver::reduce(UpperBoxType& domain, const ValidatedVectorMultivariateFunction& function, const ExactBoxType& codomain) const
{
    const ExactDouble MINIMUM_REDUCTION = 0.75_x;
    ARIADNE_ASSERT(function.argument_size()==domain.size());
    ARIADNE_ASSERT(function.result_size()==codomain.size());

    if(definitely(domain.is_empty())) { return true; }

    FloatDPUpperBound domain_magnitude={0u,dp};
    for(SizeType j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width();
    }
    FloatDPUpperBound old_domain_magnitude=domain_magnitude;

    do {
        this->hull_reduce(domain,function,codomain);
        if(definitely(domain.is_empty())) { return true; }

        for(SizeType i=0; i!=codomain.size(); ++i) {
            for(SizeType j=0; j!=domain.size(); ++j) {
                this->box_reduce(domain,function[i],codomain[i],j);
                if(definitely(domain.is_empty())) { return true; }
            }
        }
        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0u;
        for(SizeType j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude.raw() < mul(near, old_domain_magnitude.raw(), FloatDP(MINIMUM_REDUCTION,dp)));

    return false;
}

Bool ConstraintSolver::propagate(
    UpperBoxType& domain,
    const List<ValidatedConstraint>& constraints,
    ConstraintPropagationStatistics& statistics,
    Bool shaving_reduction_enabled,
    Bool hull_reduction_enabled) const
{
    for(;;) {
        UpperBoxType previous=domain;
        ++statistics.hull_rounds;
        for(SizeType i=0u; i!=constraints.size(); ++i) {
            ExactIntervalType const& bounds=constraints[i].bounds();
            auto phase_start=std::chrono::steady_clock::now();
            if(hull_reduction_enabled) {
                ValidatedProcedure procedure(constraints[i].function());
                ++statistics.hull_procedure_builds;
                statistics.hull_procedure_build_seconds+=
                    constraint_elapsed_seconds(phase_start);

                ProcedureHullReductionStatistics hull_statistics;
                phase_start=std::chrono::steady_clock::now();
                ++statistics.hull_contractor_calls;
                Bool const hull_empty=this->hull_reduce(
                    domain,procedure,bounds,hull_statistics);
                statistics.hull_contraction_seconds+=
                    constraint_elapsed_seconds(phase_start);
                statistics.hull_temporary_allocation_seconds+=
                    hull_statistics.temporary_allocation_seconds;
                statistics.hull_forward_execution_seconds+=
                    hull_statistics.forward_execution_seconds;
                statistics.hull_backward_propagation_seconds+=
                    hull_statistics.backward_propagation_seconds;
                if(hull_empty) {
                    return true;
                }
            }

            phase_start=std::chrono::steady_clock::now();
            UpperIntervalType image=apply(constraints[i].function(),domain);
            statistics.hull_direct_rejection_seconds+=
                constraint_elapsed_seconds(phase_start);
            if(definitely(disjoint(image,bounds))) {
                return true;
            }
        }
        if(not same(domain,previous)) {
            ++statistics.hull_effective;
            shaving_full_refresh=true;
        }

        if(same(domain,previous)) {
            if(not shaving_reduction_enabled) {
                return false;
            }
            UpperBoxType before_shaving=domain;
            ++statistics.shaving_rounds;
            for(SizeType i=0u; i!=constraints.size(); ++i) {
                for(SizeType variable=0u; variable!=domain.dimension(); ++variable) {
                    UpperIntervalType const before_coordinate=domain[variable];
                    auto shaving_start=std::chrono::steady_clock::now();
                    ++statistics.shaving_coordinate_attempts;
                    if(this->box_reduce(
                            domain,
                            constraints[i].function(),
                            constraints[i].bounds(),
                            variable,
                            statistics.shaving_function_evaluations)) {
                        statistics.shaving_seconds+=
                            constraint_elapsed_seconds(shaving_start);
                        return true;
                    }
                    statistics.shaving_seconds+=
                        constraint_elapsed_seconds(shaving_start);
                    Bool const same_lower=
                        before_coordinate.lower_bound().raw()
                            ==domain[variable].lower_bound().raw();
                    Bool const same_upper=
                        before_coordinate.upper_bound().raw()
                            ==domain[variable].upper_bound().raw();
                    if(not (same_lower && same_upper)) {
                        ++statistics.shaving_coordinate_effective;
                    }
                }
            }
            if(not same(domain,before_shaving)) {
                ++statistics.shaving_effective;
            }
            if(same(domain,before_shaving)) {
                return false;
            }
            continue;
        }

        for(SizeType i=0u; i!=constraints.size(); ++i) {
            UpperIntervalType image=apply(constraints[i].function(),domain);
            if(definitely(disjoint(image,constraints[i].bounds()))) {
                return true;
            }
        }
    }
}



namespace {

Bool propagation_constraint_infeasible(
    ConstraintPropagationConstraint const& constraint,
    UpperBoxType const& domain)
{
    UpperIntervalType image=apply(constraint.function,domain);
    if(constraint.strict_lower
       && definitely(image.upper_bound()<=constraint.bounds.lower_bound())) {
        return true;
    }
    if(constraint.strict_upper
       && definitely(image.lower_bound()>=constraint.bounds.upper_bound())) {
        return true;
    }
    return definitely(disjoint(image,constraint.bounds));
}

Bool propagation_monotone_coordinate_is_safe(
    ValidatedScalarMultivariateFunction const& derivative,
    UpperBoxType const& domain)
{
    UpperIntervalType derivative_image=apply(derivative,domain);
    return definitely(derivative_image.lower_bound()>0)
        || definitely(derivative_image.upper_bound()<0);
}

std::vector<SizeType> propagation_procedure_dependencies(
    ValidatedProcedure const& procedure)
{
    std::vector<SizeType> dependencies;
    std::vector<Bool> seen(procedure.argument_size(),false);
    for(SizeType i=0u;i!=procedure._instructions.size();++i) {
        auto const* variable=
            std::get_if<IndexProcedureInstruction>(
                &procedure._instructions[i].base());
        if(variable==nullptr) {
            continue;
        }
        SizeType const index=variable->_ind;
        if(index<seen.size() && not seen[index]) {
            seen[index]=true;
            dependencies.push_back(index);
        }
    }
    return dependencies;
}

std::vector<SizeType> propagation_changed_variables(
    UpperBoxType const& before,
    UpperBoxType const& after)
{
    std::vector<SizeType> changed;
    for(SizeType variable=0u;variable!=before.dimension();++variable) {
        Bool const same_lower=
            before[variable].lower_bound().raw()
                ==after[variable].lower_bound().raw();
        Bool const same_upper=
            before[variable].upper_bound().raw()
                ==after[variable].upper_bound().raw();
        if(not (same_lower && same_upper)) {
            changed.push_back(variable);
        }
    }
    return changed;
}

Bool propagation_has_cached_procedures(
    std::vector<ConstraintPropagationConstraint> const& constraints)
{
    for(auto const& constraint:constraints) {
        if(not constraint.hull_procedure) {
            return false;
        }
    }
    return true;
}

} // namespace

Bool ConstraintSolver::propagate(
    UpperBoxType& domain,
    const std::vector<ConstraintPropagationConstraint>& constraints,
    Bool monotone_reduction_enabled,
    ConstraintPropagationStatistics& statistics,
    Bool shaving_reduction_enabled,
    Bool hull_reduction_enabled) const
{
    // First incremental vertical slice: use an IBEX-style propagation agenda
    // only for the pure cached hull-contraction phase. Shaving and monotone
    // contraction retain the established full-scan implementation below until
    // their wake-up semantics are measured separately.
    if(hull_reduction_enabled
       && not shaving_reduction_enabled
       && not monotone_reduction_enabled
       && propagation_has_cached_procedures(constraints)) {
        ++statistics.hull_rounds;

        std::vector<std::vector<SizeType>> dependencies(constraints.size());
        std::vector<std::vector<SizeType>> watchers(domain.dimension());
        for(SizeType constraint_index=0u;
            constraint_index!=constraints.size();
            ++constraint_index) {
            dependencies[constraint_index]=
                propagation_procedure_dependencies(
                    *constraints[constraint_index].hull_procedure);
            for(SizeType variable:dependencies[constraint_index]) {
                if(variable<watchers.size()) {
                    watchers[variable].push_back(constraint_index);
                }
            }
        }

        std::deque<SizeType> agenda;
        std::vector<Bool> queued(constraints.size(),false);
        auto enqueue=[&](SizeType constraint_index) {
            if(not queued[constraint_index]) {
                queued[constraint_index]=true;
                agenda.push_back(constraint_index);
                ++statistics.hull_agenda_pushes;
            }
        };
        for(SizeType i=0u;i!=constraints.size();++i) {
            enqueue(i);
        }

        Bool any_effective=false;
        while(not agenda.empty()) {
            SizeType const constraint_index=agenda.front();
            agenda.pop_front();
            queued[constraint_index]=false;
            ++statistics.hull_agenda_pops;

            auto const& constraint=constraints[constraint_index];
            UpperBoxType before=domain;

            ProcedureHullReductionStatistics hull_statistics;
            auto phase_start=std::chrono::steady_clock::now();
            ++statistics.hull_contractor_calls;
            Bool const hull_empty=this->hull_reduce(
                domain,*constraint.hull_procedure,
                constraint.bounds,hull_statistics);
            statistics.hull_contraction_seconds+=
                constraint_elapsed_seconds(phase_start);
            statistics.hull_temporary_allocation_seconds+=
                hull_statistics.temporary_allocation_seconds;
            statistics.hull_forward_execution_seconds+=
                hull_statistics.forward_execution_seconds;
            statistics.hull_backward_propagation_seconds+=
                hull_statistics.backward_propagation_seconds;
            if(hull_empty) {
                return true;
            }

            phase_start=std::chrono::steady_clock::now();
            Bool const infeasible=
                propagation_constraint_infeasible(constraint,domain);
            statistics.hull_direct_rejection_seconds+=
                constraint_elapsed_seconds(phase_start);
            if(infeasible) {
                return true;
            }

            auto changed=propagation_changed_variables(before,domain);
            if(changed.empty()) {
                continue;
            }
            any_effective=true;
            ++statistics.hull_agenda_effective_calls;
            for(SizeType variable:changed) {
                for(SizeType dependent:watchers[variable]) {
                    enqueue(dependent);
                }
            }
        }

        if(any_effective) {
            ++statistics.hull_effective;
        }
        return false;
    }

    std::vector<std::vector<SizeType>> constraint_dependencies;
    constraint_dependencies.reserve(constraints.size());
    for(auto const& constraint:constraints) {
        if(constraint.hull_procedure) {
            constraint_dependencies.push_back(
                propagation_procedure_dependencies(*constraint.hull_procedure));
        } else {
            std::vector<SizeType> all_variables;
            all_variables.reserve(domain.dimension());
            for(SizeType variable=0u;variable!=domain.dimension();++variable) {
                all_variables.push_back(variable);
            }
            constraint_dependencies.push_back(std::move(all_variables));
        }
    }

    std::vector<std::vector<unsigned char>> shaving_active(
        constraints.size(),
        std::vector<unsigned char>(domain.dimension(),0u));
    for(SizeType constraint_index=0u;
        constraint_index!=constraints.size();
        ++constraint_index) {
        for(SizeType variable:constraint_dependencies[constraint_index]) {
            shaving_active[constraint_index][variable]=1u;
        }
    }
    Bool shaving_full_refresh=true;

    Bool monotone_attempted=false;
    for(;;) {
        UpperBoxType previous=domain;
        ++statistics.hull_rounds;
        for(auto const& constraint:constraints) {
            auto phase_start=std::chrono::steady_clock::now();
            if(hull_reduction_enabled) {
                std::optional<ValidatedProcedure> local_procedure;
                ValidatedProcedure const* procedure=nullptr;
                if(constraint.hull_procedure) {
                    procedure=constraint.hull_procedure.get();
                } else {
                    local_procedure.emplace(constraint.function);
                    procedure=&*local_procedure;
                    ++statistics.hull_procedure_builds;
                    statistics.hull_procedure_build_seconds+=
                        constraint_elapsed_seconds(phase_start);
                }

                ProcedureHullReductionStatistics hull_statistics;
                phase_start=std::chrono::steady_clock::now();
                ++statistics.hull_contractor_calls;
                Bool const hull_empty=this->hull_reduce(
                    domain,*procedure,constraint.bounds,hull_statistics);
                statistics.hull_contraction_seconds+=
                    constraint_elapsed_seconds(phase_start);
                statistics.hull_temporary_allocation_seconds+=
                    hull_statistics.temporary_allocation_seconds;
                statistics.hull_forward_execution_seconds+=
                    hull_statistics.forward_execution_seconds;
                statistics.hull_backward_propagation_seconds+=
                    hull_statistics.backward_propagation_seconds;
                if(hull_empty) {
                    return true;
                }
            }

            phase_start=std::chrono::steady_clock::now();
            Bool const infeasible=propagation_constraint_infeasible(
                constraint,domain);
            statistics.hull_direct_rejection_seconds+=
                constraint_elapsed_seconds(phase_start);
            if(infeasible) {
                return true;
            }
        }
        if(not same(domain,previous)) {
            ++statistics.hull_effective;
        }

        if(same(domain,previous)) {
            UpperBoxType before_shaving=domain;
            Bool shaving_refresh_without_change=false;
            if(shaving_reduction_enabled) {
                ++statistics.shaving_rounds;
                if(shaving_full_refresh) {
                    ++statistics.shaving_refresh_rounds;
                } else {
                    ++statistics.shaving_active_rounds;
                }

                std::vector<std::vector<unsigned char>> next_active(
                    constraints.size(),
                    std::vector<unsigned char>(domain.dimension(),0u));

                for(SizeType constraint_index=0u;
                    constraint_index!=constraints.size();
                    ++constraint_index) {
                    auto const& constraint=constraints[constraint_index];
                    statistics.shaving_dependency_skipped+=
                        domain.dimension()
                        -constraint_dependencies[constraint_index].size();
                    for(SizeType variable:constraint_dependencies[constraint_index]) {
                        if(not shaving_full_refresh
                           && shaving_active[constraint_index][variable]==0u) {
                            ++statistics.shaving_adaptive_skipped;
                            continue;
                        }

                        UpperIntervalType const before_coordinate=domain[variable];
                        auto shaving_start=std::chrono::steady_clock::now();
                        ++statistics.shaving_coordinate_attempts;
                        if(this->box_reduce(
                                domain,constraint.function,constraint.bounds,variable,
                                statistics.shaving_function_evaluations)) {
                            statistics.shaving_seconds+=
                                constraint_elapsed_seconds(shaving_start);
                            return true;
                        }
                        statistics.shaving_seconds+=
                            constraint_elapsed_seconds(shaving_start);
                        Bool const same_lower=
                            before_coordinate.lower_bound().raw()
                                ==domain[variable].lower_bound().raw();
                        Bool const same_upper=
                            before_coordinate.upper_bound().raw()
                                ==domain[variable].upper_bound().raw();
                        if(not (same_lower && same_upper)) {
                            ++statistics.shaving_coordinate_effective;
                            next_active[constraint_index][variable]=1u;
                        }
                    }
                }

                if(not same(domain,before_shaving)) {
                    ++statistics.shaving_effective;
                    shaving_active=std::move(next_active);
                    shaving_full_refresh=false;
                } else if(not shaving_full_refresh) {
                    // The learned working set stalled. Before declaring a
                    // fixed point, reactivate all genuine dependencies once:
                    // a previously ineffective pair may have become useful
                    // after contractions performed by other pairs.
                    shaving_full_refresh=true;
                    continue;
                } else {
                    shaving_refresh_without_change=true;
                }
            }
            if(same(domain,before_shaving)
               && (not shaving_reduction_enabled
                   || shaving_refresh_without_change)) {
                if(not monotone_reduction_enabled) {
                    return false;
                }
                if(monotone_attempted) {
                    return false;
                }
                monotone_attempted=true;
                UpperBoxType before_monotone=domain;
                ++statistics.monotone_rounds;
                for(auto const& constraint:constraints) {
                    for(SizeType variable=0u; variable!=domain.dimension(); ++variable) {
                        if(variable>=constraint.derivatives.size()) {
                            continue;
                        }
                        auto const& derivative=constraint.derivatives[variable];
                        if(not derivative.has_value()) {
                            continue;
                        }
                        if(not propagation_monotone_coordinate_is_safe(
                                *derivative,domain)) {
                            continue;
                        }
                        this->monotone_reduce(
                            domain,
                            constraint.function,
                            *derivative,
                            constraint.bounds,
                            variable);
                    }
                }
                if(not same(domain,before_monotone)) {
                    ++statistics.monotone_effective;
                    continue;
                }
                return false;
            }
            continue;
        }

        for(auto const& constraint:constraints) {
            if(propagation_constraint_infeasible(constraint,domain)) {
                return true;
            }
        }
    }
}

Bool ConstraintSolver::reduce(UpperBoxType& domain, const List<ValidatedConstraint>& constraints) const
{
    const double MINIMUM_REDUCTION = 0.75;

    if(definitely(domain.is_empty())) { return true; }

    FloatDPUpperBound domain_magnitude={0u,dp};
    for(SizeType j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width();
    }
    FloatDPUpperBound old_domain_magnitude=domain_magnitude;

    do {
        for(SizeType i=0; i!=constraints.size(); ++i) {
            this->hull_reduce(domain,constraints[i].function(),constraints[i].bounds());
        }
        if(definitely(domain.is_empty())) { return true; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0u;
        for(SizeType j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude.raw() < (old_domain_magnitude * MINIMUM_REDUCTION).raw());

    return false;
}


Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedProcedure& procedure, const ExactIntervalType& bounds) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("procedure="<<procedure);
    LOGGING_PRINTLN("bounds="<<bounds);

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(
    UpperBoxType& domain,
    const ValidatedProcedure& procedure,
    const ExactIntervalType& bounds,
    ProcedureHullReductionStatistics& statistics) const
{
    Ariadne::simple_hull_reduce(domain,procedure,bounds,statistics);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const Vector<ValidatedProcedure>& procedure, const ExactBoxType& bounds) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("procedure="<<procedure);
    LOGGING_PRINTLN("bounds="<<bounds);

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("function="<<function);
    LOGGING_PRINTLN("bounds="<<bounds);

    Procedure<ValidatedNumber> procedure(function);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedVectorMultivariateFunction& function, const ExactBoxType& bounds) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("function="<<function);
    LOGGING_PRINTLN("bounds="<<bounds);

    Vector< Procedure<ValidatedNumber> > procedure(function);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::monotone_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds, SizeType variable) const
{
    return this->monotone_reduce(
        domain,function,function.derivative(variable),bounds,variable);
}

Bool ConstraintSolver::monotone_reduce(
    UpperBoxType& domain,
    const ValidatedScalarMultivariateFunction& function,
    const ValidatedScalarMultivariateFunction& derivative,
    const ExactIntervalType& bounds,
    SizeType variable) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("function="<<function);
    LOGGING_PRINTLN("bounds="<<bounds);

    FloatDP splitpoint(dp);
    UpperIntervalType lower=domain[variable];
    UpperIntervalType upper=domain[variable];
    Box<UpperIntervalType> slice=domain;
    Box<UpperIntervalType> subdomain=domain;

    static const Int MAX_STEPS=3;
    const FloatDP threshold = div(near, lower.width().raw(), FloatDP(pow(two,MAX_STEPS),dp));
    for(Int step=0; step!=MAX_STEPS; ++step) {
        FloatDPUpperBound ub(dp); FloatDPUpperInterval ivl(-ub,+ub); FloatDP val(dp); ub=val;

        // Apply Newton contractor on lower and upper strips.
        if(lower.width().raw()>threshold) {
            splitpoint=lower.midpoint();
            slice[variable]=splitpoint;
            UpperIntervalType new_lower=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            lower=intersection(lower,new_lower);
        }
        if(upper.width().raw()>threshold) {
            splitpoint=upper.midpoint();
            slice[variable]=splitpoint;
            UpperIntervalType new_upper=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            upper=intersection(upper,new_upper);
        }
        subdomain[variable]=UpperIntervalType(lower.lower_bound(),upper.upper_bound());
        if(not (lower.width().raw()>threshold && upper.width().raw()>threshold)) {
            break;
        }
    }
    domain=subdomain;

    return definitely(domain.is_empty());
}



Bool ConstraintSolver::lyapunov_reduce(UpperBoxType& domain, const ValidatedVectorMultivariateTaylorFunctionModelDP& function, const ExactBoxType& bounds,
                                       FloatApproximationVector centre, FloatApproximationVector multipliers) const
{
    return this->lyapunov_reduce(domain,function,bounds,cast_exact(centre),cast_exact(multipliers));
}


Bool ConstraintSolver::lyapunov_reduce(UpperBoxType& domain, const ValidatedVectorMultivariateTaylorFunctionModelDP& function, const ExactBoxType& bounds,
                                       ExactFloatVector centre, ExactFloatVector multipliers) const
{
    ValidatedScalarMultivariateTaylorFunctionModelDP g(function.domain(),default_sweeper());
    UpperIntervalType C(0,0,dp);
    for(SizeType i=0; i!=function.result_size(); ++i) {
        g += cast_exact(multipliers[i]) * function[i];
        C += cast_exact(multipliers[i]) * bounds[i];
    }
    Covector<UpperIntervalType> dg = gradient_range(g,cast_vector(domain));
    C -= g(centre);

    UpperBoxType new_domain(domain);
    UpperBoxType ranges(domain.size());
    for(SizeType j=0; j!=domain.dimension(); ++j) {
        ranges[j] = dg[j]*(domain[j]-centre[j]);
    }

    // We now have sum dg(xi)[j] * (x[j]-x0[j]) in C, so we can reduce each component
    for(SizeType j=0; j!=domain.size(); ++j) {
        UpperIntervalType E = C;
        for(SizeType k=0; k!=domain.size(); ++k) {
            if(j!=k) { E-=ranges[k]; }
        }
        UpperIntervalType estimated_domain = E/dg[j]+centre[j];
        new_domain[j] = intersection(domain[j],estimated_domain);
    }

    domain=new_domain;
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::box_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds, SizeType variable) const
{
    SizeType function_evaluations=0u;
    return this->box_reduce(domain,function,bounds,variable,function_evaluations);
}

Bool ConstraintSolver::box_reduce(
    UpperBoxType& domain,
    const ValidatedScalarMultivariateFunction& function,
    const ExactIntervalType& bounds,
    SizeType variable,
    SizeType& function_evaluations) const
{
    LOGGING_SCOPE_CREATE;
    LOGGING_PRINTLN("domain="<<domain);
    LOGGING_PRINTLN("function="<<function);
    LOGGING_PRINTLN("bounds="<<bounds);

    if(definitely(domain[variable].lower_bound() >= domain[variable].upper_bound())) { return false; }

    // Try to reduce the size of the set by "shaving" off along a coordinate axis
    //
    UpperIntervalType interval=domain[variable];
    FloatDPApproximation l=interval.lower_bound();
    FloatDPApproximation u=interval.upper_bound();
    ExactIntervalType subinterval;
    UpperIntervalType new_interval(interval);
    Box<UpperIntervalType> slice=domain;

    static const Nat MAX_SLICES=(1<<3);
    const Nat n=MAX_SLICES;

    // Look for empty slices from below
    Nat imax = n;
    for(Nat i=0; i!=n; ++i) {
        subinterval=ExactIntervalType(cast_exact(affine(l,u,i,n)),cast_exact(affine(l,u,i+1,n)));
        slice[variable]=subinterval;
        ++function_evaluations;
        UpperIntervalType slice_image=apply(function,slice);
        if(definitely(intersection(slice_image,bounds).is_empty())) {
            new_interval.set_lower_bound(subinterval.upper_bound());
        } else {
            imax = i; break;
        }
    }

    // The set is proved to be empty
    if(imax==n) {
        domain[variable]=ExactIntervalType(+inf,-inf);
        return true;
    }

    // Look for empty slices from above; note that at least one nonempty slice has been found
    for(SizeType j=n-1; j!=imax; --j) {
        subinterval=ExactIntervalType(cast_exact(affine(l,u,j,n)),cast_exact(affine(l,u,j+1,n)));
        slice[variable]=subinterval;
        ++function_evaluations;
        UpperIntervalType slice_image=apply(function,slice);
        if(definitely(intersection(slice_image,bounds).is_empty())) {
            new_interval.set_upper_bound(subinterval.lower_bound());
        } else {
            break;
        }
    }

    // The set cannot be empty, since a nonempty slice has been found in the upper pass.
    // Note that the interval is an UpperIntervalType, so non-emptiness of the approximated set cannot be guaranteed,
    // but emptiness would be verified
    ARIADNE_ASSERT(not definitely(new_interval.is_empty()));

    domain[variable]=new_interval;

    return false;
}


Pair<UpperBoxType,UpperBoxType> ConstraintSolver::split(const UpperBoxType& d, const ValidatedVectorMultivariateFunction&, const ExactBoxType&) const
{
    return d.split();
}


ValidatedKleenean ConstraintSolver::check_feasibility(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& f, const ExactBoxType& c, const ExactPointType& y) const
{
    LOGGING_SCOPE_CREATE;

    for(SizeType i=0; i!=y.size(); ++i) {
        if(y[i]<d[i].lower_bound() || y[i]>d[i].upper_bound()) { return false; }
    }

    Vector<FloatDPBounds> fy=f(Vector<FloatDPBounds>(y));
    LOGGING_PRINTLN("d="<<d<<" f="<<f<<", c="<<c);
    LOGGING_PRINTLN("y="<<y<<", f(y)="<<fy);
    ValidatedKleenean result=true;
    for(SizeType j=0; j!=fy.size(); ++j) {
        if(fy[j].lower().raw()>c[j].upper_bound().raw() || fy[j].upper().raw()<c[j].lower_bound().raw()) { return false; }
        if(fy[j].upper().raw()>=c[j].upper_bound().raw() || fy[j].lower().raw()<=c[j].lower_bound().raw()) { result=indeterminate; }
    }
    return result;
}






} // namespace Ariadne
