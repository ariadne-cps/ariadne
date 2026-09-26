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

#include "utility/macros.hpp"
#include "utility/tuple.hpp"
#include "utility/tribool.hpp"
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
    CONCLOG_SCOPE_CREATE;

    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("function="<<function);
    CONCLOG_PRINTLN("codomain="<<codomain);
    ARIADNE_ASSERT(codomain.dimension()>0);

    UpperBoxType image=apply(function,domain);
    CONCLOG_PRINTLN_AT(1,"image="<<image);
    for(SizeType i=0; i!=image.size(); ++i) {
        if(definitely(disjoint(image[i],codomain[i]))) {
            CONCLOG_PRINTLN("Proved disjointness using direct evaluation");
            return make_pair(false,ExactPointType(0u,dp));
        }
    }

    NonlinearInfeasibleInteriorPointOptimiser optimiser;
    auto candidate_result=optimiser.feasible_candidate(
        domain,function,codomain);

    if(definitely(candidate_result.first)) {
        ExactPointType candidate=cast_exact(candidate_result.second);
        if(definitely(this->check_feasibility(
                domain,function,codomain,candidate))) {
            return make_pair(true,candidate);
        }
        ARIADNE_WARN(
            "NonlinearInfeasibleInteriorPointOptimiser reported a validated "
            "feasible candidate that ConstraintSolver could not certify.");
        return make_pair(indeterminate,ExactPointType());
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
        if(definitely(domain.is_empty())) { return true; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0u;
        for(SizeType j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude.raw() < mul(near, old_domain_magnitude.raw(), FloatDP(MINIMUM_REDUCTION,dp)));

    return false;
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
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("procedure="<<procedure);
    CONCLOG_PRINTLN("bounds="<<bounds);

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const Vector<ValidatedProcedure>& procedure, const ExactBoxType& bounds) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("procedure="<<procedure);
    CONCLOG_PRINTLN("bounds="<<bounds);

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("function="<<function);
    CONCLOG_PRINTLN("bounds="<<bounds);

    Procedure<ValidatedNumber> procedure(function);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedVectorMultivariateFunction& function, const ExactBoxType& bounds) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("function="<<function);
    CONCLOG_PRINTLN("bounds="<<bounds);

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
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("function="<<function);
    CONCLOG_PRINTLN("bounds="<<bounds);

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
            if(definitely(new_lower.upper_bound()<lower.lower_bound())) { lower=UpperIntervalType(lower.lower_bound().raw(),lower.lower_bound().raw()); }
            else { lower=intersection(lower,new_lower); }
        }
        if(upper.width().raw()>threshold) {
            splitpoint=upper.midpoint();
            slice[variable]=splitpoint;
            UpperIntervalType new_upper=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            if(definitely(new_upper.lower_bound()>upper.upper_bound())) { upper=UpperIntervalType(upper.upper_bound().raw(),upper.upper_bound().raw()); }
            else { upper=intersection(upper,new_upper); }
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
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("domain="<<domain);
    CONCLOG_PRINTLN("function="<<function);
    CONCLOG_PRINTLN("bounds="<<bounds);

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


Pair<UpperBoxType,UpperBoxType> ConstraintSolver::split(const UpperBoxType& d, const ValidatedVectorMultivariateFunction& f, const ExactBoxType& c) const
{
    return d.split();
}


ValidatedKleenean ConstraintSolver::check_feasibility(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& f, const ExactBoxType& c, const ExactPointType& y) const
{
    CONCLOG_SCOPE_CREATE;

    for(SizeType i=0; i!=y.size(); ++i) {
        if(y[i]<d[i].lower_bound() || y[i]>d[i].upper_bound()) { return false; }
    }

    Vector<FloatDPBounds> fy=f(Vector<FloatDPBounds>(y));
    CONCLOG_PRINTLN("d="<<d<<" f="<<f<<", c="<<c);
    CONCLOG_PRINTLN("y="<<y<<", f(y)="<<fy);
    ValidatedKleenean result=true;
    for(SizeType j=0; j!=fy.size(); ++j) {
        if(fy[j].lower().raw()>c[j].upper_bound().raw() || fy[j].upper().raw()<c[j].lower_bound().raw()) { return false; }
        if(fy[j].upper().raw()>=c[j].upper_bound().raw() || fy[j].lower().raw()<=c[j].lower_bound().raw()) { result=indeterminate; }
    }
    return result;
}






} // namespace Ariadne
