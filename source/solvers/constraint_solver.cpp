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

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../utility/tuple.hpp"
#include "../utility/tribool.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/algebra.hpp"
#include "../geometry/box.hpp"
#include "../geometry/grid_paving.hpp"
#include "../function/polynomial.hpp"
#include "../function/function.hpp"
#include "../function/formula.hpp"
#include "../function/procedure.hpp"
#include "../function/constraint.hpp"
#include "../solvers/nonlinear_programming.hpp"
#include "../function/function_mixin.hpp"
#include "../function/taylor_function.hpp"

#include "../solvers/constraint_solver.hpp"
#include "../solvers/solver.hpp"

namespace Ariadne {

typedef Vector<FloatDPApproximation> FloatApproximationVector;
typedef Vector<FloatDPValue> ExactFloatVector;

Bool has_nan(const ExactBoxType& domain);

inline Sweeper<FloatDP> default_sweeper() { return Sweeper<FloatDP>(); }

inline Sign sign(const FloatDP& x) {
    if(x>0) { return Sign::NEGATIVE; }
    else if(x<0) {  return Sign::POSITIVE; }
    else { return Sign::ZERO; }
}

inline Sign sign(const ExactIntervalType& ivl) {
    if(ivl.lower()>0) { return Sign::NEGATIVE; }
    else if(ivl.upper()<0) {  return Sign::POSITIVE; }
    else { return Sign::ZERO; }
}


inline OutputStream& operator<<(OutputStream& os, const EffectiveConstraint& c) {
    if(c.bounds().lower()==c.bounds().upper()) { return os << c.function() << "==" << c.bounds().upper(); }
    if(c.bounds().upper()==infty) { return os << c.bounds().lower() << "<=" << c.function(); }
    if(c.bounds().lower()==-infty) { return os << c.function() << "<=" << c.bounds().upper(); }
    return os << c.bounds().lower() << "<=" << c.function() << "<=" << c.bounds().upper();
}


auto ConstraintSolver::feasible(const ExactBoxType& domain,
                                const List<ValidatedConstraint>& constraints) const
    -> Pair<ValidatedKleenean,ExactPointType>
{
    if(constraints.empty()) { return make_pair(!domain.is_empty(),domain.midpoint()); }

    ValidatedVectorMultivariateFunction function(constraints.size(),constraints[0].function().domain());
    ExactBoxType bounds(constraints.size());

    for(Nat i=0; i!=constraints.size(); ++i) {
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

    static const FloatDPValue XSIGMA=0.125_exact;
    static const FloatDPValue TERR=-1.0_exact*pow(two,-10);
    static const FloatDP _inf = Ariadne::inf;

    ARIADNE_LOG(4,"domain="<<domain<<"\nfunction="<<function<<"\ncodomain="<<codomain<<"\n");
    ARIADNE_ASSERT(codomain.dimension()>0);

    // Make codomain bounded
    UpperBoxType bounds=codomain;
    UpperBoxType image=apply(function,domain);
    ARIADNE_LOG(4,"image="<<image<<"\n");
    for(Nat i=0; i!=image.size(); ++i) {
        if(definitely(disjoint(image[i],codomain[i]))) {
            ARIADNE_LOG(4,"  Proved disjointness using direct evaluation\n");
            return make_pair(false,ExactPointType());
        } else {
            bounds[i]=intersection(codomain[i],image[i]);
        }
    }


    const Nat m=domain.size(); // The total number of variables
    const Nat n=codomain.size(); // The total number of nontrivial constraints
    const Nat l=(m+n)*2; // The total number of lagrange multipliers

    DoublePrecision prec;
    FloatApproximationVector point(m); // The point in the domain which is the current test point
    FloatDPApproximation violation; // An upper bound on amount by which the constraints are violated by the test point
    FloatApproximationVector multipliers(l); // The lagrange multipliers for the constraints
    FloatApproximationVector slack(l); // The slack between the test point and the violated constraints

    FloatDPApproximation& t=violation; FloatApproximationVector& x=multipliers; FloatApproximationVector& y=point; FloatApproximationVector& z=slack; // Aliases for the main quantities used
    const ExactBoxType& d=domain; const ValidatedVectorMultivariateFunction& fn=function; const ExactBoxType& c=codomain; // Aliases for the main quantities used
    ValidatedVectorMultivariateTaylorFunctionModelDP tfn(d,fn,default_sweeper());

    point=static_cast<FloatApproximationVector>(midpoint(d));
    for(Nat k=0; k!=l; ++k) { multipliers[k]=1.0/l; }

    NonlinearInteriorPointOptimiser optimiser;
    optimiser.compute_tz(domain,function,cast_exact_box(bounds),point,violation,slack);

    ARIADNE_LOG(4,"d="<<d<<", f="<<fn<<", c="<<c<<"\n");


    // TODO: Don't use fixed number of steps
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"    t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
        optimiser.feasibility_step(d,fn,c,x,y,z,t);
        if(decide(t>=TERR)) {
            ARIADNE_LOG(4,"t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
            if(definitely(this->check_feasibility(domain,function,codomain,cast_exact(point)))) { return make_pair(true,cast_exact(point)); }
            else { ARIADNE_LOG(2,"f(y)="<<fn(cast_exact(y))<<"\n"); return make_pair(indeterminate,cast_exact(point)); }
        }
    }
    ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");

    if(decide(t<TERR)) {
        // Probably disjoint, so try to prove this
        UpperBoxType subdomain=domain;

        Vector<FloatDPValue> x_exact=cast_exact(x);
        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ValidatedScalarMultivariateTaylorFunctionModelDP txg=ValidatedScalarMultivariateTaylorFunctionModelDP::zero(d,default_sweeper());
        ValidatedNumericType cnst(0,prec);
        for(Nat j=0; j!=n; ++j) {
            txg = txg - (x_exact[j]-x_exact[n+j])*tfn[j];
            cnst += (c[j].upper()*x_exact[j]-c[j].lower()*x_exact[n+j]);
        }
        for(Nat i=0; i!=m; ++i) {
            txg = txg - (x_exact[2*n+i]-x_exact[2*n+m+i])*ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(d,i,default_sweeper());
            cnst += (d[i].upper()*x_exact[2*n+i]-d[i].lower()*x_exact[2*n+m+i]);
        }
        txg = cnst + txg;

        ARIADNE_LOG(4,"    txg="<<txg<<"\n");

        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        this->hull_reduce(subdomain,txg,ExactIntervalType(0,_inf));
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        if(definitely(subdomain.is_empty())) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return make_pair(false,ExactPointType());
        }

        for(Nat i=0; i!=m; ++i) {
            this->box_reduce(subdomain,txg,ExactIntervalType(0,_inf),i);
            ARIADNE_LOG(8,"  dom="<<subdomain<<"\n");
            if(definitely(subdomain.is_empty())) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return make_pair(false,ExactPointType()); }
        }
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");

        //Pair<ExactBoxType,ExactBoxType> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<ExactBoxType,ExactBoxType> sd=d.split();
        Vector<FloatDPApproximation> nx = FloatDPApproximation(1.0_approx-XSIGMA)*x + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size());
        Vector<FloatDPApproximation> ny = midpoint(sd.first);
        ValidatedKleenean result=this->feasible(sd.first, fn, c).first;
        nx = FloatDPApproximation(1-XSIGMA)*x + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        result = result || this->feasible(sd.second, fn, c).first;
        return make_pair(result,ExactPointType());
    }

    return make_pair(indeterminate,ExactPointType());
}


Bool ConstraintSolver::reduce(UpperBoxType& domain, const ValidatedVectorMultivariateFunction& function, const ExactBoxType& codomain) const
{
    const FloatDP MINIMUM_REDUCTION = 0.75;
    ARIADNE_ASSERT(function.argument_size()==domain.size());
    ARIADNE_ASSERT(function.result_size()==codomain.size());

    if(definitely(domain.is_empty())) { return true; }

    FloatDPUpperBound domain_magnitude={0u,dp};
    for(Nat j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width();
    }
    FloatDPUpperBound old_domain_magnitude=domain_magnitude;

    do {
        this->hull_reduce(domain,function,codomain);
        if(definitely(domain.is_empty())) { return true; }

        for(Nat i=0; i!=codomain.size(); ++i) {
            for(Nat j=0; j!=domain.size(); ++j) {
                this->box_reduce(domain,function[i],codomain[i],j);
                if(definitely(domain.is_empty())) { return true; }
            }
        }
        if(definitely(domain.is_empty())) { return true; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0u;
        for(Nat j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude.raw() < old_domain_magnitude.raw() * MINIMUM_REDUCTION);

    return false;
}

Bool has_nan(const ExactBoxType& domain) {
    for(Nat i=0; i!=domain.size(); ++i) {
        if(is_nan(domain[i].lower().raw()) || is_nan(domain[i].upper().raw())) { return true; }
    }
    return false;
}

Bool ConstraintSolver::reduce(UpperBoxType& domain, const List<ValidatedConstraint>& constraints) const
{
    static const Bool USE_BOX_REDUCE = false;

    const double MINIMUM_REDUCTION = 0.75;

    if(definitely(domain.is_empty())) { return true; }

    FloatDPUpperBound domain_magnitude={0u,dp};
    for(Nat j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width();
    }
    FloatDPUpperBound old_domain_magnitude=domain_magnitude;

    do {
        for(Nat i=0; i!=constraints.size(); ++i) {
            this->hull_reduce(domain,constraints[i].function(),constraints[i].bounds());
        }
        if(definitely(domain.is_empty())) { return true; }

        if(USE_BOX_REDUCE) {
            for(Nat i=0; i!=constraints.size(); ++i) {
                for(Nat j=0; j!=domain.size(); ++j) {
                    this->box_reduce(domain,constraints[i].function(),constraints[i].bounds(),j);
                    if(definitely(domain[j].is_empty())) { return true; }
                }
            }
        }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0u;
        for(Nat j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width();
        }
    } while(domain_magnitude.raw() < old_domain_magnitude.raw() * MINIMUM_REDUCTION);

    return false;
}


Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedProcedure& procedure, const ExactIntervalType& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBoxType domain, ValidatedProcedure procedure, ExactIntervalType bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const Vector<ValidatedProcedure>& procedure, const ExactBoxType& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBoxType domain, Vector<ValidatedProcedure> procedure, ExactBoxType bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBoxType domain, ValidatedScalarMultivariateFunction function, ExactIntervalType bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Procedure<ValidatedNumber> procedure(function);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::hull_reduce(UpperBoxType& domain, const ValidatedVectorMultivariateFunction& function, const ExactBoxType& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBoxType domain, ValidatedScalarMultivariateFunction function, ExactIntervalType bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Vector< Procedure<ValidatedNumber> > procedure(function);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::monotone_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds, Nat variable) const
{
    ValidatedScalarMultivariateFunction derivative=function.derivative(variable);

    ARIADNE_LOG(2,"ConstraintSolver::monotone_reduce(ExactBoxType domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<", derivative="<<derivative<<"\n");

    FloatDPValue splitpoint;
    UpperIntervalType lower=domain[variable];
    UpperIntervalType upper=domain[variable];
    Box<UpperIntervalType> slice=domain;
    Box<UpperIntervalType> subdomain=domain;

    static const Int MAX_STEPS=3;
    const FloatDP threshold = lower.width().raw() / (1<<MAX_STEPS);
    do {
        FloatDPUpperBound ub; FloatDPUpperInterval ivl; FloatDPValue val; ub=val;

        // Apply Newton contractor on lower and upper strips
        if(lower.width().raw()>threshold) {
            splitpoint=lower.midpoint();
            slice[variable]=splitpoint;
            UpperIntervalType new_lower=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            if(definitely(new_lower.upper()<lower.lower())) { lower=UpperIntervalType(lower.lower().raw(),lower.lower().raw()); }
            else { lower=intersection(lower,new_lower); }
        }
        if(upper.width().raw()>threshold) {
            splitpoint=upper.midpoint();
            slice[variable]=splitpoint;
            UpperIntervalType new_upper=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            if(definitely(new_upper.lower()>upper.upper())) { upper=UpperIntervalType(upper.upper().raw(),upper.upper().raw()); }
            else { upper=intersection(upper,new_upper); }
        }
        subdomain[variable]=UpperIntervalType(lower.lower(),upper.upper());
    } while(lower.width().raw()>threshold && upper.width().raw()>threshold);
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
    UpperIntervalType C(0);
    for(Nat i=0; i!=function.result_size(); ++i) {
        g += cast_exact(multipliers[i]) * function[i];
        C += cast_exact(multipliers[i]) * bounds[i];
    }
    Covector<UpperIntervalType> dg = gradient_range(g,cast_vector(domain));
    C -= g(centre);

    UpperBoxType new_domain(domain);
    UpperBoxType ranges(domain.size());
    for(Nat j=0; j!=domain.dimension(); ++j) {
        ranges[j] = dg[j]*(domain[j]-centre[j]);
    }

    // We now have sum dg(xi)[j] * (x[j]-x0[j]) in C, so we can reduce each component
    for(Nat j=0; j!=domain.size(); ++j) {
        UpperIntervalType E = C;
        for(Nat k=0; k!=domain.size(); ++k) {
            if(j!=k) { E-=ranges[k]; }
        }
        UpperIntervalType estimated_domain = E/dg[j]+centre[j];
        new_domain[j] = intersection(domain[j],estimated_domain);
    }

    domain=new_domain;
    return definitely(domain.is_empty());
}

Bool ConstraintSolver::box_reduce(UpperBoxType& domain, const ValidatedScalarMultivariateFunction& function, const ExactIntervalType& bounds, Nat variable) const
{
    ARIADNE_LOG(2,"ConstraintSolver::box_reduce(ExactBoxType domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<"\n");

    if(definitely(domain[variable].lower() >= domain[variable].upper())) { return false; }

    // Try to reduce the size of the set by "shaving" off along a coordinate axis
    //
    UpperIntervalType interval=domain[variable];
    RawFloatDP l=interval.lower().raw();
    RawFloatDP u=interval.upper().raw();
    ExactIntervalType subinterval;
    UpperIntervalType new_interval(interval);
    Box<UpperIntervalType> slice=domain;

    static const Nat MAX_SLICES=(1<<3);
    const Nat n=MAX_SLICES;

    // Look for empty slices from below
    Nat imax = n;
    for(Nat i=0; i!=n; ++i) {
        subinterval=ExactIntervalType((l*(n-i)+u*i)/n,(l*(n-i-1)+u*(i+1))/n);
        slice[variable]=subinterval;
        UpperIntervalType slice_image=apply(function,slice);
        if(definitely(intersection(slice_image,bounds).is_empty())) {
            new_interval.set_lower(subinterval.upper());
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
    for(Nat j=n-1; j!=imax; --j) {
        subinterval=ExactIntervalType((l*(n-j)+u*j)/n,(l*(n-j-1)+u*(j+1))/n);
        slice[variable]=subinterval;
        UpperIntervalType slice_image=apply(function,slice);
        if(definitely(intersection(slice_image,bounds).is_empty())) {
            new_interval.set_upper(subinterval.lower());
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
    for(Nat i=0; i!=y.size(); ++i) {
        if(y[i]<d[i].lower() || y[i]>d[i].upper()) { return false; }
    }

    Vector<FloatDPBounds> fy=f(Vector<FloatDPBounds>(y));
    ARIADNE_LOG(4,"d="<<d<<" f="<<f<<", c="<<c<<"\n  y="<<y<<", f(y)="<<fy<<"\n");
    ValidatedKleenean result=true;
    for(Nat j=0; j!=fy.size(); ++j) {
        if(fy[j].lower().raw()>c[j].upper().raw() || fy[j].upper().raw()<c[j].lower().raw()) { return false; }
        if(fy[j].upper().raw()>=c[j].upper().raw() || fy[j].lower().raw()<=c[j].lower().raw()) { result=indeterminate; }
    }
    return result;
}






} // namespace Ariadne
