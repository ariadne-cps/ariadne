/***************************************************************************
 *            constraint_solver.cc
 *
 *  Copyright 2000  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "function/functional.h"
#include "config.h"

#include "utility/macros.h"
#include "utility/tuple.h"
#include "utility/tribool.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "geometry/box.h"
#include "geometry/grid_set.h"
#include "function/polynomial.h"
#include "function/function.h"
#include "expression/formula.h"
#include "function/procedure.h"
#include "function/constraint.h"
#include "solvers/nonlinear_programming.h"
#include "function/function_mixin.h"
#include "function/taylor_function.h"

#include "solvers/constraint_solver.h"
#include "solvers/solver.h"

namespace Ariadne {

typedef Vector<ApproximateFloat> ApproximateFloatVector;
typedef Vector<ExactFloat> ExactFloatVector;

inline Sweeper default_sweeper() { return Sweeper(); }

Sign sign(const Float& x) {
    if(x>0) { return NEGATIVE; }
    else if(x<0) {  return POSITIVE; }
    else { return ZERO; }
}

Sign sign(const ExactInterval& ivl) {
    if(ivl.lower()>0) { return NEGATIVE; }
    else if(ivl.upper()<0) {  return POSITIVE; }
    else { return ZERO; }
}


OutputStream& operator<<(OutputStream& os, const EffectiveConstraint& c) {
    static const Float inf = Ariadne::inf;
    if(c.bounds().lower()==c.bounds().upper()) { return os << c.function() << "==" << c.bounds().upper(); }
    if(c.bounds().upper()==infty) { return os << c.bounds().lower() << "<=" << c.function(); }
    if(c.bounds().lower()==-infty) { return os << c.function() << "<=" << c.bounds().upper(); }
    return os << c.bounds().lower() << "<=" << c.function() << "<=" << c.bounds().upper();
}



Pair<Tribool,ExactPoint> ConstraintSolver::feasible(const ExactBox& domain, const List<ValidatedConstraint>& constraints) const
{
    if(constraints.empty()) { return make_pair(!domain.empty(),domain.centre()); }

    ValidatedVectorFunction function(constraints.size(),constraints[0].function().argument_size());
    ExactBox bounds(constraints.size());

    for(Nat i=0; i!=constraints.size(); ++i) {
        function[i]=constraints[i].function();
        bounds[i]=constraints[i].bounds();
    }
    return this->feasible(domain,function,bounds);
}


Pair<Tribool,ExactPoint> ConstraintSolver::feasible(const ExactBox& domain, const ValidatedVectorFunction& function, const ExactBox& codomain) const
{

    static const double XSIGMA=0.125;
    static const double TERR=-1.0/(1<<10);
    static const Float inf = Ariadne::inf;

    ARIADNE_LOG(4,"domain="<<domain<<"\nfunction="<<function<<"\ncodomain="<<codomain<<"\n");

    // Make codomain bounded
    UpperBox bounds=codomain;
    UpperBox image=apply(function,domain);
    ARIADNE_LOG(4,"image="<<image<<"\n");
    for(Nat i=0; i!=image.size(); ++i) {
        if(definitely(disjoint(image[i],codomain[i]))) {
            ARIADNE_LOG(4,"  Proved disjointness using direct evaluation\n");
            return make_pair(false,ExactPoint());
        } else {
            bounds[i]=intersection(codomain[i],image[i]);
        }
    }


    const Nat m=domain.size(); // The total number of variables
    const Nat n=codomain.size(); // The total number of nontrivial constraints
    const Nat l=(m+n)*2; // The total number of lagrange multipliers

    ApproximateFloatVector point(m); // The point in the domain which is the current test point
    ApproximateFloat violation; // An upper bound on amount by which the constraints are violated by the test point
    ApproximateFloatVector multipliers(l); // The lagrange multipliers for the constraints
    ApproximateFloatVector slack(l); // The slack between the test point and the violated constraints

    ApproximateFloat& t=violation; ApproximateFloatVector& x=multipliers; ApproximateFloatVector& y=point; ApproximateFloatVector& z=slack; // Aliases for the main quantities used
    const ExactBox& d=domain; const ValidatedVectorFunction& fn=function; const ExactBox& c=codomain; // Aliases for the main quantities used
    VectorTaylorFunction tfn(d,fn,default_sweeper());

    point=static_cast<ApproximateFloatVector>(midpoint(d));
    for(Nat k=0; k!=l; ++k) { multipliers[k]=1.0/l; }

    NonlinearInteriorPointOptimiser optimiser;
    optimiser.compute_tz(domain,function,make_exact_box(bounds),point,violation,slack);

    ARIADNE_LOG(4,"d="<<d<<", f="<<fn<<", c="<<c<<"\n");


    // TODO: Don't use fixed number of steps
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"    t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
        optimiser.feasibility_step(d,fn,c,x,y,z,t);
        if(t>=TERR) {
            ARIADNE_LOG(4,"t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
            if(definitely(this->check_feasibility(domain,function,codomain,make_exact(point)))) { return make_pair(true,make_exact(point)); }
            else { ARIADNE_LOG(2,"f(y)="<<fn(make_exact(y))<<"\n"); return make_pair(indeterminate,make_exact(point)); }
        }
    }
    ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        ExactBox subdomain=domain;

        Vector<ExactFloat> x_exact=make_exact(x);
        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ScalarTaylorFunction txg=ScalarTaylorFunction::zero(d,default_sweeper());
        ValidatedNumber cnst=0;
        for(Nat j=0; j!=n; ++j) {
            txg = txg - (x_exact[j]-x_exact[n+j])*tfn[j];
            cnst += (c[j].upper()*x_exact[j]-c[j].lower()*x_exact[n+j]);
        }
        for(Nat i=0; i!=m; ++i) {
            txg = txg - (x_exact[2*n+i]-x_exact[2*n+m+i])*ScalarTaylorFunction::coordinate(d,i,default_sweeper());
            cnst += (d[i].upper()*x_exact[2*n+i]-d[i].lower()*x_exact[2*n+m+i]);
        }
        txg = cnst + txg;

        ARIADNE_LOG(4,"    txg="<<txg<<"\n");

        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        this->hull_reduce(subdomain,txg,ExactInterval(0,inf));
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        if(subdomain.empty()) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return make_pair(false,ExactPoint());
        }

        for(Nat i=0; i!=m; ++i) {
            this->box_reduce(subdomain,txg,ExactInterval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<subdomain<<"\n");
            if(subdomain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return make_pair(false,ExactPoint()); }
        }
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");

        //Pair<ExactBox,ExactBox> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<ExactBox,ExactBox> sd=d.split();
        Vector<ApproximateFloat> nx = ApproximateFloat(1.0-XSIGMA)*x + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size());
        Vector<ApproximateFloat> ny = midpoint(sd.first);
        Tribool result=this->feasible(sd.first, fn, c).first;
        nx = ApproximateFloat(1.0-XSIGMA)*x + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        result = result || this->feasible(sd.second, fn, c).first;
        return make_pair(result,ExactPoint());
    }

    return make_pair(indeterminate,ExactPoint());
}


Bool ConstraintSolver::reduce(UpperBox& domain, const ValidatedVectorFunction& function, const ExactBox& codomain) const
{
    const double MINIMUM_REDUCTION = 0.75;
    ARIADNE_ASSERT(function.argument_size()==domain.size());
    ARIADNE_ASSERT(function.result_size()==codomain.size());

    if(definitely(domain.empty())) { return true; }

    Float domain_magnitude=0.0;
    for(Nat j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width().raw();
    }
    Float old_domain_magnitude=domain_magnitude;

    do {
        this->hull_reduce(domain,function,codomain);
        if(definitely(domain.empty())) { return true; }

        for(Nat i=0; i!=codomain.size(); ++i) {
            for(Nat j=0; j!=domain.size(); ++j) {
                this->box_reduce(domain,function[i],codomain[i],j);
                if(definitely(domain.empty())) { return true; }
            }
        }
        if(definitely(domain.empty())) { return true; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0.0;
        for(Nat j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width().raw();
        }
    } while(domain_magnitude/old_domain_magnitude <= MINIMUM_REDUCTION);

    return false;
}

Bool has_nan(const ExactBox& domain) {
    for(Nat i=0; i!=domain.size(); ++i) {
        if(is_nan(domain[i].lower().raw()) || is_nan(domain[i].upper().raw())) { return true; }
    }
    return false;
}

Bool ConstraintSolver::reduce(UpperBox& domain, const List<ValidatedConstraint>& constraints) const
{
    static const Bool USE_BOX_REDUCE = false;

    const double MINIMUM_REDUCTION = 0.75;

    if(definitely(domain.empty())) { return true; }

    Float domain_magnitude=0.0;
    for(Nat j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width().raw();
    }
    Float old_domain_magnitude=domain_magnitude;

    do {
        for(Nat i=0; i!=constraints.size(); ++i) {
            this->hull_reduce(domain,constraints[i].function(),constraints[i].bounds());
        }
        if(definitely(domain.empty())) { return true; }

        if(USE_BOX_REDUCE) {
            for(Nat i=0; i!=constraints.size(); ++i) {
                for(Nat j=0; j!=domain.size(); ++j) {
                    this->box_reduce(domain,constraints[i].function(),constraints[i].bounds(),j);
                    if(definitely(domain[j].empty())) { return true; }
                }
            }
        }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0.0;
        for(Nat j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width().raw();
        }
    } while(domain_magnitude/old_domain_magnitude <= MINIMUM_REDUCTION);

    return false;
}


Bool ConstraintSolver::hull_reduce(UpperBox& domain, const ValidatedProcedure& procedure, const ExactInterval& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBox domain, ValidatedProcedure procedure, ExactInterval bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.empty());
}

Bool ConstraintSolver::hull_reduce(UpperBox& domain, const Vector<ValidatedProcedure>& procedure, const ExactBox& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBox domain, Vector<ValidatedProcedure> procedure, ExactBox bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return definitely(domain.empty());
}

Bool ConstraintSolver::hull_reduce(UpperBox& domain, const ValidatedScalarFunctionInterface& function, const ExactInterval& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBox domain, ValidatedScalarFunction function, ExactInterval bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Formula<ValidatedNumber> formula=function.evaluate(Formula<ValidatedNumber>::identity(function.argument_size()));
    Procedure<ValidatedNumber> procedure(formula);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::hull_reduce(UpperBox& domain, const ValidatedVectorFunctionInterface& function, const ExactBox& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBox domain, ValidatedScalarFunction function, ExactInterval bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Vector< Formula<ValidatedNumber> > formula=function.evaluate(Formula<ValidatedNumber>::identity(function.argument_size()));
    Vector< Procedure<ValidatedNumber> > procedure(formula);
    return this->hull_reduce(domain,procedure,bounds);
}

Bool ConstraintSolver::monotone_reduce(UpperBox& domain, const ValidatedScalarFunctionInterface& function, const ExactInterval& bounds, Nat variable) const
{
    ValidatedScalarFunction derivative=function.derivative(variable);

    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(ExactBox domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<", derivative="<<derivative<<"\n");

    ExactFloat splitpoint;
    UpperInterval lower=domain[variable];
    UpperInterval upper=domain[variable];
    Vector<UpperInterval> slice=domain;
    Vector<UpperInterval> subdomain=domain;

    static const Int MAX_STEPS=3;
    const Float size = lower.width().raw() / (1<<MAX_STEPS);
    do {
        // Apply Newton contractor on lower and upper strips
        if(lower.width().raw()>size) {
            splitpoint=lower.centre();
            slice[variable]=splitpoint;
            UpperInterval new_lower=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            if(new_lower.upper_raw()<lower.lower_raw()) { lower=UpperInterval(lower.lower_raw()); }
            else { lower=intersection(lower,new_lower); }
        }
        if(upper.width().raw()>size) {
            splitpoint=upper.centre();
            slice[variable]=splitpoint;
            UpperInterval new_upper=splitpoint+(bounds-apply(function,slice))/apply(derivative,subdomain);
            if(new_upper.lower_raw()>upper.upper_raw()) { upper=UpperInterval(upper.upper_raw()); }
            else { upper=intersection(upper,new_upper); }
        }
        subdomain[variable]=UpperInterval(lower.lower(),upper.upper());
    } while(lower.width().raw()>size && upper.width().raw()>size);
    domain=subdomain;

    return definitely(domain.empty());
}



Bool ConstraintSolver::lyapunov_reduce(UpperBox& domain, const VectorTaylorFunction& function, const ExactBox& bounds,
                                       ApproximateFloatVector centre, ApproximateFloatVector multipliers) const
{
    return this->lyapunov_reduce(domain,function,bounds,make_exact(centre),make_exact(multipliers));
}


Bool ConstraintSolver::lyapunov_reduce(UpperBox& domain, const VectorTaylorFunction& function, const ExactBox& bounds,
                                       ExactFloatVector centre, ExactFloatVector multipliers) const
{
    ScalarTaylorFunction g(function.domain(),default_sweeper());
    UpperInterval C(0);
    for(Nat i=0; i!=function.result_size(); ++i) {
        g += make_exact(multipliers[i]) * function[i];
        C += make_exact(multipliers[i]) * bounds[i];
    }
    Covector<UpperInterval> dg = gradient(g,domain);
    C -= g(centre);

    UpperBox new_domain(domain);
    UpperIntervalVector ranges(domain.size());
    for(Nat j=0; j!=domain.size(); ++j) {
        ranges[j] = dg[j]*(domain[j]-centre[j]);
    }

    // We now have sum dg(xi)[j] * (x[j]-x0[j]) in C, so we can reduce each component
    for(Nat j=0; j!=domain.size(); ++j) {
        UpperInterval E = C;
        for(Nat k=0; k!=domain.size(); ++k) {
            if(j!=k) { E-=ranges[k]; }
        }
        UpperInterval estimated_domain = E/dg[j]+centre[j];
        new_domain[j] = intersection(domain[j],estimated_domain);
    }

    domain=new_domain;
    return definitely(domain.empty());
}

Bool ConstraintSolver::box_reduce(UpperBox& domain, const ValidatedScalarFunctionInterface& function, const ExactInterval& bounds, Nat variable) const
{
    ARIADNE_LOG(2,"ConstraintSolver::box_reduce(ExactBox domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<"\n");

    if(domain[variable].lower() == domain[variable].upper()) { return false; }

    // Try to reduce the size of the set by "shaving" off along a coordinate axis
    //
    UpperInterval interval=domain[variable];
    RawFloat l=interval.lower().raw();
    RawFloat u=interval.upper().raw();
    ExactInterval subinterval;
    UpperInterval new_interval(interval);
    Vector<UpperInterval> slice=domain;

    static const Nat MAX_SLICES=(1<<3);
    const Nat n=MAX_SLICES;

    // Look for empty slices from below
    Nat imax = n;
    for(Nat i=0; i!=n; ++i) {
        subinterval=ExactInterval((l*(n-i)+u*i)/n,(l*(n-i-1)+u*(i+1))/n);
        slice[variable]=subinterval;
        UpperInterval slice_image=apply(function,slice);
        if(definitely(intersection(slice_image,bounds).empty())) {
            new_interval.set_lower(subinterval.upper());
        } else {
            imax = i; break;
        }
    }

    // The set is proved to be empty
    if(imax==n) {
        domain[variable]=ExactInterval(+inf,-inf);
        return true;
    }

    // Look for empty slices from above; note that at least one nonempty slice has been found
    for(Nat j=n-1; j!=imax; --j) {
        subinterval=ExactInterval((l*(n-j)+u*j)/n,(l*(n-j-1)+u*(j+1))/n);
        slice[variable]=subinterval;
        UpperInterval slice_image=apply(function,slice);
        if(definitely(intersection(slice_image,bounds).empty())) {
            new_interval.set_upper(subinterval.lower());
        } else {
            break;
        }
    }

    // The set cannot be empty, since a nonempty slice has been found in the upper pass.
    ARIADNE_ASSERT(new_interval.upper()!=new_interval.lower());

    domain[variable]=new_interval;

    return false;
}



namespace {

Void compute_monotonicity(UpperBox& domain, const EffectiveConstraint& constraint) {
/*
    static const Nat n = domain.size();

    // Compute monotone formulae
    boost::Array<Sign,0> monotonicity(n);
    Covector<ExactInterval> grad=constraint.function().gradient(domain);
    for(Nat j=0; j!=n; ++j) {
        monotonicity[j]=sign(grad[j]);
    }
*/
}
} // namespace


Pair<UpperBox,UpperBox> ConstraintSolver::split(const UpperBox& d, const ValidatedVectorFunction& f, const ExactBox& c) const
{
    return d.split();
}


Tribool ConstraintSolver::check_feasibility(const ExactBox& d, const ValidatedVectorFunction& f, const ExactBox& c, const ExactPoint& y) const
{
    for(Nat i=0; i!=y.size(); ++i) {
        if(y[i]<d[i].lower() || y[i]>d[i].upper()) { return false; }
    }

    Vector<ValidatedFloat> fy=f(Vector<ValidatedFloat>(y));
    ARIADNE_LOG(4,"d="<<d<<" f="<<f<<", c="<<c<<"\n  y="<<y<<", f(y)="<<fy<<"\n");
    Tribool result=true;
    for(Nat j=0; j!=fy.size(); ++j) {
        if(fy[j].lower().raw()>c[j].upper().raw() || fy[j].upper().raw()<c[j].lower().raw()) { return false; }
        if(fy[j].upper().raw()>=c[j].upper().raw() || fy[j].lower().raw()<=c[j].lower().raw()) { result=indeterminate; }
    }
    return result;
}






} // namespace Ariadne
