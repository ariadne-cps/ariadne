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

#include "functional.h"
#include "config.h"

#include "macros.h"
#include "tuple.h"
#include "tribool.h"
#include "numeric.h"
#include "vector.h"
#include "box.h"
#include "grid_set.h"
#include "polynomial.h"
#include "function.h"
#include "formula.h"
#include "procedure.h"
#include "constraint.h"
#include "nonlinear_programming.h"
#include "function_mixin.h"
#include "taylor_function.h"

#include "constraint_solver.h"
#include "solver.h"

namespace Ariadne {

typedef Vector<ApproximateFloatType> ApproximateFloatVector;
typedef Vector<ExactFloatType> ExactFloatVector;

inline Sweeper default_sweeper() { return Sweeper(); }

static const double error =  1e-2;

Sign sign(const Float& x) {
    if(x>0) { return NEGATIVE; }
    else if(x<0) {  return POSITIVE; }
    else { return ZERO; }
}

Sign sign(const Interval& ivl) {
    if(ivl.lower()>0) { return NEGATIVE; }
    else if(ivl.upper()<0) {  return POSITIVE; }
    else { return ZERO; }
}


std::ostream& operator<<(std::ostream& os, const EffectiveConstraint& c) {
    static const Float inf = Ariadne::inf;
    if(c.bounds().lower()==c.bounds().upper()) { return os << c.function() << "==" << c.bounds().upper(); }
    if(c.bounds().upper()==infty) { return os << c.bounds().lower() << "<=" << c.function(); }
    if(c.bounds().lower()==-infty) { return os << c.function() << "<=" << c.bounds().upper(); }
    return os << c.bounds().lower() << "<=" << c.function() << "<=" << c.bounds().upper();
}



Pair<Tribool,ExactPoint> ConstraintSolver::feasible(const Box& domain, const List<ValidatedConstraint>& constraints) const
{
    if(constraints.empty()) { return make_pair(!domain.empty(),domain.centre()); }

    ValidatedVectorFunction function(constraints.size());
    Box bounds(constraints.size());

    for(uint i=0; i!=constraints.size(); ++i) {
        function[i]=constraints[i].function();
        bounds[i]=constraints[i].bounds();
    }
    return this->feasible(domain,function,bounds);
}


Pair<Tribool,ExactPoint> ConstraintSolver::feasible(const Box& domain, const ValidatedVectorFunction& function, const Box& codomain) const
{

    static const double XSIGMA=0.125;
    static const double TERR=-1.0/(1<<10);
    static const Float inf = Ariadne::inf;

    ARIADNE_LOG(4,"domain="<<domain<<"\nfunction="<<function<<"\ncodomain="<<codomain<<"\n");

    // Make codomain bounded
    IntervalVector bounds=codomain;
    IntervalVector image=apply(function,domain);
    ARIADNE_LOG(4,"image="<<image<<"\n");
    for(uint i=0; i!=image.size(); ++i) {
        if(disjoint(image[i],codomain[i])) { ARIADNE_LOG(4,"  Proved disjointness using direct evaluation\n");  return make_pair(false,ExactPoint()); }
        else bounds[i]=intersection(codomain[i],image[i]);
    }


    const uint m=domain.size(); // The total number of variables
    const uint n=codomain.size(); // The total number of nontrivial constraints
    const uint l=(m+n)*2; // The total number of lagrange multipliers

    ApproximateFloatVector point(m); // The point in the domain which is the current test point
    ApproximateFloat violation; // An upper bound on amount by which the constraints are violated by the test point
    ApproximateFloatVector multipliers(l); // The lagrange multipliers for the constraints
    ApproximateFloatVector slack(l); // The slack between the test point and the violated constraints

    ApproximateFloat& t=violation; ApproximateFloatVector& x=multipliers; ApproximateFloatVector& y=point; ApproximateFloatVector& z=slack; // Aliases for the main quantities used
    const Box& d=domain; const ValidatedVectorFunction& fn=function; const Box& c=codomain; // Aliases for the main quantities used
    VectorTaylorFunction tfn(d,fn,default_sweeper());

    point=static_cast<ApproximateFloatVector>(midpoint(d));
    for(uint k=0; k!=l; ++k) { multipliers[k]=1.0/l; }

    NonlinearInteriorPointOptimiser optimiser;
    optimiser.compute_tz(domain,function,bounds,point,violation,slack);

    ARIADNE_LOG(4,"d="<<d<<", f="<<fn<<", c="<<c<<"\n");


    // TODO: Don't use fixed number of steps
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"    t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
        optimiser.feasibility_step(d,fn,c,x,y,z,t);
        if(t>=TERR) {
            ARIADNE_LOG(4,"t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");
            if(this->check_feasibility(domain,function,codomain,make_exact(point))) { return make_pair(true,make_exact(point)); }
            else { ARIADNE_LOG(2,"f(y)="<<fn(make_exact(y))<<"\n"); return make_pair(indeterminate,make_exact(point)); }
        }
    }
    ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", x="<<x<<", z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        Box subdomain=domain;

        Vector<ExactFloat> x_exact=make_exact(x);
        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ScalarTaylorFunction txg=ScalarTaylorFunction::zero(d,default_sweeper());
        ValidatedNumberType cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            txg = txg - (x_exact[j]-x_exact[n+j])*tfn[j];
            cnst += (c[j].upper()*x_exact[j]-c[j].lower()*x_exact[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            txg = txg - (x_exact[2*n+i]-x_exact[2*n+m+i])*ScalarTaylorFunction::coordinate(d,i,default_sweeper());
            cnst += (d[i].upper()*x_exact[2*n+i]-d[i].lower()*x_exact[2*n+m+i]);
        }
        txg = cnst + txg;

        ARIADNE_LOG(4,"    txg="<<txg<<"\n");

        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        this->hull_reduce(subdomain,txg,Interval(0,inf));
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");
        if(subdomain.empty()) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return make_pair(false,ExactPoint());
        }

        for(uint i=0; i!=m; ++i) {
            this->box_reduce(subdomain,txg,Interval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<subdomain<<"\n");
            if(subdomain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return make_pair(false,ExactPoint()); }
        }
        ARIADNE_LOG(6,"  dom="<<subdomain<<"\n");

        //Pair<Box,Box> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=d.split();
        Vector<ApproximateFloatType> nx = ApproximateFloatType(1.0-XSIGMA)*x + Vector<ApproximateFloatType>(x.size(),XSIGMA/x.size());
        Vector<ApproximateFloatType> ny = midpoint(sd.first);
        tribool result=this->feasible(sd.first, fn, c).first;
        nx = ApproximateFloatType(1.0-XSIGMA)*x + Vector<ApproximateFloatType>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        result = result || this->feasible(sd.second, fn, c).first;
        return make_pair(result,ExactPoint());
    }

    return make_pair(indeterminate,ExactPoint());
}


bool ConstraintSolver::reduce(Box& domain, const ValidatedVectorFunction& function, const Box& codomain) const
{
    const double MINIMUM_REDUCTION = 0.75;
    ARIADNE_ASSERT(function.argument_size()==domain.size());
    ARIADNE_ASSERT(function.result_size()==codomain.size());

    if(domain.empty()) { return true; }

    Float domain_magnitude=0.0;
    for(uint j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width().raw();
    }
    Float old_domain_magnitude=domain_magnitude;

    do {
        this->hull_reduce(domain,function,codomain);
        if(domain.empty()) { return true; }

        for(uint i=0; i!=codomain.size(); ++i) {
            for(uint j=0; j!=domain.size(); ++j) {
                this->box_reduce(domain,function[i],codomain[i],j);
		        if(domain.empty()) { return true; }
            }
        }
        if(domain.empty()) { return true; }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0.0;
        for(uint j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width().raw();
        }
    } while(domain_magnitude/old_domain_magnitude <= MINIMUM_REDUCTION);

    return false;
}

inline bool is_nan(const Float& x) { return isnan(x.get_d()); }

bool has_nan(const Box& domain) {
    for(uint i=0; i!=domain.size(); ++i) {
        if(is_nan(domain[i].lower().raw()) || is_nan(domain[i].upper().raw())) { return true; }
    }
    return false;
}

bool ConstraintSolver::reduce(Box& domain, const List<ValidatedConstraint>& constraints) const
{
    static const bool USE_BOX_REDUCE = false;

    const double MINIMUM_REDUCTION = 0.75;

    if(domain.empty()) { return true; }

    Float domain_magnitude=0.0;
    for(uint j=0; j!=domain.size(); ++j) {
        domain_magnitude+=domain[j].width().raw();
    }
    Float old_domain_magnitude=domain_magnitude;

    do {
        for(uint i=0; i!=constraints.size(); ++i) {
            this->hull_reduce(domain,constraints[i].function(),constraints[i].bounds());
        }
        if(domain.empty()) { return true; }

        if(USE_BOX_REDUCE) {
            for(uint i=0; i!=constraints.size(); ++i) {
                for(uint j=0; j!=domain.size(); ++j) {
                    this->box_reduce(domain,constraints[i].function(),constraints[i].bounds(),j);
                    if(domain[j].empty()) { return true; }
                }
            }
        }

        old_domain_magnitude=domain_magnitude;
        domain_magnitude=0.0;
        for(uint j=0; j!=domain.size(); ++j) {
            domain_magnitude+=domain[j].width().raw();
        }
    } while(domain_magnitude/old_domain_magnitude <= MINIMUM_REDUCTION);

    return false;
}


bool ConstraintSolver::hull_reduce(Box& domain, const ValidatedProcedure& procedure, const Interval& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, ValidatedProcedure procedure, Interval bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return domain.empty();
}

bool ConstraintSolver::hull_reduce(Box& domain, const Vector<ValidatedProcedure>& procedure, const Box& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, Vector<ValidatedProcedure> procedure, Box bounds): "
                  "procedure="<<procedure<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Ariadne::simple_hull_reduce(domain, procedure, bounds);
    return domain.empty();
}

bool ConstraintSolver::hull_reduce(Box& domain, const ValidatedScalarFunctionInterface& function, const Interval& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, ValidatedScalarFunction function, Interval bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Formula<ValidatedNumberType> formula=function.evaluate(Formula<ValidatedNumberType>::identity(function.argument_size()));
    Procedure<ValidatedNumberType> procedure(formula);
    return this->hull_reduce(domain,procedure,bounds);
}

bool ConstraintSolver::hull_reduce(Box& domain, const ValidatedVectorFunctionInterface& function, const Box& bounds) const
{
    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain, ValidatedScalarFunction function, Interval bounds): "
                  "function="<<function<<", bounds="<<bounds<<", domain="<<domain<<"\n");

    Vector< Formula<ValidatedNumberType> > formula=function.evaluate(Formula<ValidatedNumberType>::identity(function.argument_size()));
    Vector< Procedure<ValidatedNumberType> > procedure(formula);
    return this->hull_reduce(domain,procedure,bounds);
}

bool ConstraintSolver::monotone_reduce(Box& domain, const ValidatedScalarFunctionInterface& function, const Interval& bounds, uint variable) const
{
    ValidatedScalarFunction derivative=function.derivative(variable);

    ARIADNE_LOG(2,"ConstraintSolver::hull_reduce(Box domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<", derivative="<<derivative<<"\n");

    ExactFloatType splitpoint;
    Interval lower=domain[variable];
    Interval upper=domain[variable];
    Vector<Interval> slice=domain;
    Vector<Interval> subdomain=domain;

    static const int MAX_STEPS=3;
    const Float size = lower.width().raw() / (1<<MAX_STEPS);
    do {
        if(lower.width().raw()>size) {
            splitpoint=lower.centre();
            slice[variable]=splitpoint;
            Interval new_lower=splitpoint+(bounds-apply(function,slice)).lower()/apply(derivative,subdomain);
            if(new_lower.upper().raw()<lower.lower().raw()) { lower=lower.lower().raw(); }
            else { lower=intersection(lower,new_lower); }
        }
        if(upper.width().raw()>size) {
            splitpoint=upper.centre();
            slice[variable]=splitpoint;
            Interval new_upper=splitpoint+(bounds-apply(function,slice)).upper()/apply(derivative,subdomain);
            if(new_upper.lower()>upper.upper()) { upper=upper.upper(); }
            else { upper=intersection(upper,new_upper); }
        }
        subdomain[variable]=Interval(lower.lower(),upper.upper());
    } while(lower.width().raw()>size && upper.width().raw()>size);
    domain=subdomain;

    return domain.empty();
}



bool ConstraintSolver::lyapunov_reduce(Box& domain, const VectorTaylorFunction& function, const Box& bounds,
                                       ApproximateFloatVector centre, ApproximateFloatVector multipliers) const
{
    return this->lyapunov_reduce(domain,function,bounds,make_exact(centre),make_exact(multipliers));
}


bool ConstraintSolver::lyapunov_reduce(Box& domain, const VectorTaylorFunction& function, const Box& bounds,
                                       ExactFloatVector centre, ExactFloatVector multipliers) const
{
    ScalarTaylorFunction g(function.domain(),default_sweeper());
    Interval C(0);
    for(uint i=0; i!=function.result_size(); ++i) {
        g += make_exact(multipliers[i]) * function[i];
        C += make_exact(multipliers[i]) * bounds[i];
    }
    IntervalVector dg = gradient(g,domain);
    C -= g(centre);

    Box new_domain(domain);
    IntervalVector ranges(domain.size());
    for(uint j=0; j!=domain.size(); ++j) {
        ranges[j] = dg[j]*(domain[j]-centre[j]);
    }

    // We now have sum dg(xi)[j] * (x[j]-x0[j]) in C, so we can reduce each component
    for(uint j=0; j!=domain.size(); ++j) {
        Interval E = C;
        for(uint k=0; k!=domain.size(); ++k) {
            if(j!=k) { E-=ranges[k]; }
        }
        Interval estimated_domain = E/dg[j]+centre[j];
        new_domain[j] = intersection(domain[j],estimated_domain);
    }

    domain=new_domain;
    return domain.empty();
}

bool ConstraintSolver::box_reduce(Box& domain, const ValidatedScalarFunctionInterface& function, const Interval& bounds, uint variable) const
{
    ARIADNE_LOG(2,"ConstraintSolver::box_reduce(Box domain): function="<<function<<", bounds="<<bounds<<", domain="<<domain<<", variable="<<variable<<"\n");

    if(domain[variable].lower() == domain[variable].upper()) { return false; }

    // Try to reduce the size of the set by "shaving" off along a coordinate axis
    //
    Interval interval=domain[variable];
    RawFloatType l=interval.lower().raw();
    RawFloatType u=interval.upper().raw();
    Interval subinterval;
    Interval new_interval(interval);
    Vector<Interval> slice=domain;

    static const uint MAX_SLICES=(1<<3);
    const uint n=MAX_SLICES;

    // Look for empty slices from below
    uint imax = n;
    for(uint i=0; i!=n; ++i) {
        subinterval=Interval((l*(n-i)+u*i)/n,(l*(n-i-1)+u*(i+1))/n);
        slice[variable]=subinterval;
        Interval slice_image=apply(function,slice);
        if(intersection(slice_image,bounds).empty()) {
            new_interval.set_lower(subinterval.upper());
        } else {
            imax = i; break;
        }
    }

    // The set is proved to be empty
    if(imax==n) {
        domain[variable]=Interval(+inf,-inf);
        return true;
    }

    // Look for empty slices from above; note that at least one nonempty slice has been found
    for(uint j=n-1; j!=imax; --j) {
        subinterval=Interval((l*(n-j)+u*j)/n,(l*(n-j-1)+u*(j+1))/n);
        slice[variable]=subinterval;
        Interval slice_image=apply(function,slice);
        if(intersection(slice_image,bounds).empty()) {
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

void compute_monotonicity(Box& domain, const EffectiveConstraint& constraint) {
/*
    static const uint n = domain.size();

    // Compute monotone formulae
    boost::Array<Sign,0> monotonicity(n);
    Vector<Interval> grad=constraint.function().gradient(domain);
    for(uint j=0; j!=n; ++j) {
        monotonicity[j]=sign(grad[j]);
    }
*/
}
} // namespace


Pair<Box,Box> ConstraintSolver::split(const Box& d, const ValidatedVectorFunction& f, const Box& c) const
{
    return d.split();
}


tribool ConstraintSolver::check_feasibility(const Box& d, const ValidatedVectorFunction& f, const Box& c, const ExactPoint& y) const
{
    for(uint i=0; i!=y.size(); ++i) {
        if(y[i]<d[i].lower() || y[i]>d[i].upper()) { return false; }
    }

    Vector<ValidatedFloatType> fy=f(Vector<ValidatedFloatType>(y));
    ARIADNE_LOG(4,"d="<<d<<" f="<<f<<", c="<<c<<"\n  y="<<y<<", f(y)="<<fy<<"\n");
    tribool result=true;
    for(uint j=0; j!=fy.size(); ++j) {
        if(fy[j].lower()>c[j].upper() || fy[j].upper()<c[j].lower()) { return false; }
        if(fy[j].upper()>=c[j].upper() || fy[j].lower()<=c[j].lower()) { result=indeterminate; }
    }
    return result;
}






} // namespace Ariadne
