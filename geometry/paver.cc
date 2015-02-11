/***************************************************************************
 *            paver.cc
 *
 *  Copyright 2011-12  Pieter Collins
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

#include "geometry/paver.h"

#include "utility/macros.h"
#include "utility/logging.h"
#include "function/polynomial.h"
#include "function/function.h"
#include "function/taylor_function.h"
#include "function/procedure.h"
#include "geometry/function_set.h"
#include "geometry/affine_set.h"
#include "geometry/paving_interface.h"
#include "geometry/grid_set.h"
#include "solvers/nonlinear_programming.h"
#include "solvers/constraint_solver.h"
#include "geometry/affine_set.h"

namespace Ariadne {

Pair<Nat,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const ExactBox& domain);
Pair<Nat,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const UpperBox& domain) {
    return lipschitz_index_and_error(function,make_exact_box(domain));
}

namespace {

static const Nat verbosity = 0u;

Void subdivision_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                            const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth);

Void affine_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                       const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth);

Void constraint_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                           const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth);

Void optimal_constraint_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                                   const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth);

} // namespace

namespace {

ValidatedProcedure make_procedure(const ValidatedScalarFunction& f) {
    Formula<ValidatedNumber> e=f(Formula<ValidatedNumber>::identity(f.argument_size()));
    return Procedure<ValidatedNumber>(e);
}

UpperInterval emulrng(const ExactFloatVector& x, const ExactFloatVector& z) {
    UpperInterval r=make_interval(mul(x[0],z[0]));
    for(Nat i=0; i!=x.size(); ++i) { r=hull(mul(x[i],z[i]),r); }
    return r;
}

UpperInterval emulrng(const RawFloatVector& x, const RawFloatVector& z) {
    return emulrng(reinterpret_cast<ExactFloatVector const&>(x),reinterpret_cast<ExactFloatVector const&>(z));
}

UpperFloat64 total_widths(const UpperBox& bx) {
    UpperFloat64 res=0u;
    for(Nat i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width());
    }
    return res;
}

UpperFloat64 average_width(const UpperBox& bx) {
    UpperFloat64 res=0u;
    for(Nat i=0; i!=bx.size(); ++i) {
        if(definitely(bx[i].lower()>bx[i].upper())) { return -infty; }
        res+=bx[i].width();
    }
    return res/bx.size();
}

UpperFloat64 maximum_scaled_width(const UpperBox& bx, const Vector<ExactFloat64>& sf) {
    UpperFloat64 res=0u;
    for(Nat i=0; i!=bx.size(); ++i) {
        res=max(bx[i].width()/sf[i],res);
    }
    return res;
}

UpperFloat64 average_scaled_width(const UpperBox& bx, const Vector<ExactFloat64>& sf) {
    UpperFloat64 res=0u;
    for(Nat i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width()/sf[i]);
    }
    return res/bx.size();
}

Float64 maximum_scaled_width(const UpperBox& bx, const Vector<Float64>& sf) {
    return maximum_scaled_width(bx,reinterpret_cast<Vector<ExactFloat64>const&>(sf)).raw();
}

Float64 average_scaled_width(const UpperBox& bx, const Vector<Float64>& sf) {
    return average_scaled_width(bx,reinterpret_cast<Vector<ExactFloat64>const&>(sf)).raw();
}

} // namespace

Void SubdivisionPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                           const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::subdivision_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

Void AffinePaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                           const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::affine_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

Void ConstraintPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                           const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::constraint_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

Void OptimalConstraintPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                           const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::optimal_constraint_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}


namespace {

using Ariadne::verbosity;

Void subdivision_adjoin_outer_approximation_recursion(PavingInterface& paving, const ExactBox& subdomain, const ValidatedVectorFunction& function,
                                                      const List<ValidatedConstraint>& constraints, const Int depth, const RawFloatVector& errors)
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;

    for(List<ValidatedConstraint>::ConstIterator iter=constraints.begin();
        iter!=constraints.end(); ++iter)
    {
        UpperInterval constraint_range=apply(iter->function(),subdomain);
        if(constraint_range.lower() > iter->bounds().upper() || constraint_range.upper() < iter->bounds().lower() ) { return; }
    }

    UpperBox range=apply(function,subdomain);
    Bool small=true;
    for(Nat i=0; i!=range.size(); ++i) {
        if(range[i].width().raw()>errors[i]*RELATIVE_SMALLNESS*2) {
            small=false;
            break;
        }
    }

    if(small) {
        paving.adjoin_outer_approximation(make_exact_box(range),depth);
    } else {
        ExactBox subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=Ariadne::split(subdomain);
        subdivision_adjoin_outer_approximation_recursion(paving,subdomain1,function,constraints,depth,errors);
        subdivision_adjoin_outer_approximation_recursion(paving,subdomain2,function,constraints,depth,errors);
    }
}



static Nat COUNT_TESTS=0u;

// Adjoin an over-approximation to the solution of $f(dom)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
Void procedure_constraint_adjoin_outer_approximation_recursion(
        PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& f,
        const ValidatedVectorFunction& g, const ExactBox& codomain, const GridCell& cell, Int max_dpth, Nat splt, const List<ValidatedProcedure>& procedures)
{
    const Nat m=domain.size();
    const Nat nf=f.result_size();
    const Nat ng=g.result_size();

    const ExactBox& cell_box=cell.box();
    const RawFloatVector scalings=paving.grid().lengths();

    UpperBox bbox = apply(f,domain);

    Float64 domwdth = average_width(domain).raw();
    Float64 bbxwdth = average_scaled_width(bbox,paving.grid().lengths()).raw();
    Float64 clwdth = average_scaled_width(cell_box,paving.grid().lengths()).raw();

    ARIADNE_LOG(2,"\nconstraint_adjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"   splt="<<splt<<" dpth="<<cell.tree_depth()<<" max_dpth="<<max_dpth<<"\n");
    ARIADNE_LOG(2,"     domwdth="<<domwdth<<" bbxwdth="<<bbxwdth<<" clwdth="<<clwdth<<" dom="<<domain<<" bbox="<<bbox<<" cell="<<cell.box()<<"\n");

    ConstraintSolver constraint_solver;

    if(paving.superset(cell)) {
        ARIADNE_LOG(4,"  Cell is already a subset of paving\n");
        return;
    }

    ++COUNT_TESTS;

    // Try to prove disjointness
    const ExactBox& old_domain=domain;
    ExactBox new_domain=old_domain;
    Float64 olddomwdth = average_width(domain).raw();
    Float64 newdomwdth = olddomwdth;

    static const double ACCEPTABLE_REDUCTION_FACTOR = 0.75;


    // ExactBox reduction steps
    for(Nat i=0; i!=nf; ++i) {
        for(Nat j=0; j!=m; ++j) {
            constraint_solver.box_reduce(new_domain,f[i],cell_box[i],j);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
    }
    for(Nat i=0; i!=ng; ++i) {
        for(Nat j=0; j!=m; ++j) {
            constraint_solver.box_reduce(new_domain,g[i],codomain[i],j);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
    }
    newdomwdth=average_width(new_domain).raw();
    ARIADNE_LOG(6,"     domwdth="<<newdomwdth<<" olddomwdth="<<olddomwdth<<" dom="<<new_domain<<" box reduce\n");

    // Hull reduction steps
    do {
        olddomwdth=newdomwdth;
        for(Nat i=0; i!=nf; ++i) {
            constraint_solver.hull_reduce(new_domain,procedures[i],cell_box[i]);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n"); return; }
            //constraint_solver.hull_reduce(new_domain,f[i],cell_box[i]);
        }
        for(Nat i=0; i!=ng; ++i) {
            constraint_solver.hull_reduce(new_domain,procedures[nf+i],codomain[i]);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n"); return; }
            //constraint_solver.hull_reduce(new_domain,g[i],codomain[i]);
        }
        newdomwdth=average_width(new_domain).raw();
        ARIADNE_LOG(6,"     domwdth="<<newdomwdth<<" dom="<<new_domain<<"\n");
    } while( !new_domain.empty() && (newdomwdth < ACCEPTABLE_REDUCTION_FACTOR * olddomwdth) );

    ARIADNE_LOG(6,"new_domain="<<new_domain);


    domwdth = average_scaled_width(new_domain,RawFloatVector(new_domain.size(),1.0)).raw();
    bbox=apply(f,new_domain);
    bbxwdth=average_scaled_width(bbox,paving.grid().lengths()).raw();
    if(definitely(bbox.disjoint(cell_box)) || definitely(disjoint(apply(g,new_domain),codomain))) {
        ARIADNE_LOG(4,"  Proved disjointness using image of new domain\n");
        return;
    }

    ARIADNE_LOG(4,"                 domwdth="<<domwdth<<" bbxwdth="<<bbxwdth<<" clwdth="<<clwdth<<" dom="<<new_domain<<" bbox="<<bbox<<" cell="<<cell.box()<<"\n");

    // Decide whether to split cell or split domain by comparing size of
    // bounding box with the cell and splitting the larger.
    // It seems that a more efficient algorithm results if the domain
    // is only split if the bounding box is much larger, so we preferentiably
    // split the cell unless the bounding box is 4 times as large
    Float64 bbxmaxwdth = maximum_scaled_width(bbox,scalings).raw();
    Float64 clmaxwdth = maximum_scaled_width(cell_box,scalings).raw();

    if( (bbxmaxwdth > 4.0*clmaxwdth) || (cell.tree_depth()>=max_dpth && (bbxmaxwdth > clmaxwdth)) ) {
        Pair<Nat,double> lipsch = lipschitz_index_and_error(f,new_domain);
        ARIADNE_LOG(4,"  Splitting domain on coordinate "<<lipsch.first<<"\n");
        Pair<ExactBox,ExactBox> sd=new_domain.split(lipsch.first);
        procedure_constraint_adjoin_outer_approximation_recursion(paving, sd.first, f, g, codomain, cell, max_dpth, splt+1, procedures);
        procedure_constraint_adjoin_outer_approximation_recursion(paving, sd.second, f, g, codomain, cell, max_dpth, splt+1, procedures);
    } else if(cell.tree_depth()>=max_dpth) {
        ARIADNE_LOG(4,"  Adjoining cell "<<cell_box<<"\n");
        paving.adjoin(cell);
    } else {
        ARIADNE_LOG(4,"  Splitting cell "<<cell_box<<"\n");
        Pair<GridCell,GridCell> sb = cell.split();
        procedure_constraint_adjoin_outer_approximation_recursion(paving,make_exact_box(new_domain),f,g,codomain,sb.first, max_dpth, splt, procedures);
        procedure_constraint_adjoin_outer_approximation_recursion(paving,make_exact_box(new_domain),f,g,codomain,sb.second, max_dpth, splt, procedures);
    }


}



Void hotstarted_constraint_adjoin_outer_approximation_recursion(
    PavingInterface& r, const ExactBox& d, const ValidatedVectorFunction& f,
    const ValidatedVectorFunction& g, const ExactBox& c, const GridCell& b, ExactPoint x, ExactPoint y, Int e)
{
    Nat verbosity=0;

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const double TERR = -1.0/((1<<e)*1024.0);
    static const double XZMIN = 1.0/(1<<16);
    static const Float64 inf = Ariadne::inf;

    // Set up the classes used for constraint propagation and
    // optimisation using the Kuhn-Tucker conditions
    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;
    ValidatedVectorFunction fg=join(f,g);

    const Nat m=fg.argument_size();
    const Nat n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");
    ARIADNE_LOG(2,"  x0="<<x<<", y0="<<y<<"\n");

    ExactPoint z(x.size());
    ExactFloat64 t;

    Vector<ApproximateFloat64>& ax=reinterpret_cast<Vector<ApproximateFloat64>&>(x);
    Vector<ApproximateFloat64>& ay=reinterpret_cast<Vector<ApproximateFloat64>&>(y);
    Vector<ApproximateFloat64> az=reinterpret_cast<Vector<ApproximateFloat64>&>(z);
    ApproximateFloat64 at=reinterpret_cast<ApproximateFloat64&>(t);

    if(r.superset(b)) {
        ARIADNE_LOG(2,"  Cell already in set\n");
        return;
    }

    ExactBox bx=product(static_cast<const ExactBox&>(b.box()),static_cast<const ExactBox&>(c));

    ARIADNE_LOG(2,"  fg(d)="<<apply(fg,d)<<", bx="<<bx<<"\n");
    if(definitely(disjoint(apply(fg,d),bx))) {
        ARIADNE_LOG(2,"  Proved disjointness using direct evaluation\n");
        return;
    }


    // Relax x away from boundary
    optimiser.compute_tz(d,fg,bx,ay,at,az);
    ARIADNE_LOG(2,"  z0="<<az<<", t0="<<at<<"\n");
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG(4," t="<<at);
        //optimiser.linearised_feasibility_step(d,fg,bx,x,y,z,t);
        try {
            optimiser.feasibility_step(d,fg,bx,ax,ay,az,at);
        }
        catch(NearBoundaryOfFeasibleDomainException e) {
            break;
        }
        catch(std::runtime_error e) {
            ARIADNE_ERROR(""<<e.what()<<"\n");
            break;
        }
        ARIADNE_LOG(6,", x="<<ax<<", y="<<ay<<", z="<<az<<"\n");
        ARIADNE_LOG(6,"  x.z="<<emulrng(x,z)<<"\n");
        if(t>0) { break; }
        if(emulrng(x,z).upper()<XZMIN) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");
    ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<"\n");

    if(!(t.raw()<inf)) {
        ARIADNE_WARN("feasibility failed\n");
        char c; cin >> c;
        at=0;
        ay=midpoint(d);
        ax=ApproximateFloatVector(x.size(),1.0/x.size());
    }
    ax = ApproximateFloat64(1-XSIGMA)*ax + Vector<ApproximateFloat64>(x.size(),XSIGMA/x.size());

    //assert(t>=-1000);

    if(t.raw()<TERR) {

        // Probably disjoint, so try to prove this
        UpperBox nd=d;
        const ExactBox& domain=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        TrivialSweeper sweeper;
        EffectiveScalarFunction zero_function=EffectiveScalarFunction::zero(m);
        EffectiveVectorFunction identity_function=EffectiveVectorFunction::identity(m);
        ScalarTaylorFunction txg(domain,zero_function,sweeper);
        ValidatedFloat64 cnst=0;
        for(Nat j=0; j!=n; ++j) {
            txg = txg - (ValidatedFloat64(x[j])-ValidatedFloat64(x[n+j]))*ScalarTaylorFunction(domain,ValidatedScalarFunction(fg[j]),sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(Nat i=0; i!=m; ++i) {
            txg = txg - (ValidatedFloat64(x[2*n+i])-ValidatedFloat64(x[2*n+m+i]))*ScalarTaylorFunction(domain,ValidatedScalarFunction(identity_function[i]),sweeper);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        txg = ValidatedFloat64(cnst) + txg;

        ARIADNE_LOG(6,"    txg="<<txg<<"\n");

        ValidatedConstraint constraint=(txg>=0);

        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,txg,ExactInterval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.empty())) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(Nat i=0; i!=m; ++i) {
            solver.box_reduce(nd,txg,ExactInterval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(definitely(nd.empty())) { ARIADNE_LOG(2,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        solver.hull_reduce(nd,txg,ExactInterval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.empty())) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }
    }

    if(t<=0.0_exact && UpperBox(apply(f,d)).radius()>b.box().radius()) {
        ARIADNE_LOG(2,"  Splitting domain\n");
        Pair<ExactBox,ExactBox> sd=d.split();
        ax = ApproximateFloat64(1-XSIGMA)*ax + Vector<ApproximateFloat64>(ax.size(),XSIGMA/x.size());
        ay=midpoint(sd.first);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.first, f,g, c, b, x, y, e);
        ay = midpoint(sd.second);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.second, f,g, c, b, x, y, e);
        return;
    }

    if(t>0.0_exact) {
        ARIADNE_LOG(2," Intersection point: parameter="<<y<<"\n");
    }

    if(b.tree_depth()>=e*Int(b.dimension())) {
        ARIADNE_LOG(2,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(2,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        hotstarted_constraint_adjoin_outer_approximation_recursion(r,d,f,g,c,sb.first,x,y,e);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r,d,f,g,c,sb.second,x,y,e);
    }
}


Void hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(PavingInterface& r, const ExactBox& d, const VectorTaylorFunction& fg, const ExactBox& c, const GridCell& b, ExactPoint& x, ExactPoint& y, Int e)
{
    Sweeper sweeper = fg.sweeper();

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const ExactFloat64 TERR ( -1.0/((1<<e)*1024.0) );
    static const Float64 inf = Ariadne::inf;

    const Nat m=fg.argument_size();
    const Nat n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");

    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;

    ExactFloat64 t;
    ExactPoint z(x.size());

    ApproximateFloatVector& ax=reinterpret_cast<ApproximateFloatVector&>(x);
    ApproximateFloatVector& ay=reinterpret_cast<ApproximateFloatVector&>(y);
    ApproximateFloatVector& az=reinterpret_cast<ApproximateFloatVector&>(z);
    ApproximateFloat64& at=reinterpret_cast<ApproximateFloat64&>(t);

    if(r.superset(b)) {
        return;
    }

    ExactBox bx=product(static_cast<const ExactBox&>(b.box()),static_cast<const ExactBox&>(c));

    optimiser.compute_tz(d,fg,bx,ay,at,az);
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG(4," t="<<t);
        optimiser.linearised_feasibility_step(d,fg,bx,ax,ay,az,at);
        if(t>0) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        UpperBox nd=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ScalarTaylorFunction xg=ScalarTaylorFunction::zero(d,sweeper);
        ValidatedFloat64 cnst=0;
        for(Nat j=0; j!=n; ++j) {
            xg = xg - (x[j]-x[n+j])*ScalarTaylorFunction(d,fg[j],sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(Nat i=0; i!=m; ++i) {
            xg = xg - (x[2*n+i]-x[2*n+m+i])*ScalarTaylorFunction::coordinate(d,i,sweeper);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        xg = (cnst) + xg;

        ARIADNE_LOG(4,"    xg="<<xg<<"\n");


        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,xg,ExactInterval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.empty())) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(Nat i=0; i!=m; ++i) {
            solver.box_reduce(nd,xg,ExactInterval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(definitely(nd.empty())) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        //Pair<ExactBox,ExactBox> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<ExactBox,ExactBox> sd=split(d);
        ExactPoint nx = make_exact(ApproximateFloat64(1.0-XSIGMA)*ax + Vector<ApproximateFloat64>(x.size(),XSIGMA/x.size()));
        ExactPoint ny = midpoint(sd.first);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.first, fg, c, b, nx, ny, e);
        nx = make_exact(ApproximateFloat64(1.0-XSIGMA)*x + Vector<ApproximateFloat64>(x.size(),XSIGMA/x.size()));
        ny = midpoint(sd.second);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.tree_depth()>=e*Int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        ExactPoint sx = make_exact(ApproximateFloat64(1-XSIGMA)*x + Vector<ApproximateFloat64>(x.size(),XSIGMA/x.size()));
        ExactPoint sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.first,sx,sy,e);
        sx = make_exact(ApproximateFloat64(1-XSIGMA)*x + Vector<ApproximateFloat64>(x.size(),XSIGMA/x.size()));
        sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.second,sx,sy,e);
    }


}



Void subdivision_adjoin_outer_approximation(PavingInterface& paving,
                                            const ExactBox& subdomain,
                                            const ValidatedVectorFunction& function,
                                            const ValidatedVectorFunction& constraint_functions,
                                            const ExactBox& constraint_bounds,
                                            Int depth)
{
    List<ValidatedConstraint> constraints;
    for(Nat i=0; i!=constraint_functions.result_size(); ++i) {
        constraints.append(ValidatedConstraint(constraint_bounds[i].lower(),constraint_functions[i],constraint_bounds[i].upper()));
    }

    RawFloatVector errors(paving.dimension());
    for(Nat i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }

    ::subdivision_adjoin_outer_approximation_recursion(paving,subdomain,function,constraints,depth,errors);
}

Void affine_adjoin_outer_approximation(PavingInterface& paving,
                                       const ExactBox& subdomain,
                                       const ValidatedVectorFunction& function,
                                       const ValidatedVectorFunction& constraints,
                                       const ExactBox& bounds,
                                       Int depth)
{
    ARIADNE_NOT_IMPLEMENTED;
}

Void
constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBox& d, const ValidatedVectorFunction& f,
                                      const ValidatedVectorFunction& g, const ExactBox& c, Int e)
{
    ARIADNE_ASSERT(p.dimension()==f.result_size());

    GridCell b=GridCell::smallest_enclosing_primary_cell(make_exact_box(apply(f,d)),p.grid());
    ExactBox r=make_exact_box(apply(g,d)+UpperBox(g.result_size(),UpperInterval(-1,1)));
    ExactBox rc=intersection(r,c);

    ExactPoint y=midpoint(d);
    const Nat l=(d.size()+f.result_size()+g.result_size())*2;
    ExactPoint x(l); for(Nat k=0; k!=l; ++k) { x[k]=ExactFloat64(1.0/l); }

    ::hotstarted_constraint_adjoin_outer_approximation_recursion(p,d,f,g,rc,b,x,y,e);
}

Void
procedure_constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBox& d, const ValidatedVectorFunction& f,
                                                const ValidatedVectorFunction& g, const ExactBox& c, Int e)
{
    GridCell b=p.smallest_enclosing_primary_cell(make_exact_box(apply(f,d)));

    List<ValidatedProcedure> procedures;
    procedures.reserve(f.result_size()+g.result_size());
    for(Nat i=0; i!=f.result_size(); ++i) { procedures.append(make_procedure(f[i])); }
    for(Nat i=0; i!=g.result_size(); ++i) { procedures.append(make_procedure(g[i])); }

    Ariadne::procedure_constraint_adjoin_outer_approximation_recursion(p,d,f,g,c,b,e*p.dimension(),0, procedures);
    //std::cerr<<"Computing outer approximation considered a total of "<<COUNT_TESTS<<" domains/cells\n";
    //std::cerr<<"Measure of paving is "<<p.measure()<<"\n";

    if(dynamic_cast<GridTreeSet*>(&p)) { dynamic_cast<GridTreeSet&>(p).recombine(); }
}

Void optimal_constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBox& d, const ValidatedVectorFunction& f,
                                                   const ValidatedVectorFunction& g, const ExactBox& c, Int e)
{
    GridCell b=GridCell::smallest_enclosing_primary_cell(make_exact_box(apply(g,d)),p.grid());
    ExactBox rc=intersection(make_exact_box(apply(g,d)+UpperBox(g.result_size(),UpperInterval(-1,1))),c);

    ExactPoint y=midpoint(d);
    const Nat l=(d.size()+f.result_size()+g.result_size())*2;
    ExactPoint x(l); for(Nat k=0; k!=l; ++k) { x[k]=ExactFloat64(1.0/l); }

    VectorTaylorFunction fg;
    const VectorTaylorFunction* tfptr;
    if( (tfptr=dynamic_cast<const VectorTaylorFunction*>(f.raw_pointer())) ) {
        const VectorTaylorFunction* tgptr;
        if( ( tgptr = dynamic_cast<const VectorTaylorFunction*>(g.raw_pointer()) ) ) {
            fg=join(*tfptr,*tgptr);
        } else {
            if(g.result_size()>0) {
                fg=join(*tfptr,VectorTaylorFunction(tfptr->domain(),g,tfptr->sweeper()));
            } else {
                fg=*tfptr;
            }
        }
    } else {
        ThresholdSweeper swp(1e-12);
        fg=VectorTaylorFunction(d,join(f,g),swp);
    }
    ::hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(p,d,fg,rc,b,x,y,e);
}

} // namespace

} // namespace Ariadne


