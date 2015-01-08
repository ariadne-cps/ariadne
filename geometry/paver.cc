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

Pair<uint,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const ExactBox& domain);
Pair<uint,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const UpperBox& domain) {
    return lipschitz_index_and_error(function,make_exact_box(domain));
}

namespace {

static const uint verbosity = 0u;

void subdivision_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                            const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, int depth);

void affine_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                       const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, int depth);

void constraint_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                           const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, int depth);

void optimal_constraint_adjoin_outer_approximation(PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& function,
                                                   const ValidatedVectorFunction& constraint_function, const ExactBox& constraint_bounds, int depth);

} // namespace

namespace {

ValidatedProcedure make_procedure(const ValidatedScalarFunctionInterface& f) {
    Formula<ValidatedNumber> e=f.evaluate(Formula<ValidatedNumber>::identity(f.argument_size()));
    return Procedure<ValidatedNumber>(e);
}

UpperInterval emulrng(const RawFloatVector& x, const RawFloatVector& z) {
    UpperInterval r=mul_ivl(x[0],z[0]);
    for(uint i=0; i!=x.size(); ++i) { r=hull(mul_ivl(x[i],z[i]),r); }
    return r;
}

UpperInterval emulrng(const ExactFloatVector& x, const ExactFloatVector& z) {
    return emulrng(reinterpret_cast<Vector<Float>const&>(x),reinterpret_cast<Vector<Float>const&>(z));
}

UpperFloat total_widths(const UpperBox& bx) {
    UpperFloat res=0u;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width());
    }
    return res;
}

UpperFloat average_width(const UpperBox& bx) {
    UpperFloat res=0u;
    for(uint i=0; i!=bx.size(); ++i) {
        if(bx[i].lower().raw()>bx[i].upper().raw()) { return -infty; }
        res+=bx[i].width();
    }
    return res/bx.size();
}

UpperFloat maximum_scaled_width(const UpperBox& bx, const Vector<ExactFloat>& sf) {
    UpperFloat res=0u;
    for(uint i=0; i!=bx.size(); ++i) {
        res=max(bx[i].width()/sf[i],res);
    }
    return res;
}

UpperFloat average_scaled_width(const UpperBox& bx, const Vector<ExactFloat>& sf) {
    UpperFloat res=0u;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width()/sf[i]);
    }
    return res/bx.size();
}

Float maximum_scaled_width(const UpperBox& bx, const Vector<Float>& sf) {
    return maximum_scaled_width(bx,reinterpret_cast<Vector<ExactFloat>const&>(sf)).raw();
}

Float average_scaled_width(const UpperBox& bx, const Vector<Float>& sf) {
    return average_scaled_width(bx,reinterpret_cast<Vector<ExactFloat>const&>(sf)).raw();
}

} // namespace

void SubdivisionPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::subdivision_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

void AffinePaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::affine_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

void ConstraintPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::constraint_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

void OptimalConstraintPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const ExactBox& constraint_bounds, Int depth) const
{
    return Ariadne::optimal_constraint_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}


namespace {

using Ariadne::verbosity;

void subdivision_adjoin_outer_approximation_recursion(PavingInterface& paving, const ExactBox& subdomain, const ValidatedVectorFunction& function,
                                                      const List<ValidatedConstraint>& constraints, const int depth, const RawFloatVector& errors)
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
    bool small=true;
    for(uint i=0; i!=range.size(); ++i) {
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



static uint COUNT_TESTS=0u;

// Adjoin an over-approximation to the solution of $f(dom)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
void procedure_constraint_adjoin_outer_approximation_recursion(
        PavingInterface& paving, const ExactBox& domain, const ValidatedVectorFunction& f,
        const ValidatedVectorFunction& g, const ExactBox& codomain, const GridCell& cell, int max_dpth, uint splt, const List<ValidatedProcedure>& procedures)
{
    const uint m=domain.size();
    const uint nf=f.result_size();
    const uint ng=g.result_size();

    const ExactBox& cell_box=cell.box();
    const RawFloatVector scalings=paving.grid().lengths();

    UpperBox bbox = apply(f,domain);

    Float domwdth = average_width(domain).raw();
    Float bbxwdth = average_scaled_width(bbox,paving.grid().lengths()).raw();
    Float clwdth = average_scaled_width(cell_box,paving.grid().lengths()).raw();

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
    Float olddomwdth = average_width(domain).raw();
    Float newdomwdth = olddomwdth;

    static const double ACCEPTABLE_REDUCTION_FACTOR = 0.75;


    // ExactBox reduction steps
    for(uint i=0; i!=nf; ++i) {
        for(uint j=0; j!=m; ++j) {
            constraint_solver.box_reduce(new_domain,f[i],cell_box[i],j);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
    }
    for(uint i=0; i!=ng; ++i) {
        for(uint j=0; j!=m; ++j) {
            constraint_solver.box_reduce(new_domain,g[i],codomain[i],j);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
    }
    newdomwdth=average_width(new_domain).raw();
    ARIADNE_LOG(6,"     domwdth="<<newdomwdth<<" olddomwdth="<<olddomwdth<<" dom="<<new_domain<<" box reduce\n");

    // Hull reduction steps
    do {
        olddomwdth=newdomwdth;
        for(uint i=0; i!=nf; ++i) {
            constraint_solver.hull_reduce(new_domain,procedures[i],cell_box[i]);
            if(new_domain.empty()) { ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n"); return; }
            //constraint_solver.hull_reduce(new_domain,f[i],cell_box[i]);
        }
        for(uint i=0; i!=ng; ++i) {
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
    Float bbxmaxwdth = maximum_scaled_width(bbox,scalings).raw();
    Float clmaxwdth = maximum_scaled_width(cell_box,scalings).raw();

    if( (bbxmaxwdth > 4.0*clmaxwdth) || (cell.tree_depth()>=max_dpth && (bbxmaxwdth > clmaxwdth)) ) {
        Pair<uint,double> lipsch = lipschitz_index_and_error(f,new_domain);
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



void hotstarted_constraint_adjoin_outer_approximation_recursion(
    PavingInterface& r, const ExactBox& d, const ValidatedVectorFunction& f,
    const ValidatedVectorFunction& g, const ExactBox& c, const GridCell& b, ExactPoint x, ExactPoint y, int e)
{
    uint verbosity=0;

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const double TERR = -1.0/((1<<e)*1024.0);
    static const double XZMIN = 1.0/(1<<16);
    static const Float inf = Ariadne::inf;

    // Set up the classes used for constraint propagation and
    // optimisation using the Kuhn-Tucker conditions
    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;
    ValidatedVectorFunction fg=join(f,g);

    const uint m=fg.argument_size();
    const uint n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");
    ARIADNE_LOG(2,"  x0="<<x<<", y0="<<y<<"\n");

    ExactPoint z(x.size());
    ExactFloat t;

    Vector<ApproximateFloat>& ax=reinterpret_cast<Vector<ApproximateFloat>&>(x);
    Vector<ApproximateFloat>& ay=reinterpret_cast<Vector<ApproximateFloat>&>(y);
    Vector<ApproximateFloat> az=reinterpret_cast<Vector<ApproximateFloat>&>(z);
    ApproximateFloat at=reinterpret_cast<ApproximateFloat&>(t);

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
    for(uint i=0; i!=12; ++i) {
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

    if(!(t<=1e10)) {
        ARIADNE_WARN("feasibility failed\n");
        char c; cin >> c;
        at=0;
        ay=midpoint(d);
        ax=ApproximateFloatVector(x.size(),1.0/x.size());
    }
    ax = ApproximateFloat(1-XSIGMA)*ax + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size());

    //assert(t>=-1000);

    if(t<TERR) {

        // Probably disjoint, so try to prove this
        UpperBox nd=d;
        const ExactBox& domain=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        TrivialSweeper sweeper;
        EffectiveScalarFunction zero_function=EffectiveScalarFunction::zero(m);
        EffectiveVectorFunction identity_function=EffectiveVectorFunction::identity(m);
        ScalarTaylorFunction txg(domain,zero_function,sweeper);
        ValidatedFloat cnst=0;
        for(uint j=0; j!=n; ++j) {
            txg = txg - (ValidatedFloat(x[j])-ValidatedFloat(x[n+j]))*ScalarTaylorFunction(domain,ValidatedScalarFunction(fg[j]),sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            txg = txg - (ValidatedFloat(x[2*n+i])-ValidatedFloat(x[2*n+m+i]))*ScalarTaylorFunction(domain,ValidatedScalarFunction(identity_function[i]),sweeper);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        txg = ValidatedFloat(cnst) + txg;

        ARIADNE_LOG(6,"    txg="<<txg<<"\n");

        ValidatedConstraint constraint=(txg>=0);

        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,txg,ExactInterval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.empty())) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(uint i=0; i!=m; ++i) {
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

    if(t<=0.0 && UpperBox(apply(f,d)).radius()>b.box().radius()) {
        ARIADNE_LOG(2,"  Splitting domain\n");
        Pair<ExactBox,ExactBox> sd=d.split();
        ax = ApproximateFloat(1-XSIGMA)*ax + Vector<ApproximateFloat>(ax.size(),XSIGMA/x.size());
        ay=midpoint(sd.first);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.first, f,g, c, b, x, y, e);
        ay = midpoint(sd.second);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.second, f,g, c, b, x, y, e);
        return;
    }

    if(t>0.0) {
        ARIADNE_LOG(2," Intersection point: parameter="<<y<<"\n");
    }

    if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(2,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(2,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        hotstarted_constraint_adjoin_outer_approximation_recursion(r,d,f,g,c,sb.first,x,y,e);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r,d,f,g,c,sb.second,x,y,e);
    }
}


void hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(PavingInterface& r, const ExactBox& d, const VectorTaylorFunction& fg, const ExactBox& c, const GridCell& b, ExactPoint& x, ExactPoint& y, int e)
{
    Sweeper sweeper = fg.sweeper();

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const double TERR = -1.0/((1<<e)*1024.0);
    static const Float inf = Ariadne::inf;

    const uint m=fg.argument_size();
    const uint n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");

    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;

    ExactFloat t;
    ExactPoint z(x.size());

    ApproximateFloatVector& ax=reinterpret_cast<ApproximateFloatVector&>(x);
    ApproximateFloatVector& ay=reinterpret_cast<ApproximateFloatVector&>(y);
    ApproximateFloatVector& az=reinterpret_cast<ApproximateFloatVector&>(z);
    ApproximateFloat& at=reinterpret_cast<ApproximateFloat&>(t);

    if(r.superset(b)) {
        return;
    }

    ExactBox bx=product(static_cast<const ExactBox&>(b.box()),static_cast<const ExactBox&>(c));

    optimiser.compute_tz(d,fg,bx,ay,at,az);
    for(uint i=0; i!=12; ++i) {
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
        ValidatedFloat cnst=0;
        for(uint j=0; j!=n; ++j) {
            xg = xg - (x[j]-x[n+j])*ScalarTaylorFunction(d,fg[j],sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
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

        for(uint i=0; i!=m; ++i) {
            solver.box_reduce(nd,xg,ExactInterval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(definitely(nd.empty())) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        //Pair<ExactBox,ExactBox> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<ExactBox,ExactBox> sd=split(d);
        ExactPoint nx = make_exact(ApproximateFloat(1.0-XSIGMA)*ax + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size()));
        ExactPoint ny = midpoint(sd.first);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.first, fg, c, b, nx, ny, e);
        nx = make_exact(ApproximateFloat(1.0-XSIGMA)*x + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size()));
        ny = midpoint(sd.second);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        ExactPoint sx = make_exact(ApproximateFloat(1-XSIGMA)*x + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size()));
        ExactPoint sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.first,sx,sy,e);
        sx = make_exact(ApproximateFloat(1-XSIGMA)*x + Vector<ApproximateFloat>(x.size(),XSIGMA/x.size()));
        sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.second,sx,sy,e);
    }


}



void subdivision_adjoin_outer_approximation(PavingInterface& paving,
                                            const ExactBox& subdomain,
                                            const ValidatedVectorFunction& function,
                                            const ValidatedVectorFunction& constraint_functions,
                                            const ExactBox& constraint_bounds,
                                            int depth)
{
    List<ValidatedConstraint> constraints;
    for(uint i=0; i!=constraint_functions.result_size(); ++i) {
        constraints.append(ValidatedConstraint(constraint_bounds[i].lower(),constraint_functions[i],constraint_bounds[i].upper()));
    }

    RawFloatVector errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }

    ::subdivision_adjoin_outer_approximation_recursion(paving,subdomain,function,constraints,depth,errors);
}

void affine_adjoin_outer_approximation(PavingInterface& paving,
                                       const ExactBox& subdomain,
                                       const ValidatedVectorFunction& function,
                                       const ValidatedVectorFunction& constraints,
                                       const ExactBox& bounds,
                                       int depth)
{
    ARIADNE_NOT_IMPLEMENTED;
}

void
constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBox& d, const ValidatedVectorFunction& f,
                                      const ValidatedVectorFunction& g, const ExactBox& c, int e)
{
    ARIADNE_ASSERT(p.dimension()==f.result_size());

    GridCell b=GridCell::smallest_enclosing_primary_cell(make_exact_box(apply(f,d)),p.grid());
    ExactBox r=make_exact_box(apply(g,d)+UpperBox(g.result_size(),UpperInterval(-1,1)));
    ExactBox rc=intersection(r,c);

    ExactPoint y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    ExactPoint x(l); for(uint k=0; k!=l; ++k) { x[k]=ExactFloat(1.0/l); }

    ::hotstarted_constraint_adjoin_outer_approximation_recursion(p,d,f,g,rc,b,x,y,e);
}

void
procedure_constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBox& d, const ValidatedVectorFunction& f,
                                                const ValidatedVectorFunction& g, const ExactBox& c, int e)
{
    GridCell b=p.smallest_enclosing_primary_cell(make_exact_box(apply(f,d)));

    List<ValidatedProcedure> procedures;
    procedures.reserve(f.result_size()+g.result_size());
    for(uint i=0; i!=f.result_size(); ++i) { procedures.append(make_procedure(f[i])); }
    for(uint i=0; i!=g.result_size(); ++i) { procedures.append(make_procedure(g[i])); }

    Ariadne::procedure_constraint_adjoin_outer_approximation_recursion(p,d,f,g,c,b,e*p.dimension(),0, procedures);
    //std::cerr<<"Computing outer approximation considered a total of "<<COUNT_TESTS<<" domains/cells\n";
    //std::cerr<<"Measure of paving is "<<p.measure()<<"\n";

    if(dynamic_cast<GridTreeSet*>(&p)) { dynamic_cast<GridTreeSet&>(p).recombine(); }
}

void optimal_constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBox& d, const ValidatedVectorFunction& f,
                                                   const ValidatedVectorFunction& g, const ExactBox& c, int e)
{
    GridCell b=GridCell::smallest_enclosing_primary_cell(make_exact_box(apply(g,d)),p.grid());
    ExactBox rc=intersection(make_exact_box(apply(g,d)+UpperBox(g.result_size(),UpperInterval(-1,1))),c);

    ExactPoint y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    ExactPoint x(l); for(uint k=0; k!=l; ++k) { x[k]=ExactFloat(1.0/l); }

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


