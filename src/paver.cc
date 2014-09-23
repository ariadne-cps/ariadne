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

#include "functional.h"
#include "config.h"

#include "paver.h"

#include "macros.h"
#include "logging.h"
#include "polynomial.h"
#include "function.h"
#include "taylor_function.h"
#include "procedure.h"
#include "function_set.h"
#include "affine_set.h"
#include "paving_interface.h"
#include "grid_set.h"
#include "nonlinear_programming.h"
#include "constraint_solver.h"
#include "affine_set.h"

namespace Ariadne {

Pair<uint,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const Box& domain);

namespace {

static const uint verbosity = 0u;

void subdivision_adjoin_outer_approximation(PavingInterface& paving, const Box& domain, const ValidatedVectorFunction& function,
                                            const ValidatedVectorFunction& constraint_function, const Box& constraint_bounds, int depth);

void affine_adjoin_outer_approximation(PavingInterface& paving, const Box& domain, const ValidatedVectorFunction& function,
                                       const ValidatedVectorFunction& constraint_function, const Box& constraint_bounds, int depth);

void constraint_adjoin_outer_approximation(PavingInterface& paving, const Box& domain, const ValidatedVectorFunction& function,
                                           const ValidatedVectorFunction& constraint_function, const Box& constraint_bounds, int depth);

void optimal_constraint_adjoin_outer_approximation(PavingInterface& paving, const Box& domain, const ValidatedVectorFunction& function,
                                                   const ValidatedVectorFunction& constraint_function, const Box& constraint_bounds, int depth);

} // namespace

namespace {

ValidatedProcedure make_procedure(const ValidatedScalarFunctionInterface& f) {
    Formula<Interval> e=f.evaluate(Formula<Interval>::identity(f.argument_size()));
    return Procedure<Interval>(e);
}

Interval emulrng(const FloatVector& x, const FloatVector& z) {
    Interval r=mul_ivl(x[0],z[0]);
    for(uint i=0; i!=x.size(); ++i) { r=hull(mul_ivl(x[i],z[i]),r); }
    return r;
}

Float widths(const Box& bx) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width());
    }
    return res;
}

Float maximum_scaled_width(const Box& bx, const FloatVector& sf) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res=max(bx[i].width()/sf[i],res);
    }
    return res;
}

Float average_scaled_width(const Box& bx, const FloatVector& sf) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width()/sf[i]);
    }
    return res/bx.size();
}

Float average_width(const Box& bx) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        if(bx[i].lower()>bx[i].upper()) { return -inf; }
        res+=bx[i].width();
    }
    return res/bx.size();
}

} // namespace

void SubdivisionPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const Box& constraint_bounds, Int depth) const
{
    return Ariadne::subdivision_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

void AffinePaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const Box& constraint_bounds, Int depth) const
{
    return Ariadne::affine_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

void ConstraintPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const Box& constraint_bounds, Int depth) const
{
    return Ariadne::constraint_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}

void OptimalConstraintPaver::
adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunctionInterface& space_function,
                           const ValidatedVectorFunctionInterface& constraint_function, const Box& constraint_bounds, Int depth) const
{
    return Ariadne::optimal_constraint_adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
}


namespace {

using Ariadne::verbosity;

void subdivision_adjoin_outer_approximation_recursion(PavingInterface& paving, const Box& subdomain, const ValidatedVectorFunction& function,
                                                      const List<ValidatedConstraint>& constraints, const int depth, const FloatVector& errors)
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;

    for(List<ValidatedConstraint>::const_iterator iter=constraints.begin();
        iter!=constraints.end(); ++iter)
    {
        Interval constraint_range=iter->function().evaluate(subdomain);
        if(constraint_range.lower() > iter->bounds().upper() || constraint_range.upper() < iter->bounds().lower() ) { return; }
    }

    Box range=evaluate(function,subdomain);
    bool small=true;
    for(uint i=0; i!=range.size(); ++i) {
        if(range[i].radius()>errors[i]*RELATIVE_SMALLNESS) {
            small=false;
            break;
        }
    }

    if(small) {
        paving.adjoin_outer_approximation(range,depth);
    } else {
        Box subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=Ariadne::split(subdomain);
        subdivision_adjoin_outer_approximation_recursion(paving,subdomain1,function,constraints,depth,errors);
        subdivision_adjoin_outer_approximation_recursion(paving,subdomain2,function,constraints,depth,errors);
    }
}



static uint COUNT_TESTS=0u;

// Adjoin an over-approximation to the solution of $f(dom)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
void procedure_constraint_adjoin_outer_approximation_recursion(
        PavingInterface& paving, const Box& domain, const ValidatedVectorFunction& f,
        const ValidatedVectorFunction& g, const Box& codomain, const GridCell& cell, int max_dpth, uint splt, const List<ValidatedProcedure>& procedures)
{
    const uint m=domain.size();
    const uint nf=f.result_size();
    const uint ng=g.result_size();

    const Box& cell_box=cell.box();
    const FloatVector scalings=paving.grid().lengths();

    Box bbox = f(domain);

    Float domwdth = average_width(domain);
    Float bbxwdth = average_scaled_width(bbox,paving.grid().lengths());
    Float clwdth = average_scaled_width(cell_box,paving.grid().lengths());

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
    const Box& old_domain=domain;
    Box new_domain=old_domain;
    Float olddomwdth = average_width(domain);
    Float newdomwdth = olddomwdth;

    static const double ACCEPTABLE_REDUCTION_FACTOR = 0.75;


    // Box reduction steps
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
    newdomwdth=average_width(new_domain);
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
        newdomwdth=average_width(new_domain);
        ARIADNE_LOG(6,"     domwdth="<<newdomwdth<<" dom="<<new_domain<<"\n");
    } while( !new_domain.empty() && (newdomwdth < ACCEPTABLE_REDUCTION_FACTOR * olddomwdth) );

    ARIADNE_LOG(6,"new_domain="<<new_domain);


    domwdth = average_scaled_width(new_domain,FloatVector(new_domain.size(),1.0));
    bbox=f(new_domain);
    bbxwdth=average_scaled_width(bbox,paving.grid().lengths());
    if(bbox.disjoint(cell_box) || disjoint(g(new_domain),codomain)) {
        ARIADNE_LOG(4,"  Proved disjointness using image of new domain\n");
        return;
    }

    ARIADNE_LOG(4,"                 domwdth="<<domwdth<<" bbxwdth="<<bbxwdth<<" clwdth="<<clwdth<<" dom="<<new_domain<<" bbox="<<bbox<<" cell="<<cell.box()<<"\n");

    // Decide whether to split cell or split domain by comparing size of
    // bounding box with the cell and splitting the larger.
    // It seems that a more efficient algorithm results if the domain
    // is only split if the bounding box is much larger, so we preferentiably
    // split the cell unless the bounding box is 4 times as large
    Float bbxmaxwdth = maximum_scaled_width(bbox,scalings);
    Float clmaxwdth = maximum_scaled_width(cell_box,scalings);

    if( (bbxmaxwdth > 4.0*clmaxwdth) || (cell.tree_depth()>=max_dpth && (bbxmaxwdth > clmaxwdth)) ) {
        Pair<uint,double> lipsch = lipschitz_index_and_error(f,new_domain);
        ARIADNE_LOG(4,"  Splitting domain on coordinate "<<lipsch.first<<"\n");
        Pair<Box,Box> sd=new_domain.split(lipsch.first);
        procedure_constraint_adjoin_outer_approximation_recursion(paving, sd.first, f, g, codomain, cell, max_dpth, splt+1, procedures);
        procedure_constraint_adjoin_outer_approximation_recursion(paving, sd.second, f, g, codomain, cell, max_dpth, splt+1, procedures);
    } else if(cell.tree_depth()>=max_dpth) {
        ARIADNE_LOG(4,"  Adjoining cell "<<cell_box<<"\n");
        paving.adjoin(cell);
    } else {
        ARIADNE_LOG(4,"  Splitting cell "<<cell_box<<"\n");
        Pair<GridCell,GridCell> sb = cell.split();
        procedure_constraint_adjoin_outer_approximation_recursion(paving,new_domain,f,g,codomain,sb.first, max_dpth, splt, procedures);
        procedure_constraint_adjoin_outer_approximation_recursion(paving,new_domain,f,g,codomain,sb.second, max_dpth, splt, procedures);
    }


}



void hotstarted_constraint_adjoin_outer_approximation_recursion(
    PavingInterface& r, const Box& d, const ValidatedVectorFunction& f,
    const ValidatedVectorFunction& g, const Box& c, const GridCell& b, Point x, Point y, int e)
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

    Float t;
    FloatVector z(x.size());

    if(r.superset(b)) {
        ARIADNE_LOG(2,"  Cell already in set\n");
        return;
    }

    Box bx=join(static_cast<const Box&>(b.box()),static_cast<const Box&>(c));

    ARIADNE_LOG(2,"  fg(d)="<<fg(d)<<", bx="<<bx<<"\n");
    if(disjoint(fg(d),bx)) {
        ARIADNE_LOG(2,"  Proved disjointness using direct evaluation\n");
        return;
    }


    // Relax x away from boundary
    optimiser.compute_tz(d,fg,bx,y,t,z);
    ARIADNE_LOG(2,"  z0="<<z<<", t0="<<t<<"\n");
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4," t="<<t);
        //optimiser.linearised_feasibility_step(d,fg,bx,x,y,z,t);
        try {
            optimiser.feasibility_step(d,fg,bx,x,y,z,t);
        }
        catch(NearBoundaryOfFeasibleDomainException e) {
            break;
        }
        catch(std::runtime_error e) {
            ARIADNE_ERROR(""<<e.what()<<"\n");
            break;
        }
        ARIADNE_LOG(6,", x="<<x<<", y="<<y<<", z="<<z<<"\n");
        ARIADNE_LOG(6,"  x.z="<<emulrng(x,z)<<"\n");
        if(t>0) { break; }
        if(emulrng(x,z).upper()<XZMIN) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");
    ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<"\n");

    if(!(t<=1e10)) {
        ARIADNE_WARN("feasibility failed\n");
        char c; cin >> c;
        t=0.0;
        y=midpoint(d);
        x=FloatVector(x.size(),1.0/x.size());
    }
        x = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());

    //assert(t>=-1000);

    if(t<TERR) {

        // Probably disjoint, so try to prove this
        Box nd=d;
        const Box& domain=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        TrivialSweeper sweeper;
        EffectiveScalarFunction zero_function=EffectiveScalarFunction::zero(m);
        EffectiveVectorFunction identity_function=EffectiveVectorFunction::identity(m);
        ScalarTaylorFunction txg(domain,zero_function,sweeper);
        Interval cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            txg = txg - (Interval(x[j])-Interval(x[n+j]))*ScalarTaylorFunction(domain,ValidatedScalarFunction(fg[j]),sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            txg = txg - (Interval(x[2*n+i])-Interval(x[2*n+m+i]))*ScalarTaylorFunction(domain,ValidatedScalarFunction(identity_function[i]),sweeper);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        txg = Interval(cnst) + txg;

        ARIADNE_LOG(6,"    txg="<<txg<<"\n");

        ValidatedConstraint constraint=(txg>=0.0);

        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,txg,Interval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(nd.empty()) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(uint i=0; i!=m; ++i) {
            solver.box_reduce(nd,txg,Interval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(nd.empty()) { ARIADNE_LOG(2,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        solver.hull_reduce(nd,txg,Interval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(nd.empty()) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }
    }

    if(t<=0.0 && Box(f(d)).radius()>b.box().radius()) {
        ARIADNE_LOG(2,"  Splitting domain\n");
        Pair<Box,Box> sd=d.split();
        x = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        y=midpoint(sd.first);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.first, f,g, c, b, x, y, e);
        y = midpoint(sd.second);
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


void hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(PavingInterface& r, const Box& d, const VectorTaylorFunction& fg, const Box& c, const GridCell& b, Point& x, Point& y, int e)
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

    Float t;
    Point z(x.size());

    if(r.superset(b)) {
        return;
    }

    Box bx=join(static_cast<const Box&>(b.box()),static_cast<const Box&>(c));

    optimiser.compute_tz(d,fg,bx,y,t,z);
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4," t="<<t);
        optimiser.linearised_feasibility_step(d,fg,bx,x,y,z,t);
        if(t>0) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        Box nd=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ScalarTaylorFunction xg=ScalarTaylorFunction::zero(d,sweeper);
        Interval cnst=0.0;
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
        solver.hull_reduce(nd,xg,Interval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(nd.empty()) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(uint i=0; i!=m; ++i) {
            solver.box_reduce(nd,xg,Interval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(nd.empty()) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        //Pair<Box,Box> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=split(d);
        Point nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point ny = midpoint(sd.first);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.first, fg, c, b, nx, ny, e);
        nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        Point sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.first,sx,sy,e);
        sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.second,sx,sy,e);
    }


}



void subdivision_adjoin_outer_approximation(PavingInterface& paving,
                                            const Box& subdomain,
                                            const ValidatedVectorFunction& function,
                                            const ValidatedVectorFunction& constraint_functions,
                                            const Box& constraint_bounds,
                                            int depth)
{
    List<ValidatedConstraint> constraints;
    for(uint i=0; i!=constraint_functions.result_size(); ++i) {
        constraints.append(ValidatedConstraint(constraint_bounds[i].lower(),constraint_functions[i],constraint_bounds[i].upper()));
    }

    FloatVector errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }

    ::subdivision_adjoin_outer_approximation_recursion(paving,subdomain,function,constraints,depth,errors);
}

void affine_adjoin_outer_approximation(PavingInterface& paving,
                                       const Box& subdomain,
                                       const ValidatedVectorFunction& function,
                                       const ValidatedVectorFunction& constraints,
                                       const Box& bounds,
                                       int depth)
{
    ARIADNE_NOT_IMPLEMENTED;
}

void
constraint_adjoin_outer_approximation(PavingInterface& p, const Box& d, const ValidatedVectorFunction& f,
                                      const ValidatedVectorFunction& g, const Box& c, int e)
{
    ARIADNE_ASSERT(p.dimension()==f.result_size());

    GridCell b=GridCell::smallest_enclosing_primary_cell(f(d),p.grid());
    Box r=g(d)+Box(g.result_size(),Interval(-1,1));
    Box rc=intersection(r,c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    ::hotstarted_constraint_adjoin_outer_approximation_recursion(p,d,f,g,rc,b,x,y,e);
}

void
procedure_constraint_adjoin_outer_approximation(PavingInterface& p, const Box& d, const ValidatedVectorFunction& f,
                                                const ValidatedVectorFunction& g, const Box& c, int e)
{
    GridCell b=p.smallest_enclosing_primary_cell(f(d));

    List<ValidatedProcedure> procedures;
    procedures.reserve(f.result_size()+g.result_size());
    for(uint i=0; i!=f.result_size(); ++i) { procedures.append(make_procedure(f[i])); }
    for(uint i=0; i!=g.result_size(); ++i) { procedures.append(make_procedure(g[i])); }

    Ariadne::procedure_constraint_adjoin_outer_approximation_recursion(p,d,f,g,c,b,e*p.dimension(),0, procedures);
    //std::cerr<<"Computing outer approximation considered a total of "<<COUNT_TESTS<<" domains/cells\n";
    //std::cerr<<"Measure of paving is "<<p.measure()<<"\n";

    if(dynamic_cast<GridTreeSet*>(&p)) { dynamic_cast<GridTreeSet&>(p).recombine(); }
}

void optimal_constraint_adjoin_outer_approximation(PavingInterface& p, const Box& d, const ValidatedVectorFunction& f,
                                                   const ValidatedVectorFunction& g, const Box& c, int e)
{
    GridCell b=GridCell::smallest_enclosing_primary_cell(g(d),p.grid());
    Box rc=intersection(g(d)+Box(g.result_size(),Interval(-1,1)),c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

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


