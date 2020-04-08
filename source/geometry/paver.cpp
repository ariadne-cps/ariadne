/***************************************************************************
 *            geometry/paver.cpp
 *
 *  Copyright  2011-20  Pieter Collins
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

#include "../geometry/paver.hpp"

#include "../utility/macros.hpp"
#include "../output/logging.hpp"
#include "../function/polynomial.hpp"
#include "../function/function.hpp"
#include "../function/taylor_function.hpp"
#include "../function/procedure.hpp"
#include "../solvers/nonlinear_programming.hpp"
#include "../solvers/constraint_solver.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/affine_set.hpp"
#include "../geometry/paving_interface.hpp"
#include "../geometry/grid_paving.hpp"
#include "../geometry/affine_set.hpp"
#include "../algebra/algebra.hpp"

namespace Ariadne {

Pair<Nat,FloatDP> nonlinearity_index_and_error(const ValidatedVectorMultivariateFunction& function, const ExactBoxType& domain);
Pair<Nat,FloatDP> lipschitz_index_and_error(const ValidatedVectorMultivariateFunction& function, const ExactBoxType& domain);
inline Pair<Nat,FloatDP> lipschitz_index_and_error(const ValidatedVectorMultivariateFunction& function, const UpperBoxType& domain) {
    return lipschitz_index_and_error(function,cast_exact_box(domain));
}

namespace {

UpperIntervalType emulrng(const FloatDPValueVector& x, const FloatDPValueVector& z) {
    UpperIntervalType r=make_interval(mul(x[0],z[0]));
    for(Nat i=0; i!=x.size(); ++i) { r=hull(mul(x[i],z[i]),r); }
    return r;
}


PositiveFloatDPUpperBound average_width(const UpperBoxType& bx) {
    PositiveFloatDPUpperBound res(0u,double_precision);
    for(Nat i=0; i!=bx.size(); ++i) {
        if(definitely(bx[i].lower()>bx[i].upper())) { return cast_positive(-infty); }
        res+=bx[i].width();
    }
    return res/bx.size();
}

PositiveFloatDPUpperBound maximum_scaled_width(const UpperBoxType& bx, const Vector<PositiveFloatDPValue>& sf) {
    PositiveFloatDPUpperBound res(0u,double_precision);
    for(Nat i=0; i!=bx.size(); ++i) {
        res=max(bx[i].width()/sf[i],res);
    }
    return res;
}

PositiveFloatDPUpperBound average_scaled_width(const UpperBoxType& bx, const Vector<PositiveFloatDPValue>& sf) {
    PositiveFloatDPUpperBound res(0u,double_precision);
    for(Nat i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width()/sf[i]);
    }
    return res/bx.size();
}


FloatDP average_scaled_width(const UpperBoxType& bx, const Vector<FloatDP>& sf) {
    return average_scaled_width(bx,reinterpret_cast<Vector<PositiveFloatDPValue>const&>(sf)).raw();
}

} // namespace


OutputStream& AffinePaver::_write(OutputStream& os) const { return os << "AffinePaver()"; }
OutputStream& SubdivisionPaver::_write(OutputStream& os) const { return os << "SubdivisionPaver()"; }
OutputStream& ReducePaver::_write(OutputStream& os) const { return os << "ReducePaver()"; }
OutputStream& ConstraintPaver::_write(OutputStream& os) const { return os << "ConstraintPaver()"; }
OutputStream& OptimalConstraintPaver::_write(OutputStream& os) const { return os << "OptimalConstraintPaver()"; }

Void SubdivisionPaver::adjoin_outer_approximation(PavingInterface& paving, const ValidatedConstrainedImageSet& set, Nat fineness) const
{
    Vector<FloatDPValue> max_errors(paving.dimension());
    for(Nat i=0; i!=max_errors.size(); ++i) {
        max_errors[i]=shft(static_cast<FloatDPValue>(paving.grid().lengths()[i]),-static_cast<int>(fineness));
    }

    this->adjoin_outer_approximation_recursion(paving,set,fineness,max_errors);
}

Void SubdivisionPaver::adjoin_outer_approximation_recursion(PavingInterface& paving, ValidatedConstrainedImageSet const& set, Nat fineness, const Vector<FloatDPValue>& max_errors) const
{
    // How small an over-approximating box needs to be relative to the cell size
    static const ExactDouble RELATIVE_SMALLNESS=0.5_x;

    const ExactBoxType& subdomain = set.reduced_domain();
    const ValidatedVectorMultivariateFunction& function = set.function();
    const List<ValidatedConstraint>& constraints = set.constraints();

    for(List<ValidatedConstraint>::ConstIterator iter=constraints.begin();
        iter!=constraints.end(); ++iter)
    {
        UpperIntervalType constraint_range=apply(iter->function(),subdomain);
        if( definitely(constraint_range.lower() > iter->bounds().upper() || constraint_range.upper() < iter->bounds().lower() ) ) { return; }
    }

    UpperBoxType range=apply(function,subdomain);
    Bool small=true;
    for(Nat i=0; i!=range.size(); ++i) {
        if(possibly(hlf(range[i].width())>max_errors[i]*RELATIVE_SMALLNESS)) {
            small=false;
            break;
        }
    }

    if(small) {
        paving.adjoin_outer_approximation(cast_exact_box(range),fineness);
    } else {
        Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> subsets=set.split();
        this->adjoin_outer_approximation_recursion(paving,subsets.first,fineness,max_errors);
        this->adjoin_outer_approximation_recursion(paving,subsets.second,fineness,max_errors);
    }
}



Void AffinePaver::adjoin_outer_approximation(PavingInterface& paving,
                                             const ValidatedConstrainedImageSet& set,
                                             Nat fineness) const
{
    const ExactBoxType& subdomain = set.reduced_domain();
    const ValidatedVectorMultivariateFunction& function = set.function();
    const List<ValidatedConstraint>& constraints = set.constraints();

    // Bound the maximum number of splittings allowed to draw a particular set.
    // Note that this gives rise to possibly 2^MAX_DEPTH split sets!!
    static const Nat MAXIMUM_DEPTH = 16;

    // The basic approximation error when plotting with accuracy=0
    static const double BASIC_ERROR = 0.0625;

    const double max_error=BASIC_ERROR/(1<<fineness);

    ARIADNE_DEBUG_ASSERT(function.domain()==set.constraint_function().domain());
    ValidatedVectorMultivariateFunction fg=join(function,set.constraint_function());

    List<ExactBoxType> subdomains;
    List<ExactBoxType> unsplitdomains;
    List<ExactBoxType> splitdomains;
    unsplitdomains.append(subdomain);
    ExactBoxType splitdomain1,splitdomain2;
    for(Nat i=0; i!=MAXIMUM_DEPTH; ++i) {
        //std::cerr<<"i="<<i<<"\nsubdomains="<<subdomains<<"\nunsplitdomains="<<unsplitdomains<<"\n\n";
        for(Nat n=0; n!=unsplitdomains.size(); ++n) {
            Nat k; FloatDP err;
            make_lpair(k,err)=nonlinearity_index_and_error(fg,unsplitdomains[n]);
            //std::cerr<<"  domain="<<unsplitdomains[n]<<" k="<<k<<" err="<<err<<" max_err="<<max_error<<"\n";
            if(k==subdomain.size() || err < max_error) {
                subdomains.append(unsplitdomains[n]);
            } else {
                make_lpair(splitdomain1,splitdomain2)=unsplitdomains[n].split(k);
                splitdomains.append(splitdomain1);
                splitdomains.append(splitdomain2);
            }
        }
        unsplitdomains.swap(splitdomains);
        splitdomains.clear();
        if(unsplitdomains.empty()) { break; }
    }
    subdomains.concatenate(unsplitdomains);
    if(!unsplitdomains.empty()) {
        ARIADNE_WARN("Cannot obtain desired accuracy in outer approximation without excessive splitting.");
    }

    for(Nat n=0; n!=subdomains.size(); ++n) {
        ValidatedConstrainedImageSet subset(subdomains[n],function,constraints);
        subset.affine_over_approximation().adjoin_outer_approximation_to(paving,fineness);
    }
}

Nat verbosity;

namespace {

using Ariadne::verbosity;

template<class Y, class PR> decltype(auto) make_raw_float(Y const& y, PR pr) { return RawFloat<PR>(y,pr); }

inline Bool strictly_smaller_by_factor(PositiveFloatDPUpperBound x1, PositiveFloatDPUpperBound x2, ExactDouble sf) {
    return x1.raw() < mul(down,x2.raw(),make_raw_float(sf,dp)); }
inline Bool strictly_smaller(PositiveFloatDPUpperBound x1, PositiveFloatDPUpperBound x2) {
    return x1.raw() < x2.raw(); }


// Adjoin an over-approximation to the solution of $f(dom)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
Void procedure_constraint_adjoin_outer_approximation_recursion(
        PavingInterface& paving, const ExactBoxType& domain, const ValidatedVectorMultivariateFunction& f,
        const ValidatedVectorMultivariateFunction& g, const ExactBoxType& codomain, const GridCell& cell, Int max_dpth, Nat splt, const List<ValidatedProcedure>& procedures)
{
    const Nat m=domain.size();
    const Nat nf=f.result_size();
    const Nat ng=g.result_size();

    const ExactBoxType& cell_box=cell.box();
    const Vector<PositiveFloatDPValue> scalings=Vector<PositiveFloatDPValue>(paving.grid().lengths());

    UpperBoxType bbox = apply(f,domain);

    FloatDP domwdth = average_width(domain).raw();
    FloatDP bbxwdth = average_scaled_width(bbox,paving.grid().lengths()).raw();
    FloatDP clwdth = average_scaled_width(cell_box,paving.grid().lengths()).raw();

    ARIADNE_LOG(2,"\nconstraint_adjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"   splt="<<splt<<" dpth="<<cell.depth()<<" max_dpth="<<max_dpth<<"\n");
    ARIADNE_LOG(2,"     domwdth="<<domwdth<<" bbxwdth="<<bbxwdth<<" clwdth="<<clwdth<<" dom="<<domain<<" bbox="<<bbox<<" cell="<<cell.box()<<"\n");

    ConstraintSolver constraint_solver;

    if(paving.superset(cell)) {
        ARIADNE_LOG(4,"  Cell is already a subset of paving\n");
        return;
    }

    // Try to prove disjointness
    const ExactBoxType& old_domain=domain;
    ExactBoxType exact_new_domain=old_domain;
    UpperBoxType& new_domain = reinterpret_cast<UpperBoxType&>(exact_new_domain);
    PositiveFloatDPUpperBound olddomwdth = average_width(domain);
    PositiveFloatDPUpperBound newdomwdth = olddomwdth;

    static const ExactDouble ACCEPTABLE_REDUCTION_FACTOR = 0.75_x;


    // ExactBoxType reduction steps
    for(Nat i=0; i!=nf; ++i) {
        for(Nat j=0; j!=m; ++j) {
            constraint_solver.box_reduce(new_domain,f[i],cell_box[i],j);
            if(definitely(new_domain.is_empty())) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
    }
    for(Nat i=0; i!=ng; ++i) {
        for(Nat j=0; j!=m; ++j) {
            constraint_solver.box_reduce(new_domain,g[i],codomain[i],j);
            if(definitely(new_domain.is_empty())) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
    }
    newdomwdth=average_width(new_domain);
    ARIADNE_LOG(6,"     domwdth="<<newdomwdth<<" olddomwdth="<<olddomwdth<<" dom="<<new_domain<<" box reduce\n");

    // Hull reduction steps
    do {
        olddomwdth=newdomwdth;
        for(Nat i=0; i!=nf; ++i) {
            constraint_solver.hull_reduce(new_domain,procedures[i],cell_box[i]);
            if(definitely(new_domain.is_empty())) { ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n"); return; }
            //constraint_solver.hull_reduce(new_domain,f[i],cell_box[i]);
        }
        for(Nat i=0; i!=ng; ++i) {
            constraint_solver.hull_reduce(new_domain,procedures[nf+i],codomain[i]);
            if(definitely(new_domain.is_empty())) { ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n"); return; }
            //constraint_solver.hull_reduce(new_domain,g[i],codomain[i]);
        }
        newdomwdth=average_width(new_domain);
        ARIADNE_LOG(6,"     domwdth="<<newdomwdth<<" dom="<<new_domain<<"\n");
    } while( !definitely(new_domain.is_empty()) && strictly_smaller_by_factor(newdomwdth , olddomwdth, ACCEPTABLE_REDUCTION_FACTOR) );

    ARIADNE_LOG(6,"new_domain="<<new_domain);


    domwdth = average_scaled_width(new_domain,RawFloatDPVector(new_domain.size(),1.0));
    bbox=apply(f,new_domain);
    bbxwdth=average_scaled_width(bbox,paving.grid().lengths());
    if(definitely(bbox.disjoint(cell_box)) || definitely(codomain.disjoint(apply(g,new_domain)))) {
        ARIADNE_LOG(4,"  Proved disjointness using image of new domain\n");
        return;
    }

    ARIADNE_LOG(4,"                 domwdth="<<domwdth<<" bbxwdth="<<bbxwdth<<" clwdth="<<clwdth<<" dom="<<new_domain<<" bbox="<<bbox<<" cell="<<cell.box()<<"\n");

    // Decide whether to split cell or split domain by comparing size of
    // bounding box with the cell and splitting the larger.
    // It seems that a more efficient algorithm results if the domain
    // is only split if the bounding box is much larger, so we preferentiably
    // split the cell unless the bounding box is 4 times as large
    PositiveFloatDPUpperBound bbxmaxwdth = maximum_scaled_width(bbox,scalings);
    PositiveFloatDPUpperBound clmaxwdth = maximum_scaled_width(cell_box,scalings);
    ExactDouble RELATIVE_SPLITTING_SIZE = 4.0_x;

    if( !strictly_smaller_by_factor(bbxmaxwdth, clmaxwdth, RELATIVE_SPLITTING_SIZE) || (cell.depth()>=max_dpth && strictly_smaller(clmaxwdth, bbxmaxwdth)) ) {
        Pair<Nat,FloatDP> lipsch = lipschitz_index_and_error(f,new_domain);
        ARIADNE_LOG(4,"  Splitting domain on coordinate "<<lipsch.first<<"\n");
        Pair<ExactBoxType,ExactBoxType> sd=exact_new_domain.split(lipsch.first);
        procedure_constraint_adjoin_outer_approximation_recursion(paving, sd.first, f, g, codomain, cell, max_dpth, splt+1, procedures);
        procedure_constraint_adjoin_outer_approximation_recursion(paving, sd.second, f, g, codomain, cell, max_dpth, splt+1, procedures);
    } else if(cell.depth()>=max_dpth) {
        ARIADNE_LOG(4,"  Adjoining cell "<<cell_box<<"\n");
        paving.adjoin(cell);
    } else {
        ARIADNE_LOG(4,"  Splitting cell "<<cell_box<<"\n");
        Pair<GridCell,GridCell> sb = cell.split();
        procedure_constraint_adjoin_outer_approximation_recursion(paving,cast_exact_box(new_domain),f,g,codomain,sb.first, max_dpth, splt, procedures);
        procedure_constraint_adjoin_outer_approximation_recursion(paving,cast_exact_box(new_domain),f,g,codomain,sb.second, max_dpth, splt, procedures);
    }


}



Void hotstarted_constraint_adjoin_outer_approximation_recursion(
    PavingInterface& r, const ExactBoxType& d, const ValidatedVectorMultivariateFunction& f,
    const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c, const GridCell& b, ExactPointType x, ExactPointType y, Nat e)
{
    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const FloatDPValue XSIGMA { 0.125 };
    static const FloatDPValue TERR { -1.0/((1<<e)*1024.0) };
    static const FloatDPValue XZMIN { 1.0/(1<<16) };
    static const FloatDPValue inf { Ariadne::inf };
    DoublePrecision pr;

    // Set up the classes used for constraint propagation and
    // optimisation using the Kuhn-Tucker conditions
    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;
    ValidatedVectorMultivariateFunction fg=join(f,g);

    const Nat m=fg.argument_size();
    const Nat n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.depth()<<" e="<<e<<"\n");
    ARIADNE_LOG(2,"  x0="<<x<<", y0="<<y<<"\n");

    FloatDPValuePoint z(x.size());
    FloatDPValue t;
    FloatDPApproximation one = 1.0_exact;

    Vector<FloatDPApproximation>& ax=reinterpret_cast<Vector<FloatDPApproximation>&>(x);
    Vector<FloatDPApproximation>& ay=reinterpret_cast<Vector<FloatDPApproximation>&>(y);
    Vector<FloatDPApproximation> az=reinterpret_cast<Vector<FloatDPApproximation>&>(z);
    FloatDPApproximation at=reinterpret_cast<FloatDPApproximation&>(t);

    if(r.superset(b)) {
        ARIADNE_LOG(2,"  Cell already in set\n");
        return;
    }

    ExactBoxType bx=product(static_cast<const ExactBoxType&>(b.box()),static_cast<const ExactBoxType&>(c));

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
        catch(const NearBoundaryOfFeasibleDomainException& exc) {
            break;
        }
        catch(const std::runtime_error& err) {
            ARIADNE_FAIL_MSG(""<<err.what()<<"\n");
            break;
        }
        ARIADNE_LOG(6,", x="<<ax<<", y="<<ay<<", z="<<az<<"\n");
        ARIADNE_LOG(6,"  x.z="<<emulrng(x,z)<<"\n");
        if(t>0) { break; }
        if(definitely(emulrng(x,z).upper()<XZMIN)) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");
    ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<"\n");

    if(!(t<inf)) {
        ARIADNE_WARN("feasibility failed\n");
        char ch; std::cin >> ch;
        at=0;
        ay=midpoint(d);
        ax=FloatDPApproximationVector(x.size(),one/x.size());
    }
    ax = FloatDPApproximation(1-XSIGMA)*ax + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size());

    //assert(t>=-1000);

    if(t<TERR) {

        // Probably disjoint, so try to prove this
        UpperBoxType nd=d;
        const ExactBoxType& domain=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        TrivialSweeper<FloatDP> sweeper{dp};
        EffectiveScalarMultivariateFunction zero_function=EffectiveScalarMultivariateFunction::zero(m);
        EffectiveVectorMultivariateFunction identity_function=EffectiveVectorMultivariateFunction::identity(m);
        ValidatedScalarMultivariateTaylorFunctionModelDP txg(domain,zero_function,sweeper);
        FloatDPBounds cnst = {0,pr};
        for(Nat j=0; j!=n; ++j) {
            txg = txg - (FloatDPBounds(x[j])-FloatDPBounds(x[n+j]))*ValidatedScalarMultivariateTaylorFunctionModelDP(domain,ValidatedScalarMultivariateFunction(fg[j]),sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(Nat i=0; i!=m; ++i) {
            txg = txg - (FloatDPBounds(x[2*n+i])-FloatDPBounds(x[2*n+m+i]))*ValidatedScalarMultivariateTaylorFunctionModelDP(domain,ValidatedScalarMultivariateFunction(identity_function[i]),sweeper);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        txg = FloatDPBounds(cnst) + txg;

        ARIADNE_LOG(6,"    txg="<<txg<<"\n");

        ValidatedConstraint constraint=(txg>=0.0_exact);

        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,txg,ExactIntervalType(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.is_empty())) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(Nat i=0; i!=m; ++i) {
            solver.box_reduce(nd,txg,ExactIntervalType(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(definitely(nd.is_empty())) { ARIADNE_LOG(2,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        solver.hull_reduce(nd,txg,ExactIntervalType(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.is_empty())) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }
    }

    if(decide(t<=0.0_exact) && decide(UpperBoxType(apply(f,d)).radius()>b.box().radius()) ) {
        ARIADNE_LOG(2,"  Splitting domain\n");
        Pair<ExactBoxType,ExactBoxType> sd=d.split();
        ax = FloatDPApproximation(1-XSIGMA)*ax + Vector<FloatDPApproximation>(ax.size(),XSIGMA/x.size());
        ay=midpoint(sd.first);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.first, f,g, c, b, x, y, e);
        ay = midpoint(sd.second);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r, sd.second, f,g, c, b, x, y, e);
        return;
    }

    if(t>0.0_exact) {
        ARIADNE_LOG(2," Intersection point: parameter="<<y<<"\n");
    }

    if(b.depth()>=Int(e*b.dimension())) {
        ARIADNE_LOG(2,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(2,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        hotstarted_constraint_adjoin_outer_approximation_recursion(r,d,f,g,c,sb.first,x,y,e);
        hotstarted_constraint_adjoin_outer_approximation_recursion(r,d,f,g,c,sb.second,x,y,e);
    }
}


Void hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(PavingInterface& r, const ExactBoxType& d, const ValidatedVectorMultivariateTaylorFunctionModelDP& fg, const ExactBoxType& c, const GridCell& b, ExactPointType& x, ExactPointType& y, Nat e)
{
    auto properties = fg.properties();
    auto pr = properties.precision();

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const FloatDPValue XSIGMA = {TwoExp(-3),pr};
    static const FloatDPValue  TERR = {TwoExp(-10),pr};
    static const FloatDPValue inf { Ariadne::inf };

    const Nat m=fg.argument_size();
    const Nat n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.depth()<<" e="<<e<<"\n");

    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;

    FloatDPValue t{pr};
    FloatDPValuePoint z(x.size());

    FloatDPApproximationVector& ax=reinterpret_cast<FloatDPApproximationVector&>(x);
    FloatDPApproximationVector& ay=reinterpret_cast<FloatDPApproximationVector&>(y);
    FloatDPApproximationVector& az=reinterpret_cast<FloatDPApproximationVector&>(z);
    FloatDPApproximation& at=reinterpret_cast<FloatDPApproximation&>(t);

    if(r.superset(b)) {
        return;
    }

    ExactBoxType bx=product(static_cast<const ExactBoxType&>(b.box()),static_cast<const ExactBoxType&>(c));

    optimiser.compute_tz(d,fg,bx,ay,at,az);
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG(4," t="<<t);
        optimiser.linearised_feasibility_step(d,fg,bx,ax,ay,az,at);
        if(t>0) { break; }
    }
    ARIADNE_LOG(4,"\n  t="<<t<<"\n  y="<<y<<"\n    x="<<x<<"\n    z="<<z<<"\n");

    if(t<TERR) {
        // Probably disjoint, so try to prove this
        UpperBoxType nd=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        ValidatedScalarMultivariateTaylorFunctionModelDP xg=ValidatedScalarMultivariateTaylorFunctionModelDP::zero(d,properties);
        FloatDPBounds cnst = {0,pr};
        for(Nat j=0; j!=n; ++j) {
            xg = xg - (x[j]-x[n+j])*ValidatedScalarMultivariateTaylorFunctionModelDP(d,fg[j],properties);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(Nat i=0; i!=m; ++i) {
            xg = xg - (x[2*n+i]-x[2*n+m+i])*ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(d,i,properties);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        xg = (cnst) + xg;

        ARIADNE_LOG(4,"    xg="<<xg<<"\n");


        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,xg,ExactIntervalType(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(definitely(nd.is_empty())) {
            ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(Nat i=0; i!=m; ++i) {
            solver.box_reduce(nd,xg,ExactIntervalType(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(definitely(nd.is_empty())) { ARIADNE_LOG(4,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        //Pair<ExactBoxType,ExactBoxType> sd=solver.split(List<EffectiveConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<ExactBoxType,ExactBoxType> sd=split(d);
        ExactPointType nx = cast_exact((1-XSIGMA)*ax + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size()));
        ExactPointType ny = midpoint(sd.first);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.first, fg, c, b, nx, ny, e);
        nx = cast_exact((1-XSIGMA)*x + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size()));
        ny = midpoint(sd.second);
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.depth()>=Int(e*b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        ExactPointType sx = cast_exact((1-XSIGMA)*x + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size()));
        ExactPointType sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.first,sx,sy,e);
        sx = cast_exact(FloatDPApproximation(1-XSIGMA)*x + Vector<FloatDPApproximation>(x.size(),XSIGMA/x.size()));
        sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(r,d,fg,c,sb.second,sx,sy,e);
    }


}





Void
constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBoxType& d, const ValidatedVectorMultivariateFunction& f,
                                      const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c, Nat e)
{
    ARIADNE_ASSERT(p.dimension()==f.result_size());

    GridCell b=GridCell::smallest_enclosing_primary_cell(cast_exact_box(apply(f,d)),p.grid());
    ExactBoxType r=cast_exact_box(widen(apply(g,d),1));
    ExactBoxType rc=intersection(r,c);

    Point<FloatDPValue> y=midpoint(d);
    const Nat l=(d.size()+f.result_size()+g.result_size())*2;
    Point<FloatDPValue> x(l); for(Nat k=0; k!=l; ++k) { x[k]=FloatDPValue(1.0/l); }

    Ariadne::hotstarted_constraint_adjoin_outer_approximation_recursion(p,d,f,g,rc,b,x,y,e);
}

Void
procedure_constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBoxType& d, const ValidatedVectorMultivariateFunction& f,
                                                const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c, Nat e)
{
    GridCell b=p.smallest_enclosing_primary_cell(cast_exact_box(apply(f,d)));

    List<ValidatedProcedure> procedures;
    procedures.reserve(f.result_size()+g.result_size());
    for(Nat i=0; i!=f.result_size(); ++i) { procedures.append(make_procedure(f[i])); }
    for(Nat i=0; i!=g.result_size(); ++i) { procedures.append(make_procedure(g[i])); }

    Ariadne::procedure_constraint_adjoin_outer_approximation_recursion(p,d,f,g,c,b,e*p.dimension(),0, procedures);
    //std::cerr<<"Computing outer approximation considered a total of "<<COUNT_TESTS<<" domains/cells\n";
    //std::cerr<<"Measure of paving is "<<p.measure()<<"\n";

    if(dynamic_cast<GridTreePaving*>(&p)) { dynamic_cast<GridTreePaving&>(p).recombine(); }
}

Void optimal_constraint_adjoin_outer_approximation(PavingInterface& p, const ExactBoxType& d, const ValidatedVectorMultivariateFunction& f,
                                                   const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c, Nat e)
{
    GridCell b=GridCell::smallest_enclosing_primary_cell(cast_exact_box(apply(g,d)),p.grid());
    ExactBoxType rc=intersection(cast_exact_box(widen(apply(g,d),1)),c);

    ExactPointType y=midpoint(d);
    const Nat l=(d.size()+f.result_size()+g.result_size())*2;
    ExactPointType x(l); for(Nat k=0; k!=l; ++k) { x[k]=FloatDPValue(1.0/l); }

    ValidatedVectorMultivariateTaylorFunctionModelDP fg;
    const ValidatedVectorMultivariateTaylorFunctionModelDP* tfptr;
    if( (tfptr=dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP*>(f.raw_pointer())) ) {
        const ValidatedVectorMultivariateTaylorFunctionModelDP* tgptr;
        if( ( tgptr = dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP*>(g.raw_pointer()) ) ) {
            fg=join(*tfptr,*tgptr);
        } else {
            if(g.result_size()>0) {
                fg=join(*tfptr,ValidatedVectorMultivariateTaylorFunctionModelDP(tfptr->domain(),g,tfptr->properties()));
            } else {
                fg=*tfptr;
            }
        }
    } else {
        ThresholdSweeper<FloatDP> swp(dp,1e-12);
        fg=ValidatedVectorMultivariateTaylorFunctionModelDP(d,join(f,g),swp);
    }
    Ariadne::hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(p,d,fg,rc,b,x,y,e);
}

} // namespace


Void ReducePaver::adjoin_outer_approximation(PavingInterface& paving, const ValidatedConstrainedImageSet& set, Nat fineness) const {
    return procedure_constraint_adjoin_outer_approximation(paving,set.domain(),set.function(),set.constraint_function(),set.constraint_bounds(),fineness);
}

Void ConstraintPaver::adjoin_outer_approximation(PavingInterface& paving, const ValidatedConstrainedImageSet& set, Nat fineness) const {
    return constraint_adjoin_outer_approximation(paving,set.domain(),set.function(),set.constraint_function(),set.constraint_bounds(),fineness);
}

Void OptimalConstraintPaver::adjoin_outer_approximation(PavingInterface& paving, const ValidatedConstrainedImageSet& set, Nat fineness) const {
    return optimal_constraint_adjoin_outer_approximation(paving,set.reduced_domain(),set.function(),set.constraint_function(),set.constraint_bounds(),fineness);
}

} // namespace Ariadne


