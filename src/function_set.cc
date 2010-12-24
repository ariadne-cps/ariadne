/***************************************************************************
 *            function_set.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "boost/multi_array.hpp"

#include "macros.h"
#include "logging.h"
#include "polynomial.h"
#include "function.h"
#include "function_set.h"
#include "affine_set.h"
#include "grid_set.h"
#include "nonlinear_programming.h"
#include "constraint_solver.h"
#include "affine_set.h"

namespace Ariadne {

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

Interval emulrng(const FloatVector& x, const FloatVector& z) {
    Interval r=mul_ivl(x[0],z[0]);
    for(uint i=0; i!=x.size(); ++i) { r=hull(mul_ivl(x[i],z[i]),r); }
    return r;
}


ImageSet::ImageSet()
    : _domain(), _function()
{
}

ImageSet::ImageSet(const Vector<Interval>& dom)
    : _domain(dom),
      _function(IdentityFunction(dom.size()))
{
}


ImageSet::ImageSet(const Vector<Interval>& dom, const RealVectorFunction& fn)
    : _domain(dom), _function(fn)
{
    ARIADNE_ASSERT(dom.size()==fn.argument_size());
}


ImageSet*
ImageSet::clone() const
{
    return new ImageSet(*this);
}


uint
ImageSet::dimension() const
{
    return this->_function.result_size();
}


tribool
ImageSet::empty() const
{
    return this->_domain.empty();
}


tribool
ImageSet::disjoint(const Box& bx) const
{
    return !ConstraintSolver().feasible(this->_domain,this->_function,bx).first;
}


tribool
ImageSet::overlaps(const Box& bx) const
{
    return ConstraintSolver().feasible(this->_domain,this->_function,bx).first;
}


tribool
ImageSet::inside(const Box& bx) const
{
    return this->bounding_box().inside(bx) || indeterminate;
}


Box
ImageSet::bounding_box() const
{
    return this->_function.evaluate(this->_domain);
}


void
ImageSet::draw(CanvasInterface& canvas) const
{
    return ConstrainedImageSet(this->domain(),this->function()).draw(canvas);
}


std::ostream&
ImageSet::write(std::ostream& os) const
{
    return os << "ImageSet( domain=" << this->domain() << ", function=" << this->function() << ")";
}





ConstraintSet::ConstraintSet(const Vector<Interval>& codom, const RealVectorFunction& fn)
    : _codomain(codom), _function(fn)
{
    ARIADNE_ASSERT(codom.size()==fn.result_size());
}

ConstraintSet::ConstraintSet(const List<NonlinearConstraint>& c)
    : _codomain(c.size()), _function(c.size())
{
    uint d=0u;
    if(!c.empty()) { d=c[0].function().argument_size(); }
    for(uint i=0; i!=c.size(); ++i) {
        ARIADNE_ASSERT(c[i].function().argument_size()==d);
        _function[i]=c[i].function();
        _codomain[i]=c[i].bounds();
    }
}

ConstraintSet*
ConstraintSet::clone() const
{
    return new ConstraintSet(*this);
}


uint
ConstraintSet::dimension() const
{
    return this->_function.argument_size();
}


tribool
ConstraintSet::disjoint(const Box& bx) const
{
    return ConstrainedImageSet(bx,this->function()).disjoint(this->codomain());
}


tribool
ConstraintSet::overlaps(const Box& bx) const
{
    return ConstrainedImageSet(bx,this->function()).overlaps(this->codomain());
}


tribool
ConstraintSet::covers(const Box& bx) const
{
    return Box(this->function().evaluate(bx)).inside(this->codomain());
}


std::ostream&
ConstraintSet::write(std::ostream& os) const
{
    return os << "ConstraintSet( function=" << this->function() << ", codomain=" << this->codomain() << ")";
}

BoundedConstraintSet intersection(const ConstraintSet& set, const Box& bound) {
    ARIADNE_ASSERT(set.dimension()==bound.dimension());
    return BoundedConstraintSet(bound, set.function(), set.codomain());
}



BoundedConstraintSet::BoundedConstraintSet(const Vector<Interval>& dom, const RealVectorFunction& fn, const Vector<Interval>& codom)
    : _domain(dom), _function(fn), _codomain(codom)
{
    ARIADNE_ASSERT(codom.size()==fn.result_size());
    ARIADNE_ASSERT(dom.size()==fn.argument_size());
}

BoundedConstraintSet::BoundedConstraintSet(const Vector<Interval>& dom, const List<NonlinearConstraint>& c)
    : _domain(dom), _function(c.size(),dom.size()), _codomain(c.size())
{
    for(uint i=0; i!=c.size(); ++i) {
        ARIADNE_ASSERT(_domain.size()==c[i].function().argument_size());
        _function[i]=c[i].function();
        _codomain[i]=c[i].bounds();
    }
}

BoundedConstraintSet*
BoundedConstraintSet::clone() const
{
    return new BoundedConstraintSet(*this);
}


uint
BoundedConstraintSet::dimension() const
{
    return this->_domain.size();
}


tribool
BoundedConstraintSet::disjoint(const Box& bx) const
{
    if(Ariadne::disjoint(this->domain(),bx)) { return true; }
    return ConstrainedImageSet(Ariadne::intersection(bx,Box(this->domain())),this->function()).disjoint(this->codomain());
}


tribool
BoundedConstraintSet::overlaps(const Box& bx) const
{
    if(Ariadne::disjoint(this->domain(),bx)) { return false; }
    return ConstrainedImageSet(Ariadne::intersection(bx,Box(this->domain())),this->function()).overlaps(this->codomain());
}


tribool
BoundedConstraintSet::covers(const Box& bx) const
{
    if(!Ariadne::covers(this->domain(),bx)) { return false; }
    return Box(this->function().evaluate(bx)).inside(this->codomain());
}

tribool
BoundedConstraintSet::inside(const Box& bx) const
{
    if(Ariadne::inside(this->domain(),bx)) { return true; }
    return indeterminate;
}

Box
BoundedConstraintSet::bounding_box() const
{
    Box result=this->_domain;
    result.widen();
    return result;
}


std::ostream&
BoundedConstraintSet::write(std::ostream& os) const
{
    return os << "BoundedConstraintSet( domain=" << this->domain() << ", function=" << this->function() << ", codomain=" << this->codomain() << ")";
}

void
BoundedConstraintSet::draw(CanvasInterface& os) const
{
    return ConstrainedImageSet(*this).draw(os);
}

ConstrainedImageSet image(const BoundedConstraintSet& set, const RealVectorFunction& function) {
    ARIADNE_ASSERT(set.dimension()==function.argument_size());
    ConstrainedImageSet result(set.domain(),function);
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        result.new_parameter_constraint(set.constraint(i));
    }
    return result;
}




ConstrainedImageSet::ConstrainedImageSet(const BoundedConstraintSet& set)
    : _domain(set.domain()), _function(IdentityFunction(set.dimension()))
{
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        this->new_parameter_constraint(NonlinearConstraint(set.function()[i],set.codomain()[i]));
    }
}

ConstrainedImageSet::ConstrainedImageSet(const ImageSet& set)
    : _domain(set.domain()), _function(set.function())
{
}


Box ConstrainedImageSet::bounding_box() const
{
    return this->_function(this->_domain);
}

template<class X> Vector<X> operator*(const Matrix<X>& A, const Vector<X>& b) { return prod(A,b); }

AffineSet
ConstrainedImageSet::affine_approximation() const
{
    static const Float inf=std::numeric_limits<double>::infinity();

    Vector<Float> m=midpoint(this->domain());

    const Vector<Interval> D=this->domain();
    Matrix<Float> G=this->_function.jacobian(m);
    Vector<Float> h=this->_function.evaluate(m)-G*m;
    AffineSet result(D,G,h);


    Vector<Float> a(this->number_of_parameters());
    Float b,l,u;
    for(List<NonlinearConstraint>::const_iterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        RealScalarFunction function=iter->function();
        Interval bounds=iter->bounds();
        a=function.gradient(m);
        b=function(m)-dot(a,m);
        l=bounds.lower();
        u=bounds.upper();
        if(l==u) {
            result.new_equality_constraint(a,u-b);
        } else {
            if(u< inf) { result.new_inequality_constraint( a,u-b); }
            if(l>-inf) { result.new_inequality_constraint(-a,b-l); }
        }
    }

    return result;
}


tribool ConstrainedImageSet::satisfies(const NonlinearConstraint& nc) const
{
    const Float infty = inf<Float>();

    if( subset(nc.function().evaluate(this->bounding_box()),nc.bounds()) ) {
        return true;
    }

    ConstraintSolver solver;
    const Box& domain=this->_domain;
    List<NonlinearConstraint> all_constraints=this->_constraints;
    RealScalarFunction composed_function = compose(nc.function(),this->_function);
    const Interval& bounds = nc.bounds();

    Tribool result;
    if(bounds.upper()<+infty) {
        all_constraints.append( composed_function >= bounds.upper() );
        result=solver.feasible(domain,all_constraints).first;
        all_constraints.pop_back();
        if(definitely(result)) { return false; }
    }
    if(bounds.lower()>-infty) {
        all_constraints.append(composed_function <= bounds.lower());
        result = result || solver.feasible(domain,all_constraints).first;
    }
    return !result;
}


tribool ConstrainedImageSet::disjoint(const Box& bx) const
{
    ConstraintSolver solver;
    const Box& domain=this->_domain;

/*
    // Set up constraints as f(D)\cap C\neq\emptyset
    const Box& codomain(this->dimension()+this->number_of_constraints());
    RealVectorFunction function(codomain.size(),domain.size());
    for(uint i=0; i!=this->dimension(); ++i) {
        function.set(i,this->_function[i]);
        codomain[i]=bx[i];
    }
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        function.set(this->dimension()+i,this->constraints()[i].function());
        codomain[this->dimension()+i=this->constraints()[i].function());
    }
    return solver.feasible(domain,function,codomain).first;
*/

    // Set up constraints as f_i(D)\cap C_i\neq\emptyset
    List<NonlinearConstraint> all_constraints;
    for(uint i=0; i!=this->dimension(); ++i) {
        all_constraints.append(NonlinearConstraint(this->_function[i],bx[i]));
    }
    all_constraints.append(this->_constraints);
    Pair<Tribool,FloatVector> result=ConstraintSolver().feasible(domain,all_constraints);
    return !result.first;
}


tribool ConstrainedImageSet::overlaps(const Box& bx) const
{
    return !this->disjoint(bx);
}


void
ConstrainedImageSet::adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    //paving.adjoin_outer_approximation(*this,depth);
    this->constraint_adjoin_outer_approximation_to(paving,depth);
}


Pair<ConstrainedImageSet,ConstrainedImageSet>
ConstrainedImageSet::split() const
{
    uint k=this->number_of_parameters();
    Float rmax=0.0;
    for(uint j=0; j!=this->number_of_parameters(); ++j) {
        if(this->domain()[j].radius()>rmax) {
            k=j;
            rmax=this->domain()[j].radius();
        }
    }
    return this->split(k);
}

Pair<ConstrainedImageSet,ConstrainedImageSet>
ConstrainedImageSet::split(uint j) const
{
    Pair<Box,Box> subdomains=Box(this->domain()).split(j);
    return make_pair(ConstrainedImageSet(subdomains.first,this->_function,this->_constraints),
                     ConstrainedImageSet(subdomains.second,this->_function,this->_constraints));
}

namespace {

void subdivision_adjoin_outer_approximation_to(GridTreeSet& paving,
                                               const Vector<Interval>& subdomain, const RealVectorFunction& function, const List<NonlinearConstraint>& constraints,
                                               uint depth, const Vector<Float>& errors)
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;

    for(List<NonlinearConstraint>::const_iterator iter=constraints.begin();
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
        IntervalVector subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=Ariadne::split(subdomain);
        subdivision_adjoin_outer_approximation_to(paving,subdomain1,function,constraints,depth,errors);
        subdivision_adjoin_outer_approximation_to(paving,subdomain2,function,constraints,depth,errors);
    }
}

} // namespace


void ConstrainedImageSet::
subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    IntervalVector subdomain=this->_domain;

    FloatVector errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }

    ::subdivision_adjoin_outer_approximation_to(paving,subdomain,this->_function,this->_constraints,depth,errors);
}



namespace {

void constraint_adjoin_outer_approximation_to(GridTreeSet& r, const Box& d, const RealVectorFunction& f, const RealVectorFunction& g, const Box& c, const GridCell& b, Point x, Point y, int e)
{
    uint verbosity=0;

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const double TERR = -1.0/((1<<e)*1024.0);
    static const double XZMIN = 1.0/(1<<16);
    static const Float inf = Ariadne::inf<Float>();

    // Set up the classes used for constraint propagation and
    // optimisation using the Kuhn-Tucker conditions
    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;
    RealVectorFunction fg=join(f,g);

    const uint m=fg.argument_size();
    const uint n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");
    ARIADNE_LOG(2,"  x0="<<x<<", y0="<<y<<"\n");

    Float t;
    FloatVector z(x.size());

    if(subset(b,r)) {
        ARIADNE_LOG(2,"  Cell already in set\n");
        return;
    }

    Box bx=join(static_cast<const IntervalVector&>(b.box()),static_cast<const IntervalVector&>(c));

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

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        RealScalarFunction xg=RealScalarFunction::constant(m,0);
        Interval cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            xg = xg - (Real(x[j])-Real(x[n+j]))*fg[j];
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            xg = xg - (Real(x[2*n+i])-Real(x[2*n+m+i]))*RealScalarFunction::coordinate(m,i);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        xg = Real(cnst) + xg;

        ARIADNE_LOG(6,"    xg="<<xg<<"\n");
        ScalarTaylorFunction txg(d,xg);
        ARIADNE_LOG(6,"    txg="<<txg.polynomial()<<"\n");

        xg=RealScalarFunction(txg.polynomial());
        NonlinearConstraint constraint=(xg>=0.0);

        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        solver.hull_reduce(nd,xg,Interval(0,inf));
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");
        if(nd.empty()) {
            ARIADNE_LOG(2,"  Proved disjointness using hull reduce\n");
            return;
        }

        for(uint i=0; i!=m; ++i) {
            solver.box_reduce(nd,xg,Interval(0,inf),i);
            ARIADNE_LOG(8,"  dom="<<nd<<"\n");
            if(nd.empty()) { ARIADNE_LOG(2,"  Proved disjointness using box reduce\n"); return; }
        }
        ARIADNE_LOG(6,"  dom="<<nd<<"\n");

        solver.hull_reduce(nd,xg,Interval(0,inf));
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
        constraint_adjoin_outer_approximation_to(r, sd.first, f,g, c, b, x, y, e);
        y = midpoint(sd.second);
        constraint_adjoin_outer_approximation_to(r, sd.second, f,g, c, b, x, y, e);
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
        constraint_adjoin_outer_approximation_to(r,d,f,g,c,sb.first,x,y,e);
        constraint_adjoin_outer_approximation_to(r,d,f,g,c,sb.second,x,y,e);
    }


}

} // namespace

void ConstrainedImageSet::constraint_adjoin_outer_approximation_to(GridTreeSet& p, int e) const
{
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    const RealVectorFunction& f=this->function();
    RealVectorFunction g(this->number_of_constraints(),d.size());
    Box c(this->number_of_constraints());

    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        g.set(i,this->_constraints[i].function());
        c[i]=this->_constraints[i].bounds();
    }

    GridCell b=GridCell::smallest_enclosing_primary_cell(f(d),p.grid());
    Box r=g(d)+IntervalVector(g.result_size(),Interval(-1,1));
    c=intersection(r,c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    ::constraint_adjoin_outer_approximation_to(p,d,f,g,c,b,x,y,e);
}

void draw(CanvasInterface& cnvs, const ConstrainedImageSet& set, uint depth)
{
    if( depth==0) {
        set.affine_approximation().draw(cnvs);
    } else {
        Pair<ConstrainedImageSet,ConstrainedImageSet> split=set.split();
        draw(cnvs,split.first,depth-1u);
        draw(cnvs,split.second,depth-1u);
    }
}

void
ConstrainedImageSet::draw(CanvasInterface& cnvs) const
{
    static const uint DEPTH = 0;
    Ariadne::draw(cnvs,*this,DEPTH);
}



std::ostream&
ConstrainedImageSet::write(std::ostream& os) const
{
    return os << "ConstrainedImageSet( domain=" << this->_domain
              << ", function=" << this->_function << ", constraints=" << this->_constraints << " )";
}



} // namespace Ariadne
