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
#include "procedure.h"
#include "function_set.h"
#include "affine_set.h"
#include "grid_set.h"
#include "nonlinear_programming.h"
#include "constraint_solver.h"
#include "affine_set.h"

#include "graphics_interface.h"

namespace Ariadne {

static const uint verbosity = 0u;

//! \related TaylorConstrainedImageSet \brief The possible types of method used to draw a nonlinear set.
enum DrawingMethod { CURVE_DRAW, BOX_DRAW, AFFINE_DRAW, GRID_DRAW };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DrawingMethod DRAWING_METHOD;
//! \related TaylorConstrainedImageSet \brief The accuracy used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern unsigned int DRAWING_ACCURACY;

//! \related TaylorConstrainedImageSet \brief The possible types of method used to discretise a nonlinear set.
enum DiscretisationMethod { SUBDIVISION_DISCRETISE, AFFINE_DISCRETISE, CONSTRAINT_DISCRETISE };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to discretise a nonlinear set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DiscretisationMethod DISCRETISATION_METHOD;

DrawingMethod DRAWING_METHOD=AFFINE_DRAW;
DiscretisationMethod DISCRETISATION_METHOD=SUBDIVISION_DISCRETISE;
unsigned int DRAWING_ACCURACY=1u;

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;



RealBox::RealBox(const IntervalVector& bx) : _ary(bx.size()) {
    for(uint i=0; i!=bx.size(); ++i) {
        this->_ary[i]=RealInterval(Real(bx[i].lower()),Real(bx[i].upper()));
    }
}

Box under_approximation(const RealBox& rbx) {
    Box bx(rbx.size());
    for(uint i=0; i!=bx.size(); ++i) {
        bx[i]=under_approximation(rbx[i]);
    }
    return bx;
}

Box over_approximation(const RealBox& rbx) {
    Box bx(rbx.size());
    for(uint i=0; i!=bx.size(); ++i) {
        bx[i]=over_approximation(rbx[i]);
    }
    return bx;
}

Box approximation(const RealBox& rbx) {
    Box bx(rbx.size());
    for(uint i=0; i!=bx.size(); ++i) {
        bx[i]=approximation(rbx[i]);
    }
    return bx;
}




RealBoundedConstraintSet::RealBoundedConstraintSet(const RealBox& bx)
    : _domain(bx), _function(RealVectorFunction(0u,bx.dimension())), _codomain()
{
}

RealBoundedConstraintSet::RealBoundedConstraintSet(const RealBox& dom, const RealVectorFunction& fn, const RealBox& codom)
    : _domain(dom), _function(fn), _codomain(codom)
{
    ARIADNE_ASSERT(codom.size()==fn.result_size());
    ARIADNE_ASSERT(dom.size()==fn.argument_size());
}

RealBoundedConstraintSet::RealBoundedConstraintSet(const RealBox& dom, const List<RealNonlinearConstraint>& c)
    : _domain(dom), _function(c.size(),dom.size()), _codomain(c.size())
{
    for(uint i=0; i!=c.size(); ++i) {
        ARIADNE_ASSERT(_domain.size()==c[i].function().argument_size());
        _function[i]=c[i].function();
        _codomain[i]=RealInterval(c[i].lower_bound(),c[i].upper_bound());
    }
}

RealBoundedConstraintSet*
RealBoundedConstraintSet::clone() const
{
    return new RealBoundedConstraintSet(*this);
}


uint
RealBoundedConstraintSet::dimension() const
{
    return this->_domain.size();
}


tribool
RealBoundedConstraintSet::disjoint(const Box& bx) const
{
    Box domain=over_approximation(this->domain());
    if(Ariadne::disjoint(domain,bx)) { return true; }
    Box codomain=over_approximation(this->codomain());
    return ConstrainedImageSet(Ariadne::intersection(bx,domain),this->function()).disjoint(codomain);
}


tribool
RealBoundedConstraintSet::overlaps(const Box& bx) const
{
    if(Ariadne::disjoint(over_approximation(this->domain()),bx)) { return false; }
    Box domain=under_approximation(this->domain());
    Box codomain=under_approximation(this->codomain());
    return ConstrainedImageSet(Ariadne::intersection(bx,domain),this->function()).overlaps(codomain);
}


tribool
RealBoundedConstraintSet::covers(const Box& bx) const
{
    Box domain=under_approximation(this->domain());
    Box codomain=under_approximation(this->codomain());
    if(!Ariadne::covers(domain,bx)) { return false; }
    return Box(this->function().evaluate(bx)).inside(codomain);
}

tribool
RealBoundedConstraintSet::inside(const Box& bx) const
{
    if(Ariadne::inside(over_approximation(this->domain()),bx)) { return true; }
    return indeterminate;
}

Box
RealBoundedConstraintSet::bounding_box() const
{
    Box result=over_approximation(this->_domain);
    result.widen();
    return result;
}


std::ostream&
RealBoundedConstraintSet::write(std::ostream& os) const
{
    return os << "RealBoundedConstraintSet( domain=" << this->domain() << ", function=" << this->function() << ", codomain=" << this->codomain() << ")";
}

void
RealBoundedConstraintSet::draw(CanvasInterface& os) const
{
    return ConstrainedImageSet(BoundedConstraintSet(approximation(this->domain()),this->function(),approximation(this->codomain()))).draw(os);
}






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
    return Ariadne::empty(this->_domain);
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

ConstraintSet::ConstraintSet(const List<RealNonlinearConstraint>& c)
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



BoundedConstraintSet::BoundedConstraintSet(const Box& bx)
    : _domain(bx), _function(RealVectorFunction(0u,bx.dimension())), _codomain()
{
}

BoundedConstraintSet::BoundedConstraintSet(const Vector<Interval>& dom, const RealVectorFunction& fn, const Vector<Interval>& codom)
    : _domain(dom), _function(fn), _codomain(codom)
{
    ARIADNE_ASSERT(codom.size()==fn.result_size());
    ARIADNE_ASSERT(dom.size()==fn.argument_size());
}

BoundedConstraintSet::BoundedConstraintSet(const Vector<Interval>& dom, const List<RealNonlinearConstraint>& c)
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
    return ConstrainedImageSet(Ariadne::intersection(static_cast<const IntervalVector&>(bx),this->domain()),this->function()).disjoint(this->codomain());
}


tribool
BoundedConstraintSet::overlaps(const Box& bx) const
{
    if(Ariadne::disjoint(this->domain(),bx)) { return false; }
    return ConstrainedImageSet(Ariadne::intersection(static_cast<const IntervalVector&>(bx),this->domain()),this->function()).overlaps(this->codomain());
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
        this->new_parameter_constraint(RealNonlinearConstraint(set.function()[i],set.codomain()[i]));
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
    for(List<RealNonlinearConstraint>::const_iterator iter=this->_constraints.begin();
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


tribool ConstrainedImageSet::satisfies(const RealNonlinearConstraint& nc) const
{
    const Float infty = inf<Float>();

    if( subset(nc.function().evaluate(this->bounding_box()),nc.bounds()) ) {
        return true;
    }

    ConstraintSolver solver;
    const Box& domain=this->_domain;
    List<RealNonlinearConstraint> all_constraints=this->_constraints;
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
    List<RealNonlinearConstraint> all_constraints;
    for(uint i=0; i!=this->dimension(); ++i) {
        all_constraints.append(RealNonlinearConstraint(this->_function[i],bx[i]));
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

} // namespace

IntervalProcedure make_procedure(const IntervalScalarFunctionInterface& f) {
    Formula<Interval> e=f.evaluate(Formula<Interval>::identity(f.argument_size()));
    return Procedure<Interval>(e);
}

void subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, const IntervalVector& subdomain, const IntervalVectorFunction& function, const List<IntervalNonlinearConstraint>& constraints, const int depth, const FloatVector& errors)
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;

    for(List<IntervalNonlinearConstraint>::const_iterator iter=constraints.begin();
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




void subdivision_adjoin_outer_approximation_to(GridTreeSet& paving,
                                               const Vector<Interval>& subdomain,
                                               const IntervalVectorFunction& function,
                                               const IntervalVectorFunction constraint_functions,
                                               const IntervalVector& constraint_bounds,
                                               int depth)
{
    List<IntervalNonlinearConstraint> constraints;
    for(uint i=0; i!=constraint_functions.result_size(); ++i) {
        constraints.append(IntervalNonlinearConstraint(constraint_functions[i],constraint_bounds[i]));
    }

    FloatVector errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }

    subdivision_adjoin_outer_approximation_to(paving,subdomain,function,constraints,depth,errors);
}



void affine_adjoin_outer_approximation_to(GridTreeSet& paving, const IntervalVector& subdomain, const IntervalVectorFunction& function,                                                const IntervalVectorFunction constraint_functions, const IntervalVector& constraint_bounds, const int depth)
{
    ARIADNE_NOT_IMPLEMENTED;
}



namespace {

void hotstarted_constraint_adjoin_outer_approximation_to(GridTreeSet& r, const Box& d, const IntervalVectorFunction& f, const IntervalVectorFunction& g, const Box& c, const GridCell& b, Point x, Point y, int e)
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
    IntervalVectorFunction fg=join(f,g);

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
        const Box& domain=d;

        // Use the computed dual variables to try to make a scalar function which is negative over the entire domain.
        // This should be easier than using all constraints separately
        TrivialSweeper sweeper;
        RealScalarFunction zero_function=RealScalarFunction::zero(m);
        RealVectorFunction identity_function=RealVectorFunction::identity(m);
        ScalarTaylorFunction txg(domain,zero_function,sweeper);
        Interval cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            txg = txg - (Interval(x[j])-Interval(x[n+j]))*ScalarTaylorFunction(domain,IntervalScalarFunction(fg[j]),sweeper);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            txg = txg - (Interval(x[2*n+i])-Interval(x[2*n+m+i]))*ScalarTaylorFunction(domain,IntervalScalarFunction(identity_function[i]),sweeper);
            cnst += (d[i].upper()*x[2*n+i]-d[i].lower()*x[2*n+m+i]);
        }
        txg = Interval(cnst) + txg;

        ARIADNE_LOG(6,"    txg="<<txg<<"\n");

        IntervalNonlinearConstraint constraint=(txg>=0.0);

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
        hotstarted_constraint_adjoin_outer_approximation_to(r, sd.first, f,g, c, b, x, y, e);
        y = midpoint(sd.second);
        hotstarted_constraint_adjoin_outer_approximation_to(r, sd.second, f,g, c, b, x, y, e);
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
        hotstarted_constraint_adjoin_outer_approximation_to(r,d,f,g,c,sb.first,x,y,e);
        hotstarted_constraint_adjoin_outer_approximation_to(r,d,f,g,c,sb.second,x,y,e);
    }
}


void hotstarted_optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& r, const Box& d, const VectorTaylorFunction& fg, const Box& c, const GridCell& b, Point& x, Point& y, int e)
{
    Sweeper sweeper = fg.sweeper();

    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const double TERR = -1.0/((1<<e)*1024.0);
    static const Float inf = Ariadne::inf<Float>();

    const uint m=fg.argument_size();
    const uint n=fg.result_size();
    ARIADNE_LOG(2,"\nadjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");

    ConstraintSolver solver;
    NonlinearInteriorPointOptimiser optimiser;

    Float t;
    Point z(x.size());

    if(subset(b,r)) {
        return;
    }

    Box bx=join(static_cast<const IntervalVector&>(b.box()),static_cast<const IntervalVector&>(c));

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

        //Pair<Box,Box> sd=solver.split(List<RealNonlinearConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=d.split();
        Point nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point ny = midpoint(sd.first);
        hotstarted_optimal_constraint_adjoin_outer_approximation_to(r, sd.first, fg, c, b, nx, ny, e);
        nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        hotstarted_optimal_constraint_adjoin_outer_approximation_to(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        Point sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_to(r,d,fg,c,sb.first,sx,sy,e);
        sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        sy = y;
        hotstarted_optimal_constraint_adjoin_outer_approximation_to(r,d,fg,c,sb.second,sx,sy,e);
    }


}


} // namespace


void
constraint_adjoin_outer_approximation_to(GridTreeSet& p, const Box& d, const IntervalVectorFunction& f, const IntervalVectorFunction& g, const IntervalVector& c, int e)
{
    ARIADNE_ASSERT(p.dimension()==f.result_size());

    GridCell b=GridCell::smallest_enclosing_primary_cell(f(d),p.grid());
    IntervalVector r=g(d)+IntervalVector(g.result_size(),Interval(-1,1));
    IntervalVector rc=intersection(r,c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    ::hotstarted_constraint_adjoin_outer_approximation_to(p,d,f,g,rc,b,x,y,e);
}

void optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& p, const Box& d, const VectorTaylorFunction& f, const VectorTaylorFunction& g, const IntervalVector& c, int e)
{

    GridCell b=GridCell::smallest_enclosing_primary_cell(g(d),p.grid());
    Box rc=intersection(g(d)+IntervalVector(g.result_size(),Interval(-1,1)),c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    VectorTaylorFunction fg=join(f,g);

    ::hotstarted_optimal_constraint_adjoin_outer_approximation_to(p,d,fg,rc,b,x,y,e);

}


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



void ConstrainedImageSet::
constraint_adjoin_outer_approximation_to(GridTreeSet& p, int e) const
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

    Ariadne::constraint_adjoin_outer_approximation_to(p,d,f,g,c,e);
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

#include "procedure.h"

namespace Ariadne {

typedef tribool Tribool;
typedef unsigned int Nat;
typedef std::ostream OutputStream;

template<class SF> struct FunctionTraits;
template<class X> struct FunctionTraits< ScalarFunction<X> > { typedef VectorFunction<X> VectorFunctionType; };
template<> struct FunctionTraits< ScalarTaylorFunction > { typedef VectorTaylorFunction VectorFunctionType; };

template<class SF> class TemplatedConstraintSet;
template<class SF> class TemplatedConstrainedImageSet;



const List<IntervalNonlinearConstraint> IntervalConstrainedImageSet::constraints() const
{
    List<IntervalNonlinearConstraint> result;
    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        result.push_back(this->_negative_constraints[i] <= 0.0);
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        result.push_back(this->_zero_constraints[i] == 0.0);
    }
    return result;
}

IntervalVectorFunction IntervalConstrainedImageSet::constraint_function() const
{
    IntervalVectorFunction result(this->number_of_constraints(),IntervalScalarFunction());
    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        result[i]=this->_negative_constraints[i];
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        result[i+this->_negative_constraints.size()]=this->_zero_constraints[i];
    }
    return result;
}

IntervalVector IntervalConstrainedImageSet::constraint_bounds() const
{
    IntervalVector result(this->number_of_constraints());
    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        result[i]=Interval(-inf<Float>(),0.0);
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        result[i+this->_negative_constraints.size()]=Interval(0.0);
    }
    return result;
}


Box
IntervalConstrainedImageSet::bounding_box() const
{
    return this->_function(this->_reduced_domain);
}

AffineSet
IntervalConstrainedImageSet::affine_over_approximation() const
{
    typedef List<IntervalScalarFunction>::const_iterator const_iterator;

    const uint nx=this->dimension();
    //const uint nnc=this->_negative_constraints.size();
    //const uint nzc=this->_zero_constraints.size();
    const uint np=this->number_of_parameters();

    AffineSweeper affine_sweeper;

    VectorTaylorFunction function(this->domain(),this->function(),affine_sweeper);

    // Compute the number of values with a nonzero error
    uint nerr=0;
    for(uint i=0; i!=nx; ++i) {
        if(function[i].error()>0.0) { ++nerr; }
    }

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np+nerr);
    uint ierr=0; // The index where the error bound should go
    for(uint i=0; i!=nx; ++i) {
        ScalarTaylorFunction component_function=function[i];
        h[i]=component_function.model().value();
        for(uint j=0; j!=np; ++j) {
            G[i][j]=component_function.model().gradient(j);
        }
        if(component_function.model().error()>0.0) {
            G[i][np+ierr]=component_function.model().error();
            ++ierr;
        }
    }

    AffineSet result(G,h);

    Vector<Float> a(np+nerr, 0.0);
    Float b;

    for(const_iterator iter=this->_negative_constraints.begin();
            iter!=this->_negative_constraints.end(); ++iter) {
        ScalarTaylorFunction constraint_function(this->_reduced_domain,*iter,affine_sweeper);
        b=sub_up(constraint_function.model().error(),constraint_function.model().value());
        for(uint j=0; j!=np; ++j) { a[j]=constraint_function.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=this->_zero_constraints.begin();
            iter!=this->_zero_constraints.end(); ++iter) {
        ScalarTaylorFunction constraint_function(this->_reduced_domain,*iter,affine_sweeper);
        for(uint j=0; j!=np; ++j) { a[j]=constraint_function.model().gradient(j); }
        if(constraint_function.model().error()==0.0) {
            b=-constraint_function.model().value();
            result.new_equality_constraint(a,b);
        } else {
            b=sub_up(constraint_function.model().error(),constraint_function.model().value());
            result.new_inequality_constraint(a,b);
            b=add_up(constraint_function.model().error(),constraint_function.model().value());
            result.new_inequality_constraint(-a,b);
        }
    }

    ARIADNE_LOG(2,"set="<<*this<<"\nset.affine_over_approximation()="<<result<<"\n");
    return result;
}

AffineSet IntervalConstrainedImageSet::affine_approximation() const
{
    static const Float inf=std::numeric_limits<double>::infinity();

    Vector<Float> m=midpoint(this->domain());

    const Vector<Interval> D=this->domain();
    Matrix<Float> G=this->_function.jacobian(m);
    Vector<Float> h=this->_function.evaluate(m)-G*m;
    AffineSet result(D,G,h);


    Vector<Float> a(this->number_of_parameters());
    Float b,l,u;
    List<IntervalNonlinearConstraint> constraints=this->constraints();
    for(List<IntervalNonlinearConstraint>::const_iterator iter=constraints.begin();
        iter!=constraints.end(); ++iter)
    {
        IntervalScalarFunction function=iter->function();
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


Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> IntervalConstrainedImageSet::split(uint j) const
{
    Pair<Box,Box> subdomains = Ariadne::split(this->_domain,j);
    subdomains.first=intersection(subdomains.first,this->_reduced_domain);
    subdomains.second=intersection(subdomains.second,this->_reduced_domain);

    Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> result(
        IntervalConstrainedImageSet(subdomains.first,this->_function),
        IntervalConstrainedImageSet(subdomains.second,this->_function));

    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        result.first.new_negative_parameter_constraint(Ariadne::restrict(this->_negative_constraints[i],subdomains.first));
        result.second.new_negative_parameter_constraint(Ariadne::restrict(this->_negative_constraints[i],subdomains.second));
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        result.first.new_zero_parameter_constraint(Ariadne::restrict(this->_zero_constraints[i],subdomains.first));
        result.second.new_zero_parameter_constraint(Ariadne::restrict(this->_zero_constraints[i],subdomains.second));
    }

    return result;
}

Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet>
IntervalConstrainedImageSet::split() const
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


void
IntervalConstrainedImageSet::reduce()
{
    ConstraintSolver solver;
    solver.reduce(this->_reduced_domain, this->constraint_function(), this->constraint_bounds());
}

tribool IntervalConstrainedImageSet::empty() const
{
    const_cast<IntervalConstrainedImageSet*>(this)->reduce();
    return this->_reduced_domain.empty();
}

tribool IntervalConstrainedImageSet::inside(const Box& bx) const
{
    return Ariadne::inside(this->bounding_box(),bx);
}

tribool IntervalConstrainedImageSet::disjoint(const Box& bx) const
{
    Box subdomain = this->_reduced_domain;
    IntervalVectorFunction function = join(this->_function,this->constraint_function());
    IntervalVector codomain = join(bx,this->constraint_bounds());
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return subdomain.empty() || indeterminate;
}

tribool IntervalConstrainedImageSet::overlaps(const Box& bx) const
{
    Box subdomain = this->_reduced_domain;
    IntervalVectorFunction function = join(this->_function,this->constraint_function());
    IntervalVector codomain = join(bx,this->constraint_bounds());
    NonlinearInteriorPointOptimiser optimiser;
    return optimiser.feasible(subdomain,function,codomain);
}

void IntervalConstrainedImageSet::adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    const IntervalVector subdomain=this->_reduced_domain;
    const IntervalVectorFunction function = this->function();
    const IntervalVectorFunction constraint_function = this->constraint_function();
    const IntervalVector constraint_bounds = this->constraint_bounds();

    switch(DISCRETISATION_METHOD) {
        case SUBDIVISION_DISCRETISE:
            Ariadne::subdivision_adjoin_outer_approximation_to(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        case AFFINE_DISCRETISE:
            Ariadne::affine_adjoin_outer_approximation_to(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        case CONSTRAINT_DISCRETISE:
            Ariadne::constraint_adjoin_outer_approximation_to(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown discretisation method\n");
    }

    paving.recombine();
}



tribool IntervalConstrainedImageSet::satisfies(const IntervalNonlinearConstraint& nc) const
{
    const Float infty = inf<Float>();

    if( subset(nc.function().evaluate(this->bounding_box()),nc.bounds()) ) {
        return true;
    }

    ConstraintSolver solver;
    const Box& domain=this->_domain;
    List<IntervalNonlinearConstraint> all_constraints=this->constraints();
    IntervalScalarFunction composed_function = compose(nc.function(),this->_function);
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


void draw(CanvasInterface& cnvs, const IntervalConstrainedImageSet& set, uint depth)
{
    if( depth==0) {
        set.affine_approximation().draw(cnvs);
    } else {
        Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> split=set.split();
        draw(cnvs,split.first,depth-1u);
        draw(cnvs,split.second,depth-1u);
    }
}

void
IntervalConstrainedImageSet::draw(CanvasInterface& cnvs) const
{
    static const uint DEPTH = 0;
    Ariadne::draw(cnvs,*this,DEPTH);
}



std::ostream& IntervalConstrainedImageSet::write(std::ostream& os) const
{
    return os << "IntervalConstrainedImageSet( domain=" << this->domain() << ", function="<< this->function() << ", constraints=" << this->constraints() << " )";
}

std::ostream& operator<<(std::ostream& os, const IntervalConstrainedImageSet& set) {
    return set.write(os);
}







} // namespace Ariadne;