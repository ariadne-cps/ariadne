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

#include "config.h"

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


void subdivision_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                            const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void affine_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                       const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void constraint_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                           const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void optimal_constraint_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                                   const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

Matrix<Float> nonlinearities_zeroth_order(const IntervalVectorFunction& f, const IntervalVector& dom);
Pair<uint,double> nonlinearity_index_and_error(const IntervalVectorFunction& function, const IntervalVector domain);
Pair<uint,double> lipschitz_index_and_error(const IntervalVectorFunction& function, const IntervalVector& domain);

Matrix<Float> nonlinearities_zeroth_order(const VectorTaylorFunction& f, const IntervalVector& dom)
{
    const uint m=f.result_size();
    const uint n=f.argument_size();
    VectorTaylorFunction g=restrict(f,dom);

    Matrix<Float> nonlinearities=Matrix<Float>::zero(m,n);
    MultiIndex a;
    for(uint i=0; i!=m; ++i) {
        const IntervalTaylorModel& tm=g.model(i);
        for(IntervalTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
            a=iter->key();
            if(a.degree()>1) {
                for(uint j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=mag(iter->data()); }
                }
            }
        }
    }

    return nonlinearities;
}

Matrix<Float> nonlinearities_first_order(const IntervalVectorFunctionInterface& f, const IntervalVector& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const uint m=f.result_size();
    const uint n=f.argument_size();
    Vector<IntervalDifferential> ivl_dx=IntervalDifferential::constants(m,n, 1, dom);
    MultiIndex a(n);
    for(uint i=0; i!=n; ++i) {
        Float sf=dom[i].radius();
        ++a[i];
        ivl_dx[i].expansion().append(a,Interval(sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<IntervalDifferential> df=f.evaluate(ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float> nonlinearities=Matrix<Float>::zero(m,n);
    for(uint i=0; i!=m; ++i) {
        const IntervalDifferential& d=df[i];
        for(IntervalDifferential::const_iterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==1) {
                for(uint j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=radius(iter->data()); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}

Matrix<Float> nonlinearities_second_order(const IntervalVectorFunctionInterface& f, const IntervalVector& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const uint m=f.result_size();
    const uint n=f.argument_size();
    Vector<IntervalDifferential> ivl_dx=IntervalDifferential::constants(m,n, 2, dom);
    MultiIndex a(n);
    for(uint i=0; i!=n; ++i) {
        Float sf=dom[i].radius();
        ++a[i];
        ivl_dx[i].expansion().append(a,Interval(sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<IntervalDifferential> df=f.evaluate(ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float> nonlinearities=Matrix<Float>::zero(m,n);
    for(uint i=0; i!=m; ++i) {
        const IntervalDifferential& d=df[i];
        for(IntervalDifferential::const_iterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==2) {
                for(uint j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=mag(iter->data()); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}

Pair<uint,double> nonlinearity_index_and_error(const VectorTaylorFunction& function, const IntervalVector domain) {
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(function,domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    uint imax=nonlinearities.row_size();
    uint jmax_in_row_imax=nonlinearities.column_size();
    Float max_row_sum=0.0;
    for(uint i=0; i!=nonlinearities.row_size(); ++i) {
        uint jmax=nonlinearities.column_size();
        Float row_sum=0.0;
        Float max_mag_j_in_i=0.0;
        for(uint j=0; j!=nonlinearities.column_size(); ++j) {
            row_sum+=mag(nonlinearities[i][j]);
            if(mag(nonlinearities[i][j])>max_mag_j_in_i) {
                jmax=j;
                max_mag_j_in_i=mag(nonlinearities[i][j]);
            }
        }
        if(row_sum>max_row_sum) {
            imax=i;
            max_row_sum=row_sum;
            jmax_in_row_imax=jmax;
        }
    }

    return make_pair(jmax_in_row_imax,numeric_cast<double>(max_row_sum));
}


RealBox::RealBox(const IntervalVector& bx) : _ary(bx.size()) {
    for(uint i=0; i!=bx.size(); ++i) {
        this->_ary[i]=RealInterval(ExactFloat(bx[i].lower()),ExactFloat(bx[i].upper()));
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


Float widths(const IntervalVector& bx) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width());
    }
    return res;
}

Float maximum_scaled_width(const IntervalVector& bx, const FloatVector& sf) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res=max(bx[i].width()/sf[i],res);
    }
    return res;
}

Float average_scaled_width(const IntervalVector& bx, const FloatVector& sf) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width()/sf[i]);
    }
    return res/bx.size();
}

Float average_width(const IntervalVector& bx) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        if(bx[i].lower()>bx[i].upper()) { return -inf<Float>(); }
        res+=bx[i].width();
    }
    return res/bx.size();
}



RealConstraintSet::RealConstraintSet(const RealVectorFunction& fn, const RealBox& codom)
    : _function(fn), _codomain(codom)
{
    ARIADNE_ASSERT(codom.size()==fn.result_size());
}

RealConstraintSet::RealConstraintSet(const List<RealNonlinearConstraint>& c)
    : _function(), _codomain()
{
    if(c.size()!=0) {
        _function=RealVectorFunction(c.size(),c[0].function().argument_size());
        _codomain=RealBox(c.size());
    }
    for(uint i=0; i!=c.size(); ++i) {
        _function[i]=c[i].function();
        _codomain[i]=RealInterval(c[i].lower_bound(),c[i].upper_bound());
    }
}

RealConstraintSet*
RealConstraintSet::clone() const
{
    return new RealConstraintSet(*this);
}


uint
RealConstraintSet::dimension() const
{
    return this->_function.argument_size();
}


tribool
RealConstraintSet::separated(const Box& bx) const
{
    Box codomain=over_approximation(this->codomain());
    return ConstrainedImageSet(bx,this->function()).separated(codomain);
}

tribool
RealConstraintSet::overlaps(const Box& bx) const
{
    Box codomain=under_approximation(this->codomain());
    return ConstrainedImageSet(bx,this->function()).overlaps(codomain);
}

tribool
RealConstraintSet::covers(const Box& bx) const
{
    Box codomain=under_approximation(this->codomain());
    return Box(this->function().evaluate(bx)).inside(codomain);
}


std::ostream&
RealConstraintSet::write(std::ostream& os) const
{
    return os << "RealConstraintSet( function=" << this->function() << ", codomain=" << this->codomain() << ")";
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
RealBoundedConstraintSet::separated(const Box& bx) const
{
    Box domain=over_approximation(this->domain());
    if(Ariadne::disjoint(domain,bx)) { return true; }
    Box codomain=over_approximation(this->codomain());
    return ConstrainedImageSet(Ariadne::intersection(bx,domain),this->function()).separated(codomain);
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
RealBoundedConstraintSet::draw(CanvasInterface& c, const Projection2d& p) const
{
    return RealConstrainedImageSet(*this).draw(c,p);
}






Interval emulrng(const FloatVector& x, const FloatVector& z) {
    Interval r=mul_ivl(x[0],z[0]);
    for(uint i=0; i!=x.size(); ++i) { r=hull(mul_ivl(x[i],z[i]),r); }
    return r;
}






RealConstrainedImageSet::RealConstrainedImageSet(const RealBoundedConstraintSet& set)
    : _domain(over_approximation(set.domain())), _function(RealVectorFunction::identity(set.dimension()))
{
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        this->new_parameter_constraint(RealNonlinearConstraint(ExactFloat(set.codomain()[i].lower()),set.function()[i],ExactFloat(set.codomain()[i].upper())));
    }
}


Box RealConstrainedImageSet::bounding_box() const
{
    return this->_function(over_approximation(this->_domain));
}

AffineSet
RealConstrainedImageSet::affine_approximation() const
{
    static const Float inf=std::numeric_limits<double>::infinity();

    const Vector<Interval> D=approximation(this->domain());
    Vector<Float> m=midpoint(D);
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


tribool RealConstrainedImageSet::satisfies(const RealNonlinearConstraint& nc) const
{
    if( subset(nc.function().evaluate(this->bounding_box()),nc.bounds()) ) {
        return true;
    }

    ConstraintSolver solver;
    const RealBox& domain=this->_domain;
    List<RealNonlinearConstraint> all_constraints=this->_constraints;
    RealScalarFunction composed_function = compose(nc.function(),this->_function);
    const Real& lower_bound = nc.lower_bound();
    const Real& upper_bound = nc.upper_bound();

    Tribool result;
    if(upper_bound<+infinity) {
        all_constraints.append( composed_function >= upper_bound );
        result=solver.feasible(over_approximation(domain),all_constraints).first;
        all_constraints.pop_back();
        if(definitely(result)) { return false; }
    }
    if(lower_bound>-infinity) {
        all_constraints.append(composed_function <= lower_bound);
        result = result || solver.feasible(over_approximation(domain),all_constraints).first;
    }
    return !result;
}


tribool RealConstrainedImageSet::separated(const Box& bx) const
{
    ConstraintSolver solver;
    const RealBox& domain=this->_domain;

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
        all_constraints.append(RealNonlinearConstraint(ExactFloat(bx[i].lower()),this->_function[i],ExactFloat(bx[i].upper())));
    }
    all_constraints.append(this->_constraints);
    Pair<Tribool,FloatVector> result=ConstraintSolver().feasible(over_approximation(domain),all_constraints);
    return !result.first;
}


tribool RealConstrainedImageSet::overlaps(const Box& bx) const
{
    return !this->separated(bx);
}


void
RealConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    //paving.adjoin_outer_approximation(*this,depth);
    this->constraint_adjoin_outer_approximation_to(paving,depth);
}


Pair<RealConstrainedImageSet,RealConstrainedImageSet>
RealConstrainedImageSet::split() const
{
    uint k=this->number_of_parameters();
    Float rmax=0.0;
    for(uint j=0; j!=this->number_of_parameters(); ++j) {
        if(Float(this->domain()[j].radius())>rmax) {
            k=j;
            rmax=this->domain()[j].radius();
        }
    }
    return this->split(k);
}

Pair<RealConstrainedImageSet,RealConstrainedImageSet>
RealConstrainedImageSet::split(uint j) const
{
    RealInterval interval = this->domain()[j];
    Real midpoint = interval.midpoint();
    Pair<RealBox,RealBox> subdomains(this->domain(),this->domain());
    subdomains.first[j]=RealInterval(interval.lower(),midpoint);
    subdomains.second[j]=RealInterval(midpoint,interval.upper());
    return make_pair(RealConstrainedImageSet(subdomains.first,this->_function,this->_constraints),
                     RealConstrainedImageSet(subdomains.second,this->_function,this->_constraints));
}


RealConstrainedImageSet image(const RealBoundedConstraintSet& set, const RealVectorFunction& function) {
    ARIADNE_ASSERT(set.dimension()==function.argument_size());
    RealConstrainedImageSet result(set.domain(),function);
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        result.new_parameter_constraint(set.constraint(i));
    }
    return result;
}
















namespace {

} // namespace

IntervalProcedure make_procedure(const IntervalScalarFunctionInterface& f) {
    Formula<Interval> e=f.evaluate(Formula<Interval>::identity(f.argument_size()));
    return Procedure<Interval>(e);
}




Matrix<Float> nonlinearities_zeroth_order(const VectorTaylorFunction& f, const IntervalVector& dom);


Matrix<Float> nonlinearities_zeroth_order(const IntervalVectorFunction& f, const IntervalVector& dom)
{
    ARIADNE_ASSERT(dynamic_cast<const VectorTaylorFunction*>(f.raw_pointer()));
    return nonlinearities_zeroth_order(dynamic_cast<const VectorTaylorFunction&>(*f.raw_pointer()),dom);
}

/*
Matrix<Float> nonlinearities_first_order(const IntervalVectorFunctionInterface& f, const IntervalVector& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const uint m=f.result_size();
    const uint n=f.argument_size();
    Vector<IntervalDifferential> ivl_dx=IntervalDifferential::constants(m,n, 1, dom);
    MultiIndex a(n);
    for(uint i=0; i!=n; ++i) {
        Float sf=dom[i].radius();
        ++a[i];
        ivl_dx[i].expansion().append(a,Interval(sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<IntervalDifferential> df=f.evaluate(ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float> nonlinearities=Matrix<Float>::zero(m,n);
    for(uint i=0; i!=m; ++i) {
        const IntervalDifferential& d=df[i];
        for(IntervalDifferential::const_iterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==1) {
                for(uint j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=radius(iter->data()); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}

Matrix<Float> nonlinearities_second_order(const IntervalVectorFunctionInterface& f, const IntervalVector& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const uint m=f.result_size();
    const uint n=f.argument_size();
    Vector<IntervalDifferential> ivl_dx=IntervalDifferential::constants(m,n, 2, dom);
    MultiIndex a(n);
    for(uint i=0; i!=n; ++i) {
        Float sf=dom[i].radius();
        ++a[i];
        ivl_dx[i].expansion().append(a,Interval(sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<IntervalDifferential> df=f.evaluate(ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float> nonlinearities=Matrix<Float>::zero(m,n);
    for(uint i=0; i!=m; ++i) {
        const IntervalDifferential& d=df[i];
        for(IntervalDifferential::const_iterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==2) {
                for(uint j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=mag(iter->data()); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}
*/

Pair<uint,double> lipschitz_index_and_error(const IntervalVectorFunction& function, const IntervalVector& domain)
{
    Matrix<Interval> jacobian=function.jacobian(domain);

    // Compute the column of the matrix which has the norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    uint jmax=domain.size();
    Float max_column_norm=0.0;
    for(uint j=0; j!=domain.size(); ++j) {
        Float column_norm=0.0;
        for(uint i=0; i!=function.result_size(); ++i) {
            column_norm+=mag(jacobian[i][j]);
        }
        column_norm *= domain[j].radius();
        if(column_norm>max_column_norm) {
            max_column_norm=column_norm;
            jmax=j;
        }
    }
    return make_pair(jmax,numeric_cast<double>(max_column_norm));
}

Pair<uint,double> nonlinearity_index_and_error(const IntervalVectorFunction& function, const IntervalVector domain)
{
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(function,domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    uint imax=nonlinearities.row_size();
    uint jmax_in_row_imax=nonlinearities.column_size();
    Float max_row_sum=0.0;
    for(uint i=0; i!=nonlinearities.row_size(); ++i) {
        uint jmax=nonlinearities.column_size();
        Float row_sum=0.0;
        Float max_mag_j_in_i=0.0;
        for(uint j=0; j!=nonlinearities.column_size(); ++j) {
            row_sum+=mag(nonlinearities[i][j]);
            if(mag(nonlinearities[i][j])>max_mag_j_in_i) {
                jmax=j;
                max_mag_j_in_i=mag(nonlinearities[i][j]);
            }
        }
        if(row_sum>max_row_sum) {
            imax=i;
            max_row_sum=row_sum;
            jmax_in_row_imax=jmax;
        }
    }

    return make_pair(jmax_in_row_imax,numeric_cast<double>(max_row_sum));
}





namespace {

void subdivision_adjoin_outer_approximation_recursion(PavingInterface& paving, const IntervalVector& subdomain, const IntervalVectorFunction& function,
                                                      const List<IntervalNonlinearConstraint>& constraints, const int depth, const FloatVector& errors)
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
        subdivision_adjoin_outer_approximation_recursion(paving,subdomain1,function,constraints,depth,errors);
        subdivision_adjoin_outer_approximation_recursion(paving,subdomain2,function,constraints,depth,errors);
    }
}



static uint COUNT_TESTS=0u;

// Adjoin an over-approximation to the solution of $f(dom)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
void procedure_constraint_adjoin_outer_approximation_recursion(
        PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& f,
        const IntervalVectorFunction& g, const Box& codomain, const GridCell& cell, int max_dpth, uint splt, const List<IntervalProcedure>& procedures)
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
    PavingInterface& r, const Box& d, const IntervalVectorFunction& f,
    const IntervalVectorFunction& g, const IntervalVector& c, const GridCell& b, Point x, Point y, int e)
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

    if(r.superset(b)) {
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


void hotstarted_optimal_constraint_adjoin_outer_approximation_recursion(PavingInterface& r, const IntervalVector& d, const VectorTaylorFunction& fg, const Box& c, const GridCell& b, Point& x, Point& y, int e)
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

    if(r.superset(b)) {
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


} // namespace


void subdivision_adjoin_outer_approximation(PavingInterface& paving,
                                            const IntervalVector& subdomain,
                                            const IntervalVectorFunction& function,
                                            const IntervalVectorFunction& constraint_functions,
                                            const IntervalVector& constraint_bounds,
                                            int depth)
{
    List<IntervalNonlinearConstraint> constraints;
    for(uint i=0; i!=constraint_functions.result_size(); ++i) {
        constraints.append(IntervalNonlinearConstraint(constraint_bounds[i].lower(),constraint_functions[i],constraint_bounds[i].upper()));
    }

    FloatVector errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }

    ::subdivision_adjoin_outer_approximation_recursion(paving,subdomain,function,constraints,depth,errors);
}

void affine_adjoin_outer_approximation(PavingInterface& paving,
                                       const IntervalVector& subdomain,
                                       const IntervalVectorFunction& function,
                                       const IntervalVectorFunction& constraints,
                                       const IntervalVector& bounds,
                                       int depth)
{
    ARIADNE_NOT_IMPLEMENTED;
}

void
constraint_adjoin_outer_approximation(PavingInterface& p, const IntervalVector& d, const IntervalVectorFunction& f,
                                      const IntervalVectorFunction& g, const IntervalVector& c, int e)
{
    ARIADNE_ASSERT(p.dimension()==f.result_size());

    GridCell b=GridCell::smallest_enclosing_primary_cell(f(d),p.grid());
    IntervalVector r=g(d)+IntervalVector(g.result_size(),Interval(-1,1));
    IntervalVector rc=intersection(r,c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    ::hotstarted_constraint_adjoin_outer_approximation_recursion(p,d,f,g,rc,b,x,y,e);
}

void
procedure_constraint_adjoin_outer_approximation(PavingInterface& p, const IntervalVector& d, const IntervalVectorFunction& f,
                                                const IntervalVectorFunction& g, const IntervalVector& c, int e)
{
    GridCell b=p.smallest_enclosing_primary_cell(f(d));

    List<IntervalProcedure> procedures;
    procedures.reserve(f.result_size()+g.result_size());
    for(uint i=0; i!=f.result_size(); ++i) { procedures.append(make_procedure(f[i])); }
    for(uint i=0; i!=g.result_size(); ++i) { procedures.append(make_procedure(g[i])); }

    Ariadne::procedure_constraint_adjoin_outer_approximation_recursion(p,d,f,g,c,b,e*p.dimension(),0, procedures);
    //std::cerr<<"Computing outer approximation considered a total of "<<COUNT_TESTS<<" domains/cells\n";
    //std::cerr<<"Measure of paving is "<<p.measure()<<"\n";

    if(dynamic_cast<GridTreeSet*>(&p)) { dynamic_cast<GridTreeSet&>(p).recombine(); }
}

void optimal_constraint_adjoin_outer_approximation(PavingInterface& p, const IntervalVector& d, const IntervalVectorFunction& f,
                                                   const IntervalVectorFunction& g, const IntervalVector& c, int e)
{
    GridCell b=GridCell::smallest_enclosing_primary_cell(g(d),p.grid());
    Box rc=intersection(g(d)+IntervalVector(g.result_size(),Interval(-1,1)),c);

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    std::cerr<<"Here\n";
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



void RealConstrainedImageSet::
subdivision_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    ARIADNE_ASSERT(paving.dimension()==this->dimension());
    const Box domain=over_approximation(this->domain());
    const RealVectorFunction& function=this->function();
    RealVectorFunction constraints(this->number_of_constraints(),domain.size());
    Box bounds(this->number_of_constraints());

    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        constraints.set(i,this->_constraints[i].function());
        bounds[i]=this->_constraints[i].bounds();
    }

    Ariadne::subdivision_adjoin_outer_approximation(paving,domain,function,constraints,bounds,depth);
}



void RealConstrainedImageSet::
constraint_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    ARIADNE_ASSERT(paving.dimension()==this->dimension());
    const Box domain=over_approximation(this->domain());
    const RealVectorFunction& function=this->function();
    RealVectorFunction constraints(this->number_of_constraints(),domain.size());
    Box bounds(this->number_of_constraints());

    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        constraints.set(i,this->_constraints[i].function());
        bounds[i]=this->_constraints[i].bounds();
    }

    Ariadne::constraint_adjoin_outer_approximation(paving,domain,function,constraints,bounds,depth);
}

void draw(CanvasInterface& cnvs, const Projection2d& proj, const RealConstrainedImageSet& set, uint depth)
{
    if( depth==0) {
        set.affine_approximation().draw(cnvs,proj);
    } else {
        Pair<RealConstrainedImageSet,RealConstrainedImageSet> split=set.split();
        draw(cnvs,proj,split.first,depth-1u);
        draw(cnvs,proj,split.second,depth-1u);
    }
}

void
RealConstrainedImageSet::draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    static const uint DEPTH = 0;
    Ariadne::draw(cnvs,proj,*this,DEPTH);
}



std::ostream&
RealConstrainedImageSet::write(std::ostream& os) const
{
    return os << "RealConstrainedImageSet( domain=" << this->_domain
              << ", function=" << this->_function << ", constraints=" << this->_constraints << " )";
}



} // namespace Ariadne

#include "procedure.h"
#include <include/container.h>

namespace Ariadne {

typedef tribool Tribool;
typedef unsigned int Nat;
typedef std::ostream OutputStream;

template<class SF> struct FunctionTraits;
template<class X> struct FunctionTraits< ScalarFunction<X> > { typedef VectorFunction<X> VectorFunctionType; };
template<> struct FunctionTraits< ScalarTaylorFunction > { typedef VectorTaylorFunction VectorFunctionType; };

template<class SF> class TemplatedConstraintSet;
template<class SF> class TemplatedConstrainedImageSet;



const IntervalNonlinearConstraint ConstrainedImageSet::constraint(uint i) const
{
    ARIADNE_ASSERT(i<this->_negative_constraints.size()+_zero_constraints.size());
    if(i<this->_negative_constraints.size()) {
        return this->_negative_constraints[i] <= 0.0;
    } else {
        return this->_zero_constraints[i-this->_negative_constraints.size()] == 0.0;
    }
}

const List<IntervalNonlinearConstraint> ConstrainedImageSet::constraints() const
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

IntervalVectorFunction ConstrainedImageSet::constraint_function() const
{
    IntervalVectorFunction result(this->number_of_constraints(),this->number_of_parameters());
    for(uint i=0; i!=this->_negative_constraints.size(); ++i) {
        result[i]=this->_negative_constraints[i];
    }
    for(uint i=0; i!=this->_zero_constraints.size(); ++i) {
        result[i+this->_negative_constraints.size()]=this->_zero_constraints[i];
    }
    return result;
}

IntervalVector ConstrainedImageSet::constraint_bounds() const
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
ConstrainedImageSet::bounding_box() const
{
    return this->_function(this->_reduced_domain);
}

AffineSet
ConstrainedImageSet::affine_over_approximation() const
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

AffineSet ConstrainedImageSet::affine_approximation() const
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


Pair<ConstrainedImageSet,ConstrainedImageSet> ConstrainedImageSet::split(uint j) const
{
    Pair<Box,Box> subdomains = Ariadne::split(this->_domain,j);
    subdomains.first=intersection(subdomains.first,this->_reduced_domain);
    subdomains.second=intersection(subdomains.second,this->_reduced_domain);

    Pair<ConstrainedImageSet,ConstrainedImageSet> result(
        ConstrainedImageSet(subdomains.first,this->_function),
        ConstrainedImageSet(subdomains.second,this->_function));

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


void
ConstrainedImageSet::reduce()
{
    ConstraintSolver solver;
    solver.reduce(this->_reduced_domain, this->constraint_function(), this->constraint_bounds());
}

tribool ConstrainedImageSet::empty() const
{
    const_cast<ConstrainedImageSet*>(this)->reduce();
    return this->_reduced_domain.empty();
}

tribool ConstrainedImageSet::inside(const Box& bx) const
{
    return Ariadne::inside(this->bounding_box(),bx);
}

tribool ConstrainedImageSet::separated(const Box& bx) const
{
    Box subdomain = this->_reduced_domain;
    IntervalVectorFunction function = join(this->_function,this->constraint_function());
    IntervalVector codomain = join(bx,this->constraint_bounds());
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return subdomain.empty() || indeterminate;
}

tribool ConstrainedImageSet::overlaps(const Box& bx) const
{
    Box subdomain = this->_reduced_domain;
    IntervalVectorFunction function = join(this->_function,this->constraint_function());
    IntervalVector codomain = join(bx,this->constraint_bounds());
    NonlinearInteriorPointOptimiser optimiser;
    return optimiser.feasible(subdomain,function,codomain);
}

void ConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    const IntervalVector subdomain=this->_reduced_domain;
    const IntervalVectorFunction function = this->function();
    const IntervalVectorFunction constraint_function = this->constraint_function();
    const IntervalVector constraint_bounds = this->constraint_bounds();

    switch(DISCRETISATION_METHOD) {
        case SUBDIVISION_DISCRETISE:
            Ariadne::subdivision_adjoin_outer_approximation(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        case AFFINE_DISCRETISE:
            Ariadne::affine_adjoin_outer_approximation(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        case CONSTRAINT_DISCRETISE:
            Ariadne::constraint_adjoin_outer_approximation(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown discretisation method\n");
    }

    if(dynamic_cast<GridTreeSet*>(&paving)) {
        dynamic_cast<GridTreeSet&>(paving).recombine();
    }
}



tribool ConstrainedImageSet::satisfies(const IntervalNonlinearConstraint& nc) const
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


void draw(CanvasInterface& cnvs, const Projection2d& proj, const ConstrainedImageSet& set, uint depth)
{
    if( depth==0) {
        set.affine_approximation().draw(cnvs,proj);
    } else {
        Pair<ConstrainedImageSet,ConstrainedImageSet> split=set.split();
        draw(cnvs,proj,split.first,depth-1u);
        draw(cnvs,proj,split.second,depth-1u);
    }
}

void
ConstrainedImageSet::draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    static const uint DEPTH = 0;
    Ariadne::draw(cnvs,proj,*this,DEPTH);
}



std::ostream& ConstrainedImageSet::write(std::ostream& os) const
{
    return os << "ConstrainedImageSet( domain=" << this->domain() << ", function="<< this->function() << ", constraints=" << this->constraints() << " )";
}

std::ostream& operator<<(std::ostream& os, const ConstrainedImageSet& set) {
    return set.write(os);
}







} // namespace Ariadne;
