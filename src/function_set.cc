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
#include "paver.h"
#include "affine_set.h"

#include "graphics_interface.h"
#include "drawer.h"

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


Matrix<Float> nonlinearities_zeroth_order(const IntervalVectorFunction& f, const Box& dom);
Pair<uint,double> nonlinearity_index_and_error(const IntervalVectorFunction& function, const Box domain);
Pair<uint,double> lipschitz_index_and_error(const IntervalVectorFunction& function, const Box& domain);

Matrix<Float> nonlinearities_zeroth_order(const VectorTaylorFunction& f, const Box& dom)
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

Matrix<Float> nonlinearities_first_order(const IntervalVectorFunctionInterface& f, const Box& dom)
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

Matrix<Float> nonlinearities_second_order(const IntervalVectorFunctionInterface& f, const Box& dom)
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

Pair<uint,double> nonlinearity_index_and_error(const VectorTaylorFunction& function, const Box domain) {
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(function,domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
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
            max_row_sum=row_sum;
            jmax_in_row_imax=jmax;
        }
    }

    return make_pair(jmax_in_row_imax,numeric_cast<double>(max_row_sum));
}


BoxSet::BoxSet(const IntervalVector& bx) : _ary(bx.size()) {
    for(uint i=0; i!=bx.size(); ++i) {
        this->_ary[i]=IntervalSet(ExactFloat(bx[i].lower()),ExactFloat(bx[i].upper()));
    }
}

Box under_approximation(const BoxSet& rbx) {
    Box bx(rbx.size());
    for(uint i=0; i!=bx.size(); ++i) {
        bx[i]=under_approximation(rbx[i]);
    }
    return bx;
}

Box over_approximation(const BoxSet& rbx) {
    Box bx(rbx.size());
    for(uint i=0; i!=bx.size(); ++i) {
        bx[i]=over_approximation(rbx[i]);
    }
    return bx;
}

Box approximation(const BoxSet& rbx) {
    Box bx(rbx.size());
    for(uint i=0; i!=bx.size(); ++i) {
        bx[i]=approximation(rbx[i]);
    }
    return bx;
}


namespace {

uint argument_size(const List<RealConstraint>& c) {
    uint as = ( c.size()==0 ? 0 : c[0].function().argument_size() );
    for(uint i=0; i!=c.size(); ++i) {
        ARIADNE_ASSERT_MSG(c[i].function().argument_size()==as,"c="<<c);
    }
    return as;
}

RealVectorFunction constraint_function(uint as, const List<RealConstraint>& c) {
    RealVectorFunction f(c.size(),as);
    for(uint i=0; i!=c.size(); ++i) {
        //f[i]=c[i].function();
        f.set(i,c[i].function());
    }
    return f;
}

BoxSet constraint_bounds(const List<RealConstraint>& c) {
    BoxSet b(c.size());
    for(uint i=0; i!=c.size(); ++i) {
        b[i]=IntervalSet(c[i].lower_bound(),c[i].upper_bound());
    }
    return b;
}

List<RealConstraint> constraints(const RealVectorFunction& f, const BoxSet& b) {
    ARIADNE_ASSERT(f.result_size()==b.size());
    List<RealConstraint> c; c.reserve(b.size());
    for(uint i=0; i!=b.size(); ++i) {
        c.append(RealConstraint(b[i].lower(),f[i],b[i].upper()));
    }
    return c;
}

} //namespace


ConstraintSet::ConstraintSet(const RealVectorFunction& f, const BoxSet& b)
    : _dimension(f.argument_size()), _constraints()
{
    this->_constraints=::constraints(f,b);
}

ConstraintSet::ConstraintSet(const List<RealConstraint>& c)
    : _dimension(argument_size(c)), _constraints(c)
{
}

RealVectorFunction const ConstraintSet::constraint_function() const
{
    return ::constraint_function(this->dimension(),this->constraints());
}

BoxSet const ConstraintSet::constraint_bounds() const
{
    return ::constraint_bounds(this->constraints());
}

ConstraintSet*
ConstraintSet::clone() const
{
    return new ConstraintSet(*this);
}


uint
ConstraintSet::dimension() const
{
    return this->_dimension;
}


tribool
ConstraintSet::separated(const Box& bx) const
{
    Box codomain=over_approximation(this->codomain());
    return IntervalConstrainedImageSet(bx,this->constraint_function()).separated(codomain) || indeterminate;
}

tribool
ConstraintSet::overlaps(const Box& bx) const
{
    Box codomain=under_approximation(this->codomain());
    return IntervalConstrainedImageSet(bx,this->constraint_function()).overlaps(codomain) || indeterminate;
}

tribool
ConstraintSet::covers(const Box& bx) const
{
    Box codomain=under_approximation(this->codomain());
    return Box(this->constraint_function().evaluate(bx)).inside(codomain) || indeterminate;
}


std::ostream&
ConstraintSet::write(std::ostream& os) const
{
    return os << "ConstraintSet( constraints=" << this->constraints() << " )";
}




BoundedConstraintSet::BoundedConstraintSet(const BoxSet& bx)
    : _domain(bx), _constraints()
{
}

BoundedConstraintSet::BoundedConstraintSet(const BoxSet& d, const RealVectorFunction& f, const BoxSet& b)
    : _domain(d), _constraints(::constraints(f,b))
{
    ARIADNE_ASSERT(b.size()==f.result_size());
    ARIADNE_ASSERT(d.size()==f.argument_size());
}

BoundedConstraintSet::BoundedConstraintSet(const BoxSet& d, const List<RealConstraint>& c)
    : _domain(d), _constraints(c)
{
}

RealVectorFunction const BoundedConstraintSet::constraint_function() const
{
    return ::constraint_function(this->dimension(),this->constraints());
}

BoxSet const BoundedConstraintSet::constraint_bounds() const
{
    return ::constraint_bounds(this->constraints());
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
BoundedConstraintSet::separated(const Box& bx) const
{
    Box domain=over_approximation(this->domain());
    if(Ariadne::disjoint(domain,bx)) { return true; }
    Box codomain=over_approximation(this->codomain());
    return IntervalConstrainedImageSet(Ariadne::intersection(bx,domain),this->constraint_function()).separated(codomain) || indeterminate;
}


tribool
BoundedConstraintSet::overlaps(const Box& bx) const
{
    if(Ariadne::disjoint(over_approximation(this->domain()),bx)) { return false; }
    Box domain=under_approximation(this->domain());
    Box codomain=under_approximation(this->codomain());
    return IntervalConstrainedImageSet(Ariadne::intersection(bx,domain),this->constraint_function()).overlaps(codomain) || indeterminate;
}


tribool
BoundedConstraintSet::covers(const Box& bx) const
{
    Box domain=under_approximation(this->domain());
    Box codomain=under_approximation(this->codomain());
    if(!Ariadne::covers(domain,bx)) { return false; }
    return Box(this->constraint_function().evaluate(bx)).inside(codomain) || indeterminate;
}

tribool
BoundedConstraintSet::inside(const Box& bx) const
{
    return Ariadne::inside(over_approximation(this->domain()),bx) || indeterminate;
}

Box
BoundedConstraintSet::bounding_box() const
{
    Box result=over_approximation(this->_domain);
    result.widen();
    return result;
}


std::ostream&
BoundedConstraintSet::write(std::ostream& os) const
{
    return os << "BoundedConstraintSet( domain=" << this->domain() << ", constraints=" << this->constraints() << ")";
}

void
BoundedConstraintSet::draw(CanvasInterface& c, const Projection2d& p) const
{
    return RealConstrainedImageSet(*this).draw(c,p);
}


BoundedConstraintSet
intersection(const ConstraintSet& cs,const BoxSet& bx)
{
    return BoundedConstraintSet(bx,cs.constraints());
}








RealConstrainedImageSet::RealConstrainedImageSet(const BoundedConstraintSet& set)
    : _domain(over_approximation(set.domain())), _function(RealVectorFunction::identity(set.dimension()))
{
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        this->new_parameter_constraint(set.constraint(i));
    }
}


const RealVectorFunction RealConstrainedImageSet::constraint_function() const
{
    RealVectorFunction result(this->number_of_constraints(),this->number_of_parameters());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=this->constraint(i).function();
    }
    return result;
}

const BoxSet RealConstrainedImageSet::constraint_bounds() const
{
    BoxSet result(this->number_of_constraints());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=IntervalSet(this->constraint(i).lower_bound(),this->constraint(i).upper_bound());
    }
    return result;
}


Box RealConstrainedImageSet::bounding_box() const
{
    return this->_function(over_approximation(this->_domain));
}


IntervalAffineConstrainedImageSet
RealConstrainedImageSet::affine_approximation() const
{
    const Vector<Interval> D=approximation(this->domain());
    Vector<Float> m=midpoint(D);
    Matrix<Float> G=this->_function.jacobian(m);
    Vector<Float> h=this->_function.evaluate(m)-G*m;
    IntervalAffineConstrainedImageSet result(D,G,h);


    Vector<Float> a(this->number_of_parameters());
    Float b,l,u;
    for(List<RealConstraint>::const_iterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        AffineModel<Interval> a=affine_model(D,iter->function());
        Interval b=iter->bounds();
        result.new_constraint(b.lower()<=a<=b.upper());
    }

    return result;
}


tribool RealConstrainedImageSet::satisfies(const RealConstraint& nc) const
{
    if( subset(nc.function().evaluate(this->bounding_box()),nc.bounds()) ) {
        return true;
    }

    ConstraintSolver solver;
    const BoxSet& domain=this->_domain;
    List<RealConstraint> all_constraints=this->_constraints;
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
    Box subdomain = over_approximation(this->_domain);
    RealVectorFunction function = join(this->function(),this->constraint_function());
    Box codomain = product(bx,Box(over_approximation(this->constraint_bounds())));
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return tribool(subdomain.empty()) || indeterminate;


}


tribool RealConstrainedImageSet::overlaps(const Box& bx) const
{
    return IntervalConstrainedImageSet(under_approximation(this->_domain),this->_function,this->_constraints).overlaps(bx);
    return indeterminate;
    //ARIADNE_NOT_IMPLEMENTED;
    Box subdomain = under_approximation(this->_domain);
    RealVectorFunction function = join(this->function(),this->constraint_function());
    Box codomain = product(bx,Box(under_approximation(this->constraint_bounds())));
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return !tribool(subdomain.empty()) || indeterminate;
}


void
RealConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    ARIADNE_ASSERT(paving.dimension()==this->dimension());
    const Box domain=over_approximation(this->domain());
    const RealVectorFunction& space_function=this->function();
    RealVectorFunction constraint_function(this->number_of_constraints(),domain.size());
    BoxSet codomain(this->number_of_constraints());

    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        constraint_function.set(i,this->_constraints[i].function());
        codomain[i]=IntervalSet(this->_constraints[i].lower_bound(),this->_constraints[i].upper_bound());
    }
    Box constraint_bounds=over_approximation(codomain);

    switch(DISCRETISATION_METHOD) {
        case SUBDIVISION_DISCRETISE:
            SubdivisionPaver().adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
            break;
        case AFFINE_DISCRETISE:
            AffinePaver().adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
            break;
        case CONSTRAINT_DISCRETISE:
            ConstraintPaver().adjoin_outer_approximation(paving,domain,space_function,constraint_function,constraint_bounds,depth);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown discretisation method\n");
    }

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
    IntervalSet interval = this->domain()[j];
    Real midpoint = interval.midpoint();
    Pair<BoxSet,BoxSet> subdomains(this->domain(),this->domain());
    subdomains.first[j]=IntervalSet(interval.lower(),midpoint);
    subdomains.second[j]=IntervalSet(midpoint,interval.upper());
    return make_pair(RealConstrainedImageSet(subdomains.first,this->_function,this->_constraints),
                     RealConstrainedImageSet(subdomains.second,this->_function,this->_constraints));
}


RealConstrainedImageSet image(const BoundedConstraintSet& set, const RealVectorFunction& function) {
    ARIADNE_ASSERT(set.dimension()==function.argument_size());
    RealConstrainedImageSet result(set.domain(),function);
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        result.new_parameter_constraint(set.constraint(i));
    }
    return result;
}



















Matrix<Float> nonlinearities_zeroth_order(const VectorTaylorFunction& f, const Box& dom);


Matrix<Float> nonlinearities_zeroth_order(const IntervalVectorFunction& f, const Box& dom)
{
    ARIADNE_ASSERT(dynamic_cast<const VectorTaylorFunction*>(f.raw_pointer()));
    return nonlinearities_zeroth_order(dynamic_cast<const VectorTaylorFunction&>(*f.raw_pointer()),dom);
}

/*
Matrix<Float> nonlinearities_first_order(const IntervalVectorFunctionInterface& f, const Box& dom)
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

Matrix<Float> nonlinearities_second_order(const IntervalVectorFunctionInterface& f, const Box& dom)
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

Pair<uint,double> lipschitz_index_and_error(const IntervalVectorFunction& function, const Box& domain)
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

Pair<uint,double> nonlinearity_index_and_error(const IntervalVectorFunction& function, const Box& domain)
{
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(function,domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
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
            max_row_sum=row_sum;
            jmax_in_row_imax=jmax;
        }
    }

    return make_pair(jmax_in_row_imax,numeric_cast<double>(max_row_sum));
}








void
RealConstrainedImageSet::draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    IntervalConstrainedImageSet(approximation(this->_domain),this->_function,this->_constraints).draw(cnvs,proj);
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
#include <include/vector.h>

namespace Ariadne {

typedef tribool Tribool;
typedef unsigned int Nat;
typedef std::ostream OutputStream;

template<class SF> struct FunctionTraits;
template<class X> struct FunctionTraits< ScalarFunction<X> > { typedef VectorFunction<X> VectorFunctionType; };
template<> struct FunctionTraits< ScalarTaylorFunction > { typedef VectorTaylorFunction VectorFunctionType; };

template<class SF> class TemplatedConstraintSet;
template<class SF> class TemplatedConstrainedImageSet;



IntervalVectorFunction IntervalConstrainedImageSet::constraint_function() const
{
    IntervalVectorFunction result(this->number_of_constraints(),this->number_of_parameters());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=this->constraint(i).function();
    }
    return result;
}

IntervalVector IntervalConstrainedImageSet::constraint_bounds() const
{
    IntervalVector result(this->number_of_constraints());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=Interval(this->constraint(i).lower_bound(),this->constraint(i).upper_bound());
    }
    return result;
}


Box
IntervalConstrainedImageSet::bounding_box() const
{
    return this->_function(this->_reduced_domain);
}


IntervalAffineConstrainedImageSet
IntervalConstrainedImageSet::affine_over_approximation() const
{
    typedef List<IntervalConstraint>::const_iterator const_iterator;

    Vector<Interval> domain = this->domain();
    Vector<IntervalAffineModel> space_models=affine_models(domain,this->function());
    List<IntervalAffineModelConstraint> constraint_models;
    constraint_models.reserve(this->number_of_constraints());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        const IntervalConstraint& constraint=this->constraint(i);
        constraint_models.append(IntervalAffineModelConstraint(constraint.lower_bound(),affine_model(domain,constraint.function()),constraint.upper_bound()));
    }

    return IntervalAffineConstrainedImageSet(domain,space_models,constraint_models);

/*
    const uint nx=this->dimension();
    //const uint nnc=this->_negative_constraints.size();
    //const uint nzc=this->_zero_constraints.size();
    const uint np=this->number_of_parameters();

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

    IntervalAffineConstrainedImageSet result(G,h);

    Vector<Float> a(np+nerr, 0.0);
    Float b;

    for(const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        ScalarTaylorFunction constraint_function(this->_reduced_domain,iter->function(),affine_sweeper);
        b=sub_up(constraint_function.model().error(),constraint_function.model().value());
        for(uint j=0; j!=np; ++j) { a[j]=constraint_function.model().gradient(j); }
        result.new_parameter_constraint(-inf,a,b);
    }

    ARIADNE_NOT_IMPLEMENTED;

    ARIADNE_LOG(2,"set="<<*this<<"\nset.affine_over_approximation()="<<result<<"\n");
    return result;
*/
}

IntervalAffineConstrainedImageSet IntervalConstrainedImageSet::affine_approximation() const
{
    typedef List<IntervalConstraint>::const_iterator const_iterator;

    Vector<Interval> domain = this->domain();
    Vector<IntervalAffineModel> space_models=affine_models(domain,this->function());
    List<IntervalAffineModelConstraint> constraint_models;
    constraint_models.reserve(this->number_of_constraints());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        const IntervalConstraint& constraint=this->constraint(i);
        constraint_models.append(IntervalAffineModelConstraint(constraint.lower_bound(),affine_model(domain,constraint.function()),constraint.upper_bound()));
    }

    for(uint i=0; i!=space_models.size(); ++i) { space_models[i].set_error(0.0); }
    for(uint i=0; i!=constraint_models.size(); ++i) { constraint_models[i].function().set_error(0.0); }

    return IntervalAffineConstrainedImageSet(domain,space_models,constraint_models);
}


Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> IntervalConstrainedImageSet::split(uint j) const
{
    Pair<Box,Box> subdomains = Ariadne::split(this->_domain,j);
    subdomains.first=intersection(subdomains.first,this->_reduced_domain);
    subdomains.second=intersection(subdomains.second,this->_reduced_domain);

    Pair<IntervalConstrainedImageSet,IntervalConstrainedImageSet> result(
        IntervalConstrainedImageSet(subdomains.first,this->_function),
        IntervalConstrainedImageSet(subdomains.second,this->_function));

    for(uint i=0; i!=this->_constraints.size(); ++i) {
        result.first.new_parameter_constraint(this->_constraints[i]);
        result.second.new_parameter_constraint(this->_constraints[i]);
        //result.first.new_parameter_constraint(Ariadne::restrict(this->_negative_constraints[i],subdomains.first));
        //result.second.new_parameter_constraint(Ariadne::restrict(this->_negative_constraints[i],subdomains.second));
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

tribool IntervalConstrainedImageSet::separated(const Box& bx) const
{
    Box subdomain = this->_reduced_domain;
    IntervalVectorFunction function(this->dimension()+this->number_of_constraints(),this->number_of_parameters());
    for(uint i=0; i!=this->dimension(); ++i) { function[i]=this->_function[i]; }
    for(uint i=0; i!=this->number_of_constraints(); ++i) { function[i+this->dimension()]=this->_constraints[i].function(); }
    //IntervalVectorFunction function = join(this->_function,this->constraint_function());
    IntervalVector codomain = join(bx,this->constraint_bounds());
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return tribool(subdomain.empty()) || indeterminate;
}

tribool IntervalConstrainedImageSet::overlaps(const Box& bx) const
{
    //std::cerr<<"domain="<<this->_domain<<"\n";
    //std::cerr<<"subdomain="<<this->_reduced_domain<<"\n";

    Box subdomain = this->_reduced_domain;
    IntervalVectorFunction space_function = this->_function;
    IntervalVectorFunction constraint_function = this->constraint_function();
    IntervalVectorFunction function = join(space_function,constraint_function);
    //std::cerr<<"function="<<function<<"\n";
    IntervalVector constraint_bounds = intersection(this->constraint_bounds(),this->constraint_function().evaluate(subdomain));
    IntervalVector codomain = join(bx,constraint_bounds);
    //std::cerr<<"codomain="<<codomain<<"\n";
    NonlinearInfeasibleInteriorPointOptimiser optimiser;
    optimiser.verbosity=0;

    List<Pair<Nat,Box> > subdomains;
    Nat depth(0);
    Nat MAX_DEPTH=2;
    tribool feasible = false;
    subdomains.append(make_pair(depth,subdomain));

    while(!subdomains.empty()) {
        make_lpair(depth,subdomain)=subdomains.back();
        subdomains.pop_back();
        tribool found_feasible = optimiser.feasible(subdomain,function,codomain);
        if(definitely(found_feasible)) { return true; }
        if(possibly(found_feasible)) {
            if(depth==MAX_DEPTH) {
                feasible = indeterminate;
            } else {
                Pair<Box,Box> split_subdomains=subdomain.split();
                subdomains.append(make_pair(depth+1,split_subdomains.first));
                subdomains.append(make_pair(depth+1,split_subdomains.second));
            }
        }
    }
    return false;

}

void IntervalConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    const IntervalVector subdomain=this->_reduced_domain;
    const IntervalVectorFunction function = this->function();
    const IntervalVectorFunction constraint_function = this->constraint_function();
    const IntervalVector constraint_bounds = this->constraint_bounds();

    switch(DISCRETISATION_METHOD) {
        case SUBDIVISION_DISCRETISE:
            SubdivisionPaver().adjoin_outer_approximation(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        case AFFINE_DISCRETISE:
            AffinePaver().adjoin_outer_approximation(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        case CONSTRAINT_DISCRETISE:
            ConstraintPaver().adjoin_outer_approximation(paving,subdomain,function,constraint_function,constraint_bounds,depth);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown discretisation method\n");
    }

    paving.recombine();
}



tribool IntervalConstrainedImageSet::satisfies(const IntervalConstraint& nc) const
{
    if( subset(nc.function().evaluate(this->bounding_box()),nc.bounds()) ) {
        return true;
    }

    ConstraintSolver solver;
    const Box& domain=this->_domain;
    List<IntervalConstraint> all_constraints=this->constraints();
    IntervalScalarFunction composed_function = compose(nc.function(),this->_function);
    const Interval& bounds = nc.bounds();

    Tribool result;
    if(bounds.upper()<+inf) {
        all_constraints.append( composed_function >= bounds.upper() );
        result=solver.feasible(domain,all_constraints).first;
        all_constraints.pop_back();
        if(definitely(result)) { return false; }
    }
    if(bounds.lower()>-inf) {
        all_constraints.append(composed_function <= bounds.lower());
        result = result || solver.feasible(domain,all_constraints).first;
    }
    return !result;
}


void
IntervalConstrainedImageSet::draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    AffineDrawer(Depth(8)).draw(cnvs,proj,*this);
}

IntervalConstrainedImageSet
join(const IntervalConstrainedImageSet& set1, const IntervalConstrainedImageSet& set2)
{
    ARIADNE_ASSERT(set1.dimension()==set2.dimension());
    ARIADNE_ASSERT(set1.number_of_parameters()==set2.number_of_parameters());
    ARIADNE_ASSERT(set1.number_of_constraints()==set2.number_of_constraints());

    const uint np = set1.number_of_parameters();
    const uint nc = set1.number_of_parameters();

    const Box& domain1 = set1.domain();
    const Box& domain2 = set2.domain();

    uint join_parameter=np;
    for(uint i=0; i!=np; ++i) {
        if(domain1[i]!=domain2[i]) {
            if(domain1[i].upper() < domain2[i].lower() || domain1[i].lower() > domain2[i].upper()) {
                ARIADNE_FAIL_MSG("Cannot joint sets "<<set1<<" and "<<set2<<" since domains are separated.");
            }
            if(domain1[i].upper() >= domain2[i].lower() || domain1[i].lower() <= domain2[i].upper()) {
                if(join_parameter==np) {
                    join_parameter=i;
                } else {
                    ARIADNE_FAIL_MSG("Cannot joint sets "<<set1<<" and "<<set2<<" since domains do not join to form a box.");
                }
            }
        }
    }

    Box new_domain = hull(domain1,domain2);

    IntervalVectorFunctionModel function1
        = IntervalVectorFunctionModel( dynamic_cast<VectorFunctionModelInterface<Interval> const&>(set1.function().reference()));
    Vector<Float> function_error1=function1.errors();
    function1.clobber();
    function1.restrict(new_domain);

    IntervalVectorFunctionModel function2
        = IntervalVectorFunctionModel( dynamic_cast<VectorFunctionModelInterface<Interval> const&>(set2.function().reference()));
    Vector<Float> function_error2=function2.errors();
    function2.clobber();
    function2.restrict(new_domain);

    IntervalVectorFunctionModel new_function=(function1+function2)*ExactFloat(0.5);
    new_function.clobber();
    for(uint i=0; i!=new_function.result_size(); ++i) {
        function_error1[i]=add_up(norm(new_function[i]-function1[i]),function_error1[i]);
        function_error2[i]=add_up(norm(new_function[i]-function2[i]),function_error2[i]);
        Float new_function_error = max(function_error1[i],function_error2[i]);
        new_function[i] = new_function[i] + Interval(-new_function_error,+new_function_error);
    }

    ARIADNE_ASSERT(set1.number_of_constraints()==0);
    return IntervalConstrainedImageSet(new_domain,new_function);
}


std::ostream& IntervalConstrainedImageSet::write(std::ostream& os) const
{
    return os << "IntervalConstrainedImageSet( domain=" << this->domain() << ", function="<< this->function() << ", constraints=" << this->constraints() << " )";
}

std::ostream& operator<<(std::ostream& os, const IntervalConstrainedImageSet& set) {
    return set.write(os);
}







} // namespace Ariadne;
