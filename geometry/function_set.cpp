/***************************************************************************
 *            function_set.cc
 *
 *  Copyright 2008--17  Pieter Collins
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
#include "utility/logging.h"
#include "function/polynomial.h"
#include "function/function.h"
#include "function/taylor_function.h"
#include "function/formula.h"
#include "function/procedure.h"
#include "solvers/nonlinear_programming.h"
#include "solvers/constraint_solver.h"
#include "geometry/function_set.h"
#include "geometry/affine_set.h"
#include "geometry/grid_set.h"
#include "geometry/paver.h"
#include "algebra/algebra.h"

#include "output/graphics_interface.h"
#include "output/drawer.h"

namespace Ariadne {

static const Nat verbosity = 0u;

//! \related TaylorConstrainedImageSet \brief The possible types of method used to draw a nonlinear set.
enum DrawingMethod { CURVE_DRAW, BOX_DRAW, AFFINE_DRAW, GRID_DRAW };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DrawingMethod DRAWING_METHOD;
//! \related TaylorConstrainedImageSet \brief The accuracy used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern uint DRAWING_ACCURACY;

//! \related TaylorConstrainedImageSet \brief The possible types of method used to discretise a nonlinear set.
enum DiscretisationMethod { SUBDIVISION_DISCRETISE, AFFINE_DISCRETISE, CONSTRAINT_DISCRETISE };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to discretise a nonlinear set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DiscretisationMethod DISCRETISATION_METHOD;

DrawingMethod DRAWING_METHOD=AFFINE_DRAW;
DiscretisationMethod DISCRETISATION_METHOD=SUBDIVISION_DISCRETISE;
uint DRAWING_ACCURACY=1u;

template<class T> StringType str(const T& t) { StringStream ss; ss<<t; return ss.str(); }

Matrix<Float64> nonlinearities_zeroth_order(const ValidatedVectorFunction& f, const ExactBoxType& dom);
Pair<Nat,double> nonlinearity_index_and_error(const ValidatedVectorFunction& function, const ExactBoxType domain);
Pair<Nat,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const ExactBoxType& domain);

Matrix<Float64> nonlinearities_zeroth_order(const ValidatedVectorTaylorFunctionModel64& f, const ExactBoxType& dom)
{
    const Nat m=f.result_size();
    const Nat n=f.argument_size();
    ValidatedVectorTaylorFunctionModel64 g=restriction(f,dom);

    Matrix<Float64> nonlinearities=Matrix<Float64>::zero(m,n);
    MultiIndex a;
    for(Nat i=0; i!=m; ++i) {
        const ValidatedTaylorModel64& tm=g.model(i);
        for(ValidatedTaylorModel64::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
            a=iter->key();
            if(a.degree()>1) {
                for(Nat j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=mag(iter->data()).raw(); }
                }
            }
        }
    }

    return nonlinearities;
}

Matrix<Float64> nonlinearities_first_order(const ValidatedVectorFunction& f, const ExactBoxType& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const Nat m=f.result_size();
    const Nat n=f.argument_size();
    Vector<UpperIntervalDifferentialType> ivl_dx=UpperIntervalDifferentialType::constants(m,n, 1, dom);
    MultiIndex a(n);
    for(Nat i=0; i!=n; ++i) {
        Float64 sf=dom[i].radius().upper().raw();
        ++a[i];
        ivl_dx[i].expansion().append(a,ExactIntervalType(sf,sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<UpperIntervalDifferentialType> df=derivative_range(f,ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float64> nonlinearities=Matrix<Float64>::zero(m,n);
    for(Nat i=0; i!=m; ++i) {
        const UpperIntervalDifferentialType& d=df[i];
        for(UpperIntervalDifferentialType::ConstIterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==1) {
                for(Nat j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=iter->data().radius().raw(); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}

Matrix<Float64> nonlinearities_second_order(const ValidatedVectorFunction& f, const ExactBoxType& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const Nat m=f.result_size();
    const Nat n=f.argument_size();
    Vector<UpperIntervalDifferentialType> ivl_dx=UpperIntervalDifferentialType::constants(m,n, 2, dom);
    MultiIndex a(n);
    for(Nat i=0; i!=n; ++i) {
        Float64 sf=dom[i].radius().upper().raw();
        ++a[i];
        ivl_dx[i].expansion().append(a,ExactIntervalType(sf,sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<UpperIntervalDifferentialType> df=derivative_range(f,ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float64> nonlinearities=Matrix<Float64>::zero(m,n);
    for(Nat i=0; i!=m; ++i) {
        const UpperIntervalDifferentialType& d=df[i];
        for(UpperIntervalDifferentialType::ConstIterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==2) {
                for(Nat j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=mag(iter->data()).raw(); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}

Pair<Nat,double> nonlinearity_index_and_error(const ValidatedVectorTaylorFunctionModel64& function, const ExactBoxType domain) {
    Matrix<Float64> nonlinearities=Ariadne::nonlinearities_zeroth_order(function,domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    Nat jmax_in_row_imax=nonlinearities.column_size();
    Float64 max_row_sum=0.0;
    for(Nat i=0; i!=nonlinearities.row_size(); ++i) {
        Nat jmax=nonlinearities.column_size();
        Float64 row_sum=0.0;
        Float64 max_mag_j_in_i=0.0;
        for(Nat j=0; j!=nonlinearities.column_size(); ++j) {
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


ExactIntervalType over_approximating_interval(ValidatedLowerNumber const& lb, ValidatedUpperNumber const& ub) {
    Precision64 prec;
    return cast_exact_interval(UpperIntervalType(lb,ub,prec));
}

ExactBoxType under_approximation(const EffectiveBoxType& rbx) {
    Precision64 prec;
    return cast_exact_box(LowerBoxType(rbx,prec));
}

ExactBoxType over_approximation(const EffectiveBoxType& rbx) {
    Precision64 prec;
    return cast_exact_box(UpperBoxType(rbx,prec));
}

ExactBoxType approximation(const EffectiveBoxType& rbx) {
    Precision64 prec;
    return cast_exact_box(ApproximateBoxType(rbx,prec));
}


namespace {

Nat get_argument_size(const List<EffectiveConstraint>& c) {
    Nat as = ( c.size()==0 ? 0 : c[0].function().argument_size() );
    for(Nat i=0; i!=c.size(); ++i) {
        ARIADNE_ASSERT_MSG(c[i].function().argument_size()==as,"c="<<c);
    }
    return as;
}

EffectiveVectorFunction constraint_function(Nat as, const List<EffectiveConstraint>& c) {
    EffectiveVectorFunction f(c.size(),EuclideanDomain(as));
    for(Nat i=0; i!=c.size(); ++i) {
        //f[i]=c[i].function();
        f.set(i,c[i].function());
    }
    return f;
}

EffectiveBoxType constraint_bounds(const List<EffectiveConstraint>& c) {
    EffectiveBoxType b(c.size());
    for(Nat i=0; i!=c.size(); ++i) {
        b[i]=EffectiveIntervalType(c[i].lower_bound(),c[i].upper_bound());
    }
    return b;
}

List<EffectiveConstraint> constraints(const EffectiveVectorFunction& f, const EffectiveBoxType& b) {
    ARIADNE_ASSERT(f.result_size()==b.size());
    List<EffectiveConstraint> c; c.reserve(b.size());
    for(Nat i=0; i!=b.size(); ++i) {
        c.append(EffectiveConstraint(b[i].lower(),f[i],b[i].upper()));
    }
    return c;
}

} //namespace


ConstraintSet::ConstraintSet(const EffectiveVectorFunction& f, const EffectiveBoxType& b)
    : _dimension(f.argument_size()), _constraints()
{
    this->_constraints=Ariadne::constraints(f,b);
}

ConstraintSet::ConstraintSet(const List<EffectiveConstraint>& c)
    : _dimension(get_argument_size(c)), _constraints(c)
{
}

EffectiveVectorFunction const ConstraintSet::constraint_function() const
{
    return Ariadne::constraint_function(this->dimension(),this->constraints());
}

EffectiveBoxType const ConstraintSet::constraint_bounds() const
{
    return Ariadne::constraint_bounds(this->constraints());
}

ConstraintSet*
ConstraintSet::clone() const
{
    return new ConstraintSet(*this);
}


DimensionType
ConstraintSet::dimension() const
{
    return this->_dimension;
}


ValidatedSierpinskian
ConstraintSet::separated(const ExactBoxType& bx) const
{
    ExactBoxType codomain=over_approximation(this->codomain());
    return ValidatedConstrainedImageSet(bx,this->constraint_function()).separated(codomain);
}

ValidatedSierpinskian
ConstraintSet::overlaps(const ExactBoxType& bx) const
{
    ExactBoxType codomain=under_approximation(this->codomain());
    return ValidatedConstrainedImageSet(bx,this->constraint_function()).overlaps(codomain);
}

ValidatedSierpinskian
ConstraintSet::covers(const ExactBoxType& bx) const
{
    ExactBoxType codomain=under_approximation(this->codomain());
    return UpperBoxType(apply(this->constraint_function(),bx)).inside(codomain);
}


OutputStream&
ConstraintSet::write(OutputStream& os) const
{
    return os << "ConstraintSet( constraints=" << this->constraints() << " )";
}




BoundedConstraintSet::BoundedConstraintSet(const EffectiveBoxType& bx)
    : _domain(bx), _constraints()
{
}

BoundedConstraintSet::BoundedConstraintSet(const EffectiveBoxType& d, const EffectiveVectorFunction& f, const EffectiveBoxType& b)
    : _domain(d), _constraints(Ariadne::constraints(f,b))
{
    ARIADNE_ASSERT(b.size()==f.result_size());
    ARIADNE_ASSERT(d.size()==f.argument_size());
}

BoundedConstraintSet::BoundedConstraintSet(const EffectiveBoxType& d, const List<EffectiveConstraint>& c)
    : _domain(d), _constraints(c)
{
}

EffectiveVectorFunction const BoundedConstraintSet::constraint_function() const
{
    return Ariadne::constraint_function(this->dimension(),this->constraints());
}

EffectiveBoxType const BoundedConstraintSet::constraint_bounds() const
{
    return Ariadne::constraint_bounds(this->constraints());
}

BoundedConstraintSet*
BoundedConstraintSet::clone() const
{
    return new BoundedConstraintSet(*this);
}


DimensionType
BoundedConstraintSet::dimension() const
{
    return this->_domain.size();
}


ValidatedSierpinskian
BoundedConstraintSet::separated(const ExactBoxType& bx) const
{
    ExactBoxType domain=over_approximation(this->domain());
    if(Ariadne::disjoint(domain,bx)) { return true; }
    ExactBoxType codomain=over_approximation(this->codomain());
    return ValidatedConstrainedImageSet(Ariadne::intersection(bx,domain),this->constraint_function()).separated(codomain);
}


ValidatedSierpinskian
BoundedConstraintSet::overlaps(const ExactBoxType& bx) const
{
    if(Ariadne::disjoint(over_approximation(this->domain()),bx)) { return false; }
    ExactBoxType domain=under_approximation(this->domain());
    ExactBoxType codomain=under_approximation(this->codomain());
    return ValidatedConstrainedImageSet(Ariadne::intersection(bx,domain),this->constraint_function()).overlaps(codomain);
}


ValidatedSierpinskian
BoundedConstraintSet::covers(const ExactBoxType& bx) const
{
    ExactBoxType domain=under_approximation(this->domain());
    ExactBoxType codomain=under_approximation(this->codomain());
    if(!Ariadne::covers(domain,bx)) { return false; }
    return UpperBoxType(apply(this->constraint_function(),bx)).inside(codomain);
}

ValidatedSierpinskian
BoundedConstraintSet::inside(const ExactBoxType& bx) const
{
    return Ariadne::inside(UpperBoxType(over_approximation(this->domain())),bx);
}

UpperBoxType
BoundedConstraintSet::bounding_box() const
{
    ExactBoxType result=over_approximation(this->_domain);
    return widen(result);
}


OutputStream&
BoundedConstraintSet::write(OutputStream& os) const
{
    return os << "BoundedConstraintSet( domain=" << this->domain() << ", constraints=" << this->constraints() << ")";
}

Void
BoundedConstraintSet::draw(CanvasInterface& c, const Projection2d& p) const
{
    return ConstrainedImageSet(*this).draw(c,p);
}


BoundedConstraintSet
intersection(const ConstraintSet& cs,const EffectiveBoxType& bx)
{
    return BoundedConstraintSet(bx,cs.constraints());
}








ConstrainedImageSet::ConstrainedImageSet(const BoundedConstraintSet& set)
    : _domain(over_approximation(set.domain())), _function(EffectiveVectorFunction::identity(set.dimension()))
{
    for(Nat i=0; i!=set.number_of_constraints(); ++i) {
        this->new_parameter_constraint(set.constraint(i));
    }
}


const EffectiveVectorFunction ConstrainedImageSet::constraint_function() const
{
    EffectiveVectorFunction result(this->number_of_constraints(),EuclideanDomain(this->number_of_parameters()));
    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=this->constraint(i).function();
    }
    return result;
}

const EffectiveBoxType ConstrainedImageSet::constraint_bounds() const
{
    EffectiveBoxType result(this->number_of_constraints());
    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=EffectiveIntervalType(this->constraint(i).lower_bound(),this->constraint(i).upper_bound());
    }
    return result;
}


UpperBoxType ConstrainedImageSet::bounding_box() const
{
    return Ariadne::apply(this->_function,over_approximation(this->_domain));
}

Matrix<Float64Value> cast_exact(Matrix<Float64Bounds> vA) {
    Matrix<Float64Approximation> aA=vA;
    return reinterpret_cast<Matrix<Float64Value>&>(aA);
}

ValidatedAffineConstrainedImageSet
ConstrainedImageSet::affine_approximation() const
{
    Precision64 prec;
    const Vector<ExactIntervalType> D=approximation(this->domain());
    Vector<Float64Value> m=midpoint(D);
    Matrix<Float64Value> G=cast_exact(jacobian(this->_function,m));
    Vector<Float64Value> h=cast_exact(this->_function.evaluate(m)-G*m);
    ValidatedAffineConstrainedImageSet result(D,G,h);


    Vector<Float64Approximation> a(this->number_of_parameters());
    Float64Approximation b,l,u;
    for(List<EffectiveConstraint>::ConstIterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        AffineModel<ValidatedTag,Float64> a=affine_model(D,iter->function(),prec);
        ExactIntervalType b=iter->bounds();
        result.new_constraint(b.lower()<=a<=b.upper());
    }

    return result;
}


ValidatedKleenean ConstrainedImageSet::satisfies(const EffectiveConstraint& nc) const
{
    if( definitely(subset(Ariadne::apply(nc.function(),this->bounding_box()),nc.bounds())) ) {
        return true;
    }

    ConstraintSolver solver;
    const EffectiveBoxType& domain=this->_domain;
    List<EffectiveConstraint> all_constraints=this->_constraints;
    EffectiveScalarFunction composed_function = compose(nc.function(),this->_function);
    const EffectiveNumber& lower_bound = nc.lower_bound();
    const EffectiveNumber& upper_bound = nc.upper_bound();

    ValidatedKleenean result;
    if(definitely(upper_bound<+infinity)) {
        all_constraints.append( composed_function >= upper_bound );
        result=solver.feasible(over_approximation(domain),all_constraints).first;
        all_constraints.pop_back();
        if(definitely(result)) { return false; }
    }
    if(definitely(lower_bound>-infinity)) {
        all_constraints.append(composed_function <= lower_bound);
        result = result || solver.feasible(over_approximation(domain),all_constraints).first;
    }
    return !result;
}


//! \brief Test if the set is contained in (the interior of) a box.
ValidatedSierpinskian ConstrainedImageSet::inside(const ExactBoxType& bx) const {
    return this->bounding_box().inside(bx);
}

ValidatedSierpinskian ConstrainedImageSet::separated(const ExactBoxType& bx) const
{
    UpperBoxType subdomain = over_approximation(this->_domain);
    EffectiveVectorFunction function = join(this->function(),this->constraint_function());
    ExactBoxType codomain = product(bx,ExactBoxType(over_approximation(this->constraint_bounds())));
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return subdomain.is_empty();


}


ValidatedSierpinskian ConstrainedImageSet::overlaps(const ExactBoxType& bx) const
{
    return ValidatedConstrainedImageSet(under_approximation(this->_domain),this->_function,this->_constraints).overlaps(bx);
}


Void
ConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const
{
    ARIADNE_ASSERT(paving.dimension()==this->dimension());
    const ExactBoxType domain=over_approximation(this->domain());
    const EffectiveVectorFunction& space_function=this->function();
    EffectiveVectorFunction constraint_function(this->number_of_constraints(),domain);
    EffectiveBoxType codomain(this->number_of_constraints());

    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        constraint_function.set(i,this->_constraints[i].function());
        codomain[i]=EffectiveIntervalType(this->_constraints[i].lower_bound(),this->_constraints[i].upper_bound());
    }
    ExactBoxType constraint_bounds=over_approximation(codomain);

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


Pair<ConstrainedImageSet,ConstrainedImageSet>
ConstrainedImageSet::split() const
{
    ARIADNE_ASSERT(this->number_of_parameters()>0);
    auto rmax(this->domain()[0].radius());
    Nat k=0;
    for(Nat j=1; j!=this->number_of_parameters(); ++j) {
       auto rj(this->domain()[j].radius());
       if(decide(rj>rmax)) {
            k=j;
            rmax=rj;
        }
    }
    return this->split(k);
}

Pair<ConstrainedImageSet,ConstrainedImageSet>
ConstrainedImageSet::split(Nat j) const
{
    EffectiveIntervalType interval = this->domain()[j];
    Real midpoint = interval.midpoint();
    Pair<EffectiveBoxType,EffectiveBoxType> subdomains(this->domain(),this->domain());
    subdomains.first[j]=EffectiveIntervalType(interval.lower(),midpoint);
    subdomains.second[j]=EffectiveIntervalType(midpoint,interval.upper());
    return make_pair(ConstrainedImageSet(subdomains.first,this->_function,this->_constraints),
                     ConstrainedImageSet(subdomains.second,this->_function,this->_constraints));
}


ConstrainedImageSet image(const BoundedConstraintSet& set, const EffectiveVectorFunction& function) {
    ARIADNE_ASSERT(set.dimension()==function.argument_size());
    ConstrainedImageSet result(set.domain(),function);
    for(Nat i=0; i!=set.number_of_constraints(); ++i) {
        result.new_parameter_constraint(set.constraint(i));
    }
    return result;
}



















Matrix<Float64> nonlinearities_zeroth_order(const ValidatedVectorTaylorFunctionModel64& f, const ExactBoxType& dom);


Matrix<Float64> nonlinearities_zeroth_order(const ValidatedVectorFunction& f, const ExactBoxType& dom)
{
    ARIADNE_ASSERT(dynamic_cast<const ValidatedVectorTaylorFunctionModel64*>(f.raw_pointer()));
    return nonlinearities_zeroth_order(dynamic_cast<const ValidatedVectorTaylorFunctionModel64&>(*f.raw_pointer()),dom);
}

/*
Matrix<Float64> nonlinearities_first_order(const ValidatedVectorFunction& f, const ExactBoxType& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const Nat m=f.result_size();
    const Nat n=f.argument_size();
    Vector<UpperIntervalDifferentialType> ivl_dx=UpperIntervalDifferentialType::constants(m,n, 1, dom);
    MultiIndex a(n);
    for(Nat i=0; i!=n; ++i) {
        Float64 sf=dom[i].radius();
        ++a[i];
        ivl_dx[i].expansion().append(a,ExactIntervalType(sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<UpperIntervalDifferentialType> df=f.evaluate(ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float64> nonlinearities=Matrix<Float64>::zero(m,n);
    for(Nat i=0; i!=m; ++i) {
        const UpperIntervalDifferentialType& d=df[i];
        for(UpperIntervalDifferentialType::ConstIterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==1) {
                for(Nat j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=radius(iter->data()); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}

Matrix<Float64> nonlinearities_second_order(const ValidatedVectorFunction& f, const ExactBoxType& dom)
{
    //std::cerr<<"\n\nf="<<f<<"\n";
    //std::cerr<<"dom="<<dom<<"\n";
    const Nat m=f.result_size();
    const Nat n=f.argument_size();
    Vector<UpperIntervalDifferentialType> ivl_dx=UpperIntervalDifferentialType::constants(m,n, 2, dom);
    MultiIndex a(n);
    for(Nat i=0; i!=n; ++i) {
        Float64 sf=dom[i].radius();
        ++a[i];
        ivl_dx[i].expansion().append(a,ExactIntervalType(sf));
        --a[i];
    }
    //std::cerr<<"dx="<<ivl_dx<<"\n";
    Vector<UpperIntervalDifferentialType> df=f.evaluate(ivl_dx);
    //std::cerr<<"df="<<df<<"\n";

    Matrix<Float64> nonlinearities=Matrix<Float64>::zero(m,n);
    for(Nat i=0; i!=m; ++i) {
        const UpperIntervalDifferentialType& d=df[i];
        for(UpperIntervalDifferentialType::ConstIterator iter=d.begin(); iter!=d.end(); ++iter) {
            a=iter->key();
            if(a.degree()==2) {
                for(Nat j=0; j!=n; ++j) {
                    if(a[j]>0) { nonlinearities[i][j]+=mag(iter->data()); }
                }
            }
        }
    }
    //std::cerr<<"nonlinearities="<<nonlinearities<<"\n";

    return nonlinearities;
}
*/

Pair<Nat,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const ExactBoxType& domain)
{
    Matrix<UpperIntervalType> jacobian=Ariadne::jacobian_range(function,domain);

    // Compute the column of the matrix which has the norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    Nat jmax=domain.size();
    Float64 max_column_norm=0.0;
    for(Nat j=0; j!=domain.size(); ++j) {
        Float64 column_norm=0.0;
        for(Nat i=0; i!=function.result_size(); ++i) {
            column_norm+=mag(jacobian[i][j]).raw();
        }
        column_norm *= domain[j].radius().upper().raw();
        if(column_norm>max_column_norm) {
            max_column_norm=column_norm;
            jmax=j;
        }
    }
    return make_pair(jmax,numeric_cast<double>(max_column_norm));
}

Pair<Nat,double> nonlinearity_index_and_error(const ValidatedVectorFunction& function, const ExactBoxType& domain)
{
    Matrix<Float64> nonlinearities=Ariadne::nonlinearities_zeroth_order(function,domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    Nat jmax_in_row_imax=nonlinearities.column_size();
    Float64 max_row_sum=0.0;
    for(Nat i=0; i!=nonlinearities.row_size(); ++i) {
        Nat jmax=nonlinearities.column_size();
        Float64 row_sum=0.0;
        Float64 max_mag_j_in_i=0.0;
        for(Nat j=0; j!=nonlinearities.column_size(); ++j) {
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








Void
ConstrainedImageSet::draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    ValidatedConstrainedImageSet(approximation(this->_domain),this->_function,this->_constraints).draw(cnvs,proj);
}



OutputStream&
ConstrainedImageSet::write(OutputStream& os) const
{
    return os << "ConstrainedImageSet( domain=" << this->_domain
              << ", function=" << this->_function << ", constraints=" << this->_constraints << " )";
}




template<class SF> struct FunctionTraits;
template<class X> struct FunctionTraits< ScalarFunction<X> > { typedef VectorFunction<X> VectorFunctionType; };
template<> struct FunctionTraits< ValidatedScalarTaylorFunctionModel64 > { typedef ValidatedVectorTaylorFunctionModel64 VectorFunctionType; };

template<class SF> class TemplatedConstraintSet;
template<class SF> class TemplatedConstrainedImageSet;



ValidatedVectorFunction ValidatedConstrainedImageSet::constraint_function() const
{
    ValidatedVectorFunction result(this->number_of_constraints(),EuclideanDomain(this->number_of_parameters()));
    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=this->constraint(i).function();
    }
    return result;
}

ExactBoxType ValidatedConstrainedImageSet::constraint_bounds() const
{
    Precision64 prec;
    ExactBoxType result(this->number_of_constraints());
    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        result[i]=over_approximating_interval(this->constraint(i).lower_bound(),this->constraint(i).upper_bound());
    }
    return result;
}


UpperBoxType
ValidatedConstrainedImageSet::bounding_box() const
{
    return Ariadne::apply(this->_function,this->_reduced_domain);
}


ValidatedAffineConstrainedImageSet
ValidatedConstrainedImageSet::affine_over_approximation() const
{
    typedef List<ValidatedConstraint>::ConstIterator ConstIterator;
    Precision64 prec;

    Vector<ExactIntervalType> domain = this->domain();
    Vector<ValidatedAffineModel> space_models=affine_models(domain,this->function(),prec);
    List<ValidatedAffineModelConstraint> constraint_models;
    constraint_models.reserve(this->number_of_constraints());
    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        const ValidatedConstraint& constraint=this->constraint(i);
        auto am=affine_model(domain,constraint.function(),prec);
        constraint_models.append(ValidatedAffineModelConstraint(constraint.lower_bound(),am,constraint.upper_bound()));
    }

    return ValidatedAffineConstrainedImageSet(domain,space_models,constraint_models);

/*
    const Nat nx=this->dimension();
    //const Nat nnc=this->_negative_constraints.size();
    //const Nat nzc=this->_zero_constraints.size();
    const Nat np=this->number_of_parameters();

    // Compute the number of values with a nonzero error
    Nat nerr=0;
    for(Nat i=0; i!=nx; ++i) {
        if(function[i].error()>0.0) { ++nerr; }
    }

    Vector<Float64> h(nx);
    Matrix<Float64> G(nx,np+nerr);
    Nat ierr=0; // The index where the error bound should go
    for(Nat i=0; i!=nx; ++i) {
        ValidatedScalarTaylorFunctionModel64 component_function=function[i];
        h[i]=component_function.model().value();
        for(Nat j=0; j!=np; ++j) {
            G[i][j]=component_function.model().gradient(j);
        }
        if(component_function.model().error()>0.0) {
            G[i][np+ierr]=component_function.model().error();
            ++ierr;
        }
    }

    ValidatedAffineConstrainedImageSet result(G,h);

    Vector<Float64> a(np+nerr, 0.0);
    Float64 b;

    for(ConstIterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        ValidatedScalarTaylorFunctionModel64 constraint_function(this->_reduced_domain,iter->function(),affine_sweeper);
        b=sub_up(constraint_function.model().error(),constraint_function.model().value());
        for(Nat j=0; j!=np; ++j) { a[j]=constraint_function.model().gradient(j); }
        result.new_parameter_constraint(-inf,a,b);
    }

    ARIADNE_NOT_IMPLEMENTED;

    ARIADNE_LOG(2,"set="<<*this<<"\nset.affine_over_approximation()="<<result<<"\n");
    return result;
*/
}

ValidatedAffineConstrainedImageSet ValidatedConstrainedImageSet::affine_approximation() const
{
    Precision64 prec;

    Vector<ExactIntervalType> domain = this->domain();

    Vector<ValidatedAffineModel> space_models=affine_models(domain,this->function(),prec);
    List<ValidatedAffineModelConstraint> constraint_models;
    constraint_models.reserve(this->number_of_constraints());
    for(Nat i=0; i!=this->number_of_constraints(); ++i) {
        const ValidatedConstraint& constraint=this->constraint(i);
        auto am=affine_model(domain,constraint.function(),prec);
        constraint_models.append(ValidatedAffineModelConstraint(constraint.lower_bound(),am,constraint.upper_bound()));
    }

    for(Nat i=0; i!=space_models.size(); ++i) { space_models[i].set_error(0u); }
    for(Nat i=0; i!=constraint_models.size(); ++i) { constraint_models[i].function().set_error(0u); }

    return ValidatedAffineConstrainedImageSet(domain,space_models,constraint_models);
}


Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> ValidatedConstrainedImageSet::split(Nat j) const
{
    Pair<ExactBoxType,ExactBoxType> subdomains = Ariadne::split(this->_domain,j);
    subdomains.first=intersection(subdomains.first,this->_reduced_domain);
    subdomains.second=intersection(subdomains.second,this->_reduced_domain);

    Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet> result(
        ValidatedConstrainedImageSet(subdomains.first,this->_function),
        ValidatedConstrainedImageSet(subdomains.second,this->_function));

    for(Nat i=0; i!=this->_constraints.size(); ++i) {
        result.first.new_parameter_constraint(this->_constraints[i]);
        result.second.new_parameter_constraint(this->_constraints[i]);
        //result.first.new_parameter_constraint(Ariadne::restrict(this->_negative_constraints[i],subdomains.first));
        //result.second.new_parameter_constraint(Ariadne::restrict(this->_negative_constraints[i],subdomains.second));
    }
    return result;
}

Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>
ValidatedConstrainedImageSet::split() const
{
    Nat k=this->number_of_parameters();
    Float64 rmax=0.0;
    for(Nat j=0; j!=this->number_of_parameters(); ++j) {
        Float64UpperBound rj=this->domain()[j].radius();
        if(rj.raw()>rmax) {
            k=j;
            rmax=rj.raw();
        }
    }
    return this->split(k);
}


Void
ValidatedConstrainedImageSet::reduce()
{
    ConstraintSolver solver;
    UpperBoxType& reduced_domain=reinterpret_cast<UpperBoxType&>(this->_reduced_domain);
    solver.reduce(reduced_domain, this->constraint_function(), this->constraint_bounds());
}

ValidatedKleenean ValidatedConstrainedImageSet::is_empty() const
{
    const_cast<ValidatedConstrainedImageSet*>(this)->reduce();
    return this->_reduced_domain.is_empty();
}

ValidatedSierpinskian ValidatedConstrainedImageSet::inside(const ExactBoxType& bx) const
{
    return Ariadne::inside(this->bounding_box(),bx);
}

ValidatedSierpinskian ValidatedConstrainedImageSet::separated(const ExactBoxType& bx) const
{
    UpperBoxType subdomain = this->_reduced_domain;
    ValidatedVectorFunction function(this->dimension()+this->number_of_constraints(),EuclideanDomain(this->number_of_parameters()));
    for(Nat i=0; i!=this->dimension(); ++i) { function[i]=this->_function[i]; }
    for(Nat i=0; i!=this->number_of_constraints(); ++i) { function[i+this->dimension()]=this->_constraints[i].function(); }
    //ValidatedVectorFunction function = join(this->_function,this->constraint_function());
    ExactBoxType codomain = product(bx,this->constraint_bounds());
    ConstraintSolver solver;
    solver.reduce(subdomain,function,codomain);
    return subdomain.is_empty();
}

ValidatedSierpinskian ValidatedConstrainedImageSet::overlaps(const ExactBoxType& bx) const
{
    //std::cerr<<"domain="<<this->_domain<<"\n";
    //std::cerr<<"subdomain="<<this->_reduced_domain<<"\n";

    ExactBoxType subdomain = this->_reduced_domain;
    ValidatedVectorFunction space_function = this->_function;
    ValidatedVectorFunction constraint_function = this->constraint_function();
    ValidatedVectorFunction function = join(space_function,constraint_function);
    //std::cerr<<"function="<<function<<"\n";
    ExactBoxType constraint_bounds = intersection(this->constraint_bounds(),cast_exact_box(Ariadne::apply(this->constraint_function(),subdomain)));
    ExactBoxType codomain = product(bx,constraint_bounds);
    //std::cerr<<"codomain="<<codomain<<"\n";
    NonlinearInfeasibleInteriorPointOptimiser optimiser;
    optimiser.verbosity=0;

    List<Pair<Nat,ExactBoxType> > subdomains;
    Nat depth(0);
    Nat MAX_DEPTH=2;
    ValidatedKleenean feasible = false;
    subdomains.append(make_pair(depth,subdomain));

    while(!subdomains.empty()) {
        make_lpair(depth,subdomain)=subdomains.back();
        subdomains.pop_back();
        ValidatedKleenean found_feasible = optimiser.feasible(subdomain,function,codomain);
        if(definitely(found_feasible)) { return true; }
        if(possibly(found_feasible)) {
            if(depth==MAX_DEPTH) {
                feasible = ValidatedKleenean(indeterminate);
            } else {
                Pair<ExactBoxType,ExactBoxType> split_subdomains=subdomain.split();
                subdomains.append(make_pair(depth+1,split_subdomains.first));
                subdomains.append(make_pair(depth+1,split_subdomains.second));
            }
        }
    }
    return false;

}

Void ValidatedConstrainedImageSet::adjoin_outer_approximation_to(PavingInterface& paving, Int depth) const
{
    const ExactBoxType subdomain=this->_reduced_domain;
    const ValidatedVectorFunction function = this->function();
    const ValidatedVectorFunction constraint_function = this->constraint_function();
    const ExactBoxType constraint_bounds = this->constraint_bounds();

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



ValidatedKleenean ValidatedConstrainedImageSet::satisfies(const ValidatedConstraint& nc) const
{
    if( definitely(subset(Ariadne::apply(nc.function(),this->bounding_box()),nc.bounds())) ) {
        return true;
    }

    ConstraintSolver solver;
    const ExactBoxType& domain=this->_domain;
    List<ValidatedConstraint> all_constraints=this->constraints();
    ValidatedScalarFunction composed_function = compose(nc.function(),this->_function);
    const ExactIntervalType& bounds = nc.bounds();

    ValidatedKleenean result;
    if(definitely(bounds.upper()<+infty)) {
        all_constraints.append( composed_function >= bounds.upper() );
        result=solver.feasible(domain,all_constraints).first;
        all_constraints.pop_back();
        if(definitely(result)) { return false; }
    }
    if(definitely(bounds.lower()>-infty)) {
        all_constraints.append(composed_function <= bounds.lower());
        result = result || solver.feasible(domain,all_constraints).first;
    }
    return !result;
}


Void
ValidatedConstrainedImageSet::draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    AffineDrawer(Depth(8)).draw(cnvs,proj,*this);
}

Void
ValidatedConstrainedImageSet::box_draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    BoxDrawer().draw(cnvs,proj,*this);
}

Void
ValidatedConstrainedImageSet::affine_draw(CanvasInterface& cnvs, const Projection2d& proj, Int depth) const
{
    AffineDrawer(depth).draw(cnvs,proj,*this);
}

Void
ValidatedConstrainedImageSet::grid_draw(CanvasInterface& cnvs, const Projection2d& proj) const
{
    GridDrawer().draw(cnvs,proj,*this);
}

ValidatedConstrainedImageSet
join(const ValidatedConstrainedImageSet& set1, const ValidatedConstrainedImageSet& set2)
{
    ARIADNE_ASSERT(set1.dimension()==set2.dimension());
    ARIADNE_ASSERT(set1.number_of_parameters()==set2.number_of_parameters());
    ARIADNE_ASSERT(set1.number_of_constraints()==set2.number_of_constraints());

    const Nat np = set1.number_of_parameters();
    const Nat nc = set1.number_of_parameters();

    const ExactBoxType& domain1 = set1.domain();
    const ExactBoxType& domain2 = set2.domain();

    Nat join_parameter=np;
    for(Nat i=0; i!=np; ++i) {
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

    ExactBoxType new_domain = hull(domain1,domain2);

    ValidatedVectorFunctionModel64 function1
        = ValidatedVectorFunctionModel64( dynamic_cast<VectorFunctionModel64Interface<ValidatedTag> const&>(set1.function().reference()));
    Vector<Float64Error> function_error1=function1.errors();
    function1.clobber();
    function1.restrict(new_domain);

    ValidatedVectorFunctionModel64 function2
        = ValidatedVectorFunctionModel64( dynamic_cast<VectorFunctionModel64Interface<ValidatedTag> const&>(set2.function().reference()));
    Vector<Float64Error> function_error2=function2.errors();
    function2.clobber();
    function2.restrict(new_domain);

    ValidatedVectorFunctionModel64 new_function=(function1+function2)*Float64Value(0.5);
    new_function.clobber();
    for(Nat i=0; i!=new_function.result_size(); ++i) {
        function_error1[i]=norm(new_function[i]-function1[i])+function_error1[i];
        function_error2[i]=norm(new_function[i]-function2[i])+function_error2[i];
        Float64Error new_function_error = max(function_error1[i],function_error2[i]);
        new_function[i] = new_function[i] + Float64Bounds(-new_function_error,+new_function_error);
    }

    ARIADNE_ASSERT(set1.number_of_constraints()==0);
    return ValidatedConstrainedImageSet(new_domain,new_function);
}


OutputStream& ValidatedConstrainedImageSet::write(OutputStream& os) const
{
    return os << "ValidatedConstrainedImageSet( domain=" << this->domain() << ", function="<< this->function() << ", constraints=" << this->constraints() << " )";
}

OutputStream& operator<<(OutputStream& os, const ValidatedConstrainedImageSet& set) {
    return set.write(os);
}







} // namespace Ariadne;
