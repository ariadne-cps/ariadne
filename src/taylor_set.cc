/***************************************************************************
 *            taylor_set.cc
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

#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "differential.h"
#include "polynomial.h"
#include "function.h"

#include "taylor_model.h"
#include "taylor_function.h"

#include "box.h"

#include "function_set.h"
#include "taylor_set.h"
#include "affine_set.h"

#include "list_set.h"
#include "grid_set.h"

#include "constraint_solver.h"
#include "nonlinear_programming.h"

#include "graphics_interface.h"

#include "discrete_event.h"

#include "logging.h"

#include "config.h"
#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif // HAVE_CAIRO_H


namespace Ariadne {

static const uint verbosity = 0u;

DrawingMethod DRAWING_METHOD=AFFINE_DRAW;
unsigned int DRAWING_ACCURACY=1u;

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;




IntervalTaylorModel& operator-=(IntervalTaylorModel& tm, const MultiIndex& a) {
    for(IntervalTaylorModel::iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        iter->key()-=a;
    }
    return tm;
}

// TODO: Make more efficient
inline void assign_all_but_last(MultiIndex& r, const MultiIndex& a) {
    for(uint i=0; i!=r.size(); ++i) { r[i]=a[i]; }
}


void TaylorConstrainedImageSet::_check() const {
    ARIADNE_ASSERT_MSG(this->_function.argument_size()==this->domain().size(),*this);
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->argument_size()==this->domain().size(),*this);
    }
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->argument_size()==this->domain().size(),*this);
    }
}

// FIXME: What if solving for constraint leaves domain?
void TaylorConstrainedImageSet::_solve_zero_constraints() {
    this->_check();
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ) {
        const Vector<Interval>& domain=this->domain();
        const IntervalTaylorModel& model=iter->model();
        const uint k=model.argument_size()-1u;
        IntervalTaylorModel zeroth_order(k);
        IntervalTaylorModel first_order(k);
        bool is_zeroth_order=true;
        bool is_first_order=true;
        MultiIndex r(k);
        // Try linear approach in last coefficient
        for(IntervalTaylorModel::const_iterator tmiter=model.begin(); tmiter!=model.end(); ++tmiter) {
            if(tmiter->key()[k]==0) {
                assign_all_but_last(r,tmiter->key());
                zeroth_order.expansion().append(r,tmiter->data());
            } else if(tmiter->key()[k]==1) {
                is_zeroth_order=false;
                assign_all_but_last(r,tmiter->key());
                first_order.expansion().append(r,tmiter->data());
            } else {
                is_first_order=false; break;
            }
        }
        if(is_first_order && !is_zeroth_order) {
            const Vector<Interval> new_domain=project(domain,range(0,k));
            IntervalTaylorModel substitution_model=-zeroth_order/first_order;
            this->_function=VectorTaylorFunction(new_domain,Ariadne::substitute(this->_function.models(),k,substitution_model));
            for(List<ScalarTaylorFunction>::iterator constraint_iter=this->_constraints.begin();
                    constraint_iter!=this->_constraints.end(); ++constraint_iter) {
                ScalarTaylorFunction& constraint=*constraint_iter;
                constraint=ScalarTaylorFunction(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            for(List<ScalarTaylorFunction>::iterator constraint_iter=this->_equations.begin();
                    constraint_iter!=this->_equations.end(); ++constraint_iter) {
                ScalarTaylorFunction& constraint=*constraint_iter;
                constraint=ScalarTaylorFunction(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            // Since we are using an std::vector, assign iterator to next element
            iter=this->_equations.erase(iter);
            this->_check();
        } else {
            ARIADNE_WARN("No method for solving constraint "<<*iter<<" currently implemented.");
            ++iter;
        }
    }
}


TaylorConstrainedImageSet::TaylorConstrainedImageSet()
    : _domain(), _function(), _reduced_domain()
{
}

TaylorConstrainedImageSet* TaylorConstrainedImageSet::clone() const
{
    return new TaylorConstrainedImageSet(*this);
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const Box& box)
{
    // Ensure domain elements have nonempty radius
    const float min_float=std::numeric_limits<float>::min();
    List<uint> proper_coordinates;
    proper_coordinates.reserve(box.dimension());
    for(uint i=0; i!=box.dimension(); ++i) {
        if(box[i].radius()>=min_float) {
            proper_coordinates.append(i);
        }
    }
    this->_domain=IntervalVector(proper_coordinates.size());
    for(uint j=0; j!=this->_domain.size(); ++j) {
        this->_domain[j]=box[proper_coordinates[j]];
    }

    // HACK: Make a dummy variable for the domain to avoid bugs which occur
    // with a zero-dimensional domain.
    // FIXME: Fix issues with TaylorFunction on zero-dimensional domain.
    if(proper_coordinates.size()==0) { this->_domain=IntervalVector(1u,Interval(-1,+1)); }


    this->_function=VectorTaylorFunction(box.dimension(),this->_domain);
    uint j=0;
    proper_coordinates.append(box.dimension());
    for(uint i=0; i!=box.dimension(); ++i) {
        if(proper_coordinates[j]==i) {
            this->_function[i]=ScalarTaylorFunction::coordinate(this->_domain,j);
            ++j;
        } else {
            this->_function[i]=ScalarTaylorFunction::constant(this->_domain,box[i]);
        }
    }
    this->_reduced_domain=this->_domain;
}


TaylorConstrainedImageSet::TaylorConstrainedImageSet(const IntervalVector& domain, const RealVectorFunction& function)
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    this->_domain=domain;
    this->_function=VectorTaylorFunction(this->_domain,function);
    this->_reduced_domain=this->_domain;
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const IntervalVector& domain, const RealVectorFunction& function, const List<NonlinearConstraint>& constraints)
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    const double min=std::numeric_limits<double>::min();
    this->_domain=domain;
    for(uint i=0; i!=this->_domain.size(); ++i) {
        if(this->_domain[i].radius()==0) {
            this->_domain[i]+=Interval(-min,+min);
        }
    }

    this->_function=VectorTaylorFunction(this->_domain,function);

    for(uint i=0; i!=constraints.size(); ++i) {
        ARIADNE_ASSERT_MSG(domain.size()==constraints[i].function().argument_size(),"domain="<<domain<<", constraint="<<constraints[i]);
        if(constraints[i].bounds().singleton()) {
            this->new_equality_constraint(constraints[i].function()-Real(constraints[i].bounds().midpoint()));
        } else {
            if(constraints[i].bounds().lower()>-inf<Float>()) {
                this->new_negative_constraint(Real(constraints[i].bounds().lower())-constraints[i].function());
            }
            if(constraints[i].bounds().upper()<+inf<Float>()) {
                this->new_negative_constraint(constraints[i].function()-Real(constraints[i].bounds().upper()));
            }
        }
    }

    this->_reduced_domain=domain;
    this->reduce();

}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const IntervalVector& domain, const RealVectorFunction& function, const NonlinearConstraint& constraint)
{
    *this=TaylorConstrainedImageSet(domain,function,make_list(constraint));
}



TaylorConstrainedImageSet::TaylorConstrainedImageSet(const VectorTaylorFunction& function)
{
    this->_domain=function.domain();
    this->_function=function;
    this->_reduced_domain=this->_domain;
}



// Returns true if the entire set is positive; false if entire set is negative
tribool TaylorConstrainedImageSet::satisfies(RealScalarFunction constraint) const
{
    Interval constraint_range=constraint(this->codomain());
    if(constraint_range.upper()<0.0) { return false; }
    if(constraint_range.lower()>0.0) { return true; }
    return indeterminate;
}


void TaylorConstrainedImageSet::substitute(uint j, ScalarTaylorFunction v)
{
    ARIADNE_ASSERT_MSG(v.argument_size()+1u==this->number_of_parameters(),
                       "number_of_parameters="<<this->number_of_parameters()<<", variable="<<v);
                       this->_function = Ariadne::substitute(this->_function,j,v);
                       for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }
                       for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }

                       this->_check();
}

void TaylorConstrainedImageSet::substitute(uint j, Float c)
{
    this->_function = Ariadne::partial_evaluate(this->_function,j,c);
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
    this->_check();
}

void TaylorConstrainedImageSet::apply_map(const IntervalVectorFunctionInterface& map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    VectorTaylorFunction& function=this->_function;
    function=compose(map,function);
    this->_check();
}

void TaylorConstrainedImageSet::apply_map(RealVectorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    VectorTaylorFunction& function=this->_function;
    function=compose(map,function);
    this->_check();
}

void TaylorConstrainedImageSet::apply_map(VectorTaylorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    VectorTaylorFunction& function=this->_function;
    function=compose(map,function);
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow(RealVectorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,combine(this->_function,VectorTaylorFunction::identity(Vector<Interval>(1u,time))));
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow_step(VectorTaylorFunction flow, Float time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,join(this->_function,ScalarTaylorFunction::constant(this->_function.domain(),time)));
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow(VectorTaylorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,combine(this->_function,ScalarTaylorFunction::identity(time)));
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}

void TaylorConstrainedImageSet::new_state_constraint(NonlinearConstraint constraint) {
    Float infty=+inf<Float>();
    Interval interval=constraint.bounds();
    ScalarTaylorFunction composed_function=compose(constraint.function(),this->_function);
    if(interval.lower()==0.0 && interval.upper()==0.0) {
        this->new_zero_constraint(composed_function);
    } else if(interval.lower()==0.0 && interval.upper()==infty) {
        this->new_negative_constraint(-composed_function);
    } else if(interval.lower()==-infty && interval.upper()==0.0) {
        this->new_negative_constraint(composed_function);
    }
}

void TaylorConstrainedImageSet::new_negative_constraint(RealScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_constraints.append(ScalarTaylorFunction(this->domain(),constraint));
}

void TaylorConstrainedImageSet::new_negative_constraint(ScalarTaylorFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.domain()==this->domain(),std::setprecision(17)<<"domain="<<this->domain()<<", constraint="<<constraint);
    this->_constraints.append(constraint);
}

void TaylorConstrainedImageSet::new_equality_constraint(RealScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_equations.append(ScalarTaylorFunction(this->domain(),constraint));
}

void TaylorConstrainedImageSet::new_zero_constraint(RealScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_equations.append(ScalarTaylorFunction(this->domain(),constraint));
}

List<ScalarTaylorFunction> const&
TaylorConstrainedImageSet::negative_constraints() const {
    return this->_constraints;
}

void TaylorConstrainedImageSet::new_zero_constraint(ScalarTaylorFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.domain()==this->domain(),std::setprecision(17)<<"domain="<<this->domain()<<", constraint="<<constraint);
    this->_equations.append(constraint);
}

List<ScalarTaylorFunction> const&
TaylorConstrainedImageSet::zero_constraints() const {
    return this->_equations;
}

List<NonlinearConstraint>
TaylorConstrainedImageSet::constraints() const {
    List<NonlinearConstraint> result;
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        result.append(this->_constraints[i]<=0.0);
    }
    for(uint i=0; i!=this->_equations.size(); ++i) {
        result.append(this->_equations[i]==0.0);
    }
    return result;
}

uint TaylorConstrainedImageSet::number_of_constraints() const {
    return this->_constraints.size()+this->_equations.size();
}

uint TaylorConstrainedImageSet::number_of_negative_constraints() const {
    return this->_constraints.size();
}

uint TaylorConstrainedImageSet::number_of_zero_constraints() const {
    return this->_equations.size();
}

ScalarTaylorFunction TaylorConstrainedImageSet::negative_constraint(uint i) const {
    return this->_constraints[i];
}

ScalarTaylorFunction TaylorConstrainedImageSet::zero_constraint(uint i) const {
    return this->_equations[i];
}




IntervalVector TaylorConstrainedImageSet::domain() const {
    return this->_domain;
}

IntervalVector TaylorConstrainedImageSet::reduced_domain() const {
    return this->_reduced_domain;
}

IntervalVector TaylorConstrainedImageSet::codomain() const {
    return Box(this->_function.range()).bounding_box();
}

VectorTaylorFunction const&  TaylorConstrainedImageSet::function() const {
    return this->_function;
}

VectorTaylorFunction TaylorConstrainedImageSet::taylor_function() const {
    return this->_function;
}

RealVectorFunction TaylorConstrainedImageSet::real_function() const {
    return RealVectorFunction(Vector< Polynomial<Real> >(this->_function.polynomial()));
}

uint TaylorConstrainedImageSet::dimension() const {
    return this->_function.result_size();
}

uint TaylorConstrainedImageSet::number_of_parameters() const {
    return this->_function.argument_size();
}

Box TaylorConstrainedImageSet::bounding_box() const {
    return Box(this->_function.codomain()).bounding_box();
}

Float TaylorConstrainedImageSet::radius() const {
    return this->bounding_box().radius();
}

Point TaylorConstrainedImageSet::centre() const {
    return this->bounding_box().centre();
}


tribool
TaylorConstrainedImageSet::satisfies(NonlinearConstraint c) const
{
    TaylorConstrainedImageSet copy=*this;
    copy.new_state_constraint(c);
    if(copy.empty()) { return false; }
    else { return indeterminate; }
}

tribool TaylorConstrainedImageSet::bounded() const
{
    return Box(this->domain()).bounded() || indeterminate;
}

tribool TaylorConstrainedImageSet::empty() const
{
    List<NonlinearConstraint> constraints=this->constraints();
    if(constraints.empty()) { return this->domain().empty(); }
    for(uint i=0; i!=constraints.size(); ++i) {
        if(Ariadne::disjoint(constraints[i].function().evaluate(this->_reduced_domain),constraints[i].bounds())) {
            return true;
        }
    }
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(this->_reduced_domain,constraints);
    for(uint i=0; i!=this->number_of_parameters(); ++i) {
        double l=this->_reduced_domain[i].lower().get_d();
        double u=this->_reduced_domain[i].upper().get_d();
        if(isnan(l) || isnan(u)) {
            ARIADNE_WARN("Reducing domain "<<_domain<<" yields "<<this->_reduced_domain);
            _reduced_domain[i]=_domain[i];
        }
    }
    if(this->_reduced_domain.empty()) { return true; }
    else { return indeterminate; }
}

tribool TaylorConstrainedImageSet::inside(const Box& bx) const
{
    return Ariadne::subset(this->_function.evaluate(this->_reduced_domain),bx);
}

tribool TaylorConstrainedImageSet::subset(const Box& bx) const
{
    List<NonlinearConstraint> constraints=this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(this->_reduced_domain,constraints);

    if(_reduced_domain.empty()) { return true; }

    for(uint i=0; i!=bx.dimension(); ++i) {
        const Box test_domain=this->_reduced_domain;
        constraints.append(ScalarTaylorFunction(this->_function[i]) <= bx[i].lower());
        if(possibly(contractor.feasible(test_domain,constraints).first)) { return indeterminate; }
        constraints.back()=(ScalarTaylorFunction(this->_function[i]) >= bx[i].upper());
        if(possibly(contractor.feasible(test_domain,constraints).first)) { return indeterminate; }
        constraints.pop_back();
    }
    return true;
}

tribool TaylorConstrainedImageSet::disjoint(const Box& bx) const
{
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"TaylorConstrainedImageSet::subset(Box): self="<<*this<<", box="<<bx);
    List<NonlinearConstraint> constraints=this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(this->_reduced_domain,constraints);

    if(_reduced_domain.empty()) { return true; }

    const Box test_domain=this->_reduced_domain;
    for(uint i=0; i!=bx.dimension(); ++i) {
        constraints.append(ScalarTaylorFunction(this->_function[i]) >= bx[i].lower());
        constraints.append(ScalarTaylorFunction(this->_function[i]) <= bx[i].upper());
    }
    return !contractor.feasible(test_domain,constraints).first;
}

void TaylorConstrainedImageSet::reduce() const
{
    List<NonlinearConstraint> constraints=this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(this->_reduced_domain,constraints);

    for(uint i=0; i!=this->number_of_parameters(); ++i) {
        double l=this->_reduced_domain[i].lower().get_d();
        double u=this->_reduced_domain[i].upper().get_d();
        if(isnan(l) || isnan(u)) {
            ARIADNE_WARN("Reducing domain "<<_domain<<" yields "<<this->_reduced_domain);
            _reduced_domain[i]=_domain[i];
        }
    }
}


Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split_zeroth_order() const
{
    Matrix<Interval> jacobian=this->function().jacobian(this->reduced_domain());

    // Compute the column of the matrix which has the norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    uint jmax=this->number_of_parameters();
    Float max_column_norm=0.0;
    for(uint j=0; j!=this->number_of_parameters(); ++j) {
        Float column_norm=0.0;
        for(uint i=0; i!=this->dimension(); ++i) {
            column_norm+=mag(jacobian[i][j]);
        }
        column_norm *= this->reduced_domain()[j].radius();
        if(column_norm>max_column_norm) {
            max_column_norm=column_norm;
            jmax=j;
        }
    }
    return this->split(jmax);
}


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

    // Compute the row of the nonlinearities array which has the highest norm
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


Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split_first_order() const
{
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(this->_function,this->_reduced_domain);

    // Compute the row of the nonlinearities array which has the highest norm
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

    if(jmax_in_row_imax==nonlinearities.column_size()) { ARIADNE_THROW(std::runtime_error, "split_first_order", "No need to split"); }

    //std::cerr<<"nonlinearities="<<nonlinearities<<"\nsplit_on_parameter="<<jmax_in_row_imax<<"\n\n\n";

    return this->split(jmax_in_row_imax);
}


Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split() const
{
    return this->split_zeroth_order();
}


Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split(uint d) const
{
    Vector<Interval> subdomain1,subdomain2;
    make_lpair(subdomain1,subdomain2)=Ariadne::split(this->_function.domain(),d);

    VectorTaylorFunction function1,function2;
    make_lpair(function1,function2)=Ariadne::split(this->_function,d);

    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
    result=make_pair(TaylorConstrainedImageSet(function1),TaylorConstrainedImageSet(function2));
    TaylorConstrainedImageSet& result1=result.first;
    TaylorConstrainedImageSet& result2=result.second;

    ScalarTaylorFunction constraint1,constraint2;
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        const ScalarTaylorFunction& constraint=*iter;
        make_lpair(constraint1,constraint2)=Ariadne::split(constraint,d);
        result1._constraints.append(constraint1);
        result2._constraints.append(constraint2);
    }

    ScalarTaylorFunction equation1,equation2;
    for(List<ScalarTaylorFunction>::const_iterator iter=this->_equations.begin();
        iter!=this->_equations.end(); ++iter)
    {
        const ScalarTaylorFunction& equation=*iter;
        make_lpair(equation1,equation1)=Ariadne::split(equation,d);
        result1._equations.append(equation1);
        result2._equations.append(equation1);
    }

    return make_pair(result1,result2);
}






RealScalarFunction make_function(const ScalarTaylorFunction& stf) {
    return RealScalarFunction(stf.polynomial())+Real(Interval(-stf.error(),+stf.error()));
}

void optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& r, const Box& d, const RealVectorFunction& fg, const Box& c, const GridCell& b, Point& x, Point& y, int e)
{
    // When making a new starting primal point, need to move components away from zero
    // This constant shows how far away from zero the points are
    static const double XSIGMA = 0.125;
    static const double TERR = -1.0/((1<<e)*1024.0);
    static const Float inf = Ariadne::inf<Float>();

    uint verbosity=0u;

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
        ScalarTaylorFunction xg=ScalarTaylorFunction::constant(d,0);
        Interval cnst=0.0;
        for(uint j=0; j!=n; ++j) {
            xg = xg - (x[j]-x[n+j])*ScalarTaylorFunction(d,fg[j]);
            cnst += (bx[j].upper()*x[j]-bx[j].lower()*x[n+j]);
        }
        for(uint i=0; i!=m; ++i) {
            xg = xg - (x[2*n+i]-x[2*n+m+i])*ScalarTaylorFunction::coordinate(d,i);
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

        //Pair<Box,Box> sd=solver.split(List<NonlinearConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=d.split();
        Point nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point ny = midpoint(sd.first);
        optimal_constraint_adjoin_outer_approximation_to(r, sd.first, fg, c, b, nx, ny, e);
        nx = (1.0-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        ny = midpoint(sd.second);
        optimal_constraint_adjoin_outer_approximation_to(r, sd.second, fg, c, b, x, ny, e);
    }

    if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        r.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell; t="<<t<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        Point sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        Point sy = y;
        optimal_constraint_adjoin_outer_approximation_to(r,d,fg,c,sb.first,sx,sy,e);
        sx = (1-XSIGMA)*x + Vector<Float>(x.size(),XSIGMA/x.size());
        sy = y;
        optimal_constraint_adjoin_outer_approximation_to(r,d,fg,c,sb.second,sx,sy,e);
    }


}


Float widths(const IntervalVector& bx) {
    Float res=0.0;
    for(uint i=0; i!=bx.size(); ++i) {
        res+=(bx[i].width());
    }
    return res;
}

// Adjoin an over-approximation to the solution of $f(D)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
void constraint_adjoin_outer_approximation_to(GridTreeSet& p, const Box& d, const RealVectorFunction& f, const RealVectorFunction& g, const Box& c, const GridCell& b, int e)
{
    uint verbosity=0u;

    const uint m=d.size();
    const uint nf=f.result_size();
    const uint ng=g.result_size();
    ARIADNE_LOG(2,"\nconstraint_adjoin_outer_approximation(...)\n");
    ARIADNE_LOG(2,"  dom="<<d<<" cnst="<<c<<" cell="<<b.box()<<" dpth="<<b.tree_depth()<<" e="<<e<<"\n");

    ConstraintSolver constraint_solver;

    if(subset(b,p)) {
        return;
    }

    // Try to prove disjointness
    const Box& old_domain=d;
    Box new_domain=old_domain;
    Box bx=b.box();

    ARIADNE_LOG(6,"  dom="<<old_domain<<"\n");
    for(uint i=0; i!=nf; ++i) {
        constraint_solver.hull_reduce(new_domain,f[i],bx[i]);
    }
    for(uint i=0; i!=ng; ++i) {
        constraint_solver.hull_reduce(new_domain,g[i],c[i]);
    }
    ARIADNE_LOG(6,"  dom="<<new_domain<<"\n");
    if(new_domain.empty()) {
        ARIADNE_LOG(4,"  Proved disjointness using hull reduce\n");
        return;
    }

    for(uint i=0; i!=nf; ++i) {
        constraint_solver.hull_reduce(new_domain,f[i],bx[i]);
    }
    for(uint i=0; i!=ng; ++i) {
        constraint_solver.hull_reduce(new_domain,g[i],c[i]);
    }
    ARIADNE_LOG(8,"  dom="<<new_domain<<"\n");
    if(new_domain.empty()) {
        ARIADNE_LOG(4,"  Proved disjointness using box reduce\n");
        return;
    }
    ARIADNE_LOG(6,"  dom="<<new_domain<<"\n");

    Box fd=f(new_domain);
    if(fd.disjoint(bx)) { return; }

    if(4*widths(fd)>widths(bx)) {
        //Pair<Box,Box> sd=solver.split(List<NonlinearConstraint>(1u,constraint),d);
        ARIADNE_LOG(4,"  Splitting domain\n");
        Pair<Box,Box> sd=new_domain.split();
        constraint_adjoin_outer_approximation_to(p, sd.first, f, g, c, b, e);
        constraint_adjoin_outer_approximation_to(p, sd.second, f, g, c, b, e);
    } else if(b.tree_depth()>=e*int(b.dimension())) {
        ARIADNE_LOG(4,"  Adjoining cell "<<b.box()<<"\n");
        p.adjoin(b);
    } else {
        ARIADNE_LOG(4,"  Splitting cell "<<b.box()<<"\n");
        Pair<GridCell,GridCell> sb = b.split();
        constraint_adjoin_outer_approximation_to(p,new_domain,f,g,c,sb.first,e);
        constraint_adjoin_outer_approximation_to(p,new_domain,f,g,c,sb.second,e);
    }


}


void adjoin_outer_approximation_to(GridTreeSet&, const Box& domain, const RealVectorFunction& function, const RealVectorFunction& negative_constraints, int depth);

void TaylorConstrainedImageSet::adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    //this->affine_adjoin_outer_approximation_to(paving,depth);
    this->constraint_adjoin_outer_approximation_to(paving,depth);
    //this->this->_adjoin_outer_approximation_to(paving,this->bounding_box(),depth);
}

void TaylorConstrainedImageSet::subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    Vector<Float> errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }
    this->_subdivision_adjoin_outer_approximation_to(paving,this->domain(),depth,errors);
}

void TaylorConstrainedImageSet::affine_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    this->affine_over_approximation().adjoin_outer_approximation_to(paving,depth);
}

void TaylorConstrainedImageSet::constraint_adjoin_outer_approximation_to(GridTreeSet& p, int e) const
{
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    const RealVectorFunction& f=this->real_function();

    RealVectorFunction g(this->_constraints.size()+this->_equations.size(),d.size());
    uint i=0;
    for(List<ScalarTaylorFunction>::const_iterator citer=this->_constraints.begin(); citer!=this->_constraints.end(); ++citer) {
        g.set(i,make_function(*citer));
        ++i;
    }
    for(List<ScalarTaylorFunction>::const_iterator eiter=this->_equations.begin(); eiter!=this->_equations.end(); ++eiter) {
        g.set(i,make_function(*eiter));
        ++i;
    }
    GridCell b=GridCell::smallest_enclosing_primary_cell(f(d),p.grid());
    IntervalVector cc(this->_constraints.size(),Interval(-inf<Float>(),0.0));
    IntervalVector ce(this->_equations.size(),Interval(0.0,0.0));
    Box c=intersection(Box(g(d)),Box(join(cc,ce)));

    Ariadne::constraint_adjoin_outer_approximation_to(p,d,f,g,c,b,e);
    p.recombine();
}

void TaylorConstrainedImageSet::optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& p, int e) const
{
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    RealVectorFunction f=this->real_function();

    RealVectorFunction g(this->_constraints.size(),d.size());
    uint i=0;
    for(List<ScalarTaylorFunction>::const_iterator citer=this->_constraints.begin(); citer!=this->_constraints.end(); ++citer) {
        g.set(i,make_function(*citer));
        ++i;
    }

    GridCell b=GridCell::smallest_enclosing_primary_cell(g(d),p.grid());
    Box c=intersection(g(d)+IntervalVector(g.result_size(),Interval(-1,1)),Box(g.result_size(),Interval(-inf<Float>(),0.0)));

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    RealVectorFunction fg=join(f,g);

    Ariadne::optimal_constraint_adjoin_outer_approximation_to(p,d,fg,c,b,x,y,e);

}

GridTreeSet TaylorConstrainedImageSet::outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->adjoin_outer_approximation_to(paving,depth);
    return paving;
}

GridTreeSet TaylorConstrainedImageSet::subdivision_outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->subdivision_adjoin_outer_approximation_to(paving,depth);
    return paving;
}

GridTreeSet TaylorConstrainedImageSet::affine_outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->affine_adjoin_outer_approximation_to(paving,depth);
    return paving;
}



TaylorConstrainedImageSet TaylorConstrainedImageSet::restriction(const Vector<Interval>& subdomain) const
{
    ARIADNE_ASSERT_MSG(subdomain.size()==this->number_of_parameters(),"set="<<*this<<", subdomain="<<subdomain);
    ARIADNE_ASSERT_MSG(Ariadne::subset(subdomain,this->domain()),"set.domain()="<<this->domain()<<", subdomain="<<subdomain);
    TaylorConstrainedImageSet result(*this);
    result._domain=subdomain;
    result._reduced_domain=Ariadne::intersection(static_cast<const Vector<Interval>&>(result._reduced_domain),subdomain);
    result._function=Ariadne::restrict(result._function,subdomain);
    ScalarTaylorFunction new_constraint;
    for(List<ScalarTaylorFunction>::iterator iter=result._constraints.begin();
        iter!=result._constraints.end(); ++iter)
    {
        ScalarTaylorFunction& constraint=*iter;
        constraint=Ariadne::restrict(constraint,subdomain);
    }
    for(List<ScalarTaylorFunction>::iterator iter=result._equations.begin();
        iter!=result._equations.end(); ++iter)
    {
        ScalarTaylorFunction& equation=*iter;
        equation=Ariadne::restrict(equation,subdomain);
    }
    return result;
}



void TaylorConstrainedImageSet::draw(CanvasInterface& canvas) const {
    switch(DRAWING_METHOD) {
        case BOX_DRAW:
            this->box_draw(canvas);
            break;
        case AFFINE_DRAW:
            //if(this->number_of_zero_constraints()!=0) { this->box_draw(canvas); }
            this->affine_draw(canvas,DRAWING_ACCURACY);
            break;
        case GRID_DRAW:
            this->grid_draw(canvas);
            break;
        default:
            ARIADNE_WARN("Unknown drawing method\n");
    }
}

void TaylorConstrainedImageSet::box_draw(CanvasInterface& canvas) const {
    this->reduce();
    Box(this->_function(this->_reduced_domain)).draw(canvas);
}

void TaylorConstrainedImageSet::affine_draw(CanvasInterface& canvas, uint accuracy) const {
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->_reduced_domain,this->_domain),*this);

    // Bound the maximum number of splittings allowed to draw a particular set.
    // Note that this gives rise to possibly 2^MAX_DEPTH split sets!!
    static const int MAXIMUM_DEPTH = 16;

    // The basic approximation error when plotting with accuracy=0
    static const double BASIC_ERROR = 0.0625;

    const double max_error=BASIC_ERROR/(1<<accuracy);

    VectorTaylorFunction fg(this->dimension()+this->number_of_constraints(),this->domain());
    for(uint i=0; i!=this->dimension(); ++i) { fg[i]=this->_function[i]; }
    for(uint i=0; i!=this->_constraints.size(); ++i) { fg[i+this->dimension()]=this->_constraints[i]; }
    for(uint i=0; i!=this->_equations.size(); ++i) { fg[i+this->dimension()+this->_constraints.size()]=this->_equations[i]; }

    List<Box> subdomains;
    List<Box> unsplitdomains;
    List<Box> splitdomains;
    unsplitdomains.append(this->_reduced_domain);
    Box splitdomain1,splitdomain2;
    for(int i=0; i!=MAXIMUM_DEPTH; ++i) {
        //std::cerr<<"i="<<i<<"\nsubdomains="<<subdomains<<"\nunsplitdomains="<<unsplitdomains<<"\n\n";
        for(uint n=0; n!=unsplitdomains.size(); ++n) {
            uint k; double err;
            make_lpair(k,err)=nonlinearity_index_and_error(fg,unsplitdomains[n]);
            //std::cerr<<"  domain="<<unsplitdomains[n]<<" k="<<k<<" err="<<err<<" max_err="<<max_error<<"\n";
            if(k==this->number_of_parameters() || err < max_error) {
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
        ARIADNE_WARN("Cannot obtain desired accuracy in drawing "<<*this<<" without excessive splitting.");
    }

    for(uint n=0; n!=subdomains.size(); ++n) {
        try {
            this->restriction(subdomains[n]).affine_over_approximation().draw(canvas);
        } catch(...) {
            this->restriction(subdomains[n]).box_draw(canvas);
        }
    }
};


void TaylorConstrainedImageSet::grid_draw(CanvasInterface& canvas, uint accuracy) const {
    // TODO: Project to grid first
    this->outer_approximation(Grid(this->dimension()),accuracy).draw(canvas);
}

Map<List<DiscreteEvent>,RealScalarFunction> pretty(const Map<List<DiscreteEvent>,ScalarTaylorFunction>& constraints) {
    Map<List<DiscreteEvent>,RealScalarFunction> result;
    for(Map<List<DiscreteEvent>,ScalarTaylorFunction>::const_iterator iter=constraints.begin();
    iter!=constraints.end(); ++iter) {
        result.insert(iter->first,iter->second.real_function());
    }
    return result;
}



template<class K, class V> Map<K,V> filter(const Map<K,V>& m, const Set<K>& s) {
    Map<K,V> r;
    for(typename Set<K>::const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        r.insert(*m.find(*iter));
    }
    return r;
}

template<class T> struct Repr { const T& ref; Repr(const T& r) : ref(r) { } };
template<class T> Repr<T> repr(const T& t) { return Repr<T>(t); }
template<class T> std::ostream& operator<<(std::ostream& os, const Repr<T>& repr) { repr.ref.repr(os); return os; }

template<class T> std::ostream& operator<<(std::ostream& os, const Repr< List<T> >& repr) {
    const List<T>& lst=repr.ref; os << "["; for(uint i=0; i!=lst.size(); ++i) { if(i!=0) { os << ","; } lst[i].repr(os); } os << "]"; return os; }

std::ostream& TaylorConstrainedImageSet::write(std::ostream& os) const {
    const bool LONG_FORMAT=false;

    if(LONG_FORMAT) {
        os << "TaylorConstrainedImageSet"
           << "(\n  domain=" << this->domain()
           << ",\n  range=" << this->bounding_box()
           << ",\n  function=" << this->taylor_function()
           << ",\n  negative_constraints=" << this->_constraints
           << ",\n  zero_constraints=" << this->_equations
           << "\n)\n";
    } else {
        os << "TaylorConstrainedImageSet"
           << "( domain=" << this->domain()
           << ", range=" << this->bounding_box()
           << ", function=" << repr(this->taylor_function())
           << ", negative_constraints=" << repr(this->_constraints)
           << ", zero_constraints=" << repr(this->_equations)
           << ")";

    } return os;
}



void TaylorConstrainedImageSet::
_subdivision_adjoin_outer_approximation_to(GridTreeSet& gts, const IntervalVector& subdomain,
                                           uint depth, const FloatVector& errors) const
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;

    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    for(const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0) {
            return;
        }
    }
    for(const_iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0 || constraint_range.upper() < 0.0 ) {
            return;
        }
    }

    Box range=evaluate(this->_function,subdomain);
    bool small=true;
    for(uint i=0; i!=range.size(); ++i) {
        if(range[i].radius()>errors[i]*RELATIVE_SMALLNESS) {
            small=false;
            break;
        }
    }

    if(small) {
        gts.adjoin_outer_approximation(range,depth);
    } else {
        Vector<Interval> subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=Ariadne::split(subdomain);
        this->_subdivision_adjoin_outer_approximation_to(gts,subdomain1,depth,errors);
        this->_subdivision_adjoin_outer_approximation_to(gts,subdomain2,depth,errors);
    }
}


AffineSet
TaylorConstrainedImageSet::affine_approximation() const
{
    this->_check();
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();

    TaylorConstrainedImageSet set(*this);

    //if(set._equations.size()>0) { set._solve_zero_constraints(); }
    this->_check();

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np);
    for(uint i=0; i!=nx; ++i) {
        ScalarTaylorFunction component=set._function[i];
        h[i]=component.model().value();
        for(uint j=0; j!=np; ++j) {
            G[i][j]=component.model().gradient(j);
        }
    }
    AffineSet result(G,h);

    Vector<Float> a(np);
    Float b;

    for(const_iterator iter=set._constraints.begin();
            iter!=set._constraints.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=set._equations.begin();
            iter!=set._equations.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_equality_constraint(a,b);
    }
    return result;
}

struct IntervalAffineModel {
    Float _c; Vector<Float> _g; Float _e;
    IntervalAffineModel(Float c, const Vector<Float>& g, Float e) : _c(c), _g(g), _e(e) { }
};

IntervalAffineModel _affine_model(const IntervalTaylorModel& tm) {
    IntervalAffineModel result(0.0,Vector<Float>(tm.argument_size(),0.0),tm.error());
    set_rounding_upward();
    for(IntervalTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()>=2) { result._e+=abs(iter->data()); }
        else if(iter->key().degree()==0) {result. _c=iter->data(); }
        else {
            for(uint j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]!=0) { result._g[j]=iter->data(); break; }
            }
        }
    }
    set_rounding_to_nearest();
    return result;
}

AffineSet
TaylorConstrainedImageSet::affine_over_approximation() const
{
    this->_check();
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint nc=this->_constraints.size();
    const uint neq=this->_equations.size();
    const uint np=this->number_of_parameters();

    TaylorConstrainedImageSet set(*this);
    for(uint i=0; i!=nx; ++i) {
        const_cast<IntervalTaylorModel&>(set._function.models()[i]).truncate(1u);
    }
    for(uint i=0; i!=nc; ++i) {
        const_cast<IntervalTaylorModel&>(set._constraints[i].model()).truncate(1u);
    }
    for(uint i=0; i!=neq; ++i) {
        const_cast<IntervalTaylorModel&>(set._equations[i].model()).truncate(1u);
        // Code below introduces artificial error into equality constraints to make them inequality constraints
        //const_cast<IntervalTaylorModel&>(set._equations[i].model()).error()+=std::numeric_limits<float>::epsilon();
    }

    // Compute the number of values with a nonzero error
    uint nerr=0;
    for(uint i=0; i!=nx; ++i) { if(set.function()[i].error()>0.0) { ++nerr; } }

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np+nerr);
    uint ierr=0; // The index where the error bound should go
    for(uint i=0; i!=nx; ++i) {
        ScalarTaylorFunction component=set._function[i];
        h[i]=component.model().value();
        for(uint j=0; j!=np; ++j) {
            G[i][j]=component.model().gradient(j);
        }
        if(component.model().error()>0.0) {
            G[i][np+ierr]=component.model().error();
            ++ierr;
        }
    }

    AffineSet result(G,h);

    Vector<Float> a(np+nerr, 0.0);
    Float b;

    for(const_iterator iter=set._constraints.begin();
            iter!=set._constraints.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        b=sub_up(constraint.model().error(),constraint.model().value());
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=set._equations.begin();
            iter!=set._equations.end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        if(constraint.model().error()==0.0) {
            b=-constraint.model().value();
            result.new_equality_constraint(a,b);
        } else {
            b=sub_up(constraint.model().error(),constraint.model().value());
            result.new_inequality_constraint(a,b);
            b=add_up(constraint.model().error(),constraint.model().value());
            result.new_inequality_constraint(-a,b);
        }
    }

    ARIADNE_LOG(2,"set="<<*this<<"\nset.affine_over_approximation()="<<result<<"\n");
    return result;
}

TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set, const Interval& ivl) {
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    VectorTaylorFunction new_function=combine(set.taylor_function(),ScalarTaylorFunction::identity(ivl));

    TaylorConstrainedImageSet result(new_function);
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(embed(*iter,ivl));
    }
    for(const_iterator iter=set._equations.begin(); iter!=set._equations.end(); ++iter) {
        result._equations.append(embed(*iter,ivl));
    }

    return result;
}

TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set, const Box& bx) {
    return product(set,TaylorConstrainedImageSet(bx));
}

TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set1, const TaylorConstrainedImageSet& set2) {
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    VectorTaylorFunction new_function=combine(set1.taylor_function(),set2.taylor_function());

    TaylorConstrainedImageSet result(new_function);
    for(const_iterator iter=set1._constraints.begin(); iter!=set1._constraints.end(); ++iter) {
        result._constraints.append(embed(*iter,set2.domain()));
    }
    for(const_iterator iter=set2._constraints.begin(); iter!=set2._constraints.end(); ++iter) {
        result._constraints.append(embed(set1.domain(),*iter));
    }
    for(const_iterator iter=set1._equations.begin(); iter!=set1._equations.end(); ++iter) {
        result._equations.append(embed(*iter,set2.domain()));
    }
    for(const_iterator iter=set2._equations.begin(); iter!=set2._equations.end(); ++iter) {
        result._equations.append(embed(set1.domain(),*iter));
    }

    return result;
}

TaylorConstrainedImageSet apply(const IntervalVectorFunctionInterface& function, const TaylorConstrainedImageSet& set) {
    TaylorConstrainedImageSet result(set);
    result.apply_map(function);
    return result;
}

} // namespace Ariadne


