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
#include "procedure.h"

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
#include <boost/concept_check.hpp>
#include <include/operators.h>


namespace Ariadne {

static const uint verbosity = 0u;

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

void constraint_adjoin_outer_approximation(GridTreeSet&, const Box& domain, const IntervalVectorFunction& function, const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);


Interval make_domain(const RealInterval& ivl) {
    rounding_mode_t rnd=get_rounding_mode();
    Interval dom_lower_ivl=Interval(ivl.lower());
    Interval dom_upper_ivl=Interval(ivl.upper());
    Float dom_lower=dom_lower_ivl.lower();
    Float dom_upper=dom_upper_ivl.upper();
    set_rounding_downward();
    float flt_dom_lower=numeric_cast<double>(dom_lower);
    while(double(flt_dom_lower)>dom_lower) {
        flt_dom_lower-=std::numeric_limits<float>::min();
    }
    dom_lower=flt_dom_lower;
    set_rounding_upward();
    float flt_dom_upper=numeric_cast<double>(dom_upper);
    while(double(flt_dom_upper)<dom_upper) {
        flt_dom_upper+=std::numeric_limits<float>::min();
    }
    dom_upper=flt_dom_upper;
    set_rounding_mode(rnd);
    return Interval(dom_lower,dom_upper);
}


VectorTaylorFunction make_identity(const RealBox& bx, const Sweeper& swp) {
    IntervalVector dom(bx.dimension());
    FloatVector errs(bx.dimension());

    for(uint i=0; i!=bx.dimension(); ++i) {
        Interval dom_lower_ivl=numeric_cast<Interval>(bx[i].lower());
        Interval dom_upper_ivl=numeric_cast<Interval>(bx[i].upper());
        // Convert to single-precision values
        Float dom_lower_flt=numeric_cast<float>(bx[i].lower());
        Float dom_upper_flt=numeric_cast<float>(bx[i].upper());
        set_rounding_upward();
        Float err=max( max(dom_upper_ivl.upper()-dom_upper_flt,dom_upper_flt-dom_upper_ivl.lower()),
                       max(dom_lower_ivl.upper()-dom_lower_flt,dom_lower_flt-dom_lower_ivl.lower()) );
        set_rounding_to_nearest();
        dom[i]=Interval(dom_lower_flt,dom_upper_flt);
        errs[i]=err;
    }

    VectorTaylorFunction res=VectorTaylorFunction::identity(dom,swp);
    for(uint i=0; i!=bx.dimension(); ++i) {
        res[i]=res[i]+Interval(-errs[i],+errs[i]);
    }

    return res;
};


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

Sweeper TaylorConstrainedImageSet::sweeper() const {
    return this->_function.sweeper();
}

TaylorFunctionFactory TaylorConstrainedImageSet::function_factory() const {
    return TaylorFunctionFactory(this->sweeper());
}

// FIXME: What if solving for constraint leaves domain?
void TaylorConstrainedImageSet::_solve_zero_constraints() {
    this->_check();
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ) {
        const Vector<Interval>& domain=this->domain();
        const IntervalTaylorModel& model=iter->model();
        const uint k=model.argument_size()-1u;
        IntervalTaylorModel zeroth_order(k,this->sweeper());
        IntervalTaylorModel first_order(k,this->sweeper());
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
    : _domain(), _function(), _reduced_domain(), _is_fully_reduced(true)
{
}

TaylorConstrainedImageSet* TaylorConstrainedImageSet::clone() const
{
    return new TaylorConstrainedImageSet(*this);
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const Box& box, Sweeper sweeper)
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


    this->_function=VectorTaylorFunction(box.dimension(),this->_domain,sweeper);
    uint j=0;
    proper_coordinates.append(box.dimension());
    for(uint i=0; i!=box.dimension(); ++i) {
        if(proper_coordinates[j]==i) {
            this->_function[i]=ScalarTaylorFunction::coordinate(this->_domain,j,sweeper);
            ++j;
        } else {
            this->_function[i]=ScalarTaylorFunction::constant(this->_domain,box[i],sweeper);
        }
    }
    this->_reduced_domain=this->_domain;
    this->_is_fully_reduced=true;
}


TaylorConstrainedImageSet::TaylorConstrainedImageSet(const Box& box, const TaylorFunctionFactory& factory)
{
    *this = TaylorConstrainedImageSet(box,factory.sweeper());
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const IntervalVector& domain, const IntervalVectorFunction& function, Sweeper sweeper)
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    this->_domain=domain;
    this->_function=VectorTaylorFunction(this->_domain,function,sweeper);
    this->_reduced_domain=this->_domain;
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const IntervalVector& domain, const IntervalVectorFunction& function, const List<IntervalNonlinearConstraint>& constraints, Sweeper sweeper)
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    const double min=std::numeric_limits<double>::min();
    this->_domain=domain;
    for(uint i=0; i!=this->_domain.size(); ++i) {
        if(this->_domain[i].radius()==0) {
            this->_domain[i]+=Interval(-min,+min);
        }
    }

    this->_function=VectorTaylorFunction(this->_domain,function,sweeper);

    for(uint i=0; i!=constraints.size(); ++i) {
        ARIADNE_ASSERT_MSG(domain.size()==constraints[i].function().argument_size(),"domain="<<domain<<", constraint="<<constraints[i]);
        if(constraints[i].bounds().singleton()) {
            this->new_equality_constraint(constraints[i].function()-Interval(constraints[i].bounds().midpoint()));
        } else {
            if(constraints[i].bounds().lower()>-inf<Float>()) {
                this->new_negative_constraint(Interval(constraints[i].bounds().lower())-constraints[i].function());
            }
            if(constraints[i].bounds().upper()<+inf<Float>()) {
                this->new_negative_constraint(constraints[i].function()-Interval(constraints[i].bounds().upper()));
            }
        }
    }

    this->_reduced_domain=domain;
    this->reduce();

}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const IntervalVector& domain, const IntervalVectorFunction& function, const IntervalNonlinearConstraint& constraint, Sweeper sweeper)
{
    *this=TaylorConstrainedImageSet(domain,function,make_list(constraint),sweeper);
}

TaylorConstrainedImageSet::TaylorConstrainedImageSet(const VectorTaylorFunction& taylor_function)
    : _domain(taylor_function.domain()), _function(taylor_function), _reduced_domain(_domain), _is_fully_reduced(true)
{
}





// Returns true if the entire set is positive; false if entire set is negative
tribool TaylorConstrainedImageSet::satisfies(IntervalScalarFunction constraint) const
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

void TaylorConstrainedImageSet::new_parameter(Interval ivl)
{
    this->_domain=join(this->_domain,ivl);
    this->_reduced_domain=join(this->_reduced_domain,ivl);
    this->_function=embed(this->_function,ivl);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        this->_constraints[i]=embed(this->_constraints[i],ivl);
    }
    for(uint i=0; i!=this->_equations.size(); ++i) {
        this->_equations[i]=embed(this->_equations[i],ivl);
    }
    this->_check();
}

void TaylorConstrainedImageSet::new_variable(Interval ivl)
{
    ScalarTaylorFunction variable_function = ScalarTaylorFunction::identity(ivl,this->sweeper());
    this->_domain=join(this->_domain,ivl);
    this->_reduced_domain=join(this->_reduced_domain,ivl);
    this->_function=combine(this->_function,variable_function);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        this->_constraints[i]=embed(this->_constraints[i],ivl);
    }
    for(uint i=0; i!=this->_equations.size(); ++i) {
        this->_equations[i]=embed(this->_equations[i],ivl);
    }
    this->_check();
}

void TaylorConstrainedImageSet::apply_map(IntervalVectorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    VectorTaylorFunction& function=this->_function;
    function=compose(map,function);
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow(IntervalVectorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,combine(this->_function,VectorTaylorFunction::identity(Vector<Interval>(1u,time),this->sweeper())));
    for(List<ScalarTaylorFunction>::iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ScalarTaylorFunction>::iterator iter=this->_equations.begin(); iter!=this->_equations.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}

void TaylorConstrainedImageSet::apply_flow_step(IntervalVectorFunction flow, Float time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_function=compose(flow,join(this->_function,ScalarTaylorFunction::constant(this->_function.domain(),time,this->sweeper())));
    this->_check();
}

void TaylorConstrainedImageSet::apply_state_flow_step(IntervalVectorFunction flow, IntervalScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", time="<<time);
    this->_function=compose(flow,join(this->_function,compose(time,this->_function)));
    this->_check();
}

void TaylorConstrainedImageSet::apply_parameter_flow_step(IntervalVectorFunction flow, IntervalScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->number_of_parameters(),"number_of_parameters="<<this->number_of_parameters()<<", time="<<time);
    this->_function=compose(flow,join(this->_function,ScalarTaylorFunction(this->_function.domain(),time,this->sweeper())));
    this->_check();
}

void TaylorConstrainedImageSet::new_state_constraint(IntervalNonlinearConstraint constraint) {
    this->_is_fully_reduced=false;
    Float infty=+inf<Float>();
    Interval interval=constraint.bounds();
    ScalarTaylorFunction composed_function=compose(constraint.function(),this->_function);
    if(interval.lower()==0.0 && interval.upper()==0.0) {
        this->new_zero_constraint(composed_function);
    } else if(interval.lower()==0.0 && interval.upper()==infty) {
        this->new_negative_constraint(-composed_function);
    } else if(interval.lower()==-infty && interval.upper()==0.0) {
        this->new_negative_constraint(composed_function);
    } else {
        ARIADNE_FAIL_MSG("TaylorConstrainedImageSet cannot currently handle constraints which are not of the form g(x) <=> 0");
    }
}

void TaylorConstrainedImageSet::new_parameter_constraint(IntervalNonlinearConstraint constraint) {
    ARIADNE_ASSERT(constraint.function().argument_size()==this->number_of_parameters());
    this->_is_fully_reduced=false;
    Float infty=+inf<Float>();
    Interval interval=constraint.bounds();
    if(interval.lower()==0.0 && interval.upper()==0.0) {
        this->new_zero_constraint(constraint.function());
    } else if(interval.lower()==0.0 && interval.upper()==infty) {
        this->new_negative_constraint(-constraint.function());
    } else if(interval.lower()==-infty && interval.upper()==0.0) {
        this->new_negative_constraint(constraint.function());
    } else {
        ARIADNE_FAIL_MSG("TaylorConstrainedImageSet cannot currently handle constraints which are not of the form g(x) <=> 0");
    }
}


void TaylorConstrainedImageSet::new_negative_constraint(IntervalScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_is_fully_reduced=false;
    this->_constraints.append(ScalarTaylorFunction(this->domain(),constraint,this->sweeper()));
}

void TaylorConstrainedImageSet::new_equality_constraint(IntervalScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_is_fully_reduced=false;
    this->_equations.append(ScalarTaylorFunction(this->domain(),constraint,this->sweeper()));
}

void TaylorConstrainedImageSet::new_zero_constraint(IntervalScalarFunction constraint) {
    ARIADNE_ASSERT_MSG(constraint.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint="<<constraint);
    this->_is_fully_reduced=false;
    this->_equations.append(ScalarTaylorFunction(this->domain(),constraint,this->sweeper()));
}


List<ScalarTaylorFunction> const&
TaylorConstrainedImageSet::negative_constraints() const {
    return this->_constraints;
}

List<ScalarTaylorFunction> const&
TaylorConstrainedImageSet::zero_constraints() const {
    return this->_equations;
}

List<IntervalNonlinearConstraint>
TaylorConstrainedImageSet::constraints() const {
    List<IntervalNonlinearConstraint> result;
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
TaylorConstrainedImageSet::satisfies(IntervalNonlinearConstraint c) const
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
    if(definitely(Ariadne::empty(this->_reduced_domain))) { return true; }
    if(this->_constraints.empty() && this->_equations.empty()) { return Ariadne::empty(this->domain()); }
    if(!this->_is_fully_reduced) { this->reduce(); this->reduce(); this->reduce(); }

    for(uint i=0; i!=this->_constraints.size(); ++i) {
        if(this->_constraints[i](this->_reduced_domain).lower()>0.0) {
            this->_reduced_domain[0] = Interval(1,-1);
            return true;
        }
    }
    for(uint i=0; i!=this->_equations.size(); ++i) {
        if(!contains(this->_equations[i](this->_reduced_domain),0.0)) {
            this->_reduced_domain[0] = Interval(1,-1);
            return true;
        }
    }
    if(Ariadne::empty(this->_reduced_domain)) { return true; }
    return indeterminate;
}

tribool TaylorConstrainedImageSet::inside(const Box& bx) const
{
    return Ariadne::subset(this->_function.evaluate(this->_reduced_domain),bx);
}

tribool TaylorConstrainedImageSet::subset(const Box& bx) const
{
    this->reduce();

    return Ariadne::subset(this->_function.evaluate(this->_reduced_domain),bx) || indeterminate;

}

tribool TaylorConstrainedImageSet::disjoint(const Box& bx) const
{
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"TaylorConstrainedImageSet::subset(Box): self="<<*this<<", box="<<bx);
    List<IntervalNonlinearConstraint> constraints=this->constraints();
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
    List<IntervalNonlinearConstraint> constraints=this->constraints();
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

/*
    // Remove redundant constraints
    uint j=0;
    List<ScalarTaylorFunction>& mutable_constraints=const_cast<List<ScalarTaylorFunction>&>(this->_constraints);
    for(uint i=0; i!=mutable_constraints.size(); ++i) {
        if(mutable_constraints[i](this->_reduced_domain).upper()<0.0) { redundant_constraints.append(i); }
        else { if(i>j) { mutable_constraints[j]=mutable_constraints[j]; } ++j; }
    }
    mutable_constraints.resize(j);
*/


}


Pair<uint,double> lipschitz_index_and_error(const IntervalVectorFunctionInterface& function, const IntervalVector& domain)
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

Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split_zeroth_order() const
{
    return this->split(this->splitting_index_zeroth_order());
}


List<IntervalVector>
TaylorConstrainedImageSet::splitting_subdomains_zeroth_order() const
{
    List<IntervalVector> result;
    uint k=this->splitting_index_zeroth_order();
    if(k==this->number_of_parameters()) {
        result.append(this->_reduced_domain);
    } else {
        Pair<IntervalVector,IntervalVector> subdomains = this->_reduced_domain.split(this->splitting_index_zeroth_order());
        result.append(subdomains.first);
        result.append(subdomains.second);
    }
    return result;
}


uint
TaylorConstrainedImageSet::splitting_index_zeroth_order() const
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

    return jmax;
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


Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
TaylorConstrainedImageSet::split_first_order() const
{
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(this->_function,this->_reduced_domain);

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
    ARIADNE_PRECONDITION(d<this->number_of_parameters());
    Vector<Interval> subdomain1,subdomain2;
    make_lpair(subdomain1,subdomain2)=Ariadne::split(this->_function.domain(),d);

    VectorTaylorFunction function1,function2;
    make_lpair(function1,function2)=Ariadne::split(this->_function,d);

    Pair<TaylorConstrainedImageSet,TaylorConstrainedImageSet>
    result=make_pair(TaylorConstrainedImageSet(function1.domain(),function1,function1.sweeper()),
                     TaylorConstrainedImageSet(function2.domain(),function2,function2.sweeper()));
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

void optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& r, const Box& d, const VectorTaylorFunction& fg, const Box& c, const GridCell& b, Point& x, Point& y, int e)
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

static uint COUNT_TESTS=0u;
double IMAGE_MULTIPLE_OF_CELL = 1;

typedef Procedure<Interval> IntervalProcedure;

// Adjoin an over-approximation to the solution of $f(dom)$ such that $g(D) in C$ to the paving p, looking only at solutions in b.
void constraint_adjoin_outer_approximation_to(GridTreeSet& paving, const Box& domain, const VectorTaylorFunction& f, const VectorTaylorFunction& g, const Box& codomain, const GridCell& cell, int max_dpth, uint splt, const List<IntervalProcedure>& procedures)
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

    if(subset(cell,paving)) {
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
        constraint_adjoin_outer_approximation_to(paving, sd.first, f, g, codomain, cell, max_dpth, splt+1, procedures);
        constraint_adjoin_outer_approximation_to(paving, sd.second, f, g, codomain, cell, max_dpth, splt+1, procedures);
    } else if(cell.tree_depth()>=max_dpth) {
        ARIADNE_LOG(4,"  Adjoining cell "<<cell_box<<"\n");
        paving.adjoin(cell);
    } else {
        ARIADNE_LOG(4,"  Splitting cell "<<cell_box<<"\n");
        Pair<GridCell,GridCell> sb = cell.split();
        constraint_adjoin_outer_approximation_to(paving,new_domain,f,g,codomain,sb.first, max_dpth, splt, procedures);
        constraint_adjoin_outer_approximation_to(paving,new_domain,f,g,codomain,sb.second, max_dpth, splt, procedures);
    }


}


void adjoin_outer_approximation(GridTreeSet&, const Box& domain, const IntervalVectorFunction& function, const IntervalVectorFunction& negative_constraints, const IntervalVectorFunction& equality_constraints, int depth);

void TaylorConstrainedImageSet::adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    switch(DISCRETISATION_METHOD) {
        case SUBDIVISION_DISCRETISE:
            this->subdivision_adjoin_outer_approximation_to(paving,depth);
            break;
        case AFFINE_DISCRETISE:
            this->affine_adjoin_outer_approximation_to(paving,depth);
            break;
        case CONSTRAINT_DISCRETISE:
            this->constraint_adjoin_outer_approximation_to(paving,depth);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown discretisation method\n");
    }
}

void
_subdivision_adjoin_outer_approximation_to(GridTreeSet& gts, const TaylorConstrainedImageSet& set, const IntervalVector& subdomain,
                                           uint depth, const FloatVector& errors, uint splittings);

void TaylorConstrainedImageSet::subdivision_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    Vector<Float> errors(paving.dimension());
    for(uint i=0; i!=errors.size(); ++i) {
        errors[i]=paving.grid().lengths()[i]/(1<<depth);
    }
    _subdivision_adjoin_outer_approximation_to(paving,*this,this->domain(),depth,errors,0u);
}

IntervalProcedure make_procedure(const IntervalScalarFunctionInterface& f);

void TaylorConstrainedImageSet::constraint_adjoin_outer_approximation_to(GridTreeSet& p, int acc) const
{
    ARIADNE_LOG(6,"TaylorConstrainedImageSet::constraint_adjoin_outer_approximation_to(paving, depth)\n");
    ARIADNE_LOG(7,"  set="<<*this<<", grid="<<p.grid()<<", depth="<<acc<<"\n");
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    const VectorTaylorFunction& f=this->function();

    VectorTaylorFunction g(this->_constraints.size()+this->_equations.size(),d,this->sweeper());
    uint i=0;
    for(List<ScalarTaylorFunction>::const_iterator citer=this->_constraints.begin(); citer!=this->_constraints.end(); ++citer) {
        g.set(i,*citer);
        ++i;
    }
    for(List<ScalarTaylorFunction>::const_iterator eiter=this->_equations.begin(); eiter!=this->_equations.end(); ++eiter) {
        g.set(i,*eiter);
        ++i;
    }
    GridCell b=GridCell::smallest_enclosing_primary_cell(f(d),p.grid());
    IntervalVector cc(this->_constraints.size(),Interval(-inf<Float>(),0.0));
    IntervalVector ce(this->_equations.size(),Interval(0.0,0.0));
    Box c=intersection(Box(g(d)),Box(join(cc,ce)));

    List<IntervalProcedure> procedures;
    procedures.reserve(this->dimension()+this->_constraints.size()+this->_equations.size());
    for(uint i=0; i!=this->dimension(); ++i) { procedures.append(make_procedure(this->_function[i])); }
    for(uint i=0; i!=this->_constraints.size(); ++i) { procedures.append(make_procedure(this->_constraints[i])); }
    for(uint i=0; i!=this->_equations.size(); ++i) { procedures.append(make_procedure(this->_equations[i])); }

    COUNT_TESTS=0;
    Ariadne::constraint_adjoin_outer_approximation_to(p,d,f,g,c,b,acc*p.dimension(),0, procedures);
    //std::cerr<<"Computing outer approximation considered a total of "<<COUNT_TESTS<<" domains/cells\n";
    //std::cerr<<"Measure of paving is "<<p.measure()<<"\n";
    p.recombine();
}

void TaylorConstrainedImageSet::optimal_constraint_adjoin_outer_approximation_to(GridTreeSet& p, int e) const
{
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    VectorTaylorFunction f=this->function();

    VectorTaylorFunction g(this->_constraints.size(),d,this->sweeper());
    uint i=0;
    for(List<ScalarTaylorFunction>::const_iterator citer=this->_constraints.begin(); citer!=this->_constraints.end(); ++citer) {
        g.set(i,*citer);
        ++i;
    }

    GridCell b=GridCell::smallest_enclosing_primary_cell(g(d),p.grid());
    Box c=intersection(g(d)+IntervalVector(g.result_size(),Interval(-1,1)),Box(g.result_size(),Interval(-inf<Float>(),0.0)));

    Point y=midpoint(d);
    const uint l=(d.size()+f.result_size()+g.result_size())*2;
    Point x(l); for(uint k=0; k!=l; ++k) { x[k]=1.0/l; }

    VectorTaylorFunction fg=join(f,g);

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

void TaylorConstrainedImageSet::affine_adjoin_outer_approximation_to(GridTreeSet& paving, int depth) const
{
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->_reduced_domain,this->_domain),*this);

    // Bound the maximum number of splittings allowed to draw a particular set.
    // Note that this gives rise to possibly 2^MAX_DEPTH split sets!!
    static const int MAXIMUM_DEPTH = 16;

    // The basic approximation error when plotting with accuracy=0
    static const double BASIC_ERROR = 0.0625;

    const double max_error=BASIC_ERROR/(1<<depth);

    VectorTaylorFunction fg(this->dimension()+this->number_of_constraints(),this->domain(),this->sweeper());
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
        ARIADNE_WARN("Cannot obtain desired accuracy in outer approximation "<<*this<<" without excessive splitting.");
    }

    for(uint n=0; n!=subdomains.size(); ++n) {
        this->restriction(subdomains[n]).affine_over_approximation().adjoin_outer_approximation_to(paving,depth);
    }
}





void TaylorConstrainedImageSet::
recondition()
{
    const double MAXIMUM_ERROR = std::numeric_limits<double>::epsilon() * 1024;

    List<uint> large_error_indices;

    for(uint i=0; i!=this->_function.result_size(); ++i) {
        Float error=this->_function[i].model().error();
        if(error > MAXIMUM_ERROR) {
            large_error_indices.append(i);
        }
    }

    IntervalVector error_domains(large_error_indices.size());
    for(uint i=0; i!=large_error_indices.size(); ++i) {
        Float error=this->_function[large_error_indices[i]].model().error();
        error_domains[i]=Interval(-error,+error);
    }
    error_domains=IntervalVector(large_error_indices.size(),Interval(-1,+1));
    uint k=this->number_of_parameters();

    this->_domain=join(this->_domain,error_domains);
    this->_reduced_domain=join(this->_reduced_domain,error_domains);
    this->_function=embed(this->_function,error_domains);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        this->_constraints[i]=embed(this->_constraints[i],error_domains);
    }
    for(uint i=0; i!=this->_equations.size(); ++i) {
        this->_equations[i]=embed(this->_equations[i],error_domains);
    }

    for(uint i=0; i!=large_error_indices.size(); ++i) {
        Float error=this->_function[large_error_indices[i]].model().error();
        if(error > MAXIMUM_ERROR) {
            this->_function[i].set_error(0.0);
            this->_function[i] = this->_function[i] + ScalarTaylorFunction::coordinate(this->_domain,k,this->sweeper())*error;
            ++k;
        }
    }
}





void TaylorConstrainedImageSet::restrict(const Vector<Interval>& subdomain)
{
    ARIADNE_ASSERT_MSG(subdomain.size()==this->number_of_parameters(),"set="<<*this<<", subdomain="<<subdomain);
    ARIADNE_ASSERT_MSG(Ariadne::subset(subdomain,this->domain()),"set.domain()="<<this->domain()<<", subdomain="<<subdomain);
    TaylorConstrainedImageSet& result(*this);
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
    this->reduce();
}

TaylorConstrainedImageSet TaylorConstrainedImageSet::restriction(const Vector<Interval>& subdomain) const
{
    TaylorConstrainedImageSet result(*this);
    result.restrict(subdomain);
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

    // If the reduced domain is empty, then the set is empty; abort
    if(this->_reduced_domain.empty()) {
        return;
    }

    VectorTaylorFunction fg(this->dimension()+this->number_of_constraints(),this->domain(),this->sweeper());
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

template<class T> std::ostream& operator<<(std::ostream& os, const Representation< List<T> >& repr) {
    const List<T>& lst=*repr.pointer; os << "["; for(uint i=0; i!=lst.size(); ++i) { if(i!=0) { os << ","; } lst[i].repr(os); } os << "]"; return os; }

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



void
_subdivision_adjoin_outer_approximation_to(GridTreeSet& gts, const TaylorConstrainedImageSet& set, const IntervalVector& subdomain,
                                           uint depth, const FloatVector& errors, uint splittings)
{
    // How small an over-approximating box needs to be relative to the cell size
    static const double RELATIVE_SMALLNESS=0.5;
    static const uint MAXIMUM_SPLITTINGS = 18;

    //std::cerr<<"subdomain="<<subdomain<<"; "<<"depth="<<depth<<" splittings="<<splittings<<"\n";
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    for(const_iterator iter=set.negative_constraints().begin(); iter!=set.negative_constraints().end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0) {
            return;
        }
    }
    for(const_iterator iter=set.zero_constraints().begin(); iter!=set.zero_constraints().end(); ++iter) {
        const ScalarTaylorFunction& constraint=*iter;
        Interval constraint_range=constraint.evaluate(subdomain);
        if(constraint_range.lower() > 0.0 || constraint_range.upper() < 0.0 ) {
            return;
        }
    }

    Box range=evaluate(set.function(),subdomain);

    Array<Float> radii(set.dimension());
    for(uint i=0; i!=radii.size(); ++i) { radii[i]=range[i].radius()-set.function()[i].error(); }
    bool small=true;
    for(uint i=0; i!=range.size(); ++i) {
        if(radii[i]>errors[i]*RELATIVE_SMALLNESS) {
            small=false;
            break;
        }
    }

    if(small || splittings==MAXIMUM_SPLITTINGS) {
        gts.adjoin_outer_approximation(range,depth);
    } else {
        Vector<Interval> subdomain1,subdomain2;
        make_lpair(subdomain1,subdomain2)=Ariadne::split(subdomain);
        _subdivision_adjoin_outer_approximation_to(gts,set,subdomain1,depth,errors,splittings+1u);
        _subdivision_adjoin_outer_approximation_to(gts,set,subdomain2,depth,errors,splittings+1u);
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

    AffineSweeper affine_sweeper;
    TaylorConstrainedImageSet set(*this);
    for(uint i=0; i!=nx; ++i) {
        const_cast<IntervalTaylorModel&>(set._function.models()[i]).sweep(affine_sweeper);
    }
    for(uint i=0; i!=nc; ++i) {
        const_cast<IntervalTaylorModel&>(set._constraints[i].model()).sweep(affine_sweeper);
    }
    for(uint i=0; i!=neq; ++i) {
        const_cast<IntervalTaylorModel&>(set._equations[i].model()).sweep(affine_sweeper);
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

    VectorTaylorFunction new_function=combine(set.taylor_function(),ScalarTaylorFunction::identity(ivl,set.sweeper()));

    TaylorConstrainedImageSet result(new_function.domain(),new_function,new_function.sweeper());
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(embed(*iter,ivl));
    }
    for(const_iterator iter=set._equations.begin(); iter!=set._equations.end(); ++iter) {
        result._equations.append(embed(*iter,ivl));
    }

    return result;
}

TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set, const Box& bx) {
    return product(set,TaylorConstrainedImageSet(bx,set.sweeper()));
}

TaylorConstrainedImageSet product(const TaylorConstrainedImageSet& set1, const TaylorConstrainedImageSet& set2) {
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    VectorTaylorFunction new_function=combine(set1.taylor_function(),set2.taylor_function());

    TaylorConstrainedImageSet result(new_function.domain(),new_function,new_function.sweeper());
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

TaylorConstrainedImageSet apply(const IntervalVectorFunction& function, const TaylorConstrainedImageSet& set) {
    TaylorConstrainedImageSet result(set);
    result.apply_map(function);
    return result;
}

TaylorConstrainedImageSet unchecked_apply(const VectorTaylorFunction& function, const TaylorConstrainedImageSet& set) {
    TaylorConstrainedImageSet result(set);
    const VectorTaylorFunction& space_function=result.function();
    const_cast<VectorTaylorFunction&>(space_function)=Ariadne::unchecked_compose(function,set.function());
    return result;
}

} // namespace Ariadne


