/***************************************************************************
 *            enclosure.cc
 *
 *  Copyright 2008-11  Pieter Collins
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

#include "config.h"

#include "enclosure.h"

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

#include "function_model.h"
#include "taylor_function.h"

#include "box.h"
#include "grid.h"

#include "function_set.h"
#include "affine_set.h"

#include "paving_interface.h"
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
#include <include/space.h>


namespace Ariadne {

static const uint verbosity = 0u;

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

void subdivision_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                            const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void affine_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                       const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void constraint_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                           const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void procedure_constraint_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                                     const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);

void optimal_constraint_adjoin_outer_approximation(PavingInterface& paving, const IntervalVector& domain, const IntervalVectorFunction& function,
                                                   const IntervalVectorFunction& constraint_function, const IntervalVector& constraint_bounds, int depth);


namespace {

Interval make_domain(const IntervalSet& ivl) {
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


IntervalVectorFunctionModel make_identity(const BoxSet& bx, const IntervalFunctionModelFactoryInterface& fac) {
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

    IntervalVectorFunctionModel res=fac.create_identity(dom);
    for(uint i=0; i!=bx.dimension(); ++i) {
        res[i]=res[i]+Interval(-errs[i],+errs[i]);
    }

    return res;
};

// TODO: Make more efficient
inline void assign_all_but_last(MultiIndex& r, const MultiIndex& a) {
    for(uint i=0; i!=r.size(); ++i) { r[i]=a[i]; }
}
} // namespace

Pair<IntervalScalarFunctionModel,IntervalScalarFunctionModel> split(const IntervalScalarFunctionModel& f, Nat k) {
    Pair<IntervalVector,IntervalVector> domains=split(f.domain(),k);
    return make_pair(restrict(f,domains.first),restrict(f,domains.second));
}

Pair<IntervalVectorFunctionModel,IntervalVectorFunctionModel> split(const IntervalVectorFunctionModel& f, Nat k) {
    Pair<IntervalVector,IntervalVector> domains=split(f.domain(),k);
    return make_pair(restrict(f,domains.first),restrict(f,domains.second));
}


void Enclosure::_check() const {
    ARIADNE_ASSERT_MSG(this->_space_function.argument_size()==this->domain().size(),*this);
    ARIADNE_ASSERT_MSG(this->_time_function.argument_size()==this->domain().size(),*this<<"\n\n"<<this->_domain<<"\n"<<this->_time_function<<"\n\n");
    ARIADNE_ASSERT_MSG(this->_dwell_time_function.argument_size()==this->domain().size(),*this<<"\n\n"<<this->_domain<<"\n"<<this->_dwell_time_function<<"\n\n");
    for(List<IntervalConstraintModel>::const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->function().argument_size()==this->domain().size(),*this);
    }
}

IntervalFunctionModelFactoryInterface const&
Enclosure::function_factory() const {
    return *this->_function_factory_ptr;
}

/*
// FIXME: What if solving for constraint leaves domain?
void Enclosure::_solve_zero_constraints() {
    this->_check();
    for(List<IntervalScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ) {
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
            this->_space_function=this->function_factory().create(new_domain,Ariadne::substitute(this->_space_function.models(),k,substitution_model));
            for(List<IntervalScalarFunctionModel>::iterator constraint_iter=this->_negative_constraints.begin();
                    constraint_iter!=this->_negative_constraints.end(); ++constraint_iter) {
                IntervalScalarFunctionModel& constraint=*constraint_iter;
                constraint=this->function_factory().create(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            for(List<IntervalScalarFunctionModel>::iterator constraint_iter=this->_zero_constraints.begin();
                    constraint_iter!=this->_zero_constraints.end(); ++constraint_iter) {
                IntervalScalarFunctionModel& constraint=*constraint_iter;
                constraint=this->function_factory().create(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            // Since we are using an std::vector, assign iterator to next element
            iter=this->_zero_constraints.erase(iter);
            this->_check();
        } else {
            ARIADNE_WARN("No method for solving constraint "<<*iter<<" currently implemented.");
            ++iter;
        }
    }
}
*/

Enclosure::Enclosure()
    : _domain(), _space_function(), _time_function(), _dwell_time_function(), _reduced_domain(), _is_fully_reduced(true)
{
}

Enclosure* Enclosure::clone() const
{
    return new Enclosure(*this);
}

Enclosure::Enclosure(const BoundedConstraintSet& set, const IntervalFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    this->_space_function=make_identity(set.domain(),this->function_factory());
    this->_domain=this->_space_function.domain();
    this->_time_function=this->function_factory().create_zero(this->domain());
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());
    for(uint i=0; i!=set.number_of_constraints(); ++i) {
        this->new_state_constraint(set.constraint(i));
    }
    this->_reduced_domain=this->_domain;
    this->_is_fully_reduced=true;
    this->_check();
}

Enclosure::Enclosure(const Box& box, const IntervalFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
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


    this->_space_function=this->function_factory().create_zeros(box.dimension(),this->_domain);
    this->_time_function=this->function_factory().create_zero(this->_domain);
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());
    uint j=0;
    proper_coordinates.append(box.dimension());
    for(uint i=0; i!=box.dimension(); ++i) {
        if(proper_coordinates[j]==i) {
            this->_space_function[i]=this->function_factory().create_coordinate(this->_domain,j);
            ++j;
        } else {
            this->_space_function[i]=this->function_factory().create_constant(this->_domain,box[i]);
        }
    }
    this->_reduced_domain=this->_domain;
    this->_is_fully_reduced=true;
    this->_check();
}


Enclosure::Enclosure(const IntervalVector& domain, const IntervalVectorFunction& function, const IntervalFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    this->_domain=domain;
    this->_space_function=this->function_factory().create(this->_domain,function);
    this->_time_function=this->function_factory().create_zero(this->_domain);
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());
    this->_reduced_domain=this->_domain;
    this->_check();
}

Enclosure::Enclosure(const IntervalVector& domain, const IntervalVectorFunction& function, const List<IntervalConstraint>& constraints, const IntervalFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    const double min=std::numeric_limits<double>::min();
    this->_domain=domain;
    for(uint i=0; i!=this->_domain.size(); ++i) {
        if(this->_domain[i].radius()==0) {
            this->_domain[i]+=Interval(-min,+min);
        }
    }

    this->_space_function=this->function_factory().create(this->_domain,function);
    this->_time_function=this->function_factory().create_zero(this->_domain);
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());

    for(uint i=0; i!=constraints.size(); ++i) {
        ARIADNE_ASSERT_MSG(domain.size()==constraints[i].function().argument_size(),"domain="<<domain<<", constraint="<<constraints[i]);
        this->new_parameter_constraint(constraints[i]);
    }

    this->_reduced_domain=domain;
    this->_check();
    this->reduce();
    this->_check();
}

Enclosure::Enclosure(const IntervalVector& domain, const IntervalVectorFunction& space_function, const IntervalScalarFunction& time_function, const List<IntervalConstraint>& constraints, const IntervalFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    ARIADNE_ASSERT_MSG(domain.size()==space_function.argument_size(),"domain="<<domain<<", space_function="<<space_function);
    ARIADNE_ASSERT_MSG(domain.size()==time_function.argument_size(),"domain="<<domain<<", time_function="<<time_function);
    const double min=std::numeric_limits<double>::min();
    this->_domain=domain;
    for(uint i=0; i!=this->_domain.size(); ++i) {
        if(this->_domain[i].radius()==0) {
            this->_domain[i]+=Interval(-min,+min);
        }
    }

    this->_space_function=this->function_factory().create(this->_domain,space_function);
    this->_time_function=this->function_factory().create(this->_domain,time_function);
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());

    for(uint i=0; i!=constraints.size(); ++i) {
        ARIADNE_ASSERT_MSG(domain.size()==constraints[i].function().argument_size(),"domain="<<domain<<", constraint="<<constraints[i]);
        this->new_parameter_constraint(constraints[i]);
    }

    this->_reduced_domain=domain;
    this->_check();
    this->reduce();
    this->_check();
}





// Returns true if the entire set is positive; false if entire set is negative
tribool Enclosure::satisfies(IntervalScalarFunction constraint) const
{
    Interval constraint_range=constraint(this->codomain());
    if(constraint_range.upper()<0.0) { return false; }
    if(constraint_range.lower()>0.0) { return true; }
    return indeterminate;
}


/*
void Enclosure::substitute(uint j, IntervalScalarFunctionModel v)
{
    ARIADNE_ASSERT_MSG(v.argument_size()+1u==this->number_of_parameters(),
                       "number_of_parameters="<<this->number_of_parameters()<<", variable="<<v);
                       this->_space_function = Ariadne::substitute(this->_space_function,j,v);
                       for(List<IntervalScalarFunctionModel>::iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }
                       for(List<IntervalScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }

                       this->_check();
}

void Enclosure::substitute(uint j, Float c)
{
    this->_space_function = Ariadne::partial_evaluate(this->_space_function,j,c);
    for(List<IntervalScalarFunctionModel>::iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
    for(List<IntervalScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
    this->_check();
}
*/

void Enclosure::new_parameter(Interval ivl)
{
    this->_domain=join(this->_domain,ivl);
    this->_reduced_domain=join(this->_reduced_domain,ivl);
    this->_space_function=embed(this->_space_function,ivl);
    this->_time_function=embed(this->_time_function,ivl);
    this->_dwell_time_function=embed(this->_dwell_time_function,ivl);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        IntervalConstraintModel& constraint=this->_constraints[i];
        constraint.set_function(embed(constraint.function(),ivl));
    }
    this->_check();
}

void Enclosure::new_variable(Interval ivl)
{
    IntervalScalarFunctionModel variable_function = this->function_factory().create_identity(ivl);
    this->_domain=join(this->_domain,ivl);
    this->_reduced_domain=join(this->_reduced_domain,ivl);
    this->_space_function=combine(this->_space_function,variable_function);
    this->_time_function=embed(this->_time_function,ivl);
    this->_dwell_time_function=embed(this->_dwell_time_function,ivl);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        IntervalConstraintModel& constraint=this->_constraints[i];
        constraint.set_function(embed(constraint.function(),ivl));
    }
    this->_check();
}

void Enclosure::apply_map(IntervalVectorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    this->_space_function=compose(map,this->_space_function);
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());
    this->_check();
}

/*
void Enclosure::apply_flow(IntervalVectorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_space_function=compose(flow,combine(this->_space_function,this->function_factory().create_identity(Vector<Interval>(1u,time))));
    for(List<IntervalScalarFunctionModel>::iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<IntervalScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}
*/

void Enclosure::apply_fixed_evolve_step(IntervalVectorFunction flow, Float time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    IntervalScalarFunctionModel evolve_time_function=this->function_factory().create_constant(this->domain(),ExactFloat(time));
    this->_space_function=compose(flow,join(this->_space_function,evolve_time_function));
    this->_time_function=this->_time_function + evolve_time_function;
    this->_dwell_time_function=this->_dwell_time_function + evolve_time_function;
    this->_check();
}

void Enclosure::apply_space_evolve_step(IntervalVectorFunction flow, IntervalScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", time="<<time);
    IntervalScalarFunctionModel evolve_time_function=compose(time,this->_space_function);
    this->_space_function=compose(flow,join(this->_space_function,evolve_time_function));
    this->_time_function=this->_time_function + evolve_time_function;
    this->_dwell_time_function=this->_dwell_time_function + evolve_time_function;
    this->_check();
}
void Enclosure::apply_spacetime_evolve_step(IntervalVectorFunction flow, IntervalScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", time="<<time);
    IntervalScalarFunctionModel evolve_time_function=compose(time,join(this->_space_function,this->_time_function));
    this->_space_function=compose(flow,join(this->_space_function,evolve_time_function));
    this->_time_function=this->_time_function + evolve_time_function;
    this->_dwell_time_function=this->_dwell_time_function + evolve_time_function;
    this->_check();
}

void Enclosure::apply_parameter_evolve_step(IntervalVectorFunction flow, IntervalScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->number_of_parameters(),"number_of_parameters="<<this->number_of_parameters()<<", time="<<time);
    this->_space_function=compose(flow,join(this->_space_function,this->function_factory().create(this->_space_function.domain(),time)));
    this->_time_function=this->_time_function + time;
    this->_dwell_time_function=this->_dwell_time_function + time;
    this->_check();
}

void Enclosure::apply_finishing_parameter_evolve_step(IntervalVectorFunction flow, IntervalScalarFunction finishing_time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(finishing_time.argument_size()==this->number_of_parameters(),"number_of_parameters="<<this->number_of_parameters()<<", finishing_time="<<finishing_time);
    IntervalScalarFunctionModel omega=this->function_factory().create(this->domain(),finishing_time);
    this->_space_function=compose(flow,join(this->_space_function,omega-this->_time_function));
    this->_dwell_time_function=this->_dwell_time_function + (omega-this->_time_function);
    this->_time_function=omega;
    this->_check();
}

void Enclosure::apply_full_reach_step(IntervalVectorFunctionModel phi)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    Float h=phi.domain()[phi.result_size()].upper();
    IntervalScalarFunctionModel elps=this->function_factory().create_constant(this->domain(),ExactFloat(h));
    this->apply_parameter_reach_step(phi,elps);
}

void Enclosure::apply_spacetime_reach_step(IntervalVectorFunctionModel phi, IntervalScalarFunction elps)
{
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->dimension()+1);
    this->apply_parameter_reach_step(phi,compose(elps,join(this->space_function(),this->time_function())));
}

void Enclosure::apply_parameter_reach_step(IntervalVectorFunctionModel phi, IntervalScalarFunction elps)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    Float h=phi.domain()[phi.result_size()].upper();
    IntervalVector parameter_domain=this->parameter_domain();
    Interval time_domain=Interval(0,h);
    IntervalScalarFunctionModel time_function=this->function_factory().create_identity(time_domain);
    this->new_variable(time_domain);
    ARIADNE_ASSERT(phi.argument_size()==this->dimension());
    this->apply_map(phi);
    IntervalVector new_domain=this->parameter_domain();
    IntervalScalarFunctionModel time_step_function=this->function_factory().create_coordinate(new_domain,new_domain.size()-1u);
    this->_time_function=this->_time_function+time_step_function;
    this->_dwell_time_function=this->_dwell_time_function+time_step_function;
    if(phi.domain()[phi.result_size()].lower()<time_domain.upper()) {
        this->new_negative_parameter_constraint(time_step_function-embed(elps,time_domain));
    }
    this->_check();
}

void Enclosure::new_state_constraint(IntervalConstraint constraint) {
    ARIADNE_ASSERT(constraint.function().argument_size()==this->dimension());
    this->_is_fully_reduced=false;
    IntervalScalarFunctionModel composed_function_model=compose(constraint.function(),this->_space_function);
    this->_constraints.append(IntervalConstraintModel(constraint.lower_bound(),composed_function_model,constraint.upper_bound()));
}

void Enclosure::new_parameter_constraint(IntervalConstraint constraint) {
    ARIADNE_ASSERT(constraint.function().argument_size()==this->number_of_parameters());
    this->_is_fully_reduced=false;
    IntervalScalarFunctionModel function_model=this->function_factory().create(this->domain(),constraint.function());
    this->_constraints.append(IntervalConstraintModel(constraint.lower_bound(),function_model,constraint.upper_bound()));
}


void Enclosure::new_positive_state_constraint(IntervalScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(-compose(constraint_function,this->space_function())>=0.0);
}

void Enclosure::new_negative_state_constraint(IntervalScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(compose(constraint_function,this->space_function())<=0.0);
}

void Enclosure::new_zero_state_constraint(IntervalScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(compose(constraint_function,this->space_function())==0.0);
}

void Enclosure::new_negative_parameter_constraint(IntervalScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(this->function_factory().create(this->domain(),constraint_function)<=0.0);
}

void Enclosure::new_zero_parameter_constraint(IntervalScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(this->function_factory().create(this->domain(),constraint_function)==0.0);
}





IntervalVector Enclosure::domain() const {
    return this->_domain;
}

IntervalVector Enclosure::parameter_domain() const {
    return this->_domain;
}

IntervalVector Enclosure::reduced_domain() const {
    return this->_reduced_domain;
}

IntervalVector Enclosure::codomain() const {
    return Box(this->_space_function.range()).bounding_box();
}

IntervalVectorFunctionModel const& Enclosure::function() const {
    return this->_space_function;
}

IntervalVectorFunctionModel const& Enclosure::space_function() const {
    return this->_space_function;
}

IntervalScalarFunctionModel const& Enclosure::time_function() const {
    return this->_time_function;
}

IntervalScalarFunctionModel const& Enclosure::dwell_time_function() const {
    return this->_dwell_time_function;
}

uint Enclosure::number_of_constraints() const {
    return this->_constraints.size();
}

List<IntervalConstraintModel> const& Enclosure::constraint_models() const {
    return this->_constraints;
}

List<IntervalConstraint> const Enclosure::constraints() const {
    List<IntervalConstraint> result;
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        result.append(IntervalConstraint(this->_constraints[i].lower_bound(),this->_constraints[i].function(),this->_constraints[i].upper_bound()));
    }
    return result;
}

IntervalConstraintModel const& Enclosure::constraint(uint i) const {
    return this->_constraints[i];
}

IntervalVectorFunctionModel const Enclosure::constraint_function() const {
    IntervalVectorFunctionModel g=this->function_factory().create_zeros(this->number_of_constraints(),this->domain());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        g.set(i,this->constraint(i).function());
    }
    return g;
}

IntervalVector const Enclosure::constraint_bounds() const {
    IntervalVector c(this->number_of_constraints());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        c[i]=this->constraint(i).bounds();
    }
    return c;
}

uint Enclosure::dimension() const {
    return this->_space_function.result_size();
}

uint Enclosure::number_of_parameters() const {
    return this->_space_function.argument_size();
}

Box Enclosure::bounding_box() const {
    return Box(this->_space_function.codomain()).bounding_box();
}

Float Enclosure::radius() const {
    return this->bounding_box().radius();
}

Point Enclosure::centre() const {
    return this->bounding_box().centre();
}


tribool
Enclosure::satisfies(IntervalConstraint c) const
{
    Enclosure copy=*this;
    copy.new_state_constraint(c);
    if(copy.empty()) { return false; }
    else { return indeterminate; }
}

tribool Enclosure::bounded() const
{
    return Box(this->domain()).bounded() || indeterminate;
}

tribool Enclosure::empty() const
{
    if(definitely(Ariadne::empty(this->_reduced_domain))) { return true; }
    if(this->_constraints.empty()) { return Ariadne::empty(this->domain()); }
    if(!this->_is_fully_reduced) { this->reduce(); this->reduce(); this->reduce(); }

    for(uint i=0; i!=this->_constraints.size(); ++i) {
        Interval constraint_range = this->_constraints[i].function().evaluate(this->_reduced_domain);
        if( disjoint(constraint_range,this->_constraints[i].bounds()) ) {
            if(this->_reduced_domain.size()>0) { this->_reduced_domain[0] = Interval(1,-1); }
            return true;
        }
    }
    if(Ariadne::empty(this->_reduced_domain)) { return true; }
    return indeterminate;
}

tribool Enclosure::inside(const Box& bx) const
{
    return Ariadne::subset(this->_space_function.evaluate(this->_reduced_domain),bx);
}

tribool Enclosure::subset(const Box& bx) const
{
    this->reduce();

    return Ariadne::subset(this->_space_function.evaluate(this->_reduced_domain),bx) || indeterminate;

}

tribool Enclosure::separated(const Box& bx) const
{
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"Enclosure::subset(Box): self="<<*this<<", box="<<bx);
    List<IntervalConstraint> constraints = this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(this->_reduced_domain,constraints);

    if(_reduced_domain.empty()) { return true; }

    const Box test_domain=this->_reduced_domain;
    for(uint i=0; i!=bx.dimension(); ++i) {
        constraints.append(IntervalScalarFunctionModel(this->_space_function[i]) >= bx[i].lower());
        constraints.append(IntervalScalarFunctionModel(this->_space_function[i]) <= bx[i].upper());
    }
    return !contractor.feasible(test_domain,constraints).first;
}

void Enclosure::reduce() const
{
    List<IntervalConstraint> constraints=this->constraints();
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
    List<IntervalScalarFunctionModel>& mutable_constraints=const_cast<List<IntervalScalarFunctionModel>&>(this->_negative_constraints);
    for(uint i=0; i!=mutable_constraints.size(); ++i) {
        if(mutable_constraints[i](this->_reduced_domain).upper()<0.0) { redundant_constraints.append(i); }
        else { if(i>j) { mutable_constraints[j]=mutable_constraints[j]; } ++j; }
    }
    mutable_constraints.resize(j);
*/


}



Matrix<Float> nonlinearities_zeroth_order(const IntervalVectorFunction& f, const IntervalVector& dom);
Pair<uint,double> nonlinearity_index_and_error(const IntervalVectorFunction& function, const IntervalVector domain);
Pair<uint,double> lipschitz_index_and_error(const IntervalVectorFunction& function, const IntervalVector& domain);

Pair<Enclosure,Enclosure>
Enclosure::split_zeroth_order() const
{
    return this->split(this->splitting_index_zeroth_order());
}


List<IntervalVector>
Enclosure::splitting_subdomains_zeroth_order() const
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
Enclosure::splitting_index_zeroth_order() const
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


Pair<Enclosure,Enclosure>
Enclosure::split_first_order() const
{
    Matrix<Float> nonlinearities=Ariadne::nonlinearities_zeroth_order(this->_space_function,this->_reduced_domain);

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

    return this->split(jmax_in_row_imax);
}


Pair<Enclosure,Enclosure>
Enclosure::split() const
{
    return this->split_zeroth_order();
}

Pair<Enclosure,Enclosure>
Enclosure::split(uint d) const
{
    ARIADNE_PRECONDITION(d<this->number_of_parameters());
    Vector<Interval> subdomain1,subdomain2;
    make_lpair(subdomain1,subdomain2)=Ariadne::split(this->_space_function.domain(),d);

    IntervalVectorFunctionModel function1,function2;
    make_lpair(function1,function2)=Ariadne::split(this->_space_function,d);

    Pair<Enclosure,Enclosure>
    result=make_pair(Enclosure(function1.domain(),function1,this->function_factory()),
                     Enclosure(function2.domain(),function2,this->function_factory()));
    Enclosure& result1=result.first;
    Enclosure& result2=result.second;

    IntervalScalarFunctionModel constraint_function1,constraint_function2;
    for(List<IntervalConstraintModel>::const_iterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        const IntervalConstraintModel& constraint=*iter;
        make_lpair(constraint_function1,constraint_function2)=Ariadne::split(constraint.function(),d);
        result1._constraints.append(IntervalConstraintModel(constraint.lower_bound(),constraint_function1,constraint.upper_bound()));
        result2._constraints.append(IntervalConstraintModel(constraint.lower_bound(),constraint_function2,constraint.upper_bound()));
    }

    IntervalScalarFunctionModel time_function1,time_function2;
    make_lpair(time_function1,time_function2)=Ariadne::split(this->_time_function,d);
    result1._time_function=time_function1;
    result2._time_function=time_function2;
    IntervalScalarFunctionModel dwell_time_function1,dwell_time_function2;
    make_lpair(dwell_time_function1,dwell_time_function2)=Ariadne::split(this->_dwell_time_function,d);
    result1._dwell_time_function=dwell_time_function1;
    result2._dwell_time_function=dwell_time_function2;

    result1._check();
    result2._check();
    return result;
}









typedef Procedure<Interval> IntervalProcedure;


void adjoin_outer_approximation(PavingInterface&, const Box& domain, const IntervalVectorFunction& function, const IntervalVectorFunction& negative_constraints, const IntervalVectorFunction& equality_constraints, int depth);

void Enclosure::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
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

void Enclosure::subdivision_adjoin_outer_approximation_to(PavingInterface& p, int e) const
{
    ARIADNE_LOG(6,"Enclosure::subdivision_adjoin_outer_approximation_to(paving, depth)\n");
    ARIADNE_LOG(7,"  set="<<*this<<", grid="<<p.grid()<<", depth="<<e<<"\n");
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    const IntervalVectorFunctionModel& f=this->function();

    IntervalVectorFunctionModel g=this->constraint_function();
    IntervalVector c=this->constraint_bounds();

    Ariadne::subdivision_adjoin_outer_approximation(p,d,f,g,c,e);
}

void Enclosure::constraint_adjoin_outer_approximation_to(PavingInterface& p, int e) const
{
    ARIADNE_LOG(6,"Enclosure::constraint_adjoin_outer_approximation_to(paving, depth)\n");
    ARIADNE_LOG(7,"  set="<<*this<<", grid="<<p.grid()<<", depth="<<e<<"\n");
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    const IntervalVectorFunctionModel& f=this->function();

    IntervalVectorFunctionModel g=this->constraint_function();
    IntervalVector c=this->constraint_bounds();

    Ariadne::procedure_constraint_adjoin_outer_approximation(p,d,f,g,c,e);
}

void Enclosure::optimal_constraint_adjoin_outer_approximation_to(PavingInterface& p, int e) const
{
    ARIADNE_ASSERT(p.dimension()==this->dimension());
    const Box& d=this->domain();
    IntervalVectorFunctionModel f=this->function();

    IntervalVectorFunctionModel g=this->constraint_function();
    IntervalVector c=this->constraint_bounds();

    Ariadne::optimal_constraint_adjoin_outer_approximation(p,d,f,g,c,e);
}


GridTreeSet Enclosure::outer_approximation(const Grid& grid, int depth) const
{
    GridTreeSet paving(grid);
    this->adjoin_outer_approximation_to(paving,depth);
    return paving;
}



void Enclosure::affine_adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->_reduced_domain,this->_domain),*this);

    // Bound the maximum number of splittings allowed to draw a particular set.
    // Note that this gives rise to possibly 2^MAX_DEPTH split sets!!
    static const int MAXIMUM_DEPTH = 16;

    // The basic approximation error when plotting with accuracy=0
    static const double BASIC_ERROR = 0.0625;

    const double max_error=BASIC_ERROR/(1<<depth);

    IntervalVectorFunctionModel fg=join(this->function(),this->constraint_function());

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





void
Enclosure::recondition()
{
    this->uniform_error_recondition();
    this->kuhn_recondition();
}


void Enclosure::
uniform_error_recondition()
{
    const double MAXIMUM_ERROR = std::numeric_limits<double>::epsilon() * 1024;
    uint old_number_of_parameters = this->number_of_parameters();

    List<uint> large_error_indices;

    for(uint i=0; i!=this->_space_function.result_size(); ++i) {
        Float error=this->_space_function.get(i).error();
        if(error > MAXIMUM_ERROR) {
            large_error_indices.append(i);
        }
    }

    IntervalVector error_domains(large_error_indices.size());
    for(uint i=0; i!=large_error_indices.size(); ++i) {
        Float error=this->_space_function.get(large_error_indices[i]).error();
        error_domains[i]=Interval(-error,+error);
    }
    error_domains=IntervalVector(large_error_indices.size(),Interval(-1,+1));
    uint k=this->number_of_parameters();

    this->_domain=join(this->_domain,error_domains);
    this->_reduced_domain=join(this->_reduced_domain,error_domains);
    this->_space_function=embed(this->_space_function,error_domains);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        this->_constraints[i].function()=embed(this->_constraints[i].function(),error_domains);
    }

    for(uint i=0; i!=large_error_indices.size(); ++i) {
        Float error=this->_space_function.get(large_error_indices[i]).error();
        if(error > MAXIMUM_ERROR) {
            this->_space_function[i].set_error(0.0);
            this->_space_function[i] = this->_space_function.get(i) + this->function_factory().create_coordinate(this->_domain,k)*Interval(error);
            ++k;
        }
    }

    IntervalVector new_variables = project(this->parameter_domain(),range(old_number_of_parameters,this->number_of_parameters()));
    this->_time_function = embed(this->_time_function,new_variables);
    this->_dwell_time_function = embed(this->_dwell_time_function,new_variables);

}



void
Enclosure::kuhn_recondition()
{
    if(!dynamic_cast<const VectorTaylorFunction*>(&this->space_function().reference())) {
        ARIADNE_WARN("Cannot Kuhn reduce an Enclosure which is not given by TaylorFunctions.");
    }

    static const uint NUMBER_OF_BLOCKS = 2;

    const Nat number_of_kept_parameters = (NUMBER_OF_BLOCKS-1)*this->dimension();
    const Nat number_of_discarded_parameters=this->number_of_parameters()-number_of_kept_parameters;
    const Nat number_of_error_parameters = this->dimension();

    if(this->number_of_parameters()<=number_of_kept_parameters) {
        this->uniform_error_recondition();
        return;
    }

    const VectorTaylorFunction& function=dynamic_cast<const VectorTaylorFunction&>(this->space_function().reference());
    const Vector<IntervalTaylorModel>& models = function.models();
    Matrix<Float> dependencies(this->dimension(),this->number_of_parameters());
    for(uint i=0; i!=dependencies.row_size(); ++i) {
        for(IntervalTaylorModel::const_iterator iter=models[i].begin(); iter!=models[i].end(); ++iter) {
            for(uint j=0; j!=dependencies.column_size(); ++j) {
                if(iter->key()[j]!=0) {
                    dependencies[i][j]+=abs(iter->data());
                }
            }
        }
    }
    Array< Pair<Float,Nat> > column_max_dependencies(this->number_of_parameters());
    for(uint j=0; j!=dependencies.column_size(); ++j) {
        column_max_dependencies[j] = make_pair(Float(0.0),Nat(j));
        for(uint i=0; i!=dependencies.row_size(); ++i) {
            column_max_dependencies[j].first=std::max(column_max_dependencies[j].first,dependencies[i][j]);
        }
    }
    std::sort(column_max_dependencies.begin(),column_max_dependencies.end(),std::greater< Pair<Float,Nat> >());

    Array<Nat> kept_parameters(number_of_kept_parameters);
    Array<Nat> discarded_parameters(number_of_discarded_parameters);
    for(uint j=0; j!=number_of_kept_parameters; ++j) { kept_parameters[j]=column_max_dependencies[j].second; }
    for(uint j=0; j!=number_of_discarded_parameters; ++j) { discarded_parameters[j]=column_max_dependencies[number_of_kept_parameters+j].second; }
    std::sort(kept_parameters.begin(),kept_parameters.end());
    std::sort(discarded_parameters.begin(),discarded_parameters.end());

    Vector<IntervalTaylorModel> new_models(models.size(),IntervalTaylorModel(number_of_kept_parameters+number_of_error_parameters,function.sweeper()));
    for(uint i=0; i!=this->dimension(); ++i) {
        new_models[i] = Ariadne::recondition(models[i],discarded_parameters,number_of_error_parameters,i);
    }

    Vector<Interval> new_domain(number_of_kept_parameters+number_of_error_parameters);
    Vector<Interval> new_reduced_domain(number_of_kept_parameters+number_of_error_parameters);
    for(Nat j=0; j!=number_of_kept_parameters; ++j) {
        new_domain[j]=this->parameter_domain()[kept_parameters[j]];
        new_reduced_domain[j]=this->reduced_domain()[kept_parameters[j]];
    }
    for(Nat j=number_of_kept_parameters; j!=number_of_kept_parameters+number_of_error_parameters; ++j) {
        new_domain[j]=Interval(-1,+1);
        new_reduced_domain[j]=Interval(-1,+1);
    }
    this->_domain = new_domain;
    this->_reduced_domain = new_reduced_domain;

    Enclosure new_set(new_domain,VectorTaylorFunction(new_domain,new_models),this->function_factory());
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        const IntervalConstraintModel& constraint=this->_constraints[i];
        ScalarTaylorFunction const& constraint_function=dynamic_cast<const ScalarTaylorFunction&>(constraint.function().reference());
        IntervalScalarFunctionModel new_constraint_function=ScalarTaylorFunction(new_domain,Ariadne::recondition(constraint_function.model(),discarded_parameters,number_of_error_parameters));
        new_set._constraints.append(IntervalConstraintModel(constraint.lower_bound(),new_constraint_function,constraint.upper_bound()));
    }
    ScalarTaylorFunction const& time=dynamic_cast<const ScalarTaylorFunction&>(this->_time_function.reference());
    new_set._time_function=ScalarTaylorFunction(new_domain,Ariadne::recondition(time.model(),discarded_parameters,number_of_error_parameters));
    ScalarTaylorFunction const& dwell_time=dynamic_cast<const ScalarTaylorFunction&>(this->_dwell_time_function.reference());
    new_set._dwell_time_function=ScalarTaylorFunction(new_domain,Ariadne::recondition(dwell_time.model(),discarded_parameters,number_of_error_parameters));

    (*this)=new_set;
    //ScalarTaylorFunction const& dwell_time=dynamic_cast<const ScalarTaylorFunction&>(this->_dwell_time.reference());
    //this->_dwell_time =ScalarTaylorFunction(new_domain,Ariadne::recondition(dwell_time.model(),discarded_parameters,number_of_error_parameters));

    this->_check();
}



void Enclosure::restrict(const Vector<Interval>& subdomain)
{
    ARIADNE_ASSERT_MSG(subdomain.size()==this->number_of_parameters(),"set="<<*this<<", subdomain="<<subdomain);
    ARIADNE_ASSERT_MSG(Ariadne::subset(subdomain,this->domain()),"set.domain()="<<this->domain()<<", subdomain="<<subdomain);
    Enclosure& result(*this);
    result._domain=subdomain;
    result._reduced_domain=Ariadne::intersection(static_cast<const Vector<Interval>&>(result._reduced_domain),subdomain);
    result._space_function=Ariadne::restrict(result._space_function,subdomain);
    result._time_function=Ariadne::restrict(result._time_function,subdomain);
    result._dwell_time_function=Ariadne::restrict(result._dwell_time_function,subdomain);
    IntervalScalarFunctionModel new_constraint;
    for(List<IntervalConstraintModel>::iterator iter=result._constraints.begin();
        iter!=result._constraints.end(); ++iter)
    {
        IntervalScalarFunctionModel& constraint_function=iter->function();
        constraint_function=Ariadne::restrict(constraint_function,subdomain);
    }
    this->reduce();
}

Enclosure Enclosure::restriction(const Vector<Interval>& subdomain) const
{
    Enclosure result(*this);
    result.restrict(subdomain);
    return result;
}


void Enclosure::draw(CanvasInterface& canvas, const Projection2d& projection) const {
    switch(DRAWING_METHOD) {
        case BOX_DRAW:
            this->box_draw(canvas,projection);
            break;
        case AFFINE_DRAW:
            //if(this->number_of_zero_constraints()!=0) { this->box_draw(canvas); }
            this->affine_draw(canvas,projection,DRAWING_ACCURACY);
            break;
        case GRID_DRAW:
            this->grid_draw(canvas,projection);
            break;
        default:
            ARIADNE_WARN("Unknown drawing method\n");
    }
}

void Enclosure::box_draw(CanvasInterface& canvas, const Projection2d& projection) const {
    this->reduce();
    Box(join(this->_space_function(this->_reduced_domain),this->_time_function(this->_reduced_domain))).draw(canvas,projection);
}

IntervalVectorFunctionModel join(const IntervalVectorFunctionModel& f1, const IntervalScalarFunctionModel& f2, const IntervalVectorFunctionModel& f3) {
    return join(join(f1,f2),f3);
}

void Enclosure::affine_draw(CanvasInterface& canvas, const Projection2d& projection, uint accuracy) const {
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->_reduced_domain,this->_domain),*this);

    IntervalVectorFunctionModel cached_space_function=this->_space_function;
    const_cast<Enclosure*>(this)->_space_function=join(this->_space_function,this->_time_function);

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

    IntervalVectorFunctionModel fg=this->function_factory().create_zeros(this->dimension()+1u+this->number_of_constraints(),this->domain());
    for(uint i=0; i!=this->dimension(); ++i) { fg[i]=this->_space_function[i]; }
    for(uint i=0; i!=1; ++i) { fg[i+this->dimension()]=this->_time_function; }
    for(uint i=0; i!=this->_constraints.size(); ++i) { fg[i+this->dimension()+1]=this->_constraints[i].function(); }

//    IntervalVectorFunctionModel fg=join(this->space_function(),this->time_function(),this->constraint_function());

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
            this->restriction(subdomains[n]).affine_over_approximation().draw(canvas,projection);
        } catch(std::runtime_error& e) {
            ARIADNE_WARN("Error "<<e.what()<<" in Enclosure::affine_draw(...) for "<<*this<<"\n");
            this->restriction(subdomains[n]).box_draw(canvas,projection);
        }
    }

    const_cast<Enclosure*>(this)->_space_function=cached_space_function;
};


void Enclosure::grid_draw(CanvasInterface& canvas, const Projection2d& projection, uint accuracy) const {
    // TODO: Project to grid first
    this->outer_approximation(Grid(this->dimension()),accuracy).draw(canvas,projection);
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

const IntervalScalarFunctionModel& repr(const IntervalScalarFunctionModel& f) { return f; }
const IntervalVectorFunctionModel& repr(const IntervalVectorFunctionModel& f) { return f; }
const List<IntervalScalarFunctionModel>& repr(const List<IntervalScalarFunctionModel>& f) { return f; }

std::ostream& Enclosure::write(std::ostream& os) const {
    const bool LONG_FORMAT=false;

    if(LONG_FORMAT) {
        os << "Enclosure"
           << "(\n  domain=" << this->domain()
           << ",\n  range=" << this->bounding_box()
           << ",\n  space_function=" << this->space_function()
           << ",\n  time_function=" << this->time_function()
           << ",\n  constraints=" << this->constraints()
           << "\n)\n";
    } else {
        os << "Enclosure"
           << "( domain=" << this->domain()
           << ", range=" << this->bounding_box()
           << ", space_function=" << repr(this->space_function())
           << ", time_function=" << repr(this->time_function())
           << ", negative_constraints=" << this->constraints()
           << ")";

    } return os;
}





/*
ValidatedAffineConstrainedImageSet
Enclosure::affine_approximation() const
{
    this->_check();
    typedef List<IntervalScalarFunctionModel>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();

    Enclosure set(*this);

    //if(set._zero_constraints.size()>0) { set._solve_zero_constraints(); }
    this->_check();

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np);
    for(uint i=0; i!=nx; ++i) {
        IntervalScalarFunctionModel component=set._space_function[i];
        h[i]=component.model().value();
        for(uint j=0; j!=np; ++j) {
            G[i][j]=component.model().gradient(j);
        }
    }
    ValidatedAffineConstrainedImageSet result(G,h);

    Vector<Float> a(np);
    Float b;

    for(const_iterator iter=set._negative_constraints.begin();
            iter!=set._negative_constraints.end(); ++iter) {
        const IntervalScalarFunctionModel& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=set._zero_constraints.begin();
            iter!=set._zero_constraints.end(); ++iter) {
        const IntervalScalarFunctionModel& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_equality_constraint(a,b);
    }
    return result;
}
*/

/*
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
*/


ValidatedAffineConstrainedImageSet
Enclosure::affine_over_approximation() const
{
    this->_check();
    typedef List<ScalarTaylorFunction>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint nc=this->number_of_constraints();
    const uint np=this->number_of_parameters();

    AffineSweeper affine_sweeper;
    VectorTaylorFunction space_function=dynamic_cast<const VectorTaylorFunction&>(this->_space_function.reference());
    ScalarTaylorFunction time_function=dynamic_cast<const ScalarTaylorFunction&>(this->_time_function.reference());
    List<ScalarTaylorFunction> constraint_functions;
    for(uint i=0; i!=nc; ++i) { constraint_functions.append(dynamic_cast<const ScalarTaylorFunction&>(this->_constraints[i].function().reference())); }

    //std::cerr<<"\n"<<space_function<<"\n"<<time_function<<"\n"<<constraint_functions<<"\n\n";

    Vector< IntervalAffineModel > affine_function_models(space_function.result_size());
    for(uint i=0; i!=space_function.result_size(); ++i) { affine_function_models[i]=affine_model(space_function.models()[i]); }
    //affine_function_models[space_function.result_size()]=affine_model(time_function.model());

    ValidatedAffineConstrainedImageSet result(affine_function_models);
    //std::cerr<<"\n"<<*this<<"\n"<<result<<"\n\n";

    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        IntervalTaylorModel const& constraint_model=constraint_functions[i].model();
        IntervalAffineModel affine_constraint_model=affine_model(constraint_model);
        Interval constraint_bound=this->constraint(i).bounds();
        result.new_constraint(constraint_bound.lower()<=affine_constraint_model<=constraint_bound.upper());
    }

    ARIADNE_LOG(2,"set="<<*this<<"\nset.affine_over_approximation()="<<result<<"\n");
    return result;

}


Enclosure product(const Enclosure& set, const Interval& ivl) {
    typedef List<IntervalConstraintModel>::const_iterator const_iterator;

    IntervalVectorFunctionModel new_function=combine(set.function(),set.function_factory().create_identity(ivl));

    Enclosure result(new_function.domain(),new_function,set.function_factory());
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(IntervalConstraintModel(iter->lower_bound(),embed(iter->function(),ivl),iter->upper_bound()));
    }
    result._time_function=embed(set._time_function,ivl);
    result._dwell_time_function=embed(set._dwell_time_function,ivl);

    return result;
}

Enclosure product(const Enclosure& set, const Box& bx) {
    typedef List<IntervalConstraintModel>::const_iterator const_iterator;

    IntervalVectorFunctionModel new_function=combine(set.function(),set.function_factory().create_identity(bx));

    Enclosure result(new_function.domain(),new_function,set.function_factory());
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(IntervalConstraintModel(iter->lower_bound(),embed(iter->function(),bx),iter->upper_bound()));
    }
    result._time_function=embed(set._time_function,bx);
    result._dwell_time_function=embed(set._dwell_time_function,bx);

    return result;
}

Enclosure product(const Enclosure& set1, const Enclosure& set2) {
    ARIADNE_ASSERT(set1.time_function().range() == set2.time_function().range());
    ARIADNE_ASSERT(set1.dwell_time_function().range() == set2.dwell_time_function().range());

    typedef List<IntervalConstraintModel>::const_iterator const_iterator;

    IntervalVectorFunctionModel new_function=combine(set1.function(),set2.function());

    Enclosure result(new_function.domain(),new_function,set1.function_factory());
    for(const_iterator iter=set1._constraints.begin(); iter!=set1._constraints.end(); ++iter) {
        result._constraints.append(IntervalConstraintModel(iter->lower_bound(),embed(iter->function(),set2.domain()),iter->upper_bound()));
    }
    for(const_iterator iter=set2._constraints.begin(); iter!=set2._constraints.end(); ++iter) {
        result._constraints.append(IntervalConstraintModel(iter->lower_bound(),embed(set1.domain(),iter->function()),iter->upper_bound()));
    }
    result._time_function=embed(set1.time_function(),set2.time_function().domain());
    result._dwell_time_function=embed(set1.dwell_time_function(),set2.dwell_time_function().domain());

    return result;
}

Enclosure apply(const IntervalVectorFunction& function, const Enclosure& set) {
    Enclosure result(set);
    result.apply_map(function);
    return result;
}

Enclosure unchecked_apply(const IntervalVectorFunctionModel& function, const Enclosure& set) {
    Enclosure result(set);
    const IntervalVectorFunctionModel& space_function=result.function();
    const_cast<IntervalVectorFunctionModel&>(space_function)=Ariadne::unchecked_compose(function,set.function());
    return result;
}

} // namespace Ariadne


