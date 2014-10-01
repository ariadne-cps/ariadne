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

#include "functional.h"
#include "config.h"

#include <iomanip>

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
#include "paver.h"
#include "grid_set.h"

#include "constraint_solver.h"
#include "nonlinear_programming.h"

#include "graphics_interface.h"

#include "discrete_event.h"

#include "logging.h"

#include "functional.h"

#include "config.h"
#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif // HAVE_CAIRO_H
#include <boost/concept_check.hpp>
#include "operators.h"
#include "space.h"


namespace Ariadne {

static const uint verbosity = 0u;

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

namespace {

Interval make_domain(const IntervalSet& ivl) {
    rounding_mode_t rnd=get_rounding_mode();
    Interval dom_lower_ivl=Interval(ivl.lower());
    Interval dom_upper_ivl=Interval(ivl.upper());
    Float dom_lower=dom_lower_ivl.lower().raw();
    Float dom_upper=dom_upper_ivl.upper().raw();
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


ValidatedVectorFunctionModel make_identity(const BoxSet& bx, const ValidatedFunctionModelFactoryInterface& fac) {
    IntervalVector dom(bx.dimension());
    FloatVector errs(bx.dimension());

    for(uint i=0; i!=bx.dimension(); ++i) {
        Interval dom_lower_ivl=numeric_cast<Interval>(bx[i].lower());
        Interval dom_upper_ivl=numeric_cast<Interval>(bx[i].upper());
        // Convert to single-precision values
        Float dom_lower_flt=numeric_cast<float>(bx[i].lower());
        Float dom_upper_flt=numeric_cast<float>(bx[i].upper());
        set_rounding_upward();
        Float err=max( max(dom_upper_ivl.upper().raw()-dom_upper_flt,dom_upper_flt-dom_upper_ivl.lower().raw()),
                       max(dom_lower_ivl.upper().raw()-dom_lower_flt,dom_lower_flt-dom_lower_ivl.lower().raw()) );
        set_rounding_to_nearest();
        dom[i]=Interval(dom_lower_flt,dom_upper_flt);
        errs[i]=err;
    }

    ValidatedVectorFunctionModel res=fac.create_identity(dom);
    for(uint i=0; i!=bx.dimension(); ++i) {
        res[i]=res[i]+ValidatedNumberType(-errs[i],+errs[i]);
    }

    return res;
};

// TODO: Make more efficient
inline void assign_all_but_last(MultiIndex& r, const MultiIndex& a) {
    for(uint i=0; i!=r.size(); ++i) { r[i]=a[i]; }
}
} // namespace

Pair<ValidatedScalarFunctionModel,ValidatedScalarFunctionModel> split(const ValidatedScalarFunctionModel& f, Nat k) {
    Pair<Box,Box> domains=split(f.domain(),k);
    return make_pair(restrict(f,domains.first),restrict(f,domains.second));
}

Pair<ValidatedVectorFunctionModel,ValidatedVectorFunctionModel> split(const ValidatedVectorFunctionModel& f, Nat k) {
    Pair<Box,Box> domains=split(f.domain(),k);
    return make_pair(restrict(f,domains.first),restrict(f,domains.second));
}


void Enclosure::_check() const {
    ARIADNE_ASSERT_MSG(this->_space_function.argument_size()==this->domain().size(),*this);
    ARIADNE_ASSERT_MSG(this->_time_function.argument_size()==this->domain().size(),*this<<"\n\n"<<this->_domain<<"\n"<<this->_time_function<<"\n\n");
    ARIADNE_ASSERT_MSG(this->_dwell_time_function.argument_size()==this->domain().size(),*this<<"\n\n"<<this->_domain<<"\n"<<this->_dwell_time_function<<"\n\n");
    for(List<ValidatedConstraintModel>::const_iterator iter=this->_constraints.begin(); iter!=this->_constraints.end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->function().argument_size()==this->domain().size(),*this);
    }
}

ValidatedFunctionModelFactoryInterface const&
Enclosure::function_factory() const {
    return *this->_function_factory_ptr;
}

/*
// FIXME: What if solving for constraint leaves domain?
void Enclosure::_solve_zero_constraints() {
    this->_check();
    for(List<ValidatedScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ) {
        const Box& domain=this->domain();
        const ValidatedTaylorModel& model=iter->model();
        const uint k=model.argument_size()-1u;
        ValidatedTaylorModel zeroth_order(k,this->sweeper());
        ValidatedTaylorModel first_order(k,this->sweeper());
        bool is_zeroth_order=true;
        bool is_first_order=true;
        MultiIndex r(k);
        // Try linear approach in last coefficient
        for(ValidatedTaylorModel::const_iterator tmiter=model.begin(); tmiter!=model.end(); ++tmiter) {
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
            const Box new_domain=project(domain,range(0,k));
            ValidatedTaylorModel substitution_model=-zeroth_order/first_order;
            this->_space_function=this->function_factory().create(new_domain,Ariadne::substitute(this->_space_function.models(),k,substitution_model));
            for(List<ValidatedScalarFunctionModel>::iterator constraint_iter=this->_negative_constraints.begin();
                    constraint_iter!=this->_negative_constraints.end(); ++constraint_iter) {
                ValidatedScalarFunctionModel& constraint=*constraint_iter;
                constraint=this->function_factory().create(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            for(List<ValidatedScalarFunctionModel>::iterator constraint_iter=this->_zero_constraints.begin();
                    constraint_iter!=this->_zero_constraints.end(); ++constraint_iter) {
                ValidatedScalarFunctionModel& constraint=*constraint_iter;
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

Enclosure::Enclosure(const BoundedConstraintSet& set, const ValidatedFunctionModelFactoryInterface& factory)
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

Enclosure::Enclosure(const Box& box, const ValidatedFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    // Ensure domain elements have nonempty radius
    const float min_float=std::numeric_limits<float>::min();
    List<uint> proper_coordinates;
    proper_coordinates.reserve(box.dimension());
    for(uint i=0; i!=box.dimension(); ++i) {
        if(box[i].width()>=min_float) {
            proper_coordinates.append(i);
        }
    }
    this->_domain=Box(proper_coordinates.size());
    for(uint j=0; j!=this->_domain.size(); ++j) {
        this->_domain[j]=box[proper_coordinates[j]];
    }

    // HACK: Make a dummy variable for the domain to avoid bugs which occur
    // with a zero-dimensional domain.
    // FIXME: Fix issues with TaylorFunction on zero-dimensional domain.
    if(proper_coordinates.size()==0) { this->_domain=Box(1u,Interval(-1,+1)); }


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
            this->_space_function[i]=this->function_factory().create_constant(this->_domain,make_singleton(box[i]));
        }
    }
    this->_reduced_domain=this->_domain;
    this->_is_fully_reduced=true;
    this->_check();
}


Enclosure::Enclosure(const Box& domain, const ValidatedVectorFunction& function, const ValidatedFunctionModelFactoryInterface& factory)
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

Enclosure::Enclosure(const Box& domain, const ValidatedVectorFunction& function, const List<ValidatedConstraint>& constraints, const ValidatedFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    ARIADNE_ASSERT_MSG(domain.size()==function.argument_size(),"domain="<<domain<<", function="<<function);
    const double min=std::numeric_limits<double>::min();
    this->_domain=domain;
    for(uint i=0; i!=this->_domain.size(); ++i) {
        if(this->_domain[i].width()==0) {
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

Enclosure::Enclosure(const Box& domain, const ValidatedVectorFunction& space_function, const ValidatedScalarFunction& time_function, const List<ValidatedConstraint>& constraints, const ValidatedFunctionModelFactoryInterface& factory)
    : _function_factory_ptr(factory.clone())
{
    ARIADNE_ASSERT_MSG(domain.size()==space_function.argument_size(),"domain="<<domain<<", space_function="<<space_function);
    ARIADNE_ASSERT_MSG(domain.size()==time_function.argument_size(),"domain="<<domain<<", time_function="<<time_function);
    const double min=std::numeric_limits<double>::min();
    this->_domain=domain;
    for(uint i=0; i!=this->_domain.size(); ++i) {
        if(this->_domain[i].width()==0) {
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
tribool Enclosure::satisfies(ValidatedScalarFunction constraint) const
{
    Interval constraint_range=apply(constraint,this->codomain());
    if(constraint_range.upper()<0.0) { return false; }
    if(constraint_range.lower()>0.0) { return true; }
    return indeterminate;
}


/*
void Enclosure::substitute(uint j, ValidatedScalarFunctionModel v)
{
    ARIADNE_ASSERT_MSG(v.argument_size()+1u==this->number_of_parameters(),
                       "number_of_parameters="<<this->number_of_parameters()<<", variable="<<v);
                       this->_space_function = Ariadne::substitute(this->_space_function,j,v);
                       for(List<ValidatedScalarFunctionModel>::iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }
                       for(List<ValidatedScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }

                       this->_check();
}

void Enclosure::substitute(uint j, Float c)
{
    this->_space_function = Ariadne::partial_evaluate(this->_space_function,j,c);
    for(List<ValidatedScalarFunctionModel>::iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
    for(List<ValidatedScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
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
        ValidatedConstraintModel& constraint=this->_constraints[i];
        constraint.set_function(embed(constraint.function(),ivl));
    }
    this->_check();
}

void Enclosure::new_variable(Interval ivl)
{
    ValidatedScalarFunctionModel variable_function = this->function_factory().create_identity(ivl);
    this->_domain=join(this->_domain,ivl);
    this->_reduced_domain=join(this->_reduced_domain,ivl);
    this->_space_function=combine(this->_space_function,variable_function);
    this->_time_function=embed(this->_time_function,ivl);
    this->_dwell_time_function=embed(this->_dwell_time_function,ivl);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        ValidatedConstraintModel& constraint=this->_constraints[i];
        constraint.set_function(embed(constraint.function(),ivl));
    }
    this->_check();
}

void Enclosure::apply_map(ValidatedVectorFunction map)
{
    ARIADNE_ASSERT_MSG(map.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", map="<<map);
    this->_space_function=compose(map,this->_space_function);
    this->_dwell_time_function=this->function_factory().create_zero(this->domain());
    this->_check();
}

/*
void Enclosure::apply_flow(ValidatedVectorFunction flow, Interval time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    this->_space_function=compose(flow,combine(this->_space_function,this->function_factory().create_identity(Box(1u,time))));
    for(List<ValidatedScalarFunctionModel>::iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    for(List<ValidatedScalarFunctionModel>::iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
        *iter=embed(*iter,time);
    }
    this->_check();
}
*/

void Enclosure::apply_fixed_evolve_step(ValidatedVectorFunction flow, ExactFloat time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ValidatedScalarFunctionModel evolve_time_function=this->function_factory().create_constant(this->domain(),ExactFloat(time));
    this->_space_function=compose(flow,join(this->_space_function,evolve_time_function));
    this->_time_function=this->_time_function + evolve_time_function;
    this->_dwell_time_function=this->_dwell_time_function + evolve_time_function;
    this->_check();
}

void Enclosure::apply_space_evolve_step(ValidatedVectorFunction flow, ValidatedScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", time="<<time);
    ValidatedScalarFunctionModel evolve_time_function=compose(time,this->_space_function);
    this->_space_function=compose(flow,join(this->_space_function,evolve_time_function));
    this->_time_function=this->_time_function + evolve_time_function;
    this->_dwell_time_function=this->_dwell_time_function + evolve_time_function;
    this->_check();
}
void Enclosure::apply_spacetime_evolve_step(ValidatedVectorFunction flow, ValidatedScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", time="<<time);
    ValidatedScalarFunctionModel evolve_time_function=compose(time,join(this->_space_function,this->_time_function));
    this->_space_function=compose(flow,join(this->_space_function,evolve_time_function));
    this->_time_function=this->_time_function + evolve_time_function;
    this->_dwell_time_function=this->_dwell_time_function + evolve_time_function;
    this->_check();
}

void Enclosure::apply_parameter_evolve_step(ValidatedVectorFunction flow, ValidatedScalarFunction time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(time.argument_size()==this->number_of_parameters(),"number_of_parameters="<<this->number_of_parameters()<<", time="<<time);
    this->_space_function=compose(flow,join(this->_space_function,this->function_factory().create(this->_space_function.domain(),time)));
    this->_time_function=this->_time_function + time;
    this->_dwell_time_function=this->_dwell_time_function + time;
    this->_check();
}

void Enclosure::apply_finishing_parameter_evolve_step(ValidatedVectorFunction flow, ValidatedScalarFunction finishing_time)
{
    ARIADNE_ASSERT_MSG(flow.argument_size()==this->dimension()+1u,"dimension="<<this->dimension()<<", flow="<<flow);
    ARIADNE_ASSERT_MSG(finishing_time.argument_size()==this->number_of_parameters(),"number_of_parameters="<<this->number_of_parameters()<<", finishing_time="<<finishing_time);
    ValidatedScalarFunctionModel omega=this->function_factory().create(this->domain(),finishing_time);
    this->_space_function=compose(flow,join(this->_space_function,omega-this->_time_function));
    this->_dwell_time_function=this->_dwell_time_function + (omega-this->_time_function);
    this->_time_function=omega;
    this->_check();
}

void Enclosure::apply_full_reach_step(ValidatedVectorFunctionModel phi)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    Float h=phi.domain()[phi.result_size()].upper().raw();
    ValidatedScalarFunctionModel elps=this->function_factory().create_constant(this->domain(),ExactFloat(h));
    this->apply_parameter_reach_step(phi,elps);
}

void Enclosure::apply_spacetime_reach_step(ValidatedVectorFunctionModel phi, ValidatedScalarFunction elps)
{
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->dimension()+1);
    this->apply_parameter_reach_step(phi,compose(elps,join(this->space_function(),this->time_function())));
}

void Enclosure::apply_parameter_reach_step(ValidatedVectorFunctionModel phi, ValidatedScalarFunction elps)
{
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_ASSERT(phi.result_size()==this->dimension());
    ARIADNE_ASSERT(phi.argument_size()==this->dimension()+1);
    ARIADNE_ASSERT(elps.argument_size()==this->number_of_parameters());
    Float h=phi.domain()[phi.result_size()].upper().raw();
    Box parameter_domain=this->parameter_domain();
    Interval time_domain=Interval(0,h);
    ValidatedScalarFunctionModel time_function=this->function_factory().create_identity(time_domain);
    this->new_variable(time_domain);
    ARIADNE_ASSERT(phi.argument_size()==this->dimension());
    this->apply_map(phi);
    Box new_domain=this->parameter_domain();
    ValidatedScalarFunctionModel time_step_function=this->function_factory().create_coordinate(new_domain,new_domain.size()-1u);
    this->_time_function=this->_time_function+time_step_function;
    this->_dwell_time_function=this->_dwell_time_function+time_step_function;
    if(phi.domain()[phi.result_size()].lower()<time_domain.upper()) {
        this->new_negative_parameter_constraint(time_step_function-embed(elps,time_domain));
    }
    this->_check();
}

void Enclosure::new_state_constraint(ValidatedConstraint constraint) {
    ARIADNE_ASSERT(constraint.function().argument_size()==this->dimension());
    this->_is_fully_reduced=false;
    ValidatedScalarFunctionModel composed_function_model=compose(constraint.function(),this->_space_function);
    this->_constraints.append(ValidatedConstraintModel(constraint.lower_bound(),composed_function_model,constraint.upper_bound()));
}

void Enclosure::new_parameter_constraint(ValidatedConstraint constraint) {
    ARIADNE_ASSERT(constraint.function().argument_size()==this->number_of_parameters());
    this->_is_fully_reduced=false;
    ValidatedScalarFunctionModel function_model=this->function_factory().create(this->domain(),constraint.function());
    this->_constraints.append(ValidatedConstraintModel(constraint.lower_bound(),function_model,constraint.upper_bound()));
}


void Enclosure::new_positive_state_constraint(ValidatedScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(-compose(constraint_function,this->space_function())>=0);
}

void Enclosure::new_negative_state_constraint(ValidatedScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(compose(constraint_function,this->space_function())<=0);
}

void Enclosure::new_zero_state_constraint(ValidatedScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->dimension(),"dimension="<<this->dimension()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(compose(constraint_function,this->space_function())==0);
}

void Enclosure::new_negative_parameter_constraint(ValidatedScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(this->function_factory().create(this->domain(),constraint_function)<=0);
}

void Enclosure::new_zero_parameter_constraint(ValidatedScalarFunction constraint_function) {
    ARIADNE_ASSERT_MSG(constraint_function.argument_size()==this->domain().size(),"domain="<<this->domain()<<", constraint_function="<<constraint_function);
    this->_is_fully_reduced=false;
    this->_constraints.append(this->function_factory().create(this->domain(),constraint_function)==0);
}





Box Enclosure::domain() const {
    return this->_domain;
}

Box Enclosure::parameter_domain() const {
    return this->_domain;
}

Box Enclosure::reduced_domain() const {
    return this->_reduced_domain;
}

Box Enclosure::codomain() const {
    return Box(this->_space_function.range()).bounding_box();
}

ValidatedVectorFunctionModel const& Enclosure::function() const {
    return this->_space_function;
}

ValidatedVectorFunctionModel const& Enclosure::space_function() const {
    return this->_space_function;
}

ValidatedScalarFunctionModel const& Enclosure::time_function() const {
    return this->_time_function;
}

ValidatedScalarFunctionModel const& Enclosure::dwell_time_function() const {
    return this->_dwell_time_function;
}

uint Enclosure::number_of_constraints() const {
    return this->_constraints.size();
}

List<ValidatedConstraintModel> const& Enclosure::constraint_models() const {
    return this->_constraints;
}

List<ValidatedConstraint> const Enclosure::constraints() const {
    List<ValidatedConstraint> result;
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        result.append(ValidatedConstraint(this->_constraints[i].lower_bound(),this->_constraints[i].function(),this->_constraints[i].upper_bound()));
    }
    return result;
}

ValidatedConstraintModel const& Enclosure::constraint(uint i) const {
    return this->_constraints[i];
}

ValidatedVectorFunctionModel const Enclosure::constraint_function() const {
    ValidatedVectorFunctionModel g=this->function_factory().create_zeros(this->number_of_constraints(),this->domain());
    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        g.set(i,this->constraint(i).function());
    }
    return g;
}

Box const Enclosure::constraint_bounds() const {
    Box c(this->number_of_constraints());
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

ErrorFloatType Enclosure::radius() const {
    return this->bounding_box().radius();
}

Point Enclosure::centre() const {
    return this->bounding_box().centre();
}


tribool
Enclosure::satisfies(ValidatedConstraint c) const
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
        Interval constraint_range = Ariadne::apply(this->_constraints[i].function(),this->_reduced_domain);
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
    return Ariadne::subset(Ariadne::apply(this->_space_function,this->_reduced_domain),bx);
}

tribool Enclosure::subset(const Box& bx) const
{
    this->reduce();

    return Ariadne::subset(Ariadne::apply(this->_space_function,this->_reduced_domain),bx) || indeterminate;

}

tribool Enclosure::separated(const Box& bx) const
{
    ARIADNE_ASSERT_MSG(this->dimension()==bx.dimension(),"Enclosure::subset(Box): self="<<*this<<", box="<<bx);
    List<ValidatedConstraint> constraints = this->constraints();
    ConstraintSolver contractor=ConstraintSolver();
    contractor.reduce(this->_reduced_domain,constraints);

    if(_reduced_domain.empty()) { return true; }

    const Box test_domain=this->_reduced_domain;
    for(uint i=0; i!=bx.dimension(); ++i) {
        constraints.append(ValidatedScalarFunctionModel(this->_space_function[i]) >= bx[i].lower());
        constraints.append(ValidatedScalarFunctionModel(this->_space_function[i]) <= bx[i].upper());
    }
    return !contractor.feasible(test_domain,constraints).first;
}

void Enclosure::reduce() const
{
    List<ValidatedConstraint> constraints=this->constraints();
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
    List<ValidatedScalarFunctionModel>& mutable_constraints=const_cast<List<ValidatedScalarFunctionModel>&>(this->_negative_constraints);
    for(uint i=0; i!=mutable_constraints.size(); ++i) {
        if(mutable_constraints[i](this->_reduced_domain).upper()<0.0) { redundant_constraints.append(i); }
        else { if(i>j) { mutable_constraints[j]=mutable_constraints[j]; } ++j; }
    }
    mutable_constraints.resize(j);
*/


}



Matrix<Float> nonlinearities_zeroth_order(const ValidatedVectorFunction& f, const Box& dom);
Pair<uint,double> nonlinearity_index_and_error(const ValidatedVectorFunction& function, const Box& domain);
Pair<uint,double> lipschitz_index_and_error(const ValidatedVectorFunction& function, const Box& domain);

Pair<Enclosure,Enclosure>
Enclosure::split_zeroth_order() const
{
    return this->split(this->splitting_index_zeroth_order());
}


List<Box>
Enclosure::splitting_subdomains_zeroth_order() const
{
    List<Box> result;
    uint k=this->splitting_index_zeroth_order();
    if(k==this->number_of_parameters()) {
        result.append(this->_reduced_domain);
    } else {
        Pair<Box,Box> subdomains = this->_reduced_domain.split(this->splitting_index_zeroth_order());
        result.append(subdomains.first);
        result.append(subdomains.second);
    }
    return result;
}


uint
Enclosure::splitting_index_zeroth_order() const
{
    Matrix<Interval> jacobian=Ariadne::jacobian(this->function(),this->reduced_domain());

    // Compute the column of the matrix which has the norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    uint jmax=this->number_of_parameters();
    Float max_column_norm=0.0;
    for(uint j=0; j!=this->number_of_parameters(); ++j) {
        Float column_norm=0.0;
        for(uint i=0; i!=this->dimension(); ++i) {
            column_norm+=mag(jacobian[i][j]);
        }
        column_norm *= this->reduced_domain()[j].error().raw();
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
    Box subdomain1,subdomain2;
    make_lpair(subdomain1,subdomain2)=Ariadne::split(this->_space_function.domain(),d);

    ValidatedVectorFunctionModel function1,function2;
    make_lpair(function1,function2)=Ariadne::split(this->_space_function,d);

    Pair<Enclosure,Enclosure>
    result=make_pair(Enclosure(function1.domain(),function1,this->function_factory()),
                     Enclosure(function2.domain(),function2,this->function_factory()));
    Enclosure& result1=result.first;
    Enclosure& result2=result.second;

    ValidatedScalarFunctionModel constraint_function1,constraint_function2;
    for(List<ValidatedConstraintModel>::const_iterator iter=this->_constraints.begin();
        iter!=this->_constraints.end(); ++iter)
    {
        const ValidatedConstraintModel& constraint=*iter;
        make_lpair(constraint_function1,constraint_function2)=Ariadne::split(constraint.function(),d);
        result1._constraints.append(ValidatedConstraintModel(constraint.lower_bound(),constraint_function1,constraint.upper_bound()));
        result2._constraints.append(ValidatedConstraintModel(constraint.lower_bound(),constraint_function2,constraint.upper_bound()));
    }

    ValidatedScalarFunctionModel time_function1,time_function2;
    make_lpair(time_function1,time_function2)=Ariadne::split(this->_time_function,d);
    result1._time_function=time_function1;
    result2._time_function=time_function2;
    ValidatedScalarFunctionModel dwell_time_function1,dwell_time_function2;
    make_lpair(dwell_time_function1,dwell_time_function2)=Ariadne::split(this->_dwell_time_function,d);
    result1._dwell_time_function=dwell_time_function1;
    result2._dwell_time_function=dwell_time_function2;

    result1._check();
    result2._check();
    return result;
}










void adjoin_outer_approximation(PavingInterface&, const Box& domain, const ValidatedVectorFunction& function, const ValidatedVectorFunction& negative_constraints, const ValidatedVectorFunction& equality_constraints, int depth);

void Enclosure::adjoin_outer_approximation_to(PavingInterface& paving, int depth) const
{
    PavingInterface& p=paving;
    const Box& d=this->domain();
    const ValidatedVectorFunction& f=this->function();

    ValidatedVectorFunctionModel g=this->constraint_function();
    IntervalVector c=this->constraint_bounds();
    Int e=depth;

    switch(DISCRETISATION_METHOD) {
        case SUBDIVISION_DISCRETISE:
            SubdivisionPaver().adjoin_outer_approximation(p,d,f,g,c,e);
            break;
        case AFFINE_DISCRETISE:
            AffinePaver().adjoin_outer_approximation(p,d,f,g,c,e);
            break;
        case CONSTRAINT_DISCRETISE:
            ConstraintPaver().adjoin_outer_approximation(p,d,f,g,c,e);
            break;
        default:
            ARIADNE_FAIL_MSG("Unknown discretisation method\n");
    }
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

    ValidatedVectorFunctionModel fg=join(this->function(),this->constraint_function());

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
        Float error=this->_space_function.get(i).error().raw();
        if(error > MAXIMUM_ERROR) {
            large_error_indices.append(i);
        }
    }

    Box error_domains(large_error_indices.size());
    for(uint i=0; i!=large_error_indices.size(); ++i) {
        Float error=this->_space_function.get(large_error_indices[i]).error().raw();
        error_domains[i]=Interval(-error,+error);
    }
    error_domains=Box(large_error_indices.size(),Interval(-1,+1));
    uint k=this->number_of_parameters();

    this->_domain=join(this->_domain,error_domains);
    this->_reduced_domain=join(this->_reduced_domain,error_domains);
    this->_space_function=embed(this->_space_function,error_domains);
    for(uint i=0; i!=this->_constraints.size(); ++i) {
        this->_constraints[i].function()=embed(this->_constraints[i].function(),error_domains);
    }

    for(uint i=0; i!=large_error_indices.size(); ++i) {
        Float error=this->_space_function.get(large_error_indices[i]).error().raw();
        if(error > MAXIMUM_ERROR) {
            this->_space_function[i].set_error(0.0);
            this->_space_function[i] = this->_space_function.get(i) + this->function_factory().create_coordinate(this->_domain,k)*ValidatedFloatType(+error);
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
    const Vector<ValidatedTaylorModel>& models = function.models();
    Matrix<Float> dependencies(this->dimension(),this->number_of_parameters());
    for(uint i=0; i!=dependencies.row_size(); ++i) {
        for(ValidatedTaylorModel::const_iterator iter=models[i].begin(); iter!=models[i].end(); ++iter) {
            for(uint j=0; j!=dependencies.column_size(); ++j) {
                if(iter->key()[j]!=0) {
                    dependencies[i][j]+=abs(iter->data()).raw();
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

    Vector<ValidatedTaylorModel> new_models(models.size(),ValidatedTaylorModel(number_of_kept_parameters+number_of_error_parameters,function.sweeper()));
    for(uint i=0; i!=this->dimension(); ++i) {
        new_models[i] = Ariadne::recondition(models[i],discarded_parameters,number_of_error_parameters,i);
    }

    Box new_domain(number_of_kept_parameters+number_of_error_parameters);
    Box new_reduced_domain(number_of_kept_parameters+number_of_error_parameters);
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
        const ValidatedConstraintModel& constraint=this->_constraints[i];
        ScalarTaylorFunction const& constraint_function=dynamic_cast<const ScalarTaylorFunction&>(constraint.function().reference());
        ValidatedScalarFunctionModel new_constraint_function=ScalarTaylorFunction(new_domain,Ariadne::recondition(constraint_function.model(),discarded_parameters,number_of_error_parameters));
        new_set._constraints.append(ValidatedConstraintModel(constraint.lower_bound(),new_constraint_function,constraint.upper_bound()));
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



void Enclosure::restrict(const Box& subdomain)
{
    ARIADNE_ASSERT_MSG(subdomain.size()==this->number_of_parameters(),"set="<<*this<<", subdomain="<<subdomain);
    ARIADNE_ASSERT_MSG(Ariadne::subset(subdomain,this->domain()),"set.domain()="<<this->domain()<<", subdomain="<<subdomain);
    Enclosure& result(*this);
    result._domain=subdomain;
    result._reduced_domain=Ariadne::intersection(static_cast<const Box&>(result._reduced_domain),subdomain);
    result._space_function=Ariadne::restrict(result._space_function,subdomain);
    result._time_function=Ariadne::restrict(result._time_function,subdomain);
    result._dwell_time_function=Ariadne::restrict(result._dwell_time_function,subdomain);
    ValidatedScalarFunctionModel new_constraint;
    for(List<ValidatedConstraintModel>::iterator iter=result._constraints.begin();
        iter!=result._constraints.end(); ++iter)
    {
        ValidatedScalarFunctionModel& constraint_function=iter->function();
        constraint_function=Ariadne::restrict(constraint_function,subdomain);
    }
    this->reduce();
}

Enclosure Enclosure::restriction(const Box& subdomain) const
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
    Box(join(Ariadne::apply(this->_space_function,this->_reduced_domain),Ariadne::apply(this->_time_function,this->_reduced_domain))).draw(canvas,projection);
}

ValidatedVectorFunctionModel join(const ValidatedVectorFunctionModel& f1, const ValidatedScalarFunctionModel& f2, const ValidatedVectorFunctionModel& f3) {
    return join(join(f1,f2),f3);
}

void Enclosure::affine_draw(CanvasInterface& canvas, const Projection2d& projection, uint accuracy) const {
    ARIADNE_ASSERT_MSG(Ariadne::subset(this->_reduced_domain,this->_domain),*this);

    ValidatedVectorFunctionModel cached_space_function=this->_space_function;
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

    ValidatedVectorFunctionModel fg=this->function_factory().create_zeros(this->dimension()+1u+this->number_of_constraints(),this->domain());
    for(uint i=0; i!=this->dimension(); ++i) { fg[i]=this->_space_function[i]; }
    for(uint i=0; i!=1; ++i) { fg[i+this->dimension()]=this->_time_function; }
    for(uint i=0; i!=this->_constraints.size(); ++i) { fg[i+this->dimension()+1]=this->_constraints[i].function(); }

//    ValidatedVectorFunctionModel fg=join(this->space_function(),this->time_function(),this->constraint_function());

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

const ValidatedScalarFunctionModel& repr(const ValidatedScalarFunctionModel& f) { return f; }
const ValidatedVectorFunctionModel& repr(const ValidatedVectorFunctionModel& f) { return f; }
const List<ValidatedScalarFunctionModel>& repr(const List<ValidatedScalarFunctionModel>& f) { return f; }

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
           << ", constraints=" << this->constraints()
           << ")";

    } return os;
}





/*
ValidatedAffineConstrainedImageSet
Enclosure::affine_approximation() const
{
    this->_check();
    typedef List<ValidatedScalarFunctionModel>::const_iterator const_iterator;

    const uint nx=this->dimension();
    const uint np=this->number_of_parameters();

    Enclosure set(*this);

    //if(set._zero_constraints.size()>0) { set._solve_zero_constraints(); }
    this->_check();

    Vector<Float> h(nx);
    Matrix<Float> G(nx,np);
    for(uint i=0; i!=nx; ++i) {
        ValidatedScalarFunctionModel component=set._space_function[i];
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
        const ValidatedScalarFunctionModel& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_inequality_constraint(a,b);
    }

    for(const_iterator iter=set._zero_constraints.begin();
            iter!=set._zero_constraints.end(); ++iter) {
        const ValidatedScalarFunctionModel& constraint=*iter;
        b=-constraint.model().value();
        for(uint j=0; j!=np; ++j) { a[j]=constraint.model().gradient(j); }
        result.new_equality_constraint(a,b);
    }
    return result;
}
*/

/*
struct ValidatedAffineModel {
    Float _c; Vector<Float> _g; Float _e;
    ValidatedAffineModel(Float c, const Vector<Float>& g, Float e) : _c(c), _g(g), _e(e) { }
};

ValidatedAffineModel _affine_model(const ValidatedTaylorModel& tm) {
    ValidatedAffineModel result(0.0,Vector<Float>(tm.argument_size(),0.0),tm.error());
    set_rounding_upward();
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
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

    Vector< ValidatedAffineModel > affine_function_models(space_function.result_size());
    for(uint i=0; i!=space_function.result_size(); ++i) { affine_function_models[i]=affine_model(space_function.models()[i]); }
    //affine_function_models[space_function.result_size()]=affine_model(time_function.model());

    ValidatedAffineConstrainedImageSet result(affine_function_models);
    //std::cerr<<"\n"<<*this<<"\n"<<result<<"\n\n";

    for(uint i=0; i!=this->number_of_constraints(); ++i) {
        ValidatedTaylorModel const& constraint_model=constraint_functions[i].model();
        ValidatedAffineModel affine_constraint_model=affine_model(constraint_model);
        Interval constraint_bound=this->constraint(i).bounds();
        result.new_constraint(constraint_bound.lower()<=affine_constraint_model<=constraint_bound.upper());
    }

    ARIADNE_LOG(2,"set="<<*this<<"\nset.affine_over_approximation()="<<result<<"\n");
    return result;

}


Enclosure product(const Enclosure& set, const Interval& ivl) {
    typedef List<ValidatedConstraintModel>::const_iterator const_iterator;

    ValidatedVectorFunctionModel new_function=combine(set.function(),set.function_factory().create_identity(ivl));

    Enclosure result(new_function.domain(),new_function,set.function_factory());
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(),embed(iter->function(),ivl),iter->upper_bound()));
    }
    result._time_function=embed(set._time_function,ivl);
    result._dwell_time_function=embed(set._dwell_time_function,ivl);

    return result;
}

Enclosure product(const Enclosure& set, const Box& bx) {
    typedef List<ValidatedConstraintModel>::const_iterator const_iterator;

    ValidatedVectorFunctionModel new_function=combine(set.function(),set.function_factory().create_identity(bx));

    Enclosure result(new_function.domain(),new_function,set.function_factory());
    for(const_iterator iter=set._constraints.begin(); iter!=set._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(),embed(iter->function(),bx),iter->upper_bound()));
    }
    result._time_function=embed(set._time_function,bx);
    result._dwell_time_function=embed(set._dwell_time_function,bx);

    return result;
}

Enclosure product(const Enclosure& set1, const Enclosure& set2) {
    ARIADNE_ASSERT(set1.time_function().range() == set2.time_function().range());
    ARIADNE_ASSERT(set1.dwell_time_function().range() == set2.dwell_time_function().range());

    typedef List<ValidatedConstraintModel>::const_iterator const_iterator;

    ValidatedVectorFunctionModel new_function=combine(set1.function(),set2.function());

    Enclosure result(new_function.domain(),new_function,set1.function_factory());
    for(const_iterator iter=set1._constraints.begin(); iter!=set1._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(),embed(iter->function(),set2.domain()),iter->upper_bound()));
    }
    for(const_iterator iter=set2._constraints.begin(); iter!=set2._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(),embed(set1.domain(),iter->function()),iter->upper_bound()));
    }
    result._time_function=embed(set1.time_function(),set2.time_function().domain());
    result._dwell_time_function=embed(set1.dwell_time_function(),set2.dwell_time_function().domain());

    return result;
}

Enclosure apply(const ValidatedVectorFunction& function, const Enclosure& set) {
    Enclosure result(set);
    result.apply_map(function);
    return result;
}

Enclosure unchecked_apply(const ValidatedVectorFunctionModel& function, const Enclosure& set) {
    Enclosure result(set);
    const ValidatedVectorFunctionModel& space_function=result.function();
    const_cast<ValidatedVectorFunctionModel&>(space_function)=Ariadne::unchecked_compose(function,set.function());
    return result;
}

} // namespace Ariadne


