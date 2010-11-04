/***************************************************************************
 *            hybrid_evolver-constrained.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include "numeric.h"
#include "vector.h"
#include "function.h"
#include "taylor_model.h"
#include "taylor_function.h"
#include "grid_set.h"
#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver-constrained.h"
#include "orbit.h"

#include "integrator.h"
#include "solver.h"

namespace {

ScalarTaylorFunction approximate_crossing_time(const ScalarTaylorFunction& tf, Interval normal) {
    ScalarFunction f(tf.polynomial());
    Vector<Interval> sdom=project(tf.domain(),range(0u,tf.argument_size()-1u));
    VectorTaylorFunction id=VectorTaylorFunction::identity(sdom);
    ScalarTaylorFunction h=ScalarTaylorFunction::constant(sdom,tf.domain()[f.argument_size()-1u].midpoint());
    for(unsigned int i=0; i!=8; ++i) {
        h=h-compose(f,join(id,h))/normal.midpoint();
        h.clobber();
    }
    return h;
}

} // namespace

namespace Ariadne {

DiscreteEvent final_time_event("_final_time_");
DiscreteEvent time_step_event("_step_size_");

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Interval> IntervalVector;




Orbit<HybridEnclosure>
ConstraintHybridEvolver::
orbit(const HybridAutomatonInterface& system,
      const HybridEnclosure& initial,
      const HybridTime& time,
      Semantics semantics) const
{
    EnclosureListType final, reachable, intermediate;
    this->_evolution(final,reachable,intermediate,system,initial,time,semantics,true);
    Orbit<HybridEnclosure> orbit(initial);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}


ConstraintHybridEvolver::ConstraintHybridEvolver()
    : _parameters(new EvolutionParametersType())
{ }

ConstraintHybridEvolver::ConstraintHybridEvolver(const EvolutionParametersType& parameters)
    : _parameters(new EvolutionParametersType(parameters))
{ }

void
ConstraintHybridEvolver::
_evolution(EnclosureListType& final,
           EnclosureListType& reachable,
           EnclosureListType& intermediate,
           HybridAutomatonInterface const& system,
           HybridEnclosure const& initial_set,
           HybridTime const& maximum_time,
           Semantics semantics,
           bool reach) const
{
    List<HybridEnclosure> working;
    List<HybridEnclosure> pending;

    pending.push_back(HybridEnclosure(initial_set));

    while(!working.empty() || !pending.empty()) {
        this->_upper_evolution_step(working,pending,final,reachable,intermediate,system,maximum_time);
    }
}

void
ConstraintHybridEvolver::
_upper_evolution_step(List<HybridEnclosure>& working_sets,
                      List<HybridEnclosure>& pending_sets,
                      ListSet<HybridEnclosure>& evolve_sets,
                      ListSet<HybridEnclosure>& reach_sets,
                      ListSet<HybridEnclosure>& intermediate_sets,
                      HybridAutomatonInterface const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    // working sets are those for which we do not have to check initial activation
    // pending sets are those for which we do; they are typically jump sets or initial sets
    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    TaylorIntegrator integrator(32,this->_parameters->flow_accuracy);
    Float maximum_step_size=this->_parameters->maximum_step_size;

    const Real maximum_time=maximum_hybrid_time.continuous_time();
    const uint maximum_steps=maximum_hybrid_time.discrete_time();

    // Take the first pending set and
    if(working_sets.empty()) {
        HybridEnclosure starting_set=pending_sets.back();
        pending_sets.pop_back();
        ARIADNE_LOG(4,"starting_set="<<starting_set);

        DiscreteLocation starting_location=starting_set.location();
        TaylorConstrainedImageSet starting_space_set=starting_set.continuous_state_set();

        HybridEnclosure flowable_set=starting_set;

        ARIADNE_LOG(4,"\n");
        Set<DiscreteEvent> invariant_events=system.invariant_events(starting_location);
        ARIADNE_LOG(4,"invariant_events:"<<invariant_events<<"\n");
        Set<DiscreteEvent> urgent_events=system.urgent_events(starting_location);
        ARIADNE_LOG(4,"urgent_events:"<<urgent_events<<"\n");
        Set<DiscreteEvent> permissive_events=system.permissive_events(starting_location);
        ARIADNE_LOG(4,"permissive_events:"<<permissive_events<<"\n");

        // Compute possibly initially active invariants
        for(Set<DiscreteEvent>::iterator iter=invariant_events.begin(); iter!=invariant_events.end(); ++iter) {
            DiscreteEvent event=*iter;
            ScalarFunction invariant_function=system.invariant_function(starting_location,event);
            if(possibly(starting_space_set.satisfies(invariant_function>=0))) {
                flowable_set.new_invariant(event,invariant_function);
            }
        }
        // Compute possibly initially active urgent events
        for(Set<DiscreteEvent>::iterator iter=urgent_events.begin(); iter!=urgent_events.end();  ++iter) {
            DiscreteEvent event=*iter;
            ScalarFunction guard_function=system.guard_function(starting_location,event);
            if(possibly(starting_space_set.satisfies(guard_function>=0))) {
                flowable_set.new_invariant(event,guard_function);
                HybridEnclosure jump_set=starting_set;
                jump_set.new_activation(event,guard_function);
                DiscreteLocation targent=system.target(starting_location,event);
                VectorFunction reset_function=system.reset_function(starting_location,event);
                jump_set.apply_reset(event,targent,reset_function);
                pending_sets.append(jump_set);
            }
        }
        working_sets.append(flowable_set);
    }

    HybridEnclosure starting_set=working_sets.back();
    working_sets.pop_back();

    HybridEnclosure jump_set(starting_set);

    List<DiscreteEvent> starting_events=starting_set.previous_events();
    DiscreteLocation starting_location=starting_set.location();
    TaylorConstrainedImageSet starting_space_set=starting_set.continuous_state_set();
    ScalarTaylorFunction starting_time=starting_set.time_function();

    if(verbosity==1) {
        ARIADNE_LOG(1,"\r"
                    <<"#w="<<std::setw(4)<<std::left<<working_sets.size()+1u
                    <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                    <<" s="<<std::setw(3)<<std::left<<starting_events.size()
                    //<<" t="<<std::setw(7)<<std::fixed<<starting_time.value()
                    <<" t=["<<std::setw(7)<<std::left<<std::fixed<<starting_time.range().lower()
                    <<","<<std::setw(7)<<std::left<<std::fixed<<starting_time.range().upper()<<"]"
                    <<" a="<<std::setw(2)<<std::left<<starting_set.continuous_state_set().domain().size()
                    <<" nc="<<std::setw(2)<<std::left<<starting_set.continuous_state_set().number_of_constraints()
                    <<" r="<<std::setw(7)<<starting_set.continuous_state_set().radius()
                    <<" l="<<std::left<<starting_location
                    <<" c="<<starting_set.continuous_state_set().centre()
                    <<" e="<<starting_events
                    <<"                      \n");
    }

    ARIADNE_LOG(4,"starting_events:"<<starting_events<<"\n");
    ARIADNE_LOG(4,"starting_location:"<<starting_location<<"\n");
    ARIADNE_LOG(4,"starting_space_set:"<<starting_space_set<<"\n");
    ARIADNE_LOG(4,"starting_time:"<<starting_time<<"\n");

    if(starting_time.argument_size()>5) {
        ARIADNE_ERROR("Too many independent parameters ("<<starting_time.argument_size()<<")\n"); return; }

    // Set the dimension
    const uint n=starting_set.dimension();

    // Extract mode and transitions, dynamic and constraints
    //DiscreteMode mode=system.mode(starting_location);
    //Set<DiscreteTransition> transitions=system.transitions(starting_location);

    VectorFunction dynamic=system.dynamic_function(starting_location);
    ARIADNE_LOG(4,"dynamic:"<<dynamic<<"\n");

    // Compute flow and actual time step size used
    IntervalVector flow_spacial_domain(starting_set.continuous_state_set().bounding_box());
    ARIADNE_LOG(4,"flow_spacial_domain="<<flow_spacial_domain<<"\n");
    VectorTaylorFunction flow_model=integrator.flow_step(dynamic,flow_spacial_domain,min(maximum_step_size,maximum_time));
    ARIADNE_LOG(4,"twosided_flow_model:"<<flow_model<<"\n");
    Box flow_domain=Box(flow_model.domain());
    Float step_size=flow_domain[n].upper();

    // Computed flowed set
    flow_domain[n]=Interval(0,step_size);
    flow_model=restrict(flow_model,flow_domain);
    ARIADNE_LOG(4,"flow_model:"<<flow_model<<"\n");
    ARIADNE_LOG(4,"flow_model.domain():"<<flow_model.domain()<<"\n");
    ARIADNE_LOG(4,"flow_model.range():"<<flow_model.range()<<"\n");
    Box flow_bounds=flow_model.range();
    VectorFunction flow=flow_model.function();
    ARIADNE_LOG(4,"flow:"<<flow<<"\n");
    ARIADNE_LOG(4,"starting_set:"<<starting_set<<"\n");
    HybridEnclosure reached_set=starting_set;
    reached_set.apply_flow_for(flow_model,flow_domain[n].upper());
    ARIADNE_LOG(4,"flowed_set:"<<reached_set<<"\n");

    ARIADNE_LOG(4,"\n");
    Set<DiscreteEvent> invariant_events=system.invariant_events(starting_location);
    ARIADNE_LOG(4,"invariant_events:"<<invariant_events<<"\n");
    Set<DiscreteEvent> urgent_events=system.urgent_events(starting_location);
    ARIADNE_LOG(4,"urgent_events:"<<urgent_events<<"\n");
    Set<DiscreteEvent> permissive_events=system.permissive_events(starting_location);
    ARIADNE_LOG(4,"permissive_events:"<<permissive_events<<"\n");

    // Compute possibly active invariants
    for(Set<DiscreteEvent>::iterator iter=invariant_events.begin(); iter!=invariant_events.end(); ) {
        if(compose(system.invariant_function(starting_location,*iter),flow_model).range().upper()<0.0) {
            invariant_events.erase(iter++);
        } else {
            ++iter;
        }
    }
    // Compute possibly active urgent events
    for(Set<DiscreteEvent>::iterator iter=urgent_events.begin(); iter!=urgent_events.end(); ) {
        if(compose(system.guard_function(starting_location,*iter),flow_model).range().upper()<0.0) {
            urgent_events.erase(iter++);
        } else {
            ++iter;
        }
    }
    // Compute possibly active permissive events
    for(Set<DiscreteEvent>::iterator iter=permissive_events.begin(); iter!=permissive_events.end(); ) {
        if(compose(system.guard_function(starting_location,*iter),flow_model).range().upper()<0.0) {
            permissive_events.erase(iter++);
        } else {
            ++iter;
        }
    }

    ARIADNE_LOG(4,"possibly_active_invariant_events:"<<invariant_events<<"\n");
    ARIADNE_LOG(4,"possibly_active_urgent_events:"<<urgent_events<<"\n");
    ARIADNE_LOG(4,"possibly_active_permissive_events:"<<permissive_events<<"\n");
    ARIADNE_LOG(4,"\n");

    Set<DiscreteEvent> blocking_events=join(invariant_events,urgent_events);
    Set<DiscreteEvent> transition_events=join(permissive_events,urgent_events);

    // Compute restrictions on continuous evolution due to invariants
    for(Set<DiscreteEvent>::iterator iter=invariant_events.begin(); iter!=invariant_events.end(); ++iter) {
        DiscreteEvent event=*iter;
        ARIADNE_LOG(4,"    event:"<<event<<", ");
        ScalarFunction constraint=system.invariant_function(starting_location,event);
        ARIADNE_LOG(4,"constraint:"<<constraint<<"\n");
        reached_set.new_invariant(event,constraint);
    }
    ARIADNE_LOG(4,"progress_set:"<<reached_set<<"\n");

    // Cache constraint functions for urgent events
    Map<DiscreteEvent,ScalarFunction> constraint_functions;
    for(Set<DiscreteEvent>::const_iterator iter=urgent_events.begin(); iter!=urgent_events.end(); ++iter) {
        constraint_functions.insert(*iter,system.guard_function(starting_location,*iter));
    }

    // Only compute events if discrete transitions still allowed
    if(starting_events.size()<maximum_steps) {

        // Compute crossing time interval and crossing time functions
        Map<DiscreteEvent,ScalarTaylorFunction> crossing_time_models;
        Map<DiscreteEvent,Interval> crossing_time_intervals;

        for(Set<DiscreteEvent>::const_iterator iter=urgent_events.begin(); iter!=urgent_events.end(); ++iter) {
            DiscreteEvent event=*iter;
            ScalarFunction guard=constraint_functions[event];
            // Test for transversality
            ScalarFunction guard_derivative=lie_derivative(guard,dynamic);
            Interval normal_derivative_1=guard_derivative(flow_bounds);
            Interval normal_derivative_2=ScalarTaylorFunction(flow_bounds,guard_derivative).range();
            Interval normal_derivative_3=compose(guard_derivative,flow_model).range();
            Interval normal_derivative=intersection(intersection(normal_derivative_1,normal_derivative_2),normal_derivative_3);
            if(normal_derivative.lower()>0.0 || normal_derivative.upper()<0.0) {
                // If transverse crossing, try to compute crossing time
                try {
                    // Compute crossing time as a function of the set's parameters (except for the last parameter, which is the current dwell time)
                    ScalarTaylorFunction guarded_evolution=compose(guard,jump_set.continuous_state_set().taylor_function());
                    ScalarTaylorFunction crossing_time_model=approximate_crossing_time(guarded_evolution,normal_derivative);
                    /*
                    IntervalVector evolution_domain=guarded_evolution.domain();
                    uint nsp=evolution_domain.size()-1u; // The number of parameters, excluding time
                    Vector<Interval> parameter_domain=project(evolution_domain,range(0u,nsp));
                    Interval time_domain=evolution_domain[nsp];
                    ScalarTaylorFunction crossing_time_model=
                        IntervalNewtonSolver(1e-8,12).implicit(ScalarFunction(guarded_evolution.polynomial()),flow_spacial_domain,time_domain);
                    */
                    jump_set.new_guard(event,guard,crossing_time_model);
                }
                catch(const ImplicitFunctionException& e) {
                    ARIADNE_WARN(""<<e.what()<<"\n");
                    jump_set.new_guard(event,guard);
                }
            } else {
                // If non-transverse crossing introduce guard as equation
                jump_set.new_guard(event,guard);
            }
            if(!definitely(jump_set.empty())) {
                DiscreteLocation target=system.target(starting_location,event);
                VectorFunction reset=system.reset_function(starting_location,event);
                jump_set.apply_reset(event,target,reset);
                working_sets.append(jump_set);
            }
        }

        // Compute the reached set under a single urgent event
        // TODO: Better to compute crossing time on the flow domain, and apply flow to starting set
        for(Set<DiscreteEvent>::const_iterator iter=urgent_events.begin(); iter!=urgent_events.end(); ++iter) {
            DiscreteEvent event=*iter;
            HybridEnclosure jump_set=reached_set;
            for(Map<DiscreteEvent,ScalarFunction>::const_iterator other_iter=constraint_functions.begin();
            other_iter!=constraint_functions.end(); ++other_iter)
            {
                if(other_iter->first!=event) { jump_set.new_invariant(other_iter->first,other_iter->second); }
            }
            ScalarFunction guard=constraint_functions[event];
            // Test for transversality
            ScalarFunction guard_derivative=lie_derivative(guard,dynamic);
            Interval normal_derivative_1=guard_derivative(flow_bounds);
            Interval normal_derivative_2=ScalarTaylorFunction(flow_bounds,guard_derivative).range();
            Interval normal_derivative_3=compose(guard_derivative,flow_model).range();
            Interval normal_derivative=intersection(intersection(normal_derivative_1,normal_derivative_2),normal_derivative_3);
            if(normal_derivative.lower()>0.0 || normal_derivative.upper()<0.0) {
                // If transverse crossing, try to compute crossing time
                try {
                    // Compute crossing time as a function of the set's parameters (except for the last parameter, which is the current dwell time)
                    ScalarTaylorFunction guarded_evolution=compose(guard,jump_set.continuous_state_set().taylor_function());
                    ScalarTaylorFunction crossing_time_model=approximate_crossing_time(guarded_evolution,normal_derivative);
                    /*
                    IntervalVector evolution_domain=guarded_evolution.domain();
                    uint nsp=evolution_domain.size()-1u; // The number of parameters, excluding time
                    Vector<Interval> parameter_domain=project(evolution_domain,range(0u,nsp));
                    Interval time_domain=evolution_domain[nsp];
                    ScalarTaylorFunction crossing_time_model=
                    IntervalNewtonSolver(1e-8,12).implicit(ScalarFunction(guarded_evolution.polynomial()),flow_spacial_domain,time_domain);
                    */
                    jump_set.new_guard(event,guard,crossing_time_model);
                }
                catch(const ImplicitFunctionException& e) {
                    ARIADNE_WARN(""<<e.what()<<"\n");
                    jump_set.new_guard(event,guard);
                }
            } else {
                // If non-transverse crossing introduce guard as equation
                jump_set.new_guard(event,guard);
            }
            if(!definitely(jump_set.empty())) {
                DiscreteLocation target=system.target(starting_location,event);
                VectorFunction reset=system.reset_function(starting_location,event);
                jump_set.apply_reset(event,target,reset);
                working_sets.append(jump_set);
            }
        }
    }


    // Compute restrictions on continuous evolution due to urgent events
    for(Set<DiscreteEvent>::const_iterator iter=urgent_events.begin(); iter!=urgent_events.end(); ++iter) {
        DiscreteEvent event=*iter;
        ARIADNE_LOG(4,"    event:"<<event<<", ");
        ScalarFunction constraint=system.guard_function(starting_location,event);
        ARIADNE_LOG(4,"constraint:"<<constraint<<"\n");
        reached_set.new_invariant(event,constraint);
    }
    ARIADNE_LOG(4,"reached_set:"<<reached_set<<"\n");

    if(starting_events.size()<maximum_steps) {
        // Compute the reached set under a single permissive event
        for(Set<DiscreteEvent>::const_iterator iter=permissive_events.begin(); iter!=permissive_events.end(); ++iter) {
            DiscreteEvent event=*iter;
            HybridEnclosure jump_set=reached_set;
            jump_set.new_guard(event,system.guard_function(starting_location,event));
            if(!definitely(jump_set.empty())) {
                DiscreteLocation target=system.target(starting_location,event);
                VectorFunction reset=system.reset_function(starting_location,event);
                jump_set.apply_reset(event,target,reset);
                working_sets.append(jump_set);
            }
        }
    }


    // Set the reached set, the next working set and the final set
    HybridEnclosure final_set(reached_set);
    HybridEnclosure progress_set(reached_set);

    reached_set.bound_time(maximum_time);
    ARIADNE_LOG(4,"reached_set:"<<reached_set<<"\n");
    final_set.set_time(maximum_time);
    ARIADNE_LOG(4,"final_set:"<<final_set<<"\n");
    ARIADNE_LOG(4,"final_set.empty():"<<final_set.empty()<<"\n");

    progress_set.set_step_time(step_size);
    ARIADNE_LOG(4,"progress_set:"<<progress_set<<"\n");

    reach_sets.adjoin(reached_set);
    ARIADNE_LOG(4,"reach_sets.size():"<<reach_sets.size()<<"\n");
    if(!definitely(final_set.empty())) { evolve_sets.adjoin(final_set); }
    ARIADNE_LOG(4,"evolve_sets.size():"<<evolve_sets.size()<<"\n");
    if(!definitely(progress_set.empty()) && progress_set._time.range().lower()<maximum_time) {
        working_sets.append(progress_set); }
    ARIADNE_LOG(4,"working_sets.size():"<<working_sets.size()<<"\n");

    ARIADNE_LOG(4,"DONE STEP\n"<<"\n");


    return;
}

} // namespace Ariadne
