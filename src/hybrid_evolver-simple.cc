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
#include "hybrid_evolver-simple.h"
#include "orbit.h"

#include "integrator.h"
#include "solver.h"
#include <boost/concept_check.hpp>

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

enum CrossingType { increasing=0, decreasing=1 };
inline std::ostream& operator<<(std::ostream& os, const CrossingType& cr) {
    switch(cr) { case increasing: os<<"incr"; break; case decreasing: os<<"decr"; break; } return os; }

} // namespace

namespace Ariadne {

template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Interval> IntervalVector;




Orbit<HybridEnclosure>
SimpleHybridEvolver::
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


SimpleHybridEvolver::SimpleHybridEvolver()
    : _parameters(new EvolutionParametersType())
{ }

SimpleHybridEvolver::SimpleHybridEvolver(const EvolutionParametersType& parameters)
    : _parameters(new EvolutionParametersType(parameters))
{ }

void
SimpleHybridEvolver::
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

    working.push_back(HybridEnclosure(initial_set));

    while(!working.empty()) {
        this->_upper_evolution_step(working,final,reachable,intermediate,system,maximum_time);
    }
}

struct TransitionData {
    DiscreteLocation target; ScalarFunction guard_function; VectorFunction reset_function;
    TransitionData() { }
    TransitionData(DiscreteLocation t, ScalarFunction g, VectorFunction r)
    : target(t), guard_function(g), reset_function(r) { }
};

HybridEnclosure
SimpleHybridEvolver::
_process_initial_events(List<HybridEnclosure>& working_sets,
                        HybridEnclosure const& starting_set,
                        Map<DiscreteEvent,TransitionData> const& transitions) const
{
    HybridEnclosure flowable_set=starting_set;
    Set<DiscreteEvent> events=transitions.keys();

    // Compute possibly initially active events
    for(Set<DiscreteEvent>::iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        DiscreteEvent event=*event_iter;
        ScalarFunction const& guard_function=transitions[event].guard_function;
        // Simple but inaccurate test if set is initially active
        if(possibly(starting_set.satisfies(guard_function>=0))) {
            HybridEnclosure immediate_jump_set(starting_set);
            immediate_jump_set.new_activation(event,guard_function);
            // Test if the active set is nonempty
            if(!definitely(immediate_jump_set.empty())) {
                DiscreteLocation const& target=transitions[event].target;
                VectorFunction const& reset_function=transitions[event].reset_function;
                immediate_jump_set.apply_reset(event,target,reset_function);
                ARIADNE_LOG(4,"  immediate_jump_set="<<immediate_jump_set<<"\n");
                working_sets.append(immediate_jump_set);
                flowable_set.new_invariant(event,guard_function);
            }
        }
    }
    return flowable_set;
}

VectorIntervalFunction
SimpleHybridEvolver::
_compute_flow(VectorFunction dynamic,
              Box const& initial_box,
              Float& step_size) const
{
    // Compute flow and actual time step size used
    TaylorIntegrator integrator(32,this->parameters().flow_accuracy);
    VectorIntervalFunction flow_model=integrator.flow(dynamic,initial_box,step_size);
    ARIADNE_LOG(6,"twosided_flow_model="<<flow_model<<"\n");
    IntervalVector flow_domain=flow_model.domain();
    step_size=flow_domain[flow_domain.size()-1u].upper();
    flow_domain[flow_domain.size()-1u]=Interval(0,step_size);
    flow_model=restrict(flow_model,flow_domain);
    ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
    return flow_model;
}

void
SimpleHybridEvolver::
_compute_active_events(Set<DiscreteEvent>& active_events,
                       VectorFunction const& dynamic,
                       Map<DiscreteEvent,ScalarFunction> const& guards,
                       VectorIntervalFunction const& flow,
                       HybridEnclosure const& initial_set) const
{
    Set<DiscreteEvent> events=guards.keys();
    HybridEnclosure reach_set=initial_set;
    IntervalVector flow_bounds=flow.range();
    reach_set.apply_flow(flow,flow.domain()[flow.domain().size()-1].upper());
    for(Set<DiscreteEvent>::iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        HybridEnclosure test_set=reach_set;
        test_set.new_activation(*event_iter,guards[*event_iter]);
        if(!definitely(test_set.empty())) {
            // Test direction of guard increase
            ScalarFunction flow_derivative = lie_derivative(guards[*event_iter],dynamic);
            Interval flow_derivative_range = flow_derivative(flow_bounds);
            if(flow_derivative_range.lower()<=0.0 && flow_derivative_range.upper()>=0.0) {
                ARIADNE_FAIL_MSG("Event "<<*event_iter<<" with guard "<<guards[*event_iter]<<" has normal derivative "<<flow_derivative<<" with range "<<flow_derivative_range); }
            if(flow_derivative_range.lower()>0.0) {
                // Compute crossing time
                active_events.insert(*event_iter);
            }
        }
    }
}


void
SimpleHybridEvolver::
_compute_crossing_times(Map<DiscreteEvent,ScalarIntervalFunction>& crossing_times,
                        Set<DiscreteEvent> const& active_events,
                        Map<DiscreteEvent,ScalarFunction> const& guards,
                        VectorIntervalFunction const& flow,
                       IntervalVector const& domain) const
{
    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        try {
            crossing_times[*event_iter]=implicit(compose(guards[*event_iter],flow));
        }
        catch(const ImplicitFunctionException& e) {
            ARIADNE_FAIL_MSG("Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
        }
    }

}


ScalarIntervalFunction
SimpleHybridEvolver::
_evolution_time(ScalarIntervalFunction const& maximum_evolution_time,
                Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times) const
{
    // Compute the evolution time for the current step given a maximum time compatible
    // with the flow, and crossing time functions for the transverse events
    return maximum_evolution_time;
}


void SimpleHybridEvolver::
_apply_time_step(HybridEnclosure& reach_set, HybridEnclosure& evolve_set,
                 ListSet<HybridEnclosure>& final_sets, List<HybridEnclosure>& jump_sets,
                 HybridEnclosure const& initial_set, Float const& final_time,
                 VectorTaylorFunction const& flow, ScalarIntervalFunction const& evolution_time,
                 Set<DiscreteEvent> const& active_events, Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    Interval time_domain=flow.domain()[flow.domain().size()-1];

    reach_set=initial_set;
    reach_set.apply_flow(flow,evolution_time);
    evolve_set=initial_set;
    evolve_set.apply_flow_step(flow,evolution_time);

    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        HybridEnclosure jump_set=initial_set;
        jump_set.apply_flow_step(flow,crossing_times[event]);
        for(Set<DiscreteEvent>::const_iterator other_event_iter=active_events.begin();
            other_event_iter!=active_events.end(); ++other_event_iter)
        {
            DiscreteEvent other_event=*other_event_iter;
            if(other_event!=event) {
                jump_set.new_invariant(other_event,transitions[other_event].guard_function);
            }
        }
        if(!definitely(jump_set.empty())) {
            jump_set.apply_reset(event,transitions[event].target,transitions[event].reset_function);
            jump_sets.append(jump_set);
        }
    }

    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        reach_set.new_invariant(event,transitions[event].guard_function);
        evolve_set.new_invariant(event,transitions[event].guard_function);
    }

    DiscreteEvent final_time_event("final_time");
    HybridEnclosure final_set=reach_set;
    final_set.set_time(final_time_event,final_time);
    evolve_set.new_time_bound(final_time_event,final_time);
    reach_set.new_time_bound(final_time_event,final_time);
    final_sets.adjoin(final_set);
}

void write_data(uint working_sets_size, uint reach_sets_size, const HybridEnclosure& initial_set)
{
    uint verbosity=1;
    Box initial_bounding_box=initial_set.bounding_box().continuous_state_set();
    Interval initial_time_range=initial_set.time_function().range();
    ARIADNE_LOG(1,"\r"
            <<"#w="<<std::setw(4)<<std::left<<working_sets_size+1u
            <<"#r="<<std::setw(4)<<std::left<<reach_sets_size
            <<" s="<<std::setw(3)<<std::left<<initial_set.previous_events().size()
            <<" t=["<<std::setw(7)<<std::left<<std::fixed<<initial_time_range.lower()
            <<","<<std::setw(7)<<std::left<<std::fixed<<initial_time_range.upper()<<"]"
            <<" a="<<std::setw(2)<<std::left<<initial_set.number_of_parameters()
            <<" nc="<<std::setw(2)<<std::left<<initial_set.number_of_constraints()
            <<" r="<<std::setw(7)<<initial_bounding_box.radius()
            <<" l="<<std::left<<initial_set.location()
            <<" c="<<initial_bounding_box.centre()
            <<" e="<<initial_set.previous_events()
            <<"                      \n");
}


void
SimpleHybridEvolver::
_upper_evolution_step(List<HybridEnclosure>& working_sets,
                      ListSet<HybridEnclosure>& final_sets,
                      ListSet<HybridEnclosure>& reach_sets,
                      ListSet<HybridEnclosure>& intermediate_sets,
                      HybridAutomatonInterface const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    ARIADNE_LOG(4,"\nSimpleHybridEvolver::_upper_evolution_step(...): verbosity="<<verbosity<<"\n");

    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    const Float maximum_time=maximum_hybrid_time.continuous_time();
    const uint maximum_steps=maximum_hybrid_time.discrete_time();

    // Routine check for emptiness
    if(working_sets.empty()) { return; }

    // Get the starting set for this round of evolution
    HybridEnclosure starting_set=working_sets.back(); working_sets.pop_back();
    ARIADNE_LOG(4,"\nstarting_set="<<starting_set<<"\n");

    // Test if maximum number of steps has been reached
    if(starting_set.previous_events().size()>maximum_steps) {
        final_sets.adjoin(starting_set); return;
    }

    // Extract starting location
    DiscreteLocation location=starting_set.location();

    // This evolver only works if there are only urgent events
    ARIADNE_ASSERT(system.invariant_events(location).empty());
    ARIADNE_ASSERT(system.permissive_events(location).empty());

    // Cache dynamic and constraint functions
    VectorFunction dynamic=system.dynamic_function(location);
    Set<DiscreteEvent> urgent_events=system.urgent_events(location);
    Map<DiscreteEvent,ScalarFunction> guard_functions;
    Map<DiscreteEvent,TransitionData> transitions;
    for(Set<DiscreteEvent>::const_iterator event_iter=urgent_events.begin(); event_iter!=urgent_events.end(); ++event_iter) {
        DiscreteLocation target=system.target(location,*event_iter);
        ScalarFunction guard_function=system.guard_function(location,*event_iter);
        VectorFunction reset_function=system.reset_function(location,*event_iter);
        guard_functions.insert(*event_iter,guard_function);
        transitions.insert(*event_iter,TransitionData(target,guard_function,reset_function));
    }

    ARIADNE_LOG(4,"\ndynamic="<<dynamic<<"\n");
    ARIADNE_LOG(4,"guards="<<guard_functions<<"\n");

    // Process the initially active events; cut out active points to leave initial flowable set.
    HybridEnclosure initial_set = this->_process_initial_events(working_sets,starting_set,transitions);
    ARIADNE_LOG(4,"\ninitial_set="<<initial_set<<"\n");

    while(!definitely(initial_set.empty())) {

        if(verbosity==1) { write_data(working_sets.size(),reach_sets.size(),initial_set); }

        // Compute the bounding box of the enclosure
        Box initial_bounding_box=initial_set.bounding_box().continuous_state_set();

        // Extract the initial time of the set
        ScalarIntervalFunction initial_time=initial_set.time_function();
        Interval initial_time_range=initial_time.range();

        ScalarIntervalFunction remaining_time=maximum_time-initial_time;

        // Compute flow and actual time step size used
        IntervalVector flow_spacial_domain=initial_bounding_box;
        Float step_size=this->parameters().maximum_step_size;
        VectorIntervalFunction flow_model=this->_compute_flow(dynamic,flow_spacial_domain,step_size);        ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
        Box flow_domain=flow_model.domain();
        Box flow_bounds=flow_model.range();
        ARIADNE_LOG(4,"flow_model.domain()="<<flow_domain<<" flow_model.range()="<<flow_bounds<<"\n");

        ScalarIntervalFunction evolution_time=ScalarIntervalFunction::constant(initial_bounding_box,step_size);

        ARIADNE_LOG(4,"flow_bounds="<<flow_bounds<<"\n");

        // Compute possibly active urgent events with increasing guards, and crossing times
        Set<DiscreteEvent> active_events;
        this->_compute_active_events(active_events,dynamic,guard_functions,flow_model,initial_set);

        // Compute crossing times
        Map<DiscreteEvent,ScalarTaylorFunction> crossing_times;
        this->_compute_crossing_times(crossing_times,active_events,guard_functions,flow_model,flow_domain);

        // Compute the evolution time
        evolution_time = this->_evolution_time(evolution_time,crossing_times);

        // Apply the time step
        HybridEnclosure reach_set(initial_set);
        HybridEnclosure evolve_set(initial_set);
        this->_apply_time_step(reach_set,evolve_set,final_sets,working_sets,initial_set,maximum_time,
                               flow_model,evolution_time,active_events,crossing_times,transitions);

        reach_sets.adjoin(reach_set);
        initial_set=evolve_set;
    }

}




ScalarIntervalFunction
DeterministicHybridEvolver::
_evolution_time(ScalarIntervalFunction const& maximum_evolution_time,
                Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times) const
{
    // Compute the evolution time for the current step given a maximum time compatible
    // with the flow, and crossing time functions for the transverse events
    Float step_size=maximum_evolution_time.range().upper();
    ScalarIntervalFunction evolution_time=maximum_evolution_time;
    for(Map<DiscreteEvent,ScalarIntervalFunction>::const_iterator time_iter=crossing_times.begin();
        time_iter!=crossing_times.end(); ++time_iter)
    {
        const ScalarIntervalFunction& crossing_time=time_iter->second;
        Float maximum_crossing_time=crossing_time.range().upper();
        maximum_crossing_time=max(maximum_crossing_time,2*step_size);
        ScalarTaylorFunction scaled_crossing_time=crossing_time/maximum_crossing_time;
        ScalarTaylorFunction creep_factor=scaled_crossing_time*(2-scaled_crossing_time);
        evolution_time=evolution_time*scaled_crossing_time;
    }
    return evolution_time;
}



void
DeterministicHybridEvolver::
_upper_evolution_step(List<HybridEnclosure>& working_sets,
                      ListSet<HybridEnclosure>& evolve_sets,
                      ListSet<HybridEnclosure>& reach_sets,
                      ListSet<HybridEnclosure>& intermediate_sets,
                      HybridAutomatonInterface const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    ARIADNE_LOG(4,"\nSimpleHybridEvolver::_upper_evolution_step(...): verbosity="<<verbosity<<"\n");

    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    Float maximum_step_size=this->parameters().maximum_step_size;

    const Float maximum_time=maximum_hybrid_time.continuous_time();
    const uint maximum_steps=maximum_hybrid_time.discrete_time();

    // Routine check for emptiness
    if(working_sets.empty()) { return; }

    // Get the starting set for this round of evolution
    HybridEnclosure starting_set=working_sets.back();
    working_sets.pop_back();
    ARIADNE_LOG(4,"\nstarting_set="<<starting_set<<"\n");

    // Check if maximum number of events has been reached
    if(starting_set.previous_events().size()>=maximum_steps) {
        evolve_sets.adjoin(starting_set);
        return;
    }

    // Extract starting location
    DiscreteLocation location=starting_set.location();
    List<DiscreteEvent> previous_events=starting_set.previous_events();

    // This evolver only works if there are only urgent events
    ARIADNE_ASSERT(system.invariant_events(location).empty());
    ARIADNE_ASSERT(system.permissive_events(location).empty());

    // Cache dimension, dynamic and constraint functions
    VectorFunction dynamic=system.dynamic_function(location);
    Set<DiscreteEvent> urgent_events=system.urgent_events(location);
    Map<DiscreteEvent,ScalarFunction> guard_functions;
    Map<DiscreteEvent,TransitionData> transitions;
    for(Set<DiscreteEvent>::const_iterator event_iter=urgent_events.begin(); event_iter!=urgent_events.end(); ++event_iter) {
        DiscreteLocation target=system.target(location,*event_iter);
        ScalarFunction guard_function=system.guard_function(location,*event_iter);
        VectorFunction reset_function=system.reset_function(location,*event_iter);
        guard_functions.insert(*event_iter,guard_function);
        transitions.insert(*event_iter,TransitionData(target,guard_function,reset_function));
    }

    ARIADNE_LOG(4,"\ndynamic:"<<dynamic<<"\n");
    ARIADNE_LOG(4,"guards:"<<guard_functions<<"\n");

    HybridEnclosure initial_set = this->_process_initial_events(working_sets,starting_set,transitions);
    ARIADNE_LOG(4,"\ninitial_set:"<<initial_set<<"\n");

   while(!definitely(initial_set.empty())) {
        //FIXME: Proper termination
        if(initial_set.time_function().range().lower()>=maximum_time) {
            ARIADNE_LOG(4,"\ninitial_time="<<initial_set.time_function()<<"\n");
            return;
        }

        // Compute the bounding box of the enclosure
        Box initial_bounding_box=initial_set.bounding_box().continuous_state_set();

        // Compute the bounding box of the enclosure
        ScalarIntervalFunction initial_time=initial_set.time_function();
        Interval initial_time_range=initial_time.range();

        if(verbosity==1) { write_data(working_sets.size(),reach_sets.size(),initial_set); }

        // Compute flow and actual time step size used
        IntervalVector flow_spacial_domain=initial_bounding_box;
        Float step_size=maximum_step_size;
        VectorIntervalFunction flow_model=this->_compute_flow(dynamic,flow_spacial_domain,step_size);        ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
        Box flow_domain=flow_model.domain();
        Box flow_bounds=flow_model.range();
        ARIADNE_LOG(4,"flow_model.domain()="<<flow_domain<<" flow_model.range()="<<flow_bounds<<"\n");

        ScalarIntervalFunction remaining_time=maximum_time-initial_time;
        ScalarIntervalFunction evolution_time=ScalarIntervalFunction::constant(flow_spacial_domain,step_size);
        Map<DiscreteEvent,ScalarTaylorFunction> crossing_times;

        // Apply flow to reach and evolve sets
        ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
        HybridEnclosure reach_set=initial_set;
        reach_set.apply_flow(flow_model,step_size);
        ARIADNE_LOG(4,"unconstrained_reach_set="<<reach_set<<"\n");
        HybridEnclosure evolve_set=initial_set;
        evolve_set.apply_flow_step(flow_model,step_size);
        ARIADNE_LOG(4,"unconstrained_evolve_set="<<evolve_set<<"\n");

        flow_bounds=intersection(flow_bounds,reach_set.bounding_box().continuous_state_set());
        ARIADNE_LOG(4,"flow_bounds="<<flow_bounds<<"\n");

        // Compute possibly active urgent events with increasing guards, and crossing times
        Set<DiscreteEvent> active_urgent_events;
        for(Set<DiscreteEvent>::iterator event_iter=urgent_events.begin(); event_iter!=urgent_events.end(); ++event_iter) {
            HybridEnclosure test_set=reach_set;
            test_set.new_activation(*event_iter,guard_functions[*event_iter]);
            if(!definitely(test_set.empty())) {
                // Test direction of guard increase
                ScalarFunction flow_derivative = lie_derivative(guard_functions[*event_iter],dynamic);
                Interval flow_derivative_range = flow_derivative(flow_bounds);
                if(flow_derivative_range.lower()<=0.0 && flow_derivative_range.upper()>=0.0) {
                    ARIADNE_FAIL_MSG("Event "<<*event_iter<<" with guard "<<guard_functions[*event_iter]<<" has normal derivative "<<flow_derivative<<" with range "<<flow_derivative_range); }
                if(flow_derivative_range.lower()>0.0) {
                    // Compute crossing time
                    active_urgent_events.insert(*event_iter);
                    try {
                        crossing_times[*event_iter]=implicit(compose(guard_functions[*event_iter],flow_model));
                    }
                    catch(const ImplicitFunctionException& e) {
                        ARIADNE_FAIL_MSG("Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                    }
                }
            }
        }

        ARIADNE_LOG(4,"\nevolution_time="<<evolution_time<<"\n");
        ARIADNE_LOG(4,"\ncrossing_times="<<crossing_times<<"\n");
        if(reach_set.time_function().range().upper()>=maximum_time) { ARIADNE_LOG(4,"possibly_last_time_step\n"); }
        ScalarIntervalFunction maximum_time_function=ScalarIntervalFunction::constant(initial_set.parameter_domain(),maximum_time);

        List<HybridEnclosure> jump_sets;
        if(!active_urgent_events.empty()) {
            // Compute restrictions on reachable set evolution due to invariants
            for(Set<DiscreteEvent>::iterator event_iter=active_urgent_events.begin(); event_iter!=active_urgent_events.end(); ++event_iter) {
                evolve_set.new_invariant(*event_iter,guard_functions[*event_iter]);
            }
            ARIADNE_LOG(4,"evolve_set="<<evolve_set<<"\n");
            ARIADNE_LOG(4,"evolve_set.empty()="<<evolve_set.empty()<<"\n");

            // Evolution in this location ends
            if(definitely(evolve_set.empty())) {
                this->_apply_time_step(reach_set,evolve_set,evolve_sets,working_sets,initial_set,maximum_time,
                                       flow_model,evolution_time,active_urgent_events,crossing_times,transitions);
            } else {
                // Evolution in this location has to be slowed to avoid partially active guards
                ARIADNE_LOG(1,"SLOWING_EVOLUTION\n");
                evolution_time=this->_evolution_time(evolution_time,crossing_times);
                ARIADNE_LOG(4,"evolution_time="<<evolution_time<<"\n");
                reach_set=initial_set;
                reach_set.apply_flow(flow_model,evolution_time);
                evolve_set=initial_set;
                evolve_set.apply_flow_step(flow_model,evolution_time);
            }
        }

        // Check for final time being reached
        DiscreteEvent final_time_event("end");
        HybridEnclosure final_set=reach_set;
        final_set.set_time(final_time_event,maximum_time);
        if( final_set.time_function().range().upper()>=maximum_time && !definitely(final_set.empty()) ) {
            reach_set.set_maximum_time(final_time_event,maximum_time);
            evolve_set.set_maximum_time(final_time_event,maximum_time);
        }

        // Insert computes sets into result
        reach_sets.adjoin(reach_set);
        if(!definitely(evolve_set.empty())) {
            intermediate_sets.adjoin(evolve_set);
        }
        if(!definitely(final_set.empty())) {
            evolve_sets.adjoin(final_set);
        }
        for(List<HybridEnclosure>::const_iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            working_sets.append(*jump_set_iter);
        }

        initial_set=evolve_set;
    } // End loop over initial set
}


} // namespace Ariadne
