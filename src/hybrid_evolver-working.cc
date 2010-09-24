/***************************************************************************
 *            hybrid_evolver-working.cc
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
#include "hybrid_evolver-working.h"
#include "orbit.h"

#include "integrator.h"
#include "solver.h"
#include <boost/concept_check.hpp>

namespace {

} // namespace

namespace Ariadne {

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

std::ostream& operator<<(std::ostream& os, const CrossingKind& crk) {
    switch(crk) {
        case DEGENERATE_CROSSING: os<<"degenerate"; break;
        case POSITIVE_CROSSING: os<<"positive"; break;
        case NEGATIVE_CROSSING: os<<"negative"; break;
        case INCREASING_CROSSING: os<<"increasing"; break;
        case DECREASING_CROSSING: os<<"decreasing"; break;
        case CONVEX_CROSSING: os<<"convex"; break;
        case CONCAVE_CROSSING: os<<"concave"; break;
        default: os << "unknown"; break;
    } return os;
}

std::ostream& operator<<(std::ostream& os, const StepKind& crk) {
    switch(crk) {
        case FULL_STEP: os<<"full"; break;
        case CREEP_STEP: os<<"creep"; break;
        case UNWIND_STEP: os<<"unwind"; break;
        case FINAL_STEP: os<<"final"; break;
        default: os << "unknown"; break;
    } return os;
}

std::ostream& operator<<(std::ostream& os, const TransitionData& transition) {
    return os << "kind="<<transition.event_kind<<", guard="<<transition.guard_function<<", "
                 "target="<<transition.target<<", reset="<<transition.reset_function;
}

std::ostream& operator<<(std::ostream& os, const TimingData& timing) {
    os << "step_kind="<<timing.step_kind<<", step_size="<<timing.step_size<<", "
       << "final_time="<<timing.final_time;
    if(timing.step_kind==CREEP_STEP) {
        os <<", spacial_evolution_time="<<timing.spacial_evolution_time.polynomial();
    } else if(timing.step_kind==UNWIND_STEP) {
        os << ", finishing_time="<<timing.finishing_time.polynomial();
    }
    os <<", evolution_time="<<timing.evolution_time.polynomial();
    return os;
}

std::ostream& operator<<(std::ostream& os, const CrossingData& crossing_data) {
    os << "kind="<<crossing_data.crossing_kind;
    if(crossing_data.crossing_kind==INCREASING_CROSSING) {
        os << ", crossing_time="<<crossing_data.crossing_time.polynomial();
    }
    if(crossing_data.crossing_kind==CONCAVE_CROSSING) {
        os << ", critical_time="<<crossing_data.critical_time.polynomial();
    }
    return os;
}

bool is_blocking(EventKind evk) {
    switch(evk) {
        case INVARIANT: case PROGRESS: case URGENT: case IMPACT:
            return true;
        case PERMISSIVE:
            return false;
        default:
            ARIADNE_FAIL_MSG("EventKind "<<evk<<" not recognised by is_blocking(...) predicate.");
    }
}

bool is_activating(EventKind evk) {
    switch(evk) {
        case PERMISSIVE: case URGENT: case IMPACT:
            return true;
        case INVARIANT: case PROGRESS:
            return false;
        default:
            ARIADNE_FAIL_MSG("EventKind "<<evk<<" not recognised by is_activating(...) predicate.");
    }
}

Set<DiscreteEvent> blocking_events(const Map<DiscreteEvent,TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        if(is_blocking(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}

Set<DiscreteEvent> activating_events(const Map<DiscreteEvent,TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        if(is_activating(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}





Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridAutomatonInterface& system,
      const HybridEnclosure& initial,
      const HybridTime& time,
      Semantics semantics) const
{
    ARIADNE_LOG(2,"\nHybridEvolverBase::orbit(...): verbosity="<<verbosity<<"\n");

    EvolutionData evolution_data;
    evolution_data.working_sets.push_back(HybridEnclosure(initial));
    while(!evolution_data.working_sets.empty()) {
        this->_upper_evolution_flow(evolution_data,system,time);
    }
    ARIADNE_ASSERT(evolution_data.working_sets.empty());
    ARIADNE_ASSERT(evolution_data.starting_sets.empty());

    Orbit<HybridEnclosure> orbit(initial);
    orbit.adjoin_intermediate(ListSet<HybridEnclosure>(evolution_data.intermediate_sets));
    orbit.adjoin_reach(evolution_data.reach_sets);
    orbit.adjoin_final(evolution_data.evolve_sets);
    return orbit;
}


HybridEvolverBase::HybridEvolverBase()
    : _parameters(new EvolutionParametersType())
{ }

HybridEvolverBase::HybridEvolverBase(const EvolutionParametersType& parameters)
    : _parameters(new EvolutionParametersType(parameters))
{ }

void
HybridEvolverBase::
_evolution(ListSet<HybridEnclosure>& final,
           ListSet<HybridEnclosure>& reachable,
           ListSet<HybridEnclosure>& intermediate,
           HybridAutomatonInterface const& system,
           HybridEnclosure const& initial_set,
           HybridTime const& maximum_time,
           Semantics semantics,
           bool reach) const
{
    EvolutionData evolution_data;

    evolution_data.working_sets.push_back(HybridEnclosure(initial_set));

    while(!evolution_data.working_sets.empty()) {
        this->_upper_evolution_flow(evolution_data,system,maximum_time);
    }

    final=evolution_data.evolve_sets;
    reachable=evolution_data.reach_sets;
    intermediate=evolution_data.intermediate_sets;
}

void
HybridEvolverBase::
_log_summary(uint ws, uint rs, HybridEnclosure const& starting_set) const
{
    Box starting_bounding_box=starting_set.space_bounding_box();
    Interval starting_time_range=starting_set.time_range();
    ARIADNE_LOG(1,"\r"
            <<"#w="<<std::setw(4)<<std::left<<ws+1u
            <<"#r="<<std::setw(4)<<std::left<<rs
            <<"#e="<<std::setw(3)<<std::left<<starting_set.previous_events().size()
            <<" t=["<<std::setw(7)<<std::left<<std::fixed<<starting_time_range.lower()
            <<","<<std::setw(7)<<std::left<<std::fixed<<starting_time_range.upper()<<"]"
            <<" #p="<<std::setw(2)<<std::left<<starting_set.number_of_parameters()
            <<" #c="<<std::setw(2)<<std::left<<starting_set.number_of_constraints()
            <<" r="<<std::setw(7)<<starting_bounding_box.radius()
            <<" c="<<starting_bounding_box.centre()
            <<" l="<<std::left<<starting_set.location()
            <<" e="<<starting_set.previous_events()
            <<"                      \n");
}

Map<DiscreteEvent,TransitionData>
HybridEvolverBase::
_extract_transitions(DiscreteLocation const& location,
                     HybridAutomatonInterface const& system) const
{
    Map<DiscreteEvent,TransitionData> transitions;
    Set<DiscreteEvent> events = system.urgent_events(location);
    for(Set<DiscreteEvent>::const_iterator event_iter=events.begin();
        event_iter!=events.end(); ++event_iter)
    {
        DiscreteLocation target=system.target(location,*event_iter);
        EventKind event_kind=system.event_kind(location,*event_iter);
        ScalarFunction guard_function=system.guard_function(location,*event_iter);
        VectorFunction reset_function=system.reset_function(location,*event_iter);
        TransitionData transition_data={event_kind,guard_function,target,reset_function};
        transitions.insert(*event_iter,transition_data);
    }
    return transitions;
}

void
HybridEvolverBase::
_process_initial_events(EvolutionData& evolution_data,
                        HybridEnclosure const& initial_set,
                        Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG(2,"HybridEvolverBase::_process_initial_events(...)\n");
    ARIADNE_ASSERT(evolution_data.starting_sets.empty());
    HybridEnclosure invariant_set=initial_set;

    // Apply restrictions due to invariants
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event=transition_iter->first;
        TransitionData const & transition=transition_iter->second;
        if(transition.event_kind==INVARIANT) {
            if (possibly(initial_set.satisfies(transition.guard_function>=0))) {
                invariant_set.new_invariant(event,transition.guard_function);
            }
        }
    }

    // Set the flowable set, storing the invariant set as a base for jumps
    HybridEnclosure flowable_set = invariant_set;

    // Compute possibly initially active events
    Set<DiscreteEvent> events=transitions.keys();
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event=transition_iter->first;
        TransitionData const & transition=transition_iter->second;
        if(transition.event_kind!=INVARIANT) {
            if(possibly(initial_set.satisfies(transition.guard_function>=0))) {
                if(transition.event_kind!=PERMISSIVE) {
                    HybridEnclosure immediate_jump_set=invariant_set;
                    immediate_jump_set.new_activation(event,transition.guard_function);
                    if(!definitely(immediate_jump_set.empty())) {
                        // Put the immediate jump set in the reached sets, since it does not occur in the flowable set
                        evolution_data.reach_sets.append(immediate_jump_set);
                        immediate_jump_set.apply_reset(event,transition.target,transition.reset_function);
                        ARIADNE_LOG(4,"immediate_jump_set="<<immediate_jump_set<<"\n");
                        evolution_data.intermediate_sets.append(immediate_jump_set);
                        evolution_data.working_sets.append(immediate_jump_set);
                    }
                }
                if(transition_iter->second.event_kind!=PROGRESS) {
                    flowable_set.new_invariant(event,transition.guard_function);
                }
            }
        }
    }

    // Put the flowable set in the starting sets for ordinary evolution
    if(!definitely(flowable_set.empty())) {
        evolution_data.starting_sets.append(flowable_set);
    }
}

VectorIntervalFunction
HybridEvolverBase::
_compute_flow(VectorFunction dynamic,
              Box const& initial_box,
              const Float& maximum_step_size) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_flow(...)\n");
    // Compute flow and actual time step size used
    TaylorIntegrator integrator(32,this->parameters().flow_accuracy);
    VectorIntervalFunction flow_model=integrator.flow_step(dynamic,initial_box,maximum_step_size);
    ARIADNE_LOG(6,"twosided_flow_model="<<flow_model<<"\n");
    IntervalVector flow_domain=flow_model.domain();
    Float step_size=flow_domain[flow_domain.size()-1u].upper();
    flow_domain[flow_domain.size()-1u]=Interval(0,step_size);
    flow_model=restrict(flow_model,flow_domain);
    ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
    return flow_model;
}

Set<DiscreteEvent>
HybridEvolverBase::
_compute_active_events(VectorFunction const& dynamic,
                       Map<DiscreteEvent,ScalarFunction> const& guards,
                       VectorIntervalFunction const& flow,
                       HybridEnclosure const& starting_set) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_active_events(...)\n");
    Set<DiscreteEvent> events=guards.keys();
    Set<DiscreteEvent> active_events;
    HybridEnclosure reach_set=starting_set;
    IntervalVector flow_bounds=flow.range();
    reach_set.apply_flow(flow,flow.domain()[flow.domain().size()-1].upper());
    for(Set<DiscreteEvent>::iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        const DiscreteEvent event=*event_iter;
        const ScalarFunction& guard_function=guards[event];
        ARIADNE_LOG(8,"event="<<event<<", guard="<<guard_function<<", flow_derivative="<<lie_derivative(guard_function,dynamic)<<"\n");
        // First try a simple test based on the bounding box
        if(guard_function(reach_set.space_bounding_box()).upper()>=0.0) {
            // Now make a set containing the complement of the constraint,
            // and test for emptiness. If the set is empty, then the guard is
            // not satisfied anywhere.
            HybridEnclosure test_set=reach_set;
            test_set.new_activation(event,guard_function);
            if(!definitely(test_set.empty())) {
                // Test direction of guard increase
                ScalarFunction flow_derivative = lie_derivative(guard_function,dynamic);
                Interval flow_derivative_range = flow_derivative(flow_bounds);
                if(flow_derivative_range.upper()>0.0) {
                    active_events.insert(*event_iter);
                }
            }
        }
    }
    return active_events;
}


Map<DiscreteEvent,CrossingData>
HybridEvolverBase::
_compute_crossings(Set<DiscreteEvent> const& active_events,
                   VectorFunction const& dynamic,
                   Map<DiscreteEvent,ScalarFunction> const& guards,
                   VectorIntervalFunction const& flow,
                   HybridEnclosure const& initial_set) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_crossings(...)\n");
    Map<DiscreteEvent,CrossingData> crossing_data;
    Box flow_bounds=flow.range();
    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        const DiscreteEvent event=*event_iter;
        ScalarFunction const& guard=guards[event];
        ScalarFunction derivative=lie_derivative(guard,dynamic);
        Interval derivative_range=derivative.evaluate(flow_bounds);
        if(derivative_range.lower()>0.0) {
            ScalarIntervalFunction crossing_time;
            try {
                crossing_time=implicit(compose(guard,flow));
                crossing_data[event]=CrossingData(INCREASING_CROSSING,crossing_time);
            }
            catch(const ImplicitFunctionException& e) {
                ARIADNE_LOG(0,"Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                crossing_data[event]=CrossingData(CONVEX_CROSSING);
            }
        } else if(derivative_range.upper()<0.0) {
            crossing_data[event]=CrossingData(DECREASING_CROSSING);
        } else {
            ScalarFunction second_derivative=lie_derivative(derivative,dynamic);
            Interval second_derivative_range=second_derivative.evaluate(flow_bounds);
            if(second_derivative_range.lower()>0.0) {
                crossing_data[event]=CrossingData(CONVEX_CROSSING);
            } else if(second_derivative_range.upper()<0.0) {
                try {
                    ScalarIntervalFunction critical_time=implicit(compose(derivative,flow));
                    crossing_data[event]=CrossingData(CONCAVE_CROSSING);
                    crossing_data[event].critical_time=critical_time;
                }
                catch(const ImplicitFunctionException& e) {
                    ARIADNE_LOG(0,"Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                    crossing_data[event]=CrossingData(DEGENERATE_CROSSING);
                }
            } else {
                crossing_data[event]=CrossingData(DEGENERATE_CROSSING);
            }
        }
    }
    return crossing_data;
}




TimingData
HybridEvolverBase::
_compute_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                VectorIntervalFunction const& flow,
                Map<DiscreteEvent,CrossingData> const& crossings,
                HybridEnclosure const& initial_set) const
{
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_timing(...)\n");
    TimingData result;
    result.step_size=flow.domain()[flow.domain().size()-1].upper();
    result.final_time=final_time;
    result.time_domain=Interval(0.0,result.step_size);
    result.time_coordinate=ScalarIntervalFunction::coordinate(Vector<Interval>(1u,result.time_domain),0u);
    result.remaining_time=result.final_time-initial_set.time_function();
    Interval starting_time_range=initial_set.time_function().range();
    Interval remaining_time_range=result.remaining_time.range();
    // NOTE: The time function may be negative or greater than the final time
    // over part of the parameter domain.
    if(remaining_time_range.upper()<=result.step_size) {
        result.step_kind=FINAL_STEP;
        result.finishing_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),numeric_cast<Interval>(final_time));
        result.evolution_time=result.final_time-initial_set.time_function();
    } else if(remaining_time_range.lower()<=result.step_size) {
        result.step_kind=UNWIND_STEP;
        if(remaining_time_range.width()<result.step_size) {
            Float constant_finishing_time=result.final_time-remaining_time_range.upper()+result.step_size;
            result.finishing_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),constant_finishing_time);
        } else {
            // FIXME: The finishing time may need to be adjusted
            result.finishing_time=0.5*(result.step_size+initial_set.time_function());
        }
        result.evolution_time=result.finishing_time-initial_set.time_function();
    } else {
        ARIADNE_LOG(8,"\ntesting for partially active events\n");
        bool creep=false;
        ScalarIntervalFunction evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
        for(Map<DiscreteEvent,CrossingData>::const_iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); ++crossing_iter)
        {
            if(crossing_iter->second.crossing_kind==INCREASING_CROSSING) {
                const ScalarIntervalFunction& crossing_time=crossing_iter->second.crossing_time;
                // Test crossing time incomparible with evolution time
                if(contains((evolution_time-crossing_time).range(),0.0)) {
                    // FIXME: Update evolution time better
                    evolution_time=evolution_time*(crossing_time/result.step_size)*(1.0-crossing_time/(4*result.step_size));
                    creep=true;
                }
            }
        }
        if(creep==true) {
            result.step_kind=CREEP_STEP;
            result.spacial_evolution_time=evolution_time;
            result.evolution_time=compose(result.spacial_evolution_time,initial_set.space_function());
        } else {
            result.step_kind=FULL_STEP;
            result.spacial_evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
            result.evolution_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),result.step_size);
        }
    }
    return result;
}







void
_apply_guard(HybridEnclosure& set,
             const VectorIntervalFunction& starting_state,
             const VectorIntervalFunction& flow,
             const ScalarIntervalFunction& elapsed_time,
             const DiscreteEvent event,
             const ScalarFunction& guard_function,
             const CrossingData crossing_data,
             const Semantics semantics)
{
    static const uint SUBDIVISIONS_FOR_DEGENERATE_CROSSING = 2;
    switch(crossing_data.crossing_kind) {
        case INCREASING_CROSSING:
            //set.new_state_constraint(event, guard_function <= 0.0);
            set.new_invariant(event, guard_function);
            // Alternatively:
            // set.new_parameter_constraint(event, elapsed_time <= compose(crossing_data.crossing_time,starting_state) );
            break;
        case CONVEX_CROSSING:
            //set.new_state_constraint(event, guard_function <= 0.0);
            set.new_invariant(event, guard_function);
            break;
        case CONCAVE_CROSSING: {
            ScalarIntervalFunction critical_time = compose(crossing_data.critical_time,starting_state);
            ScalarIntervalFunction final_guard
                = compose( guard_function, compose( flow, join(starting_state, elapsed_time) ) );
            ScalarIntervalFunction maximal_guard
                = compose( guard_function, compose( flow, join(starting_state, critical_time) ) );
            // If no points in the set arise from trajectories which will later leave the progress set,
            // then we only need to look at the maximum value of the guard.
            HybridEnclosure eventually_hitting_set=set;
            eventually_hitting_set.new_parameter_constraint( event, elapsed_time <= critical_time );
            eventually_hitting_set.new_parameter_constraint( event, maximal_guard >= 0.0);
            if(definitely(eventually_hitting_set.empty())) {
                set.new_parameter_constraint(event, maximal_guard <= 0.0);
                break;
            }
            // If no points in the set arise from trajectories which leave the progress set and
            // later return, then we only need to look at the guard at the final value
            HybridEnclosure returning_set=set;
            returning_set.new_parameter_constraint( event, elapsed_time >= critical_time );
            returning_set.new_parameter_constraint( event, final_guard <= 0.0 );
            if(definitely(returning_set.empty())) {
                set.new_parameter_constraint(event, final_guard <= 0.0);
                break;
            }
            // The exact guard set requires two components. Since this is not currently supported,
            // we approximate based on the semantics
            // TODO: Support splitting the set.
            switch(semantics) {
                case UPPER_SEMANTICS: set.new_parameter_constraint(event,final_guard<=0.0); break;
                case LOWER_SEMANTICS: set.new_parameter_constraint(event,maximal_guard<=0.0); break;
                default: assert(false);
            }
            break;
            // Code below is always exact, but uses two sets
            // set1.new_parameter_constraint(event,final_guard <= 0);
            // set2.new_parameter_constraint(event, elapsed_time <= critical_time);
            // set1.new_parameter_constraint(event, maximal_guard <= 0.0);
            // set2.new_parameter_constraint(event, elapsed_time >= critical_time);
        }
        case DEGENERATE_CROSSING: {
            // The crossing with the guard set is not one of the kinds handles above.
            // We obtain an over-appproximation by testing at finitely many time points
            // TODO: Handle lower semantics
            ARIADNE_ASSERT(semantics==UPPER_SEMANTICS);
            const uint n=SUBDIVISIONS_FOR_DEGENERATE_CROSSING;
            for(uint i=0; i!=n; ++i) {
                Float alpha=Float(i+1)/n;
                ScalarIntervalFunction intermediate_guard
                    = compose( guard_function, compose( flow, join(starting_state, alpha*elapsed_time) ) );
                set.new_parameter_constraint(event, intermediate_guard <= 0);
            }
            // FIXME: Lower semantics
            break;
        }
        case POSITIVE_CROSSING:
        case NEGATIVE_CROSSING:
        case DECREASING_CROSSING:
            break;
    }
}



DiscreteEvent step_time_event("step_time");

void GeneralHybridEvolver::
_apply_time_step(EvolutionData& evolution_data,
                 HybridEnclosure const& starting_set,
                 VectorIntervalFunction const& flow,
                 TimingData const& timing_data,
                 Map<DiscreteEvent,CrossingData> const& crossing_data,
                 VectorFunction const& dynamic,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    // FIXME: Make semantics function argument
    static const Semantics semantics = UPPER_SEMANTICS;

    ARIADNE_LOG(2,"GeneralHybridEvolver::_apply_time_step(...)\n");
    //ARIADNE_LOG(4,"evolution_time="<<evolution_time.range()<<" final_time="<<final_time<<"\n");

    ScalarIntervalFunction critical_function;

    // Compute events enabling transitions and events blocking continuous evolution
    Set<DiscreteEvent> active_events=crossing_data.keys();
    Set<DiscreteEvent> transition_events;
    Set<DiscreteEvent> blocking_events;
    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        switch(transitions[*event_iter].event_kind) {
            case PERMISSIVE:
                transition_events.insert(*event_iter);
                break;
            case URGENT: case IMPACT:
                transition_events.insert(*event_iter);
            case INVARIANT: case PROGRESS:
                blocking_events.insert(*event_iter);
        }
    }
    ARIADNE_LOG(6,"transition_events="<<transition_events<<"\n");
    ARIADNE_LOG(6,"blocking_events="<<blocking_events<<"\n");

    // Compute reach and evolve sets
    // TODO: This should probably be a separate function
    HybridEnclosure reach_set=starting_set;
    HybridEnclosure evolve_set=starting_set;
    HybridEnclosure final_set=starting_set;

    // Apply the flow to the reach and evolve sets
    switch(timing_data.step_kind) {
        case FULL_STEP:
            reach_set.apply_flow(flow,timing_data.step_size);
            evolve_set.apply_flow_step(flow,timing_data.step_size);
            final_set.apply_flow_and_set_time(flow,timing_data.final_time);
            break;
        case CREEP_STEP:
            // TODO: Use parameterised time here
            reach_set.apply_flow(flow,timing_data.spacial_evolution_time);
            evolve_set.apply_flow_step(flow,timing_data.spacial_evolution_time);
            final_set.apply_flow_and_set_time(flow,timing_data.final_time);
            break;
        case UNWIND_STEP:
            reach_set.apply_flow_and_bound_time(flow,timing_data.finishing_time);
            evolve_set.apply_flow_and_set_time(flow,timing_data.finishing_time);
            final_set.apply_flow_and_set_time(flow,timing_data.final_time);
            break;
        case FINAL_STEP:
            reach_set.apply_flow_and_bound_time(flow,timing_data.final_time);
            evolve_set.apply_flow_and_set_time(flow,timing_data.final_time);
            final_set.apply_flow_and_set_time(flow,timing_data.final_time);
            break;
    }
    ARIADNE_LOG(6,"flowed_set="<<reach_set<<"\n");

    // Apply evolution time constraints to reach and final sets
    DiscreteEvent step_event("_h_");
    VectorIntervalFunction starting_state=starting_set.space_function();
    VectorIntervalFunction embedded_starting_state=embed(starting_set.space_function(),timing_data.time_domain);
    ScalarIntervalFunction embedded_reach_time=embed(starting_set.parameter_domain(),timing_data.time_coordinate);
    ScalarIntervalFunction embedded_evolution_time=embed(timing_data.evolution_time,timing_data.time_domain);
    ScalarIntervalFunction reach_step_time=embedded_reach_time;
    ARIADNE_LOG(8,"embedded_reach_time="<<embedded_reach_time<<"\n")
    ARIADNE_LOG(8,"embedded_evolution_time="<<embedded_evolution_time<<"\n")
    reach_set.new_parameter_constraint(step_event,embedded_reach_time<=embedded_evolution_time);
    final_set.new_parameter_constraint(step_event,timing_data.remaining_time<=timing_data.evolution_time);

    // Apply constraints on reach, evolve and final sets due to events
    for(Set<DiscreteEvent>::const_iterator event_iter=blocking_events.begin();
        event_iter!=blocking_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        const ScalarFunction& event_guard_function=transitions[event].guard_function;
        const CrossingData& event_crossing_data=crossing_data[event];
        _apply_guard(reach_set,embedded_starting_state,flow,embedded_reach_time,
                     event,event_guard_function,event_crossing_data,semantics);
        _apply_guard(evolve_set,starting_state,flow,timing_data.evolution_time,
                     event,event_guard_function,event_crossing_data,semantics);
        _apply_guard(final_set,starting_state,flow,timing_data.remaining_time,
                     event,event_guard_function,event_crossing_data,semantics);
    }

    // Apply final timing constraints to reach and evolve sets
    DiscreteEvent final_event("_tmax_");
    reach_set.new_parameter_constraint(final_event,reach_set.time_function()<=timing_data.final_time);
    evolve_set.new_parameter_constraint(final_event,evolve_set.time_function()<=timing_data.final_time);

    // Compute final set, and apply constraints to reach and evolve set
    // if necessary
    if(timing_data.step_kind==FINAL_STEP || !definitely(final_set.empty())) {
        ARIADNE_LOG(4,"final_set="<<final_set<<"\n");
    }
    if(!definitely(final_set.empty())) {
        evolution_data.evolve_sets.append(final_set);
    }

    if(timing_data.step_kind!=FINAL_STEP) {
        ARIADNE_LOG(4,"evolve_set="<<evolve_set<<"\n");
        if(!definitely(evolve_set.empty())) {
            evolution_data.intermediate_sets.append(evolve_set);
            evolution_data.working_sets.append(evolve_set);
        }
    }

    ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
    evolution_data.reach_sets.append(reach_set);


    // Compute jump sets
    for(Set<DiscreteEvent>::const_iterator event_iter=transition_events.begin();
        event_iter!=transition_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        TransitionData const & transition = transitions[event];

        // Compute active set
        HybridEnclosure jump_set=starting_set;

        // Compute flow to guard set up to evolution time.
        ScalarIntervalFunction jump_step_time;
        VectorIntervalFunction jump_starting_state;
        switch(transitions[event].event_kind) {
            case PERMISSIVE:
                jump_set.apply_flow_and_bound_time(flow,timing_data.evolution_time);
                jump_set.new_activation(event,transition.guard_function);
                break;
            case URGENT: case IMPACT:
                switch(crossing_data[event].crossing_kind) {
                    case INCREASING_CROSSING:
                        jump_starting_state=starting_state; jump_step_time=compose(crossing_data[event].crossing_time,starting_set.space_function());
                        jump_set.new_parameter_constraint(step_event,jump_step_time<=timing_data.evolution_time);
                        jump_set.apply_flow_step(flow,crossing_data[event].crossing_time);
                        break;
                    case CONVEX_CROSSING:
                        jump_set.apply_flow_and_bound_time(flow,timing_data.evolution_time);
                        jump_set.new_guard(event,transition.guard_function);
                        jump_step_time=reach_step_time;
                        jump_starting_state=embedded_starting_state;
                        break;
                    case CONCAVE_CROSSING:
                    case DEGENERATE_CROSSING: // Just check positive derivative in this case; NOT EXACT
                        jump_set.apply_flow_and_bound_time(flow,timing_data.evolution_time);
                        jump_set.new_guard(event,transition.guard_function);
                        jump_set.new_invariant(event,lie_derivative(transition.guard_function,dynamic));
                        jump_step_time=reach_step_time;
                        jump_starting_state=embedded_starting_state;
                        break;
                    case NEGATIVE_CROSSING:
                    case POSITIVE_CROSSING:
                    case DECREASING_CROSSING:
                        ARIADNE_ASSERT(false);
                    default:
                        assert(false); // No remaining cases
                }
                break;
                // FIXME: Switch on step kind
            default:
                ARIADNE_FAIL_MSG("Invalid event kind "<<transitions[event].event_kind<<" for transition.");
        }

        // Apply maximum time
        if(jump_set.time_function().range().upper() > Float(timing_data.final_time)) {
            jump_set.bound_time(timing_data.final_time);
        }

        // Apply blocking conditions for other active events
        for(Set<DiscreteEvent>::const_iterator other_event_iter=blocking_events.begin();
            other_event_iter!=blocking_events.end(); ++other_event_iter)
        {
            DiscreteEvent other_event=*other_event_iter;
            if(other_event!=event) {
                const ScalarFunction& other_event_guard_function=transitions[other_event].guard_function;
                const CrossingData& other_event_crossing_data=crossing_data[other_event];
                _apply_guard(jump_set,jump_starting_state,flow,jump_step_time,
                             other_event,other_event_guard_function,other_event_crossing_data,semantics);
            }
        }

        // Apply reset
        if(!definitely(jump_set.empty())) {
            jump_set.apply_reset(event,transitions[event].target,transitions[event].reset_function);
            evolution_data.working_sets.append(jump_set);
        }

        ARIADNE_LOG(6, "jump_set="<<jump_set<<"\n");
    }

}


void
HybridEvolverBase::
_upper_evolution_flow(EvolutionData& evolution_data,
                      HybridAutomatonInterface const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_upper_evolution_flow\n");

    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    const Real final_time=maximum_hybrid_time.continuous_time();
    const uint maximum_steps=maximum_hybrid_time.discrete_time();

    // Routine check for emptiness
    if(evolution_data.working_sets.empty()) { return; }

    // Get the starting set for this round of evolution
    HybridEnclosure initial_set=evolution_data.working_sets.back(); evolution_data.working_sets.pop_back();
    ARIADNE_LOG(2,"initial_set="<<initial_set<<"\n\n");

    // Test if maximum number of steps has been reached
    if(initial_set.previous_events().size()>maximum_steps) {
        evolution_data.evolve_sets.append(initial_set); return;
    }

    if(initial_set.time_range().lower()>=final_time) {
        ARIADNE_WARN("starting_set.time_range()="<<initial_set.time_range()<<" which exceeds final time="<<final_time<<"\n");
        return;
    }

    // Extract starting location
    const DiscreteLocation location=initial_set.location();

    // Cache dynamic and constraint functions
    VectorFunction dynamic=system.dynamic_function(location);
    Map<DiscreteEvent,TransitionData> transitions = this->_extract_transitions(location,system);
    Set<DiscreteEvent> events = transitions.keys();

    ARIADNE_LOG(4,"\ndynamic="<<dynamic<<"\n");
    ARIADNE_LOG(4,"transitions="<<transitions<<"\n\n");

    // Process the initially active events; cut out active points to leave initial flowable set.
    this->_process_initial_events(evolution_data, initial_set,transitions);
    ARIADNE_ASSERT(evolution_data.starting_sets.size()<=1);

    while(!evolution_data.starting_sets.empty()) {
        this->_upper_evolution_step(evolution_data,dynamic,transitions,final_time);
    }
}

void
HybridEvolverBase::
_upper_evolution_step(EvolutionData& evolution_data,
                      VectorFunction const& dynamic,
                      Map<DiscreteEvent,TransitionData> const& transitions,
                      Real const& final_time) const
{
    ARIADNE_LOG(3,"HybridEvolverBase::_upper_evolution_step\n");
    HybridEnclosure starting_set=evolution_data.starting_sets.back(); evolution_data.starting_sets.pop_back();

    ARIADNE_LOG(2,"starting_set="<<starting_set<<"\n");
    ARIADNE_LOG(2,"starting_time="<<starting_set.time_function().polynomial()<<"\n");
    if(definitely(starting_set.empty())) {
        ARIADNE_LOG(4,"Empty starting_set "<<starting_set<<"\n");
        return;
    }

    if(starting_set.time_range().lower()>=static_cast<Float>(final_time)) {
        ARIADNE_WARN("starting_set.time_range()="<<starting_set.time_range()<<" which exceeds final time="<<final_time<<"\n");
        return;
    }

    if(verbosity==1) { _log_summary(evolution_data.working_sets.size(),evolution_data.reach_sets.size(),starting_set); }

    Map<DiscreteEvent,ScalarFunction> guard_functions;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        guard_functions.insert(transition_iter->first,transition_iter->second.guard_function);
    }
    ARIADNE_LOG(4,"guards="<<guard_functions<<"\n");

    // Compute the bounding box of the enclosure
    const Box starting_bounding_box=starting_set.space_bounding_box();

    // Compute flow and actual time step size used
    const FlowFunctionPatch flow_model=this->_compute_flow(dynamic,starting_bounding_box,this->parameters().maximum_step_size);
    ARIADNE_LOG(4,"flow_model.domain()="<<flow_model.domain()<<" flow_model.range()="<<flow_model.range()<<"\n");

    // Compute possibly active urgent events with increasing guards, and crossing times
    Set<DiscreteEvent> active_events =
        this->_compute_active_events(dynamic,guard_functions,flow_model,starting_set);
    ARIADNE_LOG(4,"active_events="<<active_events<<"\n");

    // Compute the kind of crossing (increasing, convex, etc);
    Map<DiscreteEvent,CrossingData> crossing_data =
        this->_compute_crossings(active_events,dynamic,guard_functions,flow_model,starting_set);
    ARIADNE_LOG(4,"crossing_data="<<crossing_data<<"\n");

    // Compute end conditions for flow
    TimingData timing_data = this->_compute_timing(active_events,Real(final_time),flow_model,crossing_data,starting_set);
    ARIADNE_LOG(4,"timing_data="<<timing_data<<"\n");

    // Apply the time step
    HybridEnclosure reach_set, evolve_set;
    this->_apply_time_step(evolution_data,starting_set,flow_model,timing_data,crossing_data,dynamic,transitions);
}



void
PureConstraintHybridEvolver::
_apply_time_step(EvolutionData& evolution_data,
                 HybridEnclosure const& starting_set,
                 VectorIntervalFunction const& flow,
                 TimingData const& timing,
                 Map<DiscreteEvent,CrossingData> const& crossings,
                 VectorFunction const& dynamic,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    // Apply time step without computing transition crossing times
    ARIADNE_LOG(3,"PureConstraintHybridEvolver::_apply_time_step(...)\n");
    ARIADNE_ASSERT(crossings.size()<=1);

    HybridEnclosure reach_set=starting_set;
    HybridEnclosure evolve_set=starting_set;
    reach_set.apply_flow(flow,timing.step_size);
    evolve_set.apply_flow_step(flow,timing.step_size);
    if(!crossings.empty()) {
        DiscreteEvent event=crossings.begin()->first;
        HybridEnclosure jump_set=reach_set;
        jump_set.new_guard(event,transitions[event].guard_function);
        jump_set.apply_reset(event,transitions[event].target,transitions[event].reset_function);
        if(jump_set.time_range().upper()>timing.final_time) {
            jump_set.bound_time(timing.final_time);
        }
        if(!definitely(jump_set.empty())) {
            evolution_data.intermediate_sets.append(jump_set);
            evolution_data.working_sets.append(jump_set);
            ARIADNE_LOG(4,"  jump_set="<<jump_set<<"\n");
        }
        reach_set.new_invariant(event,transitions[event].guard_function);
        evolve_set.new_invariant(event,transitions[event].guard_function);
    }
    if(reach_set.time_function().range().upper()>=timing.final_time) {
        HybridEnclosure final_set=reach_set;
        final_set.set_time(timing.final_time);
        if(!definitely(final_set.empty())) {
            evolution_data.evolve_sets.append(final_set);
            ARIADNE_LOG(4,"  final_set="<<final_set<<"\n");
        }
        reach_set.bound_time(timing.final_time);
        evolve_set.bound_time(timing.final_time);
    }
    evolution_data.reach_sets.append(reach_set);
    ARIADNE_LOG(4,"  reach_set="<<reach_set<<"\n");
    if(!definitely(evolve_set.empty())) {
        evolution_data.intermediate_sets.append(evolve_set);
        evolution_data.starting_sets.append(evolve_set);
        ARIADNE_LOG(4,"  evolve_set="<<evolve_set<<"\n");
    }
    //{ char c; cin.get(c); }
}

TimingData
DeterministicTransverseHybridEvolver::
_compute_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                VectorIntervalFunction const& flow,
                Map<DiscreteEvent,CrossingData> const& crossings,
                HybridEnclosure const& initial_set) const
{
    ARIADNE_LOG(7,"DeterministicTransverseHybridEvolver::_compute_timing(...)\n");
    // TODO: Implement this functions correctly
    TimingData result;
    result.step_size=flow.domain()[flow.domain().size()-1].upper();
    result.time_domain=Interval(0.0,result.step_size);
    result.time_coordinate=ScalarIntervalFunction::coordinate(Vector<Interval>(1u,result.time_domain),0u);
    result.final_time=final_time;
    result.remaining_time=result.final_time-initial_set.time_function();
    Interval remaining_time_range=result.remaining_time.range();
    // NOTE: The evolution time function may be negative or greater than the final time
    // over part of the parameter domain.
    if(remaining_time_range.upper()<=result.step_size) {
        result.step_kind=FINAL_STEP;
        result.finishing_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),numeric_cast<Interval>(final_time));
        result.evolution_time=numeric_cast<Interval>(final_time)-initial_set.time_function();
    } else if(remaining_time_range.lower()<=result.step_size) {
        result.step_kind=UNWIND_STEP;
        if(remaining_time_range.width()<result.step_size) {
            Float constant_finishing_time=result.final_time-remaining_time_range.lower()+result.step_size;
            result.finishing_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),constant_finishing_time);
        } else {
            // FIXME: The finishing time may need to be adjusted
            result.finishing_time=0.5*(result.step_size+initial_set.time_function());
        }
        result.evolution_time=numeric_cast<Interval>(final_time)-initial_set.time_function();
    } else if(false) { // Don't handle CREEP_STEP yet
        result.step_kind=CREEP_STEP;
        result.spacial_evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
        result.evolution_time=compose(result.spacial_evolution_time,initial_set.space_function());
    } else {
        result.step_kind=FULL_STEP;
        result.spacial_evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
        result.evolution_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),result.step_size);
    }
    return result;
}

void
DeterministicTransverseHybridEvolver::
_apply_time_step(EvolutionData& evolution_data,
                 HybridEnclosure const& starting_set,
                 VectorIntervalFunction const& flow,
                 TimingData const& timing,
                 Map<DiscreteEvent,CrossingData> const& crossings,
                 VectorFunction const& dynamic,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG(4,"\n");
    ARIADNE_LOG(3,"DeterministicTransverseHybridEvolver::_apply_time_step(...)\n");
    ARIADNE_LOG(4,timing<<"\n");
    ARIADNE_LOG(4,"starting_set="<<starting_set<<"\n");

    Set<DiscreteEvent> events=crossings.keys();
    typedef Set<DiscreteEvent>::const_iterator EventIterator;

    for(EventIterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        const DiscreteEvent event=*event_iter;
        ARIADNE_LOG(4,"event="<<event<<", crossing_time="<<polynomial(crossings[event].crossing_time)<<"\n");
        HybridEnclosure jump_set=starting_set;
        const ScalarIntervalFunction& crossing_time=crossings[event].crossing_time;
        jump_set.new_invariant(event,-crossing_time);  // Ensure crossing time is positive
        for(EventIterator other_event_iter=events.begin(); other_event_iter!=events.end(); ++other_event_iter) {
            const DiscreteEvent other_event=*other_event_iter;
            if(other_event!=event) {
                jump_set.new_invariant(other_event,transitions[other_event].guard_function);
            }
        }
        ARIADNE_LOG(4,"  active_set="<<jump_set<<"\n");
        switch(timing.step_kind) {
            case FULL_STEP:
                jump_set.new_invariant(event,(crossing_time-timing.step_size).function());
                jump_set.apply_flow_step(flow,crossing_time);
                if(jump_set.time_function().range().upper()>timing.final_time) { jump_set.bound_time(timing.final_time); }
                break;
            case CREEP_STEP:
                jump_set.new_invariant(event,crossing_time-timing.evolution_time);
                jump_set.apply_flow_step(flow,crossing_time);
                if(jump_set.time_function().range().upper()>timing.final_time) { jump_set.bound_time(timing.final_time); }
                break;
            case UNWIND_STEP:
                jump_set.apply_flow_step(flow,crossing_time);
                jump_set.bound_time(timing.finishing_time);
                break;
            case FINAL_STEP:
                jump_set.apply_flow_step(flow,crossing_time);
                jump_set.bound_time(timing.final_time);
                break;
        }
        ARIADNE_LOG(4,"  active_set="<<jump_set<<"\n");

        jump_set.apply_reset(event,transitions[event].target,transitions[event].reset_function);
        ARIADNE_LOG(4,"  jump_set="<<jump_set<<"\n");
        if(!definitely(jump_set.empty())) {
            evolution_data.working_sets.append(jump_set);
            evolution_data.intermediate_sets.append(jump_set);
        }
    }

    HybridEnclosure evolve_set=starting_set;
    HybridEnclosure reach_set=starting_set;
    switch(timing.step_kind) {
        case FULL_STEP:
            evolve_set.apply_flow_step(flow,timing.step_size);
            reach_set.apply_flow(flow,timing.step_size);
            break;
        case CREEP_STEP:
            evolve_set.apply_flow_step(flow,timing.evolution_time);
            reach_set.apply_flow(flow,timing.evolution_time);
            break;
        case UNWIND_STEP:
            ARIADNE_LOG(4,"finishing_time="<<timing.finishing_time<<"\n");
            evolve_set.apply_flow_and_set_time(flow,timing.finishing_time);
            reach_set.apply_flow_and_bound_time(flow,timing.finishing_time);
            break;
        case FINAL_STEP:
            evolve_set.apply_flow_and_set_time(flow,timing.final_time);
            reach_set.apply_flow_and_bound_time(flow,timing.final_time);
            break;
        default:
            ARIADNE_FAIL_MSG("DeterministicTransverseHybridEvolver cannot handle flow step kind "<<timing.step_kind<<"\n");
    }
    HybridEnclosure final_set;

    for(EventIterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        const DiscreteEvent event=*event_iter;
        evolve_set.new_invariant(event,transitions[event].guard_function);
        reach_set.new_invariant(event,transitions[event].guard_function);
    }

    if(timing.step_kind!=FINAL_STEP) {
        if(reach_set.time_function().range().upper()>timing.final_time) {
            HybridEnclosure final_set=reach_set;
            final_set.set_time(timing.final_time);
            if(!definitely(final_set.empty())) {
                evolution_data.evolve_sets.append(final_set);
            }
            ARIADNE_LOG(4,"  final_set="<<final_set<<"\n");
            reach_set.bound_time(timing.final_time);
            evolve_set.bound_time(timing.final_time);
        }
    }

    evolution_data.reach_sets.append(reach_set);
    ARIADNE_LOG(4,"  reach_set="<<reach_set<<"\n");

    switch(timing.step_kind) {
        case FINAL_STEP:
            // This is definitely the final step, so the evolve set is the final set
            ARIADNE_LOG(4,"  final_set="<<evolve_set<<"\n");
            if(!definitely(evolve_set.empty())) {
                evolution_data.evolve_sets.append(evolve_set);
            }
            break;
        case FULL_STEP: case CREEP_STEP: case UNWIND_STEP:
            ARIADNE_LOG(4,"  evolve_set="<<evolve_set<<"\n");
            if(!definitely(evolve_set.empty()) && evolve_set.time_function().range().lower()<timing.final_time) {
                evolution_data.starting_sets.append(evolve_set);
                evolution_data.intermediate_sets.append(evolve_set);
            } else {
                evolution_data.evolve_sets.append(evolve_set);
            }
    }
    ARIADNE_LOG(4,"\n");
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
_apply_time_step(EvolutionData& evolution_data,
                 HybridEnclosure const& starting_set,
                 VectorIntervalFunction const& flow,
                 TimingData const& timing_data,
                 Map<DiscreteEvent,CrossingData> const& crossings,
                 VectorFunction const& dynamic,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

/*
virtual void
TransverseHybridEvolver::
_apply_blocking(HybridEnclosure& set,
                const CrossingData& crossing,
                const TransitionData& transition,
                DiscreteEvent event)
{
    if(is_blocking(transition.event_kind)) {
        switch(crossing.kind) {
            case INCREASING_CROSSING: case CONVEX_CROSSING:
                set.new_invariant(transition.guard_function; break;
            case NEGATIVE_CROSSING: case DECREASING_CROSSING:
                break;
            case CONCAVE_CROSSING: case DEGENERATE_CROSSING:
                ARIADNE_FAIL_MESSAGE("TransverseHybridEvolver cannot handle "<<crossing.kind<<" crossing.";
            default:
                ARIADNE_FAIL_MSG("CrossingKind "<<crossing.kind<<" not recognised by TransverseHybridEvolver.");
        }
    }
}
*/

void
TransverseHybridEvolver::
_apply_time_step(EvolutionData& evolution_data,
                 HybridEnclosure const& starting_set,
                 VectorIntervalFunction const& flow,
                 TimingData const& timing_data,
                 Map<DiscreteEvent,CrossingData> const& crossings,
                 VectorFunction const& dynamic,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    HybridEnclosure reach_set=starting_set;
    for(Map<DiscreteEvent,TransitionData>::const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event = transition_iter->first;
        CrossingData const& crossing = crossings[event];
        TransitionData const& transition = transition_iter->second;
        CrossingKind crossing_kind=crossing.crossing_kind;
        if(is_blocking(transition.event_kind)) {
            switch(crossing.crossing_kind) {
                case INCREASING_CROSSING: case CONVEX_CROSSING:
                    reach_set.new_invariant(event,transition.guard_function); break;
                case NEGATIVE_CROSSING: case DECREASING_CROSSING:
                    break;
                case CONCAVE_CROSSING: case DEGENERATE_CROSSING:
                    ARIADNE_FAIL_MSG("TransverseHybridEvolver cannot handle "<<crossing_kind<<" crossing.");
                default:
                    ARIADNE_FAIL_MSG("CrossingKind "<<crossing_kind<<" not recognised by TransverseHybridEvolver.");
            }
        }
    }


}



} // namespace Ariadne
