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

bool ALLOW_FINAL = true;
bool ALLOW_UNWIND = false;
bool ALLOW_CREEP = true;

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



static const DiscreteEvent final_event("_tmax_");
static const DiscreteEvent step_event("_h_");


Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridAutomatonInterface& system,
      const HybridEnclosure& initial,
      const HybridTime& time,
      Semantics semantics) const
{
    ARIADNE_LOG(2,"\nHybridEvolverBase::orbit(...): verbosity="<<verbosity<<"\n");

    EvolutionData evolution_data;
    evolution_data.semantics=semantics;
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
            <<" t=["<<std::setw(4)<<std::setprecision(3)<<std::left<<std::fixed<<starting_time_range.lower()
            <<","<<std::setw(4)<<std::left<<std::fixed<<starting_time_range.upper()<<"]"<<std::flush
            <<" te="<<std::setw(7)<<std::scientific<<std::setprecision(1)<<starting_set.time_function().error()<<std::flush
            <<" #p="<<std::setw(2)<<std::left<<starting_set.number_of_parameters()
            <<" #c="<<std::setw(1)<<std::left<<starting_set.number_of_inequality_constraints()<<"+"<<std::setw(2)<<starting_set.number_of_equality_constraints()
            <<" r="<<std::setw(4)<<std::fixed<<std::setprecision(3)<<starting_bounding_box.radius()
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
    Set<DiscreteEvent> events = system.events(location);
    for(Set<DiscreteEvent>::const_iterator event_iter=events.begin();
        event_iter!=events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        EventKind event_kind=system.event_kind(location,event);
        ScalarFunction guard_function=system.guard_function(location,event);
        VectorFunction reset_function; DiscreteLocation target;
        if(is_activating(event_kind)) {
            reset_function=system.reset_function(location,event);
            target=system.target(location,event);
        }
        TransitionData transition_data={event,event_kind,guard_function,target,reset_function};
        transitions.insert(event,transition_data);
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
        // Assume that an impact cannot occur immediately after any other event.
        // This is true after the impact, since $L_{f}g$ changes sign.
        // After a different event, this is not so clear, though it should be
        // a modelling error for the set to intersect the interior of the guard
        // TODO: The condition for an impact to occur is $L_{f}g>0$.
        //       However, $f$ is not available to this method; maybe it should be.
        // FIXME: Need to consider impact which really can occur immediately due to mapping to boundary.
        // FIXME: Maybe need to consider impact as invariant.
        if(transition.event_kind!=INVARIANT && transition.event_kind!=IMPACT) {
            if(possibly(initial_set.satisfies(transition.guard_function>=0))) {
                if(transition.event_kind!=PROGRESS) {
                    HybridEnclosure immediate_jump_set=invariant_set;
                    immediate_jump_set.new_activation(event,transition.guard_function);
                    if(!definitely(immediate_jump_set.empty())) {
                        // Put the immediate jump set in the reached sets, since it does not occur in the flowable set
                        ARIADNE_LOG(1,"  "<<event<<": "<<transition.event_kind<<", immediate\n");
                        evolution_data.reach_sets.append(immediate_jump_set);
                        immediate_jump_set.apply_reset(event,transition.target,transition.reset_function);
                        ARIADNE_LOG(4,"immediate_jump_set="<<immediate_jump_set<<"\n");
                        evolution_data.intermediate_sets.append(immediate_jump_set);
                        evolution_data.working_sets.append(immediate_jump_set);
                    }
                }
                if(transition_iter->second.event_kind!=PERMISSIVE) {
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
    TaylorIntegrator integrator(32,this->parameters().flow_accuracy);
    // Compute flow and actual time step size used
    //
    // The Integrator classes compute the flow as a function on a symmetrical
    // time domain [-h,+h], since this means the time is centred at 0.
    // We then restrict to the time domain [0,h] since this can make evaluation
    // more accurate, and the time domain might be used explicitly for the domain
    // of the resulting set.
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
    reach_set.apply_flow_for(flow,flow.domain()[flow.domain().size()-1].upper());
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
                // FIXME: Need to allow permissive events with strictly decreasing guard.
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
    Map<DiscreteEvent,CrossingData> crossings;
    crossings.clear();

    Box flow_bounds=flow.range();
    for(Set<DiscreteEvent>::const_iterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        const DiscreteEvent event=*event_iter;
        ScalarFunction const& guard=guards[event];

        // Compute the derivative of the guard function g along flow lines of $\dot(x)=f(x)$
        // This is given by the Lie derivative at a point x, defined as $L_{f}g(x) = (\nabla g\cdot f)(x)$
        ScalarFunction derivative=lie_derivative(guard,dynamic);
        Interval derivative_range=derivative.evaluate(flow_bounds);
        if(derivative_range.lower()>0.0) {
            // If the derivative $L_{f}g$is strictly positive over the bounding box for the flow,
            // then the guard function is strictly increasing.
            // There is at most one crossing with the guard, and the time of this
            // crossing must be the time of the event along the trajectory.
            // The crossing time $\gamma(x_0)$ given the initial state can usually be computed
            // by solving the equation $g(\phi(x_0,\gamma(x_0))) = 0$
            ScalarIntervalFunction crossing_time;
            try {
                crossing_time=implicit(compose(guard,flow));
                crossings[event]=CrossingData(INCREASING_CROSSING,crossing_time);
            }
            catch(const ImplicitFunctionException& e) {
                // If the crossing time cannot be computed, then we can still
                // use the fact that the crossing occurs as soon as $g(x(t))=0$.
                // Since this is the same condition as for convex guard function
                // along flow lines (L_{f}^{2} g > 0$, we set the crossing type
                // as CONVEX_CROSSING
                ARIADNE_LOG(0,"Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                crossings[event]=CrossingData(CONVEX_CROSSING);
            }
        } else if(derivative_range.upper()<0.0) {
            // If the derivative is strictly negative over the bounding box for the flow,
            // then the guard function is strictly decreasing.
            // This means that the event is either initially active, or does not occur.
            // There is no need to compute a crossing time.
            crossings[event]=CrossingData(DECREASING_CROSSING);
        } else {
            // If the derivative of the guard function along flow lines cannot be shown
            // to have a definite sign over the entire flow box, then try to compute
            // the sign of the second derivative $L_{f}^{2}g(x)=L_{f}L_{f}g(x)$.
            ScalarFunction second_derivative=lie_derivative(derivative,dynamic);
            Interval second_derivative_range=second_derivative.evaluate(flow_bounds);
            if(second_derivative_range.lower()>0.0) {
                // If the second derivative is positive, then either
                //    (i) the event is immediately active
                //   (ii) the event is never active, or
                //  (iii) the event is initially inactive, but becomes active
                //        due to a transverse crossing.
                //   (iv) the initial state is on the boundary of the guard
                //        set, possibly with the flow tangent to this set
                // We cannot compute the crossing time, even in case (iii),
                // due to the singularity due to the tangency in (iv). However,
                // we do know that in (iii), the event occurs when $t>0$ and
                // $g(\phi(x_0,t))=0$. The crossing time is not computed.
                crossings[event]=CrossingData(CONVEX_CROSSING);
            } else if(second_derivative_range.upper()<0.0) {
                // If the second derivative is negative, then the guard
                // values $g(x(t))$ are concave along flow lines. There are
                // four main cases:
                //   (i) The event is initially active.
                //  (ii) The event is not initially active, but later becomes active.
                // (iii) The event is never active, but would become active if
                //       flowing backwards in time.
                //  (iv) The event is never active, and the maximum value of
                //       the guard along the flow lines is zero.
                // Additionally, there is the degenerate case
                //   (v) At some point in the (forward) flow, the state touches
                //       the guard set at a point of tangency.
                // Due to the presence of the tangency, the event time is
                // not a smooth function of the initial state. Further, since
                // some trajectories cross the boundary of the guard set twice,
                // the condition $g(x(t))=0$ is not sufficient for determining
                // the jump time. A necessary and sufficient condition,
                // assuming the event is not initially active, is that
                // $g(x)=0$ and $L_{f}g(x)\geq 0$. Alternatively, we can compute
                // the <em>critical time</em> $|mu(x_0) at which the guard value
                // reaches a maximum. A necessary and sufficient condition
                // is then $g(\phi(x_0,t))=0$ and $t<= \mu(x_0)$.
                //   Note that while $g(\phi(x_0,t))=0$ and $L_{f}f(\phi(x_0,t))\geq0$ is a
                // necessary and sufficient condition for a crossing, there
                // is no necessary and sufficient condition for no crossing
                // which does not involve the critical time. A necessary and
                // sufficient condition for no crossing involving the critical
                // time is $(g(\phi(x_0,t))<=0 /\ t<=\mu(x_0)) \/ g(\phi(x_0,\mu(x_0)))<=0$
                try {
                    ScalarIntervalFunction critical_time=implicit(compose(derivative,flow));
                    crossings[event]=CrossingData(CONCAVE_CROSSING);
                    crossings[event].critical_time=critical_time;
                }
                catch(const ImplicitFunctionException& e) {
                    ARIADNE_LOG(0,"Error in computing crossing time for event "<<*event_iter<<":\n  "<<e.what()<<"\n");
                    crossings[event]=CrossingData(DEGENERATE_CROSSING);
                }
            } else {
                // The crossing cannot be shown to be one of the kinds mentioned
                // above. A theoretically exact expression for the crossing
                // set is generally not available.
                crossings[event]=CrossingData(DEGENERATE_CROSSING);
            }
        }
    }
    return crossings;
}




TimingData
HybridEvolverBase::
_compute_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                VectorIntervalFunction const& flow,
                Map<DiscreteEvent,CrossingData>& crossings,
                HybridEnclosure const& initial_set) const
{
    // Compute the evolution time for the given step.
    ARIADNE_LOG(7,"HybridEvolverBase::_compute_timing(...)\n");
    const Float step_size=flow.domain()[flow.domain().size()-1].upper();
    TimingData result;
    result.step_size=step_size;
    result.final_time=final_time;
    result.time_domain=Interval(0.0,step_size);
    result.time_coordinate=ScalarIntervalFunction::identity(result.time_domain);

    // NOTE: The starting time function may be negative or greater than the final time
    // over part of the parameter domain.
    ScalarIntervalFunction const& starting_time=initial_set.time_function();
    Interval starting_time_range=initial_set.time_function().range();

    ScalarIntervalFunction remaining_time=Interval(final_time)-starting_time;
    Interval remaining_time_range=remaining_time.range();
    if(remaining_time_range.upper()<=result.step_size) {
        // The rest of the evolution can be computed within a single time step.
        // The step kind is given as FINAL_STEP so that the evolution algorithm
        // knows that the evolved set does not need to be evolved further.
        // This knowledge is required to be given combinarially, since
        // specifying the final time as a constant TaylorFunction is not
        // exact if the final_time parameter is not exactly representable as
        // a Float
        result.step_kind=FINAL_STEP;
        result.finishing_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),Interval(final_time));
        result.evolution_time=Interval(final_time)-starting_time;
    } else if(starting_time_range.upper()>final_time) {
        // Some of the points may already have reached the final time.
        // Don't try anythin fancy, just to a simple full step.
        result.step_kind=FULL_STEP;
        result.spacial_evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
        result.evolution_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),result.step_size);
    } else if(remaining_time_range.lower()<=step_size) {
        // Some of the evolved points can be evolved to the final time in a single step
        // The evolution is performed over a step size which preserves points
        // which are past the final_time, but moves other points closer.

        // Using the final_time as a guide, set the finishing time to closer to the final time.
        // This method ensures that points do not pass the final time after the transition.
        result.step_kind=UNWIND_STEP;
        Float sf=1.0;
        while(remaining_time_range.upper()*sf>step_size) { sf /= 2; }
        result.finishing_time=(1-sf)*(starting_time-Interval(final_time))+Interval(final_time);
        result.evolution_time=result.finishing_time-starting_time;
    } else if(ALLOW_UNWIND && crossings.empty() && starting_time_range.lower()<starting_time_range.upper()) {
        // Try to unwind the evolution time to a constant
        result.step_kind=UNWIND_STEP;
        if(starting_time_range.width()<step_size/2) {
            result.finishing_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),starting_time_range.lower()+step_size);
        } else {
            // Try to reduce the time interval by half the step size
            // Corresponds to setting omega(smin)=tau(smin)+h, omega(smax)=tau(smax)+h/2
            // Taking omega(s)=a tau(s) + b, we obtain
            //   a=1-h/2(tmax-tmin);  b=h(tmax-tmin/2)/(tmax-tmin) = (2tmax-tmin)a
            Float h=result.step_size; Float tmin=starting_time_range.lower(); Float tmax=starting_time_range.upper();
            Float a=1-((h/2)/(tmax-tmin)); Float b=h*(tmax-tmin/2)/(tmax-tmin);
            result.finishing_time=a*starting_time+b;
            //std::cerr<<"\nparameter_domain="<<initial_set.parameter_domain()<<"\nstarting_time="<<starting_time<<"\nfinishing_time="<<result.finishing_time<<"\n"; char c; std::cin>>c;
        }
        result.evolution_time=result.finishing_time-starting_time;
    } else if(ALLOW_CREEP && !crossings.empty()) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        ARIADNE_LOG(6,"\nPossible creep step; h="<<step_size<<"\n");
        bool creep=false;
        for(Map<DiscreteEvent,CrossingData>::const_iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); ++crossing_iter)
        {
            ARIADNE_LOG(6,"  "<<crossing_iter->first<<": "<<crossing_iter->second.crossing_kind<<", "<<crossing_iter->second.crossing_time.range()<<"\n");
            if(crossing_iter->second.crossing_kind==INCREASING_CROSSING) {
                const ScalarIntervalFunction& crossing_time=crossing_iter->second.crossing_time;
                Interval crossing_time_range=crossing_time.range();
                if(crossing_time_range.upper()<=step_size) {
                    // This event ensures that the evolve set is empty after a full step, so use this.
                    creep=false;
                    ARIADNE_LOG(6,"No creep; event "<<crossing_iter->first<<"completely taken\n");
                    break;
                } else if(crossing_time_range.lower()<=0.0) {
                    // This event is already partially active, so carry on with a full step
                    ARIADNE_LOG(6,"No creep; event "<<crossing_iter->first<<"alrealy partially active\n");
                    creep=false;
                    break;
                } else if(crossing_time_range.lower()>=step_size) {
                    //crossings.erase(crossing_iter++);
                } else {
                    creep=true;
                }
            }
        }
        if(creep==true) {
            // Compute reduced evolution time; note that for every remaining increasing crossing, we have
            // crossing_time_range.lower()>0.0 and crossing_time_range.upper()>step_size.
            ScalarIntervalFunction evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),step_size);
            for(Map<DiscreteEvent,CrossingData>::const_iterator crossing_iter=crossings.begin();
                crossing_iter!=crossings.end(); ++crossing_iter)
            {
                if(crossing_iter->second.crossing_kind==INCREASING_CROSSING) {
                    // Modify the crossing time function to be the smallest possible; this ensures that the evaluation time is
                    // essentially exact
                    ScalarIntervalFunction lower_crossing_time=crossing_iter->second.crossing_time;
                    Float crossing_time_error=lower_crossing_time.error();
                    lower_crossing_time.set_error(0.0);
                    lower_crossing_time-=crossing_time_error;

                    // One possibility is to use quadratic restrictions
                    //   If 0<=x<=2h, then x(1-x/4h)<=min(x,h)
                    //   If 0<=x<=4h, then x(1-x/8h)<=min(x,2h)
                    // Formula below works if x<=2h
                    //   evolution_time=evolution_time*(crossing_time/result.step_size)*(1.0-crossing_time/(4*result.step_size));

                    // Prefer simpler linear restrictions.
                    // Multiply evolution time by crossing_time/max_crossing_time
                    evolution_time=evolution_time*(lower_crossing_time/lower_crossing_time.range().upper());
                }
            }
            // Erase increasing crossings since these cannot occur
            for(Map<DiscreteEvent,CrossingData>::iterator crossing_iter=crossings.begin();
                crossing_iter!=crossings.end(); )
            {
                if(crossing_iter->second.crossing_kind==INCREASING_CROSSING) {
                    crossings.erase(crossing_iter++);
                } else {
                    ++crossing_iter;
                }
            }
            ARIADNE_LOG(6,"Creep step: evolution_time="<<evolution_time<<"\n");
            result.step_kind=CREEP_STEP;
            result.spacial_evolution_time=evolution_time;
            result.evolution_time=compose(result.spacial_evolution_time,initial_set.space_function());
            // TODO: Remove all definitely non-active events
        } else {
            // Perform the evolution over a full time step
            result.step_kind=FULL_STEP;
            result.spacial_evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
            result.evolution_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),result.step_size);
        }
    } else {
        result.step_kind=FULL_STEP;
        result.spacial_evolution_time=ScalarIntervalFunction::constant(initial_set.space_bounding_box(),result.step_size);
        result.evolution_time=ScalarIntervalFunction::constant(initial_set.parameter_domain(),result.step_size);
    }
    return result;
}





// Apply guard to a single set.
// In the case of concave crossings, splits the set into two, one part
// corresponding to points which actually hit the set (and stop on first crossing)
// the other part corresponding to points which miss the set.
void HybridEvolverBase::
_apply_guard(List<HybridEnclosure>& sets,
             const HybridEnclosure& starting_set,
             const VectorIntervalFunction& flow,
             const ScalarIntervalFunction& elapsed_time,
             const TransitionData& transition_data,
             const CrossingData crossing_data,
             const Semantics semantics) const
{
    static const uint SUBDIVISIONS_FOR_DEGENERATE_CROSSING = 2;
    const DiscreteEvent event=transition_data.event;
    const ScalarFunction& guard_function=transition_data.guard_function;
    VectorIntervalFunction starting_state=starting_set.space_function();
    if(elapsed_time.argument_size()>starting_state.argument_size()) {
        starting_state=embed(starting_state,elapsed_time.domain()[elapsed_time.domain().size()-1]);
    }
    ARIADNE_ASSERT(starting_state.domain()==elapsed_time.domain());

    List<HybridEnclosure>::iterator end=sets.end();
    for(List<HybridEnclosure>::iterator iter=sets.begin(); iter!=end; ++iter) {
        HybridEnclosure& set=*iter;

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
                // Split the set into two components, one corresponding to
                // points which miss the guard completely, the other to points which
                // eventually hit the guard, ensuring that the set stops at the first crossing
                HybridEnclosure extra_set=set;
                set.new_parameter_constraint(event,maximal_guard<=0.0);
                extra_set.new_parameter_constraint( event, elapsed_time <= critical_time );
                extra_set.new_state_constraint(event,guard_function<=0.0);
                sets.append(extra_set);
                break;
                // Code below is always exact, but uses two sets
                // set1.new_parameter_constraint(event,final_guard <= 0);
                // set2.new_parameter_constraint(event, elapsed_time <= critical_time);
                // set1.new_parameter_constraint(event, maximal_guard <= 0.0);
                // set2.new_parameter_constraint(event, elapsed_time >= critical_time);
            }
            case DEGENERATE_CROSSING: {
                // The crossing with the guard set is not one of the kinds handled above.
                // We obtain an over-appproximation by testing at finitely many time points
                const uint n=SUBDIVISIONS_FOR_DEGENERATE_CROSSING;
                switch(semantics) {
                    case UPPER_SEMANTICS:
                        for(uint i=0; i!=n; ++i) {
                            Float alpha=Float(i+1)/n;
                            ScalarIntervalFunction intermediate_guard
                                = compose( guard_function, compose( flow, join(starting_state, alpha*elapsed_time) ) );
                            set.new_parameter_constraint(event, intermediate_guard <= 0);
                        }
                        break;
                    case LOWER_SEMANTICS:
                        // Can't continue the evolution, so set a trivially-falsified constraint
                        set.new_parameter_constraint(event, ScalarIntervalFunction::constant(set.parameter_domain(),1.0) <= 0.0);
                        break;
                }
                break;
            }
            case POSITIVE_CROSSING:
                // No need to do anything since all points are initially
                // active and should have been handled already
            case NEGATIVE_CROSSING:
                // No points are active
            case DECREASING_CROSSING:
                // No need to do anything, since only initially active points
                // become active during the evolution, and these have been
                // handled already.
                break;
        }
    }
}

DiscreteEvent step_time_event("step_time");

void
HybridEvolverBase::
_apply_reach_step(HybridEnclosure& set,
                  VectorIntervalFunction const& flow,
                  TimingData const& timing_data) const
{
    switch(timing_data.step_kind) {
        case FULL_STEP:
            set.apply_flow_for(flow,timing_data.step_size);
            break;
        case CREEP_STEP:
            set.apply_flow_for(flow,timing_data.spacial_evolution_time);
            break;
        case PARTIAL_STEP:
            set.apply_reach_step(flow,timing_data.evolution_time);
            ARIADNE_ASSERT(false); break; // Not currently implemented
        case UNWIND_STEP:
            set.apply_flow_to(flow,timing_data.finishing_time);
            break;
        case FINAL_STEP:
            set.apply_flow_to(flow,timing_data.final_time);
            break;
    }
}

void
HybridEvolverBase::
_apply_evolve_step(HybridEnclosure& set,
                  VectorIntervalFunction const& flow,
                  TimingData const& timing_data) const
{
    switch(timing_data.step_kind) {
        case FULL_STEP:
            set.apply_evolve_step(flow,timing_data.evolution_time);
            break;
        case CREEP_STEP:
            set.apply_evolve_step(flow,timing_data.evolution_time);
            break;
        case PARTIAL_STEP:
            set.apply_evolve_step(flow,timing_data.evolution_time);
            break;
        case UNWIND_STEP:
            set.apply_evolve_step(flow,timing_data.evolution_time);
            set.set_time_function(timing_data.finishing_time);
            break;
        case FINAL_STEP:
            set.apply_evolve_step(flow,timing_data.evolution_time);
            set.set_time_function(timing_data.finishing_time);
            break;
    }
}

void
HybridEvolverBase::
_apply_guard_step(HybridEnclosure& set,
                  VectorFunction const& dynamic,
                  VectorIntervalFunction const& flow,
                  TimingData const& timing_data,
                  TransitionData const& transition_data,
                  CrossingData const& crossing_data) const
{
    // Compute flow to guard set up to evolution time.
    HybridEnclosure& jump_set=set;
    const DiscreteEvent event=transition_data.event;
    VectorIntervalFunction starting_state=set.space_function();
    VectorIntervalFunction reach_starting_state=embed(starting_state,timing_data.time_domain);
    ScalarIntervalFunction reach_step_time=embed(starting_state.domain(),timing_data.time_coordinate);
    ScalarIntervalFunction step_time;

    switch(transition_data.event_kind) {
        case PERMISSIVE:
            // The continuous evolution is just the same as a reachability step,
            // so we need to embed the starting state and the step time function into one higher dimension.
            jump_set.apply_reach_step(flow,timing_data.evolution_time);
            jump_set.new_activation(event,transition_data.guard_function);
            break;
        case URGENT: case IMPACT:
            switch(crossing_data.crossing_kind) {
                case INCREASING_CROSSING:
                    step_time=compose(crossing_data.crossing_time,starting_state);
                    // If the jump step might occur after the final evolution time, then introduce constraint that this does not happen
                    if(timing_data.step_kind!=FULL_STEP && (step_time-timing_data.evolution_time).range().upper()>0.0) {
                        jump_set.new_parameter_constraint(step_event,step_time<=timing_data.evolution_time);
                    }
                    jump_set.apply_flow_step_for(flow,crossing_data.crossing_time);
                    break;
                case CONVEX_CROSSING:
                    jump_set.apply_flow_to(flow,timing_data.evolution_time);
                    jump_set.new_guard(event,transition_data.guard_function);
                    break;
                case CONCAVE_CROSSING:
                case DEGENERATE_CROSSING: // Just check positive derivative in this case; NOT EXACT
                    jump_set.apply_flow_to(flow,timing_data.evolution_time);
                    jump_set.new_guard(event,transition_data.guard_function);
                    jump_set.new_invariant(event,lie_derivative(transition_data.guard_function,dynamic));
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
            ARIADNE_FAIL_MSG("Invalid event kind "<<transition_data.event_kind<<" for transition.");
    }

}

void GeneralHybridEvolver::
_apply_time_step(EvolutionData& evolution_data,
                 HybridEnclosure const& starting_set,
                 VectorIntervalFunction const& flow,
                 TimingData const& timing_data,
                 Map<DiscreteEvent,CrossingData> const& crossings,
                 VectorFunction const& dynamic,
                 Map<DiscreteEvent,TransitionData> const& transitions) const
{
    // FIXME: Make semantics function argument

    ARIADNE_LOG(2,"GeneralHybridEvolver::_apply_time_step(...)\n");
    //ARIADNE_LOG(4,"evolution_time="<<evolution_time.range()<<" final_time="<<final_time<<"\n");

    ARIADNE_ASSERT(!definitely(starting_set.empty()));

    ScalarIntervalFunction critical_function;

    // Compute events enabling transitions and events blocking continuous evolution
    Set<DiscreteEvent> active_events=crossings.keys();
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

    // Cache starting state and time functions
    VectorIntervalFunction const& starting_state=starting_set.space_function();
    ScalarIntervalFunction const& starting_time=starting_set.time_function();
    VectorIntervalFunction embedded_starting_state=embed(starting_set.space_function(),timing_data.time_domain);
    ScalarIntervalFunction embedded_reach_time=embed(starting_set.parameter_domain(),timing_data.time_coordinate);
    ScalarIntervalFunction embedded_evolution_time=embed(timing_data.evolution_time,timing_data.time_domain);
    ScalarIntervalFunction reach_step_time=embedded_reach_time;
    Semantics semantics = evolution_data.semantics;
    ARIADNE_LOG(8,"embedded_reach_time="<<embedded_reach_time<<"\n")
    ARIADNE_LOG(8,"embedded_evolution_time="<<embedded_evolution_time<<"\n")


    {
        // Compute the reach set, allowing for possibility of splitting due to non-transverse guards
        List<HybridEnclosure> reach_sets(starting_set);
        HybridEnclosure& reach_set=reach_sets.front();
        _apply_reach_step(reach_set,flow,timing_data);
        ARIADNE_LOG(6,"unconstrained_reach_set="<<reach_sets.front()<<"\n");

        // Apply constraints on reach set due to final time.
        if(timing_data.step_kind!=FINAL_STEP) {
            if(reach_set.time_range().upper()>timing_data.final_time) {
                //reach_set.new_parameter_constraint(step_event,embedded_reach_time<=embedded_evolution_time);
                reach_set.new_parameter_constraint(final_event,reach_set.time_function()<=timing_data.final_time);
            }
        }

        // Apply constraints on reach set due to events
        for(Set<DiscreteEvent>::const_iterator event_iter=blocking_events.begin();
            event_iter!=blocking_events.end(); ++event_iter)
        {
            const TransitionData& transition_data=transitions[*event_iter];
            const CrossingData& crossing_data=crossings[*event_iter];
            _apply_guard(reach_sets,starting_set,flow,embedded_reach_time,
                         transition_data,crossing_data,semantics);
        }

        ARIADNE_LOG(4,"reach_sets="<<reach_sets<<"\n");
        evolution_data.reach_sets.append(reach_sets);

    } // Done computing reach sets

    ARIADNE_LOG(4,"Done computing reach sets\n");

    // Compute evolve set
    // TODO: This should probably be a separate function
    List<HybridEnclosure> evolve_sets((starting_set));
    HybridEnclosure& evolve_set=evolve_sets.front();

    // Apply the flow to the evolve set
    if(timing_data.step_kind!=FULL_STEP) { ARIADNE_LOG(1,"  step_kind="<<timing_data.step_kind<<"\n"); }

    uint nc=evolve_set.number_of_constraints();
    _apply_evolve_step(evolve_sets.front(),flow,timing_data);
    ARIADNE_ASSERT(nc==evolve_set.number_of_constraints());

    if(evolve_set.time_function().range().upper()>timing_data.final_time) {
        evolve_set.bound_time(timing_data.final_time);
    }

    // Apply constraints on evolve set due to events
    // FIXME: If there is a discrete transition, the guard may be
    for(Set<DiscreteEvent>::const_iterator event_iter=blocking_events.begin();
        event_iter!=blocking_events.end(); ++event_iter)
    {
        const TransitionData& transition_data=transitions[*event_iter];
        const CrossingData& crossing_data=crossings[*event_iter];
        _apply_guard(evolve_sets,starting_set,flow,timing_data.evolution_time,
                     transition_data,crossing_data,semantics);
    }

    // Compute final set, and apply constraints to evolve set
    // if necessary
    if(timing_data.step_kind!=FINAL_STEP) {
        for(List<HybridEnclosure>::const_iterator evolve_set_iter=evolve_sets.begin(); evolve_set_iter!=evolve_sets.end(); ++evolve_set_iter) {
            HybridEnclosure const& evolve_set=*evolve_set_iter;
            if(!definitely(evolve_set.empty())) {
                ARIADNE_LOG(4,"evolve_set="<<evolve_set<<"\n");
                evolution_data.intermediate_sets.append(evolve_set);
                evolution_data.working_sets.append(evolve_set);
            }
        }
    } else {
        for(List<HybridEnclosure>::const_iterator final_set_iter=evolve_sets.begin(); final_set_iter!=evolve_sets.end(); ++final_set_iter) {
            HybridEnclosure const& final_set=*final_set_iter;
            if(!definitely(final_set.empty())) {
                ARIADNE_LOG(4,"final_set="<<final_set<<"\n");
                evolution_data.evolve_sets.append(final_set);
            }
        }
    }


    // Compute jump sets
    for(Set<DiscreteEvent>::const_iterator event_iter=transition_events.begin();
        event_iter!=transition_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;

        // Compute active set
        List<HybridEnclosure> jump_sets((starting_set));
        HybridEnclosure& jump_set=jump_sets.front();
        ARIADNE_LOG(2,"  "<<event<<": "<<transitions[event].event_kind<<", "<<crossings[event].crossing_kind<<"\n");
        _apply_guard_step(jump_set,dynamic,flow,timing_data,transitions[event],crossings[event]);

        ScalarIntervalFunction jump_step_time;
        switch(crossings[event].crossing_kind) {
            case CONVEX_CROSSING: case CONCAVE_CROSSING: case DEGENERATE_CROSSING:
                jump_step_time=reach_step_time;
                break;
            case INCREASING_CROSSING:
                jump_step_time=compose(crossings[event].crossing_time,starting_set.space_function());
                break;
            case DECREASING_CROSSING: case POSITIVE_CROSSING: case NEGATIVE_CROSSING:
                // No need to set jump set time;
                break;
            default:
                ARIADNE_FAIL_MSG("Unknown crossing kind "<<crossings[event].crossing_kind);
        }

        // Apply maximum time bound, as this will be applied after the next flow step
        for(List<HybridEnclosure>::iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& jump_set=*jump_set_iter;
            if(jump_set.time_function().range().upper() > Float(timing_data.final_time)) {
                jump_set.bound_time(timing_data.final_time);
            }
        }

        // Apply blocking conditions for other active events
        for(Set<DiscreteEvent>::const_iterator other_event_iter=blocking_events.begin();
            other_event_iter!=blocking_events.end(); ++other_event_iter)
        {
            DiscreteEvent other_event=*other_event_iter;
            if(other_event!=event) {
                const TransitionData& other_transition_data=transitions[other_event];
                const CrossingData& other_crossing_data=crossings[other_event];
                _apply_guard(jump_sets,starting_set,flow,jump_step_time,
                             other_transition_data,other_crossing_data,semantics);
            }
        }

        ARIADNE_LOG(1, "  "<<event<<": "<<transitions[event].event_kind<<", "<<crossings[event].crossing_kind<<"\n");
        // Apply reset
        for(List<HybridEnclosure>::iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& jump_set=*jump_set_iter;
            if(!definitely(jump_set.empty())) {
                jump_set.apply_reset(event,transitions[event].target,transitions[event].reset_function);
                evolution_data.working_sets.append(jump_set);
                ARIADNE_LOG(6, "jump_set="<<jump_set<<"\n");
            }
        }
    }

}


void
HybridEvolverBase::
_upper_evolution_flow(EvolutionData& evolution_data,
                      HybridAutomatonInterface const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    //  Select a working set and evolve this in the current location until either
    // all initial points have undergone a discrete transition (possibly to
    // the same location) or the final time is reached.
    //   Evolving within one location avoids having to re-extract event sets,
    // and means that initially active events are tested for only once.
    ARIADNE_LOG(3,"HybridEvolverBase::_upper_evolution_flow\n");

    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    const Real final_time=maximum_hybrid_time.continuous_time();
    const uint maximum_steps=maximum_hybrid_time.discrete_time();

    // Routine check for emptiness
    if(evolution_data.working_sets.empty()) { return; }

    // Get the starting set for this round of evolution
    HybridEnclosure starting_set=evolution_data.working_sets.back(); evolution_data.working_sets.pop_back();
    ARIADNE_LOG(2,"starting_set="<<starting_set<<"\n\n");

    // Test if maximum number of steps has been exceeded; if so, the set should be discarded.
    // NOTE: We could also place a test for the maximum number of steps being reaches which computing jump sets
    // This is not done since the maximum_steps information is not passed to the _apply_time_step(...) method.
    if(starting_set.previous_events().size()>maximum_steps) {
        ARIADNE_LOG(4,"starting_set "<<starting_set<<" has undergone more than maximum number of events "<<maximum_steps<<"\n");
        return;
    }

    // NOTE: Uncomment the lines below to stop evolution immediately after the maximum event, without further flow
    //if(starting_set.previous_events().size()>maximum_steps) {
    //    evolution_data.evolve_sets.append(starting_set);
    //    return;
    //}

    if(starting_set.time_range().lower()>=final_time) {
        ARIADNE_WARN("starting_set.time_range()="<<starting_set.time_range()<<" which exceeds final time="<<final_time<<"\n");
        return;
    }

    // Extract starting location
    const DiscreteLocation location=starting_set.location();

    // Cache dynamic and constraint functions
    VectorFunction dynamic=system.dynamic_function(location);
    Map<DiscreteEvent,TransitionData> transitions = this->_extract_transitions(location,system);
    Set<DiscreteEvent> events = transitions.keys();

    ARIADNE_LOG(4,"\ndynamic="<<dynamic<<"\n");
    ARIADNE_LOG(4,"transitions="<<transitions<<"\n\n");

    // Process the initially active events; cut out active points to leave initial flowable set.
    this->_process_initial_events(evolution_data, starting_set,transitions);
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
    Map<DiscreteEvent,CrossingData> crossings =
        this->_compute_crossings(active_events,dynamic,guard_functions,flow_model,starting_set);
    ARIADNE_LOG(4,"crossings="<<crossings<<"\n");

    // Compute end conditions for flow
    TimingData timing_data = this->_compute_timing(active_events,Real(final_time),flow_model,crossings,starting_set);
    ARIADNE_LOG(4,"timing_data="<<timing_data<<"\n");

    // Apply the time step
    HybridEnclosure reach_set, evolve_set;
    this->_apply_time_step(evolution_data,starting_set,flow_model,timing_data,crossings,dynamic,transitions);
}





} // namespace Ariadne