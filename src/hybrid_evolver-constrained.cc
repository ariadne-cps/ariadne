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

namespace Ariadne {

DiscreteEvent final_time_event("_final_time_");
DiscreteEvent time_step_event("_step_size_");



template<class T> std::string str(const T& t) { std::stringstream ss; ss<<t; return ss.str(); }

typedef Vector<Interval> IntervalVector;




Orbit<HybridEnclosure>
ConstrainedImageSetHybridEvolver::
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


void
ConstrainedImageSetHybridEvolver::
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

void
ConstrainedImageSetHybridEvolver::
_upper_evolution_step(List<HybridEnclosure>& working_sets,
                      ListSet<HybridEnclosure>& evolve_sets,
                      ListSet<HybridEnclosure>& reach_sets,
                      ListSet<HybridEnclosure>& intermediate_sets,
                      HybridAutomatonInterface const& system,
                      HybridTime const& maximum_hybrid_time) const
{
    typedef Map<DiscreteEvent,ScalarFunction>::const_iterator constraint_iterator;
    typedef Set<DiscreteEvent>::const_iterator event_iterator;

    TaylorIntegrator integrator(5);
    Float step_size=1.0;

    const Float maximum_time=maximum_hybrid_time.continuous_time();
    const uint maximum_steps=maximum_hybrid_time.discrete_time();

    HybridEnclosure starting_set=working_sets.back();
    working_sets.pop_back();

    List<DiscreteEvent>& starting_events=starting_set._events;
    DiscreteLocation starting_location=starting_set._location;
    TaylorConstrainedImageSet& starting_space_set=starting_set._set;
    ScalarTaylorFunction& starting_time=starting_set._time;

    ARIADNE_LOG(4,"starting_events:"<<starting_events<<"\n");
    ARIADNE_LOG(4,"starting_location:"<<starting_location<<"\n");
    ARIADNE_LOG(4,"starting_space_set:"<<starting_space_set<<"\n");
    ARIADNE_LOG(4,"starting_time:"<<starting_time<<"\n");

    if(starting_time.argument_size()>5) {
        std::cerr<< "ABORTING: Too many independent parameters ("<<starting_time.argument_size()<<")\n"; return; }

    // Set the dimension
    const uint n=starting_set.dimension();
    const uint m=starting_set.continuous_state_set().number_of_parameters();

    // Extract mode and transitions, dynamic and constraints
    //DiscreteMode mode=system.mode(starting_location);
    //Set<DiscreteTransition> transitions=system.transitions(starting_location);

    VectorFunction dynamic=system.dynamic_function(starting_location);
    ARIADNE_LOG(4,"dynamic:"<<dynamic<<"\n");

    Set<DiscreteEvent> invariant_events=system.invariant_events(starting_location);
    ARIADNE_LOG(4,"invariant_events:"<<invariant_events<<"\n");
    Set<DiscreteEvent> transition_events=system.transition_events(starting_location);
    ARIADNE_LOG(4,"transition_events:"<<transition_events<<"\n");
    Set<DiscreteEvent> blocking_events=system.blocking_events(starting_location);
    ARIADNE_LOG(4,"blocking_events:"<<blocking_events<<"\n");
    Set<DiscreteEvent> urgent_events=system.urgent_events(starting_location);
    ARIADNE_LOG(4,"urgent_events:"<<urgent_events<<"\n");
    Set<DiscreteEvent> permissive_events=system.permissive_events(starting_location);
    ARIADNE_LOG(4,"permissive_events:"<<permissive_events<<"\n");
    ARIADNE_LOG(4,"\n");

    // Compute flow and actual time step size used
    IntervalVector flow_spacial_domain(starting_set.continuous_state_set().bounding_box());
    ARIADNE_LOG(4,"flow_spacial_domain="<<flow_spacial_domain<<"\n");
    VectorTaylorFunction flow_model=integrator.flow(dynamic,flow_spacial_domain,maximum_time);
    ARIADNE_LOG(4,"twosided_flow_model:"<<flow_model<<"\n");
    Box flow_domain=Box(flow_model.domain());

    // Computed flowed set
    flow_domain[n]=Interval(0,step_size);
    flow_model=restrict(flow_model,flow_domain);
    ARIADNE_LOG(4,"flow_model:"<<flow_model<<"\n");
    ARIADNE_LOG(4,"flow_model.domain():"<<flow_model.domain()<<"\n");
    ARIADNE_LOG(4,"flow_model.range():"<<flow_model.range()<<"\n");
    VectorFunction flow=flow_model.function();
    ARIADNE_LOG(4,"flow:"<<flow<<"\n");
    ARIADNE_LOG(4,"starting_set:"<<starting_set<<"\n");
    HybridEnclosure reached_set=starting_set;
    reached_set.apply_flow(flow,flow_domain[n]);
    ARIADNE_LOG(4,"flowed_set:"<<reached_set<<"\n");

    // Compute restrictions on continuous evolution
    Set<DiscreteEvent> invariant_and_guard_events=join(blocking_events,urgent_events);
    for(event_iterator iter=invariant_and_guard_events.begin(); iter!=invariant_and_guard_events.end(); ++iter) {
        DiscreteEvent const& event=*iter;
        ScalarFunction const& constraint=system.invariant_function(starting_location,event);
        ScalarFunction derivative=lie_derivative(constraint,dynamic);
        reached_set.new_invariant(event,constraint,derivative);
    }
    ARIADNE_LOG(4,"reached_set:"<<reached_set<<"\n");

/*
    // Compute the set reached after one time step of the flow; corresponds to setting t=t0+delta
    // By "finishing" we mean completing one full time step without a jump; this corresponds to delta=h.
    // By "final", we mean completing the full evolution; this corresponds to t=t_max.
    ARIADNE_LOG(4,"starting_time:"<<starting_time<<"\n");
    ScalarTaylorFunction remaining_time_function=maximum_time-starting_time;
    ARIADNE_LOG(4,"remaining_time_function:"<<remaining_time_function<<"\n");
    ScalarTaylorFunction starting_time_function=embed(starting_time,1u);
    ARIADNE_LOG(4,"starting_time_function:"<<starting_time_function<<"\n");
    ScalarTaylorFunction step_time_function=ScalarTaylorFunction::constant(m,step_size);
    ARIADNE_LOG(4,"step_time_function:"<<step_time_function<<"\n");
    ScalarTaylorFunction dwell_time_function=ScalarFunction::coordinate(m+1u,m);
    ARIADNE_LOG(4,"dwell_time_function:"<<dwell_time_function<<"\n");
    ScalarTaylorFunction evolution_time_function=starting_time_function+dwell_time_function;
    ARIADNE_LOG(4,"evolution_time_function:"<<evolution_time_function<<"\n");
    ScalarTaylorFunction progress_constraint=dwell_time_function-step_size;
    ARIADNE_LOG(4,"progress_constraint:"<<progress_constraint<<"\n");
    ScalarTaylorFunction final_constraint=evolution_time_function-maximum_time;
    ARIADNE_LOG(4,"final_constraint:"<<final_constraint<<"\n");
    ARIADNE_LOG(4,""<<"\n");
*/

    // Set the reached set, the next working set and the final set
    HybridEnclosure final_set(reached_set);

    reached_set.set_maximum_time(final_time_event,maximum_time);
    ARIADNE_LOG(4,"reached_set:"<<reached_set<<"\n");
    ARIADNE_LOG(4,"final_set:"<<reached_set<<"\n");
    final_set.set_time(final_time_event,maximum_time);
    ARIADNE_LOG(4,"final_set:"<<final_set<<"\n");

    HybridEnclosure progress_set(starting_set);
    progress_set.apply_flow(flow_model,step_size);
    ARIADNE_LOG(4,"progress_set:"<<progress_set<<"\n");

    reach_sets.adjoin(reached_set);
    ARIADNE_LOG(4,"reach_sets.size():"<<reach_sets.size()<<"\n");
    if(!final_set.empty()) { evolve_sets.adjoin(final_set); }
    ARIADNE_LOG(4,"evolve_sets.size():"<<evolve_sets.size()<<"\n");
    if(!progress_set.empty() && progress_set._time.range().lower()<=maximum_time) {
        working_sets.append(progress_set); }
    ARIADNE_LOG(4,"working_sets.size():"<<working_sets.size()<<"\n");

    // Compute the reached set under a single event
    if(starting_events.size()<maximum_steps) {
        for(event_iterator iter=transition_events.begin(); iter!=transition_events.end(); ++iter) {
            DiscreteEvent const& event=*iter;
            ScalarFunction const& constraint=system.guard_function(starting_location,event);
            ARIADNE_LOG(4,"constraint:"<<constraint<<"\n");
            tribool satisfied=reached_set.continuous_state_set().satisfies(constraint);
            if(possibly(satisfied)) {
                ScalarFunction derivative=lie_derivative(constraint,dynamic);
                ARIADNE_LOG(4,"derivative:"<<derivative<<"\n");
                DiscreteLocation target=system.target(starting_location,event);
                VectorFunction reset=system.reset_function(starting_location,event);
                HybridEnclosure jump_set(reached_set);
                jump_set.new_activation(event,constraint,derivative);
                jump_set.apply_reset(event,target,reset);
                ARIADNE_LOG(4,"jump_set("<<event<<"):"<<jump_set<<"\n");
                working_sets.append(jump_set);
            }
        }
    }

    ARIADNE_LOG(4,"DONE STEP\n"<<"\n");

    if(verbosity==1) {
        ARIADNE_LOG(1,"\n\r"
                    <<"#w="<<std::setw(4)<<std::left<<working_sets.size()
                    <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                    <<" s="<<std::setw(3)<<std::left<<starting_events.size()
                    //<<" t="<<std::setw(7)<<std::fixed<<starting_time.value()
                    <<" t=["<<std::setw(7)<<std::left<<std::fixed<<starting_time.range().lower()
                    <<","<<std::setw(7)<<std::left<<std::fixed<<starting_time.range().upper()<<"]"
                    <<" a="<<std::setw(3)<<std::left<<starting_set.continuous_state_set().domain().size()
                    <<" r="<<std::setw(7)<<starting_set.continuous_state_set().radius()
                    <<" l="<<std::left<<starting_location
                    <<" c="<<starting_set.continuous_state_set().centre()
                    <<" e="<<starting_events
                    <<"                      ");
    }

    return;
}

} // namespace Ariadne
