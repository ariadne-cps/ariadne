/***************************************************************************
 *            development_hybrid_evolver.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
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

#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "vector.h"
#include "matrix.h"
#include "expression_interface.h"
#include "function_interface.h"
#include "taylor_set.h"
#include "taylor_expression.h"
#include "taylor_function.h"
#include "taylor_model.h"
#include "orbit.h"
#include "taylor_calculus.h"
#include "evolution_parameters.h"

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver.h"

namespace {

using namespace Ariadne;

template<class V, class Iter> inline
void append(V& v, Iter begin, Iter end)
{
    for(;begin!=end;++begin) {
        v.push_back(*begin);
    }
}

template<class V, class C> inline
void append(V& v, const C& c)
{
    for(typename C::const_iterator iter=c.begin(); iter!=c.end(); ++iter) {
        v.push_back(*iter);
    }
}


}



namespace Ariadne {

void wait_for_keypress() {
    std::string str;
    //getline(std::cin,str);
}


// Allow subdivisions in upper evolution
const bool ENABLE_SUBDIVISIONS = false;
// Allow premature termination of lower evolution
const bool ENABLE_PREMATURE_TERMINATION = false;

static const int BLOCKING_EVENT = -2;
using boost::shared_ptr;

class DegenerateCrossingException : public std::runtime_error {
  public:
    DegenerateCrossingException(const char* msg) : std::runtime_error(msg) { }
};

const DiscreteEvent HybridEvolver::starting_event = -1;
const DiscreteEvent HybridEvolver::finishing_event = -2;
const DiscreteEvent HybridEvolver::blocking_event = -3;
const DiscreteEvent HybridEvolver::final_time_event = -4;

HybridEvolver::HybridEvolver()
    : _parameters(new EvolutionParametersType()),
      _toolbox(new TaylorCalculus())
{
}



HybridEvolver::HybridEvolver(const EvolutionParametersType& p)
    : _parameters(new EvolutionParametersType(p)),
      _toolbox(new TaylorCalculus)
{
}


Orbit<HybridEvolver::EnclosureType>
HybridEvolver::
orbit(const SystemType& system,
      const EnclosureType& initial_set,
      const TimeType& time,
      Semantics semantics) const
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,
                     system,initial_set,time,semantics,false);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}



namespace{

enum PredicateKind { INVARIANT, ACTIVATION, GUARD, TIME, MIXED };
enum CrossingKind { POSITIVE, NEGATIVE, TRANSVERSE, TOUCHING, NONE, UNKNOWN };


struct DetectionData {
    int id;
    DiscreteEvent event;
    shared_ptr<const FunctionInterface> guard_ptr;
    PredicateKind predicate_kind;
    CrossingKind crossing_kind;
    tribool active;
    tribool initially_active;
    tribool finally_active;
    Interval touching_time_interval;
    TaylorModel crossing_time_model;

    DetectionData()
        : id(0), event(0), guard_ptr()
        , predicate_kind(), crossing_kind(UNKNOWN)
        , active(indeterminate), initially_active(indeterminate), finally_active(indeterminate)
        , touching_time_interval(+1,-1), crossing_time_model() { }
};


std::ostream& operator<<(std::ostream& os, const PredicateKind& kind) {
    switch(kind) {
    case INVARIANT: return os << "INVARIANT";
    case ACTIVATION: return os << "ACTIVATION";
    case GUARD: return os << "GUARD";
    case TIME: return os << "TIME";
    case MIXED: return os << "MIXED";
    }
    return os << "INVALID";
}

std::ostream& operator<<(std::ostream& os, const CrossingKind& kind) {
    switch(kind) {
    case TRANSVERSE: case POSITIVE: case NEGATIVE:
        return os << "TRANSVERSE";
        //case CROSSING: return os << "CROSSING";
        //case GRAZING: return os << "GRAZING";
    case TOUCHING: return os << "TOUCHING";
        //case ACTIVE: return os << "ACTIVE";
        //case MISSING: return os << "MISSING";
    case NONE: return os << "NONE";
    case UNKNOWN: return os << "UNKNOWN";
    }
    return os << "INVALID";
}

std::ostream& operator<<(std::ostream& os, const DetectionData& data) {
    os << std::boolalpha;
    os << " { id="<<data.id;
    os << ", event="<<data.event;
    os << ", kind="<<data.predicate_kind;
    os << ", crossing="<<data.crossing_kind;
    os << ", active="<<data.active;
    if(data.crossing_kind!=NONE) {
        os << ", initially_active="<<data.initially_active;
        os << ", finally_active="<<data.finally_active;
        if(data.touching_time_interval.lower()<=data.touching_time_interval.upper()) {
            os << ", touching_time_interval="<<data.touching_time_interval; }
        if(data.crossing_kind==TRANSVERSE) {
            os << ", crossing_time_model="<<data.crossing_time_model; }
    }
    return os <<" } ";
}

inline bool operator<(const DetectionData& d1, const DetectionData& d2) {
    return d1.id < d2.id;
}

} // namespace


void
HybridEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const SystemType& system,
           const EnclosureType& initial_set,
           const TimeType& maximum_hybrid_time,
           Semantics semantics,
           bool reach) const
{
    typedef boost::shared_ptr< const FunctionInterface > FunctionConstPointer;

    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    const IntegerType maximum_steps=maximum_hybrid_time.discrete_time;
    const Float maximum_time=maximum_hybrid_time.continuous_time;


    typedef tuple<DiscreteState, EventListType, SetModelType, TimeModelType> HybridTimedSetType;

    std::vector< HybridTimedSetType > working_sets;

    {
        // Set up initial timed set models
        ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
        DiscreteState initial_location;
        ContinuousEnclosureType initial_continuous_set;
        make_lpair(initial_location,initial_continuous_set)=initial_set;
        ARIADNE_LOG(6,"initial_location = "<<initial_location<<"\n");
        SetModelType initial_set_model=this->_toolbox->set_model(initial_continuous_set);
        ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
        TimeModelType initial_time_model=this->_toolbox->time_model(0.0,initial_set_model.argument_size());
        ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
        TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
        ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
        working_sets.push_back(make_tuple(initial_location,EventListType(),initial_set_model,initial_time_model));
    }


    while(!working_sets.empty()) {
        HybridTimedSetType current_set=working_sets.back();
        working_sets.pop_back();
        ARIADNE_LOG(3,"working_set="<<current_set);
        DiscreteState initial_location=current_set.first;
        EventListType initial_events=current_set.second;
        SetModelType initial_set_model=current_set.third;
        TimeModelType initial_time_model=current_set.fourth;
        RealType initial_set_radius=radius(initial_set_model.range());
        if(initial_time_model.range().lower()>=maximum_time || initial_events.size()>=uint(maximum_steps)) {
            final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else if(semantics == UPPER_SEMANTICS && ENABLE_SUBDIVISIONS
                  && (initial_set_radius>this->_parameters->maximum_enclosure_radius)) {
            // Subdivide
            uint nd=initial_set_model.dimension();
            SetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
            array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
                SetModelType subdivided_set_model=Vector<TaylorModel>(project(subdivided_timed_set_model.models(),range(0,nd)));
                TimeModelType subdivided_time_model=subdivided_timed_set_model[nd];
                working_sets.push_back(make_tuple(initial_location,initial_events,subdivided_set_model,subdivided_time_model));
            }
        } else if(semantics == LOWER_SEMANTICS && ENABLE_PREMATURE_TERMINATION && initial_set_radius>this->_parameters->maximum_enclosure_radius) {
            std::cerr << "WARNING: Terminating lower evolution at time " << initial_time_model
                      << " and set " << initial_set_model << " due to maximum radius being exceeded.";
        } else {
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  system,current_set,maximum_hybrid_time,
                                  semantics,reach);
        }

        if(verbosity==1) {
            ARIADNE_LOG(1,"\r"
                        <<"#w="<<std::setw(4)<<working_sets.size()
                        <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                        <<" s="<<std::setw(3)<<std::left<<initial_events.size()
                        <<" t="<<std::setw(7)<<std::fixed<<initial_time_model.value()
                        <<" r="<<std::setw(7)<<initial_set_model.radius()
                        <<" c="<<initial_set_model.centre()
                        <<" e="<<initial_events
                        <<"                      ");
        }

    }

}



// Old evolution step using detection data
void
HybridEvolver::
_evolution_step(std::vector< HybridTimedSetType >& working_sets,
                EnclosureListType& final_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const SystemType& system,
                const HybridTimedSetType& working_set,
                const TimeType& maximum_hybrid_time,
                Semantics semantics,
                bool reach) const
{
    // Use the following basic algorithm for computing an evolution step
    //
    // 1) Extract data about working set and location
    //
    // 2) Find all blocking events (invariants and urgent transitions) which are
    //    active at the starting time. If any events are definitely active,
    //    the corresponding reset occurs, and there is no continuous evolution.
    //    If there are possibly active events, these occur along with continuous
    //    evolution for upper semantics, and evolution terminates with lower
    //    semantics.
    //
    //    Input: Set: starting_set, map<Event,Expression> guards
    //    Output: map<Event,Tribool> initially_active
    //         (includes special "blocking event")
    //
    // 3) Compute the continuous evolution for a fixed step size h, over the
    //    time interval [-h,+h]
    //
    //    Input: Function: dynamic,
    //           Box: domain OR SetModel initial_set,
    //           Time: maximum_step_size
    //    Output: FunctionModel: flow_model
    //               OR SetModel: flow_set_model
    //            Box: flow_bounds
    //
    // 4) For each blocking event, compute the crossing time with the guard set.
    //    The computed crossing time may lie outside the flow time interval;
    //    any such crossings will be ignored. Non-transverse crossings may have
    //    large crossing time intervals.
    //
    //    If the transition is definitely not initially active and the crossing
    //    is transverse, then there are no problems. If the transition is
    //    possibly initially active, then the crossing may be in the "wrong"
    //    direction, i.e. the transition may become inactive. In this case, we
    //    have upper semantics (otherwise we would already have terminated the
    //    evolution) and the transition is considered inactive for evolution
    //    purposes.
    //
    //    If the crossing is not transverse, then lower evolution is blocking
    //    and upper evolution considers the transition as non-urgent. However,
    //    we should probably evolve close to the transition, and maybe split
    //    the evolution across the transversality boundary.
    //
    //    Input: map<Event,Expression>: resets,
    //           Box: flow_bounds OR FunctionModel: flow_model
    //             OR SetModel: flow_set_model
    //    Output: map<Event,ExpressionModel> crossing_times,
    //            set<Event> tangential_events
    //
    // 5) Compute the blocking time and blocking events. The blocking time is
    //    the minimum of the computed crossing times.
    //
    //    If there is a single  blocking event, the crossing is tranverse, and
    //    occurs between the starting and finishing times, then evolution
    //    proceeds according to this event.
    //
    //    If there are multiple blocking events, and the upper bound of the
    //    crossing time range is large, then we set the finishing time to
    //    just below the crossing time. This means that in the next step
    //    we may be better able to resolve the crossing.
    //
    //    If there are multiple blocking events and the upper bound of the
    //    blocking time range is small, then lower evolution terminates, and
    //    upper evolution proceeds according to the crossing time of each
    //    blocking event.
    //
    //    Input: map<Event,ExpressionModel>: crossing_times
    //    Output: set<Event> blocking_events, ExpressionModel: blocking_time
    //
    // 6) Compute the initial and final activation times of the non-blocking
    //    events. Tangential crossings are included in this computation,
    //    as they are treated as non-urgent.
    //
    //    Compute the maximum of the initial activation time and starting time,
    //    and the minimum of the final activation time and finishing/blocking
    //    time.
    //
    //    Input: map<Event,Expression>: guards,
    //           Box: domain OR FunctionModel: flow_model
    //               OR SetModel: flow_set_model
    //    Output: map<Event,(ExpressionModel,ExpressionModel)>: initial/final times
    //
    // 7) Apply the flows, guards and resets according to the computed
    //    event times.

    // More useful typedefs
    typedef boost::shared_ptr<const ExpressionInterface> ExpressionPtr;
    typedef boost::shared_ptr<const FunctionInterface> FunctionPtr;

    const double SMALL_RELATIVE_TIME=1./16;

    // Extract information about the working set
    DiscreteState location(0);
    IntegerType steps;
    EventListType events;
    SetModelType set_model;
    TimeModelType time_model;
    ARIADNE_LOG(9,"working_set = "<<working_set<<"\n");
    make_ltuple(location,events,set_model,time_model)=working_set;
    steps=events.size();

    // Extract information about the current location
    const DiscreteMode& mode=system.mode(location);
    const FunctionPtr dynamic_ptr=mode.dynamic_ptr();
    const std::map<DiscreteEvent,FunctionPtr> guards=system.blocking_guards(location);
    std::map<DiscreteEvent,FunctionPtr> activations=system.permissive_guards(location);
    const std::set<DiscreteTransition> transitions=system.transitions(location);

    // Check to make sure dimensions are correct
    ARIADNE_ASSERT(set_model.argument_size()==time_model.argument_size());
    ARIADNE_ASSERT(set_model.result_size()==mode.dimension());

    ARIADNE_LOG(2,"\n\nHybridEvolver::_evolution_step(...)\n");
    ARIADNE_LOG(2,"events = "<<events<<" ");
    ARIADNE_LOG(2,"time_range = "<<time_model.range()<<" ");
    ARIADNE_LOG(2,"time_model generators = "<<time_model.argument_size()<<" ");
    ARIADNE_LOG(2,"location = "<<location<<" ");
    ARIADNE_LOG(2,"box = "<<set_model.range()<<" ");
    ARIADNE_LOG(2,"generators = "<<set_model.argument_size()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(set_model.range())<<"\n\n");

    ARIADNE_LOG(2,"dynamic = "<<*dynamic_ptr<<"\n");
    ARIADNE_LOG(2,"invariants = "<<mode.invariants()<<"\n");
    ARIADNE_LOG(2,"transitions = "<<transitions<<"\n\n");
    ARIADNE_LOG(2,"guards = "<<guards<<"\n\n");
    ARIADNE_LOG(2,"activations = "<<activations<<"\n\n");



    // Compute initially active guards
    std::map<DiscreteEvent,tribool> initially_active_events;
    this->compute_initially_active_events(initially_active_events, guards, set_model);
    ARIADNE_LOG(2,"initially_active_events = "<<initially_active_events<<"\n\n");

    // Test for initially active events, and process these as requred
    if(definitely(initially_active_events[blocking_event])) {
        // No continuous evolution; just apply discrete events
        for(std::set<DiscreteTransition>::const_iterator iter=transitions.begin();
            iter!=transitions.end(); ++iter)
        {
            tribool active=this->_toolbox->active(iter->activation(),set_model);
            if(definitely(active) || (possibly(active) && UPPER_SEMANTICS)) {
                SetModelType jump_set_model=this->_toolbox->reset_step(iter->reset(),set_model);
                TimeModelType jump_time_model=time_model;
                DiscreteState jump_location=iter->target().location();
                std::vector<DiscreteEvent> jump_events=events;
                jump_events.push_back(iter->event());
                ARIADNE_LOG(2,"Pushing back "<<jump_set_model<<" into working_sets\n");
                working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,jump_time_model));
            }
        }
        // Adjoin reach and intermediate sets
        reach_sets.insert(make_pair(location,set_model));
        intermediate_sets.insert(make_pair(location,set_model));
        // No need for further work this step
        return;

    } else if(possibly(initially_active_events[blocking_event]) && semantics==LOWER_SEMANTICS) {
        ARIADNE_LOG(2,"Terminating lower evolution due to possibly initially active invariant or urgent transition.");
        reach_sets.insert(make_pair(location,set_model));
        intermediate_sets.insert(make_pair(location,set_model));
        return;
    }

    // Compute continuous evolution
    FlowSetModelType flow_set_model; BoxType flow_bounds; Float time_step;
    compute_flow_model(flow_set_model,flow_bounds,time_step,dynamic_ptr,set_model);

    ARIADNE_LOG(2,"flow_bounds = "<<flow_bounds<<"\n")
    ARIADNE_LOG(2,"time_step = "<<time_step<<"\n")
    ARIADNE_LOG(2,"flow_range = "<<flow_set_model.range()<<"\n");
    ARIADNE_LOG(2,"starting_set_range = "<<set_model.range()<<"\n");
    // Partial evaluation on flow set model to obtain final set must take scaled time equal to 1.0
    SetModelType finishing_set=partial_evaluate(flow_set_model.models(),set_model.argument_size(),1.0);
    ARIADNE_LOG(2,"finishing_set_range = "<<finishing_set.range()<<"\n")

    // Set special events and times; note that the time step is scaled to [0,1]
    TimeModelType zero_time_model = this->_toolbox->time_model(0.0,time_model.argument_size());
    TimeModelType time_step_model = this->_toolbox->time_model(1.0,time_model.argument_size());
    TimeModelType remaining_time_model = (maximum_hybrid_time.continuous_time-time_model)/time_step;
    ARIADNE_LOG(2,"remaining_time = "<<remaining_time_model.range()<<"\n\n")

    // Compute event blocking times
    std::map<DiscreteEvent, TimeModelType> event_blocking_times;
    std::set<DiscreteEvent> non_transverse_events;


    event_blocking_times[finishing_event]=time_step_model;
    event_blocking_times[final_time_event]=remaining_time_model;
    compute_blocking_events(event_blocking_times,non_transverse_events,guards,flow_set_model);
    ARIADNE_LOG(3,"event_blocking_times="<<event_blocking_times<<"\n");

    std::map<DiscreteEvent,Interval> event_blocking_time_intervals;
    for(std::map<DiscreteEvent, TimeModelType>::const_iterator iter=event_blocking_times.begin();
        iter!=event_blocking_times.end(); ++iter) { event_blocking_time_intervals[iter->first]=iter->second.range(); }

    ARIADNE_LOG(2,"event_blocking_times="<<event_blocking_time_intervals<<"\n");
    ARIADNE_LOG(2,"non_transverse_events="<<non_transverse_events<<"\n\n");

    // Compute blocking events
    std::set<DiscreteEvent> blocking_events;
    TimeModelType blocking_time_model;
    compute_blocking_time(blocking_events,blocking_time_model,event_blocking_times);
    ARIADNE_LOG(2,"blocking_events="<<blocking_events<<"\n");
    ARIADNE_LOG(2,"blocking_time="<<blocking_time_model.range()<<"\n\n");

    // If multiple blocking events and flow time is large, move closer to blocking time
    if(blocking_events.size()!=1 && blocking_time_model.range().upper()>SMALL_RELATIVE_TIME) {
        blocking_events.clear();
        blocking_events.insert(finishing_event);
        if(blocking_time_model.range().lower()>0) {
            blocking_time_model-=blocking_time_model.error();
            blocking_time_model.set_error(0.0);
        }
        blocking_time_model*=(1-SMALL_RELATIVE_TIME);
    }
    ARIADNE_LOG(2,"blocking_events="<<blocking_events<<"\n");
    ARIADNE_LOG(2,"blocking_time="<<blocking_time_model.range()<<"\n\n");


    // Treat non-transverse urgent events as non-urgent in upper semantics
    for(std::set<DiscreteEvent>::const_iterator iter=non_transverse_events.begin(); iter!=non_transverse_events.end(); ++iter) {
        if(*iter > 0) {     // If the event is a transition
            activations[*iter]=guards.find(*iter)->second;
       } 
    }

    // Compute activation times
    std::map< DiscreteEvent,tuple<TimeModelType,TimeModelType> > activation_times;
    this->compute_activation_times(activation_times,activations,flow_set_model,blocking_time_model,semantics);

    // Display activation time ranges
    std::map< DiscreteEvent,tuple<Interval> > activation_time_intervals;
    for(std::map< DiscreteEvent, tuple<TimeModelType,TimeModelType> >::const_iterator
            iter=activation_times.begin(); iter!=activation_times.end(); ++iter)
    {
        activation_time_intervals.insert(make_pair(iter->first,make_tuple(iter->second.second.range())));
    }
    ARIADNE_LOG(2,"activation_times="<<activation_time_intervals<<"\n\n");



    // Compute sets
    // TODO: Make this a function;

    ARIADNE_LOG(4,"flow_set_model="<<flow_set_model<<"\n");
    ARIADNE_LOG(4,"zero_time_model="<<zero_time_model<<"\n");
    ARIADNE_LOG(4,"blocking_time_model="<<blocking_time_model<<"\n");
    SetModelType reachable_set=this->_toolbox->reachability_step(flow_set_model,zero_time_model,blocking_time_model);
    reach_sets.adjoin(make_pair(location,reachable_set));
    ARIADNE_LOG(2,"reachable_set.argument_size()="<<reachable_set.argument_size()<<"\n");
    ARIADNE_LOG(2,"reachable_set.range()="<<reachable_set.range()<<"\n");
    if(semantics==LOWER_SEMANTICS && blocking_events.size()!=1) {
        // No further evolution
    } else {
        TimeModelType final_time_model=time_model+blocking_time_model*time_step;
        ARIADNE_LOG(2,"final_time_range="<<final_time_model.range()<<"\n");
        SetModelType evolved_set_model=this->_toolbox->integration_step(flow_set_model,blocking_time_model);
        ARIADNE_LOG(2,"evolved_set_model.argument_size()="<<evolved_set_model.argument_size()<<"\n");
        ARIADNE_LOG(2,"evolved_set_range="<<evolved_set_model.range()<<"\n");
        // Compute evolution for blocking events
        for(std::set<DiscreteEvent>::const_iterator iter=blocking_events.begin(); iter!=blocking_events.end(); ++iter) {
            const DiscreteEvent event=*iter;
            if(event==finishing_event) {
                // TODO: Better estimate to use smaller blocking time
                intermediate_sets.adjoin(make_pair(location,evolved_set_model));
                working_sets.push_back(make_tuple(location,events,evolved_set_model,final_time_model));
            } else if(event==final_time_event) {
                final_sets.adjoin(make_pair(location,evolved_set_model));
            } else if(event>=0) { // not an invariant
                intermediate_sets.adjoin(make_pair(location,evolved_set_model));
                const DiscreteTransition& transition=system.transition(event,location);
                SetModelType jump_set_model=apply(*transition.reset_ptr(),evolved_set_model);
                DiscreteState jump_location=transition.target().location();
                std::vector<DiscreteEvent> jump_events=events;
                jump_events.push_back(event);
                working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,final_time_model));
            }
        }
        // Compute evolution for non-blocking events
        for(std::map< DiscreteEvent, tuple<TimeModelType,TimeModelType> >::const_iterator 
            iter=activation_times.begin(); iter!=activation_times.end(); ++iter) 
        {
            const DiscreteEvent event=iter->first;
            const TimeModelType lower_active_time_model=iter->second.first;
            const TimeModelType upper_active_time_model=iter->second.second;
            ARIADNE_LOG(3,"Non blocking event "<<event<<":\n");
            ARIADNE_LOG(3,"  lower_active_time_model="<<lower_active_time_model.range()<<";\n");
            ARIADNE_LOG(3,"  upper_active_time_model="<<upper_active_time_model.range()<<";\n");
            SetModelType active_set_model=this->_toolbox->reachability_step(flow_set_model,lower_active_time_model,upper_active_time_model);
            ARIADNE_LOG(3,"  active_set="<<active_set_model.range()<<";\n");
            const DiscreteTransition& transition=system.transition(event,location);
            SetModelType jump_set_model=apply(*transition.reset_ptr(),active_set_model);
            ARIADNE_LOG(3,"  jump_set_model="<<active_set_model.range()<<";\n");
            const TimeModelType active_time_model = this->_toolbox->reachability_time(time_model+lower_active_time_model*time_step,time_model+upper_active_time_model*time_step);
            ARIADNE_LOG(3,"  active_time_model="<<active_time_model.range()<<".\n");
            DiscreteState jump_location=transition.target().location();
            std::vector<DiscreteEvent> jump_events=events;
            jump_events.push_back(event);
            working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,active_time_model));
        } 
    }

    wait_for_keypress();
}


tribool HybridEvolver::
active(FunctionPtr guard_ptr, const SetModelType& set) const
{
    typedef TimeModelType GuardValueModelType;
    GuardValueModelType guard_set_model = apply(*guard_ptr,set)[0];
    Interval guard_range=guard_set_model.range();
    tribool guard_initially_active=guard_range.lower()>0 ? tribool(true) : guard_range.upper()<0 ? tribool(false) : indeterminate;
    return guard_initially_active;
}


HybridEvolver::TimeModelType HybridEvolver::
crossing_time(FunctionPtr guard_ptr, const FlowSetModelType& flow_set_model) const
{
    try {
        TimeModelType crossing_time_model=this->_toolbox->scaled_crossing_time(*guard_ptr,flow_set_model);
        return crossing_time_model;
    }
    catch(DegenerateCrossingException e) {
        BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
        Interval touching_time_interval=this->_toolbox->scaled_touching_time_interval(*guard_ptr,flow_set_model);
        TimeModelType touching_time_model=this->_toolbox->time_model(touching_time_interval,space_domain);
        return touching_time_model;
    } // end non-transverse crossing
}


Interval jacobian2_range(const TaylorModel& tm);

Interval HybridEvolver::
normal_derivative(FunctionPtr guard_ptr, const FlowSetModelType& flow_set_model, const TimeModelType& crossing_time_model) const
{
    typedef TimeModelType GuardValueModelType;
    GuardValueModelType guard_flow_set_model=apply(*guard_ptr,flow_set_model)[0];
    Interval normal_derivative=jacobian2_range(guard_flow_set_model);
    return normal_derivative;
}


void HybridEvolver::
compute_initially_active_events(std::map<DiscreteEvent,tribool>& initially_active_events,
                                const std::map<DiscreteEvent,FunctionPtr>& guards,
                                const ContinuousEnclosureType& initial_set) const
{
    typedef TimeModelType GuardValueModelType;
    tribool blocking_event_initially_active=false;
    for(std::map<DiscreteEvent,FunctionPtr>::const_iterator iter=guards.begin(); iter!=guards.end(); ++iter) {
        tribool initially_active=this->active(iter->second,initial_set);
        if(possibly(initially_active)) {
            initially_active_events.insert(std::make_pair(iter->first,initially_active));
            blocking_event_initially_active = blocking_event_initially_active || initially_active;
        }
    }
    initially_active_events.insert(std::make_pair(blocking_event,blocking_event_initially_active));
    return;
}

// Compute the flow, parameterising space with the set parameters
void HybridEvolver::
compute_flow_model(FlowSetModelType& flow_set_model, BoxType& flow_bounds, Float& time_step,
                   FunctionPtr dynamic_ptr, const SetModelType& starting_set_model) const
{
    const int MAXIMUM_BOUNDS_DIAMETER_FACTOR = 8;
    const Float maximum_step_size=this->_parameters->maximum_step_size;
    const Float maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*MAXIMUM_BOUNDS_DIAMETER_FACTOR;

    BoxType starting_set_bounding_box=starting_set_model.range();

    make_lpair(time_step,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,starting_set_bounding_box,maximum_step_size,maximum_bounds_diameter);
    // Compute the flow model
    FlowModelType flow_model=this->_toolbox->flow_model(*dynamic_ptr,starting_set_bounding_box,time_step,flow_bounds);
    TaylorExpression identity_time_expression=TaylorExpression::variable(BoxType(1u,Interval(-time_step,+time_step)),0u);
    flow_set_model=unchecked_apply(flow_model,combine(starting_set_model.models(),identity_time_expression.model()));
}


void HybridEvolver::
compute_flow_model(FunctionModelType& , BoxType&,
                   FunctionPtr, const BoxType&) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


void HybridEvolver::
compute_blocking_events(std::map<DiscreteEvent,TimeModelType>& event_blocking_times,
                        std::set<DiscreteEvent>& non_transverse_events,
                        const std::map<DiscreteEvent,FunctionPtr>& guards,
                        const FlowSetModelType& flow_set_model) const
{
    uint dimension=flow_set_model.result_size();
    const double SMALL_RELATIVE_TIME = 1./16;
    FlowSetModelType positive_flow_set_model(split(flow_set_model.models(),dimension,1));

    for(std::map<DiscreteEvent,FunctionPtr>::const_iterator guard_iter=guards.begin();
        guard_iter!=guards.end(); ++guard_iter)
    {
        const DiscreteEvent event=guard_iter->first;
        const FunctionPtr guard_ptr=guard_iter->second;
        tribool active = this->_toolbox->active(*guard_ptr,positive_flow_set_model);
        if(possibly(active)) {
            TimeModelType crossing_time_model;
            Interval normal_derivative;
            try {
                crossing_time_model=this->_toolbox->scaled_crossing_time(*guard_ptr,flow_set_model);
                normal_derivative=this->normal_derivative(guard_ptr,flow_set_model,crossing_time_model);
                assert(normal_derivative.lower()>0 || normal_derivative.upper()<0);
                if(normal_derivative.lower()>0) {
                    event_blocking_times[event]=crossing_time_model;
                }
            }
            catch(DegenerateCrossingException e) {
                BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
                Interval touching_time_interval=this->_toolbox->scaled_touching_time_interval(*guard_ptr,flow_set_model);
                TimeModelType touching_time_model=this->_toolbox->time_model(touching_time_interval,space_domain);
                // Use 1.0 as upper bound above since flow set model has time interval normalised to [-1,+1]
                if(touching_time_interval.upper()>=0 && touching_time_interval.lower()<=1.0) {
                    SetModelType finishing_set_model=partial_evaluate(flow_set_model.models(),dimension,1.0);
                    tribool finishing_set_active=this->_toolbox->active(*guard_ptr,finishing_set_model);
                    if(definitely(finishing_set_active)) {
                        event_blocking_times[event]=touching_time_model;
                    } else if(possibly(finishing_set_active)) {
                        if(touching_time_interval.lower()>SMALL_RELATIVE_TIME) {
                            TaylorModel lower_touching_time_model=
                                this->_toolbox->time_model(touching_time_interval.lower(),space_domain);
                            event_blocking_times[finishing_event]=lower_touching_time_model;
                        } else {
                            // FIXME: Here we are stuck, we can't determine whether the crossing is completely finished or not.
                            // Just put this in as a blocking event and hope for the best...
                            event_blocking_times[event]=touching_time_model;
                        }
                    } else {
                        // After the flow step, the event is not active again, so the crossing was tangential
                        non_transverse_events.insert(event);
                    }
                }
            } // end non-transverse crossing
        } // end possibly active
    } // end main loop

    return;
}


void HybridEvolver::
compute_blocking_time(std::set<DiscreteEvent>& blocking_events,
                      TimeModelType& blocking_time,
                      const std::map<DiscreteEvent,TimeModelType>& event_blocking_times) const
{
    assert(!event_blocking_times.empty());
    blocking_events.insert(event_blocking_times.begin()->first);
    blocking_time=event_blocking_times.begin()->second;

    for(std::map<DiscreteEvent,TimeModelType>::const_iterator iter=++event_blocking_times.begin();
        iter!=event_blocking_times.end(); ++iter)
    {
        const DiscreteEvent event=iter->first;
        const TimeModelType& time=iter->second;
        tribool is_first=(time<blocking_time);
        if(definitely(is_first)) {
            blocking_events.clear();
            blocking_events.insert(event);
            blocking_time=time;
        } else if(possibly(is_first)) {
            std::set<DiscreteEvent> new_blocking_events;
            new_blocking_events.insert(event);
            TimeModelType new_blocking_time=time;
            for(std::set<DiscreteEvent>::const_iterator iter=blocking_events.begin();
                iter!=blocking_events.end(); ++iter)
            {
                const TimeModelType& current_event_time=event_blocking_times.find(*iter)->second;
                if(possibly(current_event_time<time)) {
                    new_blocking_events.insert(*iter);
                    blocking_time=min(current_event_time,blocking_time);
                }
            }
            blocking_time.swap(new_blocking_time);
            blocking_events.swap(new_blocking_events);
        }
    }
    return;
}


void HybridEvolver::
compute_activation_events(std::map<DiscreteEvent,tuple<tribool,TimeModelType,tribool> >& activation_events,
                          const std::map<DiscreteEvent,FunctionPtr>& activations, const FlowSetModelType& flow_set_model) const
{
    SetModelType initial_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,0.0);
    SetModelType final_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,1.0);
    for(std::map<DiscreteEvent,FunctionPtr>::const_iterator iter=activations.begin(); iter!=activations.end(); ++iter) {
        DiscreteEvent event=iter->first;
        FunctionPtr activation_ptr=iter->second;
        tribool active=this->active(activation_ptr,flow_set_model);
        if(possibly(active)) {
            tribool initially_active=this->active(activation_ptr,initial_set_model);
            tribool finally_active=this->active(activation_ptr,final_set_model);
            TimeModelType crossing_time_model=this->crossing_time(activation_ptr,flow_set_model);
            activation_events.insert(make_pair(event,make_tuple(initially_active,crossing_time_model,finally_active)));
        }
    }

}


void HybridEvolver::
compute_activation_times(std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >& activation_times,
                         const std::map<DiscreteEvent,FunctionPtr>& activations,
                         const FlowSetModelType& flow_set_model,
                         const TimeModelType& blocking_time_model,
                         const Semantics semantics) const
{
    SetModelType initial_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,0.0);
    SetModelType final_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,1.0);
    TimeModelType zero_time_model=blocking_time_model*0.0;

    for(std::map<DiscreteEvent,FunctionPtr>::const_iterator iter=activations.begin(); iter!=activations.end(); ++iter) {
        DiscreteEvent event=iter->first;
        FunctionPtr activation_ptr=iter->second;

        ARIADNE_LOG(3,"Computing activation time for event "<<event<<"...");

        // Compute whether the event might be enabled on the entire time interval
        tribool active=this->active(activation_ptr,flow_set_model);

        if(definitely(active)) {
            // The event is enabled over the entire time interval
            ARIADNE_LOG(3,"event is enabled over the entire time interval.\n");
            activation_times.insert(make_pair(event,make_tuple(zero_time_model,blocking_time_model)));
        } else if(possibly(active)) {
            ARIADNE_LOG(3,"event is possibly enabled.\n");
            // Compute whether the event is enabled at the beginning and end of the time interval
            tribool initially_active=this->active(activation_ptr,initial_set_model);
            tribool finally_active=this->active(activation_ptr,final_set_model);

            TimeModelType crossing_time_model=this->crossing_time(activation_ptr,flow_set_model);

            TimeModelType lower_crossing_time_model=crossing_time_model-crossing_time_model.error();
            TimeModelType upper_crossing_time_model=crossing_time_model+crossing_time_model.error();
            lower_crossing_time_model.set_error(0);
            upper_crossing_time_model.set_error(0);

            TimeModelType lower_active_time_model, upper_active_time_model;

            // Determine lower activation time
            if(definitely(not(initially_active))) {
                ARIADNE_LOG(3,"event is definitely not initially active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: lower_active_time_model=lower_crossing_time_model; break;
                    case LOWER_SEMANTICS: lower_active_time_model=upper_crossing_time_model; break;
                }
            } else if(definitely(initially_active)) {
                ARIADNE_LOG(3,"event is definitely initially active.\n");
                lower_active_time_model=zero_time_model;
            } else {
                ARIADNE_LOG(3,"Event is undeteriminate initially active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: lower_active_time_model=zero_time_model; break;
                    case LOWER_SEMANTICS: lower_active_time_model=upper_crossing_time_model; break;
                }
            }

            // Compute upper activation time
            if(definitely(not(finally_active))) {
                ARIADNE_LOG(3,"Event is definitely not finally active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: upper_active_time_model=upper_crossing_time_model; break;
                    case LOWER_SEMANTICS: upper_active_time_model=lower_crossing_time_model; break;
                }
            } else if(definitely(finally_active)) {
                ARIADNE_LOG(3,"Event is definitely finally active.\n");
                upper_active_time_model=blocking_time_model;
            } else {
                ARIADNE_LOG(3,"Event is undeteriminate finally active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: upper_active_time_model=blocking_time_model; break;
                    case LOWER_SEMANTICS: upper_active_time_model=upper_crossing_time_model; break;
                }
            }

            // In case of lower semantics, event may not be
            switch(semantics) {
                case UPPER_SEMANTICS:
                    activation_times.insert(make_pair(event,make_tuple(lower_active_time_model,upper_active_time_model)));
                    break;
                case LOWER_SEMANTICS:
                    if(lower_active_time_model<upper_active_time_model) {
                        activation_times.insert(make_pair(event,make_tuple(lower_active_time_model,upper_active_time_model)));
                    }
                    break;
            }
        }
    }

}

/*
{
    TimeType::ContinuousTimeType maximum_time=maximum_hybrid_time.continuous_time;
    TimeModelType maximum_integration_time_model=maximum_time-initial_time_model;

    /////////////// Main Evolution ////////////////////////////////

    const DiscreteMode& initial_mode=system.mode(initial_location);
    const FunctionType* dynamic_ptr=&initial_mode.dynamic();

    ARIADNE_LOG(6,"mode="<<initial_mode<<"\n");
    std::vector< shared_ptr<const FunctionInterface> > invariants
        =initial_mode.invariants();
    ARIADNE_LOG(7,"invariants="<<invariants<<"\n");
    const std::set< DiscreteTransition > transitions = system.transitions(initial_location);
    ARIADNE_LOG(7,"transitions="<<transitions<<"\n");



    // Set evolution parameters
    const Float maximum_step_size=this->_parameters->maximum_step_size;
    const Float maximum_bounds_diameter=this->_parameters->maximum_enclosure_radius*2;
    const Float zero_time=0.0;

    // Get bounding boxes for time and space range
    Vector<Interval> initial_set_bounds=initial_set_model.range();
    ARIADNE_LOG(4,"initial_set_range = "<<initial_set_bounds<<"\n");
    Interval initial_time_range=initial_time_model.range();
    ARIADNE_LOG(4,"initial_time_range = "<<initial_time_range<<"\n");

    //ARIADNE_ASSERT(initial_time_range.width() <= maximum_step_size);

    // Compute flow bounds and find flow bounding box
    Vector<Interval> flow_bounds;
    Float step_size;
    make_lpair(step_size,flow_bounds)=this->_toolbox->flow_bounds(*dynamic_ptr,initial_set_model.range(),maximum_step_size,maximum_bounds_diameter);
    ARIADNE_LOG(4,"step_size = "<<step_size<<"\n");
    ARIADNE_LOG(4,"flow_bounds = "<<flow_bounds<<"\n");

    // Compute the flow model
    FlowModelType flow_model=this->_toolbox->flow_model(*dynamic_ptr,initial_set_bounds,step_size,flow_bounds);
    ARIADNE_LOG(6,"flow_model = "<<flow_model<<"\n");

    //ARIADNE_ASSERT_MSG(subset(flow_model.range(),flow_bounds),flow_model<<"\n  "<<flow_model.range()<<"\n  "<<flow_bounds);

    // Compute the integration time model
    TimeModelType final_time_model=initial_time_model+step_size;
    ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
    TimeModelType integration_time_model=final_time_model-initial_time_model;

    // Compute the flow tube (reachable set) model and the final set
    SetModelType final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
    ARIADNE_LOG(6,"final_set_model = "<<final_set_model<<"\n");
    SetModelType reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    ARIADNE_LOG(4,"Done computing continuous evolution\n");









    // Compute transitions

    ARIADNE_LOG(2,"Evolution step.\n");

    ARIADNE_LOG(7,"initial_time_bounds="<<initial_time_model.range()<<"\n");
    ARIADNE_LOG(7,"flow_bounds="<<flow_bounds<<"\n");
    ARIADNE_LOG(7,"intermediate_flow_bounds="<<reach_set_model.range()<<"\n");
    ARIADNE_LOG(7,"initial_flow_bounds="<<initial_set_model.range()<<"\n");
    ARIADNE_LOG(7,"final_flow_bounds="<<final_set_model.range()<<"\n");
    ARIADNE_LOG(7,"time_interval="<<Interval(0,step_size)<<"\n");

    std::map<int,DetectionData> detection_data;
    ARIADNE_LOG(9,"invariants="<<invariants<<"\n");
    for(std::vector< shared_ptr<const FunctionInterface> >::const_iterator
            iter=invariants.begin(); iter!=invariants.end(); ++iter)
        {
            DetectionData new_data;
            new_data.id=detection_data.size();
            new_data.event=BLOCKING_EVENT;
            new_data.predicate_kind=INVARIANT;
            new_data.crossing_kind=UNKNOWN;
            new_data.guard_ptr=(*iter);
            detection_data.insert(make_pair(new_data.id,new_data));
        }

    ARIADNE_LOG(9,"transitions="<<transitions<<"\n");
    for(std::set< DiscreteTransition >::const_iterator
            iter=transitions.begin(); iter!=transitions.end(); ++iter)
        {
            DetectionData new_data;
            new_data.id=detection_data.size();
            new_data.event=iter->event();
            new_data.predicate_kind=iter->forced() ? GUARD : ACTIVATION;
            new_data.crossing_kind=UNKNOWN;
            new_data.guard_ptr=iter->activation_ptr();
            detection_data.insert(make_pair(new_data.id,new_data));
        }

    ARIADNE_LOG(8,"detection_data="<<detection_data<<"\n");

    typedef std::map<int,DetectionData>::iterator predicate_iterator;

    // Write out the ranges of the
    ARIADNE_LOG(6,"time="<<initial_time_model.value()<<"\n");
    ARIADNE_LOG(6,"centre="<<initial_set_model.centre()<<" radius="<<initial_set_model.radius()<<"\n");
    ARIADNE_LOG(6,"step_size="<<step_size<<"\n");
    ARIADNE_LOG(6,"vector="<<dynamic_ptr->evaluate(initial_set_model.centre())<<"\n");
    for(predicate_iterator iter=detection_data.begin();
        iter!=detection_data.end(); ++iter)
    {
        ARIADNE_LOG(6,"event: "<<iter->first<<
                      " guard: "<<*iter->second.guard_ptr<<
                      " guard value: "<<iter->second.guard_ptr->evaluate(initial_set_model.centre())<<"\n");
    }

    wait_for_keypress();


    for(predicate_iterator iter=detection_data.begin();
        iter!=detection_data.end(); ++iter)
        {
            DetectionData& data=iter->second;
            if(possibly(this->_toolbox->active(*data.guard_ptr,flow_bounds))) {
                data.active=this->_toolbox->active(*data.guard_ptr,reach_set_model);
                data.initially_active=data.active;
                data.finally_active=data.active;
                if(data.active) {
                    data.crossing_kind=NONE;
                } else if(!data.active) {
                    data.crossing_kind=NONE;
                } else {
                    data.initially_active=this->_toolbox->active(*data.guard_ptr,initial_set_model);
                    data.finally_active=this->_toolbox->active(*data.guard_ptr,final_set_model);
                    ConstraintModelType guard_model=this->_toolbox->predicate_model(*data.guard_ptr,flow_bounds);
                    if(data.initially_active ^ data.finally_active) {
                        ARIADNE_LOG(7," Testing crossing time for: "<<data<<"\n");
                        try {
                            data.crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);
                            data.touching_time_interval=data.crossing_time_model.range();
                            data.crossing_kind=TRANSVERSE;
                        }
                        catch(DegenerateCrossingException e) { std::cerr<<e.what()<<"\n"; ARIADNE_LOG(7," DegenerateCrossing\n"); }
                    }
                    if(data.crossing_kind!=TRANSVERSE) {
                        data.touching_time_interval=this->_toolbox->touching_time_interval(guard_model,flow_model,initial_set_model);
                        data.crossing_kind=TOUCHING;
                    }
                }
            } else {
                data.active=false;
                data.crossing_kind=NONE;
            }
        }

    ARIADNE_LOG(6,"detection_data="<<detection_data<<"\n");

    // Remove inactive transitions
    for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ) {
        if(!iter->second.active) {
            detection_data.erase(iter++);
        } else {
            ++iter;
        }
    }

    ARIADNE_LOG(6,"active_detection_data="<<detection_data<<"\n");








    // Compute the evolution time for the blocking evolution
    ARIADNE_LOG(6,"Testing for blocking predicate\n");
    DetectionData blocking_data;
    blocking_data.predicate_kind=TIME;
    blocking_data.crossing_kind=NONE;
    blocking_data.touching_time_interval=step_size;
    blocking_data.active=false;
    blocking_data.initially_active=false;
    blocking_data.finally_active=false;
    for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ++iter) {
        DetectionData& data=iter->second;
        if(data.predicate_kind==INVARIANT || data.predicate_kind==GUARD) {
            if(data.initially_active) {
                blocking_data.predicate_kind=data.predicate_kind;
                blocking_data.active=indeterminate;
                blocking_data.initially_active=true;
                blocking_data.finally_active=indeterminate;
                blocking_data.crossing_kind=NONE;
                ARIADNE_LOG(6,"No continuous evolution possible\n");
                break;
            } else if(data.finally_active) {
                blocking_data.predicate_kind=data.predicate_kind;
                if(data.touching_time_interval.lower()<blocking_data.touching_time_interval.upper()) {
                    blocking_data.active=indeterminate;
                    if(data.touching_time_interval.upper()<blocking_data.touching_time_interval.lower() && data.crossing_kind==TRANSVERSE) {
                        blocking_data.predicate_kind=data.predicate_kind;
                        blocking_data.crossing_kind=TRANSVERSE;
                        blocking_data.finally_active=true;
                        blocking_data.touching_time_interval=data.touching_time_interval;
                        blocking_data.crossing_time_model=data.crossing_time_model;
                        blocking_data.event=data.event;
                    } else {
                        blocking_data.predicate_kind=MIXED;
                        blocking_data.crossing_kind=TOUCHING;
                        blocking_data.initially_active=blocking_data.initially_active || data.initially_active;
                        blocking_data.finally_active=blocking_data.initially_active || data.finally_active;
                        blocking_data.touching_time_interval=min(data.touching_time_interval,blocking_data.touching_time_interval);
                    }
                }
            }
        }
    }

    ARIADNE_LOG(6,"blocking_data="<<blocking_data<<"\n");

    // FIXME: Take final time into account here
    Float maximum_evolution_time;
    if(semantics==UPPER_SEMANTICS) {
        maximum_evolution_time=blocking_data.touching_time_interval.upper();
    } else {
        maximum_evolution_time=blocking_data.touching_time_interval.lower();
    }


    if(maximum_evolution_time>=step_size) {
        ARIADNE_LOG(6,"No restrictions on continuous evolution\n");
        maximum_evolution_time=step_size;
    }

    // Compute blocking time for time step
    if(blocking_data.touching_time_interval.lower()>step_size/16) {
        // If the blocking event is non-transverse and occurs after too long an interval,
        // just do a time step
        blocking_data.predicate_kind=TIME;
        blocking_data.crossing_time_model=this->_toolbox->time_model(blocking_data.touching_time_interval.lower(),initial_time_model.domain());
    } else if(blocking_data.predicate_kind==TIME) {
        ARIADNE_LOG(6,"  Continuous evolution for maximum time; unwinding time differentce\n");
        Interval initial_time_range=initial_time_model.range();
        Float final_time=initial_time_range.lower()+step_size;
        ARIADNE_LOG(8,"    initial_time_model="<<initial_time_model<<"\n");
        final_time_model=this->_toolbox->time_model(final_time,initial_time_model.domain());
        TimeModelType integration_time_model=final_time_model-initial_time_model;
        ARIADNE_LOG(8,"    integration_time_model="<<integration_time_model<<"\n");
        if(integration_time_model>maximum_integration_time_model) { integration_time_model=maximum_integration_time_model; }
        blocking_data.crossing_time_model=integration_time_model;
   }


    ARIADNE_LOG(6,"maximum_evolution_time="<<maximum_evolution_time<<"\n");


    // Remove transitions which only occur after the blocking time
    for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ) {
        if(!iter->second.initially_active && iter->second.touching_time_interval.lower()>maximum_evolution_time) {
            detection_data.erase(iter++);
        } else {
            ++iter;
        }
    }


    ARIADNE_LOG(6,"active_detection_data="<<detection_data<<"\n");

    // Compute continuous evolution
    ARIADNE_LOG(6,"Computing continuous evolution\n");
    if( blocking_data.initially_active || ( possibly(blocking_data.initially_active) && semantics==LOWER_SEMANTICS ) ) {
        ARIADNE_LOG(6,"  No continuous evolution\n");
        reach_set_model=initial_set_model;
        reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
    } else if(blocking_data.crossing_kind==TRANSVERSE) {
        ARIADNE_LOG(6,"  Continuous evolution to transverse invariant/guard\n");
        reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,blocking_data.crossing_time_model);
        reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
    } else if(blocking_data.finally_active) {
        ARIADNE_LOG(6,"  Continuous evolution for finite time to invariant/guard\n");
        if(semantics==UPPER_SEMANTICS) {
            reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,blocking_data.touching_time_interval.upper());
            final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,blocking_data.touching_time_interval.upper());
            reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
        } else { // semantics==LOWER_SEMANTICS
            reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,blocking_data.touching_time_interval.lower());
            final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,blocking_data.touching_time_interval.lower());
            reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
        }
    } else {
        ARIADNE_ASSERT(blocking_data.predicate_kind==TIME);
        TimeModelType integration_time_model=blocking_data.crossing_time_model;
        reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);
        final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
        reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
        intermediate_sets.adjoin(EnclosureType(initial_location,final_set_model));
        working_sets.push_back(make_tuple(initial_location,initial_events,final_set_model,final_time_model));
    }
    ARIADNE_LOG(6,"initial_set_model="<<initial_set_model<<"\n");
    ARIADNE_LOG(6,"flow_model="<<flow_model<<"\n");
    ARIADNE_LOG(6,"integration_time_model="<<integration_time_model<<"\n");
    ARIADNE_LOG(6,"final_set_model="<<final_set_model<<"\n");

    wait_for_keypress();

    // Process discrete transitions
    ARIADNE_LOG(6,"Computing discrete transitions\n");
    for(predicate_iterator iter=detection_data.begin(); iter!=detection_data.end(); ++iter) {
        DetectionData& data=iter->second;
        ARIADNE_LOG(6,"  transition"<<data<<"\n");
        if(data.predicate_kind==ACTIVATION || data.predicate_kind==GUARD) {
            EventListType jump_events=initial_events;
            jump_events.push_back(data.event);
            const DiscreteTransition& transition=system.transition(data.event,initial_location);
            DiscreteState jump_location=transition.target().location();
            const FunctionInterface* reset_ptr=&transition.reset();
            TimeModelType active_time_model;
            TimeModelType jump_time_model;

            if( definitely(blocking_data.initially_active) || ( possibly(blocking_data.initially_active) && semantics==LOWER_SEMANTICS ) ) {
                // In this case there is no continuous evolution
                SetModelType jump_set_model
                    =this->_toolbox->reset_step(*reset_ptr,initial_set_model);
                if( ( semantics==UPPER_SEMANTICS && possibly(data.initially_active) ) || definitely(data.initially_active) ) {
                    working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,initial_time_model));
                    intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
                }
                reach_sets.adjoin(EnclosureType(initial_location,initial_set_model));
            } else {
                // In this case may be continuous evolution
                if(data.crossing_kind==TRANSVERSE) {
                    if(data.initially_active) {
                        active_time_model=this->_toolbox->reachability_time(zero_time,data.crossing_time_model);
                        jump_time_model=this->_toolbox->reachability_time(initial_time_model,initial_time_model+data.crossing_time_model);
                    } else if(data.predicate_kind==GUARD) {
                        // Not sure if this code is needed.... maybe we can use crossing time as before
                        active_time_model=this->_toolbox->reachability_time(data.crossing_time_model,maximum_evolution_time);
                        jump_time_model=initial_time_model+data.crossing_time_model;
                    } else {
                        active_time_model=this->_toolbox->reachability_time(data.crossing_time_model,maximum_evolution_time);
                        jump_time_model=this->_toolbox->
                            reachability_time(initial_time_model+data.crossing_time_model,
                                              initial_time_model+maximum_evolution_time);
                    }
                } else {
                    Float lower_active_time, upper_active_time;
                    if(semantics==UPPER_SEMANTICS) {
                        lower_active_time = data.initially_active ? zero_time : data.touching_time_interval.lower();
                        upper_active_time = data.finally_active ? step_size : data.touching_time_interval.upper();
                    } else { // semantics==LOWER_SEMANTICS
                        lower_active_time = definitely(data.initially_active) ? zero_time : data.touching_time_interval.upper();
                        upper_active_time = definitely(data.finally_active) ? step_size : data.touching_time_interval.lower();
                    }
                    TimeModelType lower_active_time_model=this->_toolbox->time_model(lower_active_time,initial_time_model.domain());
                    TimeModelType upper_active_time_model=this->_toolbox->time_model(upper_active_time,initial_time_model.domain());
                    active_time_model=this->_toolbox->reachability_time(lower_active_time_model,upper_active_time_model);
                    jump_time_model=this->_toolbox->reachability_time(initial_time_model+lower_active_time_model,initial_time_model+upper_active_time_model);
                }

                SetModelType active_set_model;
                SetModelType jump_set_model;
                if(data.crossing_kind==TRANSVERSE) {
                    active_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,data.crossing_time_model);
                    jump_set_model=this->_toolbox->reset_step(*reset_ptr,active_set_model);
                } else {
                    active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,active_time_model);
                    jump_set_model=this->_toolbox->reset_step(*reset_ptr,active_set_model);
                }

                if(semantics==UPPER_SEMANTICS || data.predicate_kind==ACTIVATION
                   || data.crossing_kind==TRANSVERSE) {
                    working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,jump_time_model));
                }
                intermediate_sets.adjoin(EnclosureType(initial_location,active_set_model));
                intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
            }
        }
    }

    ARIADNE_LOG(2,"Done evolution_step.\n\n");
    wait_for_keypress();

}
*/



}  // namespace Ariadne

