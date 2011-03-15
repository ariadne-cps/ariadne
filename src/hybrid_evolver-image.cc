/***************************************************************************
 *            hybrid_evolver-image.cc
 *
 *  Copyright  2008-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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
#include "function.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "taylor_model.h"
#include "orbit.h"
#include "taylor_calculus.h"
#include "evolution_parameters.h"

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver-image.h"

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
    getline(std::cin,str);
}


static const int BLOCKING_EVENT = -2;

class DegenerateCrossingException : public std::runtime_error {
  public:
    DegenerateCrossingException(const char* msg) : std::runtime_error(msg) { }
};

const DiscreteEvent ImageSetHybridEvolver::starting_event = -1;
const DiscreteEvent ImageSetHybridEvolver::finishing_event = -2;
const DiscreteEvent ImageSetHybridEvolver::blocking_event = -3;

ImageSetHybridEvolver::ImageSetHybridEvolver()
    : _parameters(new EvolutionParametersType()),
      _toolbox(new TaylorCalculus())
{
	
}


ImageSetHybridEvolver::ImageSetHybridEvolver(const EvolutionParametersType& p)
    : _parameters(new EvolutionParametersType(p)),
      _toolbox(new TaylorCalculus())
{
}

ImageSetHybridEvolver::ImageSetHybridEvolver(const EvolutionParametersType& p,
											 const TaylorCalculus& tc)
    : _parameters(new EvolutionParametersType(p)),
      _toolbox(new TaylorCalculus(tc))
{
}


Orbit<ImageSetHybridEvolver::EnclosureType>
ImageSetHybridEvolver::
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
                     system,initial_set,time,semantics);
    orbit.adjoin_intermediate(intermediate);
    orbit.adjoin_reach(reachable);
    orbit.adjoin_final(final);
    return orbit;
}


Orbit<ImageSetHybridEvolver::EnclosureType>
ImageSetHybridEvolver::
upper_orbit_continuous(const SystemType& system,
				       const EnclosureType& initial_set,
				       const TimeType& time)
{
    Orbit<EnclosureType> orbit(initial_set);
    EnclosureListType final;
    EnclosureListType reachable;
    EnclosureListType intermediate;

    this->_upper_evolution_continuous(final,reachable,intermediate,
                     				  system,initial_set,time);

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
    shared_ptr<const VectorFunction> guard_ptr;
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
ImageSetHybridEvolver::
_upper_evolution_continuous(EnclosureListType& final_sets,
           			  		EnclosureListType& reach_sets,
				            EnclosureListType& intermediate_sets,
           			  	 	const SystemType& system,
           			  		const EnclosureType& initial_set,
           			  		const TimeType& maximum_hybrid_time) const
{
    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n"); 
    ARIADNE_LOG(1,"Computing evolution up to "<<maximum_hybrid_time.continuous_time()<<" time units and "<<maximum_hybrid_time.discrete_time()<<" steps.\n");

	// The working sets, pushed back and popped front
    std::list<HybridTimedSetType> working_sets;
    _evolution_add_initialSet(working_sets,initial_set);

	// While there exists a working set, process it and increment the total
    uint i=0;
	while(!working_sets.empty()) {

		ARIADNE_LOG(2,"\n");
		ARIADNE_LOG(2,"Processed sets: " << i++ << ", remaining sets: " << working_sets.size() << "\n\n");

		// Get the least recent working set and pop it
		HybridTimedSetType current_set = working_sets.front(); 
		working_sets.pop_front();

		// Possibly simplifies the current set
		_set_model_reduction(current_set);

        DiscreteState initial_location=current_set.first;
        EventListType initial_events=current_set.second;  	
		SetModelType initial_set_model=current_set.third;
		TimeModelType initial_time_model=current_set.fourth;		

		bool isEnclosureTooLarge = _isEnclosureTooLarge(initial_set_model);

        if(initial_time_model.range().lower()>=maximum_hybrid_time.continuous_time()) {
            ARIADNE_LOG(3,"\n\nFinal time reached: adjoining the initial set to the final sets.\n\n");
            final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else if(isEnclosureTooLarge) {
            ARIADNE_LOG(1,"\n\nTerminating evolution at time " << initial_time_model.value()
                        << " and set centre " << initial_set_model.centre() << " due to the enclosure being too large: adjoining the initial set to the final sets.\n\n");
            final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else {
            this->_upper_evolution_continuous_step(working_sets,reach_sets,intermediate_sets,
                                  				   system,current_set,maximum_hybrid_time);
        }

        _logStepAtVerbosity1(working_sets,reach_sets,initial_events,initial_time_model,initial_set_model,initial_location);
    }

}

DisproveData
ImageSetHybridEvolver::
_lower_evolution_disprove(EnclosureListType& final_sets,
						  EnclosureListType& reach_sets,
						  EnclosureListType& intermediate_sets,
						  const SystemType& system,
						  const EnclosureType& initial_set,
						  const TimeType& maximum_hybrid_time,
						  const HybridBoxes& disprove_bounds,
						  bool enable_quick_disproving) const
{
    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");
    ARIADNE_LOG(1,"Computing evolution up to "<<maximum_hybrid_time.continuous_time()<<" time units and "<<maximum_hybrid_time.discrete_time()<<" steps.\n");

    // The disprovement information for the result
    DisproveData disproveData(system.state_space());

	// The working sets, pushed back and popped front
    std::list< HybridTimedSetType > working_sets;
    _evolution_add_initialSet(working_sets,initial_set);

	// While there exists a working set, process it and increment the total
    uint i=0;
	while(!working_sets.empty()) {

		ARIADNE_LOG(2,"\n");
		ARIADNE_LOG(2,"Processed sets: " << i++ << ", remaining sets: " << working_sets.size() << "\n\n");

		// Get the least recent working set, pop it and update the corresponding size
		HybridTimedSetType current_set = working_sets.front();
		working_sets.pop_front();

		// Possibly simplifies the current set
		_set_model_reduction(current_set);

        DiscreteState initial_location=current_set.first;
        EventListType initial_events=current_set.second;
		SetModelType initial_set_model=current_set.third;
		TimeModelType initial_time_model=current_set.fourth;

		bool isEnclosureTooLarge = _isEnclosureTooLarge(initial_set_model);
		bool subdivideOverTime = (initial_time_model.range().width() > this->_parameters->hybrid_maximum_step_size[initial_location]/2);

		if(initial_time_model.range().lower()>=maximum_hybrid_time.continuous_time() ||
		   initial_events.size()>=uint(maximum_hybrid_time.discrete_time())) {
            ARIADNE_LOG(3,"  Final time reached, adjoining result to final sets.\n");
            final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else if (subdivideOverTime && this->_parameters->enable_subdivisions) {
            ARIADNE_LOG(1,"WARNING: computed time range " << initial_time_model.range() << " width larger than half the maximum step size " << this->_parameters->hybrid_maximum_step_size[initial_location] << ", subdividing over time.\n");
            _add_models_subdivisions_time(working_sets,initial_set_model,initial_time_model,initial_location,initial_events);
		} else if(!this->_parameters->enable_subdivisions &&
                  this->_parameters->enable_premature_termination && isEnclosureTooLarge) {
            ARIADNE_LOG(1,"\n\nWARNING: Terminating evolution at time " << initial_time_model.value()
                        << " and set " << initial_set_model.centre() << " due to maximum enclosure bounds being exceeded.\n\n");
        } else {
            // Compute evolution and get the result for this working set
            DisproveData localDisproveData = this->_lower_evolution_disprove_step(working_sets,reach_sets,intermediate_sets,
																			system,current_set,maximum_hybrid_time,
																			disprove_bounds);
            disproveData.updateWith(localDisproveData);
            if (enable_quick_disproving && disproveData.getIsDisproved())
				return disproveData;
        }

		_logStepAtVerbosity1(working_sets,reach_sets,initial_events,initial_time_model,initial_set_model,initial_location);
    }

    return disproveData;
}

void
ImageSetHybridEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const SystemType& system,
           const EnclosureType& initial_set,
           const TimeType& maximum_hybrid_time,
           Semantics semantics) const
{
    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");
    ARIADNE_LOG(1,"Computing evolution up to "<<maximum_hybrid_time.continuous_time()<<" time units and "<<maximum_hybrid_time.discrete_time()<<" steps.\n");

	// The working sets, pushed back and popped front
    std::list< HybridTimedSetType > working_sets;
    _evolution_add_initialSet(working_sets,initial_set);

	// While there exists a working set, process it and increment the total
    uint i=0;
	while(!working_sets.empty()) {

		ARIADNE_LOG(2,"\n");
		ARIADNE_LOG(2,"Processed sets: " << i++ << ", remaining sets: " << working_sets.size() << "\n\n");

		// Get the least recent working set, pop it and update the corresponding size
		HybridTimedSetType current_set = working_sets.front();
		working_sets.pop_front();

		// Possibly simplifies the current set
		_set_model_reduction(current_set);

		// Get the members of the current set
        DiscreteState initial_location=current_set.first;
        EventListType initial_events=current_set.second;
		SetModelType initial_set_model=current_set.third;
		TimeModelType initial_time_model=current_set.fourth;

		bool isEnclosureTooLarge = _isEnclosureTooLarge(initial_set_model);
		bool subdivideOverTime = (initial_time_model.range().width() > this->_parameters->hybrid_maximum_step_size[initial_location]/2);

		if(initial_time_model.range().lower()>=maximum_hybrid_time.continuous_time() ||
		   initial_events.size()>=uint(maximum_hybrid_time.discrete_time())) {
            ARIADNE_LOG(3,"Final time reached, adjoining result to final sets.\n");
            final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else if (subdivideOverTime && this->_parameters->enable_subdivisions) {
            ARIADNE_LOG(1,"WARNING: computed time range " << initial_time_model.range() << " width larger than half the maximum step size " << this->_parameters->hybrid_maximum_step_size[initial_location] << ", subdividing over time.\n");
            _add_models_subdivisions_time(working_sets,initial_set_model,initial_time_model,initial_location,initial_events);
		} else if (semantics == UPPER_SEMANTICS && this->_parameters->enable_subdivisions && isEnclosureTooLarge) {
            ARIADNE_LOG(1,"WARNING: computed set range " << initial_set_model.range() << " widths larger than maximum_enclosure_cell " << this->_parameters->maximum_enclosure_cell << ", subdividing.\n");
            _add_models_subdivisions_autoselect(working_sets,initial_set_model,initial_time_model,initial_location,initial_events);
        } else if((semantics == LOWER_SEMANTICS || !this->_parameters->enable_subdivisions) &&
                  this->_parameters->enable_premature_termination && isEnclosureTooLarge) {
            ARIADNE_LOG(1,"\n\nWARNING: Terminating evolution at time " << initial_time_model.value()
                        << " and set " << initial_set_model.centre() << " due to maximum enclosure bounds being exceeded.\n\n");
            if(semantics == UPPER_SEMANTICS)
                final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else {
            this->_evolution_step(working_sets,reach_sets,intermediate_sets,system,current_set,maximum_hybrid_time,semantics);
        }

		_logStepAtVerbosity1(working_sets,reach_sets,initial_events,initial_time_model,initial_set_model,initial_location);
    }
}


// Old evolution step using detection data
void
ImageSetHybridEvolver::
_evolution_step(std::list< HybridTimedSetType >& working_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const SystemType& system,
                const HybridTimedSetType& working_set,
                const TimeType& maximum_hybrid_time,
                Semantics semantics) const
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
    const VectorFunction dynamic=get_directed_dynamic(mode.dynamic(),_parameters->direction);
    const std::map<DiscreteEvent,VectorFunction> guards=system.blocking_guards(location);
    std::map<DiscreteEvent,VectorFunction> activations=system.permissive_guards(location);
    const std::map<DiscreteEvent,VectorFunction> invariants=mode.invariants();
    const std::list<DiscreteTransition> transitions=system.transitions(location);

    // Check to make sure dimensions are correct
    ARIADNE_ASSERT(set_model.argument_size()==time_model.argument_size());
    ARIADNE_ASSERT_MSG(set_model.result_size()==mode.dimension(),"set_model="<<set_model<<", mode="<<mode);

    _logEvolutionStepInitialState(events,time_model,location,set_model,dynamic,invariants,transitions,guards,activations);

    // Compute initially active guards
    std::map<DiscreteEvent,tribool> initially_active_events;
    this->compute_initially_active_events(initially_active_events, guards, set_model);
    ARIADNE_LOG(2,"initially_active_events = "<<initially_active_events<<"\n\n");

    // Test for initially active events, and process these as required
    if(definitely(initially_active_events[blocking_event])) {
    	const std::map<DiscreteEvent,VectorFunction> invariants_superset = (semantics == UPPER_SEMANTICS ? guards : invariants);
		_processInitiallyActiveBlockingEvents(working_sets,reach_sets,intermediate_sets,location,events,
	        		   	   	   	   	   	  set_model,time_model,invariants_superset,transitions);
        return;

    } else if(possibly(initially_active_events[blocking_event]) && semantics==LOWER_SEMANTICS) {
        ARIADNE_LOG(2,"Terminating lower evolution due to possibly initially active invariant or urgent transition.");
        reach_sets.insert(make_pair(location,set_model));
        intermediate_sets.insert(make_pair(location,set_model));
        return;
    }

    // Compute continuous evolution
    FlowSetModelType flow_set_model; BoxType flow_bounds; 
    Float time_step = this->_parameters->hybrid_maximum_step_size[location];
    const Float maximum_time=maximum_hybrid_time.continuous_time();
    compute_flow_model(flow_set_model,flow_bounds,time_step,dynamic,set_model,time_model,maximum_time);

    ARIADNE_LOG(2,"flow_bounds = "<<flow_bounds<<"\n")
    ARIADNE_LOG(2,"time_step = "<<time_step<<"\n")
    ARIADNE_LOG(2,"flow_range = "<<flow_set_model.range()<<"\n");
    ARIADNE_LOG(2,"starting_set_range = "<<set_model.range()<<"\n");
    // Partial evaluation on flow set model to obtain final set must take scaled time equal to 1.0
    SetModelType finishing_set=partial_evaluate(flow_set_model.models(),set_model.argument_size(),1.0);
    ARIADNE_LOG(2,"finishing_set_range = "<<finishing_set.range()<<"\n")

    // Set special events and times; note that the time step is scaled to [0,1]
    TimeModelType zero_time_model = this->_toolbox->time_model(0.0,Box(time_model.argument_size()));
    TimeModelType time_step_model = this->_toolbox->time_model(1.0,Box(time_model.argument_size()));

    std::set<DiscreteEvent> blocking_events;
    TimeModelType blocking_time_model;
    std::set<DiscreteEvent> non_transverse_events;
    _compute_blocking_info(non_transverse_events,blocking_events,blocking_time_model,
    				  time_step_model,flow_set_model,guards,SMALL_RELATIVE_TIME);

    ActivationTimesType activation_times;
    _compute_activation_info(activations,activation_times,non_transverse_events,flow_set_model,
    						 blocking_time_model,guards,semantics);

    SetModelType reachable_set;
    _compute_and_adjoin_reachableSet(reach_sets,reachable_set,location,flow_set_model,zero_time_model,blocking_time_model);

    if(semantics!=LOWER_SEMANTICS || blocking_events.size()==1)
    	_computeEvolutionForEvents(working_sets,intermediate_sets,system,location,blocking_events,events,
    								activation_times,flow_set_model,time_model,blocking_time_model,time_step);

}


// Performs an evolution step, without the generation of the evolved set due to transitions
void
ImageSetHybridEvolver::
_upper_evolution_continuous_step(std::list< HybridTimedSetType >& working_sets,
                				 EnclosureListType& reach_sets,
                				 EnclosureListType& intermediate_sets,
                				 const SystemType& system,
                				 const HybridTimedSetType& working_set,
                				 const TimeType& maximum_hybrid_time) const
{
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
    const VectorFunction dynamic=get_directed_dynamic(mode.dynamic(),_parameters->direction);
    const std::map<DiscreteEvent,VectorFunction> guards=system.blocking_guards(location);
    std::map<DiscreteEvent,VectorFunction> activations=system.permissive_guards(location);
    const std::map<DiscreteEvent,VectorFunction> invariants=mode.invariants();
    const std::list<DiscreteTransition> transitions=system.transitions(location);

    // Check to make sure dimensions are correct
    ARIADNE_ASSERT(set_model.argument_size()==time_model.argument_size());
    ARIADNE_ASSERT_MSG(set_model.result_size()==mode.dimension(),"set_model="<<set_model<<", mode="<<mode);

    _logEvolutionStepInitialState(events,time_model,location,set_model,dynamic,invariants,transitions,guards,activations);

    // Compute initially active events
    std::map<DiscreteEvent,tribool> initially_active_events;
    this->compute_initially_active_events(initially_active_events, guards, set_model);
    ARIADNE_LOG(2,"initially_active_events = "<<initially_active_events<<"\n\n");

    if(definitely(initially_active_events[blocking_event])) {
    	_processInitiallyActiveBlockingEvents_continuous(reach_sets,intermediate_sets,guards,set_model,location);
        return;
    }

    // Compute continuous evolution
    FlowSetModelType flow_set_model; BoxType flow_bounds; 
    Float time_step = this->_parameters->hybrid_maximum_step_size[location];
    const Float maximum_time=maximum_hybrid_time.continuous_time();
    compute_flow_model(flow_set_model,flow_bounds,time_step,dynamic,set_model,time_model,maximum_time);

    ARIADNE_LOG(2,"flow_bounds = "<<flow_bounds<<"\n")
    ARIADNE_LOG(2,"time_step = "<<time_step<<"\n")
    ARIADNE_LOG(2,"flow_range = "<<flow_set_model.range()<<"\n");
    ARIADNE_LOG(2,"starting_set_range = "<<set_model.range()<<"\n");
    // Partial evaluation on flow set model to obtain final set must take scaled time equal to 1.0
    SetModelType finishing_set=partial_evaluate(flow_set_model.models(),set_model.argument_size(),1.0);
    ARIADNE_LOG(2,"finishing_set_range = "<<finishing_set.range()<<"\n")

    // Set special events and times; note that the time step is scaled to [0,1]
    TimeModelType zero_time_model = this->_toolbox->time_model(0.0,Box(time_model.argument_size()));
    TimeModelType time_step_model = this->_toolbox->time_model(1.0,Box(time_model.argument_size()));

    // Compute blocking information
    std::set<DiscreteEvent> blocking_events;
    TimeModelType blocking_time_model;
    std::set<DiscreteEvent> non_transverse_events;
    _compute_blocking_info(non_transverse_events,blocking_events,blocking_time_model,
    				  time_step_model,flow_set_model,guards,SMALL_RELATIVE_TIME);

    SetModelType reachable_set;
    _compute_and_adjoin_reachableSet(reach_sets,reachable_set,location,flow_set_model,zero_time_model,blocking_time_model);

    TimeModelType final_time_model=time_model+blocking_time_model*time_step;
    ARIADNE_LOG(2,"final_time_range="<<final_time_model.range()<<"\n");
    SetModelType evolved_set_model=this->_toolbox->integration_step(flow_set_model,blocking_time_model);
    ARIADNE_LOG(2,"evolved_set_model.argument_size()="<<evolved_set_model.argument_size()<<"\n");
    ARIADNE_LOG(2,"evolved_set_range="<<evolved_set_model.range()<<"\n");

    // Compute evolution for blocking events
    for(std::set<DiscreteEvent>::const_iterator iter=blocking_events.begin(); iter!=blocking_events.end(); ++iter) {
        const DiscreteEvent event=*iter;
        if(event==finishing_event) {
            intermediate_sets.adjoin(make_pair(location,evolved_set_model));
            working_sets.push_back(make_tuple(location,events,evolved_set_model,final_time_model));
        } else if(event.is_transition()) { // not an invariant
            intermediate_sets.adjoin(make_pair(location,evolved_set_model));
        }
    }
}


DisproveData
ImageSetHybridEvolver::
_lower_evolution_disprove_step(std::list< HybridTimedSetType >& working_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const SystemType& system,
                const HybridTimedSetType& working_set,
                const TimeType& maximum_hybrid_time,
                const HybridBoxes& disprove_bounds) const
{
    const double SMALL_RELATIVE_TIME=1./16;

    // The falsification info result
    DisproveData disproveData(system.state_space());

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
    const VectorFunction dynamic=get_directed_dynamic(mode.dynamic(),_parameters->direction);
    const std::map<DiscreteEvent,VectorFunction> guards=system.blocking_guards(location);
    std::map<DiscreteEvent,VectorFunction> activations=system.permissive_guards(location);
    const std::map<DiscreteEvent,VectorFunction> invariants=mode.invariants();
    const std::list<DiscreteTransition> transitions=system.transitions(location);

    // Check to make sure dimensions are correct
    ARIADNE_ASSERT(set_model.argument_size()==time_model.argument_size());
    ARIADNE_ASSERT_MSG(set_model.result_size()==mode.dimension(),"set_model="<<set_model<<", mode="<<mode);

    _logEvolutionStepInitialState(events,time_model,location,set_model,dynamic,invariants,transitions,guards,activations);

    // Compute initially active guards
    std::map<DiscreteEvent,tribool> initially_active_events;
    this->compute_initially_active_events(initially_active_events, guards, set_model);
    ARIADNE_LOG(2,"initially_active_events = "<<initially_active_events<<"\n\n");

    // If there is any possibly (thus including definitely) active blocking event, lower evolution must be stopped
    if (possibly(initially_active_events[blocking_event])) {
    	if(definitely(initially_active_events[blocking_event])) {
    		_processInitiallyActiveBlockingEvents(working_sets,reach_sets,intermediate_sets,location,events,
    	        		   	   	   	   	   	  set_model,time_model,invariants,transitions);
    	} else {
    		ARIADNE_LOG(2,"Terminating lower evolution due to possibly initially active invariant or urgent transition.");
            reach_sets.insert(make_pair(location,set_model));
            intermediate_sets.insert(make_pair(location,set_model));

            disproveData.updateReachBounds(location,set_model.bounding_box());
            disproveData.updateEpsilon(location,set_model.bounding_box().halfWidths());
    	}

        return disproveData;
    }

    // Compute continuous evolution
    FlowSetModelType flow_set_model; BoxType flow_bounds;
    Float time_step = this->_parameters->hybrid_maximum_step_size[location];
    const Float maximum_time=maximum_hybrid_time.continuous_time();
    compute_flow_model(flow_set_model,flow_bounds,time_step,dynamic,set_model,time_model,maximum_time);

    ARIADNE_LOG(2,"flow_bounds = "<<flow_bounds<<"\n")
    ARIADNE_LOG(2,"time_step = "<<time_step<<"\n")
    ARIADNE_LOG(2,"flow_range = "<<flow_set_model.range()<<"\n");
    ARIADNE_LOG(2,"starting_set_range = "<<set_model.range()<<"\n");
    // Partial evaluation on flow set model to obtain final set must take scaled time equal to 1.0
    SetModelType finishing_set=partial_evaluate(flow_set_model.models(),set_model.argument_size(),1.0);
    ARIADNE_LOG(2,"finishing_set_range = "<<finishing_set.range()<<"\n")

    // Set special events and times; note that the time step is scaled to [0,1]
    TimeModelType zero_time_model = this->_toolbox->time_model(0.0,Box(time_model.argument_size()));
    TimeModelType time_step_model = this->_toolbox->time_model(1.0,Box(time_model.argument_size()));

    std::set<DiscreteEvent> blocking_events;
    TimeModelType blocking_time_model;
    std::set<DiscreteEvent> non_transverse_events;
    _compute_blocking_info(non_transverse_events,blocking_events,blocking_time_model,
    				  time_step_model,flow_set_model,guards,SMALL_RELATIVE_TIME);

    ActivationTimesType activation_times;
    _compute_activation_info(activations,activation_times,non_transverse_events,flow_set_model,
    						 blocking_time_model,guards,LOWER_SEMANTICS);

    SetModelType reachable_set;
    _compute_and_adjoin_reachableSet(reach_sets,reachable_set,location,flow_set_model,zero_time_model,blocking_time_model);

    // Check the reachable set against the disprove_bounds (up to numDivisions times)
    const uint numDivisions = 10;
    HybridBoxes::const_iterator bounds_it = disprove_bounds.find(location);
    ARIADNE_LOG(2,"Disproving reachable set against box " << bounds_it->second << " ... \n");
    const bool isDisproved = definitely(_is_reachableSet_outside_disproveBounds(numDivisions,reachable_set,bounds_it->second));
    ARIADNE_LOG(2,"Result: " << isDisproved << "\n");

    // Updates all the fields of the falsification info
    disproveData.updateIsDisproved(isDisproved);
    disproveData.updateReachBounds(location,reachable_set.bounding_box());
    disproveData.updateEpsilon(location,reachable_set.bounding_box().halfWidths());

    if(blocking_events.size()==1)
    	_computeEvolutionForEvents(working_sets,intermediate_sets,system,location,blocking_events,events,
    								activation_times,flow_set_model,time_model,blocking_time_model,time_step);

    return disproveData;
}


ImageSetHybridEvolver::TimeModelType ImageSetHybridEvolver::
crossing_time(VectorFunction guard, const FlowSetModelType& flow_set_model) const
{
    try {
        TimeModelType crossing_time_model=this->_toolbox->scaled_crossing_time(guard,flow_set_model);
        return crossing_time_model;
    }
    catch(DegenerateCrossingException e) {
        BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
        Interval touching_time_interval=this->_toolbox->scaled_touching_time_interval(guard,flow_set_model);
        TimeModelType touching_time_model=this->_toolbox->time_model(touching_time_interval,space_domain);
        return touching_time_model;
    } // end non-transverse crossing
}


Interval jacobian2_range(const TaylorModel& tm);

Interval ImageSetHybridEvolver::
normal_derivative(VectorFunction guard, const FlowSetModelType& flow_set_model, const TimeModelType& crossing_time_model) const
{
    typedef TimeModelType GuardValueModelType;
    GuardValueModelType guard_flow_set_model=apply(guard,flow_set_model)[0];
    Interval normal_derivative=jacobian2_range(guard_flow_set_model);
    return normal_derivative;
}


void ImageSetHybridEvolver::
compute_initially_active_events(std::map<DiscreteEvent,tribool>& initially_active_events,
                                const std::map<DiscreteEvent,VectorFunction>& guards,
                                const ContinuousEnclosureType& initial_set) const
{
    typedef TimeModelType GuardValueModelType;
    tribool blocking_event_initially_active=false;
    for(std::map<DiscreteEvent,VectorFunction>::const_iterator iter=guards.begin(); iter!=guards.end(); ++iter) {
        VectorFunction activation=iter->second;
        tribool initially_active=this->_toolbox->active(activation,initial_set);
        if(possibly(initially_active)) {
            initially_active_events.insert(std::make_pair(iter->first,initially_active));
            blocking_event_initially_active = blocking_event_initially_active || initially_active;
        }
    }
    initially_active_events.insert(std::make_pair(blocking_event,blocking_event_initially_active));
    return;
}

// Compute the flow, parameterising space with the set parameters
void ImageSetHybridEvolver::
compute_flow_model(FlowSetModelType& flow_set_model, BoxType& flow_bounds, Float& time_step,
                   VectorFunction dynamic, const SetModelType& starting_set_model, 
                   const TimeModelType& starting_time_model, Float finishing_time) const
{
    ARIADNE_LOG(3,"compute_flow_model(....)\n");
    const int MAXIMUM_BOUNDS_DIAMETER_FACTOR = 8;
    float remaining_time = finishing_time - starting_time_model.range().lower();
    const Float maximum_step_size=min(time_step, remaining_time);
    const Float maximum_bounds_diameter=max(this->_parameters->maximum_enclosure_cell)*MAXIMUM_BOUNDS_DIAMETER_FACTOR;

    BoxType starting_set_bounding_box=starting_set_model.range();
    ARIADNE_LOG(3,"starting_set_bounding_box="<<starting_set_bounding_box<<"\n");
    make_lpair(time_step,flow_bounds)=this->_toolbox->flow_bounds(dynamic,starting_set_bounding_box,maximum_step_size,maximum_bounds_diameter);
    // Compute the flow model
    ARIADNE_LOG(3,"time_step="<<time_step<<"\n");
    ARIADNE_LOG(3,"flow_bounds="<<flow_bounds<<"\n");
    FlowModelType flow_model=this->_toolbox->flow_model(dynamic,starting_set_bounding_box,time_step,flow_bounds);
    ARIADNE_LOG(3,"flow_model="<<flow_model<<"\n");
    ScalarTaylorFunction identity_time_expression=ScalarTaylorFunction::variable(BoxType(1u,Interval(-time_step,+time_step)),0u);
    flow_set_model=unchecked_apply(flow_model,combine(starting_set_model.models(),identity_time_expression.model()));
}



void ImageSetHybridEvolver::
compute_eventBlockingTimes_and_nonTransverseEvents(std::map<DiscreteEvent,TimeModelType>& event_blocking_times,
                        std::set<DiscreteEvent>& non_transverse_events,
                        const std::map<DiscreteEvent,VectorFunction>& guards,
                        const FlowSetModelType& flow_set_model) const
{
    ARIADNE_LOG(3,"Computing blocking events.\n");

    uint dimension=flow_set_model.result_size();
    const double SMALL_RELATIVE_TIME = 1./16;
    FlowSetModelType positive_flow_set_model(split(flow_set_model.models(),dimension,1));

    for(std::map<DiscreteEvent,VectorFunction>::const_iterator guard_iter=guards.begin();
        guard_iter!=guards.end(); ++guard_iter)
    {
        const DiscreteEvent event=guard_iter->first;
        const VectorFunction guard=guard_iter->second;
        tribool active = this->_toolbox->active(guard,positive_flow_set_model);
        if(possibly(active)) {
            ARIADNE_LOG(3,"Event "<<event<<" possibly active.\n");
            TimeModelType crossing_time_model;
            Interval normal_derivative;
            try {
                crossing_time_model=this->_toolbox->scaled_crossing_time(guard,flow_set_model);
                normal_derivative=this->normal_derivative(guard,flow_set_model,crossing_time_model);
                assert(normal_derivative.lower()>0 || normal_derivative.upper()<0);
                if(normal_derivative.lower()>0) {
                    ARIADNE_LOG(3,"Event "<<event<<" inserted into blocking times.\n");
                    event_blocking_times[event]=crossing_time_model;
                }
            }
            catch(DegenerateCrossingException e) {
                ARIADNE_LOG(3,"Degenerate Crossing exception catched.\n");
                BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
                Interval touching_time_interval=this->_toolbox->scaled_touching_time_interval(guard,flow_set_model);
                TimeModelType touching_time_model=this->_toolbox->time_model(touching_time_interval,space_domain);
                // Use 1.0 as upper bound above since flow set model has time interval normalised to [-1,+1]
                ARIADNE_LOG(3,"touching_time_interval="<<touching_time_interval<<"\n");
                if(touching_time_interval.upper()>=0 && touching_time_interval.lower()<=1.0) {
                    SetModelType finishing_set_model=partial_evaluate(flow_set_model.models(),dimension,1.0);
                    tribool finishing_set_active=this->_toolbox->active(guard,finishing_set_model);
                    if(definitely(finishing_set_active)) {
                        ARIADNE_LOG(3,"Event is definitely finally active, inserting it into blocking times.\n");
                        event_blocking_times[event]=touching_time_model;
                    } else if(possibly(finishing_set_active)) {
                        ARIADNE_LOG(3,"Event is possibly finally active.\n");
                        if(touching_time_interval.lower()>SMALL_RELATIVE_TIME) {
                            ARIADNE_LOG(3,"lower touching time is greater than zero, inserting event into blocking times.\n");
                            TaylorModel lower_touching_time_model=
                                this->_toolbox->time_model(touching_time_interval.lower(),space_domain);
                            event_blocking_times[finishing_event]=lower_touching_time_model;
                        } else {
                            ARIADNE_LOG(3,"DANGER: we can't determine whether the crossing is completely finished or not..\n");
                            // FIXME: Here we are stuck, we can't determine whether the crossing is completely finished or not.
                            // Just put this in as a blocking event and hope for the best...
                            // event_blocking_times[event]=touching_time_model;
                            // DAVIDE: setting this as a blocking event is not correct,
                            //         because it causes continuous evolution to stop.
                            //         I think it is better to put it into non-transverse-events.
                            non_transverse_events.insert(event);
                        }
                    } else {
                        ARIADNE_LOG(3,"After the flow step, the event is not active again, so the crossing was tangential.\n");
                        // After the flow step, the event is not active again, so the crossing was tangential
                        non_transverse_events.insert(event);
                    }
                }
            } // end non-transverse crossing
        } // end possibly active
    } // end main loop

    return;
}


void ImageSetHybridEvolver::
compute_blockingTime_and_relatedEvents(std::set<DiscreteEvent>& blocking_events,
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


void ImageSetHybridEvolver::
compute_activationEvents(std::map<DiscreteEvent,tuple<tribool,TimeModelType,tribool> >& activation_events,
                          const std::map<DiscreteEvent,VectorFunction>& activations, const FlowSetModelType& flow_set_model) const
{
    SetModelType initial_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,0.0);
    SetModelType final_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,1.0);
    for(std::map<DiscreteEvent,VectorFunction>::const_iterator iter=activations.begin(); iter!=activations.end(); ++iter) {
        DiscreteEvent event=iter->first;
        VectorFunction activation=iter->second;
        tribool active=this->_toolbox->active(activation,flow_set_model);
        if(possibly(active)) {
            tribool initially_active=this->_toolbox->active(activation,initial_set_model);
            tribool finally_active=this->_toolbox->active(activation,final_set_model);
            TimeModelType crossing_time_model=this->crossing_time(activation,flow_set_model);
            activation_events.insert(make_pair(event,make_tuple(initially_active,crossing_time_model,finally_active)));
        }
    }

}


void ImageSetHybridEvolver::
compute_activationTimes(std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >& activation_times,
                         const std::map<DiscreteEvent,VectorFunction>& activations,
                         const FlowSetModelType& flow_set_model,
                         const TimeModelType& blocking_time_model,
                         const Semantics semantics) const
{
    SetModelType initial_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,0.0);
    SetModelType final_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,1.0);
    TimeModelType zero_time_model=blocking_time_model*0.0;

    for(std::map<DiscreteEvent,VectorFunction>::const_iterator iter=activations.begin(); iter!=activations.end(); ++iter) {
        DiscreteEvent event=iter->first;
        VectorFunction activation=iter->second;

        ARIADNE_LOG(3,"Computing activation time for event "<<event<<"...");

        // Compute whether the event might be enabled on the entire time interval
        tribool active=this->_toolbox->active(activation,flow_set_model);

        if(definitely(active)) {
            // The event is enabled over the entire time interval
            ARIADNE_LOG(3,"event is enabled over the entire time interval.\n");
            activation_times.insert(make_pair(event,make_tuple(zero_time_model,blocking_time_model)));
        } else if(possibly(active)) {
            ARIADNE_LOG(3,"event is possibly enabled.\n");
            // Compute whether the event is enabled at the beginning and end of the time interval
            tribool initially_active=this->_toolbox->active(activation,initial_set_model);
            tribool finally_active=this->_toolbox->active(activation,final_set_model);

            TimeModelType crossing_time_model=this->crossing_time(activation,flow_set_model);

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


void
ImageSetHybridEvolver::
_logStepAtVerbosity1(const std::list<HybridTimedSetType>& working_sets,
					 const EnclosureListType& reach_sets,
					 const EventListType& initial_events,
					 const TimeModelType& initial_time_model,
					 const SetModelType& initial_set_model,
					 const DiscreteState& initial_location) const
{
    if(verbosity==1) {
        ARIADNE_LOG(1,"#w="<<std::setw(4)<<working_sets.size()
                    <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                    <<" s="<<std::setw(3)<<std::left<<initial_events.size()
                    <<" t="<<std::fixed<<initial_time_model.value()
                    <<" r="<<std::setw(7)<<initial_set_model.radius()
                    <<" l="<<std::setw(3)<<std::left<<initial_location
                    <<" c="<<initial_set_model.centre()
                    <<" e="<<initial_events
                    <<"\n");
    }
}

bool ImageSetHybridEvolver::
_is_reachableSet_outside_disproveBounds(const uint numDivisions,
			  const TaylorSet& reachable_set,
			  const Box& disprove_bounds) const
{
	if (definitely(reachable_set.inside(disprove_bounds)))
		return false;

	// Otherwise, if we are allowed to split further
	if (numDivisions>0) {
		// Split on the time variable (ASSUMED TO BE THE LAST)
		std::pair<TaylorSet,TaylorSet> split_sets = reachable_set.split(reachable_set.generators_size()-1);

		if (_is_reachableSet_outside_disproveBounds(numDivisions-1,split_sets.first,disprove_bounds))
			return true;
		if (_is_reachableSet_outside_disproveBounds(numDivisions-1,split_sets.second,disprove_bounds))
			return true;

		return false;
	}

	// If we cannot split further, we check for the set being outside (possibly not inside and definitely disjoint)
	return possibly(!reachable_set.inside(disprove_bounds)) && definitely(reachable_set.disjoint(disprove_bounds));
}

void ImageSetHybridEvolver::
_computeEvolutionForEvents(std::list< HybridTimedSetType >& working_sets,
						   EnclosureListType& intermediate_sets,
						   const SystemType& system,
						   const DiscreteState& location,
						   const std::set<DiscreteEvent>& blocking_events,
						   const EventListType& events,
						   const ActivationTimesType& activation_times,
						   const SetModelType& flow_set_model,
						   const TimeModelType& time_model,
						   const TimeModelType& blocking_time_model,
						   const Float& time_step) const
{
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
        } else if(event.is_transition()) { // not an invariant
            intermediate_sets.adjoin(make_pair(location,evolved_set_model));
            const DiscreteTransition& transition=system.transition(event,location);
            SetModelType jump_set_model=apply(transition.reset(),evolved_set_model);
            DiscreteState jump_location=transition.target();
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
        SetModelType jump_set_model=apply(transition.reset(),active_set_model);
        ARIADNE_LOG(3,"  jump_set_model="<<active_set_model.range()<<";\n");
        const TimeModelType active_time_model = this->_toolbox->reachability_time(time_model+lower_active_time_model*time_step,time_model+upper_active_time_model*time_step);
        ARIADNE_LOG(3,"  active_time_model="<<active_time_model.range()<<".\n");

        DiscreteState jump_location=transition.target();
        std::vector<DiscreteEvent> jump_events=events;
        jump_events.push_back(event);
        working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,active_time_model));
    }
}

void
ImageSetHybridEvolver::
_processInitiallyActiveBlockingEvents_continuous(EnclosureListType& reach_sets,
												   EnclosureListType& intermediate_sets,
												   const std::map<DiscreteEvent,VectorFunction>& invariants,
												   const SetModelType& set_model,
												   const DiscreteState& location) const
{
	bool has_active_invariant = false;

    for(std::map<DiscreteEvent,VectorFunction>::const_iterator inv_it = invariants.begin(); inv_it != invariants.end(); inv_it++) {
        tribool active=this->_toolbox->active(inv_it->second,set_model);
        if(definitely(active)) {
			has_active_invariant = true;
			break;
        }
    }

    if (!has_active_invariant) {
    	ARIADNE_LOG(2,"A forced transition is initially active: no continuous evolution, adjoining reach set.\n");
        reach_sets.insert(make_pair(location,set_model));
        intermediate_sets.insert(make_pair(location,set_model));
    }
    else
    	ARIADNE_LOG(2,"An invariant is initially active: no continuous evolution, discarding the set.\n");
}

void
ImageSetHybridEvolver::
_compute_blocking_info(std::set<DiscreteEvent>& non_transverse_events,
				  std::set<DiscreteEvent>& blocking_events,
				  TimeModelType& blocking_time_model,
				  const TimeModelType& time_step_model,
				  const SetModelType& flow_set_model,
				  const std::map<DiscreteEvent,VectorFunction>& guards,
				  double SMALL_RELATIVE_TIME) const
{
    // Compute event blocking times
    std::map<DiscreteEvent, TimeModelType> event_blocking_times;

    event_blocking_times[finishing_event]=time_step_model;
    compute_eventBlockingTimes_and_nonTransverseEvents(event_blocking_times,non_transverse_events,guards,flow_set_model);
    ARIADNE_LOG(2,"event_blocking_times="<<event_blocking_times<<"\n");

    std::map<DiscreteEvent,Interval> event_blocking_time_intervals;
    for(std::map<DiscreteEvent, TimeModelType>::const_iterator iter=event_blocking_times.begin();
        iter!=event_blocking_times.end(); ++iter) { event_blocking_time_intervals[iter->first]=iter->second.range(); }

    ARIADNE_LOG(2,"event_blocking_time_intervals="<<event_blocking_time_intervals<<"\n");
    ARIADNE_LOG(2,"non_transverse_events="<<non_transverse_events<<"\n\n");

    // Compute blocking events
    compute_blockingTime_and_relatedEvents(blocking_events,blocking_time_model,event_blocking_times);
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
    ARIADNE_LOG(2,"(adjusted)blocking_events="<<blocking_events<<"\n");
    ARIADNE_LOG(2,"(adjusted)blocking_time="<<blocking_time_model.range()<<"\n\n");
}


void
ImageSetHybridEvolver::
_compute_activation_info(std::map<DiscreteEvent,VectorFunction>& activations,
						 ActivationTimesType& activation_times,
						 const std::set<DiscreteEvent>& non_transverse_events,
						 const SetModelType& flow_set_model,
						 const TimeModelType& blocking_time_model,
						 const std::map<DiscreteEvent,VectorFunction>& guards,
						 const Semantics semantics) const
{
    // Treat non-transverse urgent events as non-urgent in upper semantics
    for(std::set<DiscreteEvent>::const_iterator iter=non_transverse_events.begin(); iter!=non_transverse_events.end(); ++iter) {
        if(iter->is_transition()) {     // If the event is a transition
            activations[*iter]=guards.find(*iter)->second;
       }
    }

    this->compute_activationTimes(activation_times,activations,flow_set_model,blocking_time_model,semantics);

    // Display activation time ranges
    std::map< DiscreteEvent,tuple<Interval> > activation_time_intervals;
    for(ActivationTimesType::const_iterator iter=activation_times.begin(); iter!=activation_times.end(); ++iter)
        activation_time_intervals.insert(make_pair(iter->first,make_tuple(iter->second.second.range())));
    ARIADNE_LOG(2,"activation_time_intervals="<<activation_time_intervals<<"\n\n");
}


void
ImageSetHybridEvolver::
_compute_and_adjoin_reachableSet(EnclosureListType& reach_sets,
								SetModelType& reachable_set,
								const DiscreteState& location,
							    const SetModelType& flow_set_model,
							    const TimeModelType& zero_time_model,
		 	 	 	 	 	 	const TimeModelType& blocking_time_model) const
{
    ARIADNE_LOG(4,"flow_set_model="<<flow_set_model<<"\n");
    ARIADNE_LOG(4,"zero_time_model="<<zero_time_model<<"\n");
    ARIADNE_LOG(4,"blocking_time_model="<<blocking_time_model<<"\n");
    reachable_set=this->_toolbox->reachability_step(flow_set_model,zero_time_model,blocking_time_model);
    reach_sets.adjoin(make_pair(location,reachable_set));

	ARIADNE_LOG(2,"reachable_set="<<reachable_set<<"\n");
    ARIADNE_LOG(2,"reachable_set.argument_size()="<<reachable_set.argument_size()<<"\n");
    ARIADNE_LOG(2,"reachable_set.range()="<<reachable_set.range()<<"\n");
}

void
ImageSetHybridEvolver::
_logEvolutionStepInitialState(const EventListType& events,
							  const TimeModelType& time_model,
							  const DiscreteState& location,
							  const SetModelType& set_model,
							  const VectorFunction& dynamic,
							  const std::map<DiscreteEvent,VectorFunction>& invariants,
							  const std::list<DiscreteTransition>& transitions,
							  const std::map<DiscreteEvent,VectorFunction>& guards,
							  const std::map<DiscreteEvent,VectorFunction>& activations) const
{
    ARIADNE_LOG(2,"events = "<<events<<" ");
    ARIADNE_LOG(2,"time_range = "<<time_model.range()<<" ");
    ARIADNE_LOG(2,"time_model generators = "<<time_model.argument_size()<<" ");
    ARIADNE_LOG(2,"location = "<<location<<" ");
    ARIADNE_LOG(2,"box = "<<set_model.range()<<" ");
    ARIADNE_LOG(2,"generators = "<<set_model.argument_size()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(set_model.range())<<"\n\n");

    ARIADNE_LOG(2,"dynamic = "<<dynamic<<"\n");
    ARIADNE_LOG(2,"invariants = "<<invariants<<"\n");
    ARIADNE_LOG(2,"transitions = "<<transitions<<"\n\n");
    ARIADNE_LOG(2,"guards = "<<guards<<"\n\n");
    ARIADNE_LOG(2,"activations = "<<activations<<"\n\n");
}

void
ImageSetHybridEvolver::
_processInitiallyActiveBlockingEvents(std::list< HybridTimedSetType >& working_sets,
									  EnclosureListType& reach_sets,
		   	   	   	   	   	   	   	  EnclosureListType& intermediate_sets,
		   	   	   	   	   	   	   	  const DiscreteState& location,
		   	   	   	   	   	   	   	  const EventListType& events,
		   	   	   	   	   	   	   	  const SetModelType& set_model,
		   	   	   	   	   	   	   	  const TimeModelType& time_model,
		   	   	   	   	   	   	   	  const std::map<DiscreteEvent,VectorFunction>& invariants,
		   	   	   	   	   	   	   	  const std::list<DiscreteTransition>& transitions) const
{
    // Check if at least one invariant is definitely active in this mode
	bool has_active_invariant = false;
    for(std::map<DiscreteEvent,VectorFunction>::const_iterator inv_it = invariants.begin(); inv_it != invariants.end(); inv_it++) {
        tribool active=this->_toolbox->active(inv_it->second,set_model);
        if(definitely(active)) {
			has_active_invariant = true;
			break;
        }
    }

	// The set model is added, and the transitions checked, only if there is no definitely active invariant
    if (!has_active_invariant) {

        reach_sets.insert(make_pair(location,set_model));
        intermediate_sets.insert(make_pair(location,set_model));

        ARIADNE_LOG(2,"A blocking event (not an invariant) is initially active: no continuous evolution, just apply discrete events.\n");

        // No continuous evolution; just apply discrete events
        for(std::list<DiscreteTransition>::const_iterator iter=transitions.begin();
            iter!=transitions.end(); ++iter)
        {
            ARIADNE_LOG(3,"Testing transition "<<iter->event()<<".\n");
            tribool active=this->_toolbox->active(iter->activation(),set_model);
            if(definitely(active)) {
                ARIADNE_LOG(3,"Transition "<<iter->event()<<" should be activated.\n");
                SetModelType jump_set_model=this->_toolbox->reset_step(iter->reset(),set_model);
                TimeModelType jump_time_model=time_model;
                DiscreteState jump_location=iter->target();
                std::vector<DiscreteEvent> jump_events=events;
                jump_events.push_back(iter->event());
                ARIADNE_LOG(2,"Pushing back "<<jump_set_model<<" into working_sets\n");
                working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,jump_time_model));
            }
        }
    }
    else
    {
    	ARIADNE_LOG(2,"An invariant is definitely initially active: discarding the set.\n");
    }

}


void
ImageSetHybridEvolver::
_evolution_add_initialSet(std::list< HybridTimedSetType >& working_sets,
								  const EnclosureType& initial_set) const
{
    ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
    DiscreteState initial_location;
    ContinuousEnclosureType initial_continuous_set;
    make_lpair(initial_location,initial_continuous_set)=initial_set;
    ARIADNE_LOG(6,"initial_location = "<<initial_location<<"\n");
    SetModelType initial_set_model=this->_toolbox->set_model(initial_continuous_set);

	// Check for non-zero maximum step size
	ARIADNE_ASSERT_MSG(this->_parameters->hybrid_maximum_step_size[initial_location] > 0, "Error: the maximum step size for location " << initial_location.name() << " is zero.");
	// Check for match between the enclosure cell size and the set size
	ARIADNE_ASSERT_MSG(this->_parameters->maximum_enclosure_cell.size() == initial_set_model.size(), "Error: mismatch between the maximum_enclosure_cell size and the set size.");

    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    TimeModelType initial_time_model=this->_toolbox->time_model(0.0,Box(initial_set_model.argument_size()));
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push_back(make_tuple(initial_location,EventListType(),initial_set_model,initial_time_model));
}

bool
ImageSetHybridEvolver::
_isEnclosureTooLarge(const SetModelType& initial_set_model) const
{
	const Vector<Interval> initial_set_model_range = initial_set_model.range();

	for (uint i=0;i<initial_set_model_range.size();++i)
		if (initial_set_model_range[i].width() > this->_parameters->maximum_enclosure_cell[i])
			return true;

	return false;
}


void
ImageSetHybridEvolver::
_set_model_reduction(HybridTimedSetType& timed_set) const
{
	if (this->_parameters->enable_set_model_reduction &&
		!((timed_set.second.size()+1) % (this->_parameters->set_model_events_size_interleaving+1)))
	{
		SetModelType new_set_model = this->_toolbox->set_model(ContinuousEnclosureType(timed_set.third.range()));
		TimeModelType new_time_model = this->_toolbox->time_model(timed_set.fourth.range(),Box(new_set_model.argument_size()));
		timed_set = make_tuple(timed_set.first,timed_set.second, new_set_model, new_time_model);
	}
}

void
ImageSetHybridEvolver::
_add_subdivisions(std::list< HybridTimedSetType >& working_sets,
				  const array< TimedSetModelType >& subdivisions,
				  const DiscreteState& initial_location,
				  const EventListType& initial_events,
				  const uint dimension) const
{
    ARIADNE_LOG(3,"subdivisions.size()="<<subdivisions.size()<<"\n");
    for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        ARIADNE_LOG(3,"subdivided_timed_set_model.range()="<<subdivided_timed_set_model.range()<<"\n");
        SetModelType subdivided_set_model=Vector<TaylorModel>(project(subdivided_timed_set_model.models(),range(0,dimension)));
        TimeModelType subdivided_time_model=subdivided_timed_set_model[dimension];
        ARIADNE_LOG(3,"subdivided_set_model.range()="<<subdivided_set_model.range()<<"\n");
        ARIADNE_LOG(3,"subdivided_set_model.radius()*10000="<<radius(subdivided_set_model.range())*10000<<"\n");
        ARIADNE_LOG(3,"subdivided_time_model.range()="<<subdivided_time_model.range()<<"\n");
        working_sets.push_back(make_tuple(initial_location,initial_events,subdivided_set_model,subdivided_time_model));
    }
}

void
ImageSetHybridEvolver::
_add_models_subdivisions_autoselect(std::list< HybridTimedSetType >& working_sets,
		  	  	  	  	  	  	  	const SetModelType& initial_set_model,
		  	  	  	  	  	  	  	const TimeModelType& initial_time_model,
		  	  	  	  	  	  	  	const DiscreteState& initial_location,
		  	  	  	  	  	  	  	const EventListType& initial_events) const
{
    uint nd=initial_set_model.dimension();
    SetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
    _add_subdivisions(working_sets,subdivisions,initial_location,initial_events,nd);
}


void
ImageSetHybridEvolver::
_add_models_subdivisions_time(std::list< HybridTimedSetType >& working_sets,
		  	  	  	  	  	  const SetModelType& initial_set_model,
		  	  	  	  	  	  const TimeModelType& initial_time_model,
		  	  	  	  	  	  const DiscreteState& initial_location,
		  	  	  	  	  	  const EventListType& initial_events) const
{
    uint nd=initial_set_model.dimension();
    SetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model,nd);
    _add_subdivisions(working_sets,subdivisions,initial_location,initial_events,nd);
}


}  // namespace Ariadne

