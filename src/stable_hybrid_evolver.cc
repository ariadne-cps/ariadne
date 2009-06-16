/***************************************************************************
 *            stable_hybrid_evolver.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins, Davide Bresolin
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
#include "function_interface.h"
#include "taylor_model.h"
#include "taylor_set.h"
#include "taylor_expression.h"
#include "taylor_function.h"
#include "orbit.h"
#include "taylor_calculus.h"
#include "evolution_parameters.h"

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "stable_hybrid_evolver.h"

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
 
static const int BLOCKING_EVENT = -2;
using boost::shared_ptr;

class DegenerateCrossingException { };


StableHybridEvolver::StableHybridEvolver()
    : _parameters(new EvolutionParametersType()),
      _toolbox(new TaylorCalculus())
{
}



StableHybridEvolver::StableHybridEvolver(const EvolutionParametersType& p)
    : _parameters(new EvolutionParametersType(p)),
      _toolbox(new TaylorCalculus())
{
}


Orbit<StableHybridEvolver::EnclosureType> 
StableHybridEvolver::
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



namespace {


enum PredicateKind { INVARIANT, ACTIVATION, GUARD, TIME, MIXED };
enum CrossingKind { TRANSVERSE, TOUCHING, NONE, UNKNOWN };


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
    case TRANSVERSE: return os << "TRANSVERSE"; 
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
            os << ", touching_times="<<data.touching_time_interval; }
    }
    return os <<" } "; 
}

inline bool operator<(const DetectionData& d1, const DetectionData& d2) {
    return d1.id < d2.id;
}

}

void
StableHybridEvolver::
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


    typedef tuple<DiscreteState, IntegerType, SetModelType, TimeModelType> HybridTimedSetType;

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
        TimeModelType initial_time_model
            =this->_toolbox->time_model(0.0, Vector<Interval>(initial_set_model.argument_size(),Interval(-1,+1)));
        ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
        TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
        ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
        working_sets.push_back(make_tuple(initial_location,0,initial_set_model,initial_time_model));
    }


    while(!working_sets.empty()) {
        HybridTimedSetType current_set=working_sets.back();
        working_sets.pop_back();
        DiscreteState initial_location=current_set.first;
        IntegerType initial_steps=current_set.second;
        SetModelType initial_set_model=current_set.third;
        TimeModelType initial_time_model=current_set.fourth;
        RealType initial_set_radius=radius(initial_set_model.range());
        if(initial_time_model.range().lower()>=maximum_time || initial_steps>=maximum_steps) {
            final_sets.adjoin(initial_location,this->_toolbox->enclosure(initial_set_model));
        } else if(semantics == UPPER_SEMANTICS && this->_parameters->enable_subdivisions
                  && (initial_set_radius>this->_parameters->maximum_enclosure_radius)) {
            // Subdivide
            uint nd=initial_set_model.dimension();
            TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
            array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
                SetModelType subdivided_set_model=Vector<TaylorModel>(project(subdivided_timed_set_model.models(),range(0,nd)));
                TimeModelType subdivided_time_model=subdivided_timed_set_model[nd];
                working_sets.push_back(make_tuple(initial_location,initial_steps,subdivided_set_model,subdivided_time_model));
            }
        } else if((semantics == LOWER_SEMANTICS || !this->_parameters->enable_subdivisions) && 
                  this->_parameters->enable_premature_termination && 
                  initial_set_radius>this->_parameters->maximum_enclosure_radius) {
            std::cerr << "\nWARNING: Terminating evolution at time " << initial_time_model.value()
                      << " and set " << initial_set_model.centre() << " due to maximum radius being exceeded.\n";
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
                        <<" s="<<std::setw(3)<<std::left<<initial_steps
                        <<" t="<<std::setw(7)<<std::fixed<<initial_time_model.value()
                        <<" r="<<std::setw(7)<<initial_set_model.radius()
                        <<" c="<<initial_set_model.centre()
                        <<"                        ");
        }

    }

}



uint number_of_nonzero_elements(const TaylorSet& ts) {
    uint nnz=0;
    for(uint i=0; i!=ts.dimension(); ++i) {
        nnz+=ts[i].expansion().size();
    }
    return nnz;
}



void
StableHybridEvolver::
_evolution_step(std::vector< HybridTimedSetType >& working_sets,
                EnclosureListType& final_sets, 
                EnclosureListType& reach_sets, 
                EnclosureListType& intermediate_sets, 
                const SystemType& system, 
                const HybridTimedSetType& current_set,
                const TimeType& maximum_hybrid_time, 
                Semantics semantics, 
                bool reach) const
{

    ARIADNE_LOG(2,"\nStableHybridEvolver::_evolution_step(...)\n");
  
    DiscreteState initial_location(0);
    IntegerType initial_steps;
    SetModelType initial_set_model;
    TimeModelType initial_time_model;
    ARIADNE_LOG(9,"working_set = "<<current_set<<"\n");
    make_ltuple(initial_location,initial_steps,initial_set_model,initial_time_model)=current_set;

    ARIADNE_LOG(4,"initial_steps = "<<initial_steps<<"\n");
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"");
    ARIADNE_LOG(4,"initial_location = "<<initial_location<<"\n");
    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
  
    ARIADNE_LOG(2,"Starting evolution step from:\n");
    ARIADNE_LOG(2,"steps = "<<initial_steps<<" ");
    ARIADNE_LOG(2,"time_range = "<<initial_time_model.range()<<"\n");
    ARIADNE_LOG(2,"location = "<<initial_location<<" ");
    ARIADNE_LOG(2,"box = "<<initial_set_model.range()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(initial_set_model.range())<<"\n");
    // ARIADNE_LOG(2,"nnz = "<<number_of_nonzero_elements(initial_set_model)<<" ");
    ARIADNE_LOG(3,"initial_time_accuracy = "<<initial_time_model.accuracy()<<" ");
    ARIADNE_LOG(3,"initial_set_accuracy = "<<initial_set_model[0].accuracy()<<"\n");
    //const uint nd=initial_set_model.result_size();
    //const uint ng=initial_set_model.argument_size();
  

    /////////////// Main Evolution ////////////////////////////////

    const DiscreteMode& initial_mode=system.mode(initial_location);
    const FunctionType* dynamic_ptr=&initial_mode.dynamic();
  
    ARIADNE_LOG(6,"mode="<<initial_mode<<"\n");
    std::map< DiscreteEvent, shared_ptr<const FunctionInterface> > invariants
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
    ARIADNE_LOG(4,"flow_model = "<<flow_model<<"\n");
  
    // Compute the integration time model
    TimeModelType final_time_model=initial_time_model+step_size;
    ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
    TimeModelType integration_time_model=final_time_model-initial_time_model;
    ARIADNE_LOG(4, "integration_time_model="<<integration_time_model<<"\n");
  
    // Compute the flow tube (reachable set) model and the final set for the entire time step
    SetModelType final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
    ARIADNE_LOG(4,"final_set_model = "<<final_set_model<<"\n");
    SetModelType reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);         
    ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
    ARIADNE_LOG(4,"Done computing continuous evolution for the first time.\n");

    // Test for blocking.

    // Test invariants
    bool blocking=false;
    Float lower_blocking_time=step_size;
    for(std::map< DiscreteEvent, shared_ptr<const FunctionInterface> >::const_iterator
        iter=invariants.begin(); iter!=invariants.end(); ++iter) 
    {
        shared_ptr<const FunctionInterface> guard_ptr = iter->second;
        if(possibly(this->_toolbox->active(*guard_ptr,flow_bounds))) {
            ARIADNE_LOG(2,"One invariant is possibly active.\n");
            if(semantics == UPPER_SEMANTICS) {
                // For upper semantics, test if at at least one blocking event 
                // (invariant or urgent transition) is definitely finally activate.
                if(definitely(this->_toolbox->active(*guard_ptr,final_set_model))) {
                    blocking=true;
                    break;              // Skip other invariants
                }
            } else {    // LOWER_SEMANTICS
                // For lower semantics, find a time t0 such that no blocking events are 
                // (possibly) activated before t0.
                if(possibly(this->_toolbox->active(*guard_ptr,initial_set_model))) {
                    // Invariant is initially active. No continuous evolution possible
                    blocking = true;
                    lower_blocking_time = 0.0;
                    break;              // Skip remaining invariants      
                }                
                if(possibly(this->_toolbox->active(*guard_ptr,reach_set_model))) {
                    blocking = true;
                    try {
                        // Compute crossing time
                        ConstraintModelType guard_model=this->_toolbox->predicate_model(*guard_ptr,flow_bounds);
                        TimeModelType crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);                    
                        lower_blocking_time = min(lower_blocking_time, crossing_time_model.range().lower());
                    } catch(DegenerateCrossingException) { 
                        ARIADNE_LOG(3," DegenerateCrossing detected.\n"); 
                        lower_blocking_time=0.0;
                    }
                } else {
                    if(possibly(this->_toolbox->active(*guard_ptr,final_set_model))) {
                        blocking=true;                
                    }
                }
            }
        }
    }
    
    ARIADNE_LOG(3,"Blocking after checking invariants = "<<blocking<<"\n");
    if(semantics == LOWER_SEMANTICS) 
        ARIADNE_LOG(3,"Blocking time after checking invariants = "<<lower_blocking_time<<"\n");
        
    // Test transitions
    for(std::set< DiscreteTransition >::const_iterator 
        iter=transitions.begin(); iter!=transitions.end(); ++iter)
    {
        if(iter->forced()) {  // Test only forced transitions
            shared_ptr<const FunctionInterface> guard_ptr = iter->activation_ptr();
            if(possibly(this->_toolbox->active(*guard_ptr,flow_bounds))) {
                ARIADNE_LOG(2,"Forced transition "<<*iter<<" is possibly active.\n");
                if(semantics == UPPER_SEMANTICS) {
                    // For upper semantics, test if at at least one blocking event 
                    // (invariant or urgent transition) is definitely finally activate.
                    if(definitely(this->_toolbox->active(*guard_ptr,final_set_model))) {
                        blocking=true;
                        break;              // Skip other invariants
                    }
                } else {    // LOWER_SEMANTICS
                    // For lower semantics, find a time t0 such that no blocking events are 
                    // (possibly) activated before t0.
                    if(possibly(this->_toolbox->active(*guard_ptr,initial_set_model))) {
                        // Forced transition is initially active. No continuous evolution possible
                        blocking = true;
                        lower_blocking_time = 0.0;
                        break;              // Skip remaining invariants      
                    }                
                    if(possibly(this->_toolbox->active(*guard_ptr,reach_set_model))) {
                        blocking = true;
                        try {
                            // Compute crossing time
                            ConstraintModelType guard_model=this->_toolbox->predicate_model(*guard_ptr,flow_bounds);
                            TimeModelType crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);                    
                            lower_blocking_time = min(lower_blocking_time, crossing_time_model.range().lower());
                        } catch(DegenerateCrossingException) { 
                            ARIADNE_LOG(3," DegenerateCrossing detected.\n"); 
                            lower_blocking_time=0.0;
                        }
                    } else {
                        if(possibly(this->_toolbox->active(*guard_ptr,final_set_model))) {
                            blocking=true;                
                        }
                    }
                }
            }
        }
    }

    ARIADNE_LOG(3,"Blocking after checking invariants and transitions = "<<blocking<<"\n");
    if(semantics == LOWER_SEMANTICS) 
        ARIADNE_LOG(3,"Blocking time after checking invariants and transitions = "<<lower_blocking_time<<"\n");

    // For lower semantics, if blocking_time is less than step_size, final set should be recomputed
    if(lower_blocking_time < step_size && semantics == LOWER_SEMANTICS) {
        ARIADNE_LOG(3,"Blocking time smaller than step size. Recomputing final set model.\n");
        TimeModelType blocking_time_model=this->_toolbox->time_model(lower_blocking_time,initial_time_model.domain());
        final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,blocking_time_model);
        ARIADNE_LOG(3,"final_set_model = "<<final_set_model<<"\n");
    }
    
    

    // Compute discrete evolution
    bool activated=false;   // this will became true if at least one discrete transition is activated
    if(semantics == UPPER_SEMANTICS) {        
        // Upper semantics: execute every transition that is possibly active.
        ARIADNE_LOG(2,"Computing discrete transitions.\n");
        for(std::set< DiscreteTransition >::const_iterator 
            iter=transitions.begin(); iter!=transitions.end(); ++iter)
        {
            shared_ptr<const FunctionInterface> guard_ptr = iter->activation_ptr();
            if(possibly(this->_toolbox->active(*guard_ptr,flow_bounds))) {
                ARIADNE_LOG(2," "<<*iter<<" is possibly active.\n");
                if(possibly(this->_toolbox->active(*guard_ptr,reach_set_model))) {
                    Float lower_active_time = 0.0;
                    Float upper_active_time = step_size;
                    ConstraintModelType guard_model=this->_toolbox->predicate_model(*guard_ptr,flow_bounds);
                    try {
                        // Compute crossing time
                        TimeModelType crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);
                        ARIADNE_LOG(3," Crossing time: "<<crossing_time_model.range()<<"\n");
                        Float lower_crossing_time = crossing_time_model.range().lower();
                        Float upper_crossing_time = crossing_time_model.range().upper();
                        if(lower_crossing_time < 0.0 || lower_crossing_time > step_size) {
                            lower_crossing_time = 0.0;
                        }
                        if(upper_crossing_time < 0.0 || upper_crossing_time > step_size) {
                            upper_crossing_time = step_size;
                        }                        
                        tribool initially_active = this->_toolbox->active(*guard_ptr,initial_set_model);
                        tribool finally_active = this->_toolbox->active(*guard_ptr,final_set_model);
                        if(!possibly(initially_active)) {
                            // Transition is NOT initally active
                            lower_active_time = lower_crossing_time;
                        }
                        if(iter->forced()) {
                            // Forced transition
                            upper_active_time = upper_crossing_time;
                        } else  {
                            // Unforced transition
                            if(!possibly(finally_active)) {
                                // Transition is NOT finally active
                                upper_active_time = upper_crossing_time;
                            }
                        }                        
                    } catch(DegenerateCrossingException) { 
                        ARIADNE_LOG(3," DegenerateCrossing detected.\n"); 
                        ARIADNE_LOG(3," Computing active time interval using bisections...");
                        Interval active_time_interval = this->_toolbox->touching_time_interval(guard_model,flow_model,initial_set_model);
                        lower_active_time=active_time_interval.lower();
                        upper_active_time=active_time_interval.upper();
                        ARIADNE_LOG(3,"done\n");
                    }
                    ARIADNE_LOG(2,"Activation time = ["<<lower_active_time<<","<<upper_active_time<<"]\n");
                    
                    // the activation time should be included in [0.0, step_size] for the transition to be executed
                    if(upper_active_time >= 0.0 && lower_active_time <= step_size) {
                        activated = true;     // at least one transition have been activated
                        if(lower_active_time < 0.0) lower_active_time = 0.0;
                        if(upper_active_time > step_size) upper_active_time = step_size;
                        // Compute activation time models
                        TimeModelType lower_active_time_model=this->_toolbox->time_model(lower_active_time,initial_time_model.domain());
                        TimeModelType upper_active_time_model=this->_toolbox->time_model(upper_active_time,initial_time_model.domain());
                        TimeModelType active_time_model=this->_toolbox->reachability_time(lower_active_time_model,upper_active_time_model);
                        TimeModelType jump_time_model=this->_toolbox->reachability_time(initial_time_model+lower_active_time_model,initial_time_model+upper_active_time_model);
                        ARIADNE_LOG(4,"jump_time_model="<<jump_time_model<<"\n");
                        
                        // Compute activation set model and jump set model
                        shared_ptr<const FunctionInterface> reset_ptr=iter->reset_ptr();
                        DiscreteState jump_location=iter->target().location();
                        DiscreteEvent jump_event=iter->event();

                        SetModelType active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,lower_active_time,upper_active_time);
                        ARIADNE_LOG(4,"initial_location="<<initial_location<<", active_set_model="<<active_set_model<<"\n");
                        SetModelType jump_set_model=this->_toolbox->reset_step(*reset_ptr,active_set_model);
                        ARIADNE_LOG(4,"jump_location="<<jump_location<<", jump_set_model="<<jump_set_model<<"\n");
                        
                        // Adjoin intermediate and jump sets to the result
                        ARIADNE_LOG(3,"Adding jump set to the working sets.\n");                                
                        working_sets.push_back(make_tuple(jump_location,initial_steps+1,jump_set_model,jump_time_model));
                        intermediate_sets.adjoin(EnclosureType(initial_location,active_set_model));
                        intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
                    } else {
                        ARIADNE_LOG(3," Transition skipped!\n");
                    }
                }
            }
        }
    } else { // LOWER_SEMANTICS
        // Lower semantics: execute every transition that is definitely active.        
        ARIADNE_LOG(2,"Computing discrete transitions.\n");
        for(std::set< DiscreteTransition >::const_iterator 
            iter=transitions.begin(); iter!=transitions.end(); ++iter)
        {
            shared_ptr<const FunctionInterface> guard_ptr = iter->activation_ptr();
            shared_ptr<const FunctionInterface> reset_ptr=iter->reset_ptr();
            DiscreteState jump_location=iter->target().location();
            if(lower_blocking_time <= 0.0) {
                // No continuous evolution possible, execute only transitions that are initially active
                if(definitely(this->_toolbox->active(*guard_ptr,initial_set_model))) {
                    ARIADNE_LOG(3," "<<*iter<<" is initially active. Adding jump set to the working sets.\n");                                
                    activated = true;     // at least one transition have been activated
                    SetModelType jump_set_model=this->_toolbox->reset_step(*reset_ptr,initial_set_model);
                    // Adjoin intermediate and jump sets to the result
                    working_sets.push_back(make_tuple(jump_location,initial_steps+1,jump_set_model,initial_time_model));
                    intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
                }
            } else {    // lower_blocking_time > 0.0
                if(possibly(this->_toolbox->active(*guard_ptr,flow_bounds))) {
                    ARIADNE_LOG(2," "<<*iter<<" is possibly active.\n");
                    if(possibly(this->_toolbox->active(*guard_ptr,reach_set_model))) {
                        Float lower_active_time = lower_blocking_time;
                        Float upper_active_time = 0.0;
                        if(definitely(this->_toolbox->active(*guard_ptr,initial_set_model))) {
                            // Transition is definitely initially active
                            lower_active_time = 0.0;
                        }
                        if(definitely(this->_toolbox->active(*guard_ptr,final_set_model))) {
                            // Transition is definitely finally active
                            upper_active_time = lower_blocking_time;
                        }
                        if(lower_active_time > upper_active_time) {
                            // transition is not both finally and initially active, check crossing time
                            try {
                                // Compute crossing time
                                ConstraintModelType guard_model=this->_toolbox->predicate_model(*guard_ptr,flow_bounds);
                                TimeModelType crossing_time_model=this->_toolbox->crossing_time(guard_model,flow_model,initial_set_model);
                                ARIADNE_LOG(3," Crossing time: "<<crossing_time_model.range()<<"\n");
                                if(crossing_time_model.range().lower() >= 0.0) {
                                    lower_active_time = min(lower_active_time, crossing_time_model.range().lower());
                                }
                                 if(crossing_time_model.range().upper() >= 0.0) {
                                    upper_active_time = max(upper_active_time, crossing_time_model.range().upper());
                                }
                            } catch(DegenerateCrossingException) { 
                                ARIADNE_LOG(3," DegenerateCrossing detected.\n"); 
                            }
                        }
                        ARIADNE_LOG(2,"Activation time = ["<<lower_active_time<<","<<upper_active_time<<"]\n");
                    
                        // lower active time should be less or equal to upper active time for the transition to be executed
                        if(lower_active_time <= upper_active_time) {
                            activated = true;     // at least one transition have been activated
                            // Compute activation time models
                            TimeModelType lower_active_time_model=this->_toolbox->time_model(lower_active_time,initial_time_model.domain());
                            ARIADNE_LOG(4,"lower_active_time_model="<<lower_active_time_model<<"\n");                                           
                            TimeModelType upper_active_time_model=this->_toolbox->time_model(upper_active_time,initial_time_model.domain());
                            ARIADNE_LOG(4,"upper_active_time_model="<<upper_active_time_model<<"\n");
                            TimeModelType jump_time_model=this->_toolbox->reachability_time(initial_time_model+lower_active_time_model,initial_time_model+upper_active_time_model);
                            ARIADNE_LOG(4,"jump_time_model="<<jump_time_model<<"\n");
                            
                            // Compute activation set model and jump set model
                            shared_ptr<const FunctionInterface> reset_ptr=iter->reset_ptr();
                            DiscreteState jump_location=iter->target().location();
                            SetModelType active_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,lower_active_time,upper_active_time);
                            ARIADNE_LOG(4,"initial_location="<<initial_location<<", active_set_model="<<active_set_model<<"\n");
                            SetModelType jump_set_model=this->_toolbox->reset_step(*reset_ptr,active_set_model);
                            ARIADNE_LOG(4,"jump_location="<<jump_location<<", jump_set_model="<<jump_set_model<<"\n");
                        
                            // Adjoin intermediate and jump sets to the result
                            ARIADNE_LOG(3,"Adding jump set to the working sets.\n");                                
                            working_sets.push_back(make_tuple(jump_location,initial_steps+1,jump_set_model,jump_time_model));
                            intermediate_sets.adjoin(EnclosureType(initial_location,active_set_model));
                            intermediate_sets.adjoin(EnclosureType(jump_location,jump_set_model));
                        } else {
                            ARIADNE_LOG(3," Transition skipped!\n");
                        }
                    }
                }
            }
        }
    }

    // For lower semantics, if lower_blocking_time is smaller than step_size, refine continuous evolution
    // WARNING: this step is executed AFTER the activation of the transitions to improve transition detection.
    if(lower_blocking_time < step_size && lower_blocking_time > 0.0 && semantics == LOWER_SEMANTICS) {
        // recompute continuous evolution
        ARIADNE_LOG(2,"Blocking time smaller than step size. Recomputing continuous evolution.\n");
        // Compute the flow model
        flow_model=this->_toolbox->flow_model(*dynamic_ptr,initial_set_bounds,lower_blocking_time,flow_bounds);
        ARIADNE_LOG(4,"flow_model = "<<flow_model<<"\n"); 
        // Compute the integration time model
        final_time_model=initial_time_model+lower_blocking_time;
        ARIADNE_LOG(6,"final_time_model = "<<final_time_model<<"\n");
        integration_time_model=this->_toolbox->time_model(lower_blocking_time,initial_time_model.domain());
        ARIADNE_LOG(3, "integration_time_model="<<integration_time_model<<"\n");
        final_set_model=this->_toolbox->integration_step(flow_model,initial_set_model,integration_time_model);
        ARIADNE_LOG(4,"final_set_model = "<<final_set_model<<"\n");
        reach_set_model=this->_toolbox->reachability_step(flow_model,initial_set_model,zero_time,integration_time_model);         
        ARIADNE_LOG(6,"reach_set_model = "<<reach_set_model<<"\n");
        ARIADNE_LOG(4,"Done computing continuous evolution for the second time.\n");
    }    

    // Adjoin intermediate and reach set to the results
    reach_sets.adjoin(EnclosureType(initial_location,reach_set_model));
    intermediate_sets.adjoin(EnclosureType(initial_location,final_set_model));
    if(blocking && !activated) {
        // Blocking: no continuous evolution allowed, no discrete transition activated
        ARIADNE_LOG(2,"Blocking dected, evolution stopped.\n");
        final_sets.adjoin(EnclosureType(initial_location,final_set_model));
    } 
    if(!blocking) {
        // Continuous evolution should go on
        ARIADNE_LOG(4,"initial_location="<<initial_location<<", initial_steps="<<initial_steps<<"\n");
        ARIADNE_LOG(4,"final_set_model="<<final_set_model<<"\n");
        ARIADNE_LOG(4,"final_time_model="<<final_time_model<<"\n");
        ARIADNE_LOG(3,"Adding final set to the working sets.\n");        
        working_sets.push_back(make_tuple(initial_location,initial_steps,final_set_model,final_time_model));
    }
  
    ARIADNE_LOG(2,"Done evolution_step.\n\n");

}

StableHybridEvolver::TimedEnclosureListType
StableHybridEvolver::
timed_evolution(const SystemType& system, 
                const EnclosureType& initial_set, 
                const TimeType& maximum_hybrid_time, 
                Semantics semantics, 
                bool reach) const
{
    verbosity=0;
  
    typedef boost::shared_ptr< const FunctionInterface > FunctionConstPointer;

    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");

    // dummy list sets for _evolution_step
    EnclosureListType final_sets;
    EnclosureListType reach_sets;
    EnclosureListType intermediate_sets;
    
    // result set
    TimedEnclosureListType result;

    const IntegerType maximum_steps=maximum_hybrid_time.discrete_time;
    const Float maximum_time=maximum_hybrid_time.continuous_time;



    typedef tuple<DiscreteState, IntegerType, SetModelType, TimeModelType> HybridTimedSetType;

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
        TimeModelType initial_time_model
            =this->_toolbox->time_model(0.0, Vector<Interval>(initial_set_model.argument_size(),Interval(-1,+1)));
        ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
        TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
        ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
        working_sets.push_back(make_tuple(initial_location,0,initial_set_model,initial_time_model));
    }


    while(!working_sets.empty()) {
        HybridTimedSetType current_set=working_sets.back();
        working_sets.pop_back();
        DiscreteState initial_location=current_set.first;
        IntegerType initial_steps=current_set.second;
        SetModelType initial_set_model=current_set.third;
        TimeModelType initial_time_model=current_set.fourth;
        RealType initial_set_radius=radius(initial_set_model.range());
        ARIADNE_LOG(5,"initial_steps = "<<initial_steps<<"\n");
        ARIADNE_LOG(5,"initial_time_model.range = "<<initial_time_model.range() <<"\n");
        ARIADNE_LOG(5,"initial_location = "<<initial_location<<"\n");
        ARIADNE_LOG(5,"initial_set_model.range = "<<initial_set_model.range()<<"\n");

        if(initial_time_model.range().lower()>=maximum_time || initial_steps>=maximum_steps) {
            Interval final_time(initial_time_model.range());
            EnclosureType final_enclosure(initial_location,this->_toolbox->enclosure(initial_set_model));
            result.push_back(TimedEnclosureType(final_time,final_enclosure));
        } else if(semantics == UPPER_SEMANTICS && this->_parameters->enable_subdivisions
                  && (initial_set_radius>this->_parameters->maximum_enclosure_radius)) {
            // Subdivide
            uint nd=initial_set_model.dimension();
            TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
            array< TimedSetModelType > subdivisions=this->_toolbox->subdivide(initial_timed_set_model);
            for(uint i=0; i!=subdivisions.size(); ++i) {
                TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
                SetModelType subdivided_set_model=Vector<TaylorModel>(project(subdivided_timed_set_model.models(),range(0,nd)));
                TimeModelType subdivided_time_model=subdivided_timed_set_model[nd];
                working_sets.push_back(make_tuple(initial_location,initial_steps,subdivided_set_model,subdivided_time_model));
            }
        } else if((semantics == LOWER_SEMANTICS || !this->_parameters->enable_subdivisions) && 
                  this->_parameters->enable_premature_termination && 
                  initial_set_radius>this->_parameters->maximum_enclosure_radius) {
            Interval final_time(initial_time_model.range());
            EnclosureType final_enclosure(initial_location,this->_toolbox->enclosure(initial_set_model));
            result.push_back(TimedEnclosureType(final_time,final_enclosure));
            std::cerr << "WARNING: Terminating evolution at time " << initial_time_model
                      << " and set " << initial_set_model << " due to maximum radius being exceeded.";
        } else {
            // insert current working set into result
            Interval current_time(initial_time_model.range());
            EnclosureType current_enclosure(initial_location,this->_toolbox->enclosure(initial_set_model));
            result.push_back(TimedEnclosureType(current_time,current_enclosure));            
            // Compute evolution
            this->_evolution_step(working_sets,
                                  final_sets,reach_sets,intermediate_sets,
                                  system,current_set,maximum_hybrid_time,
                                  semantics,reach);
        }
    }

    return result;
}




}  // namespace Ariadne

