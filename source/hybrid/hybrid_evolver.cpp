/***************************************************************************
 *            hybrid/hybrid_evolver.cpp
 *
 *  Copyright  2009-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <typeinfo>
#include "config.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/polynomial.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"
#include "geometry/grid_paving.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_automaton_interface.hpp"
#include "hybrid/hybrid_evolver.hpp"
#include "dynamics/orbit.hpp"

#include "solvers/integrator.hpp"
#include "solvers/solver.hpp"

#include "hybrid/hybrid_evolver.hpp"

namespace {

} // namespace

namespace Ariadne {

static const FloatDPValue zero={0,dp};

inline auto operator+(Int n, FloatDPBounds x) -> decltype(FloatDPValue(n,dp)+x) { return FloatDPValue(n,dp)+x; }
inline auto operator-(Int n, FloatDPBounds x) -> decltype(FloatDPValue(n,dp)-x) { return FloatDPValue(n,dp)-x; }
inline auto operator/(FloatDPValue x, Nat n) -> decltype(x/FloatDPValue(n,dp)) { return x/FloatDPValue(n,dp); }

template<class PR> inline auto operator*(FloatUpperBound<PR> x1, PositiveFloatValue<PR> x2) -> FloatUpperBound<PR> {
    return FloatUpperBound<PR>(mul(up,x1._u,x2._v)); }
template<class PR> inline auto operator*(PositiveFloatUpperBound<PR> x1, PositiveFloatValue<PR> x2) -> PositiveFloatUpperBound<PR> {
    return PositiveFloatUpperBound<PR>(mul(up,x1._u,x2._v)); }

inline auto operator> (FloatDPLowerBound x, Real r) -> decltype(x> declval<FloatDPBounds>()) { return x> FloatDPBounds(r,dp); }
inline auto operator>=(FloatDPLowerBound x, Real r) -> decltype(x>=declval<FloatDPBounds>()) { return x>=FloatDPBounds(r,dp); }
inline auto operator> (FloatDPUpperBound x, Real r) -> decltype(x> declval<FloatDPBounds>()) { return x> FloatDPBounds(r,dp); }
inline auto operator<=(FloatDPUpperBound x, Real r) -> decltype(x<=declval<FloatDPBounds>()) { return x<=FloatDPBounds(r,dp); }

static const DiscreteEvent final_event("_tmax_");
static const DiscreteEvent step_event("_h_");

typedef Vector<FloatDPValue> ExactFloatVector;
typedef Vector<ExactIntervalType> ExactIntervalVectorType;

Set<DiscreteEvent> blocking_events(const Map<DiscreteEvent,TransitionData>& transitions);
Set<DiscreteEvent> activating_events(const Map<DiscreteEvent,TransitionData>& transitions);

OutputStream& operator<<(OutputStream& os, const CrossingKind& crk) {
    switch(crk) {
        case CrossingKind::DEGENERATE: os<<"DEGENERATE"; break;
        case CrossingKind::POSITIVE: os<<"POSITIVE"; break;
        case CrossingKind::NEGATIVE: os<<"NEGATIVE"; break;
        case CrossingKind::INCREASING: os<<"INCREASING"; break;
        case CrossingKind::DECREASING: os<<"DECREASING"; break;
        case CrossingKind::CONVEX: os<<"CONVEX"; break;
        case CrossingKind::CONCAVE: os<<"CONCAVE"; break;
        case CrossingKind::TRANSVERSE: os<<"TRANSVERSE"; break;
        case CrossingKind::GRAZING: os<<"GRAZING"; break;
        default: os << "unknown"; break;
    } return os;
}

OutputStream& operator<<(OutputStream& os, const StepKind& stpk) {
    switch(stpk) {
        case StepKind::CONSTANT_EVOLUTION_TIME: os<<"CONSTANT_EVOLUTION_TIME"; break;
        case StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME: os<<"SPACETIME_DEPENDENT_EVOLUTION_TIME"; break;
        case StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME: os<<"PARAMETER_DEPENDENT_EVOLUTION_TIME"; break;
        case StepKind::PARAMETER_DEPENDENT_FINISHING_TIME: os<<"PARAMETER_DEPENDENT_FINISHING_TIME"; break;
        case StepKind::SPACETIME_DEPENDENT_FINISHING_TIME: os<<"SPACETIME_DEPENDENT_FINISHING_TIME"; break;
        case StepKind::CONSTANT_FINISHING_TIME: os<<"CONSTANT_FINISHING_TIME"; break;
        default: os << "UNKNOWN"; break;
    } return os;
}

OutputStream& operator<<(OutputStream& os, const FinishingKind& fnshk) {
    switch(fnshk) {
        case FinishingKind::BEFORE_FINAL_TIME: os<<"BEFORE_FINAL_TIME"; break;
        case FinishingKind::AT_FINAL_TIME: os<<"AT_FINAL_TIME"; break;
        case FinishingKind::AFTER_FINAL_TIME: os<<"AFTER_FINAL_TIME"; break;
        case FinishingKind::STRADDLE_FINAL_TIME: os<<"STRADDLE_FINAL_TIME"; break;
        default: os << "UNKNOWN"; break;
    } return os;
}

OutputStream& operator<<(OutputStream& os, const TransitionData& transition) {
    os << "kind="<<transition.event_kind<<", guard="<<transition.guard_function<<", "
                 "target="<<transition.target;
    if(transition.event_kind!=EventKind::INVARIANT) { os<<", reset="<<transition.reset_function; }
    return os;
}

OutputStream& operator<<(OutputStream& os, const TimingData& timing) {
    os << "step_kind="<<timing.step_kind<<", finishing_kind="<<timing.finishing_kind<<", step_size="<<timing.step_size<<", "
       << "final_time="<<timing.final_time;
    if(timing.step_kind==StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME) {
        os <<", spacial_evolution_time="<<timing.spacetime_dependent_evolution_time;
    } else if(timing.step_kind==StepKind::PARAMETER_DEPENDENT_FINISHING_TIME) {
        os << ", parameter_dependent_finishing_time="<<timing.parameter_dependent_finishing_time;
    }
    os <<", parameter_dependent_evolution_time="<<timing.parameter_dependent_evolution_time;
    return os;
}

OutputStream& operator<<(OutputStream& os, const CrossingData& crossing_data) {
    os << "{kind="<<crossing_data.crossing_kind;
    if(crossing_data.crossing_kind==CrossingKind::TRANSVERSE) {
        os << ", crossing_time="<<crossing_data.crossing_time;
    }
    if(crossing_data.crossing_kind==CrossingKind::GRAZING) {
        os << ", critical_time="<<crossing_data.critical_time;
    }
    os << "}";
    return os;
}

template<class T> struct LogOutput { T const& _ref; operator T const& () const { return _ref; } };
template<class T> LogOutput<T> log_output(T const& t) { return LogOutput<T>{t}; }

template<class K, class V> OutputStream& operator<<(OutputStream& os, LogOutput<Map<K,V>> const& map_output) {
    Map<K,V>const& map=map_output;
    os << "{"; for(auto iter=map.begin(); iter!=map.end(); ++iter) { os << (iter==map.begin()?"":",") << iter->first << ":" << log_output(iter->second); } return os << "}"; }
template<class T> OutputStream& operator<<(OutputStream& os, LogOutput<List<T>> const& lst) {
    os << "["; for(SizeType i=0; i!=lst._ref.size(); ++i) { os<<(i==0?"":",")<<LogOutput<T>{lst._ref[i]}; } return os << "]"; }
inline OutputStream& operator<<(OutputStream& os, LogOutput<CrossingData> const& crossing_data_output) {
    CrossingData const& crossing_data = crossing_data_output._ref;
    os << "{kind="<<crossing_data.crossing_kind;
    if(crossing_data.crossing_kind==CrossingKind::TRANSVERSE) {
        os<<",crossing_time_range="<<crossing_data.crossing_time.range();
        os<<",crossing_time_error="<<crossing_data.crossing_time.error();
    } else if(crossing_data.crossing_kind==CrossingKind::GRAZING) {
        os<<",critical_time_range="<<crossing_data.critical_time.range();
        os<<",critical_time_error="<<crossing_data.critical_time.error();
        os<<",guard_range_at_critical_time="<<crossing_data.guard_range_at_critical_time;
        os<<",evolve_bounds_at_critical_time="<<crossing_data.evolve_bounds_at_critical_time;
    }

    return os;
}

// Test if an event 'blocks' continuous evolution.
Bool is_blocking(EventKind evk) {
    switch(evk) {
        case EventKind::INVARIANT: case EventKind::PROGRESS: case EventKind::URGENT: case EventKind::IMPACT:
            return true;
        case EventKind::PERMISSIVE:
            return false;
        default:
            ARIADNE_FAIL_MSG("EventKind "<<evk<<" not recognised by is_blocking(...) predicate.");
    }
}

// Test if an event 'activates' a discrete transition.
Bool is_activating(EventKind evk) {
    switch(evk) {
        case EventKind::PERMISSIVE: case EventKind::URGENT: case EventKind::IMPACT:
            return true;
        case EventKind::INVARIANT: case EventKind::PROGRESS:
            return false;
        default:
            ARIADNE_FAIL_MSG("EventKind "<<evk<<" not recognised by is_activating(...) predicate.");
    }
}

// Extract the blocking events.
Set<DiscreteEvent> blocking_events(const Map<DiscreteEvent,TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for(Map<DiscreteEvent,TransitionData>::ConstIterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        if(is_blocking(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}

// Extract the activating events.
Set<DiscreteEvent> activating_events(const Map<DiscreteEvent,TransitionData>& transitions) {
    Set<DiscreteEvent> events;
    for(Map<DiscreteEvent,TransitionData>::ConstIterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        if(is_activating(transition_iter->second.event_kind)) {
            events.insert(transition_iter->first);
        }
    }
    return events;
}




HybridEvolverBase::EnclosureListType
HybridEvolverBase::evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const
{
    EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,false);
    return final;
}

HybridEvolverBase::EnclosureListType
HybridEvolverBase::reach(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const
{
    EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,true);
    return reachable;
}

//! \brief Compute an approximation to the evolution set under the given semantics.
Pair<HybridEvolverBase::EnclosureListType,HybridEvolverBase::EnclosureListType>
HybridEvolverBase::reach_evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const
{
    EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
    this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,true);
    return make_pair(reachable,final);
}

Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridExactBoxType& initial_box,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(2,"factory="<<this->function_factory());
    ARIADNE_LOG_PRINTLN("initial_box="<<initial_box);
    HybridEnclosure initial_enclosure(initial_box,EnclosureConfiguration(this->function_factory()));
    ARIADNE_LOG_PRINTLN_AT(1,"initial_enclosure="<<initial_enclosure);
    return this->orbit(initial_enclosure,termination,semantics);
}

Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridBoxSet& initial_box,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(2,"factory="<<this->function_factory());
    ARIADNE_LOG_PRINTLN("initial_box="<<initial_box);
    HybridEnclosure initial_enclosure(initial_box,this->system().continuous_state_space(initial_box.location()),EnclosureConfiguration(this->function_factory()));
    ARIADNE_LOG_PRINTLN_AT(1,"initial_enclosure="<<initial_enclosure);
    return this->orbit(initial_enclosure,termination,semantics);
}

Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridBoundedConstraintSet& initial_set,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("initial_set="<<initial_set);
    HybridEnclosure initial_enclosure(initial_set,this->system().continuous_state_space(initial_set.location()),EnclosureConfiguration(this->function_factory()));
    ARIADNE_LOG_PRINTLN_AT(1,"initial_enclosure="<<initial_enclosure);
    return this->orbit(initial_enclosure,termination,semantics);
}


Orbit<HybridEnclosure>
HybridEvolverBase::
orbit(const HybridEnclosure& initial,
      const HybridTerminationCriterion& termination,
      Semantics semantics) const
{
    DiscreteLocation initial_location=initial.location();
    RealSpace auxiliary_space=this->system().continuous_auxiliary_space(initial_location);
    EffectiveVectorMultivariateFunction auxiliary_function=this->system().auxiliary_function(initial_location);
    const_cast<HybridEnclosure&>(initial).set_auxiliary(auxiliary_space.variables(),auxiliary_function);

    EvolutionData evolution_data;
    evolution_data.semantics=semantics;
    evolution_data.initial_sets.push_back(HybridEnclosure(initial));
    while(!evolution_data.initial_sets.empty()) {
        this->_evolution_in_mode(evolution_data,termination);
    }
    ARIADNE_ASSERT(evolution_data.initial_sets.empty());
    ARIADNE_ASSERT(evolution_data.working_sets.empty());

    Orbit<HybridEnclosure> orbit(initial);
    orbit.adjoin_intermediate(ListSet<HybridEnclosure>(evolution_data.intermediate_sets));
    orbit.adjoin_reach(evolution_data.reach_sets);
    orbit.adjoin_final(evolution_data.final_sets);
    return orbit;
}

HybridEvolverBase::HybridEvolverBase(const SystemType& system)
{
    this->_create(system,make_taylor_function_factory());
}

HybridEvolverBase::HybridEvolverBase(const SystemType& system,
                                     const FunctionFactoryType& factory)
{
    this->_create(system,factory.clone());
}

Void
HybridEvolverBase::_create(
        const SystemType& system,
        FunctionFactoryType* factory)
{
    this->_sys_ptr=std::shared_ptr<SystemType>(system.clone());
    this->_function_factory_ptr=std::shared_ptr<FunctionFactoryType>(factory);
    this->_solver_ptr=std::shared_ptr<SolverInterface>(new IntervalNewtonSolver(1e-8,12));
    this->ALLOW_CREEP=true;
    this->ALLOW_UNWIND=false;
    //this->_configuration_ptr=std::shared_ptr<ConfigurationType>(new ConfigurationType());
}


const HybridEvolverBase::SystemType&
HybridEvolverBase::system() const
{
    return *this->_sys_ptr;
}

HybridEvolverBase::ConfigurationType&
HybridEvolverBase::configuration()
{
    return *this->_configuration_ptr;
}

const HybridEvolverBase::ConfigurationType&
HybridEvolverBase::configuration() const
{
    return *this->_configuration_ptr;
}

Void
HybridEvolverBase::set_function_factory(const FunctionFactoryType& factory)
{
    this->_function_factory_ptr=std::shared_ptr<FunctionFactoryType>(factory.clone());
}

const HybridEvolverBase::FunctionFactoryType&
HybridEvolverBase::function_factory() const
{
    return*this->_function_factory_ptr;
}

Void
HybridEvolverBase::set_integrator(const IntegratorInterface& integrator)
{
    this->_integrator_ptr=std::shared_ptr<IntegratorInterface>(integrator.clone());
}

Void
HybridEvolverBase::set_solver(const SolverInterface& solver)
{
    this->_solver_ptr=std::shared_ptr<SolverInterface>(solver.clone());
}


HybridEvolverBase::EnclosureType
HybridEvolverBase::enclosure(const HybridExactBox& initial_box) const
{
    return HybridEnclosure(initial_box,EnclosureConfiguration(this->function_factory()));
}

HybridEvolverBase::EnclosureType
HybridEvolverBase::enclosure(const HybridBoundedConstraintSet& initial_set) const
{
    return HybridEnclosure(initial_set,this->system().continuous_state_space(initial_set.location()),EnclosureConfiguration(this->function_factory()));
}


Void
HybridEvolverBase::
_evolution(ListSet<HybridEnclosure>& final,
           ListSet<HybridEnclosure>& reachable,
           ListSet<HybridEnclosure>& intermediate,
           HybridEnclosure const& initial_set,
           HybridTerminationCriterion const& termination_criterion,
           Semantics semantics,
           Bool reach) const
{
    EvolutionData evolution_data;
    evolution_data.semantics=semantics;
    evolution_data.initial_sets.push_back(HybridEnclosure(initial_set));

    while(!evolution_data.initial_sets.empty()) {
        this->_evolution_in_mode(evolution_data,termination_criterion);
    }

    final=evolution_data.final_sets;
    reachable=evolution_data.reach_sets;
    intermediate=evolution_data.intermediate_sets;
}

struct EvolutionStepData {
    EvolutionStepData() : progress(false), finishing(false), events() { }
    Bool progress;
    Bool finishing;
    Set<DiscreteEvent> events;
};

Void
HybridEvolverBase::
_log_summary(const EvolutionData& evolution_data, HybridEnclosure const& starting_set) const
{
    UpperBoxType starting_bounding_box=starting_set.state_bounding_box();
    UpperIntervalType starting_time_range=starting_set.time_range();
    UpperIntervalType starting_dwell_time_range=starting_set.dwell_time_range();
    Int old_precision = std::clog.precision();
    ARIADNE_LOG_PRINTLN(
              "#i="<<std::setw(4)<<std::left<<evolution_data.initial_sets.size()+1u
            <<"#w="<<std::setw(4)<<std::left<<evolution_data.working_sets.size()+1u
            <<"#r="<<std::setw(5)<<std::left<<evolution_data.reach_sets.size()
            <<"#f="<<std::setw(4)<<std::left<<evolution_data.final_sets.size()
            <<"#e="<<std::setw(3)<<std::left<<starting_set.previous_events().size()
            <<" #p="<<std::setw(2)<<std::left<<starting_set.number_of_parameters()
            <<" #c="<<std::setw(1)<<std::left<<starting_set.number_of_constraints()
            <<" t=["<<std::setw(6)<<std::setprecision(3)<<std::left<<std::fixed<<starting_time_range.lower_bound()
            <<","<<std::setw(6)<<std::left<<std::fixed<<starting_time_range.upper_bound()<<"]"<<std::flush
            <<" dwt=["<<std::setw(6)<<std::setprecision(3)<<std::left<<std::fixed<<starting_dwell_time_range.lower_bound()
            <<","<<std::setw(6)<<std::left<<std::fixed<<starting_dwell_time_range.upper_bound()<<"]"<<std::flush
            <<" c="<<starting_bounding_box.centre()
            <<" r="<<std::setw(4)<<std::fixed<<std::setprecision(3)<<starting_bounding_box.radius()
            <<" te="<<std::setw(7)<<std::scientific<<std::setprecision(1)<<starting_set.time_function().error()<<std::flush
            <<" se="<<std::setw(7)<<std::scientific<<std::setprecision(1)<<sup_norm(starting_set.state_function().errors())<<std::fixed<<std::flush
            <<" l="<<std::left<<starting_set.location()
            <<" e="<<starting_set.previous_events()<<std::setprecision(old_precision));
}

Map<DiscreteEvent,TransitionData>
HybridEvolverBase::
_extract_transitions(DiscreteLocation const& location) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    const EffectiveVectorMultivariateFunction& dynamic=this->system().dynamic_function(location);

    Map<DiscreteEvent,TransitionData> transitions;
    Set<DiscreteEvent> events = this->system().events(location);
    for(Set<DiscreteEvent>::ConstIterator event_iter=events.begin();
        event_iter!=events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;
        EventKind event_kind=this->system().event_kind(location,event);
        ARIADNE_LOG_PRINTLN_AT(1,"event="<<event<<", kind="<<event_kind);
        EffectiveScalarMultivariateFunction constraint_function;
        if(is_activating(event_kind)) {
            constraint_function=this->system().guard_function(location,event);
        } else {
            constraint_function=this->system().invariant_function(location,event);
        }
        ARIADNE_LOG_PRINTLN_AT(1,"constraint_function="<<constraint_function);
        EffectiveScalarMultivariateFunction constraint_flow_derivative_function=lie_derivative(constraint_function,dynamic);
        EffectiveVectorMultivariateFunction reset_function; DiscreteLocation target; RealSpace target_space;
        EffectiveVectorMultivariateFunction target_auxiliary_function; RealSpace target_auxiliary_space;
        if(is_activating(event_kind)) {
            reset_function=this->system().reset_function(location,event);
            target=this->system().target(location,event);
            target_space=this->system().continuous_state_space(target);
            target_auxiliary_function=this->system().auxiliary_function(target);
            target_auxiliary_space=this->system().continuous_auxiliary_space(target);
        }
        TransitionData transition_data={event,event_kind,constraint_function,constraint_flow_derivative_function,target,reset_function,target_space,target_auxiliary_function,target_auxiliary_space};
        transitions.insert(event,transition_data);
    }
    ARIADNE_LOG_PRINTLN("transitions="<<transitions);
    return transitions;
}

Void
HybridEvolverBase::
_apply_invariants(HybridEnclosure& initial_set,
                  Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    HybridEnclosure& invariant_set=initial_set;
    ARIADNE_LOG_PRINTLN("initial_set="<<initial_set);
    ARIADNE_LOG_PRINTLN("transitions="<<transitions);

    // Apply restrictions due to invariants
    for(Map<DiscreteEvent,TransitionData>::ConstIterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event=transition_iter->first;
        ARIADNE_LOG_PRINTLN_AT(1,"event="<<event);
        TransitionData const & transition=transition_iter->second;
        if(transition.event_kind==EventKind::INVARIANT) {
            if (possibly(initial_set.satisfies(transition.guard_function>=0))) {
                invariant_set.new_invariant(event,transition.guard_function);
            }
        }
    }
}

Void
HybridEvolverBase::
_process_starting_events(EvolutionData& evolution_data,
                         HybridEnclosure const& initial_set,
                         Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_ASSERT(evolution_data.working_sets.empty());
    HybridEnclosure invariant_set=initial_set;

    // Apply restrictions due to invariants
    this->_apply_invariants(invariant_set,transitions);

    // Set the flowable set, storing the invariant set as a base for jumps
    HybridEnclosure flowable_set = invariant_set;

    // Compute possibly initially active events
    Set<DiscreteEvent> events=transitions.keys();
    for(Map<DiscreteEvent,TransitionData>::ConstIterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        DiscreteEvent event=transition_iter->first;
        TransitionData const & transition=transition_iter->second;
        // Assume that an impact cannot occur immediately after any other event.
        // This is true after the impact, since $L_{f}g$ changes sign.
        // After a different event, this is not so clear, though it should be
        // a modelling error for the set to intersect the interior of the guard
        // FIXME: Need to consider impact which really can occur immediately due to mapping to boundary.
        // TODO: The condition for an impact to occur is $L_{f}g>0$.
        //       Since this is now available, we should implement this!
        if(transition.event_kind!=EventKind::INVARIANT && transition.event_kind!=EventKind::IMPACT) {
            if(possibly(initial_set.satisfies(transition.guard_function>=0))) {
                if(transition.event_kind!=EventKind::PROGRESS) {
                    HybridEnclosure immediate_jump_set=invariant_set;
                    immediate_jump_set.new_activation(event,transition.guard_function);
                    if(!definitely(immediate_jump_set.is_empty())) {
                        // Put the immediate jump set in the reached sets, since it does not occur in the flowable set
                        ARIADNE_LOG_PRINTLN_AT(1,event<<": "<<transition.event_kind<<", immediate");
                        evolution_data.reach_sets.append(immediate_jump_set);
                        immediate_jump_set.apply_reset(event,transition.target,transition.target_space,transition.reset_function,transition.target_auxiliary_space,transition.target_auxiliary_function);
                        ARIADNE_LOG_PRINTLN_AT(2,"immediate_jump_set="<<immediate_jump_set);
                        evolution_data.intermediate_sets.append(immediate_jump_set);
                        evolution_data.initial_sets.append(immediate_jump_set);
                    }
                }
                if(transition.event_kind!=EventKind::PERMISSIVE) {
                    flowable_set.new_invariant(event,transition.guard_function);
                }
            }
        }
    }

    // Put the flowable set in the starting sets for ordinary evolution
    if(!definitely(flowable_set.is_empty())) {
        evolution_data.working_sets.append(flowable_set);
    }
}

ValidatedVectorMultivariateFunctionModelDP
HybridEvolverBase::
_compute_flow(EffectiveVectorMultivariateFunction dynamic,
              ExactBoxType const& initial_box,
              StepSizeType& previous_step_size,
              const StepSizeType& maximum_step_size) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    IntegratorInterface& integrator=*this->_integrator_ptr;

    // Compute flow and actual time step size used
    //
    // The Integrator classes compute the flow as a function on a symmetrical
    // time domain [-h,+h], since this means the time is centred at 0.
    // We then restrict to the time domain [0,h] since this can make evaluation
    // more accurate, and the time domain might be used explicitly for the domain
    // of the resulting set.
    StepSizeType step_size=maximum_step_size;
    ValidatedVectorMultivariateFunctionModelDP flow_model=integrator.flow_step(dynamic,initial_box,previous_step_size,step_size);
    previous_step_size=step_size;
    ARIADNE_LOG_PRINTLN_AT(1,"twosided_flow_model="<<flow_model);
    ExactBoxType flow_domain=flow_model.domain();
    ARIADNE_ASSERT(step_size==flow_domain[flow_domain.size()-1u].upper_bound());
    flow_domain[flow_domain.size()-1u]=ExactIntervalType(0,step_size);
    flow_model=restriction(flow_model,flow_domain);
    ARIADNE_LOG_PRINTLN_AT(1,"flow_model="<<flow_model);
    ARIADNE_LOG_PRINTLN_AT(1,"flow_model: step_size="<<step_size<<", errors="<<std::scientific<<flow_model.errors()<<", range="<<std::fixed<<flow_model.range());
    ARIADNE_LOG_PRINTLN("flow_model="<<flow_model);
    return flow_model;
}

Set<DiscreteEvent>
HybridEvolverBase::
_compute_active_events(EffectiveVectorMultivariateFunction const& dynamic,
                       Map<DiscreteEvent,EffectiveScalarMultivariateFunction> const& guards,
                       ValidatedVectorMultivariateFunctionModelDP const& flow,
                       HybridEnclosure const& starting_set) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    Set<DiscreteEvent> events=guards.keys();
    Set<DiscreteEvent> active_events;
    HybridEnclosure reach_set=starting_set;
    UpperBoxType flow_bounds=flow.range();
    ARIADNE_LOG_PRINTLN("flow_bounds="<<flow_bounds);
    reach_set.apply_full_reach_step(flow);
    ARIADNE_LOG_PRINTLN("reach_set="<<reach_set);
    for(Set<DiscreteEvent>::Iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
        const DiscreteEvent event=*event_iter;
        const EffectiveScalarMultivariateFunction& guard_function=guards[event];
        EffectiveScalarMultivariateFunction flow_derivative = lie_derivative(guard_function,dynamic);
        ARIADNE_LOG_PRINTLN_AT(1,"event="<<event<<", guard="<<guard_function<<", flow_derivative="<<flow_derivative);
        // First try a simple test based on the bounding box
        UpperIntervalType guard_range=apply(guard_function,reach_set.state_bounding_box());
        // For IMPACT events, also look at flow direction
        UpperIntervalType flow_derivative_range=apply(flow_derivative,reach_set.state_bounding_box());
        ARIADNE_LOG_PRINTLN_AT(1,"guard_range="<<guard_range<<", flow_derivative_range="<<flow_derivative_range);
        if(possibly(guard_range.upper_bound()>=zero)) {
            // Now make a set containing the complement of the constraint,
            // and test for emptiness. If the set is empty, then the guard is
            // not satisfied anywhere.
            HybridEnclosure test_set=reach_set;
            test_set.new_activation(event,guard_function);
            if(!definitely(test_set.is_empty())) {
                active_events.insert(*event_iter);
            }
        }
    }
    return active_events;
}


Map<DiscreteEvent,CrossingData>
HybridEvolverBase::
_compute_crossings(Set<DiscreteEvent> const& active_events,
                   EffectiveVectorMultivariateFunction const& dynamic,
                   Map<DiscreteEvent,EffectiveScalarMultivariateFunction> const& guards,
                   FlowFunctionModel const& flow,
                   HybridEnclosure const& initial_set) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    const SolverInterface& solver=*this->_solver_ptr;

    Map<DiscreteEvent,CrossingData> crossings;
    crossings.clear();

    ExactBoxType flow_spacial_domain=project(flow.domain(),range(0,flow.argument_size()-1u));
    ExactIntervalType flow_time_domain=flow.domain()[flow.argument_size()-1u];
    UpperBoxType flow_bounds=flow.range();
    for(Set<DiscreteEvent>::ConstIterator event_iter=active_events.begin();
        event_iter!=active_events.end(); ++event_iter)
    {
        const DiscreteEvent event=*event_iter;
        ARIADNE_LOG_PRINTLN_AT(1,"event="<<event);
        EffectiveScalarMultivariateFunction const& guard=guards[event];
        ARIADNE_LOG_PRINTLN_AT(2,"guard="<<guard);

        // Compute the value of the guard function g over the flow of $\dot(x)=f(x)$
        UpperIntervalType guard_range=apply(guard,flow_bounds);
        ARIADNE_LOG_PRINTLN_AT(2,"guard_range="<<guard_range);
        // Compute the derivative of the guard function g along flow lines of $\dot(x)=f(x)$
        // This is given by the Lie derivative at a point x, defined as $L_{f}g(x) = (\nabla g\cdot f)(x)$
        EffectiveScalarMultivariateFunction guard_derivative=lie_derivative(guard,dynamic);
        UpperIntervalType guard_derivative_range=apply(guard_derivative,flow_bounds);
        ARIADNE_LOG_PRINTLN_AT(2,"guard_derivative_range="<<guard_derivative_range);
        if(definitely(guard_derivative_range.lower_bound()>zero)) {
            // If the derivative $L_{f}g$is strictly positive over the bounding box for the flow,
            // then the guard function is strictly increasing.
            // There is at most one crossing with the guard, and the time of this
            // crossing must be the time of the event along the trajectory.
            // The crossing time $\gamma(x_0)$ given the initial state can usually be computed
            // by solving the equation $g(\phi(x_0,\gamma(x_0))) = 0$
            ValidatedScalarMultivariateFunctionModelDP crossing_time;
            try {
                crossing_time=solver.implicit(compose(guard,flow),flow_spacial_domain,flow_time_domain);
                if(decide(crossing_time.error()>1e-8)) { ARIADNE_LOG_PRINTLN_AT(2,event<<": crossing_time: error="<<crossing_time.error()<<", range="<<crossing_time.range()); }
                crossings[event]=CrossingData(CrossingKind::TRANSVERSE,crossing_time);
                ARIADNE_LOG_PRINTLN_AT(2,"crossing_time="<<crossing_time);
            }
            catch(const UnknownSolutionException& e) {
                // If the crossing time cannot be computed, then we can still
                // use the fact that the crossing occurs as soon as $g(x(t))=0$.
                ARIADNE_LOG_PRINTLN_AT(2,"Error in computing crossing time for event "<<*event_iter<<":"<<e.what());
                crossings[event]=CrossingData(CrossingKind::INCREASING);
            }
            catch(const SolverException& e) {
                // If the crossing time cannot be computed, then we can still
                // use the fact that the crossing occurs as soon as $g(x(t))=0$.
                ARIADNE_LOG_PRINTLN_AT(2,"Error in computing crossing time for event "<<*event_iter<<": "<<e.what());
                crossings[event]=CrossingData(CrossingKind::INCREASING);
            }
            catch(const std::runtime_error& e) {
                ARIADNE_LOG_PRINTLN_AT(2,"Unexpected error in computing crossing time for event "<<*event_iter<<": "<<e.what());
                ARIADNE_FAIL_MSG("ERROR!!");
                crossings[event]=CrossingData(CrossingKind::INCREASING);
            }
        } else if(definitely(guard_derivative_range.upper_bound()<zero)) {
            // If the guard derivative is strictly negative over the bounding box for the flow,
            // then the guard function is strictly decreasing.
            // This means that the event is either initially active, or does not occur.
            // There is no need to compute a crossing time.
            crossings[event]=CrossingData(CrossingKind::DECREASING);
        } else {
            // If the derivative of the guard function along flow lines cannot be shown
            // to have a definite sign over the entire flow box, then try to compute
            // the sign of the second derivative $L_{f}^{2}g(x)=L_{f}L_{f}g(x)$.
            ValidatedScalarMultivariateFunction guard_second_derivative=lie_derivative(guard_derivative,dynamic);
            UpperIntervalType guard_second_derivative_bounds_range=apply(guard_second_derivative,flow_bounds);
            UpperIntervalType guard_second_derivative_flow_range=compose(guard_second_derivative,flow).range();
            UpperIntervalType guard_second_derivative_range
                =intersection(guard_second_derivative_bounds_range,guard_second_derivative_flow_range);
            ARIADNE_LOG_PRINTLN_AT(2,"guard_second_derivative_range="<<guard_second_derivative_range);
            if(definitely(guard_second_derivative_range.lower_bound()>zero)) {
                // If the second derivative is positive, then either
                //    (i) the event is immediately active
                //   (ii) the event is never active, or
                //  (iii) the event is initially inactive, but becomes active
                //        due to a transverse crossing.
                //   (iv) the initial state is on the boundary of the guard
                //        set, possibly with the flow tangent to this set
                // It is hard compute the crossing time, even in case (iii),
                // due to the singularity due to the tangency in (iv). However,
                // we do know that in (iii), the event occurs when $t>0$ and
                // $g(\phi(x_0,t))=0$. The crossing time is not computed.
                crossings[event]=CrossingData(CrossingKind::CONVEX);
            } else if(definitely(guard_second_derivative_range.upper_bound()<zero)) {
                // If the second derivative is negative, then the guard
                // values $g(x(t))$ are concave along flow lines. There are
                // five main cases:
                //   (i) The event is initially active.
                //  (ii) The event is not initially active, but later becomes active.
                // (iii) The event is never active, but would become active if
                //       flowing backwards in time.
                //  (iv) The event is never active, and the maximum value of
                //       the guard along the flow lines is zero.
                // Additionally, there is the degenerate case
                //  (vi) At some point in the (forward) flow, the state touches
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
                    ValidatedScalarMultivariateFunctionModelDP critical_time=solver.implicit(compose(guard_derivative,flow),flow_spacial_domain,flow_time_domain);
                    UpperIntervalType critical_time_range=critical_time.range();
                    ARIADNE_LOG_PRINTLN_AT(2,"critical_time_range="<<critical_time_range);
                    if(decide(critical_time.error()>1e-8)) {
                        ARIADNE_LOG_PRINTLN_AT(2,event<<": critical_time: error="<<critical_time.error()<<", range="<<critical_time.range()); }

                    HybridEnclosure evolve_set_at_critical_time=initial_set;
                    evolve_set_at_critical_time.apply_space_evolve_step(flow,critical_time);
                    ValidatedVectorMultivariateFunctionModelDP identity = factory(critical_time).create_identity();
                    //ValidatedVectorFunctionModelDP::identity(critical_time.domain(),critical_time.sweeper());
                    UpperBoxType evolve_bounds_at_critical_time=evolve_set_at_critical_time.bounding_box().euclidean_set();
                    UpperIntervalType guard_range_at_critical_time
                        =intersection(apply(guard,evolve_bounds_at_critical_time),evolve_set_at_critical_time.range_of(guard));
                    // Less accurate version: guard_range_at_critical_time=evolve_set_at_critical_time.range_of(guard);
                    ARIADNE_LOG_PRINTLN_AT(2,"evolve_bounds_at_critical_time="<<evolve_bounds_at_critical_time);
                    ARIADNE_LOG_PRINTLN_AT(2,"guard_range_at_critical_time="<<guard_range_at_critical_time);
                    if(definitely(guard_range_at_critical_time.upper_bound()<0)) {
                        // No crossing
                        const_cast<Set<DiscreteEvent>&>(active_events).erase(event); // WARNING: Maybe removing event is unsafe
                    } else if(definitely(guard_range_at_critical_time.lower_bound()>0)) {
                        ARIADNE_LOG_PRINTLN_AT(2,"guard range eventually positive");
                        // Transverse crossing
                        // FIXME: Find a more reliable way of solving the implicit equation for the crossing time
                        //   which takes into account the fact that the derivative over the domain goes negative
                        static const Rational INTERVAL_REDUCTION_FACTOR(15,16);
                        ValidatedScalarMultivariateFunctionModelDP reduced_critical_time=INTERVAL_REDUCTION_FACTOR * critical_time;
                        HybridEnclosure evolve_set_at_reduced_critical_time=initial_set;
                        evolve_set_at_reduced_critical_time.apply_space_evolve_step(flow,reduced_critical_time);
                        HybridEnclosure evolve_set_at_upper_reduced_critical_time=initial_set;
                        evolve_set_at_upper_reduced_critical_time.apply_fixed_evolve_step(flow,static_cast<StepSizeType>(cast_exact(reduced_critical_time.range().upper_bound())));
                        ARIADNE_LOG_PRINTLN_AT(2,"guard_range_at_initial_time="<<initial_set.range_of(guard));
                        ARIADNE_LOG_PRINTLN_AT(2,"guard_range_at_reduced_critical_time="<<evolve_set_at_reduced_critical_time.range_of(guard));
                        ARIADNE_LOG_PRINTLN_AT(2,"guard_derivative_range_at_initial_time="<<initial_set.range_of(guard_derivative));
                        ARIADNE_LOG_PRINTLN_AT(2,"guard_derivative_range_at_upper_reduced_critical_time="<<evolve_set_at_upper_reduced_critical_time.range_of(guard_derivative));
                        ARIADNE_LOG_PRINTLN_AT(2,"guard_derivative_range_at_reduced_critical_time="<<evolve_set_at_reduced_critical_time.range_of(guard_derivative));
                        UpperIntervalType crossing_flow_time_range(0, critical_time_range.upper_bound());
                        ARIADNE_LOG_PRINTLN_AT(2,"crossing_flow_time_range="<<crossing_flow_time_range);
                        ExactIntervalType crossing_flow_time_domain = cast_exact_interval(INTERVAL_REDUCTION_FACTOR*UpperIntervalType(0, critical_time_range.upper_bound()));
                        ARIADNE_LOG_PRINTLN_AT(2,"crossing_flow_time_domain="<<crossing_flow_time_domain);
                        try {
                            //KrawczykSolver solver=KrawczykSolver(1e-10,20);
                            //solver.verbosity=9;//this->verbosity;
                            ValidatedScalarMultivariateFunctionModelDP crossing_time=solver.implicit(compose(guard,flow),flow_spacial_domain,crossing_flow_time_domain);
                            UpperIntervalType crossing_time_range=crossing_time.range();
                            ARIADNE_LOG_PRINTLN_AT(2,"crossing_time_range="<<crossing_time_range);
                            crossings[event]=CrossingData(CrossingKind::TRANSVERSE,crossing_time);
                        } catch(const SolverException& e) {
                            ARIADNE_LOG_PRINTLN_AT(2,"Error in computing crossing time over reduced time interval for event "<<*event_iter<<": " <<e.what());
                            // NOTE: Use GRAZING here since this will later put in a block on continuous evolution at critical time.
                            crossings[event]=CrossingData(CrossingKind::GRAZING);
                            crossings[event].critical_time=critical_time;
                            crossings[event].evolve_bounds_at_critical_time=evolve_bounds_at_critical_time;
                            crossings[event].guard_range_at_critical_time=guard_range_at_critical_time;
                        }
                    } else {
                        crossings[event]=CrossingData(CrossingKind::GRAZING);
                        crossings[event].critical_time=critical_time;
                        crossings[event].evolve_bounds_at_critical_time=evolve_bounds_at_critical_time;
                        crossings[event].guard_range_at_critical_time=guard_range_at_critical_time;
                    }
                }
                catch(const SolverException& e) {
                    ARIADNE_LOG_PRINTLN_AT(2,"Error in computing critical time for event "<<*event_iter<<":\n  "<<e.what());
                    crossings[event]=CrossingData(CrossingKind::CONCAVE);
                }
                catch(const std::runtime_error& e) {
                    ARIADNE_LOG_PRINTLN_AT(2,"Unexpected error in computing critical time for event "<<*event_iter<<":\n  "<<e.what());
                    crossings[event]=CrossingData(CrossingKind::CONCAVE);
                }
            } else {
                // The crossing cannot be shown to be one of the kinds mentioned
                // above. A theoretically exact expression for the crossing
                // set is generally not available.
                crossings[event]=CrossingData(CrossingKind::DEGENERATE);
            }
        }
    }
    if(!crossings.empty()) {
        ARIADNE_LOG_PRINTLN("crossings: "<<log_output(crossings));
    }
    return crossings;
}






Void
HybridEvolverBase::
_recondition(HybridEnclosure& set) const
{
    set.recondition();
}



Void
HybridEvolverBase::
_apply_reach_step(HybridEnclosure& set,
                  ValidatedVectorMultivariateFunctionModelDP const& flow,
                  TimingData const& timing_data) const
{
    set.apply_parameter_reach_step(flow,timing_data.parameter_dependent_evolution_time);
}

Void
HybridEvolverBase::
_apply_evolve_step(HybridEnclosure& set,
                  ValidatedVectorMultivariateFunctionModelDP const& flow,
                  TimingData const& timing_data) const
{

    switch(timing_data.step_kind) {
        case StepKind::CONSTANT_EVOLUTION_TIME:
        case StepKind::SPACE_DEPENDENT_EVOLUTION_TIME:
        case StepKind::TIME_DEPENDENT_EVOLUTION_TIME:
        case StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME:
        case StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME:
            set.apply_parameter_evolve_step(flow,timing_data.parameter_dependent_evolution_time);
            break;
        case StepKind::PARAMETER_DEPENDENT_FINISHING_TIME:
        case StepKind::SPACETIME_DEPENDENT_FINISHING_TIME:
        case StepKind::CONSTANT_FINISHING_TIME:
            set.apply_finishing_parameter_evolve_step(flow,timing_data.parameter_dependent_finishing_time);
            break;
        default:
            ARIADNE_FAIL_MSG("Unhandled step kind "<<timing_data.step_kind);
    }
}

Void
HybridEvolverBase::
_apply_guard_step(HybridEnclosure& set,
                  EffectiveVectorMultivariateFunction const& dynamic,
                  ValidatedVectorMultivariateFunctionModelDP const& flow,
                  TimingData const& timing_data,
                  TransitionData const& transition_data,
                  CrossingData const& crossing_data,
                  const Semantics semantics) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    // Compute flow to guard set up to evolution time.
    HybridEnclosure& jump_set=set;
    const DiscreteEvent event=transition_data.event;
    ValidatedVectorMultivariateFunctionModelDP starting_state=set.state_function();
    ValidatedVectorMultivariateFunctionModelDP reach_starting_state=embed(starting_state,timing_data.evolution_time_domain);
    ValidatedScalarMultivariateFunctionModelDP reach_step_time=embed(starting_state.domain(),timing_data.evolution_time_coordinate);
    ValidatedScalarMultivariateFunctionModelDP step_time, step_critical_time;

    ARIADNE_LOG_PRINTLN("transition_data.event_kind="<<transition_data.event_kind);
    switch(transition_data.event_kind) {
        case EventKind::PERMISSIVE:
            // The continuous evolution is just the same as a reachability step,
            // so we need to embed the starting state and the step time function into one higher dimension.
            jump_set.apply_parameter_reach_step(flow,timing_data.parameter_dependent_evolution_time);
            jump_set.new_activation(event,transition_data.guard_function);
            break;
        case EventKind::URGENT: case EventKind::IMPACT:
            ARIADNE_LOG_PRINTLN("crossing_data.crossing_kind="<<crossing_data.crossing_kind);
            switch(crossing_data.crossing_kind) {
                case CrossingKind::TRANSVERSE:
                    // The jump set is given by \f$\xi'(s)=\phi(\xi(s),\gamma(s)),\ \tau'(s)=\tau(s)+\gamma(s)\f$ where \f$\gamma(s)\f$ is the crossing time.
                    ARIADNE_LOG_PRINTLN_AT(1,"crossing_time="<<crossing_data.crossing_time);
                    step_time=unchecked_compose(crossing_data.crossing_time,starting_state);
                    // If the jump step might occur after the final evolution time, then introduce constraint that this does not happen
                    if(timing_data.step_kind!=StepKind::CONSTANT_EVOLUTION_TIME && possibly((step_time-timing_data.parameter_dependent_evolution_time).range().upper_bound()>zero)) {
                        jump_set.new_parameter_constraint(step_event,step_time<=timing_data.parameter_dependent_evolution_time);
                    }
                    jump_set.apply_parameter_evolve_step(flow,unchecked_compose(crossing_data.crossing_time,starting_state));
                    break;
                case CrossingKind::INCREASING: case CrossingKind::CONVEX:
                    // The jump set is given by \f$\xi'(s,t)=\phi(\xi(s),t),\ \tau'(s)=\tau(s)+t\f$ where \f$t\in[0,h]\f$ and \f$g(\xi'(s,t))=0\f$.
                    jump_set.apply_parameter_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                    jump_set.new_guard(event,transition_data.guard_function);
                    break;
                case CrossingKind::GRAZING:
                    ARIADNE_LOG_PRINTLN_AT(1,"critical_time="<<crossing_data.critical_time);
                    // Fallthrough
                case CrossingKind::CONCAVE:
                    // The jump set is given by \f$\xi'(s,t)=\phi(\xi(s),t),\ \tau'(s)=\tau(s)+t\f$ where \f$t\in[0,h]\f$, \f$g(\xi'(s,t))=0\f$
                    // and \f$L_f{g}(\xi'(s,t))>=0\f$ (equivalently, \f$t<=\gamma(\xi(s))\f$ where \f$\gamma(x)\f$ is the critical time
                    // where \f$g(x(t))\f$ reaches a maximum.
                    jump_set.apply_parameter_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                    jump_set.new_guard(event,transition_data.guard_function);
                    jump_set.new_invariant(event,-lie_derivative(transition_data.guard_function,dynamic));
                    break;
/*
                case CrossingKind::GRAZING:
                    ARIADNE_LOG(6,"critical_time="<<crossing_data.critical_time);
                    ARIADNE_LOG(9,"jump_set.domain()="<<jump_set.domain());
                    IntervalDomainType evolution_time_domain=timing_data.evolution_time_domain;
                    ValidatedScalarMultivariateFunctionModelDP embedded_space_function=embed(set.space_function(),timing_data.evolution_time_domain);
                    jump_set.apply_parameter_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                    ValidatedScalarMultivariateFunctionModelDP embedded_time_step_function=factory(embedded_space_function).create_coordinate(jump_set.number_of_parameters()-1u);
                    jump_set.new_parameter_constraint(event,embedded_time_step_function<=compose(crossing_data.critical_time,embedded_space_function));
                    jump_set.new_guard(event,transition_data.guard_function);
                    jump_set.reduce(); // Reduce the size of the parameter domain to take guards into account
                    break;
*/
                case CrossingKind::DEGENERATE: // Just check positive derivative in this case; NOT EXACT
                    if(semantics==Semantics::UPPER) {
                        jump_set.apply_parameter_reach_step(flow,timing_data.parameter_dependent_evolution_time);
                        jump_set.new_guard(event,transition_data.guard_function);
                        jump_set.new_invariant(event,-lie_derivative(transition_data.guard_function,dynamic));
                    } else {
                        // Make the empty set
                        jump_set.new_invariant(event,transition_data.guard_function*0+1);
                    }
                    break;
                case CrossingKind::NEGATIVE:
                case CrossingKind::POSITIVE:
                case CrossingKind::DECREASING:
                    {
                        ARIADNE_WARN("Crossing "<<crossing_data<<" does not introduce additional restrictions on flow");
                    }
                    break;
                default:
                    ARIADNE_FAIL_MSG("Unhandled crossing kind in "<<crossing_data);
            }
            if (transition_data.event_kind==EventKind::IMPACT) {
                // For an IMPACT event, also impose that the flow derivative is positive at the crossing.
                EffectiveScalarMultivariateFunction flow_derivative_function=lie_derivative(transition_data.guard_function,dynamic);
                ARIADNE_LOG_PRINTLN_AT(1,"flow_derivative_function="<<flow_derivative_function);
                ARIADNE_LOG_PRINTLN_AT(1,"jump_set.state_bounding_box()="<<jump_set.state_bounding_box());
                UpperIntervalType flow_derivative_range=apply(flow_derivative_function,jump_set.state_bounding_box());
                ARIADNE_LOG_PRINTLN_AT(1,"flow_derivative_range="<<flow_derivative_range);
                if(possibly(flow_derivative_range.lower_bound()<=0)) {
                    jump_set.new_activation(event,flow_derivative_function);
                }
                ARIADNE_LOG_PRINTLN_AT(1,"jump_set="<<jump_set);
                if(definitely(jump_set.is_empty())) { return; }
            }
            break;
        default:
            ARIADNE_FAIL_MSG("Invalid event kind "<<transition_data.event_kind<<" for transition.");
    }



}


// Apply guard to a single set.
// In the case of concave crossings, splits the set into two, one part
// corresponding to points which actually hit the set (and stop on first crossing)
// the other part corresponding to points which miss the set.
Void HybridEvolverBase::
_apply_guard(List<HybridEnclosure>& sets,
             const ValidatedScalarMultivariateFunctionModelDP& elapsed_time_function,
             const HybridEnclosure& starting_set,
             const ValidatedVectorMultivariateFunctionModelDP& flow,
             const TransitionData& transition_data,
             const CrossingData guard_crossing_data,
             const Semantics semantics) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("guard_crossing_data="<<log_output(guard_crossing_data));
    // Multiple sets may be input, but they must all have the parametrised elapsed time function
    static const Nat SUBDIVISIONS_FOR_DEGENERATE_CROSSING = 2;
    const DiscreteEvent event=transition_data.event;
    const ValidatedScalarMultivariateFunction& guard_function=transition_data.guard_function;

    ValidatedVectorMultivariateFunctionModelDP starting_state_function=starting_set.state_function();
    if(elapsed_time_function.domain().dimension()>starting_set.parameter_domain().dimension()) {
        IntervalDomainType elapsed_time_domain=elapsed_time_function.domain()[elapsed_time_function.argument_size()-1u];
        starting_state_function = embed(starting_state_function,elapsed_time_domain);
    }

    List<HybridEnclosure>::Iterator end=sets.end(); // Store current end set since we may adjoin new sets at end
    for(List<HybridEnclosure>::Iterator iter=sets.begin(); iter!=end; ++iter) {
        HybridEnclosure& set=*iter;
        ARIADNE_LOG_PRINTLN_AT(1,"set.parameter_domain()="<<set.parameter_domain());
        ARIADNE_ASSERT(set.parameter_domain()==elapsed_time_function.domain());
        switch(guard_crossing_data.crossing_kind) {
            case CrossingKind::TRANSVERSE:
                ARIADNE_ASSERT(set.parameter_domain()==starting_state_function.domain());
                set.new_parameter_constraint(event, elapsed_time_function <= unchecked_compose(guard_crossing_data.crossing_time,starting_state_function));
                // Alternatively:
                // set.new_invariant(event, guard_function);
                // set.new_state_constraint(event, guard_function <= zero);
                break;
            case CrossingKind::CONVEX: case CrossingKind::INCREASING:
                //set.new_state_constraint(event, guard_function <= zero);
                set.new_invariant(event, guard_function);
                break;
            case CrossingKind::GRAZING: {
                ValidatedScalarMultivariateFunctionModelDP critical_time_function = unchecked_compose(guard_crossing_data.critical_time,starting_state_function);
                ValidatedScalarMultivariateFunctionModelDP final_guard_function
                    = compose( guard_function, unchecked_compose( flow, join(starting_state_function, elapsed_time_function) ) );
                ValidatedScalarMultivariateFunctionModelDP maximal_guard_function
                    = compose( guard_function, unchecked_compose( flow, join(starting_state_function, critical_time_function) ) );
                UpperIntervalType guard_range_at_critical_time = guard_crossing_data.guard_range_at_critical_time;
                ARIADNE_ASSERT_MSG(starting_state_function.argument_size()==set.parameter_domain().size(),
                                   starting_state_function<<" "<<set);
                ARIADNE_ASSERT_MSG(critical_time_function.argument_size()==set.parameter_domain().size(),
                                   critical_time_function<<" "<<set);
                ARIADNE_ASSERT_MSG(maximal_guard_function.argument_size()==set.parameter_domain().size(),
                                   maximal_guard_function<<" "<<set);
                ARIADNE_ASSERT_MSG(final_guard_function.argument_size()==set.parameter_domain().size(),
                                   final_guard_function<<" "<<set);

                // If no points in the set arise from trajectories which will later leave the progress set,
                // then we only need to look at the maximum value of the guard.
                HybridEnclosure eventually_hitting_set=set;
                eventually_hitting_set.new_parameter_constraint( event, elapsed_time_function <= critical_time_function );
                eventually_hitting_set.new_parameter_constraint( event, maximal_guard_function >= zero);
                ARIADNE_LOG_PRINTLN_AT(1,"eventually_hitting_set.is_empty()="<<eventually_hitting_set.is_empty());
                if(definitely(eventually_hitting_set.is_empty())) {
                    // NOTE: Constraint below should always be satisfied
                    // set.new_parameter_constraint(event, maximal_guard_function <= zero);
                    break;
                }
                // If no points in the set arise from trajectories which leave the progress set and
                // later return, then we only need to look at the guard at the final value
                HybridEnclosure returning_set=set;
                returning_set.new_parameter_constraint( event, elapsed_time_function >= critical_time_function );
                returning_set.new_parameter_constraint( event, final_guard_function <= zero );
                ValidatedLowerKleenean returning_set_empty=returning_set.is_empty();
                ARIADNE_LOG_PRINTLN_AT(1,"returning_set.is_empty()="<<returning_set_empty);
                if(definitely(returning_set_empty)) {
                    set.new_parameter_constraint( event, final_guard_function <= zero );
                    break;
                }

                // Split the set into three components,
                //   hitting_set: points which hit the guard in the future
                //   missing_set: points whose trajectories miss the guard completely
                //   past_set: points which would have hit the guard in the past, but do not in the future
                HybridEnclosure& missing_set=set;
                HybridEnclosure hitting_set=set;
                HybridEnclosure past_set=set;
                missing_set.new_parameter_constraint( event, maximal_guard_function <= zero );
                hitting_set.new_parameter_constraint( event, maximal_guard_function >= zero );
                hitting_set.new_parameter_constraint( event, elapsed_time_function <= critical_time_function );
                hitting_set.new_parameter_constraint( event, final_guard_function <= zero );
                past_set.new_parameter_constraint( event, maximal_guard_function >= zero );
                past_set.new_parameter_constraint( event, critical_time_function <= zero );
                hitting_set.reduce();
                missing_set.reduce();
                past_set.reduce();
                ARIADNE_LOG_PRINTLN_AT(1,"final_guard_range="<<final_guard_function.range());
                ARIADNE_LOG_PRINTLN_AT(1,"maximal_guard_function="<<maximal_guard_function);
                ARIADNE_LOG_PRINTLN_AT(1,"final_guard_function="<<final_guard_function);
                ARIADNE_LOG_PRINTLN_AT(1,"elapsed_time_function="<<elapsed_time_function);
                ARIADNE_LOG_PRINTLN_AT(1,"critical_time_function="<<critical_time_function);
                ARIADNE_LOG_PRINTLN_AT(1,"missing_set.is_empty()="<<missing_set.is_empty()<<", hitting_set.is_empty()="<<hitting_set.is_empty()<<", past_set.is_empty()="<<past_set.is_empty());
                ARIADNE_LOG_PRINTLN_AT(1,"hitting_set="<<hitting_set);
                // FIXME: The guard range at the critical time may provide a stronger guarantee of emptiness than maximal_guard_function > zero, but we should ideally not rely on this
                if(definitely(guard_range_at_critical_time > 0 || missing_set.is_empty())) {
                    // swap out missing set
                    std::swap(hitting_set,missing_set);
                } else if(not definitely(hitting_set.is_empty())) {
                    sets.append(hitting_set);
                }
                if (not definitely(past_set.is_empty())) {
                    sets.append(past_set);
                }
                break;
            }
            case CrossingKind::DEGENERATE: case CrossingKind::CONCAVE: {
                // The crossing with the guard set is not one of the kinds handled above.
                // We obtain an over-approximation by testing at finitely many time points
                const Nat n=SUBDIVISIONS_FOR_DEGENERATE_CROSSING;
                switch(semantics) {
                    case Semantics::UPPER:
                        for(Nat i=0; i!=n; ++i) {
                            FloatDPBounds alpha=FloatDPValue(i+1,dp)/n;
                            ValidatedScalarMultivariateFunctionModelDP intermediate_guard
                                = compose( guard_function, unchecked_compose( flow, join(starting_state_function, alpha*elapsed_time_function) ) );
                            set.new_parameter_constraint(event, intermediate_guard <= zero);
                        }
                        break;
                    case Semantics::LOWER:
                        // Can't continue the evolution, so set a trivially-falsified constraint
                        set.new_parameter_constraint(event, this->function_factory().create_constant(set.parameter_domain(),1) <= zero);
                        break;
                    default:
                        ARIADNE_FAIL_MSG("Unhandled semantics.");
                }
                break;
            }
            case CrossingKind::POSITIVE:
                // No need to do anything since all points are initially
                // active and should have been handled already
            case CrossingKind::NEGATIVE:
                // No points are active
            case CrossingKind::DECREASING:
                // No need to do anything, since only initially active points
                // become active during the evolution, and these have been
                // handled already.
                break;
            default:
                ARIADNE_FAIL_MSG("Unhandled crossing "<<guard_crossing_data);
        }
    }
}



Void
HybridEvolverBase::
_evolution_in_mode(EvolutionData& evolution_data,
                   HybridTerminationCriterion const& termination_criterion) const
{

    //  Select a working set and evolve this in the current location until either
    // all initial points have undergone a discrete transition (possibly to
    // the same location) or the final time is reached.
    //   Evolving within one location avoids having to re-extract event sets,
    // and means that initially active events are tested for only once.
    ARIADNE_LOG_SCOPE_CREATE;

    const Real final_time=termination_criterion.maximum_time();
    const Integer maximum_steps=termination_criterion.maximum_steps();
    const Set<DiscreteEvent>& terminating_events=termination_criterion.terminating_events();

    // Routine check for emptiness
    if(evolution_data.initial_sets.empty()) { return; }

    // Get the initial set for this round of evolution
    HybridEnclosure initial_set=evolution_data.initial_sets.back(); evolution_data.initial_sets.pop_back();
    ARIADNE_LOG_PRINTLN_AT(2,"initial_set="<<initial_set);

    ARIADNE_LOG_RUN_MUTED(auto time_range = initial_set.time_range());
    if(definitely(time_range.lower_bound()>=final_time)) {
        ARIADNE_WARN("initial_set.time_range()="<<time_range<<" which exceeds final time="<<final_time);
        return;
    }

    // Extract starting location
    const DiscreteLocation location=initial_set.location();

    // Cache dynamic and constraint functions
    EffectiveVectorMultivariateFunction dynamic=this->system().dynamic_function(location);
    ARIADNE_LOG_RUN_AT(1,auto transitions = this->_extract_transitions(location));
    Set<DiscreteEvent> events = transitions.keys();

    ARIADNE_LOG_PRINTLN_AT(2,"dynamic="<<dynamic);
    ARIADNE_LOG_PRINTLN_AT(2,"initial_set "<<initial_set);
    ARIADNE_LOG_PRINTLN_AT(2,"transitions="<<transitions);

    // Test if maximum number of steps has been exceeded; if so, the set should be discarded.
    // NOTE: We could also place a test for the maximum number of steps being reaches which computing jump sets
    // This is not done since the maximum_steps information is not passed to the _apply_evolution_step(...) method.
    if(initial_set.previous_events().size()>=maximum_steps) {
        ARIADNE_LOG_PRINTLN_AT(2,"initial_set has undergone more than the maximum number of events ("<<maximum_steps<<")");
        ARIADNE_LOG_RUN_AT(1,this->_apply_invariants(initial_set,transitions));
        evolution_data.final_sets.append(initial_set);
        return;
    }

    // Test if a terminating event has been reached.
    if(initial_set.previous_events().size()>=1 && terminating_events.contains(initial_set.previous_events().back())) {
        ARIADNE_LOG_PRINTLN_AT(2,"initial_set has undergone event "<<initial_set.previous_events().back());
        ARIADNE_LOG_RUN_AT(1,this->_apply_invariants(initial_set,transitions));
        evolution_data.final_sets.append(initial_set);
        return;
    }

    // Process the initially active events; cut out active points to leave initial flowable set.
    ARIADNE_LOG_PRINTLN_AT(2,"initial_set has not reached the discrete termination condition");

    this->_process_starting_events(evolution_data, initial_set,transitions);
    ARIADNE_ASSERT(evolution_data.working_sets.size()<=1);

    // NOTE: Uncomment the lines below to stop evolution immediately after the maximum event, without further flow
    //if(initial_set.previous_events().size()>maximum_steps) {
    //    evolution_data.final_sets.append(initial_set);
    //    return;
    //}

    StepSizeType previous_step_size=0;
    while(!evolution_data.working_sets.empty()) {
        this->_evolution_step(evolution_data,dynamic,transitions,final_time,previous_step_size);
    }
}

Void
HybridEvolverBase::
_evolution_step(EvolutionData& evolution_data,
                EffectiveVectorMultivariateFunction const& dynamic,
                Map<DiscreteEvent,TransitionData> const& transitions,
                Real const& final_time,
                StepSizeType& previous_step_size) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    HybridEnclosure starting_set=evolution_data.working_sets.back(); evolution_data.working_sets.pop_back();

    _log_summary(evolution_data,starting_set);

    ARIADNE_LOG_PRINTLN_AT(2,"starting_set="<<starting_set);
    ARIADNE_LOG_PRINTLN_AT(2,"starting_time="<<starting_set.time_function());
    if(definitely(starting_set.is_empty())) {
        ARIADNE_LOG_PRINTLN_AT(1,"Empty starting_set "<<starting_set);
        return;
    }

    if(definitely(starting_set.time_range().lower_bound() >= final_time)) {
        ARIADNE_WARN("starting_set.time_range()="<<starting_set.time_range()<<" which exceeds final time="<<final_time);
        return;
    }

    // Compute the bounding box of the enclosure
    const ExactBoxType starting_bounding_box=cast_exact_box(starting_set.state_bounding_box());
    ARIADNE_LOG_PRINTLN_AT(2,"starting_bounding_box="<<starting_bounding_box);

    // Test to see if set requires reconditioning
    if(this->_configuration_ptr->enable_reconditioning() &&
            possibly(norm(starting_set.state_function().errors()) > this->_configuration_ptr->maximum_spacial_error())) {
        ARIADNE_LOG_PRINTLN_AT(1,"reconditioning: errors "<<starting_set.state_function().errors());
        HybridEnclosure reconditioned_set=starting_set;
        reconditioned_set.recondition();
        evolution_data.working_sets.append(reconditioned_set);
        return;
    }

    // Handle a set that is too large, based on semantics
    if (possibly(starting_bounding_box.radius() > this->_configuration_ptr->maximum_enclosure_radius())) {
        if (evolution_data.semantics == Semantics::LOWER) {
            ARIADNE_LOG_PRINTLN_AT(1,"Set too large, discarding");
            return;
        } else if (this->_configuration_ptr->enable_subdivisions()) {
            ARIADNE_LOG_PRINTLN_AT(1,"Set too large, splitting");
            List<HybridEnclosure> split_sets = starting_set.split();
            for(SizeType i=0; i!=split_sets.size(); ++i) {
                if(!definitely(split_sets[i].is_empty())) { evolution_data.working_sets.append(split_sets[i]); }
            }
            return;
        }
    }

    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guard_functions;
    for(Map<DiscreteEvent,TransitionData>::ConstIterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
        guard_functions.insert(transition_iter->first,transition_iter->second.guard_function);
    }
    ARIADNE_LOG_PRINTLN_AT(1,"guards="<<guard_functions);

    // Compute flow and actual time step size used
    const FlowFunctionModel flow_model=this->_compute_flow(dynamic,starting_bounding_box,previous_step_size,this->configuration().maximum_step_size());
    ARIADNE_LOG_PRINTLN_AT(1,"flow_model.domain()="<<flow_model.domain()<<" flow_model.range()="<<flow_model.range());

    // Compute possibly active events
    Set<DiscreteEvent> active_events =
        this->_compute_active_events(dynamic,guard_functions,flow_model,starting_set);
    if(!active_events.empty()) {
        Map<DiscreteEvent,EventKind> active_events_log_data;
        for(Set<DiscreteEvent>::ConstIterator event_iter=active_events.begin(); event_iter!=active_events.end(); ++event_iter) {
            active_events_log_data.insert(*event_iter,transitions[*event_iter].event_kind);
        }
        ARIADNE_LOG_PRINTLN_AT(1,"active_events: "<<active_events_log_data);
    }

    // Compute the kind of crossing (increasing, convex, etc);
    Map<DiscreteEvent,CrossingData> crossings =
        this->_compute_crossings(active_events,dynamic,guard_functions,flow_model,starting_set);
    ARIADNE_LOG_PRINTLN_AT(1,"crossings="<<crossings);

    // Compute end conditions for flow
    TimingData timing_data = this->_estimate_timing(active_events,Real(final_time),flow_model,crossings,transitions,starting_set);
    ARIADNE_LOG_PRINTLN_AT(1,"timing_data: "<<timing_data);

    // Apply the time step
    HybridEnclosure reach_set, evolve_set;
    this->_apply_evolution_step(evolution_data,starting_set,flow_model,timing_data,crossings,dynamic,transitions);
}


Void HybridEvolverBase::
_apply_evolution_step(EvolutionData& evolution_data,
                      HybridEnclosure const& starting_set,
                      ValidatedVectorMultivariateFunctionModelDP const& flow,
                      TimingData const& timing_data,
                      Map<DiscreteEvent,CrossingData> const& crossings,
                      EffectiveVectorMultivariateFunction const& dynamic,
                      Map<DiscreteEvent,TransitionData> const& transitions) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    EvolutionStepData _step_data;
    HybridEnclosure starting_set_copy=starting_set;
    ValidatedLowerKleenean starting_set_empty=starting_set_copy.is_empty();

    if(definitely(starting_set_empty)) {
        ExactBoxType reduced_domain=starting_set.continuous_set().reduced_domain();
        ARIADNE_WARN("empty starting_set "<<representation(starting_set));
        return;
    }
    //ARIADNE_ASSERT_MSG(!definitely(starting_set.empty()),"starting_set="<<repr(starting_set));

    Semantics semantics = evolution_data.semantics;

    // Compute events enabling transitions and events blocking continuous evolution
    Set<DiscreteEvent> active_events=crossings.keys();
    Set<DiscreteEvent> activating_events=intersection(Ariadne::activating_events(transitions),active_events);
    Set<DiscreteEvent> blocking_events=intersection(Ariadne::blocking_events(transitions),active_events);

    ARIADNE_LOG_PRINTLN("activating_events="<<activating_events);
    ARIADNE_LOG_PRINTLN("blocking_events="<<blocking_events);

    ValidatedScalarMultivariateFunctionModelDP const evolve_step_time=timing_data.parameter_dependent_evolution_time;
    ValidatedScalarMultivariateFunctionModelDP const reach_step_time=embed(starting_set.parameter_domain(),timing_data.evolution_time_coordinate);

    ARIADNE_LOG_PRINTLN_AT(1,"evolve_step_time="<<evolve_step_time<<"\n")
    ARIADNE_LOG_PRINTLN_AT(1,"reach_step_time="<<reach_step_time<<"\n")

    // Compute the reach and evolve sets, without introducing bounds due to the final time.
    List<HybridEnclosure> reach_sets={starting_set};
    List<HybridEnclosure> evolve_sets={starting_set};

    _apply_reach_step(reach_sets.front(),flow,timing_data);
    _apply_evolve_step(evolve_sets.front(),flow,timing_data);

    ARIADNE_LOG_PRINTLN_AT(1,"flow_reach_set="<<reach_sets.front()<<"\n")
    ARIADNE_LOG_PRINTLN_AT(1,"flow_evolve_set="<<evolve_sets.front()<<"\n")

    // Apply constraints on reach and evolve sets due to invariants and urgent guards
    for(Set<DiscreteEvent>::ConstIterator event_iter=blocking_events.begin();
        event_iter!=blocking_events.end(); ++event_iter)
    {
        const TransitionData& transition_data=transitions[*event_iter];
        const CrossingData& crossing_data=crossings[*event_iter];

        this->_apply_guard(reach_sets,reach_step_time,starting_set,flow,
                           transition_data,crossing_data,semantics);

        switch(crossing_data.crossing_kind) {
            case CrossingKind::INCREASING: case CrossingKind::TRANSVERSE:
                // Delay applying guard until all splittings due to non-increasing crossings have been processed.
                // We can then test emptiness on sets with crossings, but not actually apply crossing immediately,
                // since it will be introduced at the next time step anyway; this avoids duplication
                break;
            default:
                this->_apply_guard(evolve_sets,evolve_step_time,starting_set,flow,
                                   transition_data,crossing_data,semantics);
        }
    }

    // Make copy of evolve sets without transverse and increasing crossings which can be used for future evolution
    List<HybridEnclosure> next_working_sets = evolve_sets;

    for(Set<DiscreteEvent>::ConstIterator event_iter=blocking_events.begin();
        event_iter!=blocking_events.end(); ++event_iter)
    {
        const TransitionData& transition_data=transitions[*event_iter];
        const CrossingData& crossing_data=crossings[*event_iter];
        switch (crossing_data.crossing_kind) {
            case CrossingKind::INCREASING: case CrossingKind::TRANSVERSE:
                this->_apply_guard(evolve_sets,evolve_step_time,starting_set,flow,
                                   transition_data,crossing_data,semantics);
                break;
            default:
                break; // Constraint has already been handled
        }
    }

    // Compute final set depending on whether the finishing kind is exactly AT_FINAL_TIME.
    // Insert sets into evolution_data as appropriate
    if(timing_data.finishing_kind==FinishingKind::AT_FINAL_TIME) {
        for(List<HybridEnclosure>::ConstIterator evolve_set_iter=evolve_sets.begin();
            evolve_set_iter!=evolve_sets.end(); ++evolve_set_iter)
        {
            HybridEnclosure const& evolve_set=*evolve_set_iter;
            if(!definitely(evolve_set.is_empty())) {
                ARIADNE_LOG_PRINTLN("final_set="<<evolve_set);
                evolution_data.final_sets.append(evolve_set);
                _step_data.finishing = true;
            }
        }
        for(List<HybridEnclosure>::ConstIterator reach_set_iter=reach_sets.begin();
            reach_set_iter!=reach_sets.end(); ++reach_set_iter)
        {
            HybridEnclosure const& reach_set=*reach_set_iter;
            evolution_data.reach_sets.append(reach_set);
        }
    } else { // (timing_data.finishing_kind!=AT_FINAL_TIME)
        List<HybridEnclosure>::Iterator next_working_set_iter=next_working_sets.begin();
        for(List<HybridEnclosure>::Iterator evolve_set_iter=evolve_sets.begin();
            evolve_set_iter!=evolve_sets.end(); ++evolve_set_iter, ++next_working_set_iter )
        {
            HybridEnclosure& evolve_set=*evolve_set_iter;
            HybridEnclosure& next_working_set=*next_working_set_iter;
            UpperIntervalType evolve_set_time_range=evolve_set.time_range();
            if(definitely(evolve_set_time_range.lower_bound()>timing_data.final_time)) {
                // Do nothing, since evolve set is past final time is definitely empty
            } else if(definitely(evolve_set_time_range.upper_bound()<=timing_data.final_time)) {
                // No need to introduce timing constraints
                if(!definitely(evolve_set.is_empty())) {
                    ARIADNE_LOG_PRINTLN("evolve_set="<<evolve_set);
                    ARIADNE_LOG_PRINTLN("next_working_set="<<next_working_set);
                    evolution_data.working_sets.append(next_working_set);
                    evolution_data.intermediate_sets.append(evolve_set);
                    _step_data.progress = true;
                }
            } else {
                // Bound time if necessary
                if(possibly(evolve_set_time_range.upper_bound()>timing_data.final_time)) {
                    evolve_set.bound_time(timing_data.final_time);
                }
                // Only continue evolution if the time-bounded evolve set is nonempty;
                // However, continue evolution without adding time constraint,
                // since this will be introduced in the next step
                if(!definitely(evolve_set.is_empty())) {
                    ARIADNE_LOG_PRINTLN("evolve_set="<<evolve_set);
                    ARIADNE_LOG_PRINTLN("next_working_set="<<next_working_set);
                    evolution_data.working_sets.append(next_working_set);
                    evolution_data.intermediate_sets.append(evolve_set);
                    _step_data.progress = true;
                }
            }
        }
        for(List<HybridEnclosure>::Iterator reach_set_iter=reach_sets.begin();
            reach_set_iter!=reach_sets.end(); ++reach_set_iter)
        {
            HybridEnclosure& reach_set=*reach_set_iter;
            //reach_set.continuous_set().reduce();
            if(possibly(reach_set.time_range().upper_bound()>timing_data.final_time)) {
                HybridEnclosure final_set=reach_set;
                final_set.set_time(timing_data.final_time);
                if(timing_data.finishing_kind!=FinishingKind::BEFORE_FINAL_TIME && !definitely(final_set.is_empty())) {
                //if(!definitely(final_set.is_empty())) {
                    ARIADNE_LOG_PRINTLN("final_set="<<final_set);
                    evolution_data.final_sets.append(final_set);
                    _step_data.finishing = true;
                }
            }
            reach_set.bound_time(timing_data.final_time);
            evolution_data.reach_sets.append(reach_set);
        }
    }

    ARIADNE_LOG_PRINTLN("Computing jump sets");

    // Compute jump sets
    for(Set<DiscreteEvent>::ConstIterator event_iter=activating_events.begin();
        event_iter!=activating_events.end(); ++event_iter)
    {
        DiscreteEvent event=*event_iter;


        // Compute active set
        List<HybridEnclosure> jump_sets={starting_set};
        HybridEnclosure& jump_set=jump_sets.front();
        ARIADNE_LOG_PRINTLN_AT(1,event<<": "<<transitions[event].event_kind<<", "<<crossings[event].crossing_kind);
        _apply_guard_step(jump_set,dynamic,flow,timing_data,transitions[event],crossings[event],semantics);
        ValidatedScalarMultivariateFunctionModelDP jump_step_time=reach_step_time;
        if(reach_step_time.argument_size()!=jump_set.number_of_parameters()) {
            ARIADNE_ASSERT(starting_set.number_of_parameters()==jump_set.number_of_parameters());
            switch(crossings[event].crossing_kind) {
                case CrossingKind::TRANSVERSE:
                    jump_step_time=unchecked_compose(crossings[event].crossing_time,starting_set.state_function());
                    break;
                case CrossingKind::INCREASING: case CrossingKind::CONVEX: case CrossingKind::CONCAVE: case CrossingKind::DEGENERATE:
                    // Apply guard using equation
                case CrossingKind::GRAZING:
                    // Apply guard using equation; critical time not used for flow time
                case CrossingKind::DECREASING: case CrossingKind::POSITIVE: case CrossingKind::NEGATIVE:
                    // Crossing time should not be used
                    break;
                default:
                    ARIADNE_FAIL_MSG("Unknown crossing kind in "<<crossings[event]);
            }
        }

        // Apply maximum time bound, as this will be applied after the next flow step
        for(List<HybridEnclosure>::Iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& _jump_set=*jump_set_iter;
            if(possibly(_jump_set.time_range().upper_bound()>timing_data.final_time)) {
                _jump_set.bound_time(timing_data.final_time);
                if(definitely(not _jump_set.is_empty())) { ARIADNE_WARN("Explicitly bounding time in jump set"); }
            }
        }

        ARIADNE_LOG_PRINTLN_AT(1,"Applying guard to jump sets");
        // Apply blocking conditions for other active events
        for(Set<DiscreteEvent>::ConstIterator other_event_iter=blocking_events.begin();
            other_event_iter!=blocking_events.end(); ++other_event_iter)
        {
            ARIADNE_LOG_PRINTLN_AT(2,"jump_sets.size()="<<jump_sets.size());
            DiscreteEvent other_event=*other_event_iter;
            if(other_event!=event) {
                const TransitionData& other_transition_data=transitions[other_event];
                const CrossingData& other_crossing_data=crossings[other_event];
                _apply_guard(jump_sets,jump_step_time,starting_set,flow,
                             other_transition_data,other_crossing_data,semantics);
            }
        }

        ARIADNE_LOG_PRINTLN_AT(1,event<<": "<<transitions[event].event_kind<<", "<<crossings[event].crossing_kind);
        // Apply reset
        for(List<HybridEnclosure>::Iterator jump_set_iter=jump_sets.begin(); jump_set_iter!=jump_sets.end(); ++jump_set_iter) {
            HybridEnclosure& _jump_set=*jump_set_iter;
            if(!definitely(_jump_set.is_empty())) {
                DiscreteLocation const& target=transitions[event].target;
                ARIADNE_LOG_PRINTLN_AT(2,"target="<<target<<", auxiliary_space="<<this->system().continuous_auxiliary_space(target)<<", auxiliary_function="<<this->system().auxiliary_function(target));
                TransitionData const& transition=transitions[event];
                _jump_set.apply_reset(event,target,transition.target_space,transition.reset_function,transition.target_auxiliary_space,transition.target_auxiliary_function);
//                _jump_set.set_auxiliary(this->system().continuous_auxiliary_space(target).variables(),this->system().auxiliary_function(target));
                evolution_data.initial_sets.append(_jump_set);
                _step_data.events.insert(event);
                ARIADNE_LOG_PRINTLN_AT(2,"jump_set="<<_jump_set);
            }
        }
    }

    if(_step_data.finishing || !_step_data.progress || !_step_data.events.empty()) {
        ARIADNE_LOG_PRINTLN_AT(2,timing_data.finishing_kind<<" "<<_step_data.events<<(_step_data.progress?" progress":" no progress"));
    }

}





TimingData
HybridEvolverBase::
_estimate_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                FlowFunctionModel const& flow,
                Map<DiscreteEvent,CrossingData>& crossings,
                Map<DiscreteEvent,TransitionData> const& transitions,
                HybridEnclosure const& initial_set) const
{
    // Compute the evolution time for the given step.
    ARIADNE_LOG_SCOPE_CREATE;
    const StepSizeType step_size=static_cast<StepSizeType>(flow.domain()[flow.domain().size()-1].upper_bound());
    TimingData result;
    result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
    result.finishing_kind=FinishingKind::STRADDLE_FINAL_TIME;
    result.step_size=step_size;
    result.final_time=final_time;
    result.evolution_time_domain=ExactIntervalType(0,step_size);
    result.evolution_time_coordinate=this->function_factory().create_identity(result.evolution_time_domain);
    result.parameter_dependent_evolution_time=this->function_factory().create_constant(initial_set.parameter_domain(),FloatDPValue(result.step_size));
    ARIADNE_LOG_PRINTLN("timing_data="<<result);
    return result;
}


HybridEvolverBaseConfiguration::HybridEvolverBaseConfiguration(HybridEvolverBase& evolver)
    : _evolver(evolver)
{
    set_flow_accuracy(1e-5);
    set_maximum_step_size(1.0);
    set_maximum_enclosure_radius(100.0);
    set_maximum_spacial_error(1e-2);
    set_enable_reconditioning(true);
    set_enable_subdivisions(true);
}

Void
HybridEvolverBaseConfiguration::set_flow_accuracy(const ApproximateRealType value)
{
    Configuration<TaylorPicardIntegrator> config;
    config.set_step_maximum_error(value);
    _evolver._integrator_ptr=std::shared_ptr<TaylorPicardIntegrator>(new TaylorPicardIntegrator(config));
    _flow_accuracy = cast_exact(value);
}


OutputStream&
HybridEvolverBaseConfiguration::_write(OutputStream& os) const
{
    os << "HybridEvolverBaseConfiguration"
       << "(\n  flow_accuracy=" << flow_accuracy()
       << ",\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  maximum_spacial_error=" << maximum_spacial_error()
       << ",\n  enable_reconditioning=" << enable_reconditioning()
       << ",\n  enable_subdivisions=" << enable_subdivisions()
       << "\n)";
    return os;
}

GeneralHybridEvolver::GeneralHybridEvolver(const SystemType& system)
    : HybridEvolverBase(system)
{
    this->_configuration_ptr.reset(new GeneralHybridEvolverConfiguration(*this));
}


GeneralHybridEvolver::GeneralHybridEvolver(
        const SystemType& system,
        const ValidatedFunctionModelDPFactoryInterface& factory)
    : HybridEvolverBase(system,factory)
{
    this->_configuration_ptr.reset(new GeneralHybridEvolverConfiguration(*this));
}


TimingData
GeneralHybridEvolver::
_estimate_timing(Set<DiscreteEvent>& active_events,
                Real final_time,
                FlowFunctionModel const& flow,
                Map<DiscreteEvent,CrossingData>& crossings,
                Map<DiscreteEvent,TransitionData> const& transitions,
                HybridEnclosure const& initial_set) const
{
    // Compute the evolution time for the given step.
    ARIADNE_LOG_SCOPE_CREATE;

    const SizeType n = flow.result_size();
    const FloatDPValue step_size=flow.domain()[flow.domain().size()-1].upper_bound();
    const FloatDPBounds final_time_bounds(final_time,dp);

    TimingData result;
    result.step_size=flow.step_size();
    result.final_time=final_time;

    ExactBoxType state_domain = cast_exact_box(initial_set.state_bounding_box());
    ExactIntervalType time_domain = cast_exact_interval(initial_set.time_range()+ExactIntervalType(zero,step_size));
    ExactBoxType statetime_domain = product(state_domain,time_domain);

    //ValidatedVectorMultivariateFunctionModelDP space_coordinates=this->function_factory().create_identity(space_domain);
    ValidatedScalarMultivariateFunctionModelDP time_coordinate=this->function_factory().create_coordinate(statetime_domain,n);
    ValidatedScalarMultivariateFunctionModelDP time_identity=this->function_factory().create_identity(time_domain);

    result.evolution_time_domain=ExactIntervalType(zero,step_size);
    result.evolution_time_coordinate=this->function_factory().create_identity(result.evolution_time_domain);

    ExactBoxType flow_state_domain = ExactBoxType(project(flow.domain(),range(0,n)));
    if(!subset(state_domain,flow_state_domain)) {
        ARIADNE_WARN(std::setprecision(17)<<"Bounding box "<<state_domain<<" is not subset of the flow spacial domain "<<flow_state_domain);
        state_domain=hull(state_domain,flow_state_domain);
    }

    // NOTE: The starting time function may be negative or greater than the final time
    // over part of the parameter domain.
    ValidatedVectorMultivariateFunctionModelDP const& starting_state_function=initial_set.state_function();
    ValidatedScalarMultivariateFunctionModelDP const& starting_time_function=initial_set.time_function();
    UpperIntervalType starting_time_range=initial_set.time_range();
    UpperIntervalType remaining_time_range=final_time_bounds-starting_time_range;

    ARIADNE_LOG_PRINTLN(std::fixed<<"starting_time_range="<<starting_time_range<<" step_size="<<step_size<<" final_time="<<final_time);


    // The time-dependent part of the evolution time
    ValidatedScalarMultivariateFunctionModelDP temporal_evolution_time=this->function_factory().create_zero(ExactIntervalVectorType(1u,time_domain));

    ARIADNE_LOG_PRINTLN(std::fixed<<"remaining_time_range="<<remaining_time_range);
    if(possibly(remaining_time_range.lower_bound()<zero)) {
        // Some of the points may already have reached the final time.
        // Don't try anything fancy, just do a simple constant time step.
        if(definitely(remaining_time_range.upper_bound()<=step_size)) {
            if(false) {
                // Within one time step we can go beyond final time
                result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
                result.finishing_kind=FinishingKind::AFTER_FINAL_TIME;
                temporal_evolution_time=FloatDPBounds(step_size); //   remaining_time_range.upper_bound();
            } else {
                result.step_kind=StepKind::CONSTANT_FINISHING_TIME;
                result.finishing_kind=FinishingKind::AT_FINAL_TIME;
                temporal_evolution_time=final_time_bounds-time_identity;
            }
        } else {
            result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
            result.finishing_kind=FinishingKind::STRADDLE_FINAL_TIME;
            temporal_evolution_time=FloatDPBounds(step_size);
        }
    } else if(definitely(remaining_time_range.upper_bound()<=result.step_size)) {
        // The rest of the evolution can be computed within a single time step.
        // The finishing kind is given as AT_FINAL_TIME so that the evolution algorithm
        // knows that the evolved set does not need to be evolved further.
        // This knowledge is required to be given combinarially, since
        // specifying the final time as a constant Function is not
        // exact if the final_time parameter is not exactly representable as
        // a FloatDPValue
        result.step_kind=StepKind::CONSTANT_FINISHING_TIME;
        result.finishing_kind=FinishingKind::AT_FINAL_TIME;
        temporal_evolution_time=final_time_bounds-time_identity;
    } else if(possibly(remaining_time_range.lower_bound()<=step_size) && ALLOW_CREEP) {
        // Some of the evolved points can be evolved to the final time in a single step
        // The evolution is performed over a step size which moves points closer to the final time, but does not cross.

        // Using the final_time as a guide, set the finishing time to closer to the final time.
        // This method ensures that points do not pass the final time after the transition.
        result.step_kind=StepKind::SPACETIME_DEPENDENT_FINISHING_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
        PositiveFloatDPValue sf={1u,dp};
        while(possibly(remaining_time_range.upper_bound()*sf>step_size)) { sf = hlf(sf); }
        temporal_evolution_time= FloatDPBounds(sf)*(final_time_bounds-time_identity);
    } else { // remaining_time_range.lower_bound()>step_size)
        // As far as timing goes, perform the evolution over a full time step
        result.step_kind=StepKind::CONSTANT_EVOLUTION_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
        temporal_evolution_time=FloatDPValue(result.step_size);
    }

    ARIADNE_LOG_PRINTLN_AT(1,"finishing_kind="<<result.finishing_kind);
    ARIADNE_LOG_PRINTLN_AT(1,"temporal_evolution_time="<<temporal_evolution_time);


    ValidatedScalarMultivariateFunctionModelDP spacial_evolution_time=this->function_factory().create_constant(state_domain,FloatDPValue(step_size));

    // Select one of GUARD_CREEP or TIME_CREEP
    static const Bool GUARD_CREEP=true;
    static const Bool TIME_CREEP=false;

    Bool creep=false;

    // Test for creep step
    if(ALLOW_CREEP && !crossings.empty() && (result.finishing_kind == FinishingKind::BEFORE_FINAL_TIME || result.finishing_kind == FinishingKind::STRADDLE_FINAL_TIME) ) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        ARIADNE_LOG_PRINTLN("Possible creep step; h="<<step_size);

        HybridEnclosure evolve_set=initial_set;
        evolve_set.apply_fixed_evolve_step(flow,flow.step_size());

        for(Map<DiscreteEvent,CrossingData>::Iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); ++crossing_iter)
        {
            DiscreteEvent event=crossing_iter->first;
            EventKind event_kind=transitions[event].event_kind;
            CrossingKind crossing_kind=crossing_iter->second.crossing_kind;
            ARIADNE_LOG_PRINTLN_AT(1,"Event "<<event<<": "<<event_kind<<": "<<crossing_kind);
            if(event_kind!=EventKind::PERMISSIVE) {
                evolve_set.new_invariant(event,transitions[event].guard_function);
            }
            // FIXME: When using permissive crossings, jumps in the step after the crossing time are lost.
            // A hack to fix this is to only creep on non-permissive events. Check that evolution is correct in this case.
            // FIXME: What should we do on increasing but non-transverse crossings?
            if((crossing_kind==CrossingKind::TRANSVERSE ) // || crossing_kind==CrossingKind::INCREASING)
                    && event_kind!=EventKind::PERMISSIVE)
            {
                ARIADNE_LOG_PRINTLN_AT(1,"crossing_time_range="<<crossing_iter->second.crossing_time.range());
                const ValidatedScalarMultivariateFunctionModelDP& crossing_time=crossing_iter->second.crossing_time;
                UpperIntervalType crossing_time_range=crossing_time.range();
                if(Ariadne::is_blocking(event_kind) && definitely(crossing_time_range.upper_bound()<step_size)) {
                    // NOTE: Use strict comparison here so that guard is fully crossed
                    // In principle this is not necessary, but this would involve testing
                    // to ensure that the evolved set is not propagated

                    // This event ensures that the evolve set is empty after a full step, so use this.
                    ARIADNE_LOG_PRINTLN_AT(1,std::setprecision(18)<<"crossing_time_range="<<crossing_time_range<<", crossing_time_range .upper_bound()="<<crossing_time_range.upper_bound()<<", step_size="<<step_size);
                    EffectiveScalarMultivariateFunction guard=transitions[event].guard_function;
                    ValidatedVectorMultivariateFunctionModelDP identity=this->function_factory().create_identity(crossing_time.domain());
                    ValidatedScalarMultivariateFunctionModelDP step_time=crossing_time*zero+FloatDPValue(step_size);
                    ARIADNE_LOG_PRINTLN_AT(1,"full flow="<<compose(flow,join(identity,step_time)));
                    ARIADNE_LOG_PRINTLN_AT(1,"guard range at crossing time="<<compose(guard,compose(flow,join(initial_set.state_function(),compose(crossing_time,initial_set.state_function())))).range());
                    ARIADNE_LOG_PRINTLN_AT(1,"guard range at crossing time="<<compose(guard,compose(flow,join(identity,crossing_time))).range());
                    ARIADNE_LOG_PRINTLN_AT(1,"No creep; event "<<crossing_iter->first<<" completely taken");
                    creep=false;
                    break;
                } else if(possibly(crossing_time_range.lower_bound()<=zero)) {
                    // This event is already partially active, so carry on with a full step
                    ARIADNE_LOG_PRINTLN_AT(1,"No creep; event "<<crossing_iter->first<<" already partially active");
                    creep=false;
                    break;
                } else if(definitely(crossing_time_range.lower_bound()>=step_size)) {
                    ARIADNE_LOG_PRINTLN_AT(1,"Event "<<crossing_iter->first<<" is not actually reached");
                    //crossings.erase(pending_erase_crossing_iter);
                } else {
                    ARIADNE_LOG_PRINTLN_AT(1,"Event "<<crossing_iter->first<<" can be crept up to.");
                    creep=true;
                }
            }
        }
        // If the evolved set is definitely empty, no creeping occurs
        if(definitely(evolve_set.is_empty())) {
            ARIADNE_LOG_PRINTLN_AT(1,"No creep; evolve set is empty");
            creep=false;
        }
    }

    if(ALLOW_CREEP && creep==false) { ARIADNE_LOG_PRINTLN("No creep"); }

    if(TIME_CREEP && creep==true) {
        // Compute reduced evolution time; note that for every remaining increasing crossing, we have
        // crossing_time_range.lower_bound()>zero and crossing_time_range.upper_bound()>step_size.

        for(Map<DiscreteEvent,CrossingData>::ConstIterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); ++crossing_iter)
        {
            if(crossing_iter->second.crossing_kind==CrossingKind::TRANSVERSE) {
                // Modify the crossing time function to be the smallest possible; this ensures that the evaluation time is
                // essentially exact
                ValidatedScalarMultivariateFunctionModelDP lower_crossing_time=crossing_iter->second.crossing_time;
                FloatDPError crossing_time_error=lower_crossing_time.error();
                lower_crossing_time.clobber();
                lower_crossing_time-=FloatDPValue(crossing_time_error.raw());

                // One possibility is to use quadratic restrictions
                //   If 0<=x<=2h, then x(1-x/4h)<=min(x,h)
                //   If 0<=x<=4h, then x(1-x/8h)<=min(x,2h)
                // Formula below works if x<=2h
                //   evolution_time=evolution_time*(crossing_time/result.step_size)*(1.0-crossing_time/(4*result.step_size));

                // Prefer simpler linear restrictions.
                // Multiply evolution time by crossing_time/max_crossing_time
                spacial_evolution_time=spacial_evolution_time*lower_crossing_time/cast_exact(lower_crossing_time.range().upper_bound());
            }
        }
        // Erase increasing transverse crossings since these cannot occur
        for(Map<DiscreteEvent,CrossingData>::Iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); )
        {
            if(crossing_iter->second.crossing_kind==CrossingKind::TRANSVERSE) {
                crossings.erase(crossing_iter++);
            } else {
                ++crossing_iter;
            }
        }
        ARIADNE_LOG_PRINTLN("Creep step: spacial_evolution_time="<<spacial_evolution_time);
        result.step_kind=StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
    }

    if(GUARD_CREEP && creep==true) {
        // If an event is possible, but only some points reach the guard set
        // after a full time step, then terminating the evolution here will
        // cause a splitting of the enclosure set. To prevent this, attempt
        // to "creep" up to the event guard boundary, so that in the next step,
        // all points can be made to cross.
        // Test to see if a creep step is needed
        ARIADNE_LOG_PRINTLN("Possible creep step; h="<<step_size);
        ARIADNE_LOG_PRINTLN("crossings="<<crossings);

        const SolverInterface& solver=*this->_solver_ptr;

        HybridEnclosure evolve_set=initial_set;
        evolve_set.apply_fixed_evolve_step(flow,flow.step_size());
        EffectiveVectorMultivariateFunction dynamic=this->system().dynamic_function(initial_set.location());
        ExactBoxType flow_spacial_domain=project(flow.domain(),range(0,flow.argument_size()-1u));
        ExactIntervalType flow_time_domain=flow.domain()[flow.argument_size()-1u];
        ValidatedScalarMultivariateFunctionModelDP zero_function=factory(flow).create_zero();
        ValidatedVectorMultivariateFunctionModelDP identity_function=factory(flow).create_identity();
        ValidatedVectorMultivariateFunctionModelDP space_projection=flow*zero;
        for(SizeType i=0; i!=n; ++i) { space_projection[i]=space_projection[i]+identity_function[i]; }

        //static const ExactDouble CREEP_MAXIMUM=1.0_X;
        static const ExactDouble CREEP_MAXIMUM=0.9375_x;
        spacial_evolution_time=this->function_factory().create_constant(flow.space_domain(),flow.step_size()*CREEP_MAXIMUM);

        for(Map<DiscreteEvent,CrossingData>::Iterator crossing_iter=crossings.begin();
            crossing_iter!=crossings.end(); )
        {
            DiscreteEvent event=crossing_iter->first;
            EventKind event_kind=transitions[event].event_kind;
            CrossingKind crossing_kind=crossing_iter->second.crossing_kind;
            EffectiveScalarMultivariateFunction guard_function=transitions[event].guard_function;
            ARIADNE_LOG_PRINTLN_AT(1,"Event "<<event<<": "<<event_kind<<": "<<crossing_kind);
            if(event_kind!=EventKind::PERMISSIVE) {
                UpperIntervalType guard_range = compose(guard_function,flow).range();
                ARIADNE_ASSERT(decide(guard_range.lower_bound()<zero));
                UpperIntervalType guard_derivative_range = compose(lie_derivative(guard_function,dynamic),flow).range();

                FloatDPBounds alpha_val=(1+flow.step_size()*cast_exact(guard_derivative_range.lower_bound())/cast_exact(guard_range.lower_bound()));
                FloatDPValue alpha=cast_exact(alpha_val);
                ARIADNE_ASSERT(alpha_val.value()==alpha);
                ARIADNE_LOG_PRINTLN_AT(1,"step_size: "<<flow.step_size()<<", guard_range: "<<guard_range<<", guard_derivative_range: "<<guard_derivative_range<<", alpha: "<<alpha);
                if(alpha>0 && alpha<=1) {
                    ValidatedScalarMultivariateFunctionModelDP guard_creep_time;
                    Bool successfully_computed_guard_creep_time=false;
                    try {
                        guard_creep_time=solver.implicit(compose(guard_function,flow)-alpha*compose(guard_function,space_projection),
                                                        flow_spacial_domain,flow_time_domain);
                        ARIADNE_LOG_PRINTLN_AT(1,"guard_creep_time= "<<guard_creep_time);
                        ARIADNE_LOG_PRINTLN_AT(1,"guard_creep_time.range()="<<guard_creep_time.range());
                        successfully_computed_guard_creep_time=true;
                        ARIADNE_LOG_PRINTLN_AT(2,"sucessfully_computed_guard_creep_time="<<successfully_computed_guard_creep_time);
                    }
                    catch(...) {
                        ARIADNE_LOG_PRINTLN_AT(1,"Error in computing guard creep time");
                    }
                    if(successfully_computed_guard_creep_time) {
                        spacial_evolution_time = spacial_evolution_time * (guard_creep_time/flow.step_size());
                        ARIADNE_LOG_PRINTLN_AT(2,"spacial_evolution_time="<<spacial_evolution_time);
                        ARIADNE_LOG_PRINTLN_AT(2,"crossings before erasing="<<crossings);
                        crossings.erase(crossing_iter++);
                        ARIADNE_LOG_PRINTLN_AT(2,"crossings after erasing="<<crossings);
                    } else {
                      ++crossing_iter;
                    }
                } else {
                  ++crossing_iter;
                }
            } else {
              ++crossing_iter;
            }
        }
        spacial_evolution_time.clobber();


        ARIADNE_LOG_PRINTLN_AT(1,"Creep step: spacial_evolution_time="<<spacial_evolution_time);
        ARIADNE_LOG_PRINTLN_AT(1,"spacial_evolution_time.range()="<<spacial_evolution_time.range());
        ARIADNE_LOG_PRINTLN_AT(1,"remaining crossings="<<crossings);
        result.step_kind=StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
    }


    ValidatedScalarMultivariateFunctionModelDP evolution_time = embed(spacial_evolution_time,time_domain) * embed(state_domain,temporal_evolution_time/FloatDPValue(step_size));
    ValidatedScalarMultivariateFunctionModelDP finishing_time=evolution_time+time_coordinate;

    ARIADNE_LOG_PRINTLN_AT(1,"evolution_time="<<(evolution_time));
    ARIADNE_LOG_PRINTLN_AT(1,"finishing_time="<<(finishing_time));
    result.spacetime_dependent_evolution_time=evolution_time;
    result.spacetime_dependent_finishing_time=finishing_time;
    result.parameter_dependent_evolution_time=unchecked_compose(evolution_time,join(starting_state_function,starting_time_function));
    result.parameter_dependent_finishing_time=unchecked_compose(finishing_time,join(starting_state_function,starting_time_function));

    // Test to see if it is possible to unwind crossings
    if(this->ALLOW_UNWIND && crossings.empty() && (starting_time_range.lower_bound().raw()<starting_time_range.upper_bound().raw())
            && result.step_kind==StepKind::CONSTANT_EVOLUTION_TIME && result.finishing_kind==FinishingKind::BEFORE_FINAL_TIME) {
        ARIADNE_LOG_PRINTLN_AT(1,"Possible unwind step; starting_time_range="<<starting_time_range<<", step_size="<<step_size);
        // Try to unwind the evolution time to a constant
        result.step_kind=StepKind::PARAMETER_DEPENDENT_FINISHING_TIME;
        result.finishing_kind=FinishingKind::BEFORE_FINAL_TIME;
        if(decide(starting_time_range.width()*2<step_size)) {
            result.parameter_dependent_finishing_time=this->function_factory().create_constant(initial_set.parameter_domain(),cast_exact(starting_time_range.lower_bound())+step_size);
        } else {
            // Try to reduce the time interval by half the step size
            // Corresponds to setting omega(smin)=tau(smin)+h, omega(smax)=tau(smax)+h/2
            // Taking omega(s)=a tau(s) + b, we obtain
            //   a=1-h/2(tmax-tmin);  b=h(tmax-tmin/2)/(tmax-tmin) = (2tmax-tmin)a
            FloatDPValue h=result.step_size;
            FloatDPValue tmin=cast_exact(starting_time_range.lower_bound());
            FloatDPValue tmax=cast_exact(starting_time_range.upper_bound());
            FloatDPBounds a=1-(hlf(h)/(tmax-tmin));
            FloatDPBounds b=h*(tmax-hlf(tmin))/(tmax-tmin);
            result.parameter_dependent_finishing_time=a*starting_time_function+b;
        }
        ARIADNE_LOG_PRINTLN_AT(1,"Unwinding to time "<<result.parameter_dependent_finishing_time);
        result.parameter_dependent_evolution_time=result.parameter_dependent_finishing_time-starting_time_function;
    }

    ARIADNE_LOG_PRINTLN_AT(1,"step_kind="<<result.step_kind<<", finishing_kind="<<result.finishing_kind);
    ARIADNE_LOG_PRINTLN_AT(1,"parameter_dependent_evolution_time="<<result.parameter_dependent_evolution_time);
    ARIADNE_LOG_PRINTLN_AT(1,"parameter_dependent_finishing_time="<<result.parameter_dependent_finishing_time<<"\n");
    return result;
}


GeneralHybridEvolverConfiguration::GeneralHybridEvolverConfiguration(GeneralHybridEvolver& evolver)
    : HybridEvolverBaseConfiguration(evolver)
{
}

GeneralHybridEvolverFactory::GeneralHybridEvolverFactory()
    : _function_factory(make_taylor_function_factory())
{
}

GeneralHybridEvolverFactory::GeneralHybridEvolverFactory(const ValidatedFunctionModelDPFactoryInterface& factory)
    : _function_factory(factory.clone())
{
}


GeneralHybridEvolver*
GeneralHybridEvolverFactory::create(const HybridAutomatonInterface& system) const
{
    return new GeneralHybridEvolver(system,*_function_factory);
}


} // namespace Ariadne
