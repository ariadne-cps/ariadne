/***************************************************************************
 *            hybrid_automaton.cc
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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

#include <map>

#include "macros.h"
#include "stlio.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "grid_set.h"

namespace Ariadne {


typedef uint DimensionType;

class HybridSet {};


AtomicDiscreteMode::
AtomicDiscreteMode(AtomicDiscreteLocation location,
                   const List<RealAlgebraicAssignment>& equations,
                   const List<RealDifferentialAssignment>& dynamic)
    : _location(location), _algebraic_assignments(equations), _differential_assignments(dynamic)
{
}

std::ostream&
AtomicDiscreteMode::write(std::ostream& os) const
{
    const AtomicDiscreteMode& mode=*this;
    os << "AtomicDiscreteMode( "
       << "location=" << mode._location;
    if(mode._algebraic_assignments.size()>0) {
        os << ", algebraic_equations="<<mode._algebraic_assignments; }
    if(mode._differential_assignments.size()>0) {
        os << ", differential_equations="<<mode._differential_assignments; }
    return os << " )";
}




AtomicDiscreteTransition::
AtomicDiscreteTransition(DiscreteEvent event,
                   const AtomicDiscreteMode& source_mode,
                   const AtomicDiscreteMode& target_mode,
                   const List<RealUpdateAssignment>& reset,
                   const ContinuousPredicate& guard)
    : _event(event), _source(source_mode.location()), _target(target_mode.location()),
      _guard_predicate(guard), _update_assignments(reset)
{
}

AtomicDiscreteTransition::
AtomicDiscreteTransition(DiscreteEvent event,
                         AtomicDiscreteLocation source,
                         AtomicDiscreteLocation target,
                         const List<RealUpdateAssignment>& reset,
                         const ContinuousPredicate& guard)
    : _event(event), _source(source), _target(target),
      _guard_predicate(guard), _update_assignments(reset)
{
}

std::ostream&
AtomicDiscreteTransition::write(std::ostream& os) const
{
    const AtomicDiscreteTransition& transition=*this;
    return os << "AtomicDiscreteTransition( "
              << "event=" << transition._event << ", "
              << "source=" << transition._source << ", "
              << "target=" << transition._target << ", "
              << "reset=" << transition._update_assignments << ", "
              << "guard=" << transition._guard_predicate << " )";
}




AtomicHybridAutomaton::~AtomicHybridAutomaton()
{
}

AtomicHybridAutomaton::AtomicHybridAutomaton()
{
}

AtomicHybridAutomaton::AtomicHybridAutomaton(const std::string& name)
    : _name(name)
{
}





const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(AtomicDiscreteLocation location,
                          const List<RealAlgebraicAssignment>& equations,
                          const List<RealDifferentialAssignment>& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }

    RealSpace state_space;
    RealSpace auxiliary_space;
    RealSpace input_space;
    Set<String> defined_variables;
    Set<String> argument_variables;

    // Compute the auxiliary variables ordered by the given equations
    for(uint i=0; i!=equations.size(); ++i) {
        if(defined_variables.contains(equations[i].lhs.name())) {
            std::cerr<<defined_variables,equations[i].lhs.name();
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_mode",
                          "Variable "<<equations[i].lhs<<" is defined twice by the algebraic equations "<<equations<<" for mode "<<location);
        }
        auxiliary_space.append(equations[i].lhs);
        argument_variables.adjoin(equations[i].rhs.arguments());
    }

    // Compute the state variables ordered by the given differential equations
    for(uint i=0; i!=dynamic.size(); ++i) {
        if(defined_variables.contains(dynamic[i].lhs.base().name())) {
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_mode",
                          "Variable "<<dynamic[i].lhs.base()<<" is defined by the differential equations "<<dynamic<<" for mode "<<location<<" is already defined");
        }
        state_space.append(dynamic[i].lhs.base());
        argument_variables.adjoin(dynamic[i].rhs.arguments());
    }

    // Compute the input variables
    for(Set<String>::const_iterator variable_iter=argument_variables.begin();
        variable_iter!=argument_variables.end(); ++variable_iter)
    {
        if(!defined_variables.contains(*variable_iter)) {
            input_space.append(RealVariable(*variable_iter));
        }
    }

    // TODO: Compute function
    AtomicDiscreteMode new_mode=AtomicDiscreteMode(location,equations,dynamic);

    this->_modes.insert(location,new_mode);
    return this->_modes.value(location);
}

const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(AtomicDiscreteLocation location,
                          const List<RealDifferentialAssignment>& dynamic)
{
    return this->new_mode(location,List<RealAlgebraicAssignment>(),dynamic);
}

const AtomicDiscreteMode&
AtomicHybridAutomaton::new_mode(AtomicDiscreteLocation location,
                          const List<RealAlgebraicAssignment>& equations)
{
    return this->new_mode(location,equations,List<RealDifferentialAssignment>());
}



const AtomicDiscreteMode&
AtomicHybridAutomaton::new_invariant(AtomicDiscreteLocation location,
                                     DiscreteEvent action,
                                     const ContinuousPredicate& constraint)
{
    if(!this->has_mode(location)) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_invariant",
                      "mode "<<location<<" is not a location of the automaton "<<*this);
    }
    AtomicDiscreteMode& mode=this->_modes.value(location);
    mode._invariant_predicates.insert(action,constraint);
    return mode;
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::new_transition(AtomicDiscreteLocation source,
                                      DiscreteEvent event,
                                      AtomicDiscreteLocation target,
                                      const List<RealUpdateAssignment>& reset,
                                      const ContinuousPredicate& guard)
{
    AtomicDiscreteMode& source_mode=const_cast<AtomicDiscreteMode&>(this->mode(source)); // Non-constant since we may wish to update input variables
    const AtomicDiscreteMode& target_mode=this->mode(target);
    Set<RealVariable> target_state_variables(this->state_variables(target));
    for(uint i=0; i!=reset.size(); ++i) {
        if(!target_state_variables.contains(reset[i].lhs.base())) {
            ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::new_transition",
                          "reset "<<reset<<" specifies variable "<<reset[i].lhs.base().name()<<
                          " which is not a state variable of the target location "<<target);
        }
    }

    AtomicDiscreteTransition new_transition=AtomicDiscreteTransition(event,source_mode,target_mode,reset,guard);
    source_mode._transitions.insert(event,new_transition);
    return source_mode._transitions.value(event);
}


const AtomicDiscreteTransition&
AtomicHybridAutomaton::new_transition(AtomicDiscreteLocation source,
                                      DiscreteEvent event,
                                      AtomicDiscreteLocation target,
                                      const ContinuousPredicate& guard)
{
    return this->new_transition(source,event,target,List<RealUpdateAssignment>(),guard);
}




void
AtomicHybridAutomaton::set_grid(AtomicDiscreteLocation location,
                                const Grid& grid)
{
    ARIADNE_DEPRECATED("AtomicHybridAutomaton::set_grid(...)","Use RealVariable::set_resolution(...) instead.");
    ARIADNE_NOT_IMPLEMENTED;
}

void
AtomicHybridAutomaton::set_grid(const Grid& grid)
{
    ARIADNE_DEPRECATED("AtomicHybridAutomaton::set_grid(...)","Use RealVariable::set_resolution(...) instead.");
    ARIADNE_NOT_IMPLEMENTED;
}

void
AtomicHybridAutomaton::set_grid(const HybridGrid& hgrid)
{
    ARIADNE_DEPRECATED("AtomicHybridAutomaton::set_grid(...)","Use RealVariable::set_resolution(...) instead.");
    ARIADNE_NOT_IMPLEMENTED;
}






const std::string&
AtomicHybridAutomaton::name() const
{
    return this->_name;
}

Set<AtomicDiscreteLocation>
AtomicHybridAutomaton::locations() const
{
    return this->_modes.keys();
}

Set<DiscreteEvent>
AtomicHybridAutomaton::events(AtomicDiscreteLocation q) const
{
    return this->mode(q)._transitions.keys();
}


bool
AtomicHybridAutomaton::has_mode(AtomicDiscreteLocation state) const
{
    return this->_modes.has_key(state);
}


bool
AtomicHybridAutomaton::has_invariant(AtomicDiscreteLocation source, DiscreteEvent action) const
{
   return this->_modes.has_key(source) && this->_modes[source]._invariant_predicates.has_key(action);
}

bool
AtomicHybridAutomaton::has_transition(AtomicDiscreteLocation source, DiscreteEvent event) const
{
   return this->_modes.has_key(source) && this->_modes[source]._transitions.has_key(event);
}




Set<AtomicDiscreteMode>
AtomicHybridAutomaton::modes() const
{
    return Set<AtomicDiscreteMode>(this->_modes.values());
}

AtomicDiscreteMode&
AtomicHybridAutomaton::mode(AtomicDiscreteLocation state)
{
    if(!this->_modes.has_key(state)) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::mode(AtomicDiscreteLocation)",
                      state<<"is not a location of the automaton with locations "<<this->locations());
    }
    return this->_modes.find(state)->second;
}

const AtomicDiscreteMode&
AtomicHybridAutomaton::mode(AtomicDiscreteLocation state) const
{
    if(!this->_modes.has_key(state)) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::mode(AtomicDiscreteLocation)",
                      state<<"is not a location of the automaton with locations "<<this->locations());
    }
    return this->_modes[state];
}


const ContinuousPredicate&
AtomicHybridAutomaton::invariant(AtomicDiscreteLocation state, DiscreteEvent action) const
{
    if(!this->has_invariant(state,action)) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::invariant(AtomicDiscreteLocation,DiscreteEvent)",
                      "The automaton "<<*this<<" has no invariant corresponding to "<<action<<" in location "<<state);
    }
    return this->_modes[state]._invariant_predicates[action];
}

const AtomicDiscreteTransition&
AtomicHybridAutomaton::transition(AtomicDiscreteLocation source, DiscreteEvent event) const
{
    if(!this->has_transition(source,event)) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::transition(AtomicDiscreteLocation,DiscreteEvent)",
                      "The automaton "<<*this<<" has no mode transition with source "<<source<<" and event "<<event);
    }
    return this->_modes[source]._transitions[event];
}

Set<AtomicDiscreteTransition>
AtomicHybridAutomaton::transitions(AtomicDiscreteLocation source) const
{
    if(!this->has_mode(source)) {
        ARIADNE_THROW(std::runtime_error,"AtomicHybridAutomaton::transitions(AtomicDiscreteLocation)",
                      "The automaton "<<*this<<" has no mode "<<source);
    }
    return Set<AtomicDiscreteTransition>(this->_modes[source]._transitions.values());
}





/*
Map<DiscreteEvent,ScalarFunction>
AtomicHybridAutomaton::blocking_guards(AtomicDiscreteLocation source) const
{
    std::map<DiscreteEvent,ScalarFunction> result;
    const AtomicDiscreteMode& mode=this->mode(source);
    for(invariant_const_iterator invariant_iter=mode._invariants.begin();
        invariant_iter!=mode._invariants.end(); ++invariant_iter)
    {
        const DiscreteEvent event=invariant_iter->first;
        const ScalarFunction invariant=invariant_iter->second;
        result[event]=invariant;
    }

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source().location()==source && transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const ScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}


std::map<DiscreteEvent,ScalarFunction>
AtomicHybridAutomaton::permissive_guards(AtomicDiscreteLocation source) const
{
    std::map<DiscreteEvent,ScalarFunction> result;

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source().location()==source && !transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const ScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}
*/


/*
Grid
AtomicHybridAutomaton::grid(AtomicDiscreteLocation location) const
{
    ARIADNE_ASSERT(this->has_mode(location));
    const AtomicDiscreteMode& mode=this->mode(location);
    return Grid(mode.dimension());
}

HybridGrid
AtomicHybridAutomaton::grid() const
{
    HybridGrid result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        result[mode_iter->location()]=mode_iter->grid();
    }
    return result;
}
*/


std::ostream&
AtomicHybridAutomaton::write(std::ostream& os) const
{
    const AtomicHybridAutomaton& ha=*this;
    os << "\nHybridAutomaton( \n  modes=\n";
    Set<AtomicDiscreteMode> modes=ha.modes();
    for(Set<AtomicDiscreteMode>::const_iterator mode_iter=modes.begin();
            mode_iter!=modes.end(); ++mode_iter)
    {
        os << "    " <<*mode_iter<<",\n";
    }
    os << "  transitions=\n";
    for(Set<AtomicDiscreteMode>::const_iterator mode_iter=modes.begin();
        mode_iter!=modes.end(); ++mode_iter)
    {
        Set<AtomicDiscreteTransition> transitions=ha.transitions(mode_iter->location());
        for(Set<AtomicDiscreteTransition>::const_iterator transition_iter=transitions.begin();
            transition_iter!=transitions.end(); ++transition_iter)
        {
            os << "    " <<*transition_iter<<",\n";
        }
    }
    return os << ")\n";
}



List<RealVariable>
AtomicHybridAutomaton::state_variables(AtomicDiscreteLocation location) const {
    List<RealVariable> result;
    const AtomicDiscreteMode& mode=this->mode(location);
    for(uint i=0; i!=mode._differential_assignments.size(); ++i) {
        result.append(mode._differential_assignments[i].lhs.base());
    }
    return result;
}

List<RealVariable>
AtomicHybridAutomaton::auxiliary_variables(AtomicDiscreteLocation location) const {
    List<RealVariable> result;
    const AtomicDiscreteMode& mode=this->mode(location);
    for(uint i=0; i!=mode._algebraic_assignments.size(); ++i) {
        result.append(mode._algebraic_assignments[i].lhs);
    }
    return result;
}


AtomicDiscreteLocation
AtomicHybridAutomaton::target(const AtomicDiscreteLocation& source, const DiscreteEvent& event) const {
    return this->transition(source,event).target();
}

List<RealAssignment>
AtomicHybridAutomaton::algebraic_assignments(const AtomicDiscreteLocation& location) const {
    return this->mode(location)._algebraic_assignments;
}

List<DottedRealAssignment>
AtomicHybridAutomaton::differential_assignments(const AtomicDiscreteLocation& location) const {
    return this->mode(location)._differential_assignments;
}

List<PrimedRealAssignment>
AtomicHybridAutomaton::update_assignments(const AtomicDiscreteLocation& source, const DiscreteEvent& event) const {
    if(this->has_transition(source,event)) {
        return this->transition(source,event)._update_assignments;
    } else {
        List<PrimedRealAssignment> nonjump_assignments;
        List<RealVariable> state_variables=this->state_variables(source);
        for(uint i=0; i!=state_variables.size(); ++i) {
            nonjump_assignments.append(next(state_variables[i])=state_variables[i]);
        }
        return nonjump_assignments;
    }
}

Map<DiscreteEvent,ContinuousPredicate>
AtomicHybridAutomaton::invariant_predicates(const AtomicDiscreteLocation& location) const {
    return this->mode(location)._invariant_predicates;
}

ContinuousPredicate
AtomicHybridAutomaton::invariant_predicate(const AtomicDiscreteLocation& location, const DiscreteEvent& action) const {
    if(this->has_invariant(location,action)) {
        return this->mode(location)._invariant_predicates[action];
    } else {
        return ContinuousPredicate(true);
    }
}

ContinuousPredicate
AtomicHybridAutomaton::guard_predicate(const AtomicDiscreteLocation& source, const DiscreteEvent& event) const {
    if(this->has_transition(source,event)) {
        return this->transition(source,event)._guard_predicate;
    } else {
        return ContinuousPredicate(true);
    }
}





CompositeHybridAutomaton::CompositeHybridAutomaton()
    : _components() { }

CompositeHybridAutomaton::CompositeHybridAutomaton(const AtomicHybridAutomaton& automaton)
    : _components(1u,automaton) { }

CompositeHybridAutomaton::CompositeHybridAutomaton(const List<AtomicHybridAutomaton>& components)
    : _components(components) { }

uint
CompositeHybridAutomaton::number_of_components() const
{
    return this->_components.size();
}

const AtomicHybridAutomaton&
CompositeHybridAutomaton::component(uint k) const {
    return this->_components[k];
}

List<DottedRealVariable> dot(const List<RealVariable>& v) {
    List<DottedRealVariable> result;
    for(uint i=0; i!=v.size(); ++i) { result.append(dot(v[i])); }
    return result;
}

List<PrimedRealVariable> next(const List<RealVariable>& v) {
    List<PrimedRealVariable> result;
    for(uint i=0; i!=v.size(); ++i) { result.append(next(v[i])); }
    return result;
}



bool
CompositeHybridAutomaton::has_mode(const DiscreteLocation& location) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(!this->_components[i].has_mode(location[i])) {
            return false;
        }
    }
    return true;
}

bool
CompositeHybridAutomaton::has_transition(const DiscreteLocation& source, const DiscreteEvent& event) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_transition(source[i],event)) {
            return true;
        }
    }
    return false;
}

Set< DiscreteEvent >
CompositeHybridAutomaton::events(const DiscreteLocation& source) const
{
    Set< DiscreteEvent > result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.adjoin(this->_components[i].events(source[i]));
    }
    return result;
}



DiscreteLocation
CompositeHybridAutomaton::target(const DiscreteLocation& source, const DiscreteEvent& event) const {
   // std::cerr<<"CompositeHybridAutomaton::target("<<source<<","<<event<<")\n";
    DiscreteLocation result; result.reserve(source.size());
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_transition(source[i],event)) {
            AtomicDiscreteLocation target=this->_components[i].target(source[i],event);
            result.append(target);
        } else {
            result.append(source[i]);
        }
    }
    return result;
}


List<RealVariable>
CompositeHybridAutomaton::variables(const DiscreteLocation& location) const {
    return catenate(this->state_variables(location),this->auxiliary_variables(location));
}

List<RealVariable>
CompositeHybridAutomaton::state_variables(const DiscreteLocation& location) const {
    List<RealVariable> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location[i]));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].state_variables(location[i]));
    }
    return result;
}

List<RealVariable>
CompositeHybridAutomaton::auxiliary_variables(const DiscreteLocation& location) const {
    List<RealVariable> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location[i]));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].auxiliary_variables(location[i]));
    }
    return result;
}

// Find all algebraic equations valid in the location
List<RealAssignment>
CompositeHybridAutomaton::algebraic_assignments(const DiscreteLocation& location) const {
    List<RealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].algebraic_assignments(location[i]));
    }
    // TODO: Sort the result to eliminate algebraic loops
    // sort(result);
    return result;
}

List<DottedRealAssignment>
CompositeHybridAutomaton::differential_assignments(const DiscreteLocation& location) const {
    List<DottedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].differential_assignments(location[i]));
    }
    // No need to sort result since dotted variables cannot appear in the right-hand side (currently)
    return result;
}


List<PrimedRealAssignment>
CompositeHybridAutomaton::update_assignments(const DiscreteLocation& location, const DiscreteEvent& event) const {
    List<PrimedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].update_assignments(location[i],event));
    }
    return result;
}

Map<DiscreteEvent,ContinuousPredicate>
CompositeHybridAutomaton::invariant_predicates(const DiscreteLocation& location) const {
    Map<DiscreteEvent,ContinuousPredicate> result;
    Set<DiscreteEvent> actions;
    for(uint i=0; i!=this->_components.size(); ++i) {
        Map<DiscreteEvent,ContinuousPredicate> invariants=this->_components[i].invariant_predicates(location[i]);
        actions.adjoin(invariants.keys());
    }
    for(Set<DiscreteEvent>::const_iterator action_iter=actions.begin(); action_iter!=actions.end(); ++action_iter) {
        result.insert(*action_iter,this->invariant_predicate(location,*action_iter));
    }
    return result;
}

ContinuousPredicate
CompositeHybridAutomaton::invariant_predicate(const DiscreteLocation& location, const DiscreteEvent& event) const {
    ContinuousPredicate result(true);
    for(uint i=0; i!=this->_components.size(); ++i) {
        ContinuousPredicate guard=this->_components[i].invariant_predicate(location[i],event);
        if(result==true || guard==false) { result=guard; }
        else if(guard==true || result==false) { }
        else { result = result && guard; }
    }
    return result;
}


ContinuousPredicate
CompositeHybridAutomaton::guard_predicate(const DiscreteLocation& location, const DiscreteEvent& event) const {
    ContinuousPredicate result(true);
    for(uint i=0; i!=this->_components.size(); ++i) {
        ContinuousPredicate guard=this->_components[i].guard_predicate(location[i],event);
        if(result==true || guard==false) { result=guard; }
        else if(guard==true || result==false) { }
        else { result = result && guard; }
    }
    return result;
}



void
CompositeHybridAutomaton::check_mode(const DiscreteLocation& location) const {
    List<RealAssignment> equations=this->algebraic_assignments(location);
    List<DottedRealAssignment> dynamics=this->differential_assignments(location);

    ARIADNE_NOT_IMPLEMENTED;
}


VectorFunction
CompositeHybridAutomaton::output_function(const DiscreteLocation& location) const {
    return VectorFunction(auxiliary_variables(location),algebraic_assignments(location),state_variables(location));
}

VectorFunction
CompositeHybridAutomaton::dynamic_function(const DiscreteLocation& location) const {
    return VectorFunction(dot(state_variables(location)),differential_assignments(location),variables(location));
}

VectorFunction
CompositeHybridAutomaton::reset_function(const DiscreteLocation& source, const DiscreteEvent& event) const {
    DiscreteLocation target=this->target(source,event);
    return VectorFunction(next(state_variables(target)),update_assignments(source,event),variables(source));
}

ScalarFunction
CompositeHybridAutomaton::invariant_function(const DiscreteLocation& location, const DiscreteEvent& event) const {
    return ScalarFunction(indicator(invariant_predicate(location,event),negative),variables(location));
}

ScalarFunction
CompositeHybridAutomaton::guard_function(const DiscreteLocation& location, const DiscreteEvent& event) const {
    return ScalarFunction(indicator(guard_predicate(location,event)),variables(location));
}

Map<DiscreteEvent,ScalarFunction>
CompositeHybridAutomaton::invariant_functions(const DiscreteLocation& location) const {
    ARIADNE_NOT_IMPLEMENTED;
}

Map<DiscreteEvent,ScalarFunction>
CompositeHybridAutomaton::guard_functions(const DiscreteLocation& location) const {
    ARIADNE_NOT_IMPLEMENTED;
}



std::ostream&
CompositeHybridAutomaton::write(std::ostream& os) const
{
    return os << "CompositeHybridAutomaton(\n" << this->_components << "\n)\n";
}





Set<DiscreteLocation>
discrete_reachability(const CompositeHybridAutomaton& automaton, const DiscreteLocation& initial_location)
{
    typedef uint Nat;

    Set<DiscreteLocation> reached;
    Set<DiscreteLocation> working;
    Set<DiscreteLocation> found;

    Map<DiscreteLocation,Nat> steps;
    Nat step=0u;

    std::cerr<<"\n\n";
    working.insert(initial_location);
    reached.insert(initial_location);
    steps.insert(initial_location,step);
    while(!working.empty()) {
        ++step;
        for(Set<DiscreteLocation>::const_iterator source_iter=working.begin(); source_iter!=working.end(); ++source_iter) {
            std::cerr<<"new_mode\n";
            DiscreteLocation location=*source_iter;
            std::cerr<<"  mode: "<<location<<":\n";
            std::cerr<<"      output="<<automaton.algebraic_assignments(location)<<"\n";
            std::cerr<<"          function="<<automaton.output_function(location)<<"\n";
            std::cerr<<"      dynamic="<<automaton.differential_assignments(location)<<"\n";
            std::cerr<<"          function="<<automaton.dynamic_function(location)<<"\n";

            Map<DiscreteEvent,ContinuousPredicate> invariants=automaton.invariant_predicates(location);
            for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator invariant_iter=invariants.begin();
                    invariant_iter!=invariants.end(); ++invariant_iter) {
                std::cerr<<"    invariant: "<<invariant_iter->second<<"\n";
                std::cerr<<"      function: "<<automaton.invariant_function(location,invariant_iter->first)<<"\n";
                std::cerr<<"        action: "<<invariant_iter->first<<"\n";
            }

            Set<DiscreteEvent> events=automaton.events(location);
            std::cerr<<"  events: "<<events<<"\n";
            for(Set<DiscreteEvent>::const_iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
                DiscreteEvent event=*event_iter;
                std::cerr<<"    event:"<<event<<"\n";
                DiscreteLocation target=automaton.target(location,event);
                std::cerr<<"    transition: "<<event<<" -> "<<target<<"\n";
                std::cerr<<"        reset="<<automaton.update_assignments(location,event)<<"\n";
                std::cerr<<"          function="<<automaton.reset_function(location,event)<<"\n";
                std::cerr<<"        guard="<<automaton.guard_predicate(location,event)<<"\n";
                std::cerr<<"          function="<<automaton.guard_function(location,event)<<"\n";
                if(!reached.contains(target)) {
                    found.insert(target);
                    reached.insert(target);
                    steps.insert(target,step);
                }
           }
        }
        std::cerr<<"\nstep "<<step<<" found: "<<found<<"\n\n";
        working.clear();
        working.swap(found);
    }

    return reached;
}

}