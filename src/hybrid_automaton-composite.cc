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
#include "expression.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "grid_set.h"

#include "hybrid_automaton-composite.h"

namespace Ariadne {


typedef uint DimensionType;

class HybridSet {};

class SystemSpecificationError : public std::runtime_error {
  public:
    SystemSpecificationError(const std::string& what) : std::runtime_error(what) { }
};

AtomicDiscreteMode::
AtomicDiscreteMode(AtomicDiscreteLocation location,
                   const List<RealAlgebraicAssignment>& equations,
                   const List<RealDifferentialAssignment>& dynamic)
    : _location(location), _algebraic_assignments(equations), _differential_assignments(dynamic)
{
}



template<class K, class T1, class T2>
Map<K,Pair<T1,T2> > merge(const Map<K,T1>& m1, const Map<K,T2>& m2) {
    assert(m1.keys()==m2.keys());
    Map<K,Pair<T1,T2> > r;
    Set<K> keys=m1.keys();
    for(typename Set<K>::const_iterator iter=keys.begin(); iter!=keys.end(); ++iter) {
        r.insert(*iter,make_pair(m1[*iter],m2[*iter]));
    }
    return r;
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
    if(mode._invariant_predicates.size()>0) {
        os << ", invariants="<<mode._invariant_predicates; }
    if(mode._urgent_guard_predicates.size()>0) {
        os << ", guards="<<mode._urgent_guard_predicates; }
    if(mode._permissive_guard_predicates.size()>0) {
        os << ", activations="<<mode._permissive_guard_predicates; }
    if(mode._targets.size()>0) {
        os << ", transitions="<<merge(mode._targets,mode._update_assignments); }
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





void
AtomicHybridAutomaton::new_mode(AtomicDiscreteLocation location,
                          const List<RealAlgebraicAssignment>& equations,
                          const List<RealDifferentialAssignment>& dynamic)
{
    if(this->has_mode(location)) {
        throw SystemSpecificationError("The hybrid automaton already has a mode with the given id");
    }

    Set<UntypedVariable> defined_variables;
    Set<UntypedVariable> argument_variables;

    // Compute the auxiliary variables ordered by the given equations
    for(uint i=0; i!=equations.size(); ++i) {
        if(defined_variables.contains(equations[i].lhs)) {
            ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_mode",
                          "Variable "<<equations[i].lhs<<" is defined twice by the algebraic equations "<<equations<<" for mode "<<location);
        }
    }

    // Compute the state variables ordered by the given differential equations
    for(uint i=0; i!=dynamic.size(); ++i) {
        if(defined_variables.contains(dynamic[i].lhs.base())) {
            ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_mode",
                          "Variable "<<dynamic[i].lhs.base()<<" is defined by the differential equations "<<dynamic<<" for mode "<<location<<" is already defined");
        }
        argument_variables.adjoin(dynamic[i].rhs.arguments());
    }

    // TODO: Compute function
    AtomicDiscreteMode new_mode=AtomicDiscreteMode(location,equations,dynamic);

    this->_modes.insert(location,new_mode);
}

void
AtomicHybridAutomaton::new_mode(AtomicDiscreteLocation location,
                          const List<RealDifferentialAssignment>& dynamic)
{
    return this->new_mode(location,List<RealAlgebraicAssignment>(),dynamic);
}

void
AtomicHybridAutomaton::new_mode(AtomicDiscreteLocation location,
                          const List<RealAlgebraicAssignment>& equations)
{
    return this->new_mode(location,equations,List<RealDifferentialAssignment>());
}



void
AtomicHybridAutomaton::new_invariant(AtomicDiscreteLocation location,
                                     DiscreteEvent action,
                                     const ContinuousPredicate& constraint)
{
    if(!this->has_mode(location)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_invariant",
                      "mode "<<location<<" is not a location of the automaton with locations "<<this->locations());
    }
    AtomicDiscreteMode& mode=this->_modes.value(location);
    mode._invariant_predicates.insert(action,constraint);
}


void
AtomicHybridAutomaton::new_urgent_guard(AtomicDiscreteLocation location,
                                        DiscreteEvent event,
                                        const ContinuousPredicate& guard)
{
    if(!this->has_mode(location)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_urgent_guard",
                      "Source mode "<<location<<" is not a location of the automaton with locations "<<this->locations());
    }


    AtomicDiscreteMode& mode=this->mode(location); // Get non-constant mode
    mode._urgent_guard_predicates.insert(event,guard);;
}


void
AtomicHybridAutomaton::new_permissive_guard(AtomicDiscreteLocation location,
                                        DiscreteEvent event,
                                        const ContinuousPredicate& guard)
{
    if(!this->has_mode(location)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_permissive_guard",
                      "Source mode "<<location<<" is not a location of the automaton with locations "<<this->locations());
    }

    AtomicDiscreteMode& mode=this->mode(location); // Get non-constant mode
    mode._permissive_guard_predicates.insert(event,guard);
}




void
AtomicHybridAutomaton::new_transition(AtomicDiscreteLocation source,
                                      DiscreteEvent event,
                                      AtomicDiscreteLocation target,
                                      const List<RealUpdateAssignment>& reset)
{
    if(!this->has_mode(source)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_reset",
                      "Source mode "<<source<<" is not a location of the automaton with locations "<<this->locations());
    }

    if(!this->has_mode(target)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_reset",
                      "Target mode "<<target<<" is not a location of the automaton with locations "<<this->locations());
    }

    AtomicDiscreteMode& source_mode=this->mode(source); // Non-constant since we may wish to update input variables
    Set<RealVariable> target_state_variables(this->state_variables(target));

    for(uint i=0; i!=reset.size(); ++i) {
        if(!target_state_variables.contains(reset[i].lhs.base())) {
            ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::new_transition",
                          "reset "<<reset<<" for event "<<event<<" in source "<<source<<
                          " specifies variable "<<reset[i].lhs.base().name()<<
                          " which is not a state variable "<<this->state_variables(target)<<" of the target location "<<target);
        }
    }

    source_mode._targets.insert(event,target);
    source_mode._update_assignments.insert(event,reset);
}


void
AtomicHybridAutomaton::new_reset(AtomicDiscreteLocation source,
                                 DiscreteEvent event,
                                 AtomicDiscreteLocation target,
                                 const List<RealUpdateAssignment>& reset)
{
    this->new_transition(source,event,target,reset);
}


void
AtomicHybridAutomaton::new_transition(AtomicDiscreteLocation source,
                                      DiscreteEvent event,
                                      AtomicDiscreteLocation target,
                                      const List<RealUpdateAssignment>& reset,
                                      const ContinuousPredicate& guard,
                                      EventKind kind)
{
    switch(kind) {
        case urgent: this->new_urgent_guard(source,event,guard); break;
        case permissive: this->new_permissive_guard(source,event,guard); break;
        default: assert(false);
    }

    this->new_reset(source,event,target,reset);
}


void
AtomicHybridAutomaton::new_transition(AtomicDiscreteLocation source,
                                      DiscreteEvent event,
                                      const ContinuousPredicate& guard,
                                      AtomicDiscreteLocation target,
                                      const List<RealUpdateAssignment>& reset,
                                      EventKind kind)
{
    return this->new_transition(source,event,target,reset,guard);
}



void
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
AtomicHybridAutomaton::transition_events(AtomicDiscreteLocation q) const
{
    return this->mode(q)._targets.keys();
}

Set<DiscreteEvent>
AtomicHybridAutomaton::invariant_events(AtomicDiscreteLocation q) const
{
    return join(this->mode(q)._invariant_predicates.keys(),this->mode(q)._urgent_guard_predicates.keys());
}

Set<DiscreteEvent>
AtomicHybridAutomaton::blocking_events(AtomicDiscreteLocation q) const
{
    return this->mode(q)._invariant_predicates.keys();
}

Set<DiscreteEvent>
AtomicHybridAutomaton::urgent_events(AtomicDiscreteLocation q) const
{
    return this->mode(q)._urgent_guard_predicates.keys();
}

Set<DiscreteEvent>
AtomicHybridAutomaton::permissive_events(AtomicDiscreteLocation q) const
{
    return this->mode(q)._permissive_guard_predicates.keys();
}


bool
AtomicHybridAutomaton::has_mode(AtomicDiscreteLocation state) const
{
    return this->_modes.has_key(state);
}


bool
AtomicHybridAutomaton::has_transition(AtomicDiscreteLocation source, DiscreteEvent event) const
{
   return this->_modes.has_key(source) && this->_modes[source]._targets.has_key(event);
}

bool
AtomicHybridAutomaton::has_invariant(AtomicDiscreteLocation source, DiscreteEvent event) const
{
   return this->_modes.has_key(source) && this->_modes[source]._invariant_predicates.has_key(event);
}


bool
AtomicHybridAutomaton::has_guard(AtomicDiscreteLocation source, DiscreteEvent event) const
{
   return this->has_invariant(source,event) || this->has_transition(source,event);
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
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::mode(AtomicDiscreteLocation)",
                      state<<"is not a location of the automaton with locations "<<this->locations());
    }
    return this->_modes.find(state)->second;
}

const AtomicDiscreteMode&
AtomicHybridAutomaton::mode(AtomicDiscreteLocation state) const
{
    if(!this->_modes.has_key(state)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::mode(AtomicDiscreteLocation)",
                      state<<"is not a location of the automaton with locations "<<this->locations());
    }
    return this->_modes[state];
}




const ContinuousPredicate&
AtomicHybridAutomaton::invariant(AtomicDiscreteLocation state, DiscreteEvent action) const
{
    if(!this->has_invariant(state,action)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::invariant(AtomicDiscreteLocation,DiscreteEvent)",
                      "The automaton "<<*this<<" has no invariant corresponding to "<<action<<" in location "<<state);
    }
    return this->_modes[state]._invariant_predicates[action];
}

const ContinuousPredicate&
AtomicHybridAutomaton::guard(AtomicDiscreteLocation state, DiscreteEvent action) const
{
    if(!this->has_invariant(state,action)) {
        ARIADNE_THROW(SystemSpecificationError,"AtomicHybridAutomaton::invariant(AtomicDiscreteLocation,DiscreteEvent)",
                      "The automaton "<<*this<<" has no invariant corresponding to "<<action<<" in location "<<state);
    }
    return this->_modes[state]._permissive_guard_predicates[action];
}




/*
Map<DiscreteEvent,RealScalarFunction>
AtomicHybridAutomaton::blocking_guards(AtomicDiscreteLocation source) const
{
    std::map<DiscreteEvent,RealScalarFunction> result;
    const AtomicDiscreteMode& mode=this->mode(source);
    for(invariant_const_iterator invariant_iter=mode._invariants.begin();
        invariant_iter!=mode._invariants.end(); ++invariant_iter)
    {
        const DiscreteEvent event=invariant_iter->first;
        const RealScalarFunction invariant=invariant_iter->second;
        result[event]=invariant;
    }

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source().location()==source && transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const RealScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}


std::map<DiscreteEvent,RealScalarFunction>
AtomicHybridAutomaton::permissive_guards(AtomicDiscreteLocation source) const
{
    std::map<DiscreteEvent,RealScalarFunction> result;

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source().location()==source && !transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const RealScalarFunction guard=transition_iter->activation();
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
AtomicHybridAutomaton::target(AtomicDiscreteLocation source, DiscreteEvent event) const {
    if(this->has_transition(source,event)) {
        return this->_modes[source]._targets[event];
    } else {
        return source;
    }
}

List<RealAssignment>
AtomicHybridAutomaton::algebraic_assignments(AtomicDiscreteLocation location) const {
    return this->mode(location)._algebraic_assignments;
}

List<DottedRealAssignment>
AtomicHybridAutomaton::differential_assignments(AtomicDiscreteLocation location) const {
    return this->mode(location)._differential_assignments;
}

List<PrimedRealAssignment>
AtomicHybridAutomaton::update_assignments(AtomicDiscreteLocation source, DiscreteEvent event) const {
    if(this->has_transition(source,event)) {
        return this->_modes[source]._update_assignments[event];
    } else {
        List<PrimedRealAssignment> nonjump_assignments;
        List<RealVariable> state_variables=this->state_variables(source);
        for(uint i=0; i!=state_variables.size(); ++i) {
            nonjump_assignments.append(next(state_variables[i])=state_variables[i]);
        }
        return nonjump_assignments;
    }
}

ContinuousPredicate
AtomicHybridAutomaton::invariant_predicate(AtomicDiscreteLocation location, DiscreteEvent action) const {
    const AtomicDiscreteMode& mode=this->mode(location);
    if(mode._invariant_predicates.has_key(action)) {
        return mode._invariant_predicates[action];
    } else if(mode._urgent_guard_predicates.has_key(action)) {
        return !mode._urgent_guard_predicates[action];
    } else {
        return ContinuousPredicate(true);
    }
}

ContinuousPredicate
AtomicHybridAutomaton::guard_predicate(AtomicDiscreteLocation source, DiscreteEvent event) const {
    const AtomicDiscreteMode& mode=this->mode(source);
    if(mode._permissive_guard_predicates.has_key(event)) {
        return mode._permissive_guard_predicates[event];
    } else if(mode._urgent_guard_predicates.has_key(event)) {
        return mode._urgent_guard_predicates[event];
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
CompositeHybridAutomaton::has_mode(DiscreteLocation location) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(!this->_components[i].has_mode(location[i])) {
            return false;
        }
    }
    return true;
}

bool
CompositeHybridAutomaton::has_guard(DiscreteLocation location, DiscreteEvent event) const
{
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_guard(location[i],event)) {
            return true;
        }
    }
    return false;
}

bool
CompositeHybridAutomaton::has_transition(DiscreteLocation source, DiscreteEvent event) const {
    for(uint i=0; i!=this->_components.size(); ++i) {
        if(this->_components[i].has_transition(source[i],event)) {
            return true;
        }
    }
    return false;
}


Set< DiscreteEvent >
CompositeHybridAutomaton::transition_events(DiscreteLocation source) const
{
    Set< DiscreteEvent > result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.adjoin(this->_components[i].transition_events(source[i]));
    }
    return result;
}

Set< DiscreteEvent >
CompositeHybridAutomaton::invariant_events(DiscreteLocation source) const
{
    Set< DiscreteEvent > result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.adjoin(this->_components[i].invariant_events(source[i]));
    }
    return result;
}

Set<DiscreteEvent>
CompositeHybridAutomaton::events(DiscreteLocation location) const
{
    return join(this->invariant_events(location),this->transition_events(location));
}

EventKind
CompositeHybridAutomaton::event_kind(DiscreteLocation location, DiscreteEvent event) const
{
    bool is_invariant=this->invariant_events(location).contains(event);
    bool is_activation=this->transition_events(location).contains(event);
    if(is_activation  && is_invariant) { return urgent; }
    else if(is_activation) { return permissive; }
    else { return PROGRESS; }
}




DiscreteLocation
CompositeHybridAutomaton::target(DiscreteLocation source, DiscreteEvent event) const {
    DiscreteLocation result;
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


uint
CompositeHybridAutomaton::dimension(DiscreteLocation location) const {
    return this->state_variables(location).size();
}

List<RealVariable>
CompositeHybridAutomaton::variables(DiscreteLocation location) const {
    return catenate(this->state_variables(location),this->auxiliary_variables(location));
}

List<RealVariable>
CompositeHybridAutomaton::state_variables(DiscreteLocation location) const {
    List<RealVariable> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        ARIADNE_ASSERT(this->_components[i].has_mode(location[i]));
        // Append state space in mode i; note that this automatically checks whether a variable is already present
        result.append(this->_components[i].state_variables(location[i]));
    }
    return result;
}

List<RealVariable>
CompositeHybridAutomaton::auxiliary_variables(DiscreteLocation location) const {
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
CompositeHybridAutomaton::algebraic_assignments(DiscreteLocation location) const {
    List<RealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].algebraic_assignments(location[i]));
    }
    // TODO: Sort the result to eliminate algebraic loops
    // sort(result);
    return result;
}

List<DottedRealAssignment>
CompositeHybridAutomaton::differential_assignments(DiscreteLocation location) const {
    List<DottedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].differential_assignments(location[i]));
    }
    // No need to sort result since dotted variables cannot appear in the right-hand side (currently)
    return result;
}


List<PrimedRealAssignment>
CompositeHybridAutomaton::update_assignments(DiscreteLocation location, DiscreteEvent event) const {
    List<PrimedRealAssignment> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.append(this->_components[i].update_assignments(location[i],event));
    }
    return result;
}


ContinuousPredicate
CompositeHybridAutomaton::invariant_predicate(DiscreteLocation location, DiscreteEvent event) const {
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
CompositeHybridAutomaton::guard_predicate(DiscreteLocation location, DiscreteEvent event) const {
    ContinuousPredicate result(true);
    for(uint i=0; i!=this->_components.size(); ++i) {
        ContinuousPredicate guard=this->_components[i].guard_predicate(location[i],event);
        if(result==true || guard==false) { result=guard; }
        else if(guard==true || result==false) { }
        else { result = result && guard; }
    }
    return result;
}





RealVectorFunction
CompositeHybridAutomaton::output_function(DiscreteLocation location) const {
    return RealVectorFunction(auxiliary_variables(location),algebraic_assignments(location),state_variables(location));
}

RealVectorFunction
CompositeHybridAutomaton::dynamic_function(DiscreteLocation location) const {
    List<RealDifferentialAssignment> differential=this->differential_assignments(location);
    List<RealAssignment> algebraic=this->algebraic_assignments(location);
    for(uint i=0; i!=algebraic.size(); ++i) {
        for(uint j=0; j!=differential.size(); ++j) {
            differential[j].rhs=substitute(differential[j].rhs,algebraic[i].lhs,algebraic[i].rhs);
        }
    }
    return RealVectorFunction(dot(state_variables(location)),differential,state_variables(location));
}

RealVectorFunction
CompositeHybridAutomaton::reset_function(DiscreteLocation source, DiscreteEvent event) const {
    DiscreteLocation target=this->target(source,event);
    List<RealUpdateAssignment> update=this->update_assignments(source,event);
    List<RealAssignment> algebraic=this->algebraic_assignments(source);
    for(uint i=0; i!=algebraic.size(); ++i) {
        for(uint j=0; j!=update.size(); ++j) {
            update[j].rhs=substitute(update[j].rhs,algebraic[i].lhs,algebraic[i].rhs);
        }
    }
    return RealVectorFunction(next(state_variables(target)),update,state_variables(source));
}

RealScalarFunction
CompositeHybridAutomaton::constraint_function(DiscreteLocation location, DiscreteEvent event) const {
    ARIADNE_NOT_IMPLEMENTED;
}

RealScalarFunction
CompositeHybridAutomaton::invariant_function(DiscreteLocation location, DiscreteEvent event) const {
    return RealScalarFunction(indicator(substitute(invariant_predicate(location,event),algebraic_assignments(location)),NEGATIVE),state_variables(location));
}

RealScalarFunction
CompositeHybridAutomaton::guard_function(DiscreteLocation location, DiscreteEvent event) const {
    return RealScalarFunction(indicator(substitute(guard_predicate(location,event),algebraic_assignments(location))),state_variables(location));
}



Set<DiscreteEvent>
CompositeHybridAutomaton::blocking_events(DiscreteLocation location) const
{
    Set<DiscreteEvent> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.adjoin(this->_components[i].blocking_events(location[i]));
    }
    return result;
}

Set<DiscreteEvent>
CompositeHybridAutomaton::urgent_events(DiscreteLocation location) const
{
    Set<DiscreteEvent> result;
    Set<DiscreteEvent> duplicates;
    for(uint i=0; i!=this->_components.size(); ++i) {
        Set<DiscreteEvent> found=this->_components[i].urgent_events(location[i]);
        duplicates.adjoin(intersection(result,found));
        result.adjoin(found);
    }
    if(!duplicates.empty()) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::urgent_events",
                      "Events "<<duplicates<<" are urgent in more than one component in location "<<location);
    }
    return result;
}


Set<DiscreteEvent>
CompositeHybridAutomaton::permissive_events(DiscreteLocation location) const
{
    Set<DiscreteEvent> result;
    for(uint i=0; i!=this->_components.size(); ++i) {
        result.adjoin(this->_components[i].permissive_events(location[i]));
    }
    return result;
}


Grid
CompositeHybridAutomaton::grid(DiscreteLocation location) const
{
    return Grid(this->state_variables(location).size());
}

HybridGrid
CompositeHybridAutomaton::grid() const
{
    return HybridGrid(*this);
}

std::ostream&
CompositeHybridAutomaton::write(std::ostream& os) const
{
    return os << "CompositeHybridAutomaton(\n" << this->_components << "\n)\n";
}







template<class LHS, class RHS>
List<LHS> assigned(const List< Assignment<LHS,RHS> >& assignments) {
    List<LHS> result; result.reserve(assignments.size());
    for(uint i=0; i!=assignments.size(); ++i) {
        result.append(assignments[i].lhs);
    }
    return result;
}

template<class T>
Set<RealVariable> real_arguments(const Expression<T>& expression) {
    Set<RealVariable> result;
    Set<UntypedVariable> args=expression.arguments();
    for(Set<UntypedVariable>::const_iterator iter=args.begin(); iter!=args.end(); ++iter) {
        if(iter->type()==variable_type<Real>()) {
            result.insert(RealVariable(iter->name()));
        }
    }
    return result;
}

template<class LHS, class RHS>
Set<RealVariable> real_arguments(const Assignment<LHS,RHS>& assignment) {
    return real_arguments(assignment.rhs);
}

template<class LHS, class RHS>
Set<RealVariable> real_arguments(const List< Assignment<LHS,RHS> >& assignments) {
    Set<RealVariable> result;
    for(uint i=0; i!=assignments.size(); ++i) {
        result.adjoin(real_arguments(assignments[i]));
    }
    return result;
}

template<class Tp, template<class> class Dec>
List< Variable<Tp> > base(const List< Dec<Tp> >& variables) {
    List< Variable<Tp> > result;
    for(uint i=0; i!=variables.size(); ++i) {
        result.append(variables[i].base());
    }
    return result;
}

template<class T>
Set<T> make_set(const List<T>& list) {
    return Set<T>(list.begin(),list.end());
}

template<class Var>
bool unique(const List<Var>& variables) {
    Set<Var> found;
    for(uint i=0; i!=variables.size(); ++i) {
        if(found.contains(variables[i])) {
            return false;
        } else {
            found.insert(variables[i]);
        }
    }
    return true;
}

template<class Var>
Set<Var> duplicates(const List<Var>& variables) {
    Set<Var> result;
    Set<Var> found;
    for(uint i=0; i!=variables.size(); ++i) {
        if(found.contains(variables[i])) {
            result.insert(variables[i]);
        } else {
            found.insert(variables[i]);
        }
    }
    return result;
}


// Order the assignments so that each variable depends only on previous ones, assuming the given input variables.
List<RealAssignment> order(const List<RealAssignment>& assignments, const Set<RealVariable>& inputs)
{
    List<RealAssignment> result;

    Map< RealVariable, Set<RealVariable> > unresolved_dependencies;
    Map< RealVariable, uint > assignment_table;

    for(uint i=0; i!=assignments.size(); ++i) {
        unresolved_dependencies.insert(assignments[i].lhs,difference(real_arguments(assignments[i].rhs),inputs));
        assignment_table.insert(assignments[i].lhs,i);
    }

    while(!unresolved_dependencies.empty()) {
        // Check for an which has been resolved and add to list of ok equations
        bool found=false;
        for(Map< RealVariable, Set<RealVariable> >::iterator iter=unresolved_dependencies.begin();
            iter!=unresolved_dependencies.end(); ++iter)
        {
            if(iter->second.empty()) {
                RealVariable variable=iter->first;
                result.append(assignments[assignment_table[variable]]);
                for(Map< RealVariable, Set<RealVariable> >::iterator rem_iter=unresolved_dependencies.begin();
                        rem_iter!=unresolved_dependencies.end(); ++rem_iter) {
                    rem_iter->second.erase(variable);
                }
                unresolved_dependencies.erase(iter);
                found=true;
                break;
            }
        }

        if(!found) {
            ARIADNE_THROW(SystemSpecificationError,"order(List<Assignment>,Set<Variable>)","Cannot order assignments "<<assignments
                            <<" with inputs "<<inputs<<" due to algebraic loops or undefined variables.");
        }
    }

    assert(result.size()==assignments.size());
    return result;
}


void
CompositeHybridAutomaton::check_mode(DiscreteLocation location) const {
    List<RealAssignment> equations=this->algebraic_assignments(location);
    List<DottedRealAssignment> dynamics=this->differential_assignments(location);

    List<RealVariable> state_variables=base(assigned(dynamics));
    List<RealVariable> auxiliary_variables=assigned(equations);
    List<RealVariable> location_variables=catenate(state_variables,auxiliary_variables);

    if(!unique(location_variables)) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<duplicates(location_variables)<<" in location "<<location<<" with defining equations "<<equations<<" and "<<dynamics<<" have more than one defining equation.");
    }

    Set<RealVariable> result_variables=make_set(location_variables);

    Set<RealVariable> undefined_variables=difference(real_arguments(equations),result_variables);
    if(!undefined_variables.empty()) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<undefined_variables<<" are used in the definition of the algebraic equations "<<equations<<" in location "<<location<<" with state variables "<<state_variables<<", but are not defined.");
    }

    undefined_variables=difference(real_arguments(dynamics),result_variables);
    if(!undefined_variables.empty()) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<undefined_variables<<" are used in the definition of the differential equations "<<dynamics<<" in location "<<location<<" with state variables "<<state_variables<<", but are not defined.");
    }

    List<RealAssignment> ordered_equations=Ariadne::order(equations,make_set(state_variables));


/*
    Map<DiscreteEvent,ContinuousPredicate> const& invariants=this->invariant_predicates(location);
    for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator invariant_iter=invariants.begin();
        invariant_iter!=invariants.end(); ++invariant_iter)
    {
        DiscreteEvent action=invariant_iter->first;
        const ContinuousPredicate& invariant=invariant_iter->second;
        if(!subset(real_arguments(invariant),result_variables)) {
            undefined_variables=difference(real_arguments(invariant),result_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the invariant "<<invariant<<" with label \""<<action<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }
    }

    Set<DiscreteEvent> events=this->events(location);
    for(Set<DiscreteEvent>::const_iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {

        // Get transition information
        DiscreteEvent event=*event_iter;
        DiscreteLocation target=this->target(location,event);
        ContinuousPredicate guard=this->guard_predicate(location,event);
        List<RealUpdateAssignment> reset=this->update_assignments(location,event);

        Set<RealVariable> target_state_variables=make_set(this->state_variables(target));
        List<RealVariable> reset_variables=base(assigned(reset));

        // Check arguments of guard
        if(!subset(real_arguments(guard),result_variables)) {
            undefined_variables=difference(real_arguments(guard),result_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the guard "<<guard<<" for event \""<<event<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }

        // Check arguments of reset
        if(!subset(real_arguments(reset),result_variables)) {
            undefined_variables=difference(real_arguments(reset),result_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are used in the reset "<<reset<<" for event \""<<event<<"\" in location "<<location<<" with variables "<<location_variables<<", but are not defined.");
        }

        // Check duplication of reset
        if(!unique(reset_variables)) {
        ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                      "Variables "<<duplicates(reset_variables)<<" have more than one defining equation in resets "<<reset<<" for event \""<<event<<"\" in location "<<location<<".");
        }

        Set<RealVariable> updated_variables=make_set(reset_variables);

        if(!subset(updated_variables,target_state_variables)) {
            Set<RealVariable> extra_variables=difference(updated_variables,target_state_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<extra_variables<<" are defined in the reset "<<reset<<" for event \""<<event<<"\" in location "<<location<<", but are not state variables "<<target_state_variables<<" in the target "<<target<<".");
        }

        if(!subset(target_state_variables,updated_variables)) {
            undefined_variables=difference(target_state_variables,updated_variables);
            ARIADNE_THROW(SystemSpecificationError,"CompositeHybridAutomaton::check_mode(DiscreteLocation)",
                          "Variables "<<undefined_variables<<" are state variables in location "<<target<<", but are not defined in the reset "<<reset<<" for the transition \""<<event<<"\" from location "<<location<<".");
        }

    } // Finished checking transitions
*/
    return;
}


void
CompositeHybridAutomaton::check_reachable_modes(DiscreteLocation initial_location) const
{
    Set<DiscreteLocation> initial_locations;
    initial_locations.insert(initial_location);
    this->check_reachable_modes(initial_locations);
}


void
CompositeHybridAutomaton::check_reachable_modes(const Set<DiscreteLocation>& initial_locations) const
{
    Set<DiscreteLocation> reachable_locations=this->discrete_reachability(initial_locations);
    for(Set<DiscreteLocation>::const_iterator iter=reachable_locations.begin(); iter!=reachable_locations.end(); ++iter) {
        this->check_mode(*iter);
    }
}

Set<DiscreteLocation>
CompositeHybridAutomaton::discrete_reachability(DiscreteLocation initial_location) const
{
    Set<DiscreteLocation> initial_locations;
    initial_locations.insert(initial_location);
    return this->discrete_reachability(initial_locations);
}

Set<DiscreteLocation>
CompositeHybridAutomaton::discrete_reachability(const Set<DiscreteLocation>& initial_locations) const
{
    const CompositeHybridAutomaton& automaton=*this;

    typedef uint Nat;

    Set<DiscreteLocation> reached=initial_locations;
    Set<DiscreteLocation> working=initial_locations;
    Set<DiscreteLocation> found;

    Map<DiscreteLocation,Nat> steps;
    Nat step=0u;

    for(Set<DiscreteLocation>::const_iterator initial_iter=initial_locations.begin(); initial_iter!=initial_locations.end(); ++initial_iter) {
        steps.insert(*initial_iter,step);
    }

    while(!working.empty()) {
        ++step;
        for(Set<DiscreteLocation>::const_iterator source_iter=working.begin(); source_iter!=working.end(); ++source_iter) {
            ARIADNE_LOG(5,"new_mode\n");
            DiscreteLocation location=*source_iter;
            ARIADNE_LOG(5,"  mode: "<<location<<":\n");
            ARIADNE_LOG(7,"      output="<<automaton.algebraic_assignments(location)<<"\n");
            ARIADNE_LOG(7,"          function="<<automaton.output_function(location)<<"\n");
            ARIADNE_LOG(7,"      dynamic="<<automaton.differential_assignments(location)<<"\n");
            ARIADNE_LOG(7,"          function="<<automaton.dynamic_function(location)<<"\n");

/*
            Map<DiscreteEvent,ContinuousPredicate> invariants=automaton.invariant_predicates(location);
            for(Map<DiscreteEvent,ContinuousPredicate>::const_iterator invariant_iter=invariants.begin();
                    invariant_iter!=invariants.end(); ++invariant_iter) {
                ARIADNE_LOG(7,"    invariant: "<<invariant_iter->second<<"\n");
                ARIADNE_LOG(7,"      function: "<<automaton.invariant_function(location,invariant_iter->first)<<"\n");
                ARIADNE_LOG(7,"        action: "<<invariant_iter->first<<"\n");
            }
*/
            Set<DiscreteEvent> events=automaton.transition_events(location);
            ARIADNE_LOG(5,"  events: "<<events<<"\n");
            for(Set<DiscreteEvent>::const_iterator event_iter=events.begin(); event_iter!=events.end(); ++event_iter) {
                DiscreteEvent event=*event_iter;
                ARIADNE_LOG(7,"    event:"<<event<<"\n");
                DiscreteLocation target=automaton.target(location,event);
                ARIADNE_LOG(7,"    transition: "<<event<<" -> "<<target<<"\n");
                ARIADNE_LOG(7,"        reset="<<automaton.update_assignments(location,event)<<"\n");
                ARIADNE_LOG(7,"          function="<<automaton.reset_function(location,event)<<"\n");
                ARIADNE_LOG(7,"        guard="<<automaton.guard_predicate(location,event)<<"\n");
                ARIADNE_LOG(7,"          function="<<automaton.guard_function(location,event)<<"\n");
                if(!reached.contains(target)) {
                    found.insert(target);
                    reached.insert(target);
                    steps.insert(target,step);
                }
           }

        }
        ARIADNE_LOG(5,"\nstep "<<step<<" found: "<<found<<"\n\n");
        working.clear();
        working.swap(found);
    }

    return reached;
}

CompositeHybridAutomaton parallel_composition(const List<AtomicHybridAutomaton>& components)
{
    return CompositeHybridAutomaton(components);
}



}
