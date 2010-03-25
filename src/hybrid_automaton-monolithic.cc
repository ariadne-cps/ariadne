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


uint
DiscreteMode::
dimension() const
{
    return this->_dynamic.argument_size();
}


DiscreteMode::
DiscreteMode(DiscreteLocation location,
             const VectorFunction& dynamic)
    :  _location(location), _dynamic(dynamic), _invariants(), _grid(new Grid(dynamic.argument_size()))
{
}

std::ostream&
operator<<(std::ostream& os, const DiscreteMode& mode)
{
    return os << "DiscreteMode( "
              << "location=" << mode.location() << ", "
              << "dynamic=" << mode.dynamic() << ", "
              << "invariants=" << mode.invariants() << " )";
}





DiscreteTransition::
DiscreteTransition(DiscreteEvent event,
                   const DiscreteMode& source,
                   const DiscreteMode& target,
                   const VectorFunction& reset,
                   const ScalarFunction& activation,
                   bool forced)
    : _event(event), _source(&source), _target(&target),
      _activation(activation), _reset(reset), _forced(forced)
{
    ARIADNE_ASSERT(activation.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.result_size()==target.dimension());
}

std::ostream&
operator<<(std::ostream& os, const DiscreteTransition& transition)
{
    return os << "DiscreteTransition( "
              << "event=" << transition.event() << ", "
              << "source=" << transition.source() << ", "
              << "target=" << transition.target() << ", "
              << "reset=" << transition.reset() << ", "
              << "activation=" << transition.activation() << " )";
}






MonolithicHybridAutomaton::~MonolithicHybridAutomaton()
{
}

MonolithicHybridAutomaton::MonolithicHybridAutomaton()
{
}

MonolithicHybridAutomaton::MonolithicHybridAutomaton(const std::string& name)
    : _name(name)
{
}





const DiscreteMode&
MonolithicHybridAutomaton::new_mode(DiscreteLocation location,
                          const VectorFunction& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }
    if(dynamic.result_size()!=dynamic.argument_size()) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_mode(location,dynamic)",
            "The dynamic has argument size " << dynamic.argument_size()
                << " and result size " << dynamic.result_size() << ", so does not define a vector field.");
    }
    this->_modes.insert(DiscreteMode(location,dynamic));
    return this->mode(location);
}


const DiscreteMode&
MonolithicHybridAutomaton::new_invariant(DiscreteLocation location,
                               const ScalarFunction& invariant)
{
    if(!this->has_mode(location)) {
        throw std::runtime_error("The location of the invariant must be in the automaton.");
    }
    DiscreteMode& mode=const_cast<DiscreteMode&>(this->mode(location));
    if(invariant.argument_size()!=mode.dimension()) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_invariant(location,invariant)",
            "The invariant has argument size " << invariant.argument_size()
                << " but the mode has state-space dimension " << mode.dimension());
    }
    DiscreteEvent invariant_event(-8-mode._invariants.size());
    mode._invariants[invariant_event]=invariant;
    return mode;
}


const DiscreteMode&
MonolithicHybridAutomaton::new_invariant(DiscreteLocation location,
                               const VectorFunction& invariant)
{
    if(invariant.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_invariant(location,invariant)",
            "The invariant has result size " << invariant.result_size()
                << " but only scalar constraints are currently supported.");
    }
    return this->new_invariant(location, invariant[0]);
}


const DiscreteTransition&
MonolithicHybridAutomaton::new_transition(DiscreteEvent event,
                                DiscreteLocation source,
                                DiscreteLocation target,
                                const VectorFunction &reset,
                                const ScalarFunction &activation,
                                bool forced)
{
    ARIADNE_ASSERT_MSG(event>0, "Transition event should be positive.");
    if(this->has_transition(event,source)) {
        throw std::runtime_error("The automaton already has a transition with the given id and source id.");
    }
    if(!this->has_mode(source)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given source id");
    }
    if(!this->has_mode(target)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
    }

    const DiscreteMode& source_mode=this->mode(source);
    const DiscreteMode& target_mode=this->mode(target);
    this->_transitions.insert(DiscreteTransition(event,source_mode,target_mode,reset,activation,forced));
    return this->transition(event,source);
}


const DiscreteTransition&
MonolithicHybridAutomaton::new_transition(DiscreteEvent event,
                                DiscreteLocation source,
                                DiscreteLocation target,
                                const VectorFunction &reset,
                                const VectorFunction &activation,
                                bool forced)
{
    if(activation.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"MonolithicHybridAutomaton::new_transition(...)",
            "The activation " << activation << " has result size " << activation.result_size()
                << " but only scalar constraints are currently supported.");
    }
    return this->new_transition(event,source,target,reset,activation[0],forced);
}

const DiscreteTransition&
MonolithicHybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      DiscreteLocation source,
                      DiscreteLocation target,
                      const VectorFunction &reset,
                      const VectorFunction &activation)
{
    return new_transition(event,source,target,reset,activation,true);
}


const DiscreteTransition&
MonolithicHybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        DiscreteLocation source,
                        DiscreteLocation target,
                        const VectorFunction &reset,
                        const VectorFunction &activation)
{
    return new_transition(event,source,target,reset,activation,false);
}





void
MonolithicHybridAutomaton::set_grid(DiscreteLocation location,
                                    const Grid& grid)
{
    if(!this->has_mode(location)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given location id");
    }
    DiscreteMode& mode=const_cast<DiscreteMode&>(this->mode(location));
    if(grid.dimension()!=mode.dimension()) {
        throw std::runtime_error("The mode of the automaton has a different dimension to the grid.");
    }
    mode._grid=shared_ptr<Grid>(new Grid(grid));
}

void
MonolithicHybridAutomaton::set_grid(const Grid& grid)
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        DiscreteMode& mode=const_cast<DiscreteMode&>(*mode_iter);
        if(grid.dimension()!=mode.dimension()) {
            throw std::runtime_error("The automaton has a different dimension to the grid.");
        }
        mode._grid=shared_ptr<Grid>(new Grid(grid));
    }
}

void
MonolithicHybridAutomaton::set_grid(const HybridGrid& hgrid)
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        DiscreteMode& mode=const_cast<DiscreteMode&>(*mode_iter);
        DiscreteLocation loc = mode.location();
        if(!hgrid.has_location(loc)) {
            throw std::runtime_error("The automaton does not contain a mode with this given location id");
        }
        if(hgrid[loc].dimension()!=mode.dimension()) {
            throw std::runtime_error("The mode of the automaton has a different dimension to the grid.");
        }
        mode._grid=shared_ptr<Grid>(new Grid(hgrid[loc]));
    }
}







bool
MonolithicHybridAutomaton::has_mode(DiscreteLocation location) const
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==location) {
                return true;
            }
        }
    return false;
}


bool
MonolithicHybridAutomaton::has_transition(DiscreteLocation source, DiscreteEvent event) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return true;
            }
        }
    return false;
}

bool
MonolithicHybridAutomaton::has_transition(DiscreteEvent event, DiscreteLocation source) const
{
    return this->has_transition(source,event);
}

const DiscreteMode&
MonolithicHybridAutomaton::mode(DiscreteLocation location) const
{
    for(Set<DiscreteMode>::const_iterator iter=this->_modes.begin(); iter!=this->_modes.end(); ++iter) {
        if(iter->location()==location) {
            return *iter;
        }
    }
    ARIADNE_FAIL_MSG("MonolithicHybridAutomaton "<<*this<<" has no mode "<<location<<"\n");
}


const DiscreteTransition&
MonolithicHybridAutomaton::transition(DiscreteLocation source, DiscreteEvent event) const
{
    return this->transition(event,source);
}

Set<DiscreteTransition>
MonolithicHybridAutomaton::transitions(DiscreteLocation source) const
{
    ARIADNE_DEPRECATED("MonolithicHybridAutomaton::transitions(DiscreteLocation)","Use guards(DiscreteLocation) to get events instead.");
    Set<DiscreteTransition> result;
    for(Set<DiscreteTransition>::const_iterator iter=this->_transitions.begin();
        iter!=this->_transitions.end(); ++iter)
    {
        if(iter->source()==source) {
            result.insert(*iter);
        }
    }
    return result;
}


uint
MonolithicHybridAutomaton::dimension(DiscreteLocation location) const
{
    return this->mode(location).dynamic().result_size();
}

VectorFunction
MonolithicHybridAutomaton::dynamic(DiscreteLocation location) const
{
    return this->mode(location).dynamic();
}

Map<DiscreteEvent,ScalarFunction>
MonolithicHybridAutomaton::invariants(DiscreteLocation location) const
{
    return this->mode(location).invariants();
}

Map<DiscreteEvent,ScalarFunction>
MonolithicHybridAutomaton::guards(DiscreteLocation location) const
{
    Map<DiscreteEvent,ScalarFunction> result;
    for(Set<DiscreteTransition>::const_iterator iter=this->_transitions.begin();
        iter!=this->_transitions.end(); ++iter)
    {
        if(iter->source()==location && iter->urgency()==urgent) {
            result.insert(iter->event(),iter->guard());
        }
    }
    return result;
}

Map<DiscreteEvent,ScalarFunction>
MonolithicHybridAutomaton::activations(DiscreteLocation location) const
{
    Map<DiscreteEvent,ScalarFunction> result;
    for(Set<DiscreteTransition>::const_iterator iter=this->_transitions.begin();
        iter!=this->_transitions.end(); ++iter)
    {
        if(iter->source()==location && iter->urgency()==permissive) {
            result.insert(iter->event(),iter->activation());
        }
    }
    return result;
}

DiscreteLocation
MonolithicHybridAutomaton::target(DiscreteLocation location, DiscreteEvent event) const
{
    return this->transition(location,event).target();
}

VectorFunction
MonolithicHybridAutomaton::reset(DiscreteLocation location, DiscreteEvent event) const
{
    return this->transition(location,event).reset();
}


Map<DiscreteEvent,DiscreteLocation>
MonolithicHybridAutomaton::targets(DiscreteLocation location) const
{
    ARIADNE_DEPRECATED("MonolithicHybridAutomaton::targets(DiscreteLocation)","Use guards(DiscreteLocation) to get events and target(DiscreteLocation,DiscreteEvent) instead.");
    Map<DiscreteEvent,DiscreteLocation> result;
    Set<DiscreteTransition> transitions=this->transitions(location);
    for(Set<DiscreteTransition>::const_iterator iter=transitions.begin(); iter!=transitions.end(); ++iter) {
        result.insert(iter->event(),iter->target());
    }
    return result;
}

Map<DiscreteEvent,VectorFunction>
MonolithicHybridAutomaton::resets(DiscreteLocation location) const
{
    ARIADNE_DEPRECATED("MonolithicHybridAutomaton::resets(DiscreteLocation)","Use guards(DiscreteLocation) to get events and reset(DiscreteLocation,DiscreteEvent) instead.");
    Map<DiscreteEvent,VectorFunction> result;
    Set<DiscreteTransition> transitions=this->transitions(location);
    for(Set<DiscreteTransition>::const_iterator iter=transitions.begin(); iter!=transitions.end(); ++iter) {
        result.insert(iter->event(),iter->reset());
    }
    return result;
}


HybridSet
MonolithicHybridAutomaton::invariant() const
{
    ARIADNE_NOT_IMPLEMENTED;
}



HybridSpace
MonolithicHybridAutomaton::state_space() const
{
    HybridSpace result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            result[mode_iter->location()]=mode_iter->dimension();
        }
    return result;
}


const Set< DiscreteMode >&
MonolithicHybridAutomaton::modes() const
{
    return this->_modes;
}



const Set< DiscreteTransition >&
MonolithicHybridAutomaton::transitions() const
{
    return this->_transitions;
}




Map<DiscreteEvent,ScalarFunction>
MonolithicHybridAutomaton::blocking_guards(DiscreteLocation source) const
{
    Map<DiscreteEvent,ScalarFunction> result;
    const DiscreteMode& mode=this->mode(source);
    for(Map<DiscreteEvent,ScalarFunction>::const_iterator invariant_iter=mode._invariants.begin();
        invariant_iter!=mode._invariants.end(); ++invariant_iter)
    {
        const DiscreteEvent event=invariant_iter->first;
        const ScalarFunction invariant=invariant_iter->second;
        result[event]=invariant;
    }

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source()==source && transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const ScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}


Map<DiscreteEvent,ScalarFunction>
MonolithicHybridAutomaton::permissive_guards(DiscreteLocation source) const
{
    Map<DiscreteEvent,ScalarFunction> result;

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source()==source && !transition_iter->forced()) {
            const DiscreteEvent event=transition_iter->event();
            const ScalarFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}



const DiscreteTransition&
MonolithicHybridAutomaton::transition(DiscreteEvent event, DiscreteLocation source) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return *transition_iter;
            }
        }
    throw std::runtime_error("The hybrid automaton does not have a transition with the given event and source.");
}


const std::string&
MonolithicHybridAutomaton::name() const
{
    return this->_name;
}

Grid
MonolithicHybridAutomaton::grid(DiscreteLocation location) const
{
    ARIADNE_ASSERT(this->has_mode(location));
    const DiscreteMode& mode=this->mode(location);
    return Grid(mode.dimension());
}

HybridGrid
MonolithicHybridAutomaton::grid() const
{
    HybridGrid result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
    {
        result[mode_iter->location()]=mode_iter->grid();
    }
    return result;
}

Set<DiscreteEvent>
MonolithicHybridAutomaton::urgent_events(DiscreteLocation source) const
{
    Set<DiscreteEvent> result;
    for(Set<DiscreteTransition>::const_iterator iter=this->_transitions.begin();
        iter!=this->_transitions.end(); ++iter)
    {
        if(iter->source()==source && iter->urgency()==urgent) {
            result.insert(iter->event());
        }
    }
    return result;
}

Set<DiscreteEvent>
MonolithicHybridAutomaton::permissive_events(DiscreteLocation source) const
{
    Set<DiscreteEvent> result;
    for(Set<DiscreteTransition>::const_iterator iter=this->_transitions.begin();
        iter!=this->_transitions.end(); ++iter)
    {
        if(iter->source()==source && iter->urgency()==permissive) {
            result.insert(iter->event());
        }
    }
    return result;
}

Set<DiscreteEvent>
MonolithicHybridAutomaton::invariant_events(DiscreteLocation location) const
{
    return this->mode(location)._invariants.keys();
}


Set<DiscreteEvent> MonolithicHybridAutomaton::blocking_events(DiscreteLocation location) const
{
    return join(this->urgent_events(location),this->invariant_events(location));
}

Set<DiscreteEvent> MonolithicHybridAutomaton::transition_events(DiscreteLocation source) const
{
    Set<DiscreteEvent> result;
    for(Set<DiscreteTransition>::const_iterator iter=this->_transitions.begin();
        iter!=this->_transitions.end(); ++iter)
    {
        if(iter->source()==source) {
            result.insert(iter->event());
        }
    }
    return result;
}

VectorFunction
MonolithicHybridAutomaton::dynamic_function(DiscreteLocation location) const
{
    return this->dynamic(location);
}

VectorFunction
MonolithicHybridAutomaton::reset_function(DiscreteLocation source, DiscreteEvent event) const
{
    return this->reset(source,event);
}

ScalarFunction
MonolithicHybridAutomaton::guard_function(DiscreteLocation source, DiscreteEvent event) const
{
    return this->guards(source)[event];
}

ScalarFunction
MonolithicHybridAutomaton::invariant_function(DiscreteLocation source, DiscreteEvent event) const
{
    return this->invariants(source)[event];
}


std::ostream&
operator<<(std::ostream& os, const MonolithicHybridAutomaton& ha)
{
    return os << "MonolithicHybridAutomaton( modes=" << ha.modes() << ", transitions=" << ha.transitions() << ")";
}




}
