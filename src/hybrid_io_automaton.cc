/***************************************************************************
 *            hybrid_io_automaton.cc
 *
 *  Copyright  2010  Davide Bresolin
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
#include "hybrid_io_automaton.h"
#include "hybrid_automaton.h"
#include "assignment.h"
#include "operators.h"
#include "space.h"
#include "grid.h"


namespace Ariadne {

//
//  Auxiliary functions
//

// Checks whether a set is disjoint from a map (i.e., no common keys)
template< class X >
bool disjoint(const std::set<RealVariable>& s, const std::map< RealVariable, X >& m) {
    // Scan all elements of s
    for(std::set<RealVariable>::iterator iter=s.begin(); iter != s.end(); iter++) {
        if(m.find(*iter) != m.end()) {
            return false;   // Found a common element, set is not disjoint from the map
        }
    }
    return true;    // No common elements, set and map are disjoint
}

// Checks whether a map is a subset of a set (i.e., all keys belongs to the set)
template< class X >
bool subset(const std::map< RealVariable, X >& m, const std::set<RealVariable>& s) {
    // Scan all elements of m
    for(typename std::map< RealVariable, X >::const_iterator iter=m.begin(); iter != m.end(); iter++) {
        if(s.find(iter->first) == s.end()) {
            return false;   // Found an element of m that is not in s, map is not a subset of the set
        }
    }
    return true;    // All elements of m are in s, map is a subset of the set
}


//
//  DiscreteIOMode
//
const RealExpression&
DiscreteIOMode::dynamics(const RealVariable& var) 
{
    std::map< RealVariable, RealExpression >::const_iterator iter;
    iter = this->_dynamics.find(var);
    if(iter == this->_dynamics.end()) {     // Dynamics not defined for var
        ARIADNE_FAIL_MSG("Dynamics for variable " << var << " not defined in location " << this->_location << ".");
    }
    return iter->second;
}

DiscreteIOMode::
DiscreteIOMode(DiscreteState location)
    :  _location(location), _dynamics(), _invariants()
{
}

DiscreteIOMode::
DiscreteIOMode(DiscreteState location,
               const std::map< RealVariable, RealExpression >& dynamics)
    :  _location(location), _dynamics(dynamics), _invariants()
{
}

DiscreteIOMode::
DiscreteIOMode(DiscreteState location,
               const std::map< RealVariable, RealExpression >& dynamics,
               const std::list< RealExpression >& invariants)
    :  _location(location), _dynamics(dynamics), _invariants(invariants)
{
}

void
DiscreteIOMode::set_dynamics(const RealVariable& var,
                            const RealExpression& dyn)
{
    this->_dynamics[var] = dyn;
}

void
DiscreteIOMode::add_invariant(const RealExpression& inv)
{
    this->_invariants.push_back(inv);
}

void 
DiscreteIOMode::substitute(const Constant<Real>& con, const Real& c) 
{
    for(std::map<RealVariable, RealExpression>::iterator it=this->_dynamics.begin();it!=this->_dynamics.end();it++)
        it->second.substitute(con,c);
    for(std::list<RealExpression>::iterator it=this->_invariants.begin();it!=this->_invariants.end();it++)
        it->substitute(con,c);
}


std::ostream&
operator<<(std::ostream& os, const DiscreteIOMode& mode)
{
    return os << "DiscreteIOMode( "
              << "location=" << mode.location() << ", "
              << "dynamics=" << mode.dynamics() << ", "
              << "invariants=" << mode.invariants() << " )";
}


//
// DiscreteIOTransition
//
const RealExpression&
DiscreteIOTransition::reset(const RealVariable& var) 
{
    std::map< RealVariable, RealExpression >::const_iterator iter;
    iter = this->_reset.find(var);
    if(iter == this->_reset.end()) {     // Reset not defined for var
        ARIADNE_FAIL_MSG("Reset for variable " << var << " not defined in transition " << this->_event 
                << " from mode " << this->_source << ".");
    }
    return iter->second;
}

DiscreteIOTransition::
DiscreteIOTransition()
    : _event(), _source(), _target(),
      _activation(),
      _reset(), _forced(false)
{
}


DiscreteIOTransition::
DiscreteIOTransition(DiscreteEvent event,
                     DiscreteState source,
                     DiscreteState target,
                     bool forced)
    : _event(event), _source(source), _target(target),
      _activation(RealExpression(1.0)),     // We assume the transition to be always active
      _reset(), _forced(forced)
{
}

DiscreteIOTransition::
DiscreteIOTransition(DiscreteEvent event,
                     DiscreteState source,
                     DiscreteState target,
                     const RealExpression& activation,
                     bool forced)
    : _event(event), _source(source), _target(target),
      _activation(activation), _reset(), _forced(forced)
{
}

DiscreteIOTransition::
DiscreteIOTransition(DiscreteEvent event,
                     DiscreteState source,
                     DiscreteState target,
                     const std::map< RealVariable, RealExpression >& reset,
                     bool forced)
    : _event(event), _source(source), _target(target),
      _activation(RealExpression(1.0)),     // We assume the transition to be always active
      _reset(reset), _forced(forced)
{
}

DiscreteIOTransition::
DiscreteIOTransition(DiscreteEvent event,
                     DiscreteState source,
                     DiscreteState target,
                     const std::map< RealVariable, RealExpression >& reset,
                     const RealExpression& activation,
                     bool forced)
    : _event(event), _source(source), _target(target),
      _activation(activation), _reset(reset), _forced(forced)
{
}

void
DiscreteIOTransition::set_event(DiscreteEvent event)
{
    this->_event = event;
}

void
DiscreteIOTransition::set_source(DiscreteState source)
{
    this->_source = source;
}

void
DiscreteIOTransition::set_target(DiscreteState target)
{
    this->_target = target;
}

void
DiscreteIOTransition::set_reset(const RealVariable& var,
                                const RealExpression& res)
{
    this->_reset[var] = res;
}

void
DiscreteIOTransition::set_reset(const std::map< RealVariable, RealExpression >& reset)
{
    this->_reset = reset;
}


void
DiscreteIOTransition::set_activation(const RealExpression& inv)
{
    this->_activation = inv;
}

void 
DiscreteIOTransition::substitute(const Constant<Real>& con, const Real& c) 
{
    this->_activation.substitute(con,c);
    for(std::map<RealVariable, RealExpression>::iterator it=this->_reset.begin(); it != this->_reset.end(); it++)
        it->second.substitute(con,c);
}


std::ostream&
operator<<(std::ostream& os, const DiscreteIOTransition& transition)
{
    return os << "DiscreteIOTransition( "
              << "event=" << transition.event() << ", "
              << "source=" << transition.source() << ", "
              << "target=" << transition.target() << ", "
              << "reset=" << transition.reset() << ", "
              << "activation=" << transition.activation() << ", "
              << (transition.forced() ? "forced" : "unforced") << " )";
}


//
// HybridIOAutomaton
//
HybridIOAutomaton::
HybridIOAutomaton()
    : _name(""), 
      _input_vars(), _output_vars(), _internal_vars(),
      _input_events(), _output_events(), _internal_events(),
      _grid()
{
}

HybridIOAutomaton::
HybridIOAutomaton(const std::string& name)
    : _name(name), 
      _input_vars(), _output_vars(), _internal_vars(),
      _input_events(), _output_events(), _internal_events(),
      _grid()
{
}

HybridIOAutomaton::
HybridIOAutomaton(const std::string& name,
                  const std::set< RealVariable >& input_vars,
                  const std::set< RealVariable >& output_vars,
                  const std::set< RealVariable >& internal_vars)
    : _name(name), 
      _input_vars(input_vars), _output_vars(output_vars), _internal_vars(internal_vars),
      _input_events(), _output_events(), _internal_events(),
      _grid()
{
    // Input, output, and internal events and variables must be disjoint
    ARIADNE_ASSERT_MSG(disjoint<RealVariable>(input_vars, output_vars), 
        "Input and output variables are not disjoint in the definition of automaton " << name);
    ARIADNE_ASSERT_MSG(disjoint<RealVariable>(input_vars, internal_vars), 
        "Input and internal variables are not disjoint in the definition of automaton " << name);
    ARIADNE_ASSERT_MSG(disjoint<RealVariable>(output_vars, internal_vars), 
        "Output and internal variables are not disjoint in the definition of automaton " << name);
}


HybridIOAutomaton::
HybridIOAutomaton(const std::string& name,
                  const std::set< RealVariable >& input_vars,
                  const std::set< RealVariable >& output_vars,
                  const std::set< RealVariable >& internal_vars,
                  const std::set< DiscreteEvent >& input_events,
                  const std::set< DiscreteEvent >& output_events,
                  const std::set< DiscreteEvent >& internal_events)
    : _name(name), 
      _input_vars(input_vars), _output_vars(output_vars), _internal_vars(internal_vars),
      _input_events(input_events), _output_events(output_events), _internal_events(internal_events),
      _grid()
{
    // Input, output, and internal events and variables must be disjoint
    ARIADNE_ASSERT_MSG(disjoint<RealVariable>(input_vars, output_vars), 
        "Input and output variables are not disjoint in the definition of automaton " << name);
    ARIADNE_ASSERT_MSG(disjoint<RealVariable>(input_vars, internal_vars), 
        "Input and internal variables are not disjoint in the definition of automaton " << name);
    ARIADNE_ASSERT_MSG(disjoint<RealVariable>(output_vars, internal_vars), 
        "Output and internal variables are not disjoint in the definition of automaton " << name);

    ARIADNE_ASSERT_MSG(disjoint<DiscreteEvent>(input_events, output_events), 
        "Input and output variables are not disjoint in the definition of automaton " << name);
    ARIADNE_ASSERT_MSG(disjoint<DiscreteEvent>(input_events, internal_events), 
        "Input and internal variables are not disjoint in the definition of automaton " << name);
    ARIADNE_ASSERT_MSG(disjoint<DiscreteEvent>(output_events, internal_events), 
        "Output and internal variables are not disjoint in the definition of automaton " << name);
}

const std::string&
HybridIOAutomaton::set_name(const std::string& name)
{
    ARIADNE_ASSERT_MSG(name != "",
        "Automaton's name must be nonempty.")
    
    this->_name = name;
    return this->_name;
}

const std::set< RealVariable >& 
HybridIOAutomaton::add_input_var(const RealVariable& u)
{
    ARIADNE_ASSERT_MSG(!this->has_output_var(u),
        "Variable " << u << " is an output variable for automaton " << this->_name);
    ARIADNE_ASSERT_MSG(!this->has_internal_var(u),
        "Variable " << u << " is an internal variable for automaton " << this->_name);
    this->_input_vars.insert(u);
    return this->_input_vars;
}

const std::set< RealVariable >& 
HybridIOAutomaton::add_output_var(const RealVariable& y)
{
    ARIADNE_ASSERT_MSG(!this->has_input_var(y),
        "Variable " << y << " is an input variable for automaton " << this->_name);
    ARIADNE_ASSERT_MSG(!this->has_internal_var(y),
        "Variable " << y << " is an internal variable for automaton " << this->_name);
    this->_output_vars.insert(y);
    return this->_output_vars;
}

const std::set< RealVariable >& 
HybridIOAutomaton::add_internal_var(const RealVariable& x)
{
    ARIADNE_ASSERT_MSG(!this->has_output_var(x),
        "Variable " << x << " is an output variable for automaton " << this->_name);
    ARIADNE_ASSERT_MSG(!this->has_input_var(x),
        "Variable " << x << " is an input variable for automaton " << this->_name);
    this->_internal_vars.insert(x);
    return this->_internal_vars;
}

const std::set< DiscreteEvent >&  
HybridIOAutomaton::add_input_event(const DiscreteEvent& e)
{
    ARIADNE_ASSERT_MSG(!this->has_output_event(e),
        "Event " << e << " is an output event for automaton " << this->_name);
    ARIADNE_ASSERT_MSG(!this->has_internal_event(e),
        "Event " << e << " is an internal event for automaton " << this->_name);
    this->_input_events.insert(e);
    return this->_input_events;
}

const std::set< DiscreteEvent >& 
HybridIOAutomaton::add_output_event(const DiscreteEvent& e)
{
    ARIADNE_ASSERT_MSG(!this->has_input_event(e),
        "Event " << e << " is an input event for automaton " << this->_name);
    ARIADNE_ASSERT_MSG(!this->has_internal_event(e),
        "Event " << e << " is an internal event for automaton " << this->_name);
    this->_output_events.insert(e);
    return this->_output_events;
}

const std::set< DiscreteEvent >&  
HybridIOAutomaton::add_internal_event(const DiscreteEvent& e)
{
    ARIADNE_ASSERT_MSG(!this->has_output_event(e),
        "Event " << e << " is an output event for automaton " << this->_name);
    ARIADNE_ASSERT_MSG(!this->has_input_event(e),
        "Event " << e << " is an input event for automaton " << this->_name);
    this->_internal_events.insert(e);
    return this->_internal_events;
}

const DiscreteIOMode&
HybridIOAutomaton::new_mode(DiscreteState location)
{
    if(this->has_mode(location)) {
        ARIADNE_FAIL_MSG("The hybrid automaton " << this->_name << " already has a mode with id " << location << ".");
    }
    this->_modes.push_back(DiscreteIOMode(location));
    return this->mode(location);
}

const DiscreteIOMode&
HybridIOAutomaton::new_mode(const DiscreteIOMode& mode)
{
    if(this->has_mode(mode.location())) {
        ARIADNE_FAIL_MSG("The hybrid automaton " << this->_name << " already has a mode with id " << mode.location() << ".");
    }
    this->_modes.push_back(mode);
    return mode;
}

const DiscreteIOMode&
HybridIOAutomaton::new_mode(DiscreteState location,
                            const std::map< RealVariable, RealExpression >& dynamics)
{
    if(this->has_mode(location)) {
        ARIADNE_FAIL_MSG("The hybrid automaton " << this->_name << " already has a mode with id " << location << ".");
    }
    if(!disjoint(this->_input_vars, dynamics)) {
        ARIADNE_FAIL_MSG("Dynamics of mode " << location << " of automaton " << this->_name << " cannot control an input variable.");
    }
    
    this->_modes.push_back(DiscreteIOMode(location, dynamics));
    return this->mode(location);
}

const DiscreteIOMode&
HybridIOAutomaton::set_dynamics(DiscreteState location,
                               const RealVariable& var,
                               const RealExpression& dyn)
{
    if(!this->has_mode(location)) {
        ARIADNE_FAIL_MSG("The hybrid automaton " << this->_name << " has no mode with id " << location << ".");
    }
    if(contains(this->_input_vars,var)) {
        ARIADNE_FAIL_MSG("Variable " << var << " is an input variable in automaton " << this->_name << ": dynamics cannot be specified.");
    }    
    if(!contains(this->_output_vars,var) && !contains(this->_internal_vars,var)) {
        ARIADNE_FAIL_MSG("Variable " << var << " is neither an output nor an internal variable in automaton " << this->_name << ": dynamics cannot be specified.");
    }        
    DiscreteIOMode& mode=this->_mode(location);
    mode.set_dynamics(var, dyn);
    return mode;
}

const DiscreteIOMode&
HybridIOAutomaton::new_invariant(DiscreteState location,
                                 const RealExpression& inv)
{
    if(!this->has_mode(location)) {
        ARIADNE_FAIL_MSG("The hybrid automaton " << this->_name << " has no mode with id " << location << ".");
    }
    DiscreteIOMode& mode=this->_mode(location);
    mode.add_invariant(inv);
    return mode;
}

const DiscreteIOTransition&
HybridIOAutomaton::new_transition(DiscreteEvent event,
                                  DiscreteState source,
                                  DiscreteState target,
                                  const std::map< RealVariable, RealExpression >& reset,
                                  const RealExpression& activation,
                                  bool forced)
{
    ARIADNE_ASSERT_MSG(!this->has_input_event(event) || identical(activation,1.0),
        "Error in transition " << event << " from " << source << " of automaton " << this->name() <<
        ": transitions labelled with input events cannot have an activation different from true.");
    this->new_transition(event, source, target, forced);
    this->set_reset(event, source, reset);
    this->set_activation(event, source, activation);
    return this->transition(event,source);
}

const DiscreteIOTransition&
HybridIOAutomaton::new_transition(DiscreteEvent event,
                                  DiscreteState source,
                                  DiscreteState target,
                                  const std::map< RealVariable, RealExpression >& reset,
                                  bool forced)
{
    this->new_transition(event, source, target, forced);
    this->set_reset(event, source, reset);
    return this->transition(event,source);
}


const DiscreteIOTransition&
HybridIOAutomaton::new_transition(DiscreteEvent event,
                                  DiscreteState source,
                                  DiscreteState target,
                                  const RealExpression& activation,
                                  bool forced)
{
    ARIADNE_ASSERT_MSG(!this->has_input_event(event) || identical(activation,1.0),
        "Error in transition " << event << " from " << source << " of automaton " << this->name() <<
        ": transitions labelled with input events cannot have an activation different from true.");
    this->new_transition(event, source, target, forced);
    this->set_activation(event, source, activation);
    return this->transition(event,source);
}


const DiscreteIOTransition&
HybridIOAutomaton::new_transition(DiscreteEvent event,
                                  DiscreteState source,
                                  DiscreteState target,
                                  bool forced)
{
    if(!event.is_transition()) {
        ARIADNE_FAIL_MSG("Transition event names cannot start with \"invariant\".");    
    }
    if(forced && this->has_input_event(event)) {
        ARIADNE_FAIL_MSG("Event " << event << " is an input event in automaton " << this->_name << 
            ": transition cannot be forced.");
    }    
    if(this->has_transition(event,source)) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " already has a transition with id "
            << event << " and source " << source << ".");
    }
    if(!this->has_mode(source)) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " does not contain a source mode with id " << source);
    }
    if(!this->has_mode(target)) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " does not contain a target mode with id " << target);
    }
          
    this->_transitions.push_back(DiscreteIOTransition(event,source,target,forced));
    return this->transition(event,source);
}

const DiscreteIOTransition&
HybridIOAutomaton::new_transition(const DiscreteIOTransition& trans)
{
    if(!trans.event().is_transition()) {
        ARIADNE_FAIL_MSG("Transition event names cannot start with \"invariant\".");    
    }
    if(this->has_transition(trans.event(),trans.source())) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " already has a transition with id "
            << trans.event() << " and source " << trans.source() << ".");
    }
    if(!this->has_mode(trans.source())) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " does not contain a source mode with id " << trans.source());
    }
    if(!this->has_mode(trans.target())) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " does not contain a target mode with id " << trans.target());
    }
          
    this->_transitions.push_back(trans);
    return trans;
}


const DiscreteIOTransition&
HybridIOAutomaton::set_reset(DiscreteEvent event,
                             DiscreteState source,
                             const std::map< RealVariable, RealExpression >& reset)
{
    if(!this->has_transition(event,source)) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " has no transition with event "
            << event << " and source " << source << ".");
    }
    if(!disjoint(this->_input_vars, reset)) {
        ARIADNE_FAIL_MSG("Reset function in transition " << event << " from mode " << source << 
            " cannot affect an input variable.");
    }
    
    DiscreteIOTransition& trans = this->_transition(event,source);
    trans.set_reset(reset);
    return trans;
}

const DiscreteIOTransition&
HybridIOAutomaton::set_reset(DiscreteEvent event,
                             DiscreteState source,
                             const RealVariable& var,
                             const RealExpression& reset)
{
    if(!this->has_transition(event,source)) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " has no transition with event "
            << event << " and source " << source << ".");
    }
    if(contains(this->_input_vars,var)) {
        ARIADNE_FAIL_MSG("Reset function in transition " << event << " from mode " << source << 
            " cannot affect input variable " << var << ".");
    }
    
    DiscreteIOTransition& trans = this->_transition(event,source);
    trans.set_reset(var, reset);
    return trans;
}

const DiscreteIOTransition&
HybridIOAutomaton::set_activation(DiscreteEvent event,
                                  DiscreteState source,
                                  const RealExpression& activation)
{
    if(!this->has_transition(event,source)) {
        ARIADNE_FAIL_MSG("The automaton " << this->_name << " has no transition with event "
            << event << " and source " << source << ".");
    }

    ARIADNE_ASSERT_MSG(!this->has_input_event(event) || identical(activation,1.0),
        "Error in transition " << event << " from " << source << " of automaton " << this->name() <<
        ": transitions labelled with input events cannot have an activation different from true.");
    
    DiscreteIOTransition& trans = this->_transition(event,source);
    trans.set_activation(activation);
    return trans;
}

void 
HybridIOAutomaton::set_grid(const std::map< RealVariable, Float>& grid) 
{
    // Check if grid is consistent with the variables of the automaton
    // grid spacing of input variables cannot be specified 
    ARIADNE_ASSERT_MSG(disjoint(this->_input_vars, grid),
        "Grid " << grid << " is not disjoint from the set of input variables of automaton " << this->name());
    // check if all variables of the grid are controlled by the automaton
    ARIADNE_ASSERT_MSG(subset(grid, this->controlled_vars()),
        "Grid " << grid << " contains variables that do not belongs to automaton " << this->name());
    
    this->_grid = grid;
}

void 
HybridIOAutomaton::set_grid(const RealVariable& var, Float scaling)
{
    // Variable var should be a controlled variable of the automaton.
    ARIADNE_ASSERT_MSG(!contains(this->_input_vars, var),
        "Variable " << var << " is an input variable for automaton " << this->name() <<
        ": grid scaling cannot be specified.");
    ARIADNE_ASSERT_MSG(contains(this->controlled_vars(), var),
        "Variable " << var << " does not belong to automaton " << this->name() <<
        ": grid scaling cannot be specified.");
    
    this->_grid[var] = scaling;
}

const std::map< RealVariable, Float >&
HybridIOAutomaton::grid() const 
{
    return this->_grid;
}
    
Float 
HybridIOAutomaton::grid(const RealVariable& var) const
{
    std::map< RealVariable, Float >::const_iterator iter = this->_grid.find(var);
    
    if(iter == this->_grid.end())   // If spacing is not specified, return the default value of 1.0
        return 1.0;
        
    return iter->second;
}

bool 
HybridIOAutomaton::has_input_var(const RealVariable& u) const
{
    return (this->_input_vars.find(u) != this->_input_vars.end());
}

bool 
HybridIOAutomaton::has_output_var(const RealVariable& y) const
{
    return (this->_output_vars.find(y) != this->_output_vars.end());
}

bool 
HybridIOAutomaton::has_internal_var(const RealVariable& x) const
{
    return (this->_internal_vars.find(x) != this->_internal_vars.end());
}

bool 
HybridIOAutomaton::has_input_event(const DiscreteEvent& e) const
{
    return (this->_input_events.find(e) != this->_input_events.end());
}

bool 
HybridIOAutomaton::has_output_event(const DiscreteEvent& e) const
{
    return (this->_output_events.find(e) != this->_output_events.end());
}


bool 
HybridIOAutomaton::has_internal_event(const DiscreteEvent& e) const
{
    return (this->_internal_events.find(e) != this->_internal_events.end());
}



bool
HybridIOAutomaton::has_mode(DiscreteState state) const
{
    // FIXME: This is a hack since we use std::list which cannot be searched by id.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return true;
            }
        }
    return false;
}


bool
HybridIOAutomaton::has_transition(DiscreteEvent event, DiscreteState source) const
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


const std::list< DiscreteIOMode >&
HybridIOAutomaton::modes() const
{
    return this->_modes;
}

DiscreteIOMode&
HybridIOAutomaton::_mode(DiscreteState state)
{
    // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
    for(std::list< DiscreteIOMode >::iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return *mode_iter;
            }
        }
    ARIADNE_FAIL_MSG("The automaton " << this->name() << " does not have a mode with id " << state);
}



const DiscreteIOMode&
HybridIOAutomaton::mode(DiscreteState state) const
{
    // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==state) {
                return *mode_iter;
            }
        }
    ARIADNE_FAIL_MSG("The automaton " << this->name() << " does not have a mode with id " << state);
}


const std::list< DiscreteIOTransition >&
HybridIOAutomaton::transitions() const
{
    return this->_transitions;
}



std::list< DiscreteIOTransition >
HybridIOAutomaton::transitions(DiscreteState source) const
{
    std::list< DiscreteIOTransition > result;
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->source()==source) {
                result.push_back(*transition_iter);
            }
        }
    return result;
}


const DiscreteIOTransition&
HybridIOAutomaton::transition(DiscreteEvent event, DiscreteState source) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return *transition_iter;
            }
        }
    ARIADNE_THROW(std::runtime_error, "HybridIOAutomaton::transition(event, source)",
        "The hybrid automaton does not have a transition with event \"" << event << "\" and source \"" << source << "\".");
}

DiscreteIOTransition&
HybridIOAutomaton::_transition(DiscreteEvent event, DiscreteState source)
{
    for(std::list< DiscreteIOTransition >::iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return *transition_iter;
            }
        }
    ARIADNE_THROW(std::runtime_error, "HybridIOAutomaton::transition(event, source)",
        "The hybrid automaton does not have a transition with event \"" << event << "\" and source \"" << source << "\".");
}


const std::string&
HybridIOAutomaton::name() const
{
    return this->_name;
}

const std::set< RealVariable >& 
HybridIOAutomaton::input_vars() const
{
    return this->_input_vars;
}

const std::set< RealVariable >& 
HybridIOAutomaton::output_vars() const
{
    return this->_output_vars;
}

const std::set< RealVariable >& 
HybridIOAutomaton::internal_vars() const
{
    return this->_internal_vars;
}

std::set< RealVariable >
HybridIOAutomaton::controlled_vars() const
{
    std::set< RealVariable > result = this->_internal_vars;
    result.insert(this->_output_vars.begin(), this->_output_vars.end());
    
    return result;
}


const std::set< DiscreteEvent >& 
HybridIOAutomaton::input_events() const
{
    return this->_input_events;
}

const std::set< DiscreteEvent >& 
HybridIOAutomaton::output_events() const
{
    return this->_output_events;
}

const std::set< DiscreteEvent >& 
HybridIOAutomaton::internal_events() const
{
    return this->_internal_events;
}

std::set< DiscreteEvent >
HybridIOAutomaton::controlled_events() const
{
    std::set< DiscreteEvent > result = this->_internal_events;
    result.insert(this->_output_events.begin(), this->_output_events.end());
    
    return result;
}


std::pair< HybridAutomaton, RealSpace > make_monolithic_automaton(const HybridIOAutomaton& hioa) 
{
    // Check if the input hioa is closed: no input variables or events
    if(!hioa.input_vars().empty()) {
        ARIADNE_FAIL_MSG("The set of input variables of I/O automaton " << hioa.name() << 
            " is not empty: cannot convert to a monolithic automaton.");
    }
    if(!hioa.input_events().empty()) {
        ARIADNE_FAIL_MSG("The set of input events of I/O automaton " << hioa.name() << 
            " is not empty: cannot convert to a monolithic automaton.");
    }
    
    // Make the list of all RealVariable of the automaton
    List<RealVariable> varlist;
    // first the internal variables...
    std::set< RealVariable >::const_iterator variter;
    for(variter = hioa.internal_vars().begin(); variter != hioa.internal_vars().end(); variter++)
    {
        varlist.append(*variter);
    }
    // ... then the output variables
    for(variter = hioa.output_vars().begin(); variter != hioa.output_vars().end(); variter++)
    {
        varlist.append(*variter);
    }
    // create the RealSpace corresponding to varlist
    RealSpace spc(varlist);
    
    // Create a monolithic automaton with the same name.
    HybridAutomaton ha(hioa.name());
    
    // Scan all modes of the hioa
    std::list< DiscreteIOMode > mode_list = hioa.modes();
    for(std::list< DiscreteIOMode >::iterator modeiter=mode_list.begin();
        modeiter != mode_list.end(); modeiter++)
    {
        // Create the expression list for the dynamics
        List<RealExpression> exprlist;
        // for each variable in varlist, add the expression defining the dynamics
        for(List<RealVariable>::const_iterator viter=varlist.begin(); viter != varlist.end(); viter++)
        {
            exprlist.append(modeiter->dynamics(*viter));
        }
        // Create a VectorFunction for the dynamics
        VectorFunction dyn(exprlist,spc);
        DiscreteState loc = modeiter->location();
        ha.new_mode(loc,dyn);
        // List all invariants
        std::list< RealExpression > invlist = modeiter->invariants();
        // Add all invariants to the mode in the monolithic automaton
        for(std::list< RealExpression >::const_iterator inviter=invlist.begin(); 
            inviter != invlist.end(); inviter++)
        {
            ha.new_invariant(loc, ScalarFunction(*inviter, spc));
        }
    }

    // Scan all transitions of the hioa
    std::list< DiscreteIOTransition > tr_list = hioa.transitions();
    for(std::list< DiscreteIOTransition >::iterator triter=tr_list.begin();
        triter != tr_list.end(); triter++)
    {
        // Create the expression list for the reset function
        List<RealExpression> exprlist;
        // for each variable in varlist, add the expression defining the reset
        for(List<RealVariable>::const_iterator viter=varlist.begin(); viter != varlist.end(); viter++)
        {
            exprlist.append(triter->reset(*viter));
        }
        // Create a VectorFunction for the reset
        VectorFunction res(exprlist,spc);
        // Add the transition to the monolithic automaton
        ha.new_transition(triter->event(), triter->source(), triter->target(),
                          res, ScalarFunction(triter->activation(), spc), triter->forced());
    }
    
    // Set the grid
    Vector<Float> lengths(spc.size());
    for(uint i = 0 ; i < spc.size() ; i++) {
        lengths[i] = hioa.grid(spc[i]);
    }
    Grid grid(lengths);
    ha.set_grid(grid);
    
    return make_pair(ha, spc);
}

//
// Auxiliary function that recursively compose ha1 and ha2.
//
DiscreteState _recursive_composition(HybridIOAutomaton& ha, 
                                     const HybridIOAutomaton& ha1, 
                                     const HybridIOAutomaton& ha2,
                                     const DiscreteState& init1,
                                     const DiscreteState& init2)
{
    // The composed location is named init1,init2
    DiscreteState newloc(init1.name() + "," + init2.name());

    if(ha.has_mode(newloc))  // The composed location init1,init2 already exists, return.
        return newloc;    
        
    // Location init1,init2 do not exist, create it.
    const DiscreteIOMode& loc1 = ha1.mode(init1);
    const DiscreteIOMode& loc2 = ha2.mode(init2);    
    ha.new_mode(newloc);
    
    // Add the dynamics of loc1.
    for(std::map< RealVariable, RealExpression >::const_iterator iter = loc1.dynamics().begin();
        iter != loc1.dynamics().end() ; iter++)
    {
        ha.set_dynamics(newloc, iter->first, iter->second);
    }    
    // Add the dynamics of loc2.
    for(std::map< RealVariable, RealExpression >::const_iterator iter = loc2.dynamics().begin();
        iter != loc2.dynamics().end() ; iter++)
    {
        ha.set_dynamics(newloc, iter->first, iter->second);
    }

    // Add the invariants of loc1.
    for(std::list< RealExpression >::const_iterator iter = loc1.invariants().begin();
        iter != loc1.invariants().end() ; iter++)
    {
        ha.new_invariant(newloc, *iter);
    }    
    // Add the invariants of loc2.
    for(std::list< RealExpression >::const_iterator iter = loc2.invariants().begin();
        iter != loc2.invariants().end() ; iter++)
    {
        ha.new_invariant(newloc, *iter);
    }    
    
    // Scan all transition exiting from loc1
    std::list< DiscreteIOTransition > transitions = ha1.transitions(init1);
    for(std::list< DiscreteIOTransition >::const_iterator iter = transitions.begin();
        iter != transitions.end() ; iter++)
    {
        DiscreteState target2 = init2;
        RealExpression act = iter->activation();
        std::map< RealVariable, RealExpression > res = iter->reset();
        bool forced=iter->forced();

        // std::cout << "Checking transition " << *iter << "..." << std::endl;
        // We should distinguish between different cases
        if( ha1.has_input_event(iter->event()) &&        // input event shared with ha2
            (ha2.has_input_event(iter->event()) || ha2.has_output_event(iter->event())) )
        {   // event is shared with ha2: components should synchronize
           // std::cout << "Event is input and shared with ha2." << std::endl;
           DiscreteIOTransition tr2;
            try {
                tr2 = ha2.transition(iter->event(), init2);
            }   
            catch(std::runtime_error e) {
                // std::cout << "No transition in ha2 with the same event." << std::endl;                    
                continue;
            }
            // the target state in the composed automaton is iter->target(),tr2.target()
            target2 = tr2.target();
            // activation depends on the second component
            act = tr2.activation();
            // join the reset functions
            res.insert(tr2.reset().begin(), tr2.reset().end());
            // the transition is forced only if it is forced in the second component
            forced = tr2.forced();
        } 
        else if( ha1.has_output_event(iter->event()) &&     // output event shared with ha2
                 ha2.has_input_event(iter->event())) 
        {   // event is shared with ha2: components should synchronize
            // std::cout << "Event is output and shared with ha2." << std::endl;
            DiscreteIOTransition tr2;
            try {
                tr2 = ha2.transition(iter->event(), init2);
            }   
            catch(std::runtime_error e) {
                // std::cout << "No transition in ha2 with the same event." << std::endl;                    
                continue;
            }
            // the target state in the composed automaton is iter->target(),tr2.target()
            target2 = tr2.target();
            // join the reset functions
            res.insert(tr2.reset().begin(), tr2.reset().end());
        } else {    // In all other cases the event is not shared with ha2
            // The reset function for the variables controlled by ha2 is the identity
            for(std::set<RealVariable>::const_iterator viter=ha2.internal_vars().begin() ;
                viter != ha2.internal_vars().end() ; viter++ )
            {
                res[*viter] = *viter;
            }
            for(std::set<RealVariable>::const_iterator viter=ha2.output_vars().begin() ;
                viter != ha2.output_vars().end() ; viter++ )
            {
                res[*viter] = *viter;
            }
        }
        // recursively create the target mode (if it does not exists)
        DiscreteState newtarget = _recursive_composition(ha, ha1, ha2, iter->target(), target2);
        // add the new transition to the automaton
        ha.new_transition(iter->event(), newloc, newtarget, res, act, forced);
        
    }   // end scanning transitions from loc1
        
    // Scan all remaining transition exiting from loc2: 
    transitions = ha2.transitions(init2);
    for(std::list< DiscreteIOTransition >::const_iterator iter = transitions.begin();
        iter != transitions.end() ; iter++)
    {
        // std::cout << "Checking transition " << *iter << "..." << std::endl;
        // remember that transition with shared events have already been added in the previous loop
        if(ha1.has_input_event(iter->event()) || ha1.has_output_event(iter->event()))
        {   // event is shared with ha1: transition have already been added, skip
            continue;
        }
        // The event is not shared with ha1: the first component should not react
        RealExpression act = iter->activation();
        std::map< RealVariable, RealExpression > res = iter->reset();
        // The reset function for the variables controlled by ha1 is the identity
        for(std::set<RealVariable>::const_iterator viter=ha1.internal_vars().begin() ;
            viter != ha1.internal_vars().end() ; viter++ )
        {
            res[*viter] = *viter;
        }
        for(std::set<RealVariable>::const_iterator viter=ha1.output_vars().begin() ;
            viter != ha1.output_vars().end() ; viter++ )
        {
            res[*viter] = *viter;
        }
        bool forced=iter->forced();
        // recursively create the target mode (if it does not exists)
        DiscreteState newtarget = _recursive_composition(ha, ha1, ha2, init1, iter->target());
        // add the new transition to the automaton
        ha.new_transition(iter->event(), newloc, newtarget, res, act, forced);
    }   // end scanning transitions from loc2    

    return newloc;
}

HybridIOAutomaton compose(const std::string& name, 
                          const HybridIOAutomaton& ha1, 
                          const HybridIOAutomaton& ha2,
                          const DiscreteState& init1,
                          const DiscreteState& init2)
{
    // 
    // FIRST STEP: check if the two HIOAs are compatible
    //
    // Internal vars of ha1 should be disjoint from
    //  the internal vars of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.internal_vars(),ha2.internal_vars()),
        "Cannot compose " << ha1.name() << " with " << ha2.name() << 
        ": internal variables are not disjoint.");
    //  the input vars of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.internal_vars(),ha2.input_vars()),
        "Cannot compose: internal variables of " << ha1.name() << 
        " are not disjoint from the input variables of " << ha2.name() << ".");
    //  the output vars of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.internal_vars(),ha2.output_vars()),
        "Cannot compose: internal variables of " << ha1.name() << 
        " are not disjoint from the output variables of " << ha2.name() << ".");
    // Output vars of ha1 should be disjoint from
    //  the output vars of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.output_vars(),ha2.output_vars()),
        "Cannot compose " << ha1.name() << " with " << ha2.name() << 
        ": output variables are not disjoint.");
    //  the internal vars of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.output_vars(),ha2.internal_vars()),
        "Cannot compose: output variables of " << ha1.name() << 
        " are not disjoint from the internal variables of " << ha2.name() << ".");
     // Input vars of ha1 should be disjoint from the internal vars of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.input_vars(),ha2.internal_vars()),
        "Cannot compose: input variables of " << ha1.name() << 
        " are not disjoint from the internal variables of " << ha2.name() << ".");
    //
    // Internal events of ha1 should be disjoint from
    //  the internal events of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.internal_events(),ha2.internal_events()),
        "Cannot compose " << ha1.name() << " with " << ha2.name() << 
        ": internal events are not disjoint.");
    //  the input events of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.internal_events(),ha2.input_events()),
        "Cannot compose: internal events of " << ha1.name() << 
        " are not disjoint from the input events of " << ha2.name() << ".");
    //  the output events of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.internal_events(),ha2.output_events()),
        "Cannot compose: internal events of " << ha1.name() << 
        " are not disjoint from the output events of " << ha2.name() << ".");
    // Output events of ha1 should be disjoint from
    //  the output events of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.output_events(),ha2.output_events()),
        "Cannot compose " << ha1.name() << " with " << ha2.name() << 
        ": output events are not disjoint.");
    //  the internal events of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.output_events(),ha2.internal_events()),
        "Cannot compose: output events of " << ha1.name() << 
        " are not disjoint from the internal events of " << ha2.name() << ".");
     // Input events of ha1 should be disjoint from the internal events of ha2
    ARIADNE_ASSERT_MSG(disjoint(ha1.input_events(),ha2.internal_events()),
        "Cannot compose: input events of " << ha1.name() << 
        " are not disjoint from the internal events of " << ha2.name() << ".");
    
    // The internal variables of the composed automaton is the union of 
    // the internal variables of the components
    std::set< RealVariable > internal_vars;
    internal_vars.insert(ha1.internal_vars().begin(), ha1.internal_vars().end());
    internal_vars.insert(ha2.internal_vars().begin(), ha2.internal_vars().end());
    // The output variables of the composed automaton is the union of 
    // the output variables of the components
    std::set< RealVariable > output_vars;
    output_vars.insert(ha1.output_vars().begin(), ha1.output_vars().end());
    output_vars.insert(ha2.output_vars().begin(), ha2.output_vars().end());
    // The input variables of the composed automaton is the union of 
    // the input variables of the components minus the output variables.
    std::set< RealVariable > input_vars1, input_vars2, input_vars;
    std::set_difference(ha1.input_vars().begin(), ha1.input_vars().end(), 
                        ha2.output_vars().begin(), ha2.output_vars().end(), 
                        std::inserter(input_vars1, input_vars1.begin()));
    std::set_difference(ha2.input_vars().begin(), ha2.input_vars().end(), 
                        ha1.output_vars().begin(), ha1.output_vars().end(), 
                        std::inserter(input_vars2, input_vars2.begin()));
    input_vars.insert(input_vars1.begin(), input_vars1.end());
    input_vars.insert(input_vars2.begin(), input_vars2.end());

    // The internal events of the composed automaton is the union of 
    // the internal events of the components
    std::set< DiscreteEvent > internal_events;
    internal_events.insert(ha1.internal_events().begin(), ha1.internal_events().end());
    internal_events.insert(ha2.internal_events().begin(), ha2.internal_events().end());
    // The output events of the composed automaton is the union of 
    // the output events of the components
    std::set< DiscreteEvent > output_events;
    output_events.insert(ha1.output_events().begin(), ha1.output_events().end());
    output_events.insert(ha2.output_events().begin(), ha2.output_events().end());
    // The input events of the composed automaton is the union of 
    // the input events of the components minus the output events.
    std::set< DiscreteEvent > input_events1, input_events2, input_events;
    std::set_difference(ha1.input_events().begin(), ha1.input_events().end(), 
                        ha2.output_events().begin(), ha2.output_events().end(), 
                        std::inserter(input_events1, input_events1.begin()));
    std::set_difference(ha2.input_events().begin(), ha2.input_events().end(), 
                        ha1.output_events().begin(), ha1.output_events().end(), 
                        std::inserter(input_events2, input_events2.begin()));
    input_events.insert(input_events1.begin(), input_events1.end());
    input_events.insert(input_events2.begin(), input_events2.end());

    // Create an HybridIOAutomaton with the given name
    HybridIOAutomaton ha(name, 
                         input_vars, output_vars, internal_vars,
                         input_events, output_events, internal_events);
    
    // Set the grid as the union of the grids of the two components
    std::map< RealVariable, Float > grid = ha1.grid();
    grid.insert(ha2.grid().begin(), ha2.grid().end());
    ha.set_grid(grid);
    
    // Start the recursive composition from the initial states init1 and init2.
    _recursive_composition(ha, ha1, ha2, init1, init2);
    
    return ha;

}

HybridIOAutomaton aasap_relaxation(const HybridIOAutomaton& hioa)
{
	typedef std::map<DiscreteEvent,bool> ReceivedEventsMap;
	typedef std::list<std::pair<DiscreteState,ReceivedEventsMap> > LocationsWithReceivedEvents;

	// Create the Delta constant with arbitrary value 0
	RealConstant Delta("Delta",0.0);
	// Create the d variable
	RealVariable d("d");

	// Create the target automaton with a copy of the original variables and events
	HybridIOAutomaton aasap(hioa.name()+"_aasap",
							hioa.input_vars(),
							hioa.output_vars(),
							hioa.internal_vars(),
							hioa.input_events(),
							hioa.output_events(),
							hioa.internal_events());

	// Add the AASAP internal variables and events
		
		// For each input event
		for (std::set<DiscreteEvent>::const_iterator input_event_it = aasap.input_events().begin(); input_event_it != aasap.input_events().end(); input_event_it++)
		{
			// Clock for waiting for acknowledgement of each input event
			aasap.add_internal_var(RealVariable("y_" + input_event_it->name()));
			// The corresponding input acknowledgment event
			aasap.add_internal_event(DiscreteEvent(input_event_it->name() + "_ack"));			
		}
		// Variable for remaining in the same location
		aasap.add_internal_var(d);

	// Get the input events, internal events and output events of the automaton
	const std::set<DiscreteEvent>& input_events = aasap.input_events();
	const std::set<DiscreteEvent>& internal_events = aasap.internal_events();
	const std::set<DiscreteEvent>& output_events = aasap.output_events();

	// Create the destination locations with the information about received events, and copy the original modes into the aasap ones

		// Initialize the list of locations with received events information
		LocationsWithReceivedEvents lwre;

		// For each mode of the original automaton
		for (std::list<DiscreteIOMode>::const_iterator mode_it = hioa.modes().begin(); mode_it != hioa.modes().end(); mode_it++)
		{
			// We must construct 2^num_input_events locations
			for (int i=0;i < (1<<input_events.size());i++)
			{
				// Initialize the state name
				std::string statename = mode_it->location().name();

				// Initialize the received events map
				ReceivedEventsMap rem;

				// Initialize the iterator and counter
				std::set<DiscreteEvent>::const_iterator input_event_it = input_events.begin();
				int j=0; 

				// For each event
				while (input_event_it != input_events.end())
				{
					// Stores the information of being received
					bool is_received = (((1<<j) & i) > 0);

					// If received, adds the input event name to the state name
					if (is_received)
						statename += "_" + input_event_it->name();

					// Insert the information about whether the event is considered received
					rem.insert(make_pair<DiscreteEvent,bool>(*input_event_it,is_received));

					// Increase the iterator and the counter
					input_event_it++; 
					j++;
				}

				// Add the corresponding mode to the AASAP automaton, along with the original dynamics
				aasap.new_mode(DiscreteState(statename),mode_it->dynamics());

				// Add the original location with received events
				lwre.push_back(make_pair<DiscreteState,ReceivedEventsMap>(mode_it->location(),rem));
			}
		}

	// For each location of the AASAP automaton
	for (LocationsWithReceivedEvents::const_iterator lwre_it = lwre.begin(); lwre_it != lwre.end(); lwre_it++)
	{
		// Get the AASAP name of the location
			
			// Initial value
			std::string location_aasap_name = lwre_it->first.name();
			// For each input event, if received 
			for (ReceivedEventsMap::const_iterator rem_it = lwre_it->second.begin(); rem_it != lwre_it->second.end(); rem_it++)
				if (rem_it->second)
					location_aasap_name += "_" + rem_it->first.name();

		// Hold the aasap location

			DiscreteState aasap_location(location_aasap_name);

		// Hold the received events map

			ReceivedEventsMap rem = lwre_it->second;

		// Insert the dynamics of the AASAP internal variables

			// The expression for a running clock
			RealExpression clock_running = 1.0;
			// The expression for an idle clock
			RealExpression clock_idle = 0.0;

			// For d, the dynamics is always to run
			aasap.set_dynamics(aasap_location, d, clock_running);

			// For the y variables, it depends whether the corresponding input event has been received or not
			for (ReceivedEventsMap::const_iterator rem_it = rem.begin(); rem_it != rem.end(); rem_it++)
				aasap.set_dynamics(aasap_location, RealVariable("y_"+rem_it->first.name()), rem_it->second ? clock_running : clock_idle);

		// Insert the transitions due to input events
			
			// For each input event
			for (std::set<DiscreteEvent>::const_iterator input_event_it = input_events.begin(); input_event_it != input_events.end(); input_event_it++)
			{
				// Identify the target location: the same if the event has already been received, the one with the received event otherwise; 
				// in the second case, a reset of the corresponding clock is required
				std::string target_location;
			
				// If the event has already been received, create a loop transition, otherwise create a transition to the proper location, resetting the related clock
				if (rem[*input_event_it])
				{
					// Create the reset (the identity)
					std::map< RealVariable, RealExpression> reset; 
					for (std::set<RealVariable>::const_iterator var_it = aasap.internal_vars().begin(); var_it != aasap.internal_vars().end(); var_it++)
						reset[*var_it] = *var_it;
					
					aasap.new_unforced_transition(*input_event_it,aasap_location,aasap_location,reset);					
				}
				else
				{
					// Create the reset (reset to zero only the variable corresponding to the event)
					std::map< RealVariable, RealExpression> reset; 
					for (std::set<RealVariable>::const_iterator var_it = aasap.internal_vars().begin(); var_it != aasap.internal_vars().end(); var_it++)
					{
						if (var_it->name() != "y_"+input_event_it->name())
							reset[*var_it] = *var_it;
						else
							reset[*var_it] = 0.0;
					}

					// Initialize the target location name
					std::string target_location_name = lwre_it->first.name();
					// Copy the received events map information
					ReceivedEventsMap target_rem = rem;
					// Set the received event flag as true
					target_rem[*input_event_it] = true;
					// Build the name by appending the received events
					for (ReceivedEventsMap::const_iterator rem_it = target_rem.begin(); rem_it != target_rem.end(); rem_it++)
						if (rem_it->second)
							target_location_name += "_" + rem_it->first.name();

					// Create the transition
					aasap.new_unforced_transition(*input_event_it,aasap_location,DiscreteState(target_location_name),reset);					
				}
			}

		// Insert the transitions due to internal and output events

			// Join the two sets
			std::set<DiscreteEvent> internal_output_events = internal_events;
			internal_output_events.insert(output_events.begin(),output_events.end());

			// For each event
			for (std::set<DiscreteEvent>::const_iterator event_it = internal_output_events.begin(); event_it != internal_output_events.end(); event_it++)
			{
				// Initialize the flags that inform whether this is an input acknowledgement event and whether the corresponding input event has been received
				bool is_input_acknowledgement_event = false;
				bool has_input_event_been_received = false;

				// Initialize the base event (used for storing the base event for input acknowledgement events)
				DiscreteEvent base_event = *event_it;
	
				// Verify if this is an input acknowledgement event
				for (std::set<DiscreteEvent>::const_iterator input_event_it = input_events.begin(); input_event_it != input_events.end(); input_event_it++)
					if (event_it->name() == input_event_it->name()+"_ack")
					{
						// Set the flag
						is_input_acknowledgement_event = true;
						// Set if the related input event has been received
						has_input_event_been_received = rem[*input_event_it];
						// Store the base event
						base_event = *input_event_it;
						// No more event checking is necessary
						break;
					}

				// Consider it only if it IS NOT an input acknowledgment event, or (if it IS) if the corresponding input event HAS been received
				if (!is_input_acknowledgement_event || has_input_event_been_received)
				{
					// For each transition of the original automaton... 
					for (std::list<DiscreteIOTransition>::const_iterator trans_it = hioa.transitions().begin(); trans_it != hioa.transitions().end(); trans_it++)
					{
						// ...featuring the same base event and the same base source location
						if (trans_it->event().name() == base_event.name() && trans_it->source().name() == lwre_it->first.name())
						{
							// Initialize the target location name as the original target name
							std::string target_location_name = trans_it->target().name();
							// Copy the received events map information
							ReceivedEventsMap target_rem = rem;
							// If this is an input acknowledgment event, set the corresponding target received information as false
							if (is_input_acknowledgement_event)
								target_rem[base_event] = false;
							// Build the name by appending the received events
							for (ReceivedEventsMap::const_iterator rem_it = target_rem.begin(); rem_it != target_rem.end(); rem_it++)
								if (rem_it->second)
									target_location_name += "_" + rem_it->first.name();

							// Create the reset (equal to the original reset, plus resetting d and keeping the previous value for the variables associated to input events)
							std::map<RealVariable,RealExpression> reset = trans_it->reset(); 
							reset[d] = 0.0;
							for (std::set<DiscreteEvent>::const_iterator input_event_it = input_events.begin(); input_event_it != input_events.end(); input_event_it++)
								reset[RealVariable("y_" + input_event_it->name())] = RealVariable("y_" + input_event_it->name());

							// If the activation is always true, do not explicitate it, otherwise enlarge it by Delta
							if (identical(trans_it->activation(),RealExpression(1.0)))							
								aasap.new_unforced_transition(*event_it,aasap_location,DiscreteState(target_location_name),reset);
							else
								aasap.new_unforced_transition(*event_it,aasap_location,DiscreteState(target_location_name),reset,trans_it->activation()+Delta);
						}
					}
				}
			}

		// Insert the invariants

			// There is an invariant for all the original transitions related to either:
			// a) input events that have been received in the considered location
			// b) internal and output events

			// For each transition of the original automaton...
			for (std::list<DiscreteIOTransition>::const_iterator trans_it = hioa.transitions().begin(); trans_it != hioa.transitions().end(); trans_it++)
			{
				// ...featuring the same source location as the base of the considered location
				if (trans_it->source().name() == lwre_it->first.name())
				{
					// Get the information of the event being present in the ReceivedEventMap (i.e. being an input event) of the considered location
					std::map<DiscreteEvent,bool>::iterator rem_it = rem.find(trans_it->event());

					// If the transition event is an input event
					if (rem_it != rem.end())
					{
						// If it has been received, add the related invariant
						if (rem_it->second)
						{
							// Prepare the expression
							RealExpression invariant = min(d-Delta,RealVariable("y_"+rem_it->first.name()) -Delta);
							// Add the invariant
							aasap.new_invariant(aasap_location,invariant);						
						}
					}
					// Otherwise it is an internal or output event
					else
					{
							// Prepare the expression
							RealExpression invariant = min(d-Delta,-left_delta_restrict(trans_it->activation(),true));
							// Add the invariant
							aasap.new_invariant(aasap_location,invariant);						
					}
				}
			}
	}

	// Returns
	return aasap;
}

bool is_elastic_controller(const HybridIOAutomaton& hioa)
{
	/*
		1) Only clocks
		2) Only internal variables
		3) All input events have to be handled in any location
		4) Only forced transitions for internal and output events
		5) Resets only to zero
		6) Guards on a single clock and no >=0.0 guards
	*/

	// If it passed all checks, return true
	return true;
}

std::ostream&
operator<<(std::ostream& os, const HybridIOAutomaton& ha)
{
    return os << "HybridIOAutomaton( name=" << ha.name() << 
        ", input vars=" << ha.input_vars() << ", output vars=" << ha.output_vars() <<
        ", internal vars=" << ha.internal_vars() <<
        ", input events=" << ha.input_events() << ", output events=" << ha.output_events() <<
        ", internal events=" << ha.internal_events() <<        
        ", modes=" << ha.modes() << ", transitions=" << ha.transitions() <<
        ", grid=" << ha.grid() << ")";
}

void 
HybridIOAutomaton::substitute(const Constant<Real>& con, const Real& c)
{
	for (std::list<DiscreteIOMode>::iterator modes_it=this->_modes.begin();modes_it!=this->_modes.end();modes_it++)
		modes_it->substitute(con,c);
	for (std::list<DiscreteIOTransition>::iterator trans_it=this->_transitions.begin();trans_it!=this->_transitions.end();trans_it++)
		trans_it->substitute(con,c);
}


}   // namespace Ariadne
