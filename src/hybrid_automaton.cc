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
#include "function_interface.h"
#include "hybrid_automaton.h"

namespace Ariadne {  

typedef uint DimensionType;

class HybridSet {};

class HybridSpace : public std::map<DiscreteState,DimensionType> {};

uint 
DiscreteMode::
dimension() const 
{ 
  return this->_dynamic->argument_size(); 
}
    

DiscreteMode::
DiscreteMode(DiscreteState location,
             const FunctionInterface& dynamic)
  :  _location(location), _dynamic(dynamic.clone()), _invariants() 
{
  ARIADNE_ASSERT(dynamic.result_size()==dynamic.argument_size());
}


DiscreteMode::
DiscreteMode(DiscreteState location,
             const boost::shared_ptr< const FunctionInterface > dynamic, 
             const std::vector< boost::shared_ptr< const FunctionInterface > >& invariants)
  :  _location(location), _dynamic(dynamic), _invariants(invariants) 
{
  ARIADNE_ASSERT(dynamic->result_size()==dynamic->argument_size());
  //ARIADNE_ASSERT(dynamic->argument_size()==invariants->argument_size());
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
                   const FunctionInterface& reset,
                   const FunctionInterface& activation)
  : _event(event), 
    _source(&source), 
    _target(&target), 
    _activation(activation.clone()), 
    _reset(reset.clone()) 
{ 
  ARIADNE_ASSERT(activation.result_size()==1);
  ARIADNE_ASSERT(activation.argument_size()==source.dimension());
  ARIADNE_ASSERT(reset.argument_size()==source.dimension());
  ARIADNE_ASSERT(reset.result_size()==target.dimension());
}


DiscreteTransition::
DiscreteTransition(DiscreteEvent event,
                   const DiscreteMode& source, 
                   const DiscreteMode& target,
                   const boost::shared_ptr< FunctionInterface > reset,
                   const boost::shared_ptr< FunctionInterface > activation) 
  : _event(event), _source(&source), _target(&target), 
    _activation(activation), _reset(reset) 
{ 
  ARIADNE_ASSERT(activation->result_size()==1);
  ARIADNE_ASSERT(activation->argument_size()==source.dimension());
  ARIADNE_ASSERT(reset->argument_size()==source.dimension());
  ARIADNE_ASSERT(reset->result_size()==target.dimension());
}


DiscreteTransition::
DiscreteTransition(DiscreteEvent event,
                   const boost::shared_ptr< DiscreteMode > source, 
                   const boost::shared_ptr< DiscreteMode > target,
                   const boost::shared_ptr< FunctionInterface > reset,
                   const boost::shared_ptr< FunctionInterface > activation) 
  : _event(event), _source(&*source), _target(&*target), 
    _activation(activation), _reset(reset) 
{ 
  ARIADNE_ASSERT(activation->result_size()==1);
  ARIADNE_ASSERT(activation->argument_size()==source->dimension());
  ARIADNE_ASSERT(reset->argument_size()==source->dimension());
  ARIADNE_ASSERT(reset->result_size()==target->dimension());
}


std::ostream& 
operator<<(std::ostream& os, const DiscreteTransition& transition)  
{ 
  return os << "DiscreteTransition( "
            << "event=" << transition.event() << ", " 
            << "source=" << transition.source().location() << ", "
            << "target=" << transition.target().location() << ", "
            << "reset=" << transition.reset() << ", " 
            << "activation=" << transition.activation() << " )";
}




HybridAutomaton::HybridAutomaton(const std::string& name)
  : _name(name) 
{ 
}





const DiscreteMode& 
HybridAutomaton::new_mode(DiscreteState location,
                          const FunctionInterface& dynamic) 
{
  ARIADNE_ASSERT(location>0);
  if(this->has_mode(location)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  this->_modes.insert(DiscreteMode(location,dynamic));
  return this->mode(location);
}


const DiscreteTransition& 
HybridAutomaton::new_transition(DiscreteEvent event,
               const DiscreteMode &source, 
               const DiscreteMode &target,
               const FunctionInterface &reset,
               const FunctionInterface &activation) 
{
  ARIADNE_ASSERT(event>0);
  DiscreteEvent event_id=event;
  DiscreteState source_id=source.location();
  DiscreteState target_id=target.location();
  if(this->has_transition(event_id,source_id)) {
    throw std::runtime_error("The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    throw std::runtime_error("The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(target_id)) {
    throw std::runtime_error("The desitination mode of the transition must be in the automaton");
  }
  
  const DiscreteMode& this_source=this->mode(source_id);
  const DiscreteMode& this_target=this->mode(target_id);
  if(&source!=&this_source) {
    throw std::runtime_error("The source mode of the transition is not equal to the mode in the automaton with the same id.");
  }
  if(&target!=&this_target) {
    throw std::runtime_error("The target mode of the transition is not equal to the mode in the automaton with the same id.");
  }
  
  this->_transitions.insert(DiscreteTransition(event_id,this_source,this_target,reset,activation));
  return this->transition(event_id,source_id);
}


const DiscreteTransition& 
HybridAutomaton::new_transition(DiscreteEvent event,
               DiscreteState source, 
               DiscreteState target,
               const FunctionInterface &reset,
               const FunctionInterface &activation) 
{
  ARIADNE_ASSERT(event>0);
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
  this->_transitions.insert(DiscreteTransition(event,source_mode,target_mode,reset,activation));
  return this->transition(event,source);
}


bool
HybridAutomaton::has_mode(DiscreteState state) const 
{
  // FIXME: This is a hack since we use std::set which cannot be searched by id.
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
HybridAutomaton::has_transition(DiscreteEvent event, DiscreteState source) const 
{
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->event()==event && transition_iter->source().location()==source) {
      return true;
    }
  }
  return false;
}



HybridSet
HybridAutomaton::invariant() const 
{
  ARIADNE_NOT_IMPLEMENTED;
}



HybridSpace
HybridAutomaton::state_space() const 
{
  HybridSpace result;
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin(); 
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    result[mode_iter->location()]=mode_iter->dimension();
  } 
  return result;
}


const std::set< DiscreteMode >& 
HybridAutomaton::modes() const 
{
  return this->_modes;
}


const DiscreteMode& 
HybridAutomaton::mode(DiscreteState state) const 
{
  // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->location()==state) {
      return *mode_iter;
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


const std::set< DiscreteTransition >& 
HybridAutomaton::transitions() const 
{
  return this->_transitions;
}



std::set< DiscreteTransition >
HybridAutomaton::transitions(DiscreteState source) const
{
  std::set< DiscreteTransition > result;
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->source().location()==source) {
      result.insert(*transition_iter);
    }
  }
  return result;
}



const DiscreteTransition& 
HybridAutomaton::transition(DiscreteEvent event, DiscreteState source) const 
{
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->event()==event && transition_iter->source().location()==source) {
      return *transition_iter;
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a transition with the given event and source.");
}


const std::string&
HybridAutomaton::name() const
{ 
  return this->_name; 
}


std::ostream& 
operator<<(std::ostream& os, const HybridAutomaton& ha) 
{
  return os << "HybridAutomaton( modes=" << ha.modes() << ", transitions=" << ha.transitions() << ")"; 
}




}
