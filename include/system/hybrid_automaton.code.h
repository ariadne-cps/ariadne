/***************************************************************************
 *            hybrid_automaton.code.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>

#include "base/stlio.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"

#include "hybrid_automaton.h"

namespace Ariadne {  
  

template<class R>
HybridAutomaton<R>::HybridAutomaton(const std::string &name)
  : _name(name) 
{ 
}


template<class R>
HybridAutomaton<R>*
HybridAutomaton<R>::clone() const {
  throw NotImplemented(__PRETTY_FUNCTION__);
  return new HybridAutomaton<R>(*this);
}

template<class R>
HybridAutomaton<R>::~HybridAutomaton() {
  this->_transitions.clear();
}


template<class R>
const DiscreteMode<R>& 
HybridAutomaton<R>::new_mode(DiscreteState id,
         const VectorField<R>& dynamic,
         const ConstraintSet<R>& invariant) 
{
  if(this->has_mode(id)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  this->_modes.insert(DiscreteMode<R>(id,dynamic,invariant));
  return this->mode(id);
}


template<class R>
const DiscreteTransition<R>& 
HybridAutomaton<R>::new_transition(DiscreteEvent event,
               const DiscreteMode<R> &source, 
               const DiscreteMode<R> &destination,
               const Map<R> &reset,
               const ConstraintSet<R> &activation) 
{
  DiscreteEvent event_id=event;
  DiscreteState source_id=source.discrete_state();
  DiscreteState destination_id=destination.discrete_state();
  if(this->has_transition(event_id,source_id)) {
    throw std::runtime_error("The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    throw std::runtime_error("The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(destination_id)) {
    throw std::runtime_error("The desitination mode of the transition must be in the automaton");
  }
  
  const DiscreteMode<R>& this_source=this->mode(source_id);
  const DiscreteMode<R>& this_destination=this->mode(destination_id);
  if(&source!=&this_source) {
    throw std::runtime_error("The source mode of the transition is not equal to the mode in the automaton with the same id.");
  }
  if(&destination!=&this_destination) {
    throw std::runtime_error("The destination mode of the transition is not equal to the mode in the automaton with the same id.");
  }
  
  this->_transitions.insert(DiscreteTransition<R>(event_id,this_source,this_destination,reset,activation));
  return this->transition(event_id,source_id);
}


template<class R>
const DiscreteTransition<R>& 
HybridAutomaton<R>::new_transition(DiscreteEvent event,
               DiscreteState source, 
               DiscreteState destination,
               const Map<R> &reset,
               const ConstraintSet<R> &activation) 
{
  if(this->has_transition(event,source)) {
    throw std::runtime_error("The automaton already has a transition with the given id and source id.");
  }
  if(!this->has_mode(source)) {
    throw std::runtime_error("The automaton does not contain a mode with ths given source id");
  }
  if(!this->has_mode(destination)) {
    throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
  }
  
  const DiscreteMode<R>& source_mode=this->mode(source);
  const DiscreteMode<R>& destination_mode=this->mode(destination);
  this->_transitions.insert(DiscreteTransition<R>(event,source_mode,destination_mode,reset,activation));
  return this->transition(event,source);
}


template<class R>
bool
HybridAutomaton<R>::has_mode(DiscreteState state) const 
{
  // FIXME: This is a hack since we use std::set which cannot be searched by id.
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->discrete_state()==state) {
      return true;
    }
  }
  return false;
}


template<class R>
bool 
HybridAutomaton<R>::has_transition(DiscreteEvent event, DiscreteState source) const 
{
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->discrete_event()==event && transition_iter->source().discrete_state()==source) {
      return true;
    }
  }
  return false;
}



template<class R>
HybridSet<R>
HybridAutomaton<R>::invariant() const 
{
  HybridSet<R> result;
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin(); 
      mode_iter!=this->_modes.end(); ++mode_iter)
  {
    result.new_location(mode_iter->discrete_state(),mode_iter->invariant());
  }
  return result;
}


template<class R>
HybridSpace
HybridAutomaton<R>::locations() const 
{
  HybridSpace result;
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin(); 
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    result.new_location(mode_iter->discrete_state(),mode_iter->dimension());
  } 
  return result;
}

template<class R>
HybridSpace
HybridAutomaton<R>::state_space() const 
{
  return this->locations();
}


template<class R>
const std::set< DiscreteMode<R> >& 
HybridAutomaton<R>::modes() const 
{
  return this->_modes;
}


template<class R>
const DiscreteMode<R>& 
HybridAutomaton<R>::mode(DiscreteState state) const 
{
  // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->discrete_state()==state) {
      return *mode_iter;
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


template<class R>
const std::set< DiscreteTransition<R> >& 
HybridAutomaton<R>::transitions() const 
{
  return this->_transitions;
}


template<class R>
reference_vector< const DiscreteTransition<R> >
HybridAutomaton<R>::transitions(DiscreteState source) const
{
  reference_vector< const DiscreteTransition<R> > result;
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->source().discrete_state()==source) {
      result.push_back(*transition_iter);
    }
  }
  return result;
}


template<class R>
const DiscreteTransition<R>& 
HybridAutomaton<R>::transition(DiscreteEvent event, DiscreteState source) const 
{
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->discrete_event()==event && transition_iter->source().discrete_state()==source) {
      return *transition_iter;
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a transition with the given event and source.");
}


template<class R>
const std::string&
HybridAutomaton<R>::name() const{ 
  return this->_name; 
}


template<class R>
std::ostream& 
HybridAutomaton<R>::write(std::ostream& os) const {
  return os << "HybridAutomaton( modes=" << this->_modes << ", transitions=" << this->_transitions << ")"; 
}




}
