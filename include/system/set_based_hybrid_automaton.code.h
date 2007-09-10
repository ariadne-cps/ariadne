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

#include "../base/stlio.h"
#include "../geometry/set_reference.h"
#include "../geometry/hybrid_set.h"

#include "set_based_hybrid_automaton.h"

namespace Ariadne {  
  

template<class R>
System::SetBasedHybridAutomaton<R>::SetBasedHybridAutomaton(const std::string &name)
  : _name(name) 
{ 
}


template<class R>
System::SetBasedHybridAutomaton<R>::~SetBasedHybridAutomaton() {
  this->_transitions.clear();
}


template<class R>
const System::SetBasedDiscreteMode<R>& 
System::SetBasedHybridAutomaton<R>::new_mode(id_type id,
         const VectorFieldInterface<R>& dynamic,
         const Geometry::SetInterface<R>& invariant) 
{
  if(this->has_mode(id)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  this->_modes.insert(SetBasedDiscreteMode<R>(id,dynamic,invariant));
  return this->mode(id);
}


template<class R>
const System::SetBasedDiscreteTransition<R>& 
System::SetBasedHybridAutomaton<R>::new_transition(id_type event_id,
               const SetBasedDiscreteMode<R> &source, 
               const SetBasedDiscreteMode<R> &destination,
               const MapInterface<R> &reset,
               const Geometry::SetInterface<R> &activation) 
{
  id_type source_id=source.id();
  id_type destination_id=destination.id();
  if(this->has_transition(event_id,source_id)) {
    throw std::runtime_error("The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    throw std::runtime_error("The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(destination_id)) {
    throw std::runtime_error("The desitination mode of the transition must be in the automaton");
  }
  
  const SetBasedDiscreteMode<R>& this_source=this->mode(source_id);
  const SetBasedDiscreteMode<R>& this_destination=this->mode(destination_id);
  if(&source!=&this_source) {
    throw std::runtime_error("The source mode of the transition is not equal to the mode in the automaton with the same id.");
  }
  if(&destination!=&this_destination) {
    throw std::runtime_error("The destination mode of the transition is not equal to the mode in the automaton with the same id.");
  }
  
  this->_transitions.insert(SetBasedDiscreteTransition<R>(event_id,this_source,this_destination,reset,activation));
  return this->transition(event_id,source_id);
}


template<class R>
const System::SetBasedDiscreteTransition<R>& 
System::SetBasedHybridAutomaton<R>::new_transition(id_type event_id,
               id_type source_id, 
               id_type destination_id,
               const MapInterface<R> &reset,
               const Geometry::SetInterface<R> &activation) 
{
  if(this->has_transition(event_id,source_id)) {
    throw std::runtime_error("The automaton already has a transition with the given id and source id.");
  }
  if(!this->has_mode(source_id)) {
    throw std::runtime_error("The automaton does not contain a mode with ths given source id");
  }
  if(!this->has_mode(destination_id)) {
    throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
  }
  
  const SetBasedDiscreteMode<R>& source=this->mode(source_id);
  const SetBasedDiscreteMode<R>& destination=this->mode(destination_id);
  this->_transitions.insert(SetBasedDiscreteTransition<R>(event_id,source,destination,reset,activation));
  return this->transition(event_id,source_id);
}


template<class R>
bool
System::SetBasedHybridAutomaton<R>::has_mode(id_type id) const 
{
  // FIXME: This is a hack since we use std::set which cannot be searched by id.
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->id()==id) {
      return true;
    }
  }
  return false;
}


template<class R>
bool 
System::SetBasedHybridAutomaton<R>::has_transition(id_type event_id, id_type source_id) const 
{
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->id()==event_id && transition_iter->source().id()==source_id) {
      return true;
    }
  }
  return false;
}



template<class R>
Geometry::HybridSet<R>
System::SetBasedHybridAutomaton<R>::invariant() const 
{
  Geometry::HybridSet<R> result;
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin(); 
      mode_iter!=this->_modes.end(); ++mode_iter)
  {
    result.new_location(mode_iter->id(),mode_iter->invariant());
  }
  return result;
}


template<class R>
Geometry::HybridSpace
System::SetBasedHybridAutomaton<R>::locations() const 
{
  Geometry::HybridSpace result;
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin(); 
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    result.new_location(mode_iter->id(),mode_iter->dimension());
  } 
  return result;
}


template<class R>
const std::set< System::SetBasedDiscreteMode<R> >& 
System::SetBasedHybridAutomaton<R>::modes() const 
{
  return this->_modes;
}


template<class R>
const System::SetBasedDiscreteMode<R>& 
System::SetBasedHybridAutomaton<R>::mode(id_type id) const 
{
  // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given id.
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->id()==id) {
      return *mode_iter;
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


template<class R>
const std::set< System::SetBasedDiscreteTransition<R> >& 
System::SetBasedHybridAutomaton<R>::transitions() const 
{
  return this->_transitions;
}


template<class R>
const System::SetBasedDiscreteTransition<R>& 
System::SetBasedHybridAutomaton<R>::transition(id_type event_id, id_type source_id) const 
{
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->id()==event_id && transition_iter->source().id()==source_id) {
      return *transition_iter;
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a transition with the given event_id and source_id.");
}


template<class R>
const std::string&
System::SetBasedHybridAutomaton<R>::name() const{ 
  return this->_name; 
}


template<class R>
std::ostream& 
System::SetBasedHybridAutomaton<R>::write(std::ostream& os) const {
  return os << "SetBasedHybridAutomaton( modes=" << this->_modes << ", transitions=" << this->_transitions << ")"; 
}




}
