/***************************************************************************
 *            constraint_hybrid_automaton.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
#include "../system/map.h"
#include "../system/vector_field.h"

#include "constraint_hybrid_automaton.h"



namespace Ariadne {  
  

template<class R>
System::ConstraintDiscreteMode<R>::ConstraintDiscreteMode(id_type id, 
                                                          const VectorFieldInterface<R>& dynamic, 
                                                          const std::vector< boost::shared_ptr< const ConstraintInterface<R> > >& invariant)
  : _id(id), _dynamic(dynamic.clone()), _invariant(invariant)
{
  for(typename std::vector< boost::shared_ptr< const ConstraintInterface<R> > >::const_iterator constraint_iter=invariant.begin();
      constraint_iter!=invariant.end(); ++constraint_iter)
  { 
    boost::shared_ptr< const ConstraintInterface<R> > constraint_ptr=*constraint_iter;
    const ConstraintInterface<R>& constraint=*constraint_ptr;
    ARIADNE_CHECK_EQUAL_DIMENSIONS(dynamic,constraint,"ConstraintDiscreteMode::ConstraintDiscreteMode(...)");
  }
}


template<class R>
id_type
System::ConstraintDiscreteMode<R>::id() const
{
  return this->_id;
}


template<class R>
dimension_type
System::ConstraintDiscreteMode<R>::dimension() const
{
  return this->_dynamic->dimension();
}


template<class R>
const System::VectorFieldInterface<R>&
System::ConstraintDiscreteMode<R>::dynamic() const
{
  return *this->_dynamic;
}


template<class R>
const System::ConstraintInterface<R>&
System::ConstraintDiscreteMode<R>::invariant(size_type k) const
{
  return *this->_invariant.at(k);
}


template<class R>
std::ostream&
System::ConstraintDiscreteMode<R>::write(std::ostream& os) const
{
  typedef typename std::set< boost::shared_ptr< const ConstraintInterface<R> > >::const_iterator constraint_iterator;
  return os << "DiscreteMode( id=" << this->id() 
            << ", dynamic=" << this->dynamic() << ", invariant={...} )";
}





template<class R>
System::ConstraintDiscreteTransition<R>::ConstraintDiscreteTransition(id_type id, 
                                                                      const ConstraintDiscreteMode<R>& source, 
                                                                      const ConstraintDiscreteMode<R>& destination,
                                                                      const MapInterface<R>& reset, 
                                                                      const ConstraintInterface<R>& activation,
                                                                      bool forced)
  : _event_id(id), _source(&source), _destination(&destination), _reset(reset.clone()), _activation(activation.clone()), _forced(forced)
{
}


template<class R>
id_type
System::ConstraintDiscreteTransition<R>::id() const
{
  return this->_event_id;
}


template<class R>
id_type
System::ConstraintDiscreteTransition<R>::source_id() const
{
  return this->_source->id();
}


template<class R>
const System::ConstraintDiscreteMode<R>&
System::ConstraintDiscreteTransition<R>::source() const
{
  return *this->_source;
}


template<class R>
const System::ConstraintDiscreteMode<R>&
System::ConstraintDiscreteTransition<R>::destination() const
{
  return *this->_source;
}


template<class R>
const System::MapInterface<R>&
System::ConstraintDiscreteTransition<R>::reset() const
{
  return *this->_reset;
}


template<class R>
const System::ConstraintInterface<R>&
System::ConstraintDiscreteTransition<R>::activation() const
{
  return *this->_activation;
}


template<class R>
bool
System::ConstraintDiscreteTransition<R>::forced() const
{
  return this->_forced;
}



template<class R>
std::ostream&
System::ConstraintDiscreteTransition<R>::write(std::ostream& os) const
{
  return os << "DiscreteTransition" 
            << "( id=" << this->id() 
            << ", source_id=" << this->source().id() 
            << ", destination_id=" << this->destination().id() 
            << ", reset=" << this->reset()
            << ", activation=" << this->activation()
            << ", forced=" << this->forced() << " )";
  
}





template<class R>
System::ConstraintHybridAutomaton<R>::ConstraintHybridAutomaton(const std::string& name)
  : _name(name) 
{ 
}


template<class R>
System::ConstraintHybridAutomaton<R>::~ConstraintHybridAutomaton() {
  this->_transitions.clear();
}


template<class R>
const System::ConstraintDiscreteMode<R>& 
System::ConstraintHybridAutomaton<R>::new_mode(id_type id,
                                               const VectorFieldInterface<R>& dynamic) 
{
  if(this->has_mode(id)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  std::vector< boost::shared_ptr< const ConstraintInterface<R> > > invariant;
  this->_modes.insert(mode_type(id,dynamic,invariant));
  return this->mode(id);
}


template<class R>
const typename System::ConstraintHybridAutomaton<R>::mode_type& 
System::ConstraintHybridAutomaton<R>::new_mode(id_type id,
                                               const VectorFieldInterface<R>& dynamic,
                                               const ConstraintInterface<R>& constraint) 
{
  if(this->has_mode(id)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  std::vector< boost::shared_ptr< const ConstraintInterface<R> > > invariant;
  invariant.push_back(boost::shared_ptr< const ConstraintInterface<R> >(constraint.clone()));
  this->_modes.insert(mode_type(id,dynamic,invariant));
  return this->mode(id);
}


template<class R>
const typename System::ConstraintHybridAutomaton<R>::transition_type& 
System::ConstraintHybridAutomaton<R>::new_transition(id_type event_id,
                                                      id_type source_id, 
                                                      id_type destination_id,
                                                      const MapInterface<R>& reset,
                                                      const ConstraintInterface<R>& activation) 
{
  if(this->has_transition(event_id,source_id)) {
    throw std::runtime_error("The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    throw std::runtime_error("The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(destination_id)) {
    throw std::runtime_error("The desitination mode of the transition must be in the automaton");
  }
  
  const mode_type& source=this->mode(source_id);
  const mode_type& destination=this->mode(destination_id);
  
  this->_transitions.insert(transition_type(event_id,source,destination,reset,activation,false));
  this->_mode(source_id)._activations.insert(std::make_pair(event_id,activation.clone()));
  return this->transition(event_id,source_id);
}



template<class R>
const typename System::ConstraintHybridAutomaton<R>::transition_type& 
System::ConstraintHybridAutomaton<R>::new_forced_transition(id_type event_id,
                                                            id_type source_id, 
                                                            id_type destination_id,
                                                            const MapInterface<R>& reset,
                                                            const ConstraintInterface<R>& guard) 
{
  if(this->has_transition(event_id,source_id)) {
    throw std::runtime_error("The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    throw std::runtime_error("The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(destination_id)) {
    throw std::runtime_error("The desitination mode of the transition must be in the automaton");
  }
  
  const mode_type& source=this->mode(source_id);
  const mode_type& destination=this->mode(destination_id);
  
  this->_transitions.insert(transition_type(event_id,source,destination,reset,guard,true));
  this->_mode(source_id)._guards.insert(std::make_pair(event_id,guard.clone()));
  return this->transition(event_id,source_id);
}




template<class R>
bool
System::ConstraintHybridAutomaton<R>::has_mode(id_type id) const 
{
  // FIXME: This is a hack since we use std::set which cannot be searched by id.
  for(discrete_mode_iterator mode_iter=this->_modes.begin();
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
System::ConstraintHybridAutomaton<R>::has_transition(id_type event_id, id_type source_id) const 
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
Geometry::HybridSpace
System::ConstraintHybridAutomaton<R>::locations() const 
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
const std::set<typename System::ConstraintHybridAutomaton<R>::mode_type>& 
System::ConstraintHybridAutomaton<R>::modes() const 
{
  return this->_modes;
}


template<class R>
const typename System::ConstraintHybridAutomaton<R>::mode_type& 
System::ConstraintHybridAutomaton<R>::mode(id_type id) const 
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
typename System::ConstraintHybridAutomaton<R>::mode_type& 
System::ConstraintHybridAutomaton<R>::_mode(id_type id) 
{
  // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given id.
  for(discrete_mode_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->id()==id) {
      return const_cast<mode_type&>(*mode_iter);
    }
  }
  throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


template<class R>
const std::set<typename System::ConstraintHybridAutomaton<R>::transition_type>& 
System::ConstraintHybridAutomaton<R>::transitions() const 
{
  return this->_transitions;
}


template<class R>
const typename System::ConstraintHybridAutomaton<R>::transition_type&
System::ConstraintHybridAutomaton<R>::transition(id_type event_id, id_type source_id) const 
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
const std::vector<typename System::ConstraintHybridAutomaton<R>::constraint_pointer>&
System::ConstraintHybridAutomaton<R>::invariants(id_type mode_id) const 
{
  if(!this->has_mode(mode_id)) {
    throw std::runtime_error("The automaton does not contain a mode with ths given id");
  }
  return this->mode(mode_id)._invariant;
}


template<class R>
const std::map<id_type, typename System::ConstraintHybridAutomaton<R>::constraint_pointer>&
System::ConstraintHybridAutomaton<R>::activations(id_type source_id) const 
{
  return this->mode(source_id)._activations;
}


template<class R>
const std::map<id_type, typename System::ConstraintHybridAutomaton<R>::constraint_pointer>&
System::ConstraintHybridAutomaton<R>::guards(id_type source_id) const 
{
  return this->mode(source_id)._guards;
}


template<class R>
typename System::ConstraintHybridAutomaton<R>::constraint_reference
System::ConstraintHybridAutomaton<R>::activation(id_type event_id, id_type source_id) const 
{
  return *this->_activations.find(source_id)->second.find(event_id)->second;
}


template<class R>
typename System::ConstraintHybridAutomaton<R>::constraint_reference
System::ConstraintHybridAutomaton<R>::guard(id_type event_id, id_type source_id) const 
{
  return *this->_guards.find(source_id)->second.find(event_id)->second;
}





template<class R>
const std::string&
System::ConstraintHybridAutomaton<R>::name() const{ 
  return this->_name; 
}


template<class R>
std::ostream& 
System::ConstraintHybridAutomaton<R>::write(std::ostream& os) const 
{
   return os << "System::ConstraintHybridAutomaton( modes=" << this->_modes << ", transitions=" << this->_transitions << ")"; 
}





}

