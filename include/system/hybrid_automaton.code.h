/***************************************************************************
 *            hybrid_automaton.code.h
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
#include "../geometry/set_interface.h"
#include "../geometry/set_reference.h"
#include "../geometry/set_constraint.h"
#include "../geometry/hybrid_set.h"
#include "../system/map_interface.h"
#include "../system/vector_field_interface.h"
#include "../output/logging.h"

#include "hybrid_automaton.h"


namespace Ariadne {  
  

template<class R>
System::DiscreteMode<R>::DiscreteMode(id_type id, 
                                      const VectorFieldInterface<R>& dynamic)
  : _id(id), _dynamic(dynamic.clone()), _invariants()
{
}


template<class R>
id_type
System::DiscreteMode<R>::id() const
{
  return this->_id;
}


template<class R>
dimension_type
System::DiscreteMode<R>::dimension() const
{
  return this->_dynamic->dimension();
}


template<class R>
const System::VectorFieldInterface<R>&
System::DiscreteMode<R>::dynamic() const
{
  return *this->_dynamic;
}


/*
template<class R>
const System::ConstraintInterface<R>&
System::DiscreteMode<R>::invariant(id_type event_id) const
{
  typename std::map< id_type, constraint_const_pointer >::const_iterator inv_iter=this->_invariants.find(event_id);
  if(inv_iter==this->_invariants.end()) {
    throw std::runtime_error("The discrete mode has no invariant with the given event id");
  }
  return *inv_iter->second;
}


template<class R>
const System::ConstraintInterface<R>&
System::DiscreteMode<R>::activation(id_type event_id) const
{
  typename std::map< id_type, constraint_const_pointer >::const_iterator act_iter=this->_activations.find(event_id);
  if(act_iter==this->_activations.end()) {
    throw std::runtime_error("The discrete mode has no unforced transition with the given event id");
  }
  return *act_iter->second;
}

*/

template<class R>
const reference_set< const System::DiscreteTransition<R> >&
System::DiscreteMode<R>::transitions() const
{
  return this->_transitions;
}


/*
template<class R>
const reference_map< id_type, const Geometry::ConstraintInterface<R> >&
System::DiscreteMode<R>::constraints() const
{
  return this->_constraints;
}
*/


 /*
template<class R>
const System::ConstraintInterface<R>&
System::DiscreteMode<R>::guard(id_type event_id) const
{
  typename std::map< id_type, constraint_const_pointer >::const_iterator grd_iter=this->_guards.find(event_id);
  if(grd_iter==this->_guards.end()) {
    throw std::runtime_error("The discrete mode has no forced transition with the given event id");
  }
  return *grd_iter->second;
}
 */

template<class R>
std::ostream&
System::DiscreteMode<R>::write(std::ostream& os) const
{
  typedef typename std::set< boost::shared_ptr< const ConstraintInterface<R> > >::const_iterator constraint_iterator;
  os << "DiscreteMode( id=" << this->id() 
            << ", dynamic=" << this->dynamic() << ", invariants=\{";
  for(typename std::map<id_type,constraint_const_pointer>::const_iterator iter=this->_invariants.begin();
      iter!=this->_invariants.end(); ++iter) 
  {
    os << *iter->second << ","; 
  }
  os << "}";
  return os << " )";
}





template<class R>
System::DiscreteTransition<R>::DiscreteTransition(id_type id, 
                                                  const DiscreteMode<R>& source,
                                                  const ConstraintInterface<R>& invariant)
  : _event_id(id), _source(&source), _destination(), _reset(), _constraint(invariant.clone()), _event_kind(invariant_tag)
{
}


template<class R>
System::DiscreteTransition<R>::DiscreteTransition(id_type id, 
                                                  const DiscreteMode<R>& source, 
                                                  const DiscreteMode<R>& destination,
                                                  const MapInterface<R>& reset, 
                                                  const ConstraintInterface<R>& activation,
                                                  bool forced)
  : _event_id(id), _source(&source), _destination(&destination), 
    _reset(reset.clone()), _constraint(activation.clone()), 
    _event_kind(forced?guard_tag:activation_tag)
{
}


template<class R>
id_type
System::DiscreteTransition<R>::id() const
{
  return this->_event_id;
}


template<class R>
id_type
System::DiscreteTransition<R>::source_id() const
{
  return this->_source->id();
}


template<class R>
const System::DiscreteMode<R>&
System::DiscreteTransition<R>::source() const
{
  return *this->_source;
}


template<class R>
id_type
System::DiscreteTransition<R>::destination_id() const
{
  return this->_destination->id();
}


template<class R>
const System::DiscreteMode<R>&
System::DiscreteTransition<R>::destination() const
{
  return *this->_destination;
}


template<class R>
const System::MapInterface<R>&
System::DiscreteTransition<R>::reset() const
{
  return *this->_reset;
}


template<class R>
const System::ConstraintInterface<R>&
System::DiscreteTransition<R>::constraint() const
{
  return *this->_constraint;
}


template<class R>
System::EventKind
System::DiscreteTransition<R>::kind() const
{
  return this->_event_kind;
}


template<class R>
bool
System::DiscreteTransition<R>::forced() const
{
  return this->_event_kind!=activation_tag;
}



template<class R>
std::ostream&
System::DiscreteTransition<R>::write(std::ostream& os) const
{
  os << "DiscreteTransition" 
     << "( id=" << this->id() 
     << ", source_id=" << this->source().id();
  if(this->kind()!=invariant_tag) {
    os << ", destination_id=" << this->destination().id() 
       << ", reset=" << this->reset();
  }
  os << ", constraint=" << this->constraint()
     << ", forced=" << this->forced() << " )";
  return os;
}





template<class R>
System::HybridAutomaton<R>::HybridAutomaton(const std::string& name)
  : _name(name) 
{ 
}


template<class R>
System::HybridAutomaton<R>::~HybridAutomaton() {
  this->_transitions.clear();
}


template<class R>
const System::VectorFieldInterface<R>&
System::HybridAutomaton<R>::new_dynamic(id_type id,
                                        const System::VectorFieldInterface<R>& dynamic) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
const System::MapInterface<R>&
System::HybridAutomaton<R>::new_reset(id_type id,
                                      const System::MapInterface<R>& reset) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
const Geometry::ConstraintInterface<R>&
System::HybridAutomaton<R>::new_constraint(id_type id,
                                           const Geometry::ConstraintInterface<R>& constraint) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
const Geometry::SetInterface<R>&
System::HybridAutomaton<R>::new_domain(id_type id,
                                        const Geometry::SetInterface<R>& domain) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
const System::DiscreteMode<R>& 
System::HybridAutomaton<R>::new_mode(id_type mode_id, id_type dynamic_id)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
const System::DiscreteTransition<R>& 
System::HybridAutomaton<R>::new_invariant(id_type event_id, id_type mode_id, id_type constraint_id)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
const System::DiscreteTransition<R>& 
System::HybridAutomaton<R>::new_transition(id_type event_id, id_type source_id, id_type destination_id, id_type reset_id, id_type constraint_id, bool forced)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
const System::DiscreteTransition<R>& 
System::HybridAutomaton<R>::new_forced_transition(id_type event_id, id_type source_id, id_type destination_id, id_type reset_id, id_type constraint_id)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
const System::DiscreteTransition<R>& 
System::HybridAutomaton<R>::new_unforced_transition(id_type event_id, id_type source_id, id_type destination_id, id_type reset_id, id_type constraint_id)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




template<class R>
const System::DiscreteMode<R>& 
System::HybridAutomaton<R>::new_mode(id_type id,
                                               const VectorFieldInterface<R>& dynamic) 
{
  if(this->has_mode(id)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  std::map< id_type,  boost::shared_ptr< const ConstraintInterface<R> > > invariant;
  this->_modes.insert(mode_type(id,dynamic));
  return this->mode(id);
}



template<class R>
const typename System::HybridAutomaton<R>::transition_type& 
System::HybridAutomaton<R>::new_invariant(id_type event_id,
                                                    id_type mode_id,
                                                    const ConstraintInterface<R>& constraint) 
{
  ARIADNE_LOG(4,"  new_invariant(event_id="<<event_id<<", mode_id="<<mode_id<<", constraint="<<constraint<<"\n");

  if(this->has_transition(event_id,mode_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant","with event_id="<<event_id<<", mode_id="<<mode_id<<": The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(mode_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant","with event_id="<<event_id<<", mode_id="<<mode_id<<": The source mode of the transition must be in the automaton");
  }
  
  mode_type& mode=this->_mode(mode_id);
  this->_transitions.insert(transition_type(event_id,mode,constraint));
  const transition_type& transition=this->transition(event_id,mode_id);;
  mode._transitions.insert(transition);
  this->_mode(mode_id)._invariants.insert(std::make_pair(event_id,constraint.clone()));
  this->_mode(mode_id)._invariants.insert(std::make_pair(event_id,constraint.clone()));
  this->_mode(mode_id)._constraints.insert(event_id,transition.constraint());
  return transition;
}



template<class R>
const typename System::HybridAutomaton<R>::transition_type& 
System::HybridAutomaton<R>::new_transition(id_type event_id,
                                                     id_type source_id, 
                                                     id_type destination_id,
                                                     const MapInterface<R>& reset,
                                                     const ConstraintInterface<R>& constraint,
                                                     bool forced) 
{
  if(forced) {
    return this->new_forced_transition(event_id,source_id,destination_id,reset,constraint);
  } else {
    return this->new_unforced_transition(event_id,source_id,destination_id,reset,constraint);
  }
}


template<class R>
const typename System::HybridAutomaton<R>::transition_type& 
System::HybridAutomaton<R>::new_unforced_transition(id_type event_id,
                                                              id_type source_id, 
                                                              id_type destination_id,
                                                              const MapInterface<R>& reset,
                                                              const ConstraintInterface<R>& activation) 
{
  ARIADNE_LOG(4,"  new_unforced_transition(event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<", reset, guard)\n");

  if(this->has_transition(event_id,source_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_unforced_transition","with event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<": The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_unforced_transition","with event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<": The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(destination_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_unforced_transition","with event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<": The destination mode of the transition must be in the automaton");
  }
  
  mode_type& source=this->_mode(source_id);
  mode_type& destination=this->_mode(destination_id);
  
  this->_transitions.insert(transition_type(event_id,source,destination,reset,activation,false));
  const transition_type& transition=this->transition(event_id,source_id);;
  source._transitions.insert(transition);
  this->_mode(source_id)._activations.insert(std::make_pair(event_id,activation.clone()));
  this->_mode(source_id)._constraints.insert(event_id,transition.constraint());
  return transition;
}



template<class R>
const typename System::HybridAutomaton<R>::transition_type& 
System::HybridAutomaton<R>::new_forced_transition(id_type event_id,
                                                            id_type source_id, 
                                                            id_type destination_id,
                                                            const MapInterface<R>& reset,
                                                            const ConstraintInterface<R>& guard) 
{
  ARIADNE_LOG(4,"  new_forced_transition(event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<", reset, guard)\n");

  if(this->has_transition(event_id,source_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_forced_transition(...)","with event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<": The automaton already has a transition with the given event_id and source id.");
  }
  if(!this->has_mode(source_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_forced_transition","with event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<": The source mode of the transition must be in the automaton");
  }
  if(!this->has_mode(destination_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_forced_transition","with event_id="<<event_id<<", source_id="<<source_id<<", destination_id="<<destination_id<<": The destination mode of the transition must be in the automaton");
  }
  
  mode_type& source=this->_mode(source_id);
  mode_type& destination=this->_mode(destination_id);
  
  this->_transitions.insert(transition_type(event_id,source,destination,reset,guard,true));
  const transition_type& transition=this->transition(event_id,source_id);;
  source._transitions.insert(transition);
  this->_mode(source_id)._guards.insert(std::make_pair(event_id,guard.clone()));
  this->_mode(source_id)._constraints.insert(event_id,transition.constraint());
  ARIADNE_LOG(5,"  transition="<<transition<<"\n");
  return transition;
}


template<class R>
const System::DiscreteMode<R>& 
System::HybridAutomaton<R>::new_mode(id_type id,
                                     const System::VectorFieldInterface<R>& dynamic,
                                     const Geometry::SetInterface<R>& invariant) 
{
  if(this->has_mode(id)) {
    throw std::runtime_error("The hybrid automaton already has a mode with the given id");
  }
  if(this->has_transition(id,id)) {
    throw std::runtime_error("The hybrid automaton already has an event with the given id; cannot make new invariant.");
  }
  Geometry::SetConstraint<R> invariant_constraint(invariant,true);
  this->new_mode(id,dynamic);
  this->new_invariant(id,id,invariant_constraint);
  return this->mode(id);
}


template<class R>
const System::DiscreteTransition<R>& 
System::HybridAutomaton<R>::new_transition(id_type event_id,
                                           id_type source_id, 
                                           id_type destination_id,
                                           const System::MapInterface<R>& reset,
                                           const Geometry::SetInterface<R>& activation)
{
  Geometry::SetConstraint<R> activation_constraint(activation,false);
  this->new_unforced_transition(event_id,source_id,destination_id,reset,activation_constraint);
  return this->transition(event_id,source_id);
}


template<class R>
bool
System::HybridAutomaton<R>::has_mode(id_type id) const 
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
System::HybridAutomaton<R>::has_transition(id_type event_id, id_type source_id) const 
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
System::HybridAutomaton<R>::locations() const 
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
const std::set<typename System::HybridAutomaton<R>::mode_type>& 
System::HybridAutomaton<R>::modes() const 
{
  return this->_modes;
}


template<class R>
const typename System::HybridAutomaton<R>::mode_type& 
System::HybridAutomaton<R>::mode(id_type id) const 
{
  // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given id.
  for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->id()==id) {
      return *mode_iter;
    }
  }
  ARIADNE_THROW(std::runtime_error,"HybridAutomaton::mode(...)"," with id="<<id<<": The automaton does not have a mode with the given id.");
}


template<class R>
typename System::HybridAutomaton<R>::mode_type& 
System::HybridAutomaton<R>::_mode(id_type id) 
{
  ARIADNE_LOG(6,"HybridAutomaton::mode(id="<<id<<")\n");
  // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given id.
  for(discrete_mode_iterator mode_iter=this->_modes.begin();
      mode_iter!=this->_modes.end(); ++mode_iter) 
  {
    if(mode_iter->id()==id) {
      return const_cast<mode_type&>(*mode_iter);
    }
  }
  ARIADNE_THROW(std::runtime_error,"HybridAutomaton::mode(...)","with id="<<id<<": The automaton does not have a mode with the given id.");
}


template<class R>
const std::set<typename System::HybridAutomaton<R>::transition_type>& 
System::HybridAutomaton<R>::transitions() const 
{
  return this->_transitions;
}


template<class R>
const typename System::HybridAutomaton<R>::transition_type&
System::HybridAutomaton<R>::transition(id_type event_id, id_type source_id) const 
{
  ARIADNE_LOG(6,"HybridAutomaton::transition(event_id="<<event_id<<", source_id="<<source_id<<")\n");
  for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
      transition_iter!=this->_transitions.end(); ++transition_iter) 
  {
    if(transition_iter->id()==event_id && transition_iter->source().id()==source_id) {
      return *transition_iter;
    } 
  }
  ARIADNE_THROW(std::runtime_error,"HybridAutomaton::transition(...)","with event_id="<<event_id<<", source_id="<<source_id<<": The automaton does not have a transition with the given event_id and source_id.");
} 


template<class R>
const reference_set<const typename System::HybridAutomaton<R>::transition_type>&
System::HybridAutomaton<R>::transitions(id_type mode_id) const 
{
  if(!this->has_mode(mode_id)) {
    ARIADNE_THROW(std::runtime_error,"HybridAutomaton::invariants(...)","with mode_id="<<mode_id<<": The automaton does not have a mode with the given id.");
  }
  return this->mode(mode_id)._transitions;
}




template<class R>
const std::string&
System::HybridAutomaton<R>::name() const{ 
  return this->_name; 
}


template<class R>
std::ostream& 
System::HybridAutomaton<R>::write(std::ostream& os) const 
{
   return os << "System::HybridAutomaton( modes=" << this->_modes << ", transitions=" << this->_transitions << ")"; 
}





}

