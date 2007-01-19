/***************************************************************************
 *            hybrid_automaton.h
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
 
#ifndef _ARIADNE_HYBRID_AUTOMATON_H
#define _ARIADNE_HYBRID_AUTOMATON_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../base/stlio.h"
#include "../system/discrete_mode.h"
#include "../system/discrete_transition.h"



namespace Ariadne {  
namespace System {
  
template< class R > class HybridAutomaton;
template< class R> void dot_print(const HybridAutomaton< R >& A);

  
/*! \ingroup HybridTime
 *  \brief The class of hybrid automata, (dynamic systems without
 *  input, which have both discrete and continuous dynamics).
 * 
 * It represents automata and its methods evaluate both reachability
 * properties and region.
 */
template< class R >
class HybridAutomaton
{
 public:
  typedef R real_type;
 
  typedef typename std::vector< DiscreteTransition<R> >::const_iterator discrete_transition_iterator;
  typedef typename std::vector< DiscreteMode<R> >::const_iterator discrete_mode_iterator;
 private: 
  static std::list<std::string> system_names;
  
 protected:
  /*! \brief The hybrid automaton's name. */
  std::string _name;
  
  /*! \brief The hybrid automaton's identificator. */
  id_type _id;
  
  /*! \brief The list of the hybrid automaton's discrete modes. */
  std::vector< DiscreteMode<R> > _modes;
  
  /*! \brief The hybrid automaton's transitions. */
  std::vector< DiscreteTransition<R> > _transitions;
  
 public:
  
  /*! \brief This is a hybrid automaton class constructor.
   *  
   * This constructor initializes the object of the
   * hybrid automaton class.
   * \param name is the name of the hybrid automaton.
   */
  HybridAutomaton(const std::string &name): _name(name) { }
  
  /*! \brief  This is the destructor of the class hybrid 
   * automaton.
   *
   * This destructor deletes in a safe way an object of the
   * hybrid automaton class.
   */
  ~HybridAutomaton() {
    //for (size_t i=0; i< (this->_transitions).size(); i++) {
    //  this->_transitions[i].clear();
    //}
    this->_transitions.clear();
  }
  
  /*! \brief Adds a discrete mode.
   *
   * This method adds a discrete mode to automaton definition.
   * \param dynamic is the discrete mode's vector field.
   * \param invariant is the discrete mode's invariant.
   */
  inline DiscreteMode<R>& new_mode(const VectorField<R>& dynamic,
                                   const Geometry::Set<R>& invariant) 
  {
    this->_modes.push_back(DiscreteMode<R>(this->_modes.size(),dynamic,invariant));
    return this->_modes.back();
  }
    
  /*! \brief Adds a discrete transition.
   *
   * This method adds an arc and the relavite reset between 
   * two discrete modes.
   * \param reset is the discrete transition's reset.
   * \param activation is the discrete transition's activation region.
   * \param source is the discrete transition's source.
   * \param destination is the discrete transition's destination.
   */
  inline DiscreteTransition<R>& new_transition(const Map<R> &reset,
                                               const Geometry::Set<R> &activation,
                                               const DiscreteMode<R> &source, 
                                               const DiscreteMode<R> &destination) 
  {
    if(&this->_modes[source.id()] != &source) {
      throw std::runtime_error("The source mode of the transition must be in the automaton");
    }
    if(&this->_modes[source.id()] != &destination) {
      throw std::runtime_error("The desitination mode of the transition must be in the automaton");
    }
    
    this->_transitions.push_back(DiscreteTransition<R>(this->_transitions.size(),reset,activation,source,destination));
    return this->_transitions.back();
  }
  
  /*! \brief Adds a discrete transition.
   *
   * This method adds an arc and the relavite reset between 
   * two discrete modes.
   * \param reset is the discrete transition's reset.
   * \param activation is the discrete transition's activation region.
   * \param source_id is the identifyer of the discrete transition's source mode.
   * \param destination_id is the identifyer of the discrete transition's destination mode.
   */
  inline DiscreteTransition<R>& new_transition(const Map<R> &reset,
                                               const Geometry::Set<R> &activation,
                                               const id_type &source_id, 
                                               const id_type &destination_id) 
  {
    if(source_id >= this->modes().size()) {
      throw std::runtime_error("The automaton does not contain a mode with ths given id");
    }
    if(destination_id >= this->modes().size()) {
      throw std::runtime_error("The automaton does not contain a mode with ths given id");
    }
    
    DiscreteMode<R>& source=this->_modes[source_id];
    DiscreteMode<R>& destination=this->_modes[destination_id];
    this->_transitions.push_back(DiscreteTransition<R>(this->_transitions.size(),reset,activation,source,destination));
    return this->_transitions.back();
  }
  
  /*! \brief The list of discrete modes. */
  inline const std::vector< DiscreteMode<R> >& modes() const 
  {
    return this->_modes;
  }
  
  /*! \brief The discrete mode with given id. */
  inline const DiscreteMode<R>& mode(id_type id) const 
  {
    return this->_modes[id];
  }
  
  /*! \brief The list of discrete transitions. */
  inline const std::vector< DiscreteTransition<R> >& transitions() const 
  {
    return this->_transitions;
  }
  
  /*! \brief The discrete transition with given id. */
  inline const DiscreteTransition<R>& transition(id_type id) const 
  {
    return this->_transitions[id];
  }
  
  /*! \brief Returns the hybrid automaton's name. */
  inline const std::string &name() const{ 
    return this->_name; 
  }
  
  inline std::ostream& write(std::ostream& os) const {
    return os << "HybridAutomaton( modes=" << this->_modes << ", transitions=" << this->_transitions << ")"; 
  }
};

template<class R> inline 
std::ostream& operator<<(std::ostream& os, const HybridAutomaton<R>& ha) {
  return ha.write(os);
}

/*

template< class R>
inline void dot_print(const HybridAutomaton< R >& A) 
{          
  std::ofstream fos;
  
  std::string f_name=A.name();
  
  f_name+=".dot";
  
  fos.open(f_name.c_str() , std::ios::out);
  
  size_t arc_number=0;
  
  fos << "digraph \""<< A.name()<<"\" {" << std::endl
      << " rankdir=LR; "<< std::endl
      << " node [shape = circle]; "<< std::endl;
  
  for (size_t i=0; i<(A._modes).size(); i++) {
    std::string l_name=A._modes[i].name();
    for (size_t j=0; j<(A._automaton[i]).size(); j++) {
      fos << "\"" <<  l_name << "\" -> \"" 
          << (((A._automaton[i])[j]).destination()).name() 
          << "\" [ label=\"a_" << arc_number++ << "\" ]; " << std::endl;
    }      
  }
  
  fos << "}" << std::endl;
  
  fos.close();
}
*/

}
}

#endif /* _ARIADNE_HYBRID_AUTOMATON_H */
