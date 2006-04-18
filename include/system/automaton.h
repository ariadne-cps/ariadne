/***************************************************************************
 *            automaton.h
 *
 *  Tue Mar 23 14:12:31 2004
 *  Copyright  2004  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _AUTOMATON_H
#define _AUTOMATON_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../system/discrete_transition.h"
#include "../system/discrete_location.h"

#include "../evaluation/solver.h"


namespace Ariadne {  
namespace HybridSystem {
  
template < typename LDT >
class HybridAutomaton;
  
template < typename LDT>
inline void dot_print(const HybridAutomaton< LDT >& A){          
  std::ofstream fos;
  
  std::string f_name=A.name();
  
  f_name+=".dot";
  
  fos.open(f_name.c_str() , std::ios::out);
  
  size_t arc_number=0;
  
  fos << "digraph \""<< A.name()<<"\" {" << std::endl
      << " rankdir=LR; "<< std::endl
      << " node [shape = circle]; "<< std::endl;
  
  for (size_t i=0; i<(A._locations).size(); i++) {
    std::string l_name=A._locations[i].name();
    for (size_t j=0; j<(A._automaton[i]).size(); j++) {
      fos << "\"" <<  l_name << "\" -> \"" 
          << (((A._automaton[i])[j]).destination()).name() 
          << "\" [ label=\"a_" << arc_number++ << "\" ]; " << std::endl;
    }      
  }
  
  fos << "}" << std::endl;
  
  fos.close();
}


/*! \class HybridAutomaton
 *  \brief This is the class of automata.
 * 
 * It represents automata and its methods evaluate both reachability
 * properties and region.
 */
template < typename LDT >
class HybridAutomaton
{
 private: 
  static std::list<std::string> system_names;
  
 public:
  typedef LDT LeavingTrans;
  typedef typename LeavingTrans::DiscreteLocation DiscreteLocation;
  typedef typename LeavingTrans::ResetMap ResetMap;
  typedef typename DiscreteLocation::Vector_typeField Vector_typeField;
  typedef typename ResetMap::DenotableSet DenotableSet;
  typedef typename DenotableSet::BasicSet BasicSet;
  typedef typename BasicSet::state_type state_type;
  typedef typename state_type::real_type real_type;
  
 protected:
  /*! \brief HybridAutomaton's name. */
  std::string _name;
  
  /*! \brief HybridAutomaton's identificator. */
  HybridAutomatonID _id;
  
  /*! \brief The list of the automaton's locations. */
  std::vector<DiscreteLocation> _locations;
  
  /*! \brief The discrete automaton. */
  std::vector< std::vector<LeavingTrans> > _automaton;
  
  /*! \brief The automaton space dimension. */
  size_t _dimension;
  
 public:
  
  /*! \brief This is a hybrid automaton class constructor.
   *  
   * This constructor initializes the object of the
   * hybrid automaton class.
   * \param name is the name of the hybrid automaton.
   */
  HybridAutomaton(const std::string &name): _name(name), _dimension(0) {}
  
  /*! \brief  This is the destructor of the class hybrid 
   * automaton.
   *
   * This destructor deletes in a safe way an object of the
   * hybrid automaton class.
   */
  ~HybridAutomaton() {
    
    for (size_t i=0; i< (this->_automaton).size(); i++) {
      this->_automaton[i].clear();
    }
    
    this->_automaton.clear();
    
  }
  
  /*! \brief Adds a location.
   *
   * This method adds a location to automaton definition.
   * \param A is the new location.
   */
  inline void add_location(DiscreteLocation &A) {
    
    if ((this->dimension()!=0)&&((A.Vector_field()).dimension()!=this->dimension())) {
      throw std::invalid_argument("The automaton has a different space dimension with respect to the discrete location's continuous objects.");
    }
    
    this->_dimension=(A.Vector_field()).dimension();
    
    /* If this location is already present into the 
     * location <vector> there is nothing to do. */
    for (size_t i=0; i<(this->_locations).size(); i++) {
      if (A.name()== this->_locations[i].name()) {
        
        A._set_id(i);
        
        return;
      }
    }
    
    /* if it is a new location, insert it into the locations
     * <vector> and upgrade its id. */
    A._set_id((this->_locations).size());
    
    (this->_locations).push_back(A);
    
    /* upgrade the automaton */
    std::vector<LeavingTrans> vec_leaving_A;
    (this->_automaton).push_back(vec_leaving_A);
  }
  
  /*! \brief Adds an arc.
   *
   * This method adds an arc and the relavite reset between 
   * two locations.
   * \param source is the discrete arc's source.
   * \param dest is the discrete arc's destination.
   * \param act is the arc's activation region.
   * \param reset is the discrete arc's reset.
   */
  inline void add_arc(DiscreteLocation &source, 
                      DiscreteLocation &dest, 
                      const DenotableSet &act,
                      const ResetMap &reset) {
    
    this->add_location(source);
    this->add_location(dest);
    
    LeavingTrans arc(dest,act,reset);
    
    this->_automaton[source.id()].push_back(arc);
  }
  
  /*! \brief Adds an arc.
   *
   * This method adds an arc and the relavite reset between 
   * two locations.
   * \param source is the discrete arc's source.
   * \param dest is the discrete arc's destination.
   * \param act is the arc's activation region.
   * \param reset is the discrete arc's reset.
   */
  inline void add_arc(DiscreteLocation &source, 
                      DiscreteLocation &dest, 
                      const BasicSet &act,
                      const ResetMap &reset) {
    
    this->add_location(source);
    this->add_location(dest);
    
    DenotableSet DS_act(act);
    
    LeavingTrans arc(dest, DS_act , reset);
    
    this->_automaton[source.id()].push_back(arc);
  }
  
  /*! \brief Returns the hybrid automaton's name.
   *
   * \return The name of the hybrid automaton.
   */
  inline const std::string &name() const{ 
    return this->_name; 
  }
  
  /*! \brief Returns the hybrid automaton's space dimension. (Deprecated) */
  inline size_t dim() const { 
    return this->_dimension; 
  }
   
  /*! \brief Returns the hybrid automaton's space dimension.
   *
   * \return The automaton's space dimension.
   */
  inline size_t dimension() const { 
    return this->_dimension; 
  }
  
  template <typename AUTO , typename HDS ,typename TRACE, 
            typename MAINTAIN, typename INT>
  friend class Ariadne::Evaluation::Solver;
  
  template <typename LeavingDT>
  friend void dot_print(const HybridAutomaton< LeavingDT >& A);
};
  
}
}

#endif /* _AUTOMATON_H */
