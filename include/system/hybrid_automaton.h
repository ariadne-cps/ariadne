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
 
#ifndef ARIADNE_HYBRID_AUTOMATON_H
#define ARIADNE_HYBRID_AUTOMATON_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../base/stlio.h"
#include "../geometry/hybrid_space.h"
#include "../system/discrete_mode.h"
#include "../system/discrete_transition.h"



namespace Ariadne {  

namespace Geometry {
template<class S> class HybridSet;
}

namespace System {
  
  
/*! \ingroup HybridTime
 *  \brief A hybrid automaton, comprising continuous-time behaviour
 *  at each DiscreteMode, coupled by instantaneous DiscreteTransition events.
 *  The state space is given by a Geometry::HybridSet.  
 *
 * A hybrid automaton is a dynamic system with evolution in both
 * continuous time and discrete time. 
 * The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
 * where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
 * is the <em>continuous state space</em> of corresponding to
 * each discrete state.
 *
 * For each %DiscreteMode, the dynamics is given by a 
 * %VectorField describing the continuous dynamics,
 * and a %Set giving an invariant which must be satisified at
 * all times.
 *
 * The discrete time behaviour is specified by %DiscreteTransition
 * objects. 
 * Each discrete transition represents an jump from a \a source
 * mode to a \a destination mode. 
 * There can be at most one discrete transition in an automaton
 * with the same event_id and source_id.
 */
template< class R >
class HybridAutomaton
{
 public:
  typedef R real_type;
 
  typedef typename std::set< DiscreteTransition<R> >::const_iterator discrete_transition_iterator;
  typedef typename std::set< DiscreteMode<R> >::const_iterator discrete_mode_iterator;
 private: 
  static std::list<std::string> system_names;
  
 protected:
  /*! \brief The hybrid automaton's name. */
  std::string _name;
  
  /*! \brief The hybrid automaton's identificator. */
  id_type _id;
  
  /*! \brief The list of the hybrid automaton's discrete modes. */
  std::set< DiscreteMode<R> > _modes;
  
  /*! \brief The hybrid automaton's transitions. */
  std::set< DiscreteTransition<R> > _transitions;
  
 public:
  
  /*! \brief This is a hybrid automaton class constructor.
   *  
   * This constructor initializes the object of the
   * hybrid automaton class.
   * \param name is the name of the hybrid automaton.
   */
  HybridAutomaton(const std::string &name);
  
  /*! \brief  This is the destructor of the class hybrid 
   * automaton.
   *
   * This destructor deletes in a safe way an object of the
   * hybrid automaton class.
   */
  ~HybridAutomaton();
  
  /*! \brief Adds a discrete mode.
   *
   * This method adds a discrete mode to automaton definition.
   * \param id is the unique key or identifyer of the discrete mode.
   * \param dynamic is the discrete mode's vector field.
   * \param invariant is the discrete mode's invariant.
   */
  const DiscreteMode<R>& new_mode(id_type id,
                                  const VectorField<R>& dynamic,
                                  const Geometry::SetInterface<R>& invariant);
    
  /*! \brief Add a discrete transition.
   *
   * This method creates a new discrete transition from the source mode to the 
   * destination mode.
   * \param event_id is the unique identifyer of the discrete event. 
   * \param source is the discrete transition's source.
   * \param destination is the discrete transition's destination.
   * \param reset is the discrete transition's reset.
   * \param activation is the discrete transition's activation region.
   */
  const DiscreteTransition<R>& new_transition(id_type event_id,
                                              const DiscreteMode<R> &source, 
                                              const DiscreteMode<R> &destination,
                                              const Map<R> &reset,
                                              const Geometry::SetInterface<R> &activation);
  
  /*! \brief Adds a discrete transition.
   *
   * This method creates a new discrete transition from the mode with  mode to the 
   * destination mode.
   * \param event_id is the identifier of the discrete transition's event.
   * \param source_id is the identifier of the discrete transition's source mode.
   * \param destination_id is the identifier of the discrete transition's destination mode.
   * \param reset is the discrete transition's reset.
   * \param activation is the discrete transition's activation region.
   */
  const DiscreteTransition<R>& new_transition(id_type event_id,
                                              id_type source_id, 
                                              id_type destination_id,
                                              const Map<R> &reset,
                                              const Geometry::SetInterface<R> &activation);
  
  /*! \brief Test if the hybrid automaton has a discrete mode with key id. */
  bool has_mode(id_type id) const;
  
  /*! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id. */
  bool has_transition(id_type event_id, id_type source_id) const;
  
  /*! \brief A set giving the dimension of the state space for each location identifier. */
  Geometry::HybridSpace locations() const;
  
  /*! \brief The hybrid set giving the invariant for each discrete location. */
  Geometry::HybridSet<R> invariant() const;
  
  /*! \brief The set of discrete modes. */
  const std::set< DiscreteMode<R> >& modes() const;
  
  /*! \brief The discrete mode with given id. */
  const DiscreteMode<R>& mode(id_type id) const;
  
  /*! \brief The set of discrete transitions. */
  const std::set< DiscreteTransition<R> >& transitions() const;
  
  /*! \brief The discrete transition with given \a event_id and \a source id. */
  const DiscreteTransition<R>& transition(id_type event_id, id_type source_id) const;
  
  /*! \brief Returns the hybrid automaton's name. */
  const std::string &name() const;
  
  std::ostream& write(std::ostream& os) const;
};

template<class R> inline 
std::ostream& operator<<(std::ostream& os, const HybridAutomaton<R>& ha) {
  return ha.write(os);
}

template< class R>
void dot_print(const HybridAutomaton< R >& A);

}
}

#endif /* ARIADNE_HYBRID_AUTOMATON_H */
