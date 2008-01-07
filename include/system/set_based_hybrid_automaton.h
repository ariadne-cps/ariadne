/***************************************************************************
 *            set_based_hybrid_automaton.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
#ifndef ARIADNE_SET_BASED_HYBRID_AUTOMATON_H
#define ARIADNE_SET_BASED_HYBRID_AUTOMATON_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "geometry/declarations.h"
#include "system/declarations.h"

#include "base/reference_container.h"
#include "base/stlio.h"
#include "function/function_interface.h"
#include "geometry/exceptions.h"
#include "geometry/discrete_state.h"
#include "geometry/set_interface.h"
#include "system/exceptions.h"
#include "system/map.h"
#include "system/vector_field.h"

#include "discrete_event.h"

namespace Ariadne {  

namespace Geometry {
template<class S> class HybridSet;
}

namespace System {

template<class R> class SetBasedDiscreteMode;  
template<class R> class SetBasedDiscreteTransition;  
template<class R> class SetBasedHybridAutomaton;  

/*!
 * \brief A discrete mode of a SetBasedHybridAutomaton, comprising continuous evolution given by a VectorField
 * within and invariant Geometry::ConstraintSet. 
 *
 * A %SetBasedDiscreteMode can only be created using the new_mode() method in
 * the %SetBasedHybridAutomaton class.
 */
template<class R>
class SetBasedDiscreteMode {
  friend class SetBasedHybridAutomaton<R>;
  public:
      /*!\brief The type of denotable real number used to describe to discrete mode. */
    typedef R real_type;
      
  private:
    
    // The discrete mode's discrete state.
    Geometry::DiscreteState _state;
  
    // The discrete mode's vector field.
    boost::shared_ptr< const VectorField<R> > _dynamic;
  
    // The discrete mode's invariant.
    boost::shared_ptr< const Geometry::ConstraintSet<R> > _invariant;
  
  public:
    /*! \brief Copy constructor. */
    SetBasedDiscreteMode(const SetBasedDiscreteMode<R>& original)
      : _state(original._state), _dynamic(original._dynamic),
        _invariant(original._invariant) {}
    
    /*! \brief Copy assignment operator. */
    SetBasedDiscreteMode<R>& operator=(const SetBasedDiscreteMode<R>& original) {
      if(this!=&original) {
        this->_state=original._state;
        this->_dynamic=original._dynamic;
        this->_invariant=original._invariant;
      }
      return *this;      
    }

    /*! \brief The discrete mode's identifier. (%Deprecated) */
    id_type id() const {
      return this->_state.id();
    }
    
    /*! \brief The mode's discrete state. */
    Geometry::DiscreteState discrete_state() const {
      return this->_state;
    }
    
    /*! \brief The discrete mode's dynamic (a vector field). */
    const VectorField<R>& dynamic() const {
      return *this->_dynamic;  
    }
    
    /*! \brief The discrete mode's invariant. */
    const Geometry::ConstraintSet<R>& invariant() const{
      return *this->_invariant;  
    }
    
    /*! \brief The dimension of the discrete mode. */
    dimension_type dimension() const { return this->_invariant->dimension(); }
    
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
    
   private:
    /* Construct discrete mode.
     *
     * \param id is the identifier of the mode.
     * \param dynamic is the mode's vector field.
     * \param invariant is the mode's invariant.
     */
    SetBasedDiscreteMode(Geometry::DiscreteState state,
                 const VectorField<R>& dynamic, 
                 const Geometry::ConstraintSet<R>& invariant)
      :  _state(state), _dynamic(dynamic.clone()), _invariant(invariant.clone()) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(dynamic,invariant,"SetBasedDiscreteMode::SetBasedDiscreteMode(...)");
    }
    
    /* Construct from objects managed by shared pointers (for internal use) */
    SetBasedDiscreteMode(Geometry::DiscreteState state,
                 const boost::shared_ptr< VectorField<R> > dynamic, 
                 const boost::shared_ptr< Geometry::ConstraintSet<R> > invariant)
      :  _state(state), _dynamic(dynamic), _invariant(invariant) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*dynamic,*invariant,"SetBasedDiscreteMode::SetBasedDiscreteMode(...)");
    }
    
};
  
template<class R> inline
std::ostream& 
SetBasedDiscreteMode<R>::write(std::ostream& os) const 
{ 
  return os << "SetBasedDiscreteMode( "
            << "discrete_state=" << this->discrete_state() << ", " 
            << "invariant=" << this->invariant() << ", "
            << "dynamic=" << this->dynamic() << " )"; 
}

template<class R> inline
std::ostream& operator<<(std::ostream& os, const SetBasedDiscreteMode<R>& dm) 
{
  return dm.write(os);
}

template<class R> inline
bool operator<(const SetBasedDiscreteMode<R>& mode1, const SetBasedDiscreteMode<R>& mode2) 
{
  return mode1.id() < mode2.id();
}
  



/*!
 * \brief A discrete transition of a SetBasedHybridAutomaton, representing an instantaneous
 * jump from one SetBasedDiscreteMode to another, governed by an activation Geometry::ConstraintSet and a reset MapInterface.
 *
 * A %SetBasedDiscreteTransition can only be created using the new_transition() method in
 * the %SetBasedHybridAutomaton class.
 */
template< class R >
class SetBasedDiscreteTransition
{
  friend class SetBasedHybridAutomaton<R>;
  private:
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;
  
    // \brief The source of the discrete transition.
    const SetBasedDiscreteMode<R>* _source;   
  
    // \brief The destination of the discrete transition.
    const SetBasedDiscreteMode<R>* _destination;   
  
    // \brief The activation region of the discrete transition.
    boost::shared_ptr< const Geometry::ConstraintSet<R> > _activation; 

    // \brief The reset of the discrete transition.
    boost::shared_ptr< const Map<R> > _reset;  
    
  public:
    /*! \brief Copy constructor. */
    SetBasedDiscreteTransition(const SetBasedDiscreteTransition<R>& original)
      : _event(original._event), _source(original._source), _destination(original._destination), 
        _activation(original._activation), _reset(original._reset) { 
    }
                
    /*! \brief Copy assignment operator. */
    SetBasedDiscreteTransition<R>& operator=(const SetBasedDiscreteTransition<R>& original) {
      this->_event=original._event;
      this->_source=original._source;
      this->_destination=original._destination;
      this->_activation=original._activation;
      this->_reset=original._reset;
      return *this;
    }
    
    /*! \brief The unique identifier of the discrete transition. (%Deprecated) */
    id_type id() const {
      return this->_event.id();
    }      

    /*! \brief The discrete event associated with the discrete transition. */
    DiscreteEvent discrete_event() const {
      return this->_event;
    }      

    /*! \brief The source mode of the discrete transition. */
    const SetBasedDiscreteMode<R>& source() const {
      return *this->_source;
    }      

    /*! \brief The destination of the discrete transition. */
    const SetBasedDiscreteMode<R>& destination() const { 
      return *this->_destination;
    }
  
    /*! \brief The activation region of the discrete transition. */
    const Geometry::ConstraintSet<R>& activation() const { 
      return *this->_activation;
    }

    /*! \brief The reset map of the discrete transition. */
    const Map<R>& reset() const { 
      return *this->_reset;
    }
    
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
  private:
 
    /* Constructor.
     * \param event is the discrete event.
     * \param source is the source mode of the discrete 
     * transition.
     * \param source is the source mode of the discrete 
     * transition.
     * \param destination is the destination mode of the discrete 
     * transition.
     * \param reset is the reset relation of the discrete 
     * transition.
     * \param activation is the activation region of the 
     * discrete transition.
      */
    SetBasedDiscreteTransition(DiscreteEvent event, 
                               const SetBasedDiscreteMode<R>& source, 
                               const SetBasedDiscreteMode<R>& destination,
                               const Map<R>& reset,
                               const Geometry::ConstraintSet<R>& activation)
      : _event(event), _source(&source), _destination(&destination), 
        _activation(activation.clone()), _reset(reset.clone()) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(activation,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(reset,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(reset,destination,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
    }

    /* Construct from shared pointers (for internal use). */
    SetBasedDiscreteTransition(DiscreteEvent event,
                       const SetBasedDiscreteMode<R>& source, 
                       const SetBasedDiscreteMode<R>& destination,
                       const boost::shared_ptr< Map<R> > reset,
                       const boost::shared_ptr< Geometry::ConstraintSet<R> > activation) 
      : _event(event), _source(&source), _destination(&destination), 
        _activation(activation), _reset(reset) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*activation,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(*reset,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(*reset,destination,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
    }

    /* Construct from shared pointers (for internal use). */
    SetBasedDiscreteTransition(DiscreteEvent event,
                       const boost::shared_ptr< SetBasedDiscreteMode<R> > source, 
                       const boost::shared_ptr< SetBasedDiscreteMode<R> > destination,
                       const boost::shared_ptr< Map<R> > reset,
                       const boost::shared_ptr< Geometry::ConstraintSet<R> > activation) 
      : _event(event), _source(&*source), _destination(&*destination), 
        _activation(activation), _reset(reset) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*activation,*source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(*reset,*source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(*reset,*destination,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
    }

};


template<class R> inline
std::ostream& 
SetBasedDiscreteTransition<R>::write(std::ostream& os) const 
{ 
  return os << "SetBasedDiscreteTransition( "
            << "event=" << this->discrete_event() << ", " 
            << "source=" << this->source().discrete_state() << ", "
            << "destination=" << this->destination().discrete_state() << ", "
            << "activation=" << this->activation() << ", "
            << "reset=" << this->reset() << " )"; 
}


template<class R> inline
std::ostream& operator<<(std::ostream& os, const SetBasedDiscreteTransition<R>& dt) 
{
  return dt.write(os);
}

template<class R> inline
bool operator<(const SetBasedDiscreteTransition<R>& transition1, const SetBasedDiscreteTransition<R>& transition2) 
{
  return transition1.id() < transition2.id()
    || (transition1.id() == transition2.id() 
          && transition1.source().id() < transition2.source().id());
}





/*! \ingroup HybridTime
 *  \brief A hybrid automaton, comprising continuous-time behaviour
 *  at each SetBasedDiscreteMode, coupled by instantaneous SetBasedDiscreteTransition events.
 *  The state space is given by a Geometry::HybridSet.  
 *
 * A hybrid automaton is a dynamic system with evolution in both
 * continuous time and discrete time. 
 * The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
 * where \f$q\f$ is the <em>discrete state</em> and \f$X_q\f$
 * is the <em>continuous state space</em> of corresponding to
 * each discrete state.
 *
 * For each %SetBasedDiscreteMode, the dynamics is given by a 
 * %VectorField describing the continuous dynamics,
 * and a %Set giving an invariant which must be satisified at
 * all times.
 *
 * The discrete time behaviour is specified by %SetBasedDiscreteTransition
 * objects. 
 * Each discrete transition represents an jump from a \a source
 * mode to a \a destination mode. 
 * There can be at most one discrete transition in an automaton
 * with the same event and source.
 */
template< class R >
class SetBasedHybridAutomaton
{
 public:
  typedef R real_type;
 
  typedef typename std::set< SetBasedDiscreteTransition<R> >::const_iterator discrete_transition_const_iterator;
  typedef typename std::set< SetBasedDiscreteMode<R> >::const_iterator discrete_mode_const_iterator;
 private: 
  static std::list<std::string> system_names;
  
 private:
  /*! \brief The hybrid automaton's name. */
  std::string _name;
  
  /*! \brief The hybrid automaton's identificator. */
  id_type _id;
  
  /*! \brief The list of the hybrid automaton's discrete modes. */
  std::set< SetBasedDiscreteMode<R> > _modes;
  
  /*! \brief The hybrid automaton's transitions. */
  std::set< SetBasedDiscreteTransition<R> > _transitions;
  
 public:
  //@{
  //! \name Constructors and destructors 
  
  /*! \brief Construct a named automaton. */
  SetBasedHybridAutomaton(const std::string& name);
  
  /*! \brief  Destructor. */
  ~SetBasedHybridAutomaton();
  //@}

  //@{ 
  //! \name Methods for building the automaton.

  /*! \brief Adds a discrete mode to the automaton.
   *
   * \param state is the mode's discrete state.
   * \param dynamic is the mode's vector field.
   * \param invariant is the mode's invariant.
   */
  const SetBasedDiscreteMode<R>& new_mode(Geometry::DiscreteState state,
                                          const System::VectorField<R>& dynamic,
                                          const Geometry::ConstraintSet<R>& invariant);
    
  /*! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and destination modes.
   *
   * \param event is the transition's event.
   * \param source is the transition's source location.
   * \param destination is the transition's destination location.
   * \param reset is the transition's reset.
   * \param activation is the transition's activation region.
   */
  const SetBasedDiscreteTransition<R>& new_transition(System::DiscreteEvent event,
                                                      Geometry::DiscreteState source, 
                                                      Geometry::DiscreteState destination,
                                                      const System::Map<R>& reset,
                                                      const Geometry::ConstraintSet<R>& activation);

  /*! \brief Adds a discrete transition to the automaton using the discrete modes to specify the source and destination.
   *
   * \param event is the discrete transition's discrete event. 
   * \param source is the discrete transition's source mode.
   * \param destination is the discrete transition's destination mode.
   * \param reset is the discrete transition's reset.
   * \param activation is the discrete transition's activation region.
   */
  const SetBasedDiscreteTransition<R>& new_transition(System::DiscreteEvent event,
                                                      const SetBasedDiscreteMode<R>& source, 
                                                      const SetBasedDiscreteMode<R>& destination,
                                                      const Map<R>& reset,
                                                      const Geometry::ConstraintSet<R>& activation);
  
  //@}
  
  //@{ 
  //! \name Data access and queries. 
  /*! \brief Returns the hybrid automaton's name. */
  const std::string& name() const;

  /*! \brief Test if the hybrid automaton has a discrete mode with discrete state \a state. */
  bool has_mode(Geometry::DiscreteState state) const;
  
  /*! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id. */
  bool has_transition(System::DiscreteEvent event, Geometry::DiscreteState source) const;
  
  /*! \brief The discrete mode with given discrete state. */
  const SetBasedDiscreteMode<R>& mode(Geometry::DiscreteState state) const;
  
  /*! \brief The discrete transition with given \a event and \a source location. */
  const SetBasedDiscreteTransition<R>& transition(System::DiscreteEvent event, Geometry::DiscreteState source) const;

  /*! \brief The discrete transitions from location \a source. */
  reference_vector< const SetBasedDiscreteTransition<R> > transitions(Geometry::DiscreteState source) const;

  /*! \brief A set giving the dimension of the state space for each location identifier. */
  Geometry::HybridSpace locations() const;
  
  /*! \brief The hybrid set giving the invariant for each discrete location. */
  Geometry::HybridSet<R> invariant() const;
  
  /*! \brief The set of discrete modes. (Not available in Python interface) */
  const std::set< SetBasedDiscreteMode<R> >& modes() const;
  
  /*! \brief The set of discrete transitions. (Not available in Python interface) */
  const std::set< SetBasedDiscreteTransition<R> >& transitions() const;
  
  //@}
 
  /* \brief Write to an output stream. */
  std::ostream& write(std::ostream& os) const;
};

template<class R> inline 
std::ostream& operator<<(std::ostream& os, const SetBasedHybridAutomaton<R>& ha) {
  return ha.write(os);
}

template< class R>
void dot_print(const SetBasedHybridAutomaton< R >& A);

}
}

#endif /* ARIADNE_SET_BASED_HYBRID_AUTOMATON_H */
