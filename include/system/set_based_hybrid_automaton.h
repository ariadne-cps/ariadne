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

#include "../geometry/declarations.h"
#include "../system/declarations.h"

#include "../base/stlio.h"
#include "../geometry/exceptions.h"
#include "../geometry/set_interface.h"
#include "../system/exceptions.h"
#include "../system/map_interface.h"
#include "../system/vector_field_interface.h"


namespace Ariadne {  

namespace Geometry {
template<class S> class HybridSet;
}

namespace System {

template<class R> class SetBasedDiscreteMode;  
template<class R> class SetBasedDiscreteTransition;  
template<class R> class SetBasedHybridAutomaton;  

/*!\ingroup HybridTime
 * \brief A discrete mode of a SetBasedHybridAutomaton, comprising continuous evolution given by a VectorField
 * within and invariant Geometry::SetInterface. 
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
    
    // The discrete mode's identificator.
    id_type _id;
  
    // The discrete mode's vector field.
    boost::shared_ptr< const VectorFieldInterface<R> > _dynamic;
  
    // The discrete mode's invariant.
    boost::shared_ptr< const Geometry::SetInterface<R> > _invariant;
  
  public:
    
    /*! \brief Copy constructor. */
    SetBasedDiscreteMode(const SetBasedDiscreteMode<R>& original)
      : _id(original._id), _dynamic(original._dynamic),
        _invariant(original._invariant) {}
    
    /*! \brief Copy assignment operator. */
    SetBasedDiscreteMode<R>& operator=(const SetBasedDiscreteMode<R> &original) {
      if(this!=&original) {
        this->_id=original._id;
        this->_dynamic=original._dynamic;
        this->_invariant=original._invariant;
      }
      return *this;      
    }

    /*! \brief The discrete mode's identifier. */
    id_type id() const {
      return this->_id;
    }
    
    /*! \brief The discrete mode's dynamic (a vector field). */
    const VectorFieldInterface<R>& dynamic() const {
      return *this->_dynamic;  
    }
    
    /*! \brief The discrete mode's invariant. */
    const Geometry::SetInterface<R>& invariant() const{
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
    SetBasedDiscreteMode(id_type id,
                 const VectorFieldInterface<R> &dynamic, 
                 const Geometry::SetInterface<R> &invariant)
      :  _id(id), _dynamic(dynamic.clone()), _invariant(invariant.clone()) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(dynamic,invariant,"SetBasedDiscreteMode::SetBasedDiscreteMode(...)");
    }
    
    /* Construct from objects managed by shared pointers (for internal use) */
    SetBasedDiscreteMode(id_type id,
                 const boost::shared_ptr< VectorFieldInterface<R> > dynamic, 
                 const boost::shared_ptr< Geometry::SetInterface<R> > invariant)
      :  _id(id), _dynamic(dynamic), _invariant(invariant) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*dynamic,*invariant,"SetBasedDiscreteMode::SetBasedDiscreteMode(...)");
    }
    
};
  
template<class R> inline
std::ostream& 
SetBasedDiscreteMode<R>::write(std::ostream& os) const 
{ 
  return os << "SetBasedDiscreteMode( "
            << "id=" << this->id() << ", " 
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
  



/*! \ingroup HybridTime
 * \brief A discrete transition of a SetBasedHybridAutomaton, representing an instantaneous
 * jump from one SetBasedDiscreteMode to another, governed by an activation Geometry::SetInterface and a reset MapInterface.
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
    id_type _event_id;
  
    // \brief The source of the discrete transition.
    const SetBasedDiscreteMode<R>* _source;   
  
    // \brief The destination of the discrete transition.
    const SetBasedDiscreteMode<R>* _destination;   
  
    // \brief The activation region of the discrete transition.
    boost::shared_ptr< const Geometry::SetInterface<R> > _activation; 

    // \brief The reset of the discrete transition.
    boost::shared_ptr< const MapInterface<R> > _reset;  
    
  public:
    /*! \brief Copy constructor. */
    SetBasedDiscreteTransition(const SetBasedDiscreteTransition<R> &original)
      : _event_id(original._event_id), _source(original._source), _destination(original._destination), 
        _activation(original._activation), _reset(original._reset) { 
    }
                
    /*! \brief Copy assignment operator. */
    SetBasedDiscreteTransition<R>& operator=(const SetBasedDiscreteTransition<R> &original) {
      this->_event_id=original._event_id;
      this->_source=original._source;
      this->_destination=original._destination;
      this->_activation=original._activation;
      this->_reset=original._reset;
      return *this;
    }
    
    /*! \brief The unique identifier of the discrete transition. */
    id_type id() const {
      return this->_event_id;
    }      

    /*! \brief The source mode of the discrete transition. */
    const SetBasedDiscreteMode<R> &source() const {
      return *this->_source;
    }      

    /*! \brief The destination of the discrete transition. */
    const SetBasedDiscreteMode<R> &destination() const { 
      return *this->_destination;
    }
  
    /*! \brief The activation region of the discrete transition. */
    const Geometry::SetInterface<R> &activation() const { 
      return *this->_activation;
    }

    /*! \brief The reset map of the discrete transition. */
    const MapInterface<R> &reset() const { 
      return *this->_reset;
    }
    
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
  private:
 
    /* Constructor.
     * \param event_id is the identifier of the discrete event.
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
    SetBasedDiscreteTransition(id_type event_id, 
                       const SetBasedDiscreteMode<R> &source, 
                       const SetBasedDiscreteMode<R> &destination,
                       const MapInterface<R> &reset,
                       const Geometry::SetInterface<R> &activation)
      : _event_id(event_id), _source(&source), _destination(&destination), 
        _activation(activation.clone()), _reset(reset.clone()) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(activation,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(reset,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(reset,destination,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
    }

    /* Construct from shared pointers (for internal use). */
    SetBasedDiscreteTransition(id_type event_id,
                       const SetBasedDiscreteMode<R> &source, 
                       const SetBasedDiscreteMode<R> &destination,
                       const boost::shared_ptr< MapInterface<R> > reset,
                       const boost::shared_ptr< Geometry::SetInterface<R> > activation) 
      : _event_id(event_id), _source(&source), _destination(&destination), 
        _activation(activation), _reset(reset) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*activation,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(*reset,source,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(*reset,destination,"SetBasedDiscreteTransition::SetBasedDiscreteTransition(...)");
    }

    /* Construct from shared pointers (for internal use). */
    SetBasedDiscreteTransition(id_type event_id,
                       const boost::shared_ptr< SetBasedDiscreteMode<R> > source, 
                       const boost::shared_ptr< SetBasedDiscreteMode<R> > destination,
                       const boost::shared_ptr< MapInterface<R> > reset,
                       const boost::shared_ptr< Geometry::SetInterface<R> > activation) 
      : _event_id(event_id), _source(&*source), _destination(&*destination), 
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
            << "id=" << this->id() << ", " 
            << "source_id=" << this->source().id() << ", "
            << "destination_id=" << this->destination().id() << ", "
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
 * with the same event_id and source_id.
 */
template< class R >
class SetBasedHybridAutomaton
{
 public:
  typedef R real_type;
 
  typedef typename std::set< SetBasedDiscreteTransition<R> >::const_iterator discrete_transition_iterator;
  typedef typename std::set< SetBasedDiscreteMode<R> >::const_iterator discrete_mode_iterator;
 private: 
  static std::list<std::string> system_names;
  
 protected:
  /*! \brief The hybrid automaton's name. */
  std::string _name;
  
  /*! \brief The hybrid automaton's identificator. */
  id_type _id;
  
  /*! \brief The list of the hybrid automaton's discrete modes. */
  std::set< SetBasedDiscreteMode<R> > _modes;
  
  /*! \brief The hybrid automaton's transitions. */
  std::set< SetBasedDiscreteTransition<R> > _transitions;
  
 public:
  
  /*! \brief This is a hybrid automaton class constructor.
   *  
   * This constructor initializes the object of the
   * hybrid automaton class.
   * \param name is the name of the hybrid automaton.
   */
  SetBasedHybridAutomaton(const std::string &name);
  
  /*! \brief  This is the destructor of the class hybrid 
   * automaton.
   *
   * This destructor deletes in a safe way an object of the
   * hybrid automaton class.
   */
  ~SetBasedHybridAutomaton();
  
  /*! \brief Adds a discrete mode.
   *
   * This method adds a discrete mode to automaton definition.
   * \param id is the unique key or identifyer of the discrete mode.
   * \param dynamic is the discrete mode's vector field.
   * \param invariant is the discrete mode's invariant.
   */
  const SetBasedDiscreteMode<R>& new_mode(id_type id,
                                  const VectorFieldInterface<R>& dynamic,
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
  const SetBasedDiscreteTransition<R>& new_transition(id_type event_id,
                                              const SetBasedDiscreteMode<R> &source, 
                                              const SetBasedDiscreteMode<R> &destination,
                                              const MapInterface<R> &reset,
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
  const SetBasedDiscreteTransition<R>& new_transition(Base::id_type event_id,
                                              id_type source_id, 
                                              id_type destination_id,
                                              const System::MapInterface<R> &reset,
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
  const std::set< SetBasedDiscreteMode<R> >& modes() const;
  
  /*! \brief The discrete mode with given id. */
  const SetBasedDiscreteMode<R>& mode(id_type id) const;
  
  /*! \brief The set of discrete transitions. */
  const std::set< SetBasedDiscreteTransition<R> >& transitions() const;
  
  /*! \brief The discrete transition with given \a event_id and \a source id. */
  const SetBasedDiscreteTransition<R>& transition(id_type event_id, id_type source_id) const;
  
  /*! \brief Returns the hybrid automaton's name. */
  const std::string &name() const;
  
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
