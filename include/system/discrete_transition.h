/****************************************************************************
 *            discrete_transition.h
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

#ifndef ARIADNE_DISCRETE_TRANSITION_H
#define ARIADNE_DISCRETE_TRANSITION_H

#include <stdexcept>
#include <memory>

#include <boost/smart_ptr.hpp>

#include "../geometry/set_interface.h"
#include "../system/exceptions.h"
#include "../system/discrete_mode.h"
#include "../system/vector_field.h"
#include "../system/map.h"

namespace Ariadne {
namespace System {
    
template< class R > class HybridAutomaton;


/*! \ingroup HybridTime
 * \brief A discrete transition of a HybridAutomaton, representing an instantaneous
 * jump from one DiscreteMode to another, governed by an activation Geometry::SetInterface and a reset MapInterface.
 *
 * A %DiscreteTransition can only be created using the new_transition() method in
 * the %HybridAutomaton class.
 */
template< class R >
class DiscreteTransition
{
  friend class HybridAutomaton<R>;
  private:
    // \brief The discrete transition's identificator.
    id_type _event_id;
  
    // \brief The source of the discrete transition.
    const DiscreteMode<R>* _source;   
  
    // \brief The destination of the discrete transition.
    const DiscreteMode<R>* _destination;   
  
    // \brief The activation region of the discrete transition.
    boost::shared_ptr< const Geometry::SetInterface<R> > _activation; 

    // \brief The reset of the discrete transition.
    boost::shared_ptr< const MapInterface<R> > _reset;  
    
  public:
    /*! \brief Copy constructor. */
    DiscreteTransition(const DiscreteTransition<R> &original)
      : _event_id(original._event_id), _source(original._source), _destination(original._destination), 
        _activation(original._activation), _reset(original._reset) { 
    }
                
    /*! \brief Copy assignment operator. */
    DiscreteTransition<R>& operator=(const DiscreteTransition<R> &original) {
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
    const DiscreteMode<R> &source() const {
      return *this->_source;
    }      

    /*! \brief The destination of the discrete transition. */
    const DiscreteMode<R> &destination() const { 
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
    DiscreteTransition(id_type event_id, 
                       const DiscreteMode<R> &source, 
                       const DiscreteMode<R> &destination,
                       const MapInterface<R> &reset,
                       const Geometry::SetInterface<R> &activation)
      : _event_id(event_id), _source(&source), _destination(&destination), 
        _activation(activation.clone()), _reset(reset.clone()) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(activation,source,"DiscreteTransition::DiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(reset,source,"DiscreteTransition::DiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(reset,destination,"DiscreteTransition::DiscreteTransition(...)");
    }

    /* Construct from shared pointers (for internal use). */
    DiscreteTransition(id_type event_id,
                       const DiscreteMode<R> &source, 
                       const DiscreteMode<R> &destination,
                       const boost::shared_ptr< MapInterface<R> > reset,
                       const boost::shared_ptr< Geometry::SetInterface<R> > activation) 
      : _event_id(event_id), _source(&source), _destination(&destination), 
        _activation(activation), _reset(reset) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*activation,source,"DiscreteTransition::DiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(*reset,source,"DiscreteTransition::DiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(*reset,destination,"DiscreteTransition::DiscreteTransition(...)");
    }

    /* Construct from shared pointers (for internal use). */
    DiscreteTransition(id_type event_id,
                       const boost::shared_ptr< DiscreteMode<R> > source, 
                       const boost::shared_ptr< DiscreteMode<R> > destination,
                       const boost::shared_ptr< MapInterface<R> > reset,
                       const boost::shared_ptr< Geometry::SetInterface<R> > activation) 
      : _event_id(event_id), _source(&*source), _destination(&*destination), 
        _activation(activation), _reset(reset) 
    { 
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*activation,*source,"DiscreteTransition::DiscreteTransition(...)");
      ARIADNE_CHECK_ARGUMENT_DIMENSION(*reset,*source,"DiscreteTransition::DiscreteTransition(...)");
      ARIADNE_CHECK_RESULT_DIMENSION(*reset,*destination,"DiscreteTransition::DiscreteTransition(...)");
    }

};


template<class R> inline
std::ostream& 
DiscreteTransition<R>::write(std::ostream& os) const 
{ 
  return os << "DiscreteTransition( "
            << "id=" << this->id() << ", " 
            << "source_id=" << this->source().id() << ", "
            << "destination_id=" << this->destination().id() << ", "
            << "activation=" << this->activation() << ", "
            << "reset=" << this->reset() << " )"; 
}


template<class R> inline
std::ostream& operator<<(std::ostream& os, const DiscreteTransition<R>& dt) 
{
  return dt.write(os);
}

template<class R> inline
bool operator<(const DiscreteTransition<R>& transition1, const DiscreteTransition<R>& transition2) 
{
  return transition1.id() < transition2.id()
    || (transition1.id() == transition2.id() 
          && transition1.source().id() < transition2.source().id());
}


}
}
 
#endif /* ARIADNE_DISCRETE_TRANSITION_H */
