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

#ifndef _ARIADNE_DISCRETE_TRANSITION_H
#define _ARIADNE_DISCRETE_TRANSITION_H

#include <stdexcept>
#include <memory>

#include <boost/smart_ptr.hpp>

#include "../exceptions.h"
#include "../geometry/set.h"
#include "../system/map.h"
#include "../system/vector_field.h"
#include "../system/discrete_mode.h"

namespace Ariadne {
namespace System {
    
template< class R > class HybridAutomaton;


/*! \ingroup HybridTime
 *  \brief A discrete transition of a HybridAutomaton.
 */
template< class R >
class DiscreteTransition
{
  friend class HybridAutomaton<R>;
  private:
    static id_type _next_transition_id;
    
    // \brief The discrete transition's identificator.
    id_type _id;
  
    // \brief The source of the discrete transition.
    const DiscreteMode<R>* _source;   
  
    // \brief The destination of the discrete transition.
    const DiscreteMode<R>* _destination;   
  
    // \brief The activation region of the discrete transition.
    boost::shared_ptr< const Geometry::Set<R> > _activation; 

    // \brief The reset of the discrete transition.
    boost::shared_ptr< const Map<R> > _reset;  
    
  public:
    /*! \brief Copy constructor.
     *
     * This constructor initializes the object of the  
     * leaving discrete transition class.
     * \param orig is the original copy of the leaving arc.
     */
    DiscreteTransition(const DiscreteTransition<R> &orig)
      : _id(orig._id), _source(orig._source), _destination(orig._destination), 
        _activation(orig._activation), _reset(orig._reset) { 
    }
                
    /*! \brief Copy assignment operator for discrete transition.
     *
     * This method copies a leaving arc.
     * \param orig is the original copy of the leaving arc.
     */
    inline DiscreteTransition<R>& operator=(const DiscreteTransition<R> &orig) {
      this->_id=orig._id;
      this->_activation=orig._activation;
      this->_destination=orig._destination;
      this->_reset=orig._reset;
      this->_source=orig._source;
      return *this;
    }
    
    /*! \brief Return the unique identifyer of the discrete transition. */
    inline const id_type &id() const {
      return this->_id;
    }      

    /*! \brief Return the source of the discrete transition.
     *
     * This method return the source of the discrete transition.
     * \return The source of the discrete transition.
     */
    inline const DiscreteMode<R> &source() const {
      return *this->_source;
    }      

      /*! \brief Return the destination of the discrete transition.
     * 
     * This method return the destination of the discrete 
     * transition.
     * \return The destination of the discrete transition.
     */
    inline const DiscreteMode<R> &destination() const { 
      return *this->_destination;
    }
  
    /*! \brief Return the activation region of the discrete 
     * transition.
     * 
     * This method return the activation region of the leaving 
     * discrete transition.
     * \return The activation region of the discrete transition.
     */
    inline const Geometry::Set<R> &activation() const { 
      return *this->_activation;
    }

    /*! \brief Return the reset.
     * 
     * This method return the reset function of the leaving 
     * discrete transition.
     * \return The reset function of the leaving discrete 
     * transition.
     */
    inline const Map<R> &reset() const { 
      return *this->_reset;
    }
    
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
  private:
    /*! \brief This is a discrete transition class constructor.
     *
     * This constructor initializes the object of the discrete 
     * transition class.
     * @see LeavingDiscreteTransition()
     * \param reset is the reset relation of the discrete 
     * transition.
     * \param act is the activation region of the 
     * discrete transition.
     * \param source is the source mode of the discrete 
     * transition.
     * \param dest is the destination mode of the discrete 
     * transition.
      */
    DiscreteTransition(id_type id,
                       const Map<R> &reset,
                       const Geometry::Set<R> &act, 
                       const DiscreteMode<R> &source, 
                       const DiscreteMode<R> &dest)
      : _id(id), _source(&source), _destination(&dest), 
        _activation(act.clone()), _reset(reset.clone()) 
    { 
      check_equal_dimensions(act,source);
      check_argument_dimension(reset,source);
      check_result_dimension(reset,dest);
    }

    // Construct from shared pointers (for internal use)
    DiscreteTransition(id_type id,
                       const boost::shared_ptr< Map<R> > reset,
                       const boost::shared_ptr< Geometry::Set<R> > act, 
                       const DiscreteMode<R> &source, 
                       const DiscreteMode<R> &dest)
      : _id(id), _source(&source), _destination(&dest), 
        _activation(act), _reset(reset) 
    { 
      check_dimension(act,source);
      check_argument_dimension(reset,source);
      check_result_dimension(reset,dest);
    }

};


template<class R> inline
std::ostream& 
DiscreteTransition<R>::write(std::ostream& os) const 
{ 
  return os << "DiscreteTransition( "
            << "id=" << this->id() << ", " 
            << "source=" << this->source().id() << ", "
            << "destination=" << this->destination().id() << ", "
            << "activation=" << this->activation() << ", "
            << "reset=" << this->reset() << " )"; 
}


template<class R> inline
std::ostream& operator<<(std::ostream& os, DiscreteTransition<R>& dt) {
  return dt.write(os);
}

}
}
 
#endif /* _ARIADNE_DISCRETE_TRANSITION_H */
