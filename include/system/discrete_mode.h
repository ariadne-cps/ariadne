/***************************************************************************
 *            discrete_mode.h
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
 
#ifndef _ARIADNE_DISCRETE_MODE_H
#define _ARIADNE_DISCRETE_MODE_H

#include <string>
#include <stdexcept>
#include <boost/smart_ptr.hpp>

#include "vector_field.h"
#include "../geometry/set.h"

namespace Ariadne {

namespace Geometry {
template<class R> class Set;
}

namespace System {

template< class R > class VectorField;
template< class R > class HybridAutomaton;
  

/*!\ingroup HybridTime
 * \brief A discrete mode or mode of a HybridAutomaton. 
 */
template<class R>
class DiscreteMode {
    friend class HybridAutomaton<R>;
  public:
      /*!\brief The type of denotable real number used to describe to discrete mode. */
    typedef R real_type;
      
  private:

    static id_type _next_mode_id;
    
    // The discrete mode's identificator.
    id_type _id;
  
    // The discrete mode's name.
    std::string _name;
  
    // The discrete mode's vector field.
    boost::shared_ptr< const VectorField<R> > _dynamic;
  
    // The discrete mode's invariant.
    boost::shared_ptr< const Geometry::Set<R> > _invariant;
  
  public:  
  
    /*! \brief Construct a discrete mode.
     *  
     * This constructor initializes the object of the 
     * discrete mode class.
     * \param name is the name of the discrete mode.
     * \param dynamic is the mode's vector field.
     * \param invariant is the mode's invariant.
     */
    DiscreteMode(const std::string &name, 
                     const VectorField<R> &dynamic, 
                     const Geometry::Set<R> &invariant)
      : _name(name), _dynamic(dynamic.clone()), _invariant(invariant.clone()) 
    {
      check_dimension(dynamic,invariant);
      this->_set_id();
    }
      
    /*! \brief Copy constructor. */
    DiscreteMode(const DiscreteMode<R>& orig)
      : _id(orig._id), _name(orig._name), _dynamic(orig._dynamic),
        _invariant(orig._invariant) {}
    
    /*! \brief Copy assignment operator. */
    inline DiscreteMode<R>& operator=(const DiscreteMode<R> &orig) {
      if(this!=&orig) {
        this->_id=orig._id;
        this->_name=orig._name;
        this->_dynamic=orig._dynamic;
        this->_invariant=orig._invariant;
      }
      return *this;      
    }

    /*! \brief The discrete mode's identifier. */
    inline const id_type& id() const {
      return this->_id;
    }
    
    /*! \brief The discrete mode's name. */
    inline const std::string& name() const {
      return this->_name;  
    }
    
    /*! \brief The discrete mode's dynamic (a vector field). */
    inline const VectorField<R>& dynamic() const {
      return *this->_dynamic;  
    }
    
    /*! \brief The discrete mode's invariant. */
    inline const Geometry::Set<R>& invariant() const{
      return *this->_invariant;  
    }
    
    /*! \brief The dimension of the discrete mode. */
    dimension_type dimension() const { return this->_invariant->dimension(); }
    

    
   private:
    // Set the identifier of the mode.
    void _set_id() { this->_id=_next_mode_id; ++_next_mode_id; }
};

template<class R> id_type DiscreteMode<R>::_next_mode_id=0;
  
}
}

#endif /* _ARIADNE_DISCRETE_MODE_H */
