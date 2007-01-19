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
    inline id_type id() const {
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
    
    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;
    
   private:
    /*! \brief Construct a discrete mode.
     *  
     * This constructor initializes the object of the 
     * discrete mode class.
     * \param name is the name of the discrete mode.
     * \param id is the identifier of the mode.
     * \param dynamic is the mode's vector field.
     * \param invariant is the mode's invariant.
     */
    DiscreteMode(const std::string &name, 
                 const id_type& id,
                 const VectorField<R> &dynamic, 
                 const Geometry::Set<R> &invariant)
      : _name(name), _id(id), _dynamic(dynamic.clone()), _invariant(invariant.clone()) 
    {
      check_dimension(dynamic,invariant);
    }
      
    /*! \brief Construct an anonymous discrete mode.
     *  
     * This constructor initializes the object of the 
     * discrete mode class.
     * \param id is the identifier of the mode.
     * \param dynamic is the mode's vector field.
     * \param invariant is the mode's invariant.
     */
    DiscreteMode(id_type id,
                 const VectorField<R> &dynamic, 
                 const Geometry::Set<R> &invariant)
      :  _id(id), _name(), _dynamic(dynamic.clone()), _invariant(invariant.clone()) 
    {
      check_equal_dimensions(dynamic,invariant);
    }
    
    // Construct from objects managed by shared pointers (for internal use)
    DiscreteMode(id_type id,
                 const boost::shared_ptr< VectorField<R> > dynamic, 
                 const boost::shared_ptr< Geometry::Set<R> > invariant)
      :  _id(id), _name(),_dynamic(dynamic), _invariant(invariant) 
    {
      check_dimension(dynamic,invariant);
    }
    
};
  
template<class R> inline
std::ostream& 
DiscreteMode<R>::write(std::ostream& os) const 
{ 
  return os << "DiscreteMode( "
            << "id=" << this->id() << ", " 
            << "invariant=" << this->invariant() << ", "
            << "dynamic=" << this->dynamic() << " )"; 
}

template<class R> inline
std::ostream& operator<<(std::ostream& os, const DiscreteMode<R>& dm) {
  return dm.write(os);
}
    
}
}

#endif /* _ARIADNE_DISCRETE_MODE_H */
