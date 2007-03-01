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
 
#ifndef ARIADNE_DISCRETE_MODE_H
#define ARIADNE_DISCRETE_MODE_H

#include <string>
#include <stdexcept>
#include <boost/smart_ptr.hpp>

#include "vector_field.h"
#include "../exceptions.h"
#include "../geometry/set_interface.h"

namespace Ariadne {
namespace System {

template< class R > class VectorField;
template< class R > class HybridAutomaton;
  

/*!\ingroup HybridTime
 * \brief A discrete mode of a HybridAutomaton, comprising continuous evolution given by a VectorField
 * within and invariant Geometry::SetInterface. 
 *
 * A %DiscreteMode can only be created using the new_mode() method in
 * the %HybridAutomaton class.
 */
template<class R>
class DiscreteMode {
    friend class HybridAutomaton<R>;
  public:
      /*!\brief The type of denotable real number used to describe to discrete mode. */
    typedef R real_type;
      
  private:
    
    // The discrete mode's identificator.
    id_type _id;
  
    // The discrete mode's vector field.
    boost::shared_ptr< const VectorField<R> > _dynamic;
  
    // The discrete mode's invariant.
    boost::shared_ptr< const Geometry::SetInterface<R> > _invariant;
  
  public:
    
    /*! \brief Copy constructor. */
    DiscreteMode(const DiscreteMode<R>& original)
      : _id(original._id), _dynamic(original._dynamic),
        _invariant(original._invariant) {}
    
    /*! \brief Copy assignment operator. */
    DiscreteMode<R>& operator=(const DiscreteMode<R> &original) {
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
    const VectorField<R>& dynamic() const {
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
    DiscreteMode(id_type id,
                 const VectorField<R> &dynamic, 
                 const Geometry::SetInterface<R> &invariant)
      :  _id(id), _dynamic(dynamic.clone()), _invariant(invariant.clone()) 
    {
      check_equal_dimensions(dynamic,invariant);
    }
    
    /* Construct from objects managed by shared pointers (for internal use) */
    DiscreteMode(id_type id,
                 const boost::shared_ptr< VectorField<R> > dynamic, 
                 const boost::shared_ptr< Geometry::SetInterface<R> > invariant)
      :  _id(id), _dynamic(dynamic), _invariant(invariant) 
    {
      check_equal_dimensions(*dynamic,*invariant);
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
std::ostream& operator<<(std::ostream& os, const DiscreteMode<R>& dm) 
{
  return dm.write(os);
}

template<class R> inline
bool operator<(const DiscreteMode<R>& mode1, const DiscreteMode<R>& mode2) 
{
  return mode1.id() < mode2.id();
}


}
}

#endif /* ARIADNE_DISCRETE_MODE_H */
