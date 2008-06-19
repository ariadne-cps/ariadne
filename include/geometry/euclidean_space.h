/***************************************************************************
 *            euclidean_space.h
 *
 *  Copyright  2008  Pieter Collins
 *
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
 
#ifndef ARIADNE_EUCLIDEAN_SPACE_H
#define ARIADNE_EUCLIDEAN_SPACE_H

#include <set>
#include <map>

#include <string>
#include <iostream>
#include <sstream>

#include "base/types.h"

namespace Ariadne {  
  


    /*! \ingroup EuclideanSet
     *  \brief A space of a euclidean system, giving the dimensions.
     */
    class EuclideanSpace
    {
     public:
      /*! \brief Constructor. */
      EuclideanSpace(const dimension_type& d)
        : _dimension(d) { }
      /*! \brief Equality operator. */
      bool operator==(const EuclideanSpace& espc) const {
        return this->_dimension==espc._dimension; }
      bool operator!=(const EuclideanSpace& espc) const {
        return !(*this==espc); }
      /*! The dimension of the space. */
      dimension_type dimension() const { return this->_dimension; }
     private:
      dimension_type _dimension;
    };


    inline std::ostream& operator<<(std::ostream& os, const EuclideanSpace& espc) {
      return os << "EuclideanSpace( dimension=" << espc.dimension() << " )" << std::endl;
    }
    

} // namespace Ariadne

#endif /* ARIADNE_EUCLIDEAN_SPACE_H */
