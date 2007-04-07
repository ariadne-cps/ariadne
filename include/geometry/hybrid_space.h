/***************************************************************************
 *            hybrid_space.h
 *
 *  Copyright  2007  Alberto Casagrande,  Pieter Collins
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
 
#ifndef ARIADNE_HYBRID_SPACE_H
#define ARIADNE_HYBRID_SPACE_H

#include <set>
#include <map>

#include <string>
#include <iostream>
#include <sstream>

#include "../base/types.h"
#include "exceptions.h"

namespace Ariadne {  
  namespace Geometry {

    /*! \brief The type identifying a discrete locatation of a hybrid system. */
    typedef id_type location_type;


    /*! \ingroup HybridSet
     *  \brief A location of a hybrid system, with an identifier and a dimension.
     */
    class HybridLocation
    {
     public:
      /*! \brief Construct from an identifier and a dimension. */
      HybridLocation(id_type id, dimension_type d)
        : _id(id), _dimension(d) { }
      /*! \brief The indentifier of the location. */
      id_type id() const { return this->_id; }
      /*! \brief The dimension of the location. */
      dimension_type dimension() const { return this->_dimension; }
      /*! \brief Equality operator. 
       * Note that it is an error to have two locations with the same identifier in a given system. */
      bool operator==(const HybridLocation& loc) const {
        return this->_id==loc._id && this->_dimension==loc._dimension;
      }
      /*! \brief Inequality operator. 
       * Note that it is an error to have two locations with the same identifier in a given system. */
      bool operator!=(const HybridLocation& loc) const {
        return !(*this==loc);
      }
      /*! \brief Comparison operator. */
      bool operator<(const HybridLocation& loc) const {
        return this->_id<loc._id;
      }
     private:
      id_type _id;
      dimension_type _dimension;
    };

    inline std::ostream& operator<<(std::ostream& os, const HybridLocation& hloc) {
      return os << "HybridLocation( id=" << hloc.id() << ", dimension=" << hloc.dimension() << " )";
    }


    /*! \ingroup HybridSet
     *  \brief A space of a hybrid system, giving the number of dimensions in each mode.
     */
    class HybridSpace
    {
     public:
      typedef std::set<HybridLocation>::const_iterator iterator;
      typedef std::set<HybridLocation>::const_iterator const_iterator;
      /*! \brief Default constructor. */
      HybridSpace() : _locations() { }
      /*! \brief Copy constructor. */
      HybridSpace(const std::map<id_type,dimension_type>& hyspcmap) : _locations() { 
        for(std::map<id_type,dimension_type>::const_iterator iter=hyspcmap.begin();
            iter!=hyspcmap.end(); ++iter)
        {
          this->_locations.insert(HybridLocation(iter->first,iter->second));
        }
      }
      /*! \brief Copy constructor. */
      HybridSpace(const HybridSpace& hspc) : _locations(hspc._locations) { }
      /*! \brief Equality operator. */
      bool operator==(const HybridSpace& hspc) const {
        return this->_locations==hspc._locations;
      }
      bool operator!=(const HybridSpace& hspc) const {
        return !(*this==hspc);
      }
      /*! The number of locations. */
      size_type number_of_locations() const { return this->_locations.size(); }
      /*! Add a new location. */
      void new_location(id_type id, dimension_type d);
      /*! Add a new location. */
      void new_location(const HybridLocation& loc);
      /*! \brief Check if the space has a location. */
      bool has_location(id_type id) const;
      /*! \brief The location with identifier \a id. */
      HybridLocation operator[](id_type id) const;
      /*! \brief The dimension of the location. */
      dimension_type dimension(id_type id) const;
      /*! \brief The set of locations (returns self). */
      const HybridSpace& locations() const {
        return *this;
      }
      /*! \brief A constant iterator to the beginning of the set of locations. */
      const_iterator begin() const {
        return this->_locations.begin();
      }
      /*! \brief A constant iterator to the end of the set of locations. */
      const_iterator end() const {
        return this->_locations.end();
      }
     private:
      std::set<HybridLocation> _locations;
    };


    inline std::ostream& operator<<(std::ostream& os, const HybridSpace& hspc) {
      os << "HybridSpace( ";
      HybridSpace::const_iterator iter=hspc.begin();
      if(iter!=hspc.end()) {
        os << iter->id() << ":" << iter->dimension();
        ++iter;
        for( ; iter!=hspc.end(); ++iter) {
          os << ", " << iter->id() << ":" << iter->dimension();
        }
      }
      os << " )";
      return os;
    }


    inline void HybridSpace::new_location(const HybridLocation& loc) {
      ARIADNE_CHECK_NEW_LOCATION(*this,loc.id(),"void HybridSpace::new_location(HybridLocation loc)");
      this->_locations.insert(loc);
    }

    inline void HybridSpace::new_location(id_type id, dimension_type d) {
      ARIADNE_CHECK_NEW_LOCATION(*this,id,"void HybridSpace::new_location(id_type id, dimension_type d)");
      this->_locations.insert(HybridLocation(id,d));
    }

    inline bool HybridSpace::has_location(id_type id) const
    {
      const_iterator iter=this->_locations.find(HybridLocation(id,0u));
      return iter!=this->_locations.end();
    }
  
    inline HybridLocation HybridSpace::operator[](id_type id) const {
      ARIADNE_CHECK_LOCATION(*this,id,"HybridLocation HybridSpace::operator[](id_type id)");
      const_iterator iter=this->_locations.find(HybridLocation(id,0u));
      return *iter;
    }
  
    inline dimension_type HybridSpace::dimension(id_type id) const {
      ARIADNE_CHECK_LOCATION(*this,id,"dimension_type HybridSpace::dimension(id_type id)");
      const_iterator iter=this->_locations.find(HybridLocation(id,0u));
      return iter->dimension();
    }

    
  }
}

#endif /* ARIADNE_HYBRID_SPACE_H */
