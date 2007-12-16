/***************************************************************************
 *            geometry/exceptions.h
 *
 *  Copyright  2005-7  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
 
/*! \file geometry/exceptions.h
 *  \brief Exceptions, error handling and assertions for the Geometry module.
 */

#ifndef ARIADNE_GEOMETRY_EXCEPTIONS_H
#define ARIADNE_GEOMETRY_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "macros/throw.h"
#include "base/types.h"

namespace Ariadne {

  namespace Geometry {
  
    //@{ \name Exceptions

    /*! \brief The dimensions of two geometric objects are incompatible. */
    struct IncompatibleDimensions : public std::runtime_error {
      IncompatibleDimensions(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The coordinate in  a geometric object was invalid. */
    struct InvalidCoordinate : public std::runtime_error {
      InvalidCoordinate(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The value of an element of a geometric object was invalid. */
    struct InvalidValue : public std::runtime_error {
      InvalidValue(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief An invalid index to the set of vertices. */
    struct InvalidVertex : public std::runtime_error {
      InvalidVertex(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The generaors of a zonotope or polytope were of the wrong size. */
    struct InvalidGenerators : public std::runtime_error {
      InvalidGenerators(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief Attempting to perform a binary operation on two objects on incompatible grids. */
    struct IncompatibleGrids : public std::runtime_error {
      IncompatibleGrids(const std::string& str) : std::runtime_error(str) { }
    };
      
    /*! \brief Attempting to perform a binary operation on two objects on incompatible grids. */
    struct InvalidGridPosition : public std::runtime_error {
      InvalidGridPosition(const std::string& str) : std::runtime_error(str) { }
    };
      
    /*! \brief A location of a hybrid set or system was invalid. */
    struct InvalidLocation : public std::runtime_error {
      InvalidLocation(const std::string& str) : std::runtime_error(str) { }
    };
      
    /*! \brief The locations of two hybrid sets/systems are incompatible. */
    struct IncompatibleLocations : public std::runtime_error {
      IncompatibleLocations(const std::string& str) : std::runtime_error(str) { }
    };
      
    /*! \brief Attempting to perform an operation on an unbounded set. */
    struct UnboundedSet : public std::runtime_error {
      UnboundedSet(const std::string& str) : std::runtime_error(str) { }
    };
     
    /*! \brief Attempting to perform an invalid operation on a set with empty interior. */
    struct EmptyInterior : public std::runtime_error {
      EmptyInterior(const std::string& str) : std::runtime_error(str) { }
    };
     
  }
}

#define ARIADNE_CHECK_DIMENSION(obj,dim,func)                          \
  { if((obj).dimension()!=dim) { using namespace Geometry; ARIADNE_THROW(IncompatibleDimensions,func,#obj"="<<obj<<", "#dim"="<<dim); } }
        
#define ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(obj1,obj2,func)                          \
  { if((obj1).dimension()!=(obj2).size()) { using namespace Geometry; ARIADNE_THROW(IncompatibleDimensions,func,#obj1"="<<obj1<<", "#obj2"="<<obj2); } }
        
#define ARIADNE_CHECK_EQUAL_DIMENSIONS(obj1,obj2,func)                  \
  { if((obj1).dimension()!=(obj2).dimension()) { using namespace Geometry; ARIADNE_THROW(IncompatibleDimensions,func,#obj1"="<<obj1<<", "#obj2"="<<obj2); } }
        
#define ARIADNE_CHECK_COORDINATE(obj,ind,func)                          \
  { if((obj).dimension()<=ind) { using namespace Geometry; ARIADNE_THROW(InvalidCoordinate,func,#obj"="<<obj<<", "#ind"="<<ind); } }

#define ARIADNE_CHECK_VERTEX_INDEX(poly,ind,func)                       \
  { if((poly).number_of_vertices()<=ind) { using namespace Geometry; ARIADNE_THROW(InvalidVertex,func,#ind"="<<ind<<" but "#poly"has "<<(poly).number_of_vertices()<<" vertices"); } }

#define ARIADNE_CHECK_SAME_GRID(set1,set2,func)                       \
  { if((set1).grid()!=(set2).grid()) { using namespace Geometry; ARIADNE_THROW(IncompatibleGrids,func,#set1"="<<set1<<", "#set2"="<<set2); } }

#define ARIADNE_CHECK_BOUNDED(set,func)                                 \
  { if(!(set).bounded()) { using namespace Geometry; ARIADNE_THROW(UnboundedSet,func,#set"="<<set); } }


#define ARIADNE_CHECK_NEW_LOCATION(hybr,loc,func)                              \
  { if((hybr).has_location(loc)) { using namespace Geometry; ARIADNE_THROW(InvalidLocation,func,#hybr".locations()="<<(hybr).locations()<<" already contains loc="<<loc); } }

#define ARIADNE_CHECK_LOCATION(hybr,loc,func)                              \
  { if(!(hybr).has_location(loc)) { using namespace Geometry; ARIADNE_THROW(InvalidLocation,func,#hybr".locations()="<<(hybr).locations()<<", "#loc"="<<loc); } }

#define ARIADNE_CHECK_SAME_LOCATIONS(hybr1,hybr2,func)                 \
  { if((hybr1).locations()!=(hybr2).locations()) { using namespace Geometry; ARIADNE_THROW(IncompatibleLocations,func,#hybr1".locations()="<<(hybr1).locations()<<", "#hybr2".locations()="<<(hybr2).locations()); } }


#endif /* ARIADNE_GEOMETRY_EXCEPTIONS_H */
