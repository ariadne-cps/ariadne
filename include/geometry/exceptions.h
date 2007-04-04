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

#include "../throw.h"
#include "../base/types.h"

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
      
    /*! \brief Attempting to perform an operation on an unbounded set. */
    struct UnboundedSet : public std::runtime_error {
      UnboundedSet(const std::string& str) : std::runtime_error(str) { }
    };
     
    template<class S1> inline
    void check_dimension(const S1& s1, const dimension_type& d2, const char* where="") {
      if(s1.dimension()!=d2) { throw IncompatibleDimensions(where); }
    }
    
    template<class S1, class S2> inline
    void check_equal_dimensions(const S1& s1, const S2& s2, const char* where="") {
      if(s1.dimension()!=s2.dimension()) { throw IncompatibleDimensions(where); }
    }
    
    template<class S1> inline
    void check_coordinate(const S1& s1, const dimension_type& i2, const char* where="") {
      if(i2>=s1.dimension()) { throw InvalidCoordinate(where); }
    }

    template<class S1> inline
    void check_vertex_index(const S1& s1, const dimension_type& i2, const char* where="") {
      if(i2>=s1.number_of_vertices()) { throw InvalidVertex(where); }
    }

    template<class S1, class S2> inline
    void check_same_grid(const S1& s1, const S2& s2, const char* where="") {
      if(s1.grid()!=s2.grid()) { throw IncompatibleGrids(where); }
    }
    
    template<class S> inline
    void check_bounded(const S& s, const char* where="") {
      if(!s.bounded()) { throw UnboundedSet(where); }
    }



  }
}

#define ARIADNE_CHECK_DIMENSION(obj,dim,func)                          \
  { if((obj).dimension()!=dim) { ARIADNE_THROW(IncompatibleDimensions,func,#obj"="<<obj<<", "#dim"="<<dim); } }
        
#define ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(obj1,obj2,func)                          \
  { if((obj1).dimension()!=obj2.size()) { ARIADNE_THROW(IncompatibleDimensions,func,#obj1"="<<obj1<<", "#obj2"="<<obj2); } }
        
#define ARIADNE_CHECK_EQUAL_DIMENSIONS(obj1,obj2,func)                  \
  { if((obj1).dimension()!=obj2.dimension()) { ARIADNE_THROW(IncompatibleDimensions,func,#obj1"="<<obj1<<", "#obj2"="<<obj2); } }
        
#define ARIADNE_CHECK_COORDINATE(obj,ind,func)                          \
  { if((obj).dimension()>=ind) { ARIADNE_THROW(InvalidCoordinate,func,#obj"="<<obj<<", "#ind"="<<ind); } }

#define ARIADNE_CHECK_VERTEX_INDEX(poly,ind,func)                       \
  { if((poly).number_of_vertices()<=ind) { ARIADNE_THROW(InvalidVertex,func,#ind"="<<ind<<" but "#poly"has "<<poly.number_of_vertices()<<" vertices"); } }

#define ARIADNE_CHECK_SAME_GRID(poly,ind,func)                          \
  { if((poly).number_of_vertices()<=ind) { ARIADNE_THROW(InvalidVertex,func,#ind"="<<ind<<" but "#poly"has "<<poly.number_of_vertices()<<" vertices"); } }

#define ARIADNE_CHECK_BOUNDED(set,func)                                 \
  { if(!(set).bounded()) { ARIADNE_THROW(UnboundedSet,func,#set); } }


#endif /* ARIADNE_GEOMETRY_EXCEPTIONS_H */
