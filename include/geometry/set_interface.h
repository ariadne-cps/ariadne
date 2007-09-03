/***************************************************************************
 *            set_interface.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file set_interface.h
 *  \brief Interface for general sets.
 */

#ifndef ARIADNE_SET_INTERFACE_H
#define ARIADNE_SET_INTERFACE_H

#include <iosfwd>

#include "../base/types.h"
#include "../base/tribool.h"


namespace Ariadne {
  namespace Geometry {

    template<class R> class Point;
    template<class R> class Rectangle;
      
    //! \ingroup SetInterface
    /*! \brief An abstract base class for general sets. */
    template<class R>
    class SetInterface {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      /*! \brief Virtual destructor. */
      virtual ~SetInterface();
      /*! \brief A dynamically-allocated copy of the set. */
      virtual SetInterface<R>* clone() const = 0;
     
      /*! \brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const = 0;
      
      /*! \brief Tests if the set contains a point. */
      virtual tribool contains(const Point<R>&) const = 0;
     
      /*! \brief Tests if the set is a superset of a rectangle. */
      virtual tribool superset(const Rectangle<R>&) const = 0;
      /*! \brief Tests if the set intersects a rectangle. */
      virtual tribool intersects(const Rectangle<R>&) const = 0;
      /*! \brief Tests if the set is disjoint from a rectangle. */
      virtual tribool disjoint(const Rectangle<R>&) const = 0;
      /*! \brief Tests if the set is a subset of a rectangle. */
      virtual tribool subset(const Rectangle<R>&) const = 0;
      
      /*! \brief Tests if the set is bounded. */
      virtual tribool bounded() const = 0;
      /*! \brief A rectangle containing the set. Throws Geometry::UnboundedSet exception if the set is unbounded. */
      virtual Rectangle<R> bounding_box() const = 0;

      /*! \brief Write to an output stream. 
       *  Called by operator<<(std::ostream&, const SetInterface<R>&) to dynamically dispatch stream output. 
       */
      virtual std::ostream& write(std::ostream& os) const = 0;

#ifdef DOXYGEN
      //@{ 
      //! \name Binary geometric predicates
      /*! \brief Tests if the set \a s is disjoint from the rectangle \a r. */
      friend tribool disjoint(const SetInterface<R>& s, const Rectangle<R>& r);
      /*! \brief Tests if the rectangle \a r is disjoint from the set \a s. */
      friend tribool disjoint(const Rectangle<R>& r, const SetInterface<R>& s);
    
      /*! \brief Tests if the rectangle \a r is a subset of the set \a s. */
      friend tribool subset(const Rectangle<R>& r, const SetInterface<R>& s);
    
      /*! \brief Tests if the set \a s is a subset of the rectangle \a r. */
      friend tribool subset(const SetInterface<R>& s, const Rectangle<R>& r);
      //@}

      //@{ 
      //! \name Stream input/output 
      /*! \brief Write to an output stream. */
      friend std::ostream& operator<<(std::ostream& os, const SetInterface<R>& s);
      //@}
#endif
    };
 

    template<class R> 
    SetInterface<R>::~SetInterface() 
    {
    }


    template<class R> inline tribool disjoint(const SetInterface<R>& s, const Rectangle<R>& r) {
      return s.disjoint(r);
    }
    
    template<class R> inline tribool disjoint(const Rectangle<R>& r, const SetInterface<R>& s) {
      return s.disjoint(r);
    }
    
    template<class R> inline tribool subset(const Rectangle<R>& r, const SetInterface<R>& s) {
      return s.superset(r);
    }
    
    template<class R> inline tribool subset(const SetInterface<R>& s, const Rectangle<R>& r) {
      return s.subset(r);
    }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const SetInterface<R>& s) {
      return s.write(os);
    }
    
  }
}


#endif /* ARIADNE_SET_INTERFACE_H */
