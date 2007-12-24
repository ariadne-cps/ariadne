/***************************************************************************
 *            empty_set.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file empty_set.h
 *  \brief An empty set corresponding to the Set interface.
 */

#ifndef ARIADNE_EMPTY_SET_H
#define ARIADNE_EMPTY_SET_H

#include <iosfwd>

#include "base/tribool.h"
#include "numeric/interval.h"
#include "geometry/set_interface.h"
#include "geometry/rectangle.h"


namespace Ariadne {
  namespace Geometry {

    //! \ingroup ExactSet
    /*! \brief A class for empty sets of a given dimension conforming to the Set interface. 
     *
     *  \internal This class is experimental and may be removed at some point.
     */
    template<class R>
    class EmptySet
      : public SetInterface<R> 
    {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
      typedef Box<R> basic_set_type;
     
      /*! \brief An empty set in \a d dimensions. */
      EmptySet(dimension_type d) : _dimension(d) { }
     
      /*! \brief A dynamically-allocated copy of the set. */
      virtual EmptySet<R>* clone() const { return new EmptySet(this->dimension()); }
     
      /*! \brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const { return this->_dimension; }
      
      /*! \brief Tests if the set contains a point. */
      virtual tribool contains(const Point<R>& pt) const { return false; }
     
      /*! \brief Tests if the set is a superset of a rectangle. */
      virtual tribool superset(const Box<R>& r) const { return r.empty(); }
      /*! \brief Tests if the set intersects a rectangle. */
      virtual tribool intersects(const Box<R>&) const { return false; }
      /*! \brief Tests if the set is disjoint from a rectangle. */
      virtual tribool disjoint(const Box<R>&) const { return true; }
      /*! \brief Tests if the set is a subset of a rectangle. */
      virtual tribool subset(const Box<R>&) const { return true; }
      
      /*! \brief A rectangle containing the set. Returns the default rectangle of the given dimension. */
      virtual tribool bounded() const { return true; }
      /*! \brief A rectangle containing the set. Returns the default rectangle of the given dimension. */
      virtual Box<R> bounding_box() const { return Box<R>(array< Numeric::Interval<R> >(this->dimension())); }


      /*! \brief Write to an output stream. 
       *  Called by operator<<(std::ostream&, const SetInterface<R>&) to dynamically dispatch stream output. 
       */
      virtual std::ostream& write(std::ostream& os) const { return os << "EmptySet( dimension=" << this->dimension() << " )"; }
     private:
      dimension_type _dimension;
    };
    
  }
}


#endif /* ARIADNE_EMPTY_SET_H */
