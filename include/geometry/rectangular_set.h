/***************************************************************************
 *            rectangular_set.h
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
 
/*! \file rectangular.h
 *  \brief Wrapper for %Rectangle class conforming to the set interface.
 */

#ifndef ARIADNE_RECTANGULAR_SET_H
#define ARIADNE_RECTANGULAR_SET_H

#include <iosfwd>

#include "../base/tribool.h"

#include "../geometry/set_interface.h"
#include "../geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {
 
    
    //! \ingroup ExactSet
    /*! \brief An adaptor for the Rectangle class conforming to the SetInterface interface. */
    template<class R>
    class RectangularSet : public SetInterface<R>, public Rectangle<R>
    {
     public:
      RectangularSet(const Rectangle<R>& r)
        : SetInterface<R>(), Rectangle<R>(r) { }
      
      virtual ~RectangularSet<R>() { }
      virtual RectangularSet<R>* clone() const { return new RectangularSet<R>(*this); }
      virtual dimension_type dimension() const { return Rectangle<R>::dimension(); }
      virtual tribool contains(const Point<R>& pt) const { return Rectangle<R>::contains(pt); }
      virtual Rectangle<R> bounding_box() const { return Rectangle<R>::bounding_box(); }      
      virtual tribool disjoint(const Rectangle<R>& r) const { 
        return Geometry::disjoint(r,static_cast<const Rectangle<R>&>(*this)); }
      virtual tribool superset(const Rectangle<R>& r) const { 
        return Geometry::subset(r,static_cast<const Rectangle<R>&>(*this)); }
      virtual tribool subset(const Rectangle<R>& r) const { 
        return Geometry::subset(static_cast<const Rectangle<R>&>(*this),r); }
      virtual std::ostream& write(std::ostream& os) const {
        return Rectangle<R>::write(os);
      }
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const RectangularSet<R>& rset) {
      return rset.write(os);
    }
    
  }
}

#endif /* ARIADNE_RECTANGULAR_SET_H */
