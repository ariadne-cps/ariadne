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
    class RectangularSet : public SetInterface<R>
    {
     public:
      /*! \brief */
      RectangularSet(const std::string& str)
        : SetInterface<R>(), _rectangle(str) { }
      /*! \brief */
      template<class R1> RectangularSet(const Rectangle<R1>& r)
        : SetInterface<R>(), _rectangle(r) { }
      /*! \brief */
      template<class R1> RectangularSet(const RectangularSet<R>& rs)
        : SetInterface<R>(), _rectangle(static_cast<const Rectangle<R1>&>(rs)) { }
      /*! \brief */
      operator const Rectangle<R>& () const { return this->_rectangle; }

      /*! \brief */
      virtual ~RectangularSet<R>() { }
      /*! \brief */
      virtual RectangularSet<R>* clone() const { return new RectangularSet<R>(this->_rectangle); }
      /*! \brief */
      virtual dimension_type dimension() const { return this->_rectangle.dimension(); }
      /*! \brief */
      virtual tribool contains(const Point<R>& pt) const { return this->_rectangle.contains(pt); }
      /*! \brief */
      virtual tribool superset(const Rectangle<R>& r) const { 
        return Geometry::subset(r,this->_rectangle); }
      /*! \brief */
      virtual tribool intersects(const Rectangle<R>& r) const { 
        return !Geometry::disjoint(r,this->_rectangle); }
      /*! \brief */
      virtual tribool disjoint(const Rectangle<R>& r) const { 
        return Geometry::disjoint(r,this->_rectangle); }
      /*! \brief */
      virtual tribool subset(const Rectangle<R>& r) const { 
        return Geometry::subset(this->_rectangle,r); }
      /*! \brief */
      virtual tribool bounded() const { return this->_rectangle.bounded(); }      
      /*! \brief */
      virtual Rectangle<R> bounding_box() const { return this->_rectangle.bounding_box(); }      
      /*! \brief */
      virtual std::ostream& write(std::ostream& os) const {
        return os << "RectangularSet(\"" << this->_rectangle << "\")";
      }
     private:
      Rectangle<R> _rectangle;
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const RectangularSet<R>& rset) {
      return rset.write(os);
    }
    
  }
}

#endif /* ARIADNE_RECTANGULAR_SET_H */
