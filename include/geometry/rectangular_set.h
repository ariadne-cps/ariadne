/***************************************************************************
 *            rectangular_set.h
 *
 *  Copyright  2007-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file rectangular.h
 *  \brief Wrapper for %Box class conforming to the set interface.
 */

#ifndef ARIADNE_RECTANGULAR_SET_H
#define ARIADNE_RECTANGULAR_SET_H

#include <iosfwd>

#include "base/tribool.h"

#include "function/identity_function.h"
#include "geometry/image_set.h"
#include "geometry/constraint_set.h"
#include "geometry/box.h"

namespace Ariadne {
  
 
    
    //! \ingroup ExactSet
    /*! \brief A adaptor for the Rectangle class conforming to the SetInterface interface. */
    template<class R>
    class RectangularSet
      : public ConstraintSet<R>
    {
     public:
      /*! \brief */
      RectangularSet(const std::string& str) {
        Box<R> codomain(str); 
        *this=ConstraintSet<R>(IdentityFunction<R>(codomain.dimension()),codomain); }
      /*! \brief */
      template<class R1> RectangularSet(const dimension_type& d, const R1 a[][2])
        : ConstraintSet<R>(IdentityFunction<R>(d),Box<R>(d,a)) { }
      /*! \brief */
      template<class R1> RectangularSet(const Point<R1>& pt)
        : ConstraintSet<R>(IdentityFunction<R>(pt.dimension()),Box<R>(pt)) { }
      /*! \brief */
      template<class R1> RectangularSet(const Box<R1>& bx)
        : ConstraintSet<R>(IdentityFunction<R>(bx.dimension()),Box<R>(bx)) { }
      /*! \brief */
      operator Box<R> () const { return this->codomain(); }      

      /*! \brief */
      virtual ~RectangularSet<R>() { }
      /*! \brief */
      virtual RectangularSet<R>* clone() const { return new RectangularSet<R>(*this); }
      /*! \brief */
      virtual dimension_type dimension() const { return this->codomain().dimension(); }
      /*! \brief */
      virtual tribool contains(const Point<R>& pt) const { return this->codomain().contains(pt); }
      /*! \brief */
      virtual tribool superset(const Box<R>& r) const { 
        return Ariadne::superset(this->codomain(),r); }
      /*! \brief */
      virtual tribool intersects(const Box<R>& r) const { 
        return Ariadne::intersect(this->codomain(),r); }
      /*! \brief */
      virtual tribool disjoint(const Box<R>& r) const { 
        return Ariadne::disjoint(this->codomain(),r); }
      /*! \brief */
      virtual tribool subset(const Box<R>& r) const { 
        return Ariadne::subset(this->codomain(),r); }
      /*! \brief */
      virtual tribool bounded() const { return this->codomain().bounded(); }      
      /*! \brief */
      virtual Box<R> bounding_box() const { return this->codomain(); }      
      /*! \brief */
      virtual std::ostream& write(std::ostream& os) const {
        return os << "RectangularSet(\"" << this->codomain() << "\")";
      }
    };
    

  
} // namespace Ariadne

#endif /* ARIADNE_RECTANGULAR_SET_H */
