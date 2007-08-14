/***************************************************************************
 *            curve_interface.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file curve_interface.h
 *  \brief An interface for curves in Euclidean space.
  */

#ifndef ARIADNE_CURVE_INTERFACE_H
#define ARIADNE_CURVE_INTERFACE_H

#include <iosfwd>

namespace Ariadne {
  namespace Geometry {
    
    template<class R> class Point;
    template<class R> class Rectangle;
    template<class R> class Polyhedron;
    template<class R0,class R1> class Zonotope;

  
    // Forward declarations for friends
    template<class R> class CurveInterface;

    //! \ingroup ExactSet
    /*! \brief A curve in Euclidean space
     */
    template<class R>
    class CurveInterface
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
     public:
      /*! \brief Destructor. */
      virtual ~CurveInterface();
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual CurveInterface<R>* clone() const = 0;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The smoothness of the constraint function. */
      virtual size_type smoothness() const = 0;

      /*! \brief The value at a point. */
      virtual Point<A> value(const A&& s) const = 0;
      /*! \brief The tangent at a point. */
      virtual LinearAlgebra::Vector<A> tangent(const A& s) const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
    };
    
    template<class R> std::ostream& operator<<(std::ostream& os, const CurveInterface<R>& c);
    
  }
}

#include "constraint.inline.h"

#endif /* ARIADNE_CONSTRAINT_H */
