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

#include "base/types.h"
#include "numeric/traits.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"

namespace Ariadne {
  namespace Geometry {
    
    template<class R> class Point;
    template<class R> class Box;
    template<class R> class Polyhedron;
    template<class R0,class R1> class Zonotope;

  
    // Forward declarations for friends
    template<class R> class CurveInterface;

    //! \ingroup SetInterface
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
      /*! \brief Return a new dynamically-allocated copy of the curve. */
      virtual CurveInterface<R>* clone() const = 0;
      /*! \brief The dimension of the space the curve lies in. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The smoothness of the curve. */
      virtual smoothness_type smoothness() const = 0;

      /*! \brief The point on the curve at a parameter value. */
      virtual Point<A> value(const A& s) const = 0;
      /*! \brief The tangent vector to the curve at a parameter value. */
      virtual LinearAlgebra::Vector<A> tangent(const A& s) const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
    };
    

    template<class R> CurveInterface<R>::~CurveInterface() {
    }
    
    template<class R> inline std::ostream& operator<<(std::ostream& os, const CurveInterface<R>& c) {
      return c.write(os);
    }
    
  }
}


#endif /* ARIADNE_CURVE_INTERFACE_H */
