/***************************************************************************
 *            constraint_interface.h
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
 
/*! \file constraint_interface.h
 *  \brief A constraint defining an invariant, activation or guard set.
  */

#ifndef ARIADNE_CONSTRAINT_INTERFACE_H
#define ARIADNE_CONSTRAINT_INTERFACE_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "../base/tribool.h"

#include "../system/function_interface.h"

namespace Ariadne {
  namespace Geometry {
    
    template<class R> class Point;
    template<class R> class Rectangle;
    template<class R0,class R1> class Zonotope;
    template<class R> class Polyhedron;

  
    // Forward declarations for friends
    template<class R> class ConstraintInterface;
    template<class R> class Constraint;
    template<class R> bool equal(const Constraint<R>& c1, const Constraint<R>& c2);
    template<class R> bool opposite(const Constraint<R>& c1, const Constraint<R>& c2);
  
    /*! \brief The type of comparison used to define the constraint. */
    enum Comparison { less, greater };


    //! \ingroup SetInterface
    /*! \brief A constraint on the state
     */
    template<class R>
    class ConstraintInterface
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
     public:
      /*! \brief Destructor. */
      virtual ~ConstraintInterface();
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual ConstraintInterface<R>* clone() const = 0;
      /*! \brief The dimension of the underlying state space. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The operation used for comparison. */
      virtual Comparison comparison() const = 0;
      /*! \brief The smoothness of the constraint function. */
      virtual smoothness_type smoothness() const = 0;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;

      /*! \brief The value of the constraint function at a point. */
      virtual A value(const Point<A>& pt) const = 0;
      /*! \brief The gradient of the constraint function at a point. */
      virtual LinearAlgebra::Vector<A> gradient(const Point<A>& pt) const = 0;
#ifdef DOXYGEN
      /*! \brief Test if the constraint is satisfied over a rectangle. */
      friend tribool satisfies(const Rectangle<R>& r, const ConstraintInterface<R>& c);
      /*! \brief Test if the constraint is satisfied over a zonotope. */
      friend tribool satisfies(const Zonotope<R,R>& z), const ConstraintInterface<R>& c);
      /*! \brief Test if the constraint is satisfied over a zonotope. */
      friend tribool satisfies(const Zonotope<I,R>& z), const ConstraintInterface<R>& c);
      /*! \brief Test if the constraint is satisfied over a zonotope. */
      friend tribool satisfies(const Zonotope<I,I>& z), const ConstraintInterface<R>& c);
#endif
    };
    
    template<class R> ConstraintInterface<R>::~ConstraintInterface() { 
    }
    
    template<class R> inline std::ostream& operator<<(std::ostream& os, const ConstraintInterface<R>& c) {
      return c.write(os);
    }
    
    template<class R> Numeric::Interval<R> value(const ConstraintInterface<R>& c, const Rectangle<R>& r);
    template<class R> Numeric::Interval<R> value(const ConstraintInterface<R>& c, const Zonotope< Numeric::Interval<R>, Numeric::Interval<R> >& z);

    template<class R> tribool satisfies(const Rectangle<R>& r, const ConstraintInterface<R>& c);
    template<class R> tribool satisfies(const Zonotope<R,R>& z, const ConstraintInterface<R>& c);
    template<class R> tribool satisfies(const Zonotope<Numeric::Interval<R>,R>& z, const ConstraintInterface<R>& c);
    template<class R> tribool satisfies(const Zonotope< Numeric::Interval<R>,Numeric::Interval<R> >& z, const ConstraintInterface<R>& c);

  }
}

#endif /* ARIADNE_CONSTRAINT_INTERFACE_H */
