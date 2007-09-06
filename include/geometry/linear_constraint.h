/***************************************************************************
 *            linear_constraint.h
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
 
/*! \file linear_constraint.h
 *  \brief A linear constraint defining an invariant, activation or guard set.
  */

#ifndef ARIADNE_LINEAR_CONSTRAINT_H
#define ARIADNE_LINEAR_CONSTRAINT_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "../base/tribool.h"

#include "../function/function_interface.h"
#include "constraint_interface.h"

namespace Ariadne {
  namespace Geometry {
    
    template<class R> class Point;
    template<class R> class Rectangle;
    template<class R0,class R1> class Zonotope;
    template<class R> class Polyhedron;

  

    /*! \brief A linear inequality constraint. */
    template<class R>
    class LinearConstraint
      : public ConstraintInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Construct the constraint \f$a\cdot x \lessgtr b\f$. */
      LinearConstraint(const LinearAlgebra::Vector<R> a, Comparison cmp, const R& b);

      /*! \brief Destructor. */
      virtual ~LinearConstraint();
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual LinearConstraint<R>* clone() const;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const;
      /*! \brief The smoothness of the constraint function. */
      virtual smoothness_type smoothness() const;
      /*! \brief The operation used for comparison. */
      virtual Comparison comparison() const;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;

      /*! \brief The value at a point. */
      virtual A value(const Point<A>& pt) const;
      /*! \brief The gradient at a point. */
      virtual LinearAlgebra::Vector<A> gradient(const Point<A>& pt) const;

      /*! \brief Convert to a polyhedron. */
      virtual Polyhedron<R> polyhedron() const;
     private:
      static void instantiate();
     private:
      LinearAlgebra::Vector<R> _a;
      R _b;
      Comparison _c;
    };


  }
}

#endif /* ARIADNE_LINEAR_CONSTRAINT_H */
