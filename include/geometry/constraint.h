/***************************************************************************
 *            constraint.h
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
 
/*! \file constraint.h
 *  \brief A constraint defining an invariant, activation or guard set.
  */

#ifndef ARIADNE_CONSTRAINT_H
#define ARIADNE_CONSTRAINT_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "base/tribool.h"

#include "function/function_interface.h"
#include "constraint_interface.h"

namespace Ariadne {
  namespace Geometry {
    
    template<class R> class Point;
    template<class R> class Box;
    template<class R0,class R1> class Zonotope;
    template<class R> class Polyhedron;

  
  /*
    // Forward declarations for friends
    template<class R> class DifferentiableConstraint;
    template<class R> bool equal(const DifferentiableConstraint<R>& c1, const DifferentiableConstraint<R>& c2);
    template<class R> bool opposite(const DifferentiableConstraint<R>& c1, const DifferentiableConstraint<R>& c2);
  */  


    //! \ingroup ExactSet
    /*! \brief A constraint on the state
     */
    template<class R>
    class DifferentiableConstraint
      : public DifferentiableConstraintInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Construct the set \f$f(x) \lessgtr 0\f$ from the function \f$f\f$. */
      DifferentiableConstraint(const Function::DifferentiableFunctionInterface<R>& f, const Comparison cmp=greater);

      /*! \brief Destructor. */
      virtual ~DifferentiableConstraint();
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual DifferentiableConstraint<R>* clone() const;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const;
      /*! \brief The smoothness of the constraint function. */
      virtual smoothness_type smoothness() const;
      /*! \brief The operation used for comparison. */
      virtual Comparison comparison() const;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;

      /*! \brief The function defining the constraint. */
      const Function::DifferentiableFunctionInterface<R>& function() const;
      /*! \brief The value at a point. */
      A value(const Point<A>& pt) const;
      /*! \brief The gradient at a point. */
      LinearAlgebra::Vector<A> gradient(const Point<A>& pt) const;
     public:
      /*! \brief Test for equality as reference. */
      //      friend bool equal<>(const Constraint<R>& c1, const Constraint<R>& c2);
      /*! \brief Test for equality as reference, but with different sign. */
      //      friend bool opposite<>(const Constraint<R>& c1, const Constraint<R>& c2);
#ifdef DOXYGEN
      /*! \brief Test if the constraint is satisfied over a rectangle. */
      friend tribool satisfies(const Box<R>& r, const Constraint<R>& c);
      /*! \brief Test if the constraint is satisfied over a zonotope. */
      friend tribool satisfies(const Zonotope<R,R>& z), const Constraint<R>& c);
      /*! \brief Test if the constraint is satisfied over a zonotope. */
      friend tribool satisfies(const Zonotope<I,R>& z), const Constraint<R>& c);
      /*! \brief Test if the constraint is satisfied over a zonotope. */
      friend tribool satisfies(const Zonotope<I,I>& z), const Constraint<R>& c);
#endif
     private:
      static void instantiate();
     private:
      boost::shared_ptr< const Function::DifferentiableFunctionInterface<R> > _function_ptr;
      Comparison _comparison;
    };
    
/*
    template<class R> bool equal(const DifferentiableConstraint<R>& c1, const DifferentiableConstraint<R>& c2);
    template<class R> bool opposite(const DifferentiableConstraint<R>& c1, const DifferentiableConstraint<R>& c2);

    template<class R> tribool satisfies(const Box<R>& r, const DifferentiableConstraint<R>& c);
    template<class R> tribool satisfies(const Zonotope<R,R>& z, const DifferentiableConstraint<R>& c);
    template<class R> tribool satisfies(const Zonotope<Numeric::Interval<R>,R>& z, const DifferentiableConstraint<R>& c);
    template<class R> tribool satisfies(const Zonotope< Numeric::Interval<R>,Numeric::Interval<R> >& z, const DifferentiableConstraint<R>& c);
*/
  }
}

#include "constraint.inline.h"

#endif /* ARIADNE_CONSTRAINT_H */
