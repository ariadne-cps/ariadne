/***************************************************************************
 *            set_constraint.h
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
 
/*! \file set_constraint.h
 *  \brief A constraint defined by a set.
  */

#ifndef ARIADNE_SET_CONSTRAINT_H
#define ARIADNE_SET_CONSTRAINT_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "../base/tribool.h"

#include "../geometry/set_interface.h"
#include "constraint_interface.h"

namespace Ariadne {
  namespace Geometry {
    

    //! \ingroup ExactSet
    /*! \brief A constraint on the state
     */
    template<class R>
    class SetConstraint
      : public ConstraintInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Construct the set \f$f(x) \lessgtr 0\f$ from the function \f$f\f$. */
      SetConstraint(const Geometry::SetInterface<R>& s, bool i=true);

      /*! \brief Destructor. */
      virtual ~SetConstraint();
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual SetConstraint<R>* clone() const;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const;
      /*! \brief The smoothness of the constraint function. */
      virtual smoothness_type smoothness() const;
      /*! \brief The operation used for comparison. */
      virtual Comparison comparison() const;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;

      /*! \brief The value of the constraint function at a point. */
      virtual A value(const Point<A>& pt) const;
      /*! \brief The gradient of the constraint function at a point. */
      virtual LinearAlgebra::Vector<A> gradient(const Point<A>& pt) const;

      /*! \brief The function defining the constraint. */
      const Geometry::SetInterface<R>& set() const;
     private:
      static void instantiate();
     private:
      boost::shared_ptr< const Geometry::SetInterface<R> > _set_ptr;
      bool _inside;
    };

  }
}

#endif /* ARIADNE_SET_CONSTRAINT_H */
