/***************************************************************************
 *            constraint_set.h
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
 
/*! \file constraint_set.h
 *  \brief Sets defined by function and constraints.
  */

#ifndef ARIADNE_CONSTRAINT_SET_H
#define ARIADNE_CONSTRAINT_SET_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "base/tribool.h"

#include "function/function_interface.h"
#include "geometry/set_interface.h"

namespace Ariadne {
  namespace Geometry {
    
    //! \ingroup ExactSet
    /*! \brief A set defined by the conditions \f$f(x)\geq0\f$ for some function \f$f\f$. 
     *   Satisfies the conditions of the RegularSetInterface, 
     *   which means that only superset(Box) and disjoint(Box) need be meaningfully defined, 
     *   and only outer- and inner-approximations can be computed.
     */
    template<class R>
    class ConstraintSet
      : public SetInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
     public:
      /*! \brief Construct the set \f$x\in D\f$ from the function \f$f\f$ and the box \a D. */
      ConstraintSet(const Geometry::Box<R>& D);
      /*! \brief Construct the set \f$f(x)\in D\f$ from the function \f$f\f$ and the box \a D. */
      ConstraintSet(const Function::FunctionInterface<R>& f, const Geometry::Box<R>& D);

      /*! \brief Destructor. */
      virtual ~ConstraintSet();
      /*! \brief Return a new dynamically-allocated copy of the set. */
      virtual ConstraintSet<R>* clone() const;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const;
      /*! \brief Test if the set contains a point. */
      virtual tribool contains(const Point<R>& pt) const;
      /*! \brief */
      virtual tribool superset(const Box<R>& r) const;
      /*! \brief */
      virtual tribool intersects(const Box<R>& r) const;
      /*! \brief */
      virtual tribool disjoint(const Box<R>& r) const;
      /*! \brief */
      virtual tribool subset(const Box<R>& r) const;
      /*! \brief */
      virtual tribool bounded() const;
      /*! \brief */
      virtual Box<R> bounding_box() const;
      /*! \brief */
      virtual std::ostream& write(std::ostream& os) const;

      /*! \brief The number of independed inequality constraints used to define the set. */
      size_type number_of_constraints() const;
      /*! \brief The codomain given the allowable values of the constraint function. */
      const Geometry::Box<R>& codomain() const;
      /*! \brief The function describing the constraints. */
      const Function::FunctionInterface<R>& function() const;

      /*! \brief Test if the set contains a point described approximately. */
      tribool contains(const Point<A>& pt) const;
     private:
      boost::shared_ptr< const Function::FunctionInterface<R> > _function_ptr;
      Geometry::Box<R> _codomain;
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const ConstraintSet<R>& cset) {
      return cset.write(os);
    }
    




  }
}

#endif /* ARIADNE_CONSTRAINT_SET_H */
