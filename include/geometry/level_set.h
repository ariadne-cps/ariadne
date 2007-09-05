/***************************************************************************
 *            level_set.h
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
 
/*! \file level_set.h
 *  \brief Sets defined by function and constraints.
  */

#ifndef ARIADNE_LEVEL_SET_H
#define ARIADNE_LEVEL_SET_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>

#include "../base/tribool.h"

#include "../function/function_interface.h"
#include "../geometry/set_interface.h"

namespace Ariadne {
  namespace Geometry {
    


    //! \ingroup ExactSet
    /*! \brief A set defined by the conditions \f$f(x)=0\f$ for some function \f$f\f$. 
     */
    template<class R>
    class LevelSet
      : public SetInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
     public:
      /*! \brief Construct the set \f$f(x)\geq0\f$ from the function \f$f\f$. */
      LevelSet(const Function::FunctionInterface<R>& f);

      /*! \brief Destructor. */
      virtual ~LevelSet();
      /*! \brief Return a new dynamically-allocated copy of the set. */
      virtual LevelSet<R>* clone() const;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const;
      /*! \brief */
      virtual tribool contains(const Point<R>& pt) const;
      /*! \brief */
      virtual tribool superset(const Rectangle<R>& r) const;
      /*! \brief */
      virtual tribool intersects(const Rectangle<R>& r) const;
      /*! \brief */
      virtual tribool disjoint(const Rectangle<R>& r) const;
      /*! \brief */
      virtual tribool subset(const Rectangle<R>& r) const;
      /*! \brief */
      virtual tribool bounded() const;
      /*! \brief */
      virtual Rectangle<R> bounding_box() const;
      /*! \brief */
      virtual std::ostream& write(std::ostream& os) const;

      /*! \brief The number of independed inequality constraints used to define the set. */
      size_type number_of_constraints() const;
      /*! \brief The countour function defining the level set. */
      const Function::FunctionInterface<R>& function() const;
      /*! \brief The countour function applied to a point. */
      Point<A> function(const Point<A>& pt) const;

      /*! \brief Test if the two points lie on opposite components of the constaint set. */
      tribool separates(const Point<A>& pt1, const Point<A>& pt2) const;
     private:
      boost::shared_ptr< const Function::FunctionInterface<R> > _function_ptr;
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const LevelSet<R>& lset) {
      return lset.write(os);
    }


  }
}

#endif /* ARIADNE_LEVEL_SET_H */
