/***************************************************************************
 *            satisfier_interface.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file satisfier_interface.h
 *  \brief Interface for computing the image of a basic set under a map.
 */

#ifndef ARIADNE_SATISFIER_INTERFACE_H
#define ARIADNE_SATISFIER_INTERFACE_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a basic set under a map. 
     *  \ingroup Satisfiers
     */
    template<class BS>
    class SatisfierInterface
    {
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ~SatisfierInterface() { }

      /*! \brief Make a dynamically-allocated copy. */
      virtual SatisfierInterface<BS>* clone() const = 0;
      
      /*! \brief Test whether a set is a subset of a given constraint set. */
      virtual tribool subset(const BS& bs, const Geometry::ConstraintSet<R>& cs) const = 0;
      
      /*! \brief Test whether a set is disjoint from a given constraint set. */
      virtual tribool disjoint(const BS& bs, const Geometry::ConstraintSet<R>& cs) const = 0;
      
      /*! \brief Test whether a set intersects from a given constraint set. */
      tribool intersects(const BS& bs, const Geometry::ConstraintSet<R>& cs) const {
        return !this->disjoint(bs,cs); }
      
    };


  }
}

#endif /* ARIADNE_SATISFIER_INTERFACE_H */
