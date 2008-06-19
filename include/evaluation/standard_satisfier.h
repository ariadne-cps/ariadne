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

#ifndef ARIADNE_STANDARD_SATISFIER_H
#define ARIADNE_STANDARD_SATISFIER_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "satisfier_interface.h"

namespace Ariadne {
  

    /*! \ingroup EvaluatorInterfaces \ingroup Approximators
     *  \brief Interface for testing whether a set satisfies a constraint. 
     */
    template<class ES>
    class StandardSatisfier
      : public SatisfierInterface<ES>
    {
      typedef typename ES::real_type R;
     public:
      /*! \brief Make a dynamically-allocated copy. */
      virtual StandardSatisfier<ES>* clone() const {
        return new StandardSatisfier<ES>(*this); }
      
      /*! \brief Test whether a set is a subset of a given constraint set. */
      virtual tribool subset(const ES& es, const ConstraintSet<R>& cs) const {
        return subset(es,cs); }
      
      /*! \brief Test whether a set is disjoint from a given constraint set. */
      virtual tribool disjoint(const ES& es, const ConstraintSet<R>& cs) const {
        return disjoint(es,cs); }
      
      /*! \brief Test whether a set intersects from a given constraint set. */
      virtual tribool intersects(const ES& es, const ConstraintSet<R>& cs) const {
        return intersects(es,cs); }
      
    };


  
} // namespace Ariadne

#endif /* ARIADNE_STANDARD_SATISFIER_H */
