/***************************************************************************
 *            standard_subdivider.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file standard_subdivider.h
 *  \brief Methods for subdividing basic sets
 */

#ifndef ARIADNE_STANDARD_SUBDIVIDER_H
#define ARIADNE_STANDARD_SUBDIVIDER_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

#include "subdivider_interface.h"
#include "geometry/zonotope.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Method for subdividing a basic set.
     *  \ingroup Approximation
     */
    template<class BS>
    class StandardSubdivider : public SubdividerInterface<BS> {
     public:
      virtual StandardSubdivider<BS>* clone() const {
        return new StandardSubdivider<BS>(*this); }
      virtual Geometry::ListSet<BS> split(const BS& bs) const {
        return Geometry::split(bs); }
    };

    template<class R>
    class OrthogonalSubdivider
      : public SubdividerInterface< Geometry::Zonotope<R> > {
      typedef Geometry::Zonotope<R> BS;
     public:
      virtual StandardSubdivider< BS >* clone() const {
        return new StandardSubdivider< BS >(*this); }
      virtual BS reduce(const BS& bs) const {
        return Geometry::orthogonal_over_approximation(bs); }
      virtual Geometry::ListSet<BS> split(const BS& bs) const {
        return Geometry::split(bs); }
    };



  }
}

#endif /* ARIADNE_STANDARD_SUBDIVIDER_H */
