/***************************************************************************
 *            fast_approximator.h
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
 
/*! \file approximator.h
 *  \brief Methods for approximating basic sets
 */

#ifndef ARIADNE_FAST_APPROXIMATOR_H
#define ARIADNE_FAST_APPROXIMATOR_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

#include "approximator_interface.h"

namespace Ariadne {
  namespace Evaluation {


    /*! \brief Geomerical approximation schemes.
     *  \ingroup Faster but less accurate approximation methods.
     */
    template<class BS>
    class FastApproximator
      : public ApproximatorInterface<BS>
    {
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      FastApproximator();
      FastApproximator(const FastApproximator<BS>& approx);
      virtual FastApproximator<BS>* clone() const;
      virtual BS basic_set(const Geometry::Box<R>&  bx) const;
      virtual R radius(const BS& bs) const;
      virtual Geometry::Box<R> bounding_box(const BS& bs) const;
      virtual Geometry::GridCellListSet<R> outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const;
    };

  }
}

#endif /* ARIADNE_FAST_APPROXIMATOR_H */
