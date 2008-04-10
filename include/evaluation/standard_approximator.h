/***************************************************************************
 *            standard_approximator.h
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
 
/*! \file standard_approximator.h
 *  \brief Methods for approximating basic sets
 */

#ifndef ARIADNE_STANDARD_APPROXIMATOR_H
#define ARIADNE_STANDARD_APPROXIMATOR_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

#include "geometry/grid_approximation_scheme.h"
#include "approximator_base.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Geomerical approximation schemes.
     *  \ingroup Approximators
     */
    template<class ES>
    class StandardApproximator
      : public ApproximatorBase< Geometry::GridApproximationScheme<typename ES::real_type>,ES>
    {
      typedef typename ES::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      StandardApproximator() { }
      virtual StandardApproximator<ES>* clone() const { return new StandardApproximator<ES>(*this); }
      virtual ES enclosure_set(const Geometry::Box<R>& r) const;
      virtual R radius(const ES& bs) const;
      virtual Geometry::Box<R> bounding_box(const ES& bs) const;
      virtual Geometry::BoxListSet<R> lower_approximation(const ES& es) const;
      virtual Geometry::GridCellListSet<R> inner_approximation(const ES& es, const Geometry::Grid<R>& g) const;
      virtual Geometry::GridCellListSet<R> outer_approximation(const ES& es, const Geometry::Grid<R>& g) const;
    };



  }
}

#endif /* ARIADNE_STANDARD_APPROXIMATOR_H */
