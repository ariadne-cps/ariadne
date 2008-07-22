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

#include "geometry/grid_approximation_scheme.h"
#include "approximator_base.h"

namespace Ariadne {
  


    /*! \ingroup Approximators
     *  \brief Fast geomerical approximation scheme.
     *  
     *  Faster but less accurate approximation methods.
     */
    template<class ES>
    class FastApproximator
      : public ApproximatorBase<GridApproximationScheme<typename ES::real_type>,ES>
    {
      typedef typename ES::real_type R;
      typedef Interval<R> I;
      typedef GridApproximationScheme<R> GAS;
     public:
      FastApproximator() : ApproximatorBase<GAS,ES>() { };
      virtual FastApproximator<ES>* clone() const;
      virtual ES enclosure_set(const Box<R>&  bx) const;
      virtual R radius(const ES& bs) const;
      virtual tribool disjoint(const ES& es, const Box<R>& bx) const;
      virtual tribool superset(const ES& es, const Box<R>& bx) const;
      virtual Box<R> bounding_box(const ES& es) const;
      virtual BoxListSet<R> over_approximation(const ES& es) const;
      virtual void adjoin_outer_approximation(GridCellListSet<R>& gcls, const ES& bs) const;
      virtual void adjoin_outer_approximation(GridMaskSet<R>& gms, const ES& bs) const;
      virtual std::ostream& write(std::ostream& os) const;
    };

  
} // namespace Ariadne

#endif /* ARIADNE_FAST_APPROXIMATOR_H */
