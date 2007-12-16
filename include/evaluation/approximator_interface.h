/***************************************************************************
 *            approximator_interface.h
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
 
/*! \file approximator_interface.h
 *  \brief Interface for approximating basic sets
 */

#ifndef ARIADNE_APPROXIMATOR_INTERFACE_H
#define ARIADNE_APPROXIMATOR_INTERFACE_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    template<class BS> 
    class ApproximatorInterface 
    { 
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Virtual destructor. */
      virtual ~ApproximatorInterface() { };

      /*! \brief Make a dynamically-allocated copy. */
      virtual ApproximatorInterface<BS>* clone() const = 0;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      virtual BS over_approximation(const Geometry::Rectangle<R>& r) const = 0;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      virtual Geometry::GridCellListSet<R> outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const = 0;

      /*! \brief Adjoins an outer approximation to a basic set to a grid mask set. */
      void adjoin_outer_approximation(Geometry::GridMaskSet<R>& gms, const BS& bs) const;

      /*! \brief Computets and over-approximation of a set from a rectangle. */
      void adjoin_outer_approximation(Geometry::GridMaskSet<R>& gms, const Geometry::ListSet<BS>& ls) const;

    };


  }
}

#include "approximator_interface.inline.h"

#endif /* ARIADNE_APPROXIMATOR_H */
