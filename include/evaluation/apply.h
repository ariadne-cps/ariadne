/***************************************************************************
 *            apply.h
 *
 *  17 January 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file apply.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef _ARIADNE_APPLY_H
#define _ARIADNE_APPLY_H

#include "../geometry/geometry_declarations.h"

#include "../evaluation/evaluation_declarations.h"
#include "../evaluation/map.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Compute the image of a rectangle under a continuous function. */
    template<typename R>
    Geometry::Rectangle<R> 
    apply(const Map<R>& f, const Geometry::Rectangle<R>& p);

    /*! \brief Compute the image of a parallelotope under a differentiable function. */
    template<typename R>
    Geometry::Parallelotope<R> 
    apply(const Map<R>& f, const Geometry::Parallelotope<R>& p);

    /*! \brief Compute the image of a list set under a map. */
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    apply(const Map<R>& f, const Geometry::ListSet<R,BS>& ds);
     
    /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set on the grid \a grid while staying within \a bounds. */
    template<typename R>
    Geometry::GridMaskSet<R> 
    chainreach(const Map<R>& map, 
               const Geometry::ListSet<R,Geometry::Rectangle>& initial_set, 
               const Geometry::FiniteGrid<R>& grid, 
               const Geometry::Rectangle<R>& bounds);
    
  }
}

#endif /* _ARIADNE_APPLY_H */
