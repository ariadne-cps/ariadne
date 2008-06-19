/***************************************************************************
 *            subdivider_interface.h
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
 
/*! \file subdivider_interface.h
 *  \brief Interface for approximating basic sets
 */

#ifndef ARIADNE_SUBDIVIDER_INTERFACE_H
#define ARIADNE_SUBDIVIDER_INTERFACE_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {
  

    /*! \brief Interface for methods subdividing enclosure sets into smaller pieces.
     *  \ingroup EvaluatorInterfaces \ingroup Approximators
     */
    template<class ES> 
    class SubdividerInterface 
    { 
      typedef typename ES::real_type R;
      typedef ListSet<ES> ESL;
     public:
      /*! \brief Virtual destructor. */
      virtual ~SubdividerInterface() { }

      /*! \brief Create a dynamically-allocated copy. */
      virtual SubdividerInterface<ES>* clone() const = 0;

      /*! \brief Computes the radius of a basic set. */
      virtual R radius(const ES& es) const = 0;

      /*! \brief Subdivide the set \a bs into two smaller pieces. */
      virtual ESL split(const ES& es) const = 0;

      /*! \brief Subdivide the set \a bs into smaller pieces, each with radius at most \a r. */
      virtual ESL subdivide(const ES& es, const R& r) const = 0;
    };

  
} // namespace Ariadne

#endif /* ARIADNE_SUBDIVIDER_INTERFACE_H */
