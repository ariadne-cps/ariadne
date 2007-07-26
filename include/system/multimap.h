/***************************************************************************
 *            multimap.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it
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
 
/*! \file multimap.h
 *  \brief Mulitvalued map interface.
 */

#ifndef ARIADNE_MULTIMAP_H
#define ARIADNE_MULTIMAP_H

#include <string>

#include "../base/types.h"

namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Abstract base class for multivalued functions.
     * 
     * A multivalued function is specified by operator()(const Geometry::Point<R>& pt) const.
     * This must return a list set \f$\hat{f}(A)\f$ such that either
     *   - <em>lower-semicontinuity:</em> all sets of \f$\hat{f}(A)\f$ intersect \f$f(A)\f$, or
     *   - <em>upper-semicontinuity:</em> the true image \f$\hat{f}(A)\f$ is a subset of \f$\overline{f}(A)\f$.
     *
     * Further, the image must converge monotonically as the set \a A converges to a point.
     *
     * Derivatives of multivalued maps are not supported.
     */
     template<class P, class S>
    class MultiMap {
     public:
      /*! \brief The real number type. */
      typedef typename P::real_type real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef P state_type;
      /*! \brief The type of set the system value takes at every point. */
      typedef S set_type;

      /*! \brief Virtual destructor. */
      virtual ~MultiMap();
      
      /*! \brief Call the function at a point. */
      S operator() (const P& x) const { return this->image(x); }
    
      /*! \brief The image of the map on a point. */
      virtual S image(const P& x) const = 0;
      
      /*! \brief The dimension of the domain space. */
      virtual dimension_type argument_dimension() const = 0;
      /*! \brief The dimension of the range space. */
      virtual dimension_type result_dimension() const = 0;
    
      /*! \brief The name of the map. */
      virtual std::string name() const = 0;
    };
  
    
  }
}

#endif /* ARIADNE_MULTIMAP_H */
