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
  namespace Evaluation {

    template<class BS> 
    class SubdividerInterface 
    { 
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
      typedef Geometry::ListSet<BS> BSL;
     public:
      /*! \brief Virtual destructor. */
      virtual ~SubdividerInterface() { }

      /*! \brief Create a dynamically-allocated copy. */
      virtual SubdividerInterface<BS>* clone() const = 0;

      /*! \brief Returns an over-approximation which simplifies the description of the set. */
      virtual BS reduce(const BS& bs) const;

      /*! \brief Subdivide the set \a bs into two smaller pieces. */
      virtual Geometry::ListSet<BS> split(const BS& bs) const = 0;

      /*! \brief Subdivide the set \a bs into smaller pieces, each with radius at most \a r. */
      virtual Geometry::ListSet<BS> subdivide(const BS& bs, const R& r) const;
    };

  }
}



namespace Ariadne {

template<class BS> inline
BS 
Evaluation::SubdividerInterface<BS>::reduce(const BS& bs) const 
{
  return bs;
}

template<class BS>
Geometry::ListSet<BS> 
Evaluation::SubdividerInterface<BS>::subdivide(const BS& bs, const R& r) const
{
  BSL result; 
  BS set=this->reduce(bs); 
  BSL working(set);
  while(working.size()!=0) { 
    set=working.pop(); 
    if(set.radius()<r) {
      result.adjoin(set);
    } else {
      working.adjoin(this->split(set)); 
    }
  }
  return result;
}


}

#endif /* ARIADNE_SUBDIVIDER_INTERFACE_H */
