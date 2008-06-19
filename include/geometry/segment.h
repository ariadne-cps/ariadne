/***************************************************************************
 *            segment.h
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
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

/*! \file segment.h
 *  \brief Line segments in Euclidean space.
 */

#ifndef ARIADNE_SEGMENT_H
#define ARIADNE_SEGMENT_H

#include <iosfwd>
#include <stdexcept>

#include "geometry/euclidean_space.h"
#include "geometry/point.h"
#include "exceptions.h"

namespace Ariadne {
     
    /*!\brief A line segment in Euclidean space. */
    template<class X>
    class Segment {
     public:
      /*! \brief The line segment from \a pt0 to \a pt1. */
      template<class X1, class X2> Segment(const Point<X1>& pt0, const Point<X2>& pt1) 
        : _initial_point(pt0), _final_point(pt1) { ARIADNE_ASSERT(pt0.dimension()==pt1.dimension()); }
       
      /*! \brief The dimension of the Euclidean space the line segment lies in. */
      dimension_type dimension() const { return this->_initial_point.dimension(); }
      /*! \brief The initial point of the line segment. */
      const Point<X>& initial_point() const { return this->_initial_point; }
      /*! \brief The final point of the line segment. */
      const Point<X>& final_point() const { return this->_final_point; }
     private:
      Point<X> _initial_point;
      Point<X> _final_point;
    };

    template<class X> 
    std::ostream& operator<<(std::ostream& os, const Segment<X>& seg) {
      return os << "Segment( initial_point=" << seg.initial_point() << ", final_point=" << seg.final_point() << " )";
    }

} // namespace Ariadne

#endif /* ARIADNE_SEGMENT_H */
