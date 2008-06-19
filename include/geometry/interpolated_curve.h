/***************************************************************************
 *            interpolated_curve.h
 *
 *  Copyright  2008  Pieter Collins
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

/*! \file interpolated_curve.h
 *  \brief Interpolated curves in Euclidean space.
 */

#ifndef ARIADNE_INTERPOLATED_CURVE_H
#define ARIADNE_INTERPOLATED_CURVE_H

#include <iosfwd>
#include <map>

#include "geometry/point.h"
#include "geometry/segment.h"
#include "exceptions.h"

namespace Ariadne {
     
    /*!\brief A line segment in Euclidean space. */
    template<class R, class X=R>
    class InterpolatedCurve {
     public:
      typedef typename std::map< R, Point<X> >::const_iterator const_iterator;

     public:
      /*! \brief Create a curve with a single point \a pt at parameter value 0. */
      template<class XX> InterpolatedCurve(const Point<XX>& pt) 
        : _points() { this->insert(0,pt); }
      /*! \brief Create a curve with a single point \a pt at parameter value \a s. */
      template<class RR, class XX> InterpolatedCurve(const RR& s, const Point<XX>& pt) 
        : _points() { this->insert(s,pt); }
      /*! \brief Create a segment from \a pt0 at parameter value 0 to \a pt1 at parameter value 1. */
      template<class XX0, class XX1> InterpolatedCurve(const Point<XX0>& pt0, const Point<XX1>& pt1) 
        : _points() { this->insert(0,pt0); this->insert(1,pt1); }
      /*! \brief Create a segment from \a pt0 at parameter value 0 to \a pt1 at parameter value 1. */
      template<class XX> InterpolatedCurve(const Segment<XX>& seg) 
        : _points() { this->insert(0,seg.initial_point()); this->insert(1,seg.final_point()); }
      /*! \brief Insert a point with parameter value \a s and spacial value \a pt. */
      template<class RR, class XX> void insert(const RR& s, const Point<XX>& pt) {
        if(!this->_points.empty()) { ARIADNE_ASSERT(pt.dimension()==this->dimension()); }
        this->_points.insert(std::pair< R, Point<X> >(s,pt)); }
       
      /*! \brief The number of segments in the curve. */
      size_type size() const { return this->_points.size(); }
      /*! \brief The dimension of the Euclidean space the line segment lies in. */
      dimension_type dimension() const { return this->_points.begin()->second.dimension(); }
      /*! \brief An iterator to the first point in the curve. */
      const_iterator begin() const { return this->_points.begin(); }
      /*! \brief An iterator to the end point in the curve, NOT the one-past-the-end! */
      const_iterator end() const { return --this->_points.end(); }
     private:
      std::map< X, Point<X> > _points;
    };

    template<class R, class X> 
    std::ostream& operator<<(std::ostream& os, const InterpolatedCurve<R,X>& curve) {
      typename InterpolatedCurve<R,X>::const_iterator iter=curve.begin();
      os << "InterpolatedCurve( size=" << curve.size() << ", [ " << iter->first << ":" << iter->second; 
      while(iter!=curve.end()) { ++iter; os << ", " << iter->first << ":" << iter->second; }
      return os << " ] )"; 
    }

} // namespace Ariadne

#endif /* ARIADNE_INTERPOLATED_CURVE_H */
