/****************************************************************************
 *            polytope.h
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_POLYTOPE_H
#define ARIADNE_POLYTOPE_H

/*! \file polytope.h
 *  \brief Polytope class for geometry output.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "utility/macros.h"
#include "geometry/point.h"

#include "geometry/set_interface.h"
#include "output/graphics_interface.h"

namespace Ariadne {

class Zonotope;
class Polytope;
class InterpolatedCurve;

Polytope polytope(const ExactBoxType& bx);
Polytope polytope(const Zonotope& z);
Polytope polytope(const Polytope& p);

ExactPoint baricentre(const Polytope& p);

Polytope& reduce2d(Polytope& p);
Float64 slope2d(const ExactPoint& pt1, const ExactPoint& pt2);


//! A polytope in Euclidean space, defined by a list of extreme points.
class Polytope
    : public LocatedSetInterface
    , public DrawableInterface
{
  public:
    typedef std::vector<ExactPoint>::const_iterator ConstIterator;

    //! \brief Default constructor constructs an empty polytope in zero dimensions.
    Polytope() { }
    //! \brief Construct an empty polytope in \a d dimensions.
    Polytope(Nat d) { }
    //! \brief Construct polytope with vertices in \a v.
    Polytope(const std::vector<ExactPoint>& v) : _vertices(v) { }

    //! \brief Add a new vertex \a v to the current vertices.
    Void new_vertex(const ExactPoint& v) {
        ARIADNE_ASSERT(this->_vertices.size()==0 || v.dimension()==this->_vertices.front().dimension());
        this->_vertices.push_back(v); }
    //! \brief The number of points defining the polytope. Note that an interior point is counted by this method.
    SizeType number_of_vertices() const { return this->_vertices.size(); }
    //! \brief The \a i<sup>th</sup> vertex.
    const ExactPoint& vertex(SizeType i) const { return this->_vertices[i]; }
    const std::vector<ExactPoint>& vertices() const { return this->_vertices; }
    std::vector<ExactPoint>& vertices() { return this->_vertices; }
    //! \brief A constant Iterator pointing to the first vertex.
    ConstIterator vertices_begin() const { return this->_vertices.begin(); }
    //! \brief A constant Iterator pointing to the past-the-end vertex.
    ConstIterator vertices_end() const { return this->_vertices.end(); }
    //! \brief Reduce the description of the polytope by removing all non-extreme points from the list of vertices. (Not currently implemented)
    Void reduce() { ARIADNE_NOT_IMPLEMENTED; }

    virtual Polytope* clone() const { return new Polytope(*this); }
    virtual DimensionType dimension() const { if(this->_vertices.size()==0) { return 0; } return this->_vertices.front().dimension(); }
    virtual ValidatedKleenean separated(const ExactBoxType& bx) const;
    virtual ValidatedKleenean overlaps(const ExactBoxType& bx) const;
    virtual ValidatedKleenean inside(const ExactBoxType& bx) const;
    virtual UpperBoxType bounding_box() const;
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;
    virtual OutputStream& write(OutputStream& os) const;

    //! \brief The convex hull of two polytopes.
    friend Polytope convex_hull(const Polytope& p1, const Polytope& p2);
    //! \brief The baricentre of a polytope.
    friend ExactPoint baricentre(const Polytope& p);
  private:
    std::vector<ExactPoint> _vertices;
};
inline OutputStream& operator<<(OutputStream& os, const Polytope& p) { return p.write(os); }




} // namespace Ariadne


#endif /* ARIADNE_POLYTOPE_H */
