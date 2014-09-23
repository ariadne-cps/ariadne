/***************************************************************************
 *            grid.h
 *
 *  Copyright  2008-9  Ivan S. Zapreev, Pieter Collins
 *            ivan.zapreev@gmail.com, pieter.collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file grid.h
 *  \brief Coordinate-aligned grids.
 */

#ifndef ARIADNE_GRID_H
#define ARIADNE_GRID_H

#include "array.h"
#include "numeric.h"
#include "vector.h"

#include "point.h"
#include "box.h"

namespace Ariadne {

class Grid;
class GridCell;

/*! \brief An infinite, uniform grid of rectangles in Euclidean space.
 *
 *  \internal Maybe a Grid should be a type of Paving or Cover.
 *  Then rather than having GridXXX classes, we can have classes such that cells of
 *  some type are mapped into concrete sets by a Paving or Cover.
 *  This should be more general, and will unify the concepts of Paving and Cover,
 *  as well as different types of covers.
 */
class Grid {
    typedef double DyadicType;
    typedef int IntegerType;
  private:
    // Structure containing actual data values
    struct Data;
  public:
    //! Destructor.
    ~Grid();

    //! Default constructor constructs a grid from a null pointer. Needed for some iterators.
    explicit Grid();

    //! Construct from a dimension and a spacing in each direction.
    explicit Grid(uint d);

    //! Construct from a dimension and a spacing in each direction.
    explicit Grid(uint d, RawFloatType l);

    //! Construct from a vector of offsets.
    explicit Grid(const Vector<RawFloatType>& lengths);

    //! Construct from a centre point and a vector of offsets.
    explicit Grid(const Vector<RawFloatType>& origin, const Vector<RawFloatType>& lengths);

    //! Copy constructor. Copies a reference to the grid data.
    Grid(const Grid& g);

    //! The underlying dimension of the grid.
    uint dimension() const;

    //! Tests equality of two grids. Tests equality of references first.
    bool operator==(const Grid& g) const;

    //! Tests inequality of two grids.
    bool operator!=(const Grid& g) const;

    //! The origin of the grid.
    const Vector<RawFloatType>& origin() const;

    //! The strides between successive integer points.
    const Vector<RawFloatType>& lengths() const;

    //! Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Grid& g);

    ExactNumberType coordinate(uint d, DyadicType x) const;
    ExactNumberType subdivision_coordinate(uint d, DyadicType x) const;
    ExactNumberType subdivision_coordinate(uint d, IntegerType n) const;

    int subdivision_index(uint d, const ExactNumberType& x) const;
    int subdivision_lower_index(uint d, const LowerNumberType& x) const;
    int subdivision_upper_index(uint d, const UpperNumberType& x) const;

    Array<DyadicType> index(const Point& pt) const;
    Array<DyadicType> lower_index(const Box& bx) const;
    Array<DyadicType> upper_index(const Box& bx) const;

    Point point(const Array<IntegerType>& a) const;
    Point point(const Array<DyadicType>& a) const;
    Box box(const Array<DyadicType>& l, const Array<DyadicType>& u) const;
    Box box(const GridCell& cell) const;
  private:
    // Create new data
    void _create(const Vector<RawFloatType>& o, const Vector<RawFloatType>& l);
  private:
    // Pointer to data. We can test grids for equality using reference semantics since data is a constant.
    std::shared_ptr<Data> _data;
};


} // namespace Ariadne

#endif /* ARIADNE_GRID_H */

