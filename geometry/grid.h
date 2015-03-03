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

#include "utility/array.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"

#include "geometry/point.h"
#include "geometry/box.h"

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
    typedef Int IntegerType;
  private:
    // Structure containing actual data values
    struct Data;
  public:
    //! Destructor.
    ~Grid();

    //! Default constructor constructs a grid from a null pointer. Needed for some iterators.
    explicit Grid();

    //! Construct from a dimension and a spacing in each direction.
    explicit Grid(Nat d);

    //! Construct from a dimension and a spacing in each direction.
    explicit Grid(Nat d, RawFloat64 l);

    //! Construct from a vector of offsets.
    explicit Grid(const Vector<RawFloat64>& lengths);

    //! Construct from a centre point and a vector of offsets.
    explicit Grid(const Vector<RawFloat64>& origin, const Vector<RawFloat64>& lengths);

    //! Copy constructor. Copies a reference to the grid data.
    Grid(const Grid& g);

    //! The underlying dimension of the grid.
    DimensionType dimension() const;

    //! Tests equality of two grids. Tests equality of references first.
    Bool operator==(const Grid& g) const;

    //! Tests inequality of two grids.
    Bool operator!=(const Grid& g) const;

    //! The origin of the grid.
    const Vector<RawFloat64>& origin() const;

    //! The strides between successive integer points.
    const Vector<RawFloat64>& lengths() const;

    //! Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const Grid& g);

    ExactNumericType coordinate(Nat d, DyadicType x) const;
    ExactNumericType subdivision_coordinate(Nat d, DyadicType x) const;
    ExactNumericType subdivision_coordinate(Nat d, IntegerType n) const;

    Int subdivision_index(Nat d, const ExactNumericType& x) const;
    Int subdivision_lower_index(Nat d, const LowerNumericType& x) const;
    Int subdivision_upper_index(Nat d, const UpperNumericType& x) const;

    Array<DyadicType> index(const ExactPoint& pt) const;
    Array<DyadicType> lower_index(const ExactBox& bx) const;
    Array<DyadicType> upper_index(const ExactBox& bx) const;

    ExactPoint point(const Array<IntegerType>& a) const;
    ExactPoint point(const Array<DyadicType>& a) const;
    ExactBox box(const Array<DyadicType>& l, const Array<DyadicType>& u) const;
    ExactBox box(const GridCell& cell) const;
  private:
    // Create new data
    Void _create(const Vector<RawFloat64>& o, const Vector<RawFloat64>& l);
  private:
    // Pointer to data. We can test grids for equality using reference semantics since data is a constant.
    std::shared_ptr<Data> _data;
};


} // namespace Ariadne

#endif /* ARIADNE_GRID_H */

