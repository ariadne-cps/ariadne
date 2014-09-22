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
    typedef double dyadic_type;
    typedef int integer_type;
    typedef Float real_type;
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
    explicit Grid(uint d, Float l);

    //! Construct from a vector of offsets.
    explicit Grid(const Vector<Float>& lengths);

    //! Construct from a centre point and a vector of offsets.
    explicit Grid(const Vector<Float>& origin, const Vector<Float>& lengths);

    //! Copy constructor. Copies a reference to the grid data.
    Grid(const Grid& g);

    //! The underlying dimension of the grid.
    uint dimension() const;

    //! Tests equality of two grids. Tests equality of references first.
    bool operator==(const Grid& g) const;

    //! Tests inequality of two grids.
    bool operator!=(const Grid& g) const;

    //! The origin of the grid.
    const Vector<Float>& origin() const;

    //! The strides between successive integer points.
    const Vector<Float>& lengths() const;

    //! Write to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Grid& g);

    Float coordinate(uint d, dyadic_type x) const;
    Float subdivision_coordinate(uint d, dyadic_type x) const;
    Float subdivision_coordinate(uint d, integer_type n) const;

    int subdivision_index(uint d, const Float& x) const;
    int subdivision_lower_index(uint d, const Float& x) const;
    int subdivision_upper_index(uint d, const Float& x) const;

    Array<double> index(const Vector<Float>& pt) const;
    Array<double> lower_index(const Vector<Interval>& bx) const;
    Array<double> upper_index(const Vector<Interval>& bx) const;

    Vector<Float> point(const Array<int>& a) const;
    Vector<Float> point(const Array<double>& a) const;
    Box box(const Array<double>& l, const Array<double>& u) const;
    Box box(const GridCell& cell) const;
  private:
    // Create new data
    void _create(const Vector<Float>& o, const Vector<Float>& l);
  private:
    // Pointer to data. We can test grids for equality using reference semantics since data is a constant.
    std::shared_ptr<Data> _data;
};


} // namespace Ariadne

#endif /* ARIADNE_GRID_H */

