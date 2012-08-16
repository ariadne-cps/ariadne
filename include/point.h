/***************************************************************************
 *            point.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

/*! \file point.h
 *  \brief Points in Euclidean space.
 */

#ifndef ARIADNE_POINT_H
#define ARIADNE_POINT_H

#include "numeric.h"
#include "vector.h"

#include "graphics_interface.h"

namespace Ariadne {

//! A point in Euclidean space.
class Point
    : public Vector<Float>
    , public DrawableInterface
{
  public:
    typedef Float real_type;
    //! Default constructor contructs the singleton point in zero dimensions.
    Point() : Vector<Float>() { }
    //! The origin in \a n dimensions.
    explicit Point(uint n) : Vector<Float>(n) { }
    //! Construct from a string literal of the form "(x1,x2,...,xd)".
    explicit Point(const std::string& str);
    Point(const Vector<Float>& v) : Vector<Float>(v) { }
    template<class E> Point(const VectorExpression<E>& ve) : Vector<Float>(ve) { }
    template<class T1, class T2> Point(const T1& t1, const T2& t2) : Vector<Float>(t1,t2) { }
    //! Construct from an initializer list of floating-point values.
    explicit Point(std::initializer_list<Float> lst);
    //! The origin in \a n dimensions.
    static Point origin(uint n) { return Point(n,0.0); }
    //! A dynamically-allocated copy.
    virtual Point* clone() const;
    //! The dimension of the point.
    uint dimension() const { return this->size(); }
    //! An explicit cast to a float vector. Useful to prevent ambiguous function overloads.
    const Vector<Float>& vector() const { return *this; }

    Vector<Float> centre() const { return *this; }

    //! Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const {
        return os << static_cast<const Vector<Float>&>(*this); }

    virtual void draw(CanvasInterface& c, const Projection2d& p) const;
    virtual Box bounding_box() const;
};

Point make_point(const std::string&);

} // namespace Ariadne

#endif // ARIADNE_POINT_H
