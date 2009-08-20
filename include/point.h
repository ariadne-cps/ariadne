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
    //! Construct from a string literal of the form "(x1,x2,...,xd)".
    explicit Point(const std::string& str);
    template<class T> Point(const T& t) : Vector<Float>(t) { }
    template<class T1, class T2> Point(const T1& t1, const T2& t2) : Vector<Float>(t1,t2) { }
    //! Construct from an integer giving the dimension and a list of floating-point values 
    //! giving the values.
    explicit Point(uint d, const Float& x0, ...);
    //! The origin in \a n dimensions.
    static Point origin(uint n) { return Point(n,0.0); }
    //! A dynamically-allocated copy.
    virtual Point* clone() const;
    //! The dimension of the point.
    uint dimension() const { return this->size(); }
    Vector<Float> centre() const { return *this; }

    //! Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const {
        return os << static_cast<const Vector<Float>&>(*this); }

    virtual void draw(CanvasInterface& c) const;
    virtual Box bounding_box() const;
};

Point make_point(const std::string&);

} // namespace Ariadne

#endif // ARIADNE_POINT_H
