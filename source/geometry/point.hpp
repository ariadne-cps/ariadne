/***************************************************************************
 *            geometry/point.hpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file geometry/point.hpp
 *  \brief Points in Euclidean space.
 */

#ifndef ARIADNE_POINT_HPP
#define ARIADNE_POINT_HPP

#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"

#include "../output/graphics_interface.hpp"

namespace Ariadne {

template<class X> class Point;

//@{
//! \relates Point
//! \name Type synonyms
using DyadicPoint = Point<Dyadic>; //!< .
using RationalPoint = Point<Rational>; //!< .
using RealPoint = Point<Real>; //!< .

template<class F> using ExactPoint = Point<Value<F>>;
template<class F> using ValidatedPoint = Point<Bounds<F>>;
template<class F> using ApproximatePoint = Point<Approximation<F>>;

using FloatDPValuePoint = Point<FloatDPValue>;
using FloatDPBoundsPoint = Point<FloatDPBounds>;
using FloatDPApproximationPoint = Point<FloatDPApproximation>;

typedef ExactPoint<FloatDP> ExactPointType;
//@}

//! A point in Euclidean space.
template<class X>
class Point
    : public Vector<X>
    , public DrawableInterface
{
  public:
    typedef X RealType;
    //! Default constructor contructs the singleton point in zero dimensions.
    Point() : Vector<RealType>() { }
    //! The origin in \a n dimensions.
    explicit Point(Nat n) : Vector<RealType>(n) { }
    Point(const Vector<RealType>& v) : Vector<RealType>(v) { }
    template<class T, EnableIf<IsConvertible<T,X>> =dummy> Point(const Point<T>& pt) : Vector<RealType>(pt.vector()) { }
    template<class Y, class PR, EnableIf<IsConstructible<X,Y,PR>> =dummy> Point(const Point<Y>& pt, PR pr) : Vector<RealType>(pt.vector(),pr) { }
    //! Construct from an initializer list of floating-point values.
    template<class T, EnableIf<IsConvertible<T,X>> =dummy> Point(SizeType n, const T& t) : Vector<RealType>(n,RealType(t)) { }
    //! Construct from an initializer list of floating-point values.
    explicit Point(InitializerList<double> lst);
    //! The origin in \a n dimensions.
    static Point origin(Nat n) { return Point(n,RealType(0)); }
    //! A dynamically-allocated copy.
    virtual Point<X>* clone() const;
    //! The dimension of the point.
    DimensionType dimension() const { return this->size(); }
    //! An explicit cast to a float vector. Useful to prevent ambiguous function overloads.
    const Vector<RealType>& vector() const { return *this; }

    Vector<RealType> centre() const { return *this; }

    friend Point<X> product(Point<X> const& pt1, Point<X> const& pt2) {
        return Point<X>(join(pt1.vector(),pt2.vector())); }

    //! Write to an output stream.
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<const Vector<RealType>&>(*this); }

    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;
    virtual FloatDPUpperBox bounding_box() const;
};

template<class X> Point(Vector<X>) -> Point<X>;

//Point<Real> make_point(const StringType&);

} // namespace Ariadne

#endif // ARIADNE_POINT_HPP
