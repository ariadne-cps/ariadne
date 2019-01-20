/***************************************************************************
 *            point.hpp
 *
 *  Copyright 2008-17  Alberto Casagrande, Pieter Collins
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

/*! \file point.hpp
 *  \brief Points in Euclidean space.
 */

#ifndef ARIADNE_POINT_HPP
#define ARIADNE_POINT_HPP

#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"

#include "../output/graphics_interface.hpp"

namespace Ariadne {

template<class X> class Point;
typedef Point<ExactNumericType> ExactPoint;
typedef Point<EffectiveNumericType> EffectivePoint;
typedef Point<ValidatedNumericType> ValidatedPoint;
typedef Point<ApproximateNumericType> ApproximatePoint;
typedef Point<Real> RealPoint;

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

    //! Write to an output stream.
    virtual OutputStream& write(OutputStream& os) const {
        return os << static_cast<const Vector<RealType>&>(*this); }

    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;
    virtual ExactBoxType bounding_box() const;
};

template<class X> Point(Vector<X>) -> Point<X>;

ExactPoint make_point(const StringType&);

} // namespace Ariadne

#endif // ARIADNE_POINT_HPP
