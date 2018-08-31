/***************************************************************************
 *            point.cpp
 *
 *  Copyright 2008--17  Alberto Casagrande, Pieter Collins
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

#include "../utility/standard.hpp"
#include "../config.hpp"

#include <cstdarg>
#include <sstream>
#include <string>
#include <vector>

#include "../geometry/point.hpp"
#include "../geometry/box.hpp"
#include "../utility/stlio.hpp"

namespace Ariadne {

inline ExactIntervalType make_interval(ApproximateNumericType x) { return ExactIntervalType(x.raw(),x.raw()); }
inline ExactIntervalType make_interval(ExactNumericType x) { return ExactIntervalType(x.raw(),x.raw()); }

template<class X, EnableIf<IsConstructible<X,ExactDouble>> =dummy> Vector<X> make_vector(InitializerList<double> lst) {
    return Vector<X>(Array<ExactDouble>(Array<double>(lst))); }

template<class X, EnableIf<IsConstructible<X,ExactDouble,DoublePrecision>> =dummy> Vector<X> make_vector(InitializerList<double> lst) {
    return Vector<X>(lst,dp); }

template<class X> Point<X>::Point(InitializerList<double> lst)
    : Vector<X>(make_vector<X>(lst))
{
}

template<class X> Point<X>* Point<X>::clone() const {
    return new Point<X>(*this);
}

template<class X> ExactBoxType Point<X>::bounding_box() const {
    typedef decltype(declval<X>()+declval<FloatDPError>()) U;
    FloatDPError eps(FloatDP::eps(dp));
    ExactBoxType r(this->dimension());
    UpperIntervalType v;
    for(Nat i=0; i!=this->dimension(); ++i) {
        X const& pti=(*this)[i];
        r[i]=cast_exact_interval(Interval<U>(pti-eps,pti+eps)); }
    return r;
}

template<class X> Void Point<X>::draw(CanvasInterface& canv, const Projection2d& proj) const {
    canv.dot(numeric_cast<double>((*this)[proj.x_coordinate()]),numeric_cast<double>((*this)[proj.y_coordinate()]));
    canv.stroke();
}

template<> String class_name<Point<FloatDPApproximation>>() { return "Point<FloatDPApproximation>"; }
template<> String class_name<Point<FloatDPBounds>>() { return "Point<FloatDPBounds>"; }
template<> String class_name<Point<FloatDPValue>>() { return "Point<FloatDPValue>"; }
template<> String class_name<Point<Real>>() { return "Point<Real>"; }

template class Point<ExactNumericType>;
template class Point<EffectiveNumericType>;
template class Point<ValidatedNumericType>;
template class Point<ApproximateNumericType>;


ExactPoint make_point(const StringType& str)
{
    std::vector<FloatDP> lst;
    StringStream ss(str);
    read_sequence(ss,lst,'(',')',',');
    Vector<FloatDP> vec(lst);

    return ExactPoint(Vector<FloatDPValue>(vec));
}

} //namespace Ariadne
