/***************************************************************************
 *            point.cc
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

#include "utility/standard.h"
#include "config.h"

#include <cstdarg>
#include <sstream>
#include <string>
#include <vector>

#include "geometry/point.h"
#include "geometry/box.h"
#include "utility/stlio.h"

namespace Ariadne {

inline ExactIntervalType make_interval(ApproximateNumericType x) { return ExactIntervalType(x.raw(),x.raw()); }
inline ExactIntervalType make_interval(ExactNumericType x) { return ExactIntervalType(x.raw(),x.raw()); }

template<class X> Point<X>::Point(InitializerList<double> lst)
    : Vector<X>(Vector<Float64Value>(Vector<Float64>(lst)))
{
}

template<class X> Point<X>* Point<X>::clone() const {
    return new Point<X>(*this);
}

template<class X> ExactBoxType Point<X>::bounding_box() const {
    Float64 eps=Float64::eps(Precision64());
    ExactBoxType r(this->dimension());
    UpperIntervalType e(-eps,+eps);
    for(Nat i=0; i!=this->dimension(); ++i) {
        r[i]=cast_exact_interval(cast_exact((*this)[i])+e); }
    return r;
}

template<class X> Void Point<X>::draw(CanvasInterface& canv, const Projection2d& proj) const {
    canv.dot(numeric_cast<double>((*this)[proj.x_coordinate()]),numeric_cast<double>((*this)[proj.y_coordinate()]));
    canv.stroke();
}

template class Point<ExactNumericType>;
template class Point<EffectiveNumericType>;
template class Point<ValidatedNumericType>;
template class Point<ApproximateNumericType>;


ExactPoint make_point(const StringType& str)
{
    std::vector<Float64> lst;
    StringStream ss(str);
    read_sequence(ss,lst,'(',')',',');
    Vector<Float64> vec(lst);

    return ExactPoint(Vector<Float64Value>(vec));
}

} //namespace Ariadne
