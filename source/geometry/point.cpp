/***************************************************************************
 *            point.cpp
 *
 *  Copyright 2008--17  Alberto Casagrande, Pieter Collins
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

#include "../config.hpp"

#include "../geometry/point.tpl.hpp"

namespace Ariadne {

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
