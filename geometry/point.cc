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

#include "standard.h"
#include "config.h"

#include <cstdarg>
#include <sstream>
#include <string>
#include <vector>

#include "point.h"
#include "box.h"
#include "stlio.h"

namespace Ariadne {

inline ExactInterval make_interval(ApproximateNumber x) { return ExactInterval(x.raw()); }
inline ExactInterval make_interval(ExactNumber x) { return ExactInterval(x.raw()); }

template<class X> Point<X>::Point(std::initializer_list<double> lst)
    : Vector<X>(Vector<Float>(Vector<double>(lst)))
{
}

template<class X> Point<X>::Point(const std::string& str)
{
    *this=make_point(str);
}

template<class X> Point<X>* Point<X>::clone() const {
    return new Point<X>(*this);
}

UpperInterval operator+(ApproximateNumber x, UpperInterval y) { return UpperInterval(x.raw())+y; }

template<class X> ExactBox Point<X>::bounding_box() const {
    ExactBox r(this->dimension());
    Float e=eps();
    for(uint i=0; i!=this->dimension(); ++i) {
        r[i]=make_exact_interval((*this)[i]+UpperInterval(-e,+e)); }
    return r;
}

template<class X> void Point<X>::draw(CanvasInterface& canv, const Projection2d& proj) const {
    canv.dot(numeric_cast<double>((*this)[proj.x_coordinate()]),numeric_cast<double>((*this)[proj.y_coordinate()]));
    canv.stroke();
}

template class Point<ExactNumber>;
template class Point<ApproximateNumber>;


ExactPoint make_point(const std::string& str)
{
    std::vector<Float> lst;
    std::stringstream ss(str);
    read_sequence(ss,lst,'(',')',',');
    Vector<Float> vec(lst);

    return ExactPoint(Vector<ExactFloat>(lst));
}

} //namespace Ariadne
