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

#include <cstdarg>
#include <sstream>
#include <string>
#include <vector>

#include "config.h"

#include "point.h"
#include "box.h"
#include "stlio.h"

namespace Ariadne {


Point::Point(std::initializer_list<Float> lst)
    : Vector<Float>(lst)
{
}

Point::Point(const std::string& str)
{
    *this=make_point(str);
}

Point* Point::clone() const {
    return new Point(*this);
}

Box Point::bounding_box() const {
    Box r(this->dimension());
    Float e=eps();
    for(uint i=0; i!=this->dimension(); ++i) {
        r[i]=(*this)[i]+Interval(-e,+e); }
    return r;
}

void Point::draw(CanvasInterface& canv, const Projection2d& proj) const {
    canv.dot(numeric_cast<double>((*this)[proj.x_coordinate()]),numeric_cast<double>((*this)[proj.y_coordinate()]));
    canv.stroke();
}

Point make_point(const std::string& str)
{
    std::vector<Float> vec;
    std::stringstream ss(str);
    read_sequence(ss,vec,'(',')',',');
    return Point(vec.size(),&vec[0]);
}

} //namespace Ariadne
