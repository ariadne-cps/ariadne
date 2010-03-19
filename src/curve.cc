/***************************************************************************
 *            curve.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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

#include <cassert>

#include "vector.h"
#include "matrix.h"
#include "point.h"
#include "box.h"

#include "curve.h"


namespace Ariadne {


Curve::~Curve()
{
}



Curve::Curve(const VectorFunction& f)
    : _function(f)
{
    assert(this->_function.argument_size()==1);
}


Curve::Curve(const Curve& c)
    : _function(c._function)
{
}


Curve*
Curve::clone() const
{
    return new Curve(*this);
}


uint
Curve::dimension() const
{
    return this->_function.result_size();
}


ushort
Curve::smoothness() const
{
    return SMOOTH;
}



Point
Curve::value(const Float& s) const
{
    Vector<Float> v(1,&s);
    return Point(this->_function.evaluate(v));
}


Vector< Float >
Curve::tangent(const Float& s) const
{
    Vector<Float> v(1,&s);
    return column(this->_function.jacobian(v),0);
}



std::ostream&
Curve::write(std::ostream& os) const
{
    return os << "Curve( function=" << this->_function << " )";
}



InterpolatedCurve*
InterpolatedCurve::clone() const
{
    return new InterpolatedCurve(*this);
}

Box
InterpolatedCurve::bounding_box() const
{
    Box bx(this->_points.begin()->second);
    for(const_iterator iter=this->_points.begin(); iter!=this->_points.end(); ++iter) {
        bx=hull(bx,iter->second);
    }
    return bx;
}

void
InterpolatedCurve::draw(CanvasInterface& c) const
{
    uint xi=c.x_coordinate(); uint yi=c.y_coordinate();
    const_iterator iter=this->begin();
    const Point& pt=iter->second;
    c.move_to(pt[xi],pt[yi]);
    while(iter!=this->end()) {
        ++iter;
        const Point& pt=iter->second;
        c.line_to(pt[xi],pt[yi]);
    }
    if(this->begin()->second==(this->end())->second) {
        c.fill();
    } else {
        c.stroke();
    }
}

std::ostream&
InterpolatedCurve::write(std::ostream& os) const
{
    return os << (*this);
}


} // namespace Ariadne