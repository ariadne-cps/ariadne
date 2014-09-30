/***************************************************************************
 *            curve.cc
 *
 *  Copyright  2007  Pieter Collins
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

#include "functional.h"

#include "config.h"

#include <cassert>

#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "point.h"
#include "box.h"

#include "curve.h"


namespace Ariadne {


Curve::~Curve()
{
}



Curve::Curve(const EffectiveVectorFunction& f)
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



Curve::PointType
Curve::value(const ParameterType& s) const
{
    Vector<ApproximateFloatType> v(1,ApproximateFloatType(s));
    Vector<ApproximateFloatType> fv=this->_function.evaluate(v);
    return PointType(fv);
}


Vector< Float >
Curve::tangent(const Float& s) const
{
    Vector<ExactFloatType> v(1,ExactFloatType(s));
    auto col=column(this->_function.jacobian(v),0);
    return reinterpret_cast<Vector<Float>const&>(col);
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

void
InterpolatedCurve::insert(const ParameterType& s, const PointType& pt) {
    if(!this->_points.empty()) { ARIADNE_ASSERT(pt.size()==this->dimension()); }
    this->_points.insert(std::pair< ParameterType, PointType >(s,pt));
}

void
InterpolatedCurve::insert(const ExactFloatType& s, const ExactPointType& pt) {
    if(!this->_points.empty()) { ARIADNE_ASSERT(pt.size()==this->dimension()); }
    this->_points.insert(std::pair< ParameterType, PointType >(reinterpret_cast<ParameterType const&>(s),reinterpret_cast<PointType const&>(pt)));
}

Box
InterpolatedCurve::bounding_box() const
{
    Point pt=this->_points.begin()->second;
    Box bx(pt);
    for(const_iterator iter=this->_points.begin(); iter!=this->_points.end(); ++iter) {
        pt=iter->second;
        bx=hull(bx,pt);
    }
    return bx;
}

void
InterpolatedCurve::draw(CanvasInterface& c, const Projection2d& p) const
{
    uint xi=p.x_coordinate(); uint yi=p.y_coordinate();
    const_iterator iter=this->begin();
    Point pt=join(iter->second,iter->first);
    c.move_to(pt[xi],pt[yi]);
    while(iter!=this->end()) {
        ++iter;
        pt=join(iter->second,iter->first);
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
