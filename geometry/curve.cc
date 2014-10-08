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

#include "function/functional.h"

#include "config.h"

#include <cassert>

#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/differential.h"
#include "geometry/point.h"
#include "geometry/box.h"

#include "geometry/curve.h"


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
    Vector<ApproximateFloat> v(1,ApproximateFloat(s));
    Vector<ApproximateFloat> fv=this->_function.evaluate(v);
    return PointType(fv);
}


Curve::TangentVectorType
Curve::tangent(const ParameterType& s) const
{
    TangentVectorType v(1,s);
    auto col=column(this->_function.jacobian(v),0);
    return col;
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
InterpolatedCurve::insert(const RawFloat& s, const Vector<RawFloat>& pt) {
    if(!this->_points.empty()) { ARIADNE_ASSERT(pt.size()==this->dimension()); }
    this->insert(ParameterType(s),PointType(reinterpret_cast<Vector<ExactFloat>const&>(pt)));
}

UpperBox
InterpolatedCurve::bounding_box() const
{
    ExactPoint pt=make_exact(this->_points.begin()->second);
    ExactBox bx(pt);
    for(const_iterator iter=this->_points.begin(); iter!=this->_points.end(); ++iter) {
        pt=make_exact(iter->second);
        bx=hull(bx,pt);
    }
    return bx;
}

void
InterpolatedCurve::draw(CanvasInterface& c, const Projection2d& p) const
{
    uint xi=p.x_coordinate(); uint yi=p.y_coordinate();
    const_iterator iter=this->begin();
    ApproximatePoint pt=join(iter->second,iter->first);
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
