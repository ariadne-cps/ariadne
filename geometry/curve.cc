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



Curve::Curve(const Function<Effective,IntervalDomain,BoxDomain>& f)
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


DimensionType
Curve::dimension() const
{
    return this->_function.result_size();
}


DegreeType
Curve::smoothness() const
{
    return SMOOTH;
}



Curve::PointType
Curve::value(const ParameterType& s) const
{
    Vector<ApproximateFloat64> fv=this->_function.evaluate(s);
    return PointType(fv);
}


Curve::TangentVectorType
Curve::tangent(const ParameterType& s) const
{
    return Ariadne::tangent(this->_function,s);
}



OutputStream&
Curve::write(OutputStream& os) const
{
    return os << "Curve( function=" << this->_function << " )";
}



InterpolatedCurve*
InterpolatedCurve::clone() const
{
    return new InterpolatedCurve(*this);
}

Void
InterpolatedCurve::insert(const ParameterType& s, const PointType& pt) {
    if(!this->_points.empty()) { ARIADNE_ASSERT(pt.size()==this->dimension()); }
    this->_points.insert(Pair< ParameterType, PointType >(s,pt));
}

Void
InterpolatedCurve::insert(const RawFloat64& s, const Vector<RawFloat64>& pt) {
    if(!this->_points.empty()) { ARIADNE_ASSERT(pt.size()==this->dimension()); }
    this->insert(ParameterType(s),PointType(reinterpret_cast<Vector<ExactFloat64>const&>(pt)));
}

UpperBox
InterpolatedCurve::bounding_box() const
{
    ExactPoint pt=cast_exact(this->_points.begin()->second);
    ExactBox bx(pt);
    for(ConstIterator iter=this->_points.begin(); iter!=this->_points.end(); ++iter) {
        pt=cast_exact(iter->second);
        bx=hull(bx,pt);
    }
    return bx;
}

Void
InterpolatedCurve::draw(CanvasInterface& c, const Projection2d& p) const
{
    Nat xi=p.x_coordinate(); Nat yi=p.y_coordinate();
    ConstIterator iter=this->begin();
    ApproximatePoint pt=join(iter->second,iter->first);
    c.move_to(pt[xi],pt[yi]);
    while(iter!=this->end()) {
        ++iter;
        pt=join(iter->second,iter->first);
        c.line_to(pt[xi],pt[yi]);
    }
    if( decide(this->begin()->second==(this->end())->second) ) {
        c.fill();
    } else {
        c.stroke();
    }
}

OutputStream&
InterpolatedCurve::write(OutputStream& os) const
{
    return os << (*this);
}


} // namespace Ariadne
