/***************************************************************************
 *            hybrid/hybrid_set.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "config.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "hybrid/hybrid_graphics.hpp"
#include "drawing/geometry_draw.hpp"
#include "numeric/casts.hpp"

namespace Ariadne {

Orbit<HybridApproximatePoint>::Orbit(const HybridApproximatePoint& hpt)
    : _curves_ptr(new List<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(0,hpt.point()))))
{ }

Orbit<HybridApproximatePoint>::Orbit(List<HybridInterpolatedCurve> hcrvs)
    : _curves_ptr(new List<HybridInterpolatedCurve>(std::move(hcrvs)))
{ }

SizeType
Orbit<HybridApproximatePoint>::size() const
{
    return this->_curves_ptr->size();
}

const InterpolatedCurve&
Orbit<HybridApproximatePoint>::curve(Nat m) const
{
    return (*this->_curves_ptr)[m].euclidean_set();
}

Void
Orbit<HybridApproximatePoint>::insert(HybridTime ht, const HybridApproximatePoint& hpt)
{
    ARIADNE_ASSERT(ht.discrete_time()<=this->size());
    Real time=ht.continuous_time();
    FloatDP flt_time=cast_exact(time.get(dp));
    if(this->size()==ht.discrete_time()) {
        this->_curves_ptr->push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(flt_time,hpt.point())));
    } else {
        (*this->_curves_ptr)[static_cast<unsigned int>(ht.discrete_time().get_si())].euclidean_set().insert(flt_time,hpt.point());
    }
}

Void Orbit<HybridApproximatePoint>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<HybridApproximatePoint>& orbit=*this;
    for(SizeType i=0; i!=orbit._curves_ptr->size(); ++i) {
        HybridInterpolatedCurve const& hcurve=this->_curves_ptr->at(i);
        if(locations.empty() || locations.contains(hcurve.location())) {
            RealSpace const& space=hcurve.space();
            if(valid_axis_variables(space,axes)) {
                Ariadne::draw(canvas,projection(space,axes),hcurve.euclidean_set());
            }
        }
    }
}

Orbit<Vector<HybridApproximatePoint>>::Orbit(const Vector<HybridApproximatePoint>& vhpt)
    : _curves_ptr(new Vector(vhpt.size(), List<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(vhpt.at(0).location(),vhpt.at(0).space(),InterpolatedCurve(0,vhpt.at(0).point())))))
{
    for(SizeType i=1; i<vhpt.size(); i++)
    {
        this->_curves_ptr->at(i) = List<HybridInterpolatedCurve>(1u,HybridInterpolatedCurve(vhpt.at(i).location(),vhpt.at(i).space(),InterpolatedCurve(0,vhpt.at(i).point())));
    }
}

SizeType
Orbit<Vector<HybridApproximatePoint>>::size(SizeType setNumber)
{
    return this->_curves_ptr->at(setNumber).size();
}

Void
Orbit<Vector<HybridApproximatePoint>>::insert(HybridTime ht, const HybridApproximatePoint& hpt, SizeType setNumber)
{
    ARIADNE_ASSERT(ht.discrete_time()<=this->size(setNumber));
    Real time=ht.continuous_time();
    FloatDP flt_time=cast_exact(time.get(dp));
    if(this->size(setNumber)==ht.discrete_time()) {
        this->_curves_ptr->at(setNumber).push_back(HybridInterpolatedCurve(hpt.location(),hpt.space(),InterpolatedCurve(flt_time,hpt.point())));
    } else {
        this->_curves_ptr->at(setNumber)[static_cast<unsigned int>(ht.discrete_time().get_si())].euclidean_set().insert(flt_time,hpt.point());
    }
}

const InterpolatedCurve&
Orbit<Vector<HybridApproximatePoint>>::curve(Nat setNumber, Nat step)
{
    return this->_curves_ptr->at(setNumber)[step].euclidean_set();
}

Void
Orbit<Vector<HybridApproximatePoint>>::draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes) const {
    const Orbit<Vector<HybridApproximatePoint>>& orbit=*this;
    for(SizeType i=0; i<orbit._curves_ptr->size(); i++)
    {
        for(SizeType j=0; j!=orbit._curves_ptr->at(i).size(); ++j) {
            HybridInterpolatedCurve const& hcurve=this->_curves_ptr->at(i).at(j);
            if(locations.empty() || locations.contains(hcurve.location())) {
                RealSpace const& space=hcurve.space();
                if(valid_axis_variables(space,axes)) {
                    Ariadne::draw(canvas,projection(space,axes),hcurve.euclidean_set());
                }
            }
        }
    }
}



template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridApproximatePoint >& orb)
{
    return os << orb.curves();
}




Void Orbit<HybridEnclosure>::draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const {
    this->reach().draw(c,l,v);
}

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridEnclosure >& orb)
{
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}


} // namespace Ariadne
