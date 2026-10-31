/***************************************************************************
 *            dynamics/orbit.cpp
 *
 *  Copyright  2007-20  Pieter Collins
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

#include "function/functional.hpp"


#include <utility>

#include "evolution/orbit.hpp"

#include "geometry/box.hpp"
#include "geometry/point.hpp"
#include "set/curve.hpp"
#include "set/function_set.hpp"
#include "set/list_set.hpp"
#include "rendering/figure.hpp"

namespace Ariadne {


template<class F> Orbit<ExactPoint<F>>::Orbit(const ExactPoint<F>& pt)
    : _curve(new InterpolatedCurve(0,pt))
{ }

template<class F> Void
Orbit<ExactPoint<F>>::insert(F t, const ExactPoint<F>& pt)
{
    _curve->insert(t,pt);
}


template<class F> Orbit<List<Point<Approximation<F>>>>::Orbit(List<Point<Approximation<F>>> const& ptLst)
    : _curves(new List<InterpolatedCurve>())
{
    for (SizeType i=0; i<ptLst.size(); ++i) {
        _curves->push_back(InterpolatedCurve(ptLst.at(i)));
    }
}

template<class F> Void Orbit<List<Point<Approximation<F>>>>::insert(F t, const Point<Approximation<F>>& pt, SizeType curveNumber)
{
    _curves->at(curveNumber).insert(t,pt);
}


template<class F> Orbit<List<LabelledPoint<Approximation<F>>>>::Orbit(const List<LabelledPoint<Approximation<F>>>& ptLst)
    : _curves(new List<LabelledInterpolatedCurve>())
{
    for (SizeType i=0; i<ptLst.size(); ++i) {
        _curves->push_back(LabelledInterpolatedCurve(ptLst.at(i)));
    }
}

template<class F> Void Orbit<List<LabelledPoint<Approximation<F>>>>::insert(F t, const LabelledPoint<Approximation<F>>& pt, SizeType curveNumber) {
    _curves->at(curveNumber).insert(static_cast<const RawFloatDP>(t), pt);
}


template class Orbit<Point<Approximation<FloatDP>>>;
template class Orbit<List<Point<Approximation<FloatDP>>>>;
template class Orbit<LabelledPoint<Approximation<FloatDP>>>;
template class Orbit<List<LabelledPoint<Approximation<FloatDP>>>>;




} // namespace Ariadne
