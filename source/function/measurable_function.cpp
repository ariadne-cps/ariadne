/***************************************************************************
 *            measurable_function.cpp
 *
 *  Copyright  2020-21  Pieter Collins
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

#include "measurable_function.hpp"
#include "measurable_function.tpl.hpp"

#include "../geometry/union_of_intervals.tpl.hpp"

namespace Ariadne {

template ValidatedOpenSet<Real> preimage(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedOpenSet<Real> const& ops, IntervalDomainType dom, Accuracy acc);

template UnionOfIntervals<Value<FloatDP>> preimage_intervals(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedRegularSet<Real> rgps, Interval<Value<FloatDP>> ivl, Accuracy acc);
template UnionOfIntervals<Value<FloatDP>> preimage_intervals(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedOpenSet<Real> ops, Interval<Value<FloatDP>> ivl, Accuracy acc);
template UnionOfIntervals<LowerBound<FloatDP>> preimage_intervals(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedOpenSet<Real> ops, Interval<LowerBound<FloatDP>> ivl, Accuracy acc);


} // namespace Ariadne
