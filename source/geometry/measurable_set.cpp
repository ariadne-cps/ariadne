/***************************************************************************
 *            measurable_set.cpp
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

#include "measurable_set.tpl.hpp"
#include "union_of_intervals.tpl.hpp"

namespace Ariadne {

template class UnionOfIntervals<Dyadic>;
template class UnionOfIntervals<Value<FloatDP>>;
template class UnionOfIntervals<LowerBound<FloatDP>>;

template class LowerMeasurableSetModel<UnionOfIntervals<Dyadic>,Dyadic>;
template class LowerMeasurableSetModel<UnionOfIntervals<Value<FloatDP>>,Error<FloatDP>>;
template class LowerMeasurableSetModel<UnionOfIntervals<LowerBound<FloatDP>>,Error<FloatDP>>;

} // namespace Ariadne
