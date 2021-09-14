/***************************************************************************
 *            measurable_function.tpl.hpp
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
#include "../geometry/union_of_intervals.hpp"

namespace Ariadne {

template<class PR> UnionOfIntervals<Value<PR>> preimage_intervals(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedRegularSet<Real> rgs, Interval<Value<PR>> ivl, Accuracy acc) {
    Interval<Value<PR>> rng=cast_exact(image(ivl,f));
    if (definitely(rgs.covers(rng))) {
        return {cast_exact_interval(ivl)};
    } else if (definitely(rgs.separated(rng))) {
        return {{},ivl.precision()};
    } else if (definitely(ivl.width()<acc.error())) {
        return {{},ivl.precision()};
    } else {
        Pair<Interval<Value<PR>>,Interval<Value<PR>>> subivls=split(ivl);
        return join(preimage_intervals(f,rgs,subivls.first,acc),preimage_intervals(f,rgs,subivls.second,acc));
    }
}

template<class PR> UnionOfIntervals<Value<PR>> preimage_intervals(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedOpenSet<Real> ops, Interval<Value<PR>> ivl, Accuracy acc) {
    if (definitely(ops.covers(cast_exact(image(ivl,f))))) {
        return {ivl};
    } else if (definitely(ivl.width()<acc.error())) {
        return {{},ivl.precision()};
    } else {
        Pair<Interval<Value<PR>>,Interval<Value<PR>>> subivls=split(ivl);
        return join(preimage_intervals(f,ops,subivls.first,acc),preimage_intervals(f,ops,subivls.second,acc));
    }
}

template<class PR> UnionOfIntervals<LowerBound<PR>> preimage_intervals(ValidatedContinuousFunction<Real(Real)> const& f, ValidatedOpenSet<Real> ops, Interval<LowerBound<PR>> ivl, Accuracy acc) {
    if (definitely(ops.covers(cast_exact(image(cast_exact(ivl),f))))) {
        return {ivl};
    } else if (definitely(cast_exact(ivl.width())<acc.error())) {
        return {{},cast_exact(ivl).precision()};
    } else {
        Pair<Interval<LowerBound<PR>>,Interval<LowerBound<PR>>> subivls=split(ivl);
        return join(preimage_intervals(f,ops,subivls.first,acc),preimage_intervals(f,ops,subivls.second,acc));
    }
}


template<class ARG, class RES> ValidatedOpenSet<ARG> preimage(ValidatedContinuousFunction<RES(ARG)> const& f, ValidatedOpenSet<RES> const& ops, DomainOfType<ARG> dom, Accuracy acc) {
    DoublePrecision pr;
    Interval<Value<FloatDP>> ivl=cast_exact(Interval<UpperBound<FloatDP>>(dom,pr));
    auto rgsp=dynamic_pointer_cast<ValidatedRegularSet<Real>::Interface>(ops.managed_pointer());
    if (false and rgsp) {
        auto rgs=ValidatedRegularSet<Real>(rgsp);
        return ValidatedOpenSet<Real>(wrap_open(preimage_intervals(f,rgs,ivl,acc)));
    } else {
        return ValidatedOpenSet<Real>(wrap_open(preimage_intervals(f,ops,ivl,acc)));
    }
}


} // namespace Ariadne
