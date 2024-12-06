/***************************************************************************
 *            function/taylor_function.cpp
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

#include "function/functional.hpp"
#include "config.hpp"

#include <iostream>
#include <iomanip>

#include "utility/macros.hpp"
#include "utility/exceptions.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "algebra/multi_index.hpp"
#include "function/polynomial.hpp"
#include "algebra/differential.hpp"
#include "algebra/evaluate.hpp"
#include "function/taylor_model.hpp"

#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/function_patch.hpp"
#include "function/scaled_function_patch.hpp"

#include "function/taylor_function.hpp"

#include "taylor_model.tpl.hpp"
#include "scaled_function_patch.tpl.hpp"
#include "function_mixin.tpl.hpp"

#define VOLATILE ;

namespace Ariadne {

//template<class S1, class S2> decltype(auto) superset(S1&& s1, S2&& s2) {
//    return subset(std::forward<S2>(s2),std::forward<S1>(s1)); }

#warning Move to FunctionPatch. Throw exception
template<class X> Void check_domain(IntervalDomainType const& dom, X const& x);
template<> Void check_domain(IntervalDomainType const& dom, ApproximateNumber const& x) {
    ARIADNE_ASSERT(probably(contains(static_cast<DyadicInterval>(dom),x))); }
template<> Void check_domain(IntervalDomainType const& dom, ValidatedNumber const& x) {
    ARIADNE_ASSERT(definitely(contains(static_cast<DyadicInterval>(dom),x))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatDPApproximation const& x) {
    ARIADNE_ASSERT(probably(contains(dom,x))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatMPApproximation const& x) {
    ARIADNE_ASSERT(probably(contains(dom,x))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatDPBounds const& x) {
    ARIADNE_ASSERT(definitely(contains(dom,x))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatMPBounds const& x) {
    ARIADNE_ASSERT(definitely(contains(dom,x))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatDPApproximationDifferential const& x) {
    ARIADNE_ASSERT(probably(contains(dom,x.value()))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatMPApproximationDifferential const& x) {
    ARIADNE_ASSERT(probably(contains(dom,x.value()))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatDPBoundsDifferential const& x) {
    ARIADNE_ASSERT(definitely(contains(dom,x.value()))); }
template<> Void check_domain(IntervalDomainType const& dom, FloatMPBoundsDifferential const& x) {
    ARIADNE_ASSERT(definitely(contains(dom,x.value()))); }
template<> Void check_domain(IntervalDomainType const& dom, ApproximateTaylorModelDP const& x) {
    ARIADNE_ASSERT(probably(subset(x.range(),dom))); }
template<> Void check_domain(IntervalDomainType const& dom, ApproximateTaylorModelMP const& x) {
    ARIADNE_ASSERT(probably(subset(Interval<FloatDPApproximation>(x.range()),dom))); }
template<> Void check_domain(IntervalDomainType const& dom, ValidatedTaylorModelDP const& tm) {
    if (not definitely(subset(tm.range(),dom))) {
        ARIADNE_THROW(DomainException,"check_domain(dom,tm) with dom="<<dom<<", tm.range()="<<tm.range()<<", tm=-"<<tm,
                      "tm.range() is not definitely a subset of dom"); } }
template<> Void check_domain(IntervalDomainType const& dom, ValidatedTaylorModelMP const& x) {
    ARIADNE_ASSERT(definitely(subset(Interval<FloatDPUpperBound>(x.range()),dom))); }
template<> Void check_domain(IntervalDomainType const& dom, ValidatedTaylorModel<FloatDPUpperInterval> const& x) {
    ARIADNE_ASSERT(definitely(subset(x.range(),dom))); }
template<> Void check_domain(IntervalDomainType const& dom, ValidatedIntervalTaylorModelMP const& x) {
    ARIADNE_ASSERT(definitely(subset(IntervalRangeType(x.range()),dom))); }
template<> Void check_domain(IntervalDomainType const& dom, ElementaryAlgebra<ApproximateNumber> const& x) {
    ARIADNE_NOT_IMPLEMENTED; }
template<> Void check_domain(IntervalDomainType const& dom, ElementaryAlgebra<ValidatedNumber> const& x) {
    ARIADNE_NOT_IMPLEMENTED; }
template<> Void check_domain(IntervalDomainType const& dom, Formula<ApproximateNumber> const& x) {
    ARIADNE_NOT_IMPLEMENTED; }
template<> Void check_domain(IntervalDomainType const& dom, Formula<ValidatedNumber> const& x) {
    ARIADNE_NOT_IMPLEMENTED; }


#warning Should not need check_domain for an unbounded function
template<> Void check_domain(IntervalDomainType const& ivl, ValidatedScalarMultivariateFunction const& sf) {
    ARIADNE_THROW(std::runtime_error,"check_domain(IntervalDomainType ivl, ValidatedScalarMultivariateFunction const& sf)",
                  "Domain of function (sf=" << sf <<") cannot be computed as a bounded set."); }

template class ScaledFunctionPatchFactory<ValidatedTaylorModelDP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelDP>,RealVector>;

template class ScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,RealScalar(RealVector)>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,RealScalar(RealVector)>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelDP>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ApproximateTag,RealVector(RealVector)>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedTaylorModelDP>,ValidatedTag,RealVector(RealVector)>;


template class ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelDP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelDP>,RealVector>;

template class ScaledFunctionPatch<ValidatedBoundsTaylorModelDP>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedBoundsTaylorModelDP>,ApproximateTag,RealScalar(RealVector)>;
template class FunctionMixin<ScaledFunctionPatch<ValidatedBoundsTaylorModelDP>,ValidatedTag,RealScalar(RealVector)>;
template class VectorScaledFunctionPatch<ValidatedBoundsTaylorModelDP>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedBoundsTaylorModelDP>,ApproximateTag,RealVector(RealVector)>;
template class FunctionMixin<VectorScaledFunctionPatch<ValidatedBoundsTaylorModelDP>,ValidatedTag,RealVector(RealVector)>;


template class ScaledFunctionPatchFactory<ValidatedTaylorModelMP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedTaylorModelMP>,RealVector>;

template class ScaledFunctionPatch<ValidatedTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedTaylorModelMP>;


template class ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelMP>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelMP>,RealVector>;

template class ScaledFunctionPatch<ValidatedBoundsTaylorModelMP>;
template class VectorScaledFunctionPatch<ValidatedBoundsTaylorModelMP>;



FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory(Sweeper<FloatDP> const& sweeper) {
    return new TaylorFunctionFactory(sweeper);
}
FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>* make_taylor_function_factory() {
    return make_taylor_function_factory(Sweeper<FloatDP>());
}
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory(Sweeper<FloatDP> const& sweeper) {
    return new TaylorFunctionFactory(sweeper);
}
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory() {
    return make_taylor_function_patch_factory(Sweeper<FloatDP>());
}

} // namespace Ariadne
