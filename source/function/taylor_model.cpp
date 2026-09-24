/***************************************************************************
 *            function/taylor_model.cpp
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

#include "numeric/numeric.hpp"
#include "function/taylor_model.tpl.hpp"

namespace Ariadne {

namespace {
Bool g_taylor_model_product_profile_enabled=false;
TaylorModelProductProfileContext g_taylor_model_product_profile_context=
    TaylorModelProductProfileContext::GENERAL;
TaylorModelProductProfileSnapshot g_taylor_model_product_profile;

Void accumulate_product_profile(
    TaylorModelProductProfileCounters& counters,
    unsigned long long product_pairs,
    unsigned long long sweep_passes,
    unsigned long long sweep_input_terms,
    unsigned long long sweep_output_terms,
    unsigned long long maximum_sweep_input_terms,
    unsigned long long maximum_sweep_output_terms,
    double elapsed_seconds)
{
    ++counters.calls;
    counters.product_pairs+=product_pairs;
    counters.sweep_passes+=sweep_passes;
    counters.sweep_input_terms+=sweep_input_terms;
    counters.sweep_output_terms+=sweep_output_terms;
    counters.swept_terms+=sweep_input_terms-sweep_output_terms;
    counters.maximum_sweep_input_terms=
        std::max(counters.maximum_sweep_input_terms,maximum_sweep_input_terms);
    counters.maximum_sweep_output_terms=
        std::max(counters.maximum_sweep_output_terms,maximum_sweep_output_terms);
    counters.elapsed_seconds+=elapsed_seconds;
}
} // namespace

Bool taylor_model_product_profile_enabled() {
    return g_taylor_model_product_profile_enabled;
}

Void set_taylor_model_product_profile_enabled(Bool enabled) {
    g_taylor_model_product_profile_enabled=enabled;
}

Void reset_taylor_model_product_profile() {
    g_taylor_model_product_profile=TaylorModelProductProfileSnapshot();
    g_taylor_model_product_profile_context=TaylorModelProductProfileContext::GENERAL;
}

TaylorModelProductProfileSnapshot taylor_model_product_profile_snapshot() {
    return g_taylor_model_product_profile;
}

TaylorModelProductProfileContext taylor_model_product_profile_context() {
    return g_taylor_model_product_profile_context;
}

Void set_taylor_model_product_profile_context(TaylorModelProductProfileContext context) {
    g_taylor_model_product_profile_context=context;
}

Void record_taylor_model_product_profile(
    TaylorModelProductProfileContext context,
    unsigned long long product_pairs,
    unsigned long long sweep_passes,
    unsigned long long sweep_input_terms,
    unsigned long long sweep_output_terms,
    unsigned long long maximum_sweep_input_terms,
    unsigned long long maximum_sweep_output_terms,
    double elapsed_seconds)
{
    auto& counters=(context==TaylorModelProductProfileContext::COMPOSE)
        ? g_taylor_model_product_profile.compose
        : g_taylor_model_product_profile.general;
    accumulate_product_profile(
        counters,product_pairs,sweep_passes,
        sweep_input_terms,sweep_output_terms,
        maximum_sweep_input_terms,maximum_sweep_output_terms,
        elapsed_seconds);
}

template<> String class_name<UnknownError<FloatDP>>() { return "UnknownError<FloatDP>"; }
template<> String class_name<UnknownError<FloatMP>>() { return "UnknownError<FloatMP>"; }


template class Series<FloatDPBounds>;
template class Series<FloatMPBounds>;

template class SweeperBase<FloatDP>;
template class SweeperBase<FloatMP>;
template class RelativeSweeperBase<FloatDP>;
template class RelativeSweeperBase<FloatMP>;

template class TaylorModel<ValidatedTag,FloatDP>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatDP>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatDP>>;

template class TaylorModel<ValidatedTag,FloatDPBounds>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatDPBounds>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatDPBounds>>;

template class TaylorModel<ApproximateTag,FloatDP>;
template struct AlgebraOperations<TaylorModel<ApproximateTag,FloatDP>>;
template class NormedAlgebraOperations<TaylorModel<ApproximateTag,FloatDP>>;


template class TaylorModel<ValidatedTag,FloatMP>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatMP>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatMP>>;

template class TaylorModel<ValidatedTag,FloatMPBounds>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatMPBounds>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatMPBounds>>;

template class TaylorModel<ApproximateTag,FloatMP>;
template struct AlgebraOperations<TaylorModel<ApproximateTag,FloatMP>>;
template class NormedAlgebraOperations<TaylorModel<ApproximateTag,FloatMP>>;


template class TaylorModel<ValidatedTag,FloatDPUpperInterval>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatDPUpperInterval>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatDPUpperInterval>>;

template class TaylorModel<ValidatedTag,FloatMPUpperInterval>;
template struct AlgebraOperations<TaylorModel<ValidatedTag,FloatMPUpperInterval>>;
template class NormedAlgebraOperations<TaylorModel<ValidatedTag,FloatMPUpperInterval>>;

//template TaylorModel<ValidatedTag,FloatDP> value_coefficients(TaylorModel<ValidatedTag,FloatDPBounds> const&);
//template TaylorModel<ValidatedTag,FloatDPBounds> exact_coefficients(TaylorModel<ValidatedTag,FloatDPBounds> const&);
//template TaylorModel<ValidatedTag,FloatMP> value_coefficients(TaylorModel<ValidatedTag,FloatMPBounds> const&);
//template TaylorModel<ValidatedTag,FloatMPBounds> exact_coefficients(TaylorModel<ValidatedTag,FloatMPBounds> const&);



} // namespace Ariadne
