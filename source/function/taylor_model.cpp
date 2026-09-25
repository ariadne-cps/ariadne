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
Bool g_taylor_model_early_discard_enabled=false;
Bool g_taylor_model_incremental_sweep_enabled=true;
Bool g_taylor_model_product_accumulator_enabled=false;
Bool g_taylor_model_dense_accumulator_enabled=false;
Bool g_taylor_model_dense_batched_rounding_enabled=true;
TaylorModelDenseWorkspaceStats g_taylor_model_dense_workspace_stats;
TaylorModelDenseHotLoopProfile g_taylor_model_dense_hot_loop_profile;
TaylorModelAccumulatorProfile g_taylor_model_accumulator_profile;
TaylorModelProductProfileContext g_taylor_model_product_profile_context=
    TaylorModelProductProfileContext::GENERAL;
TaylorModelProductProfileSnapshot g_taylor_model_product_profile;

Void accumulate_product_profile(
    TaylorModelProductProfileCounters& counters,
    unsigned long long product_pairs,
    unsigned long long individual_products_below_threshold,
    unsigned long long individual_products_below_threshold_collision,
    unsigned long long individual_products_below_threshold_new_term,
    unsigned long long individual_products_below_threshold_trailing,
    unsigned long long individual_products_above_threshold,
    double individual_products_below_threshold_abs_mass,
    double individual_products_above_threshold_abs_mass,
    unsigned long long sweep_passes,
    unsigned long long sweep_input_terms,
    unsigned long long sweep_output_terms,
    unsigned long long maximum_sweep_input_terms,
    unsigned long long maximum_sweep_output_terms,
    double elapsed_seconds)
{
    ++counters.calls;
    counters.product_pairs+=product_pairs;
    counters.individual_products_below_threshold+=
        individual_products_below_threshold;
    counters.individual_products_below_threshold_collision+=
        individual_products_below_threshold_collision;
    counters.individual_products_below_threshold_new_term+=
        individual_products_below_threshold_new_term;
    counters.individual_products_below_threshold_trailing+=
        individual_products_below_threshold_trailing;
    counters.individual_products_above_threshold+=
        individual_products_above_threshold;
    counters.individual_products_below_threshold_abs_mass+=
        individual_products_below_threshold_abs_mass;
    counters.individual_products_above_threshold_abs_mass+=
        individual_products_above_threshold_abs_mass;
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

Bool taylor_model_early_discard_enabled() {
    return g_taylor_model_early_discard_enabled;
}

Void set_taylor_model_early_discard_enabled(Bool enabled) {
    g_taylor_model_early_discard_enabled=enabled;
}

Bool taylor_model_incremental_sweep_enabled() {
    return g_taylor_model_incremental_sweep_enabled;
}

Void set_taylor_model_incremental_sweep_enabled(Bool enabled) {
    g_taylor_model_incremental_sweep_enabled=enabled;
}

Bool taylor_model_product_accumulator_enabled() {
    return g_taylor_model_product_accumulator_enabled;
}

Void set_taylor_model_product_accumulator_enabled(Bool enabled) {
    g_taylor_model_product_accumulator_enabled=enabled;
}

Bool taylor_model_dense_accumulator_enabled() {
    return g_taylor_model_dense_accumulator_enabled;
}

Void set_taylor_model_dense_accumulator_enabled(Bool enabled) {
    g_taylor_model_dense_accumulator_enabled=enabled;
}

Bool taylor_model_dense_batched_rounding_enabled() {
    return g_taylor_model_dense_batched_rounding_enabled;
}

Void set_taylor_model_dense_batched_rounding_enabled(Bool enabled) {
    g_taylor_model_dense_batched_rounding_enabled=enabled;
}

Void reset_taylor_model_dense_workspace_stats() {
    g_taylor_model_dense_workspace_stats=TaylorModelDenseWorkspaceStats();
}

TaylorModelDenseWorkspaceStats taylor_model_dense_workspace_stats() {
    return g_taylor_model_dense_workspace_stats;
}

Void reset_taylor_model_dense_hot_loop_profile() {
    g_taylor_model_dense_hot_loop_profile=TaylorModelDenseHotLoopProfile();
}

TaylorModelDenseHotLoopProfile taylor_model_dense_hot_loop_profile() {
    return g_taylor_model_dense_hot_loop_profile;
}

Void record_taylor_model_dense_hot_loop_profile(
    unsigned long long product_pairs,
    unsigned long long new_slots,
    unsigned long long collision_slots,
    double prepare_seconds,
    double prerank_seconds,
    double pair_loop_seconds,
    double emit_sweep_seconds)
{
    ++g_taylor_model_dense_hot_loop_profile.calls;
    g_taylor_model_dense_hot_loop_profile.product_pairs+=product_pairs;
    g_taylor_model_dense_hot_loop_profile.new_slots+=new_slots;
    g_taylor_model_dense_hot_loop_profile.collision_slots+=collision_slots;
    g_taylor_model_dense_hot_loop_profile.prepare_seconds+=prepare_seconds;
    g_taylor_model_dense_hot_loop_profile.prerank_seconds+=prerank_seconds;
    g_taylor_model_dense_hot_loop_profile.pair_loop_seconds+=pair_loop_seconds;
    g_taylor_model_dense_hot_loop_profile.emit_sweep_seconds+=emit_sweep_seconds;
}

Void record_taylor_model_dense_workspace_prepare(
    Bool slot_resize,
    Bool capacity_grow,
    unsigned long long slot_count)
{
    ++g_taylor_model_dense_workspace_stats.calls;
    if(slot_resize) {
        ++g_taylor_model_dense_workspace_stats.slot_resizes;
    }
    if(capacity_grow) {
        ++g_taylor_model_dense_workspace_stats.coefficient_capacity_grows;
    }
    g_taylor_model_dense_workspace_stats.maximum_slot_count=
        std::max(g_taylor_model_dense_workspace_stats.maximum_slot_count,
                 slot_count);
}

Void record_taylor_model_dense_workspace_touched(
    unsigned long long touched_count)
{
    g_taylor_model_dense_workspace_stats.maximum_touched_count=
        std::max(g_taylor_model_dense_workspace_stats.maximum_touched_count,
                 touched_count);
}

Void reset_taylor_model_accumulator_profile() {
    g_taylor_model_accumulator_profile=TaylorModelAccumulatorProfile();
}

TaylorModelAccumulatorProfile taylor_model_accumulator_profile() {
    return g_taylor_model_accumulator_profile;
}

Void record_taylor_model_accumulator_profile(
    unsigned long long product_pairs,
    unsigned long long temporary_entries,
    unsigned long long unique_entries,
    unsigned long long num_variables,
    unsigned long long x_degree,
    unsigned long long y_degree,
    unsigned long long product_degree,
    unsigned long long dense_slots)
{
    ++g_taylor_model_accumulator_profile.calls;
    g_taylor_model_accumulator_profile.product_pairs+=product_pairs;
    g_taylor_model_accumulator_profile.temporary_entries+=temporary_entries;
    g_taylor_model_accumulator_profile.unique_entries+=unique_entries;
    g_taylor_model_accumulator_profile.maximum_temporary_entries=
        std::max(g_taylor_model_accumulator_profile.maximum_temporary_entries,
                 temporary_entries);
    g_taylor_model_accumulator_profile.maximum_unique_entries=
        std::max(g_taylor_model_accumulator_profile.maximum_unique_entries,
                 unique_entries);
    g_taylor_model_accumulator_profile.maximum_argument_size=
        std::max(g_taylor_model_accumulator_profile.maximum_argument_size,
                 num_variables);
    g_taylor_model_accumulator_profile.maximum_x_degree=
        std::max(g_taylor_model_accumulator_profile.maximum_x_degree,x_degree);
    g_taylor_model_accumulator_profile.maximum_y_degree=
        std::max(g_taylor_model_accumulator_profile.maximum_y_degree,y_degree);
    g_taylor_model_accumulator_profile.maximum_product_degree=
        std::max(g_taylor_model_accumulator_profile.maximum_product_degree,
                 product_degree);
    g_taylor_model_accumulator_profile.maximum_dense_slots=
        std::max(g_taylor_model_accumulator_profile.maximum_dense_slots,
                 dense_slots);
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
    unsigned long long individual_products_below_threshold,
    unsigned long long individual_products_below_threshold_collision,
    unsigned long long individual_products_below_threshold_new_term,
    unsigned long long individual_products_below_threshold_trailing,
    unsigned long long individual_products_above_threshold,
    double individual_products_below_threshold_abs_mass,
    double individual_products_above_threshold_abs_mass,
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
        counters,product_pairs,
        individual_products_below_threshold,
        individual_products_below_threshold_collision,
        individual_products_below_threshold_new_term,
        individual_products_below_threshold_trailing,
        individual_products_above_threshold,
        individual_products_below_threshold_abs_mass,
        individual_products_above_threshold_abs_mass,
        sweep_passes,sweep_input_terms,sweep_output_terms,
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
