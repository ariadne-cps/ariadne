/***************************************************************************
 *            dynamics/vector_field_evolver_task.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/*! \file dynamics/vector_field_evolver_task.hpp
 *  \brief Tasks specification for vector field evolver.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_TASK_HPP
#define ARIADNE_VECTOR_FIELD_EVOLVER_TASK_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>


#include "../utility/tuple.hpp"

#include "../function/function_interface.hpp"
#include "../solvers/configuration_interface.hpp"
#include "../solvers/integrator_interface.hpp"
#include "../solvers/integrator.hpp"

#include "../concurrency/task_runner.tpl.hpp"

namespace Ariadne {

using std::min, std::max;

struct FlowStepInput {
    FlowStepInput(EffectiveVectorMultivariateFunction const& dynamic_, IntegratorInterface const& integrator_, LabelledEnclosure const& current_set_,
                  FloatDPExactBox const& current_set_bounds_, Dyadic const& current_time_, Dyadic const& previous_step_size_, Dyadic const& maximum_step_size_) :
            dynamic(dynamic_), integrator(integrator_), current_set(current_set_), current_set_bounds(current_set_bounds_),
            current_time(current_time_), previous_step_size(previous_step_size_), maximum_step_size(maximum_step_size_) { }
    EffectiveVectorMultivariateFunction const& dynamic;
    IntegratorInterface const& integrator;
    LabelledEnclosure const& current_set;
    FloatDPExactBox const& current_set_bounds;
    Dyadic const& current_time;
    Dyadic const& previous_step_size;
    Dyadic const& maximum_step_size;
};

struct FlowStepOutput {
    FlowStepOutput(LabelledEnclosure const& evolve_, LabelledEnclosure const& reach_, Dyadic const& time_, Dyadic const& step_size_used_) :
            evolve(evolve_), reach(reach_), time(time_), step_size_used(step_size_used_) { }
    LabelledEnclosure const evolve;
    LabelledEnclosure const reach;
    Dyadic const time;
    Dyadic const step_size_used;
};

struct FlowStepConfiguration {
    FlowStepConfiguration(SharedPointer<TaylorPicardIntegrator> const& integrator_) : integrator(integrator_){ }
    SharedPointer<TaylorPicardIntegrator> integrator;
};

class FlowStepTask final: public TaskInterface<FlowStepInput,FlowStepOutput,FlowStepConfiguration> {
  private:
    TaskSearchSpace const _space = TaskSearchSpace(
            {MetricSearchParameter("starting_step_size_num_refinements", 0, 5, 2),
             MetricSearchParameter("sweep_threshold", exp(-RealVariable("sweep_threshold") * log(RealConstant(10))), 6, 12, 9),
             MetricSearchParameter("maximum_temporal_order", 2, 15, 12)
            });
  public:
    std::string name() const override { return "stp"; }
    TaskSearchSpace const& search_space() const override { return _space; }

    FlowStepConfiguration
    to_configuration(FlowStepInput const& in, TaskSearchPoint const& p) const override {
        TaylorPicardIntegrator const& default_integrator = static_cast<TaylorPicardIntegrator const&>(in.integrator);
        SharedPointer<TaylorPicardIntegrator> integrator(new TaylorPicardIntegrator(
                MaximumError(default_integrator.maximum_error()),
                ThresholdSweeper<FloatDP>(DoublePrecision(),p.value("sweep_threshold")),
                LipschitzConstant(default_integrator.lipschitz_tolerance()),
                StartingStepSizeNumRefinements(p.value("starting_step_size_num_refinements").get_d()),
                StepMaximumError(default_integrator.step_maximum_error()),
                MinimumTemporalOrder(default_integrator.minimum_temporal_order()),
                MaximumTemporalOrder(p.value("maximum_temporal_order").get_d())
        ));
        return FlowStepConfiguration(integrator);
    }

    FlowStepOutput
    run_task(FlowStepInput const& in, FlowStepConfiguration const& cfg) const override {
        LabelledEnclosure next_set = in.current_set;
        LabelledEnclosure reach_set = in.current_set;
        Dyadic next_time = in.current_time;
        Dyadic chosen_step_size = in.maximum_step_size;
        FlowStepModelType flow_model = cfg.integrator->flow_step(in.dynamic, in.current_set_bounds, in.previous_step_size,chosen_step_size);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, chosen_step_size);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, flow_model);
        next_time += chosen_step_size;
        ARIADNE_LOG_PRINTLN_VAR_AT(1, next_time);
        reach_set.apply_full_reach_step(flow_model);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, reach_set);
        next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
        ARIADNE_LOG_PRINTLN_VAR_AT(1, next_set);
        return FlowStepOutput(next_set, reach_set, next_time, chosen_step_size);
    }

    Set<TaskSearchPointCost>
    appraise(Map<TaskSearchPoint,Pair<FlowStepOutput,DurationType>> const& data, FlowStepInput const& input) const override {
        Set<TaskSearchPointCost> result;

        auto dim = input.current_set_bounds.dimension();
        auto initial_widths = input.current_set_bounds.widths();

        auto data_iter = data.cbegin();

        // Get the minimum and maximum execution time
        Nat min_x((Nat)data_iter->second.second.count()), max_x((Nat)data_iter->second.second.count());
        // Get the minimum and maximum step size
        StepSizeType min_step_size(data_iter->second.first.step_size_used), max_step_size(data_iter->second.first.step_size_used);
        // Get the minimum and maximum differences of widths
        Vector<CostType> min_width_diffs(dim), max_width_diffs(dim);
        auto final_widths = data_iter->second.first.evolve.bounding_box().continuous_set().widths();
        for (SizeType i=0; i<dim; ++i) {
            min_width_diffs[i] = (final_widths[i]-initial_widths[i]).get_d();
            max_width_diffs[i] = (final_widths[i]-initial_widths[i]).get_d();
        }

        ++data_iter;
        while (data_iter != data.cend()) {
            min_x = min(min_x,(Nat)data_iter->second.second.count());
            max_x = max(max_x,(Nat)data_iter->second.second.count());
            min_step_size = min(min_step_size,data_iter->second.first.step_size_used);
            max_step_size = max(max_step_size,data_iter->second.first.step_size_used);
            final_widths = data_iter->second.first.evolve.bounding_box().continuous_set().widths();
            for (SizeType i=0; i<dim; ++i) {
                min_width_diffs[i] = min(min_width_diffs[i],(final_widths[i]-initial_widths[i]).get_d());
                max_width_diffs[i] = max(max_width_diffs[i],(final_widths[i]-initial_widths[i]).get_d());
            }
            ++data_iter;
        }

        for (auto entry : data) {
            final_widths = entry.second.first.evolve.bounding_box().continuous_set().widths();
            CostType a = 0;
            for (SizeType i=0; i<dim; ++i) {
                a += (max_width_diffs[i] != min_width_diffs[i] ? CostType(((final_widths[i]-initial_widths[i]-min_width_diffs[i])/(max_width_diffs[i]-min_width_diffs[i])).get_d()) : 0);
            }
            a/=dim;
            auto x = (max_x != min_x ? (CostType(entry.second.second.count())-min_x)/(max_x-min_x) : 0);
            auto p = (max_step_size != min_step_size ? ((entry.second.first.step_size_used-min_step_size)/(max_step_size-min_step_size)).get_d() : 0.0);
            result.insert(TaskSearchPointCost(entry.first, 1*a+1*x-2*p));
        }
        return result;
    }
};

template class SequentialRunner<FlowStepTask>;
template class DetachedRunner<FlowStepTask>;
template class ParameterSearchRunner<FlowStepTask>;

using FlowStepSequentialRunner = SequentialRunner<FlowStepTask>;
using FlowStepDetachedRunner = DetachedRunner<FlowStepTask>;
using FlowStepParameterSearchRunner = ParameterSearchRunner<FlowStepTask>;

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_TASK_HPP
