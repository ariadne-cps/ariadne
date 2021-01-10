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

struct FlowStepTask final: public TaskInterface<FlowStepInput,FlowStepOutput,FlowStepConfiguration> {
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
        ARIADNE_LOG_PRINTLN_VAR(chosen_step_size);
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
    appraise(Map<TaskSearchPoint,TaskIOData<FlowStepInput,FlowStepOutput>> const& data) const override {
        Set<TaskSearchPointCost> result;

        Nat max_x = 0;
        for (auto entry : data) max_x = std::max(max_x,(Nat)entry.second.execution_time().count());
        for (auto entry : data) {
            CostType x = CostType(entry.second.execution_time().count())/max_x;
            CostType p = entry.second.output().step_size_used.get_d();
            result.insert(TaskSearchPointCost(entry.first, x/p));
        }
        return result;
    }

  private:
    TaskSearchSpace make_flow_step_runner_space() {
        RealVariable sssnr("starting_step_size_num_refinements"), st("sweep_threshold"), mto("maximum_temporal_order");
        return TaskSearchSpace({MetricSearchParameter(sssnr, 5, 2),
                                MetricSearchParameter(st, exp(-st * log(RealConstant(10))), 12, 9),
                                MetricSearchParameter(mto, 15, 12)
                               },(st*mto)/sssnr); }
    TaskSearchSpace const _space = make_flow_step_runner_space();
};

using FlowStepSequentialRunner = SequentialRunner<FlowStepTask>;
using FlowStepDetachedRunner = DetachedRunner<FlowStepTask>;
using FlowStepParameterSearchRunner = ParameterSearchRunner<FlowStepTask>;

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_EVOLVER_TASK_HPP
