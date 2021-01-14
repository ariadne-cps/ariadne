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
#include "../concurrency/task_appraisal_space.hpp"

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

class FlowStepTask final: public ParameterSearchTaskBase<FlowStepInput,FlowStepOutput,FlowStepConfiguration> {
    typedef FlowStepInput I;
    typedef FlowStepOutput O;
    typedef FlowStepConfiguration C;
  public:
    FlowStepTask() : ParameterSearchTaskBase<I,O,C>(
        "stp",
        TaskSearchSpace({
            MetricSearchParameter("starting_step_size_num_refinements", 2, 5),
            MetricSearchParameter("sweep_threshold", exp(-RealVariable("sweep_threshold") * log(RealConstant(10))), 8, 12),
            MetricSearchParameter("maximum_temporal_order", 9, 15)
        }),
        TaskAppraisalSpace<I,O>({
            execution_time_appraisal_parameter<I,O>,
            ScalarAppraisalParameter<I,O>("step_size_used",TaskAppraisalParameterOptimisation::MAXIMISE,[](I const& i,O const& o,DurationType const& d) { return o.step_size_used.get_d(); },2),
            VectorAppraisalParameter<I,O>("final_set_width_increases",TaskAppraisalParameterOptimisation::MINIMISE,
                                      [](I const& i,O const& o,DurationType const& d,SizeType const& idx) { return (o.evolve.euclidean_set().bounding_box()[idx].width() - i.current_set_bounds[idx].width()).get_d(); },
                                      [](I const& i){ return i.current_set_bounds.dimension(); })
        })) { }
  public:

    C to_configuration(I const& in, TaskSearchPoint const& p) const override {
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
        return C(integrator);
    }

    O run_task(I const& in, C const& cfg) const override {
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
        return O(next_set, reach_set, next_time, chosen_step_size);
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
