/***************************************************************************
 *            concurrency/parameter_exploration_runner.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file concurrency/parameter_exploration_runner.hpp
 *  \brief A runner for concurrent exploration of parameters.
 */

#ifndef ARIADNE_PARAMETER_EXPLORATION_RUNNER_HPP
#define ARIADNE_PARAMETER_EXPLORATION_RUNNER_HPP

#include <utility>
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#include <string>

#include "../concurrency/loggable_smart_thread.hpp"
#include "../concurrency/buffer.hpp"
#include "../concurrency/task_parameter_point.hpp"
#include "../concurrency/task_parameter_space.hpp"

#include "../dynamics/vector_field_evolver.hpp"

namespace Ariadne {

class FlowStepInput {
  public:
    FlowStepInput(LabelledEnclosure const& current_set, FloatDPExactBox const& current_set_bounds, Dyadic const& current_time, Dyadic const& previous_step_size) :
            _current_set(current_set), _current_set_bounds(current_set_bounds), _current_time(current_time), _previous_step_size(previous_step_size) { }
    LabelledEnclosure const& _current_set;
    FloatDPExactBox const& _current_set_bounds;
    Dyadic const& _current_time;
    Dyadic const& _previous_step_size;
};

struct FlowStepConfiguration {
    FlowStepConfiguration(SharedPointer<TaylorPicardIntegrator> const& integrator) : integrator(integrator){ }
    SharedPointer<TaylorPicardIntegrator> integrator;
};

class FlowStepOutput {
  public:
    FlowStepOutput(LabelledEnclosure const& evolve, LabelledEnclosure const& reach, Dyadic const& time, Dyadic const& step_size_used) :
            _evolve(evolve), _reach(reach), _time(time), _step_size_used(step_size_used) { }
    LabelledEnclosure const _evolve;
    LabelledEnclosure const _reach;
    Dyadic const _time;
    Dyadic const _step_size_used;
};

class ParameterExplorationRunner {
    typedef Buffer<Pair<FlowStepInput,TaskParameterPoint>> InputBufferType;
    typedef Buffer<Pair<FlowStepOutput,TaskParameterPoint>> OutputBufferType;
  public:
    ParameterExplorationRunner(TaskParameterSpace const& space, VectorFieldEvolver const& evolver);
    ParameterExplorationRunner(TaskParameterSpace const& space, EffectiveVectorMultivariateFunction const& dynamic, TaylorPicardIntegrator const& integrator, Dyadic const& maximum_step_size);

    ~ParameterExplorationRunner();

    //! \brief Activates the runner, activating and threads and consequently properly setting the verbosity level
    Void activate();

    Void push(FlowStepInput const& input);

    FlowStepOutput pull();

  private:

    Void _loop() {
        ARIADNE_LOG_SCOPE_CREATE;
        while(true) {
            std::unique_lock<std::mutex> locker(_input_mutex);
            _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
            if (_terminate) break;
            auto pkg = _input_buffer.pop();
            _output_buffer.push({_task(pkg.first,_to_configuration(pkg.second)),pkg.second});
            _output_availability.notify_all();
        }
    }

    FlowStepConfiguration _to_configuration(TaskParameterPoint const& p) {
        SharedPointer<TaylorPicardIntegrator> integrator(new TaylorPicardIntegrator(
                MaximumError(_integrator->maximum_error()),
                ThresholdSweeper<FloatDP>(DoublePrecision(),p.value("sweep_threshold")),
                LipschitzConstant(_integrator->lipschitz_tolerance()),
                StartingStepSizeNumRefinements(p.value("starting_step_size_num_refinements").get_d()),
                StepMaximumError(_integrator->step_maximum_error()),
                MinimumTemporalOrder(_integrator->minimum_temporal_order()),
                MaximumTemporalOrder(p.value("maximum_temporal_order").get_d())
        ));
        return FlowStepConfiguration(integrator);
    }

    FlowStepOutput _task(FlowStepInput const& in, FlowStepConfiguration const& cfg) {

        LabelledEnclosure next_set = in._current_set;
        LabelledEnclosure reach_set = in._current_set;
        Dyadic next_time = in._current_time;
        Dyadic chosen_step_size = _maximum_step_size;
        FlowStepModelType flow_model = cfg.integrator->flow_step(_dynamic, in._current_set_bounds, in._previous_step_size,chosen_step_size);
        ARIADNE_LOG_PRINTLN("step_size = " << chosen_step_size);
        ARIADNE_LOG_PRINTLN_AT(1, "flow_model = " << flow_model);
        next_time += chosen_step_size;
        ARIADNE_LOG_PRINTLN_AT(1, "next_time = " << next_time)
        reach_set.apply_full_reach_step(flow_model);
        ARIADNE_LOG_PRINTLN_AT(1, "reach_set = " << reach_set);
        next_set.apply_fixed_evolve_step(flow_model, chosen_step_size);
        ARIADNE_LOG_PRINTLN_AT(1, "next_set = " << next_set);

        return FlowStepOutput(next_set,reach_set,next_time,chosen_step_size);
    }

  private:
    // Parameter space
    TaskParameterSpace const& _parameter_space;
    // Constants
    EffectiveVectorMultivariateFunction const _dynamic;
    SharedPointer<TaylorPicardIntegrator> const _integrator;
    Dyadic const _maximum_step_size;
    // Synchronization
    LoggableSmartThread _thread;
    InputBufferType _input_buffer;
    OutputBufferType _output_buffer;
    std::atomic<bool> _terminate;
    std::mutex _input_mutex;
    std::condition_variable _input_availability;
    std::mutex _output_mutex;
    std::condition_variable _output_availability;
};

ParameterExplorationRunner::ParameterExplorationRunner(TaskParameterSpace const& space, VectorFieldEvolver const& evolver)
        :  _parameter_space(space), _dynamic(evolver.system().dynamic_function()), _integrator(static_cast<const TaylorPicardIntegrator*>(evolver.integrator())->clone()), _maximum_step_size(evolver.configuration().maximum_step_size()), _thread("step", [this]() { _loop(); }),
           _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
           _terminate(false) { }

ParameterExplorationRunner::ParameterExplorationRunner(TaskParameterSpace const& space, EffectiveVectorMultivariateFunction const& dynamic, TaylorPicardIntegrator const& integrator, Dyadic const& maximum_step_size)
        :  _parameter_space(space), _dynamic(dynamic), _integrator(integrator.clone()), _maximum_step_size(maximum_step_size), _thread("step", [this]() { _loop(); }),
           _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
           _terminate(false) { }

ParameterExplorationRunner::~ParameterExplorationRunner() {
    _terminate = true;
    _input_availability.notify_all();
}

Void
ParameterExplorationRunner::activate()
{
    _thread.activate();
}

Void
ParameterExplorationRunner::push(FlowStepInput const& input)
{
    _input_buffer.push({input,_parameter_space.initial_point()});
    _input_availability.notify_all();
}

FlowStepOutput
ParameterExplorationRunner::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>0; });
    auto result = _output_buffer.pop();
    return result.first;
}

} // namespace Ariadne

#endif // ARIADNE_PARAMETER_EXPLORATION_RUNNER_HPP
