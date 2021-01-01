/***************************************************************************
 *            concurrency/task_runner_interface.hpp
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

/*! \file concurrency/task_runner_interface.hpp
 *  \brief The interfaces for generic runners.
 */

#ifndef ARIADNE_TASK_RUNNER_INTERFACE_HPP
#define ARIADNE_TASK_RUNNER_INTERFACE_HPP

#include "../utility/tuple.hpp"
#include "../utility/pointer.hpp"
#include "../concurrency/loggable_smart_thread.hpp"
#include "../concurrency/buffer.hpp"
#include "../concurrency/task_parameter_point.hpp"
#include "../concurrency/task_parameter_space.hpp"

namespace Ariadne {

//! \brief Interface for the runner of a task.
//! \details This will usually be first implemented into an abstract base class.
template<class I, class O, class C>
class TaskRunnerInterface {
  public:
    typedef I InputType;
    typedef O OutputType;
    typedef C ConfigurationType;

    //! \brief Activates the runner, activating and threads and consequently properly setting the verbosity level
    virtual void activate() = 0;
    //! \brief Push input to the runner
    virtual void push(InputType const& input) = 0;
    //! \brief Pull output from the runner
    virtual OutputType pull() = 0;
    //! \brief Convert a task parameter point into a configuration of values for the task, possibly using \a in for values
    virtual ConfigurationType to_configuration(InputType const& in, TaskParameterPoint const& p) const = 0;
    //! \brief The task to be performed, taking \a in as input and \a cfg as a configuration of the parameters
    virtual OutputType run_task(InputType const& in, ConfigurationType const& cfg) const = 0;
};

//! \brief Interface for a class that supports a runnable task.
template<class I, class O, class C>
class TaskRunnableInterface {
  public:
    typedef I InputType;
    typedef O OutputType;
    typedef C ConfigurationType;

    //! \brief Set a new runner, useful to override the default runner
    virtual void set_runner(SharedPointer<TaskRunnerInterface<InputType,OutputType,ConfigurationType>> runner) = 0;
};

class TaskParameterSpace;

template<class I, class O, class C>
class SerialRunnerBase : public TaskRunnerInterface<I,O,C> {
  public:
    typedef typename TaskRunnerInterface<I,O,C>::InputType InputType;
    typedef typename TaskRunnerInterface<I,O,C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<I,O,C>::ConfigurationType ConfigurationType;

    SerialRunnerBase(TaskParameterSpace const& space);
    virtual ~SerialRunnerBase() = default;

    Void activate() override final;
    Void push(InputType const& input) override final;
    OutputType pull() override final;

    virtual ConfigurationType to_configuration(InputType const& in, TaskParameterPoint const& p) const override = 0;
    virtual OutputType run_task(InputType const& in, ConfigurationType const& cfg) const override = 0;

private:
    SharedPointer<OutputType> _last_output;
    // Parameter space
    SharedPointer<TaskParameterSpace> const _parameter_space;
};

template<class I, class O, class C>
SerialRunnerBase<I,O,C>::SerialRunnerBase(TaskParameterSpace const& space)
        :  _parameter_space(space.clone()) { }

template<class I, class O, class C>
Void
SerialRunnerBase<I,O,C>::activate()
{
    ARIADNE_LOG_SCOPE_CREATE;
}

    template<class I, class O, class C>
Void
SerialRunnerBase<I,O,C>::push(InputType const& input)
{
    _last_output.reset(new OutputType(run_task(input,to_configuration(input,_parameter_space->initial_point()))));
}

template<class RB, class I, class O>
typename SerialRunnerBase<RB,I,O>::OutputType
SerialRunnerBase<RB,I,O>::pull() {
    return *_last_output;
}

template<class I, class O, class C>
class ConcurrentRunnerBase : public TaskRunnerInterface<I,O,C> {
  public:
    typedef typename TaskRunnerInterface<I,O,C>::InputType InputType;
    typedef typename TaskRunnerInterface<I,O,C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<I,O,C>::ConfigurationType ConfigurationType;
    typedef Buffer<Pair<InputType,TaskParameterPoint>> InputBufferType;
    typedef Buffer<Pair<OutputType,TaskParameterPoint>> OutputBufferType;

    ConcurrentRunnerBase(String const& thread_name, TaskParameterSpace const& space);
    virtual ~ConcurrentRunnerBase();

    Void activate() override final;
    Void push(InputType const& input) override final;
    OutputType pull() override final;

    virtual ConfigurationType to_configuration(InputType const& in, TaskParameterPoint const& p) const override = 0;
    virtual OutputType run_task(InputType const& in, ConfigurationType const& cfg) const override = 0;

private:

    Void _loop() {
        ARIADNE_LOG_SCOPE_CREATE;
        while(true) {
            std::unique_lock<std::mutex> locker(_input_mutex);
            _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
            if (_terminate) break;
            auto pkg = _input_buffer.pop();
            _output_buffer.push({run_task(pkg.first,to_configuration(pkg.first,pkg.second)),pkg.second});
            _output_availability.notify_all();
        }
    }

private:
    // Parameter space
    SharedPointer<TaskParameterSpace> const _parameter_space;
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

template<class I, class O, class C>
ConcurrentRunnerBase<I,O,C>::ConcurrentRunnerBase(String const& thread_name, TaskParameterSpace const& space)
        :  _parameter_space(space.clone()), _thread(thread_name, [this]() { _loop(); }),
           _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
           _terminate(false) { }

template<class I, class O, class C>
ConcurrentRunnerBase<I,O,C>::~ConcurrentRunnerBase() {
    _terminate = true;
    _input_availability.notify_all();
}

template<class I, class O, class C>
Void
ConcurrentRunnerBase<I,O,C>::activate()
{
    _thread.activate();
}

template<class I, class O, class C>
Void
ConcurrentRunnerBase<I,O,C>::push(InputType const& input)
{
    _input_buffer.push({input,_parameter_space->initial_point()});
    _input_availability.notify_all();
}

template<class I, class O, class C>
typename ConcurrentRunnerBase<I,O,C>::OutputType
ConcurrentRunnerBase<I,O,C>::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>0; });
    auto result = _output_buffer.pop();
    return result.first;
}

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_INTERFACE_HPP
