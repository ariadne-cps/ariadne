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

    //! \brief Activates the runner, activating the threads and consequently properly setting the verbosity level
    //! \details This is necessary since the object may need to be visible at a higher level of abstraction compared
    //! to the use. Multiple activations are possible due to multiple calls of the runner, though they will not have any effect,
    //! since a thread can be activated only once, with no exception thrown.
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

//! \brief Run a task sequentially.
//! \details Used to provide a sequential alternative to any thread-based implementation.
template<class I, class O, class C>
class SequentialRunnerBase : public TaskRunnerInterface<I,O,C> {
  public:
    typedef typename TaskRunnerInterface<I,O,C>::InputType InputType;
    typedef typename TaskRunnerInterface<I,O,C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<I,O,C>::ConfigurationType ConfigurationType;

    SequentialRunnerBase(TaskParameterSpace const& space);
    virtual ~SequentialRunnerBase() = default;

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
SequentialRunnerBase<I,O,C>::SequentialRunnerBase(TaskParameterSpace const& space)
        :  _parameter_space(space.clone()) { }

template<class I, class O, class C>
Void
SequentialRunnerBase<I,O,C>::activate()
{
    ARIADNE_LOG_SCOPE_CREATE;
}

    template<class I, class O, class C>
Void
SequentialRunnerBase<I,O,C>::push(InputType const& input)
{
    _last_output.reset(new OutputType(run_task(input,to_configuration(input,_parameter_space->initial_point()))));
}

template<class RB, class I, class O>
typename SequentialRunnerBase<RB,I,O>::OutputType
SequentialRunnerBase<RB,I,O>::pull() {
    return *_last_output;
}

//! \brief Run a task in a detached thread, allowing concurrent processing.
template<class I, class O, class C>
class DetachedRunnerBase : public TaskRunnerInterface<I,O,C> {
  public:
    typedef typename TaskRunnerInterface<I,O,C>::InputType InputType;
    typedef typename TaskRunnerInterface<I,O,C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<I,O,C>::ConfigurationType ConfigurationType;
    typedef Buffer<Pair<InputType,TaskParameterPoint>> InputBufferType;
    typedef Buffer<OutputType> OutputBufferType;

    DetachedRunnerBase(String const& thread_name, TaskParameterSpace const& space);
    virtual ~DetachedRunnerBase();

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
            _output_buffer.push(run_task(pkg.first,to_configuration(pkg.first,pkg.second)));
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
DetachedRunnerBase<I,O,C>::DetachedRunnerBase(String const& thread_name, TaskParameterSpace const& space)
        :  _parameter_space(space.clone()), _thread(thread_name, [this]() { _loop(); }),
           _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
           _terminate(false) { }

template<class I, class O, class C>
DetachedRunnerBase<I,O,C>::~DetachedRunnerBase() {
    _terminate = true;
    _input_availability.notify_all();
}

template<class I, class O, class C>
Void
DetachedRunnerBase<I,O,C>::activate()
{
    _thread.activate();
}

template<class I, class O, class C>
Void
DetachedRunnerBase<I,O,C>::push(InputType const& input)
{
    _input_buffer.push({input,_parameter_space->initial_point()});
    _input_availability.notify_all();
}

template<class I, class O, class C>
typename DetachedRunnerBase<I,O,C>::OutputType
DetachedRunnerBase<I,O,C>::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>0; });
    return _output_buffer.pop();
}

typedef std::chrono::microseconds DurationType;

template<class D>
class ParameterSearchOutput {
  public:
    ParameterSearchOutput(D const& data, DurationType const& duration, TaskParameterPoint const& point) : _data(data), _point(point), _duration(duration) { }
    ParameterSearchOutput& operator=(ParameterSearchOutput<D> const& p) {
        _data = p._data;
        _point = p._point;
        _duration = p._duration;
        return *this;
    };
    D const& data() const { return _data; }
    TaskParameterPoint const& point() const { return _point; }
    DurationType const& duration() const { return _duration; }
  private:
    D _data;
    TaskParameterPoint _point;
    DurationType _duration;
};

//! \brief Run a task by concurrent search into the parameter space.
template<class I, class O, class C>
class ParameterSearchRunnerBase : public TaskRunnerInterface<I,O,C> {
public:
    typedef typename TaskRunnerInterface<I,O,C>::InputType InputType;
    typedef typename TaskRunnerInterface<I,O,C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<I,O,C>::ConfigurationType ConfigurationType;
    typedef Pair<InputType,TaskParameterPoint> InputBufferContentType;
    typedef ParameterSearchOutput<O> OutputBufferContentType;
    typedef Buffer<InputBufferContentType> InputBufferType;
    typedef Buffer<OutputBufferContentType> OutputBufferType;

    ParameterSearchRunnerBase(String const& thread_base_name, TaskParameterSpace const& space, Nat concurrency);
    virtual ~ParameterSearchRunnerBase();

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
            auto cfg = to_configuration(pkg.first,pkg.second);
            auto start = std::chrono::high_resolution_clock::now();
            auto result = run_task(pkg.first,cfg);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<DurationType>(end-start);
            ARIADNE_LOG_PRINTLN("duration: " << duration.count() << " us");
            _output_buffer.push(OutputBufferContentType(result,duration,pkg.second));
            _output_availability.notify_all();
        }
    }

private:
    // Parameter space
    SharedPointer<TaskParameterSpace> const _parameter_space;
    // Concurrency
    Nat const _concurrency;
    std::queue<TaskParameterPoint> _points;
    List<TaskParameterPoint> _best_points;
    // Synchronization
    List<SharedPointer<LoggableSmartThread>> _threads;
    InputBufferType _input_buffer;
    OutputBufferType _output_buffer;
    std::atomic<bool> _terminate;
    std::mutex _input_mutex;
    std::condition_variable _input_availability;
    std::mutex _output_mutex;
    std::condition_variable _output_availability;
};

template<class I, class O, class C>
ParameterSearchRunnerBase<I,O,C>::ParameterSearchRunnerBase(String const& thread_base_name, TaskParameterSpace const& space, Nat concurrency)
        :  _parameter_space(space.clone()), _concurrency(concurrency),
           _input_buffer(InputBufferType(concurrency)),_output_buffer(OutputBufferType(concurrency)),
           _terminate(false)
{
    for (Nat i=0; i<concurrency; ++i)
        _threads.append(SharedPointer<LoggableSmartThread>(new LoggableSmartThread(thread_base_name + (concurrency>=10 and i<10 ? "0" : "") + to_string(i), [this]() { _loop(); })));

    auto initial = _parameter_space->initial_point();
    _points.push(initial);
    if (_concurrency>1) {
        auto shifted = initial.make_random_shifted(_concurrency-1);
        for (auto point : shifted)
            _points.push(point);
    }
}

template<class I, class O, class C>
ParameterSearchRunnerBase<I,O,C>::~ParameterSearchRunnerBase() {
    ARIADNE_LOG_PRINTLN("Best points: " << _best_points);
    _terminate = true;
    _input_availability.notify_all();
}

template<class I, class O, class C>
Void
ParameterSearchRunnerBase<I,O,C>::activate()
{
    for (auto thread : _threads)
        thread->activate();
}

template<class I, class O, class C>
Void
ParameterSearchRunnerBase<I,O,C>::push(InputType const& input)
{
    for (SizeType i=0; i<_concurrency; ++i) {
        _input_buffer.push({input,_points.front()});
        _points.pop();
    }
    _input_availability.notify_all();
}

template<class I, class O, class C>
typename ParameterSearchRunnerBase<I,O,C>::OutputType
ParameterSearchRunnerBase<I,O,C>::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()==_concurrency; });

    List<OutputBufferContentType> result;
    result.push_back(_output_buffer.pop());
    while (_output_buffer.size() > 0) {
        auto current = _output_buffer.pop();
        auto execution_time = current.duration();
        if (execution_time < result.back().duration()) {
            result.pop_back();
            result.push_back(current);
        }
    }
    _best_points.push_back(result.back().point());

    _points.push(result.back().point());
    if (_concurrency>1) {
        auto shifted = result.back().point().make_adjacent_shifted(_concurrency-1);
        for (auto point : shifted)
            _points.push(point);
    }

    return result.back().data();
}

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_INTERFACE_HPP
