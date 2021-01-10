/***************************************************************************
 *            concurrency/task_runner.hpp
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

/*! \file concurrency/task_runner.hpp
 *  \brief Runner classes.
 */

#ifndef ARIADNE_TASK_RUNNER_HPP
#define ARIADNE_TASK_RUNNER_HPP

#include "../concurrency/task_runner_interface.hpp"
#include "../concurrency/loggable_smart_thread.hpp"
#include "../concurrency/buffer.hpp"
#include "../concurrency/task_search_point.hpp"
#include "../concurrency/task_search_space.hpp"

namespace Ariadne {

template<class T> class TaskRunnerBase;

//! \brief Run a task sequentially.
//! \details Used to provide a sequential alternative to any thread-based implementation.
template<class T>
class SequentialRunner final : public TaskRunnerBase<T> {
  public:
    typedef typename TaskRunnerBase<T>::InputType InputType;
    typedef typename TaskRunnerBase<T>::OutputType OutputType;
    typedef typename TaskRunnerBase<T>::ConfigurationType ConfigurationType;

    Void activate() override final;
    Void push(InputType const& input) override final;
    OutputType pull() override final;

private:
    SharedPointer<OutputType> _last_output;
};

//! \brief Run a task in a detached thread, allowing concurrent processing.
template<class T>
class DetachedRunner final : public TaskRunnerBase<T> {
  public:
    typedef typename TaskRunnerBase<T>::InputType InputType;
    typedef typename TaskRunnerBase<T>::OutputType OutputType;
    typedef typename TaskRunnerBase<T>::ConfigurationType ConfigurationType;
    typedef Buffer<Pair<InputType,TaskSearchPoint>> InputBufferType;
    typedef Buffer<OutputType> OutputBufferType;

    DetachedRunner();
    virtual ~DetachedRunner();

    Void activate() override final;
    Void push(InputType const& input) override final;
    OutputType pull() override final;

private:
    Void _loop();
private:
    LoggableSmartThread _thread;
    InputBufferType _input_buffer;
    OutputBufferType _output_buffer;
    std::atomic<bool> _terminate;
    std::mutex _input_mutex;
    std::condition_variable _input_availability;
    std::mutex _output_mutex;
    std::condition_variable _output_availability;
};

typedef std::chrono::microseconds DurationType;

template<class I, class O> class TaskIOData;
template<class I, class O> class ParameterSearchOutputBufferData;

//! \brief Run a task by concurrent search into the parameter space.
template<class T>
class ParameterSearchRunner final : public TaskRunnerBase<T> {
public:
    typedef typename TaskRunnerBase<T>::InputType InputType;
    typedef typename TaskRunnerBase<T>::OutputType OutputType;
    typedef typename TaskRunnerBase<T>::ConfigurationType ConfigurationType;
    typedef Pair<InputType,TaskSearchPoint> InputBufferContentType;
    typedef ParameterSearchOutputBufferData<InputType,OutputType> OutputBufferContentType;
    typedef Buffer<InputBufferContentType> InputBufferType;
    typedef Buffer<OutputBufferContentType> OutputBufferType;

    ParameterSearchRunner(Nat concurrency);
    virtual ~ParameterSearchRunner();

    Void activate() override final;
    Void push(InputType const& input) override final;
    OutputType pull() override final;

private:
    Void _loop();
private:
    Nat const _concurrency;
    std::queue<TaskSearchPoint> _points;
    List<TaskSearchPoint> _best_points;
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

} // namespace Ariadne

#endif // ARIADNE_TASK_RUNNER_HPP
