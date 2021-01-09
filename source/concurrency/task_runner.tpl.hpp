/***************************************************************************
 *            concurrency/task_runner.tpl.hpp
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

/*! \file concurrency/task_runner.tpl.hpp
 *  \brief Implementation code for runner classes.
 */

#include "../concurrency/task_runner.hpp"

namespace Ariadne {

template<class I, class O, class C>
SequentialRunnerBase<I,O,C>::SequentialRunnerBase(TaskSearchSpace const& space)
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

template<class I, class O, class C>
Void
DetachedRunnerBase<I,O,C>::_loop() {
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

template<class I, class O, class C>
DetachedRunnerBase<I,O,C>::DetachedRunnerBase(String const& thread_name, TaskSearchSpace const& space)
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

template<class I, class O>
class TaskIOData {
  public:
    TaskIOData(I const& input, O const& output, DurationType const& execution_time) : _input(input), _output(output), _execution_time(execution_time) { }
    TaskIOData& operator=(TaskIOData<I,O> const& p) {
        _input = p._input;
        _output = p._output;
        _execution_time = p._execution_time;
        return *this;
    };
    I const& input() const { return _input; }
    O const& output() const { return _output; }
    DurationType const& execution_time() const { return _execution_time; }
  private:
    I _input;
    O _output;
    DurationType _execution_time;
};

template<class I, class O>
class ParameterSearchOutputBufferData {
  public:
    ParameterSearchOutputBufferData(I const& input, O const& output, DurationType const& execution_time, TaskSearchPoint const& point) : _input(input), _output(output), _execution_time(execution_time), _point(point) { }
    ParameterSearchOutputBufferData& operator=(ParameterSearchOutputBufferData<I,O> const& p) {
        _input = p._input;
        _output = p._output;
        _execution_time = p._execution_time;
        _point = p._point;
        return *this;
    };
    I const& input() const { return _input; }
    O const& output() const { return _output; }
    DurationType const& execution_time() const { return _execution_time; }
    TaskSearchPoint const& point() const { return _point; }
  private:
    I _input;
    O _output;
    DurationType _execution_time;
    TaskSearchPoint _point;
};

template<class I, class O, class C>
Void
ParameterSearchRunnerBase<I,O,C>::_loop() {
    ARIADNE_LOG_SCOPE_CREATE;
    while(true) {
        std::unique_lock<std::mutex> locker(_input_mutex);
        _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
        locker.unlock();
        if (_terminate) break;
        auto pkg = _input_buffer.pop();
        auto cfg = to_configuration(pkg.first,pkg.second);
        auto start = std::chrono::high_resolution_clock::now();
        auto result = run_task(pkg.first,cfg);
        auto end = std::chrono::high_resolution_clock::now();
        auto execution_time = std::chrono::duration_cast<DurationType>(end-start);
        ARIADNE_LOG_PRINTLN("execution_time: " << execution_time.count() << " us");
        _output_buffer.push(OutputBufferContentType(pkg.first,result,execution_time,pkg.second));
        _output_availability.notify_all();
    }
}

template<class I, class O, class C>
ParameterSearchRunnerBase<I,O,C>::ParameterSearchRunnerBase(String const& thread_base_name, TaskSearchSpace const& space, Nat concurrency)
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

    Map<TaskSearchPoint,TaskIOData<I,O>> outputs;
    while (_output_buffer.size() > 0) {
        auto iodata = _output_buffer.pop();
        outputs.insert(Pair<TaskSearchPoint,TaskIOData<I,O>>(iodata.point(),TaskIOData<I,O>(iodata.input(),iodata.output(),
                                                                                            iodata.execution_time())));
    }
    auto appraisals = appraise(outputs);
    ARIADNE_LOG_PRINTLN("previous points:");
    for (auto a : appraisals) {
        ARIADNE_LOG_PRINTLN(a.point() << ", cost: " << a.cost());
    }

    Set<TaskSearchPoint> new_points;
    SizeType cnt = 0;
    for (auto a : appraisals) {
        new_points.insert(a.point());
        ++cnt;
        if (cnt >= _concurrency/2) break;
    }
    new_points = make_extended_set_by_shifting(new_points, _concurrency);
    for (auto p : new_points) _points.push(p);
    ARIADNE_LOG_PRINTLN("new points: " << new_points);

    auto best = appraisals.begin()->point();
    _best_points.push_back(best);
    return outputs.get(best).output();
}

} // namespace Ariadne
