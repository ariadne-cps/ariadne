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
#include "../concurrency/task_interface.hpp"
#include "../concurrency/concurrency_manager.hpp"

namespace Ariadne {

template<class C>
class TaskRunnerBase : public TaskRunnerInterface<C> {
  public:
    typedef typename TaskRunnerInterface<C>::TaskType TaskType;
    typedef typename TaskRunnerInterface<C>::InputType InputType;
    typedef typename TaskRunnerInterface<C>::OutputType OutputType;
    typedef typename TaskRunnerInterface<C>::ConfigurationType ConfigurationType;

    TaskRunnerBase() : _task(new TaskType()) { }

    TaskType& task() override { return *_task; };
    TaskType const& task() const override { return *_task; };

    virtual ~TaskRunnerBase() = default;
  protected:
    SharedPointer<TaskType> const _task;
};

template<class C>
Void
SequentialRunner<C>::activate()
{
    ARIADNE_LOG_SCOPE_CREATE;
}

template<class C>
Void
SequentialRunner<C>::dump_statistics()
{ }

template<class C>
Void
SequentialRunner<C>::push(InputType const& input, ConfigurationType const& cfg)
{
    _last_output.reset(new OutputType(this->_task->run_task(input,this->_task->singleton_configuration(cfg,this->_task->search_space().initial_point()))));
}

template<class C>
typename SequentialRunner<C>::OutputType
SequentialRunner<C>::pull() {
    return *_last_output;
}

template<class R>
Void
DetachedRunner<R>::_loop() {
    ARIADNE_LOG_SCOPE_CREATE;
    while(true) {
        std::unique_lock<std::mutex> locker(_input_mutex);
        _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
        if (_terminate) break;
        auto pkg = _input_buffer.pop();
        _output_buffer.push(this->_task->run_task(pkg.first,this->_task->singleton_configuration(pkg.second.first,pkg.second.second)));
        _output_availability.notify_all();
    }
}

template<class T>
DetachedRunner<T>::DetachedRunner()
        :  _thread(this->_task->name(), [this]() { _loop(); }),
           _input_buffer(InputBufferType(1)),_output_buffer(OutputBufferType(1)),
           _terminate(false) { }

template<class T>
DetachedRunner<T>::~DetachedRunner() {
    _terminate = true;
    _input_availability.notify_all();
}

template<class T>
Void
DetachedRunner<T>::activate()
{
    _thread.activate();
}

template<class T>
Void
DetachedRunner<T>::dump_statistics() {
}

template<class T>
Void
DetachedRunner<T>::push(InputType const& input, ConfigurationType const& cfg)
{
    _input_buffer.push({input,{cfg,this->_task->search_space().initial_point()}});
    _input_availability.notify_all();
}

template<class T>
typename DetachedRunner<T>::OutputType
DetachedRunner<T>::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>0; });
    return _output_buffer.pop();
}

template<class O>
class ParameterSearchOutputBufferData {
  public:
    ParameterSearchOutputBufferData(O const& output, DurationType const& execution_time, TaskSearchPoint const& point) : _output(output), _execution_time(execution_time), _point(point) { }
    ParameterSearchOutputBufferData(ParameterSearchOutputBufferData<O> const& p) : _output(p._output), _execution_time(p._execution_time), _point(p._point) { }
    ParameterSearchOutputBufferData& operator=(ParameterSearchOutputBufferData<O> const& p) {
        _output = p._output;
        _execution_time = p._execution_time;
        _point = p._point;
        return *this;
    };
    O const& output() const { return _output; }
    DurationType const& execution_time() const { return _execution_time; }
    TaskSearchPoint const& point() const { return _point; }
  private:
    O _output;
    DurationType _execution_time;
    TaskSearchPoint _point;
};

template<class T>
Void
ParameterSearchRunner<T>::_loop() {
    ARIADNE_LOG_SCOPE_CREATE;
    while(true) {
        std::unique_lock<std::mutex> locker(_input_mutex);
        _input_availability.wait(locker, [this]() { return _input_buffer.size()>0 || _terminate; });
        locker.unlock();
        if (_terminate) break;
        auto pkg = _input_buffer.pop();
        auto cfg = this->_task->singleton_configuration(pkg.second.first,pkg.second.second);
        auto start = std::chrono::high_resolution_clock::now();
        try {
            auto result = this->_task->run_task(pkg.first,cfg);
            auto end = std::chrono::high_resolution_clock::now();
            auto execution_time = std::chrono::duration_cast<DurationType>(end-start);
            ARIADNE_LOG_PRINTLN("task for " << pkg.second << " completed in " << execution_time.count() << " us");
            _output_buffer.push(OutputBufferContentType(result,execution_time,pkg.second.second));
        } catch (std::exception& e) {
            ++_failures;
            ARIADNE_LOG_PRINTLN("task failed: " << e.what());
        }
        _output_availability.notify_all();
    }
}

template<class T>
ParameterSearchRunner<T>::ParameterSearchRunner(Nat concurrency)
        :  _concurrency(concurrency), _failures(0), _last_used_input(1),
           _input_buffer(InputBufferType(concurrency)),_output_buffer(OutputBufferType(concurrency)),
           _terminate(false) {
    for (Nat i=0; i<concurrency; ++i)
        _threads.append(SharedPointer<LoggableSmartThread>(new LoggableSmartThread(this->_task->name() + (concurrency>=10 and i<10 ? "0" : "") + to_string(i), [this]() { _loop(); })));

    auto initial = this->_task->search_space().initial_point();
    _points.push(initial);
    if (_concurrency>1) {
        auto shifted = initial.make_random_shifted(_concurrency-1);
        for (auto point : shifted)
            _points.push(point);
    }
}

template<class T>
ParameterSearchRunner<T>::~ParameterSearchRunner() {
    _terminate = true;
    _input_availability.notify_all();
}

template<class T>
Void
ParameterSearchRunner<T>::activate() {
    for (auto thread : _threads)
        thread->activate();
}

template<class T>
Void
ParameterSearchRunner<T>::dump_statistics() {
    ConcurrencyManager::instance().set_last_search_best_points(_best_points);
    _best_points.clear();
}

template<class T>
Void
ParameterSearchRunner<T>::push(InputType const& input, ConfigurationType const& cfg) {
    for (SizeType i=0; i<_concurrency; ++i) {
        _input_buffer.push({input,{cfg,_points.front()}});
        _points.pop();
    }
    _last_used_input.push(input);
    _input_availability.notify_all();
}

template<class T>
typename ParameterSearchRunner<T>::OutputType
ParameterSearchRunner<  T>::pull() {
    std::unique_lock<std::mutex> locker(_output_mutex);
    _output_availability.wait(locker, [this]() { return _output_buffer.size()>=_concurrency-_failures; });
    ARIADNE_LOG_PRINTLN("received " << _concurrency-_failures << " completed tasks");
    _failures=0;

    InputType input = _last_used_input.pop();
    Map<TaskSearchPoint,Pair<OutputType,DurationType>> outputs;
    while (_output_buffer.size() > 0) {
        auto io_data = _output_buffer.pop();
        outputs.insert(Pair<TaskSearchPoint,Pair<OutputType,DurationType>>(
                io_data.point(),{io_data.output(),io_data.execution_time()}));
    }
    auto appraisals = this->_task->appraise(outputs,input);
    ARIADNE_LOG_PRINTLN_VAR(appraisals);

    Set<TaskSearchPoint> new_points;
    SizeType cnt = 0;
    for (auto a : appraisals) {
        new_points.insert(a.point());
        ++cnt;
        if (cnt >= _concurrency/2) break;
    }
    new_points = make_extended_set_by_shifting(new_points, _concurrency);
    for (auto p : new_points) _points.push(p);
    ARIADNE_LOG_PRINTLN_VAR(new_points);

    auto best = appraisals.begin()->point();
    if (appraisals.begin()->critical_failures() > 0)
        throw CriticalAppraisalFailureException(appraisals);
    _best_points.push_back(*appraisals.begin());
    return outputs.get(best).first;
}

} // namespace Ariadne
