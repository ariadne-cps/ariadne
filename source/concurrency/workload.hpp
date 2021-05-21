/***************************************************************************
 *            concurrency/workload.hpp
 *
 *  Copyright  2007-21  Luca Geretti
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

/*! \file concurrency/workload.hpp
 *  \brief A stack of objects to work on, supplied with a function to process them
 */

#ifndef ARIADNE_WORKLOAD_HPP
#define ARIADNE_WORKLOAD_HPP

#include <functional>
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "io/progress_indicator.hpp"
#include "concurrency_typedefs.hpp"
#include "workload_interface.hpp"
#include "task_manager.hpp"
#include "workload_advancement.hpp"

namespace Ariadne {

//! \brief Base class implementation
template<class E, class... AS>
class WorkloadBase : public WorkloadInterface<E,AS...> {
  protected:
    WorkloadBase() : _progress_acknowledge_func(std::bind_front(&WorkloadBase::_default_progress_acknowledge, this)), _advancement(0), _logger_level(0), _progress_indicator(new ProgressIndicator(0)) { }
  public:
    using TaskFunctionType = std::function<Void(E const &)>;
    using ProgressAcknowledgeFunctionType = std::function<Void(E const &, SharedPointer<ProgressIndicator>)>;
    using CompletelyBoundFunctionType = VoidFunction;

    Void process() override {
        _log_scope_manager.reset(new LogScopeManager(ARIADNE_PRETTY_FUNCTION,0));
        _logger_level = Logger::instance().current_level();
        while (true) {
            UniqueLock<Mutex> lock(_element_availability_mutex);
            _element_availability_condition.wait(lock, [=,this] { return _advancement.has_finished() or not _sequential_queue.empty() or _exception != nullptr; });
            if (_exception != nullptr) rethrow_exception(_exception);
            if (_advancement.has_finished()) { _log_scope_manager.reset(); return; }

            CompletelyBoundFunctionType task, progress_acknowledge;
            make_lpair(task,progress_acknowledge) = _sequential_queue.front();
            _sequential_queue.pop();
            if (_using_concurrency()) {
                TaskManager::instance().enqueue([this, task, progress_acknowledge] { _concurrent_task_wrapper(task,progress_acknowledge); });
            } else {
                _advancement.add_to_processing();
                if (not Logger::instance().is_muted_at(0)) { progress_acknowledge(); _print_hold(); }
                task();
                _advancement.add_to_completed();
            }
        }
    }

    SizeType size() const override { return _sequential_queue.size(); }

    WorkloadInterface<E,AS...>& append(E const& e) override {
        _advancement.add_to_waiting();
        _sequential_queue.push(std::make_pair(std::bind(std::forward<TaskFunctionType const>(_task_func), std::forward<E const&>(e)),
                                                   std::bind(std::forward<ProgressAcknowledgeFunctionType const>(_progress_acknowledge_func), std::forward<E const&>(e),
                                                  _progress_indicator)
                         ));
        return *this;
    }

    WorkloadInterface<E,AS...>& append(List<E> const& es) override { for (auto e : es) append(e); return *this; }

  private:

    Bool _using_concurrency() const { return TaskManager::instance().concurrency() > 0; }

    Void _concurrent_task_wrapper(CompletelyBoundFunctionType const& task, CompletelyBoundFunctionType const& progress_acknowledge) {
        _advancement.add_to_processing();
        Logger::instance().set_level(_logger_level);
        if (not Logger::instance().is_muted_at(0)) { progress_acknowledge(); _print_hold(); }
        try {
            task();
        } catch (...) {
            {
                LockGuard<Mutex> lock(_element_availability_mutex);
                _exception = std::current_exception();
            }
            _element_availability_condition.notify_one();
        }

        {
            LockGuard<Mutex> lock(_element_availability_mutex);
            _advancement.add_to_completed();
        }
        if (_advancement.has_finished()) { _element_availability_condition.notify_one(); }
    }

    void _default_progress_acknowledge(E const& e, SharedPointer<ProgressIndicator> indicator) {
        indicator->update_current(_advancement.completed());
        indicator->update_final(_advancement.total());
    }

    void _print_hold() {
        std::ostringstream logger_stream;
        logger_stream << "[" << _progress_indicator->symbol() << "] " << _progress_indicator->percentage() << "% ";
        logger_stream << " (w="<<std::setw(2)<<std::left<<_advancement.waiting()
                      << " p="<<std::setw(2)<<std::left<<_advancement.processing()
                      << " c="<<std::setw(3)<<std::left<<_advancement.completed()
                      << ")";
        Logger::instance().hold(_log_scope_manager->scope(),logger_stream.str());
    }

  protected:

    Void _enqueue(E const& e) {
        if (_using_concurrency()) {
            _advancement.add_to_waiting();
            auto task = std::bind(std::forward<TaskFunctionType const>(_task_func), std::forward<E const&>(e));
            auto progress_acknowledge = std::bind(std::forward<ProgressAcknowledgeFunctionType const>(_progress_acknowledge_func),
                                                  std::forward<E const&>(e), _progress_indicator);
            TaskManager::instance().enqueue([this,task,progress_acknowledge]{ _concurrent_task_wrapper(task,progress_acknowledge); });
        } else {
            // Locking and notification are used when concurrency is set to zero during processing, in order to avoid a race and to resume processing _sequential_queue respectively
            {
                LockGuard<Mutex> lock(_element_appending_mutex);
                append(e);
            }
            _element_availability_condition.notify_one();
        }
    }

  protected:

    TaskFunctionType _task_func;
    ProgressAcknowledgeFunctionType _progress_acknowledge_func;
    WorkloadAdvancement _advancement;

  private:

    // Queue of task-progress_acknowledge pairs for initial consumption and for consumption when using no concurrency
    std::queue<Pair<CompletelyBoundFunctionType,CompletelyBoundFunctionType>> _sequential_queue;

    Nat _logger_level; // The logger level to impose to the running threads
    SharedPointer<LogScopeManager> _log_scope_manager; // The scope manager required to properly hold print
    SharedPointer<ProgressIndicator> _progress_indicator; // The progress indicator to hold print

    Mutex _element_appending_mutex;
    Mutex _element_availability_mutex;
    ConditionVariable _element_availability_condition;

    ExceptionPtr _exception;
};

//! \brief A basic static workload where all elements are appended and then processed
template<class E, class... AS>
class StaticWorkload : public WorkloadBase<E,AS...> {
public:
    using TaskFunctionType = std::function<Void(E const&, AS...)>;

    StaticWorkload(TaskFunctionType f, AS... as) : WorkloadBase<E, AS...>() {
        this->_task_func = std::bind(std::forward<TaskFunctionType const>(f), std::placeholders::_1, std::forward<AS>(as)...);
    }

    //StaticWorkload& append(E const& e) { WorkloadBase<E,AS...>::append(e); return *this; }
    //StaticWorkload& append(List<E> const& es) { WorkloadBase<E,AS...>::append(es); return *this; }
};

//! \brief A dynamic workload in which it is possible to append new elements from the called function
template<class E, class... AS>
class DynamicWorkload : public WorkloadBase<E,AS...> {
  public:
    //! \brief Reduced interface to be used by the processing function (and any function called by it)
    class Access {
        friend DynamicWorkload;
    protected:
        Access(DynamicWorkload& parent) : _load(parent) { }
    public:
        Void append(E const &e) { _load._enqueue(e); }
    private:
        DynamicWorkload& _load;
    };
  public:
    using TaskFunctionType = std::function<Void(DynamicWorkload<E,AS...>::Access&, E const&, AS...)>;
    using ProgressAcknowledgeFunctionType = std::function<Void(E const&, SharedPointer<ProgressIndicator>)>;

    DynamicWorkload(ProgressAcknowledgeFunctionType p, TaskFunctionType t, AS... as) : WorkloadBase<E, AS...>(), _access(Access(*this)) {
        this->_task_func = std::bind(std::forward<TaskFunctionType const>(t),
                                     std::forward<DynamicWorkload<E,AS...>::Access const&>(_access),
                                     std::placeholders::_1,
                                     std::forward<AS>(as)...);
        this->_progress_acknowledge_func = p;
    }

  private:
    Access const _access;
};

}

#endif // ARIADNE_WORKLOAD_HPP