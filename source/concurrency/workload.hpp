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
  public:
    using PartiallyBoundFunctionType = std::function<Void(E const &)>;
    using CompletelyBoundFunctionType = VoidFunction;

    WorkloadBase() : _advancement(0), _logger_level(0) { }

    Void process() override {
        _log_scope_manager.reset(new LogScopeManager(ARIADNE_PRETTY_FUNCTION,0));
        _progress_indicator.reset(new ProgressIndicator(_tasks.size()));
        _logger_level = Logger::instance().current_level();
        while (true) {
            UniqueLock<Mutex> lock(_element_availability_mutex);
            _element_availability_condition.wait(lock, [=,this] { return _advancement.has_finished() or not _tasks.empty() or _exception != nullptr; });
            if (_exception != nullptr) rethrow_exception(_exception);
            if (_advancement.has_finished()) { _log_scope_manager.reset(); return; }

            auto task = _tasks.back();
            _tasks.pop_back();
            if (_using_concurrency()) {
                TaskManager::instance().enqueue([this, task] { _concurrent_task_wrapper(task); });
            } else {
                _advancement.add_to_processing();
                _print_hold();
                task();
                _advancement.add_to_completed();
            }
        }
    }

    SizeType size() const override { return _tasks.size(); }

    Void append(E const& e) {
        _advancement.add_to_waiting();
        _tasks.push_back(std::bind(std::forward<PartiallyBoundFunctionType const>(_f), std::forward<E const&>(e)));
    }

    Void append(List<E> const& es) override { for (auto e : es) append(e); }

  private:

    Bool _using_concurrency() const { return TaskManager::instance().concurrency() > 0; }

    Void _concurrent_task_wrapper(CompletelyBoundFunctionType const& task) {
        _advancement.add_to_processing();
        _print_hold();
        try {
            Logger::instance().set_level(_logger_level);
            task();
        } catch (...) {
            {
                LockGuard<Mutex> exc_lock(_exception_mutex);
                _exception = std::current_exception();
            }
            _element_availability_condition.notify_one();
        }

        _advancement.add_to_completed();
        if (_advancement.has_finished()) _element_availability_condition.notify_one();
    }

    void _print_hold() {
        _progress_indicator->update_current(_advancement.completed());
        _progress_indicator->update_final(_advancement.total());
        if (not Logger::instance().is_muted_at(0)) {
            std::ostringstream logger_stream;
            logger_stream << "[" << _progress_indicator->symbol() << "] " << _progress_indicator->percentage() << "% ";
            logger_stream << " (w="<<std::setw(2)<<std::left<<_advancement.waiting()
                          << " p="<<std::setw(2)<<std::left<<_advancement.processing()
                          << " c="<<std::setw(3)<<std::left<<_advancement.completed()
                          << ")";
            Logger::instance().hold(_log_scope_manager->scope(),logger_stream.str());
        }
    }

  protected:

    Void _enqueue(E const& e) {
        if (_using_concurrency()) {
            _advancement.add_to_waiting();
            auto task = std::bind(std::forward<PartiallyBoundFunctionType const>(_f), std::forward<E const&>(e));
            TaskManager::instance().enqueue([this,task]{ _concurrent_task_wrapper(task); });
        } else {
            // Locking and notification are used when concurrency is set to zero during processing, in order to avoid a race and to resume processing _tasks respectively
            {
                LockGuard<Mutex> lock(_element_appending_mutex);
                append(e);
            }
            _element_availability_condition.notify_one();
        }
    }

  protected:

    WorkloadAdvancement _advancement;
    PartiallyBoundFunctionType _f;

  private:

    List<CompletelyBoundFunctionType> _tasks; // Used for initial consumption and for consumption when using no concurrency

    Nat _logger_level; // The logger level to impose to the running threads
    SharedPointer<LogScopeManager> _log_scope_manager; // The scope manager required to properly hold print
    SharedPointer<ProgressIndicator> _progress_indicator; // The progress indicator to hold print

    Mutex _element_appending_mutex;
    Mutex _element_availability_mutex;
    Mutex _exception_mutex;
    ConditionVariable _element_availability_condition;

    ExceptionPtr _exception;
};

//! \brief A basic static workload where all elements are appended and then processed
template<class E, class... AS>
class StaticWorkload : public WorkloadBase<E,AS...> {
public:
    using FunctionType = std::function<Void(E const&, AS...)>;

    StaticWorkload(FunctionType f, AS... as) : WorkloadBase<E, AS...>() {
        this->_f = std::bind(std::forward<FunctionType const>(f), std::placeholders::_1, std::forward<AS>(as)...);
    }
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
    using FunctionType = std::function<Void(DynamicWorkload<E,AS...>::Access&, E const&, AS...)>;

    DynamicWorkload(FunctionType f, AS... as) : WorkloadBase<E, AS...>(), _access(Access(*this)) {
        this->_f = std::bind(std::forward<FunctionType const>(f),
                             std::forward<DynamicWorkload<E,AS...>::Access const&>(_access),
                             std::placeholders::_1,
                             std::forward<AS>(as)...);
    }

  private:
    Access const _access;
};

}

#endif // ARIADNE_WORKLOAD_HPP