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
#include "concurrency_typedefs.hpp"
#include "task_manager.hpp"
#include "workload_advancement.hpp"

namespace Ariadne {

//! \brief A workload expressed as a stack of elements to work on, supplied with a function to process them
//! \details E: stack element type
//!          AS: optional input arguments for processing the elements; if used as output, their synchronisation
//!              in the concurrent case is up to the designer
//!          The workload handles the non-concurrent case separately, in order to unroll the tasks breadth-first: if
//!          tasks were instead enqueued to the TaskManager, a depth-first execution would be performed.
//!          A workload object has a lifetime of multiple processings, if necessary.
template<class E, class... AS>
class Workload {
public:
    //! \brief Interface for appending elements while processing the workload in the function
    class Appender {
        friend Workload;
    protected:
        Appender(Workload& parent) : _load(parent) { }
    public:
        Void append(E const &e) { _load._enqueue(e); }
    private:
        Workload& _load;
    };
public:
    using FunctionType = std::function<Void(Workload<E,AS...>::Appender&, E, AS...)>;
    using BoundFunctionType = std::function<Void(E)>;

    Workload(FunctionType f, AS... as) :
            _appender(Appender(*this)),
            _f(std::bind(
                std::forward<FunctionType const>(f),
                std::forward<Workload<E,AS...>::Appender const&>(_appender),
                std::placeholders::_1,
                std::forward<AS>(as)...)),
            _advancement(0) { }

    //! \brief Process the given elements until completion
    Void process() {
        ARIADNE_PRECONDITION(not _elements.empty());
        while (true) {
            UniqueLock<Mutex> lock(_element_availability_mutex);
            _element_availability_condition.wait(lock, [=,this] { return _advancement.has_finished() or not _elements.empty(); });
            if (_advancement.has_finished()) return;

            auto e = _elements.back();
            _elements.pop_back();
            if (_using_concurrency()) {
                TaskManager::instance().enqueue([this, e] {
                    _advancement.add_to_processing();
                    (*e)();
                    _advancement.add_to_completed();
                    if (_advancement.has_finished()) _element_availability_condition.notify_one();
                });
            } else {
                _advancement.add_to_processing();
                (*e)();
                _advancement.add_to_completed();
            }
        }
    }

    //! \brief The size of the workload, i.e., the number of elements to process
    SizeType size() const { return _elements.size(); }

    //! \brief Append one element to process
    Void append(E const& e) {
        _advancement.add_to_waiting();
        _elements.push_back(std::make_shared<PackagedTask<Void()>>(std::bind(
                std::forward<BoundFunctionType const>(_f),
                std::forward<E const&>(e))));
    }

    //! \brief Append a list of elements to process
    Void append(List<E> const& es) { for (auto e : es) append(e); }

  private:

    Bool _using_concurrency() const { return TaskManager::instance().concurrency() > 0; }

    Void _enqueue(E const& e) {
        if (_using_concurrency()) {
            _advancement.add_to_waiting();
            TaskManager::instance().enqueue([=,this](E const& e){
                    _advancement.add_to_processing();
                    _f(e);
                    _advancement.add_to_completed();
                    if (_advancement.has_finished()) _element_availability_condition.notify_one();
                }, e);
        } else {
            // Locking and notification are used when concurrency is set to zero during processing, in order to avoid a race and to resume processing _elements
            {
                LockGuard<Mutex> lock(_element_appending_mutex);
                append(e);
            }
            _element_availability_condition.notify_one();
        }
    }

  private:
    List<SharedPointer<PackagedTask<Void()>>> _elements; // Used for initial consumption and for consumption when using no concurrency
    Appender const _appender;
    BoundFunctionType const _f;

    WorkloadAdvancement _advancement;

    Mutex _element_appending_mutex;
    Mutex _element_availability_mutex;
    ConditionVariable _element_availability_condition;
};

}

#endif // ARIADNE_WORKLOAD_HPP