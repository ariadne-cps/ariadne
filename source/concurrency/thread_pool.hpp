/***************************************************************************
 *            concurrency/thread_pool.hpp
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

/*! \file concurrency/smart_thread_pool.hpp
 *  \brief A pool of smart threads
 */

#ifndef ARIADNE_THREAD_POOL_HPP
#define ARIADNE_THREAD_POOL_HPP

#include <queue>
#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "conclog/include/logging.hpp"
#include "thread.hpp"

using namespace ConcLog;

namespace Ariadne {

const String THREAD_POOL_DEFAULT_NAME = "thr";

//! \brief Exception for stopping a thread pool
class StoppedThreadPoolException : public std::exception { };

//! \brief A pool of Thread objects managed internally given a (variable) number of threads
//! \details Differently from managing a single BufferedThread, the task queue for a pool is not upper-bounded, i.e., BufferedThread
//! objects use a buffer of one element, which receives once the wrapped task that consumes elements from the task queue.
class ThreadPool {
  public:
    //! \brief Construct from a given number of threads and possibly a name
    ThreadPool(SizeType num_threads, String name = THREAD_POOL_DEFAULT_NAME);

    //! \brief Enqueue a task for execution, returning the future handler
    //! \details The is no limits on the number of tasks to enqueue
    template<class F, class... AS> auto enqueue(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>>;

    //! \brief The name of the pool
    String name() const;

    //! \brief The size of the tasks queue
    SizeType queue_size() const;

    //! \brief The number of threads
    SizeType num_threads() const;

    //! \brief Set the number of threads
    //! \details If reducing the current number, this method will block until
    //! all the previous tasks are completed, previous threads are destroyed
    //! and new threads are spawned
    Void set_num_threads(SizeType number);

    ~ThreadPool();

  private:

    //! \brief The function wrapper handling the extraction from the queue
    //! \details Takes \a i as the index of the thread in the list, for identification when stopping selectively
    VoidFunction _task_wrapper_function(SizeType i);
    //! \brief Append threads in the given range
    Void _append_thread_range(SizeType lower, SizeType upper);

  private:
    const String _name;
    List<SharedPointer<Thread>> _threads;
    std::queue<VoidFunction> _tasks;

    mutable Mutex _task_availability_mutex;
    ConditionVariable _task_availability_condition;
    Bool _finish_all_and_stop; // Wait till the queue is empty before stopping the thread, used for destruction
    Nat _num_active_threads; // Down-counter for checking whether all the threads to stop have been stopped
    Nat _num_threads_to_use; // Reference on the number of threads to use: if lower than the threads size, the last threads will stop
    mutable Mutex _num_active_threads_mutex;
    mutable Mutex _num_threads_mutex;
    Promise<Void> _all_unused_threads_stopped_promise;
    Future<Void> _all_unused_threads_stopped_future;
};

template<class F, class... AS>
auto ThreadPool::enqueue(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>> {
    using ReturnType = ResultOf<F(AS...)>;

    auto task = std::make_shared<PackagedTask<ReturnType()> >(std::bind(std::forward<F>(f), std::forward<AS>(args)...));
    Future<ReturnType> result = task->get_future();
    {
        UniqueLock<Mutex> lock(_task_availability_mutex);
        if (_finish_all_and_stop) throw StoppedThreadPoolException();
        _tasks.emplace([task]{ (*task)(); });
    }
    _task_availability_condition.notify_one();
    return result;
}

//! \brief Utility function to construct a thread name from a \a prefix and a \a number,
//! accounting for a maximum number of threads given by \a max_number
String construct_thread_name(String prefix, SizeType number, SizeType max_number);

}

#endif // ARIADNE_THREAD_POOL_HPP
