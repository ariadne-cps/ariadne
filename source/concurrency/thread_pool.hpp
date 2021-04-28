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
#include "thread.hpp"

namespace Ariadne {

//! \brief Exception for stopping a thread pool
class StoppedThreadPoolException : public std::exception { };

//! \brief A pool of Thread objects managed internally given a (variable) number of threads
//! \details Differently from managing a single BufferedThread, the task queue for a pool is not upper-bounded, i.e., BufferedThread
//! objects use a buffer of one element, which receives once the wrapped task that consumes elements from the task queue.
class ThreadPool {
  public:
    //! \brief Construct from a given number of threads
    ThreadPool(SizeType num_threads);

    //! \brief Enqueue a task for execution, returning the future handler
    //! \details The is no limits on the number of tasks to enqueue
    template<class F, class... AS> auto enqueue(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>>;

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
    List<SharedPointer<Thread>> _threads;
    std::queue<VoidFunction> _tasks;

    //! \brief The function wrapper handling the extraction from the queue
    //! \details Takes \a i as the index of the thread in the list, for identification when stopping selectively
    VoidFunction _task_wrapper_function(SizeType i);

    mutable std::mutex _task_availability_mutex;
    std::condition_variable _task_availability_condition;
    Bool _finish_all_and_stop; // Wait till the queue is empty before stopping the thread, used for destruction
    Nat _num_active_threads; // Down-counter for checking whether all the threads to stop have been stopped
    Nat _num_threads_to_use; // Reference on the number of threads to use: if lower than the threads size, the last threads will stop
    mutable std::mutex _num_active_threads_mutex;
    mutable std::mutex _num_threads_mutex;
    Promise<Void> _all_unused_threads_stopped_promise;
    Future<Void> _all_unused_threads_stopped_future;
};

template<class F, class... AS>
auto ThreadPool::enqueue(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>> {
    using ReturnType = ResultOf<F(AS...)>;

    auto task = std::make_shared<PackagedTask<ReturnType()> >(std::bind(std::forward<F>(f), std::forward<AS>(args)...));
    Future<ReturnType> result = task->get_future();
    {
        std::unique_lock<std::mutex> lock(_task_availability_mutex);
        if (_finish_all_and_stop) throw StoppedThreadPoolException();
        _tasks.emplace([task]() { (*task)(); });
    }
    _task_availability_condition.notify_one();
    return result;
}

}

#endif // ARIADNE_THREAD_POOL_HPP
