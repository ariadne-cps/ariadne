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
#include "loggable_thread.hpp"

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

    //! \brief Add the \a number of threads
    Void add_threads(SizeType number);

    //! \brief Schedule the threads to stop
    //! \details The threads will stop as soon as their last task is completed;
    //! threads in the waiting state will not be able to stop, i.e., this method is
    //! not meant to be used on an inactive pool but an active one.
    Void schedule_stop_threads();

    //! \brief Force the threads to stop and remove the threads
    Void stop_and_remove_threads();

    //! \brief The number of threads already stopped but not removed
    SizeType num_stopped_threads() const;

    //! \brief Remove the stopped threads
    //! \details Throws an exception if there are still threads to stop
    Void remove_threads();

    ~ThreadPool();

  private:
    List<SharedPointer<LoggableThread>> _threads;
    std::queue<VoidFunction> _tasks;

    VoidFunction _task_wrapper_function();

    mutable std::mutex _task_availability_mutex;
    std::condition_variable _task_availability_condition;
    Bool _finish_all_and_stop; // Wait till the queue is empty before stopping the thread, used for destruction
    Bool _finish_current_and_stop; // Wait till the current one is completed before stopping the thread
    std::atomic<Nat> _num_stopped_threads;
    mutable std::mutex _num_threads_mutex;
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
