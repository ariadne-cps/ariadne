/***************************************************************************
 *            concurrency/smart_thread_pool.hpp
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

#ifndef ARIADNE_SMART_THREAD_POOL_HPP
#define ARIADNE_SMART_THREAD_POOL_HPP

#include "utility/container.hpp"
#include "utility/pointer.hpp"
#include "smart_thread.hpp"

namespace Ariadne {

//! \brief Exception for stopping a thread pool
class StoppedThreadPoolException : public std::exception { };

//! \brief A pool of SmartThread objects managed internally given a capacity
class SmartThreadPool {
  public:
    SmartThreadPool(SizeType size);

    template<class F, class... AS>
    auto execute(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>>;

    //! \brief The number of threads
    SizeType num_threads();
    //! \brief The size of the tasks queue
    SizeType queue_size();

    ~SmartThreadPool();

  private:
    List<SharedPointer<SmartThread>> _threads;
    std::queue<VoidFunction> _tasks;

    std::mutex _mutex;
    std::condition_variable _availability_condition;
    Bool _stop;
};

template<class F, class... AS>
auto SmartThreadPool::execute(F &&f, AS &&... args) -> Future<ResultOf<F(AS...)>> {
    using ReturnType = ResultOf<F(AS...)>;

    auto task = std::make_shared<PackagedTask<ReturnType()> >(std::bind(std::forward<F>(f), std::forward<AS>(args)...));
    Future<ReturnType> result = task->get_future();
    {
        std::unique_lock<std::mutex> lock(_mutex);
        if (_stop) throw StoppedThreadPoolException();
        _tasks.emplace([task]() { (*task)(); });
    }
    _availability_condition.notify_one();
    return result;
}

}

#endif // ARIADNE_SMART_THREAD_POOL_HPP
