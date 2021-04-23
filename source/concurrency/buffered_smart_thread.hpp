/***************************************************************************
 *            concurrency/buffered_smart_thread.hpp
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

/*! \file concurrency/buffered_smart_thread.hpp
 *  \brief A wrapper for smart handling of a thread with a buffer for incoming tasks
 */

#ifndef ARIADNE_BUFFERED_SMART_THREAD_HPP
#define ARIADNE_BUFFERED_SMART_THREAD_HPP

#include <utility>
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#include <functional>
#include "utility/typedefs.hpp"
#include "utility/string.hpp"
#include "utility/metaprogramming.hpp"
#include "concurrency/buffer.hpp"
#include "concurrency/concurrency_typedefs.hpp"

namespace Ariadne {

//! \brief A class for handling a thread in a smarter way.
//! \details It allows to wait for the start of the \a task before extracting the thread id, which is held along with
//! a readable \a name for logging purposes. The thread can execute only one task at a time, but multiple tasks can
//! be enqueued based on the internal buffer capacity.
class BufferedSmartThread {
  public:

    //! \brief Construct with a name.
    //! \details The thread will start and store the id, then register to the Logger
    BufferedSmartThread(String name);
    //! \brief Construct using the thread id as the name.
    BufferedSmartThread();

    //! \brief Enqueue a task for execution, returning the future handler
    //! \details If the buffer is full, successive calls will block until an execution is started.
    template<class F, class... AS>
    auto enqueue(F&& f, AS&&... args) -> Future<ResultOf<F(AS...)>>;

    //! \brief Get the thread id
    ThreadId id() const;
    //! \brief Get the readable name
    String name() const;

    //! \brief The current size of the queue
    SizeType queue_size() const;
    //! \brief The capacity of the tasks to execute
    SizeType queue_capacity() const;
    //! \brief Change the queue capacity
    //! \details Capacity cannot be changed to a value lower than the current size
    Void set_queue_capacity(SizeType capacity);

    //! \brief Destroy the instance, also unregistering from the Logger
    ~BufferedSmartThread();

  private:
    String _name;
    ThreadId _id;
    Thread _thread;
    Buffer<VoidFunction> _task_buffer;
    Promise<void> _got_id_promise;
    Future<void> _got_id_future;
};

template<class F, class... AS> auto BufferedSmartThread::enqueue(F&& f, AS&&... args) -> Future<ResultOf<F(AS...)>>
{
    using ReturnType = ResultOf<F(AS...)>;

    auto task = std::make_shared<PackagedTask<ReturnType()>>(std::bind(std::forward<F>(f), std::forward<AS>(args)...));
    Future<ReturnType> result = task->get_future();
    _task_buffer.push([task](){ (*task)(); });
    return result;
}

} // namespace Ariadne

#endif // ARIADNE_BUFFERED_SMART_THREAD_HPP
