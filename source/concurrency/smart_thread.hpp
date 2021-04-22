/***************************************************************************
 *            concurrency/smart_thread.hpp
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

/*! \file concurrency/smart_thread.hpp
 *  \brief A wrapper for smart handling of a thread
 */

#ifndef ARIADNE_SMART_THREAD_HPP
#define ARIADNE_SMART_THREAD_HPP

#include <utility>
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#include "utility/typedefs.hpp"
#include "utility/string.hpp"
#include "output/logging.hpp"

namespace Ariadne {

class String;

typedef std::thread::id ThreadId;

//! \brief A class for handling a thread in a smarter way.
//! \details It allows to wait for the start of the \a task before extracting the thread id, which is held along with
//! a readable \a name for logging purposes. The thread does not start immediately.
class SmartThread {
  public:

    template<class Function, class... Args>
    explicit SmartThread(String name, Function&& f, Args&&... args);

    //! \brief Get the thread id
    ThreadId id() const;
    //! \brief Get the readable name
    String name() const;

    //! \brief Start the thread. Can be done once for the object; if already started, this does not do anything.
    Void start();

    //! \brief Whether the thread is running the task.
    Bool has_started() const;

    //! \brief Whether the thread has finished the task.
    Bool has_finished() const;

    ~SmartThread();

  private:
    String _name;
    ThreadId _id;
    std::atomic<bool> _has_started = false;
    std::atomic<bool> _has_finished = false;
    std::thread _thread;
    std::promise<void> _has_started_promise;
    std::future<void> _has_started_future = _has_started_promise.get_future();
    std::promise<void> _got_id_promise;
    std::future<void> _got_id_future = _got_id_promise.get_future();
};

template<class Function, class... Args> SmartThread::SmartThread(String name, Function&& f, Args&&... args)
    : _name(name)
{
    _thread = std::thread([=,this]() {
        _id = std::this_thread::get_id();
        _got_id_promise.set_value();
        _has_started_future.get();
        f(args...);
        _has_finished = true;
    });
    _got_id_future.get();
    Logger::instance().register_thread(_id,_name);
}

} // namespace Ariadne

#endif // ARIADNE_SMART_THREAD_HPP
