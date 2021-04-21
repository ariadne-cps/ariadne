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

namespace Ariadne {

class String;

typedef std::thread::id ThreadId;
typedef std::function<Void(Void)> VoidFunction;

//! \brief A class for handling a thread in a smarter way.
//! \details It allows to wait for the start of the \a task before extracting the thread id, which is held along with
//! a readable \a name for logging purposes. An \a entry function can be specified to be run just before the task starts running, and
//! an \a exit function to be run just after the thread joins. The thread joins automatically when exiting scope.
//! The thread task can be started immediately, as soon as the id is retrieved, or not (true by default).
class SmartThread {
  public:
    SmartThread(String name, VoidFunction task, Bool start_immediately=true);
    SmartThread(String name, VoidFunction task, VoidFunction entry, VoidFunction exit, Bool start_immediately=true);

    //! \brief Get the thread id
    ThreadId id() const;
    //! \brief Get the readable name
    String name() const;

    //! \brief Start the thread. Can be done once for the object; if already started, this does not do anything.
    Void start();

    //! \brief Whether the thread is running the task.
    Bool has_started();

    ~SmartThread();

  private:
    String _name;
    VoidFunction _entry;
    VoidFunction _exit;
    ThreadId _id;
    std::atomic<bool> _has_started = false;
    std::thread _thread;
    std::promise<void> _has_started_promise;
    std::future<void> _has_started_future = _has_started_promise.get_future();
    std::promise<void> _got_id_promise;
    std::future<void> _got_id_future = _got_id_promise.get_future();
};

} // namespace Ariadne

#endif // ARIADNE_SMART_THREAD_HPP
