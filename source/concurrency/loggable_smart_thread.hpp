/***************************************************************************
 *            concurrency/loggable_smart_thread.hpp
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

/*! \file concurrency/loggable_smart_thread.hpp
 *  \brief A wrapper for smart handling of a thread
 */

#ifndef ARIADNE_LOGGABLE_SMART_THREAD_HPP
#define ARIADNE_LOGGABLE_SMART_THREAD_HPP

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

class LoggableSmartThread {
  public:
    LoggableSmartThread(String name, std::function<Void(Void)> task);

    ThreadId id() const;
    String name() const;

    Void activate();

    ~LoggableSmartThread();

  private:
    ThreadId _id;
    String _name;
    std::atomic<bool> _active = false;
    std::thread _thread;
    std::promise<void> _activate_promise;
    std::future<void> _activate_future = _activate_promise.get_future();
    std::promise<void> _start_promise;
    std::future<void> _start_future = _start_promise.get_future();
};

} // namespace Ariadne

#endif // ARIADNE_LOGGABLE_SMART_THREAD_HPP
