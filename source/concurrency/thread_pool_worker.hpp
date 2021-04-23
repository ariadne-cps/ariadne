/***************************************************************************
 *            concurrency/thread_pool_worker.hpp
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

/*! \file concurrency/thread_pool_worker.hpp
 *  \brief A wrapper for smart handling of a thread
 */

#ifndef ARIADNE_THREAD_POOL_WORKER_HPP
#define ARIADNE_THREAD_POOL_WORKER_HPP

#include <utility>
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#include <functional>
#include "utility/typedefs.hpp"
#include "utility/string.hpp"
#include "utility/metaprogramming.hpp"
#include "concurrency/concurrency_typedefs.hpp"

namespace Ariadne {

//! \brief A class for handling a thread for a pool in a smarter way.
//! \details It allows to wait for the start of the \a task before extracting the thread id, which is held along with
//! a readable \a name for logging purposes.
class ThreadPoolWorker {
  public:

    //! \brief Construct with an optional name.
    //! \details The thread will start and store the id, then register to the Logger
    ThreadPoolWorker(VoidFunction task, String name = String());

    //! \brief Get the thread id
    ThreadId id() const;
    //! \brief Get the readable name
    String name() const;

    //! \brief Destroy the instance, also unregistering from the Logger
    ~ThreadPoolWorker();

  private:
    String _name;
    ThreadId _id;
    Thread _thread;
    Promise<void> _got_id_promise;
    Future<void> _got_id_future;
};

} // namespace Ariadne

#endif // ARIADNE_THREAD_POOL_WORKER_HPP
