/***************************************************************************
 *            concurrency/thread.hpp
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

/*! \file concurrency/thread.hpp
 *  \brief A wrapper for smart handling of a thread
 */

#ifndef ARIADNE_THREAD_HPP
#define ARIADNE_THREAD_HPP

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

using ExceptionPtr = std::exception_ptr;

//! \brief A class for handling a thread for a pool in a smarter way.
//! \details It allows to wait for the start of the \a task before extracting the thread id, which is held along with
//! a readable \a name.
class Thread {
  public:

    //! \brief Construct with an optional name.
    //! \details The thread will start and store the id
    Thread(VoidFunction task, String name = String());

    //! \brief Get the thread id
    ThreadId id() const;
    //! \brief Get the readable name
    String name() const;

    //! \brief The exception, if it exists
    ExceptionPtr const& exception() const;

    //! \brief Destroy the instance
    ~Thread();

  private:
    String _name;
    ThreadId _id;
    std::thread _thread;
    Promise<Void> _got_id_promise;
    Future<Void> _got_id_future;
    Promise<Void> _registered_thread_promise;
    Future<Void> _registered_thread_future;
    ExceptionPtr _exception;
};

} // namespace Ariadne

#endif // ARIADNE_THREAD_HPP
