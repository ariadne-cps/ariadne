/***************************************************************************
 *            concurrency/concurrency_typedefs.hpp
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

/*! \file concurrency/concurrency_typedefs.hpp
 *  \brief Typedefs for the module.
 */

#ifndef ARIADNE_CONCURRENCY_TYPEDEFS_HPP
#define ARIADNE_CONCURRENCY_TYPEDEFS_HPP

#include <utility>
#include <mutex>
#include <thread>
#include <future>
#include "utility/typedefs.hpp"

namespace Ariadne {

using ConditionVariable = std::condition_variable;
using Mutex = std::mutex;
template<class T> using LockGuard = std::lock_guard<T>;
template<class T> using UniqueLock = std::unique_lock<T>;
using ThreadId = std::thread::id;
using VoidFunction = std::function<Void()>;
template<class T> using Future = std::future<T>;
template<class T> using Promise = std::promise<T>;
template<class T> using PackagedTask = std::packaged_task<T>;


} // namespace Ariadne

#endif // ARIADNE_CONCURRENCY_TYPEDEFS_HPP
