/***************************************************************************
 *            concurrency/thread_pool.cpp
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

#include "thread_pool.hpp"

namespace Ariadne {

String construct_thread_name(String prefix, SizeType number, SizeType max_number) {
    std::ostringstream ss;
    ss << prefix;
    if (max_number > 9 and number <= 9) ss << "0";
    ss << number;
    return ss.str();
}

VoidFunction ThreadPool::_task_wrapper_function(SizeType i) {
    return [i, this] {
        while (true) {
            VoidFunction task;
            {
                UniqueLock<Mutex> lock(_task_availability_mutex);
                _task_availability_condition.wait(lock, [=, this] { return _finish_all_and_stop or not _tasks.empty(); });
                if (_finish_all_and_stop and _tasks.empty()) return;
                task = std::move(_tasks.front());
                _tasks.pop();
            }
            task();
            if (i>=_num_threads_to_use) {
                UniqueLock<Mutex> lock(_num_active_threads_mutex);
                _num_active_threads--;
                if (_num_active_threads == _num_threads_to_use) _all_unused_threads_stopped_promise.set_value();
                return;
            }
        }
    };
}

Void ThreadPool::_append_thread_range(SizeType lower, SizeType upper) {
    for (SizeType i=lower; i<upper; ++i) {
        _threads.append(make_shared<Thread>(ThreadPool::_task_wrapper_function(i), construct_thread_name(_name,i,upper)));
    }
}

ThreadPool::ThreadPool(SizeType size, String name)
        : _name(name), _finish_all_and_stop(false), _num_active_threads(size), _num_threads_to_use(size),
          _all_unused_threads_stopped_future(_all_unused_threads_stopped_promise.get_future())
{
    _append_thread_range(0,size);
}

String ThreadPool::name() const {
    return _name;
}

SizeType ThreadPool::num_threads() const {
    LockGuard<Mutex> lock(_num_threads_mutex);
    return _threads.size();
}

Void ThreadPool::set_num_threads(SizeType number) {
    LockGuard<Mutex> lock(_num_threads_mutex);
    auto old_size = _threads.size();
    _num_threads_to_use = number;
    if (number > old_size) {
        _num_active_threads = number;
        _append_thread_range(old_size,number);
    } else if (number < old_size) {
        _all_unused_threads_stopped_future.get();
        _threads.resize(number);
        _all_unused_threads_stopped_promise = Promise<Void>();
        _all_unused_threads_stopped_future = _all_unused_threads_stopped_promise.get_future();
    }
}

SizeType ThreadPool::queue_size() const {
    LockGuard<Mutex> lock(_task_availability_mutex);
    return _tasks.size();
}

ThreadPool::~ThreadPool() {
    {
        LockGuard<Mutex> task_availability_lock(_task_availability_mutex);
        _finish_all_and_stop = true;
    }
    _task_availability_condition.notify_all();
    _threads.clear();
}

}
