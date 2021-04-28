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

const String THREAD_NAME_PREFIX = "thr";

VoidFunction ThreadPool::_task_wrapper_function(SizeType i) {
    return [i, this] {
        while (true) {
            VoidFunction task;
            {
                std::unique_lock<std::mutex> lock(_task_availability_mutex);
                _task_availability_condition.wait(lock, [=, this] { return _finish_all_and_stop or not _tasks.empty(); });
                if (_finish_all_and_stop and _tasks.empty()) return;
                task = std::move(_tasks.front());
                _tasks.pop();
            }
            task();
            if (i>=_num_threads_to_use) {
                std::unique_lock<std::mutex> lock(_num_active_threads_mutex);
                _num_active_threads--;
                if (_num_active_threads == _num_threads_to_use) _all_unused_threads_stopped_promise.set_value();
                return;
            }
        }
    };
}

ThreadPool::ThreadPool(SizeType size)
        : _finish_all_and_stop(false), _num_active_threads(size), _num_threads_to_use(size),
          _all_unused_threads_stopped_future(_all_unused_threads_stopped_promise.get_future())
{
    for (SizeType i = 0; i < size; ++i) {
        _threads.append(make_shared<Thread>(ThreadPool::_task_wrapper_function(i), THREAD_NAME_PREFIX + to_string(i)));
    }
}

SizeType ThreadPool::num_threads() const {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    return _threads.size();
}

Void ThreadPool::set_num_threads(SizeType number) {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    auto old_size = _threads.size();
    _num_threads_to_use = number;
    if (number > old_size) {
        _threads.resize(number);
        _num_active_threads = number;
        for (SizeType i=old_size; i<number; ++i) {
            _threads.at(i) = make_shared<Thread>(ThreadPool::_task_wrapper_function(i), THREAD_NAME_PREFIX + to_string(i));
        }
    } else if (number < old_size) {
        _all_unused_threads_stopped_future.get();
        _threads.resize(number);
        _all_unused_threads_stopped_promise = Promise<Void>();
        _all_unused_threads_stopped_future = _all_unused_threads_stopped_promise.get_future();
    }
}

SizeType ThreadPool::queue_size() const {
    std::lock_guard<std::mutex> lock(_task_availability_mutex);
    return _tasks.size();
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> task_availability_lock(_task_availability_mutex);
        _finish_all_and_stop = true;
    }
    _task_availability_condition.notify_all();
    _threads.clear();
}

}
