/***************************************************************************
 *            concurrency/smart_thread_pool.cpp
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

#include "smart_thread_pool.hpp"

namespace Ariadne {

const String THREAD_NAME_PREFIX = "thr";

VoidFunction SmartThreadPool::_task_wrapper_function() {
    return [=, this] {
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
            if (_finish_current_and_stop) {
                _num_stopped_threads++;
                return;
            }
        }
    };
}

SmartThreadPool::SmartThreadPool(SizeType size)
        : _finish_all_and_stop(false), _finish_current_and_stop(false), _num_stopped_threads(0) {
    for (SizeType i = 0; i < size; ++i) {
        _threads.append(make_shared<SmartThread>(SmartThreadPool::_task_wrapper_function(), THREAD_NAME_PREFIX + to_string(i)));
    }
}

SizeType SmartThreadPool::num_threads() const {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    return _threads.size();
}

Void SmartThreadPool::add_threads(SizeType number) {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    auto old_size = _threads.size();
    if (old_size == 0) {
        _finish_current_and_stop = false;
    }
    _threads.resize(old_size+number);
    for (SizeType i=old_size; i<old_size+number; ++i) {
        _threads.at(i) = make_shared<SmartThread>(SmartThreadPool::_task_wrapper_function(), THREAD_NAME_PREFIX + to_string(i));
    }
}

Void SmartThreadPool::schedule_stop_threads() {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    _finish_current_and_stop = true;
}

SizeType SmartThreadPool::num_stopped_threads() const {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    return _num_stopped_threads;
}

Void SmartThreadPool::remove_threads() {
    std::lock_guard<std::mutex> lock(_num_threads_mutex);
    ARIADNE_PRECONDITION(_num_stopped_threads == _threads.size());
    _threads.clear();
    _num_stopped_threads = 0;
}

SizeType SmartThreadPool::queue_size() const {
    std::lock_guard<std::mutex> lock(_task_availability_mutex);
    return _tasks.size();
}

SmartThreadPool::~SmartThreadPool() {
    {
        std::lock_guard<std::mutex> task_availability_lock(_task_availability_mutex);
        _finish_all_and_stop = true;
    }
    _task_availability_condition.notify_all();
    _threads.clear();
}

}
