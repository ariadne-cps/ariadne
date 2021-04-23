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

SmartThreadPool::SmartThreadPool(SizeType size)
        : _stop(false) {
    ARIADNE_PRECONDITION(size > 0);
    for (SizeType i = 0; i < size; ++i) {
        _threads.append(make_shared<SmartThread>(THREAD_NAME_PREFIX + to_string(i)));
        _threads.at(i)->enqueue(
                [=, this] {
                    while (true) {
                        VoidFunction task;
                        {
                            std::unique_lock<std::mutex> lock(_mutex);
                            _availability_condition.wait(lock, [=, this] { return _stop or not _tasks.empty(); });
                            if (_stop and _tasks.empty()) return;
                            task = std::move(_tasks.front());
                            _tasks.pop();
                        }
                        task();
                    }
                }
        );
    }
}

SizeType SmartThreadPool::num_threads() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return _threads.size();
}

Void SmartThreadPool::set_num_threads(SizeType number) {
    ARIADNE_PRECONDITION(number > 0);
    const SizeType previous_number = _threads.size();
    if (number > previous_number) {
        _threads.resize(number);
        for (SizeType i=previous_number; i<number; ++i)
            _threads.at(i) = make_shared<SmartThread>(THREAD_NAME_PREFIX + to_string(i));
    } else if (number < previous_number) {

    }
}

SizeType SmartThreadPool::queue_size() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return _tasks.size();
}

SmartThreadPool::~SmartThreadPool() {
    {
        std::lock_guard<std::mutex> lock(_mutex);
        _stop = true;
    }
    _availability_condition.notify_all();
    _threads.clear();
}

}
