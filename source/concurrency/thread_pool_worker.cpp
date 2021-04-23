/***************************************************************************
 *            concurrency/thread_pool_worker.cpp
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

#include "output/logging.hpp"
#include "thread_pool_worker.hpp"

namespace Ariadne {

ThreadPoolWorker::ThreadPoolWorker(VoidFunction task, String name)
        : _name(name), _got_id_future(_got_id_promise.get_future())
{
    _thread = std::thread([=,this]() {
        _id = std::this_thread::get_id();
        _got_id_promise.set_value();
        task();
    });
    _got_id_future.get();
    Logger::instance().register_thread(_id,_name);
}

ThreadId ThreadPoolWorker::id() const {
    return _id;
}

String ThreadPoolWorker::name() const {
    return _name;
}

ThreadPoolWorker::~ThreadPoolWorker() {
    _thread.join();
    Logger::instance().unregister_thread(_id);
}

} // namespace Ariadne
