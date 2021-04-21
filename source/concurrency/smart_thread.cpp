/***************************************************************************
 *            concurrency/smart_thread.cpp
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

#include "smart_thread.hpp"

namespace Ariadne {

SmartThread::SmartThread(String name, std::function<Void(Void)> task, Bool start_immediately)
    : SmartThread(name, task, [](){}, [](){}, start_immediately)
{
}

SmartThread::SmartThread(String name, std::function<Void(Void)> task, std::function<Void(Void)> entry, std::function<Void(Void)> exit, Bool start_immediately)
    : _name(name), _entry(std::move(entry)), _exit(std::move(exit))
{
    _thread = std::thread([=,this]() {
        _id = std::this_thread::get_id();
        _got_id_promise.set_value();
        _has_started_future.get();
        task();
    });
    _got_id_future.get();
    if (start_immediately)
        start();
}

ThreadId SmartThread::id() const {
    return _id;
}

String SmartThread::name() const {
    return _name;
}

Void SmartThread::start()  {
    if (!_has_started) {
        _has_started = true;
        _entry();
        _has_started_promise.set_value();
    }
}

Bool SmartThread::has_started() {
    return _has_started;
}

SmartThread::~SmartThread() {
    if(!_has_started) _has_started_promise.set_value();
    _thread.join();
    _exit();
}

} // namespace Ariadne
