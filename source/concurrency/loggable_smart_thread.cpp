/***************************************************************************
 *            concurrency/loggable_smart_thread.cpp
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

#include "loggable_smart_thread.hpp"
#include "output/logging.hpp"

namespace Ariadne {

LoggableSmartThread::LoggableSmartThread(String name, std::function<Void(Void)> task) {
    _name = name;
    _thread = std::thread([=,this]() {
                    _id = std::this_thread::get_id();
                    _start_promise.set_value();
                    _activate_future.get();
                    if(_active) task();
             });
    _start_future.get();
}

ThreadId LoggableSmartThread::id() const {
    return _id;
}

String LoggableSmartThread::name() const {
    return _name;
}

Void LoggableSmartThread::activate()  {
    if (!_active) {
        _active = true;
        Logger::instance().register_thread(_id,_name);
        _activate_promise.set_value();
    }
}

LoggableSmartThread::~LoggableSmartThread() {
    if(!_active) _activate_promise.set_value();
    _thread.join();
    Logger::instance().unregister_thread(_id);
}

} // namespace Ariadne
