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


ThreadId SmartThread::id() const {
    return _id;
}

String SmartThread::name() const {
    return _name;
}

Void SmartThread::start()  {
    if (!_has_started) {
        _has_started = true;
        _has_started_promise.set_value();
    }
}

Bool SmartThread::has_finished() const {
    return _has_finished;
}

Bool SmartThread::has_started() const {
    return _has_started;
}

SmartThread::~SmartThread() {
    if(!_has_started) _has_started_promise.set_value();
    _thread.join();
    Logger::instance().unregister_thread(this->id());
}

} // namespace Ariadne
