/***************************************************************************
 *            concurrency/thread.cpp
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
#include "thread.hpp"

namespace Ariadne {

Thread::Thread(VoidFunction task, String name)
        : _name(name), _got_id_future(_got_id_promise.get_future())
{
    _thread = std::thread([=,this]() {
        _id = std::this_thread::get_id();
        _got_id_promise.set_value();
        task();
    });
    _got_id_future.get();
    if (_name == String()) _name = to_string(_id);
    Logger::instance().register_thread(this->id(),this->name());
}

ThreadId Thread::id() const {
    return _id;
}

String Thread::name() const {
    return _name;
}

Thread::~Thread() {
    _thread.join();
    Logger::instance().unregister_thread(this->id());
}

} // namespace Ariadne