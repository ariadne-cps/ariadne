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

#include "output/logging.hpp"
#include "smart_thread.hpp"

namespace Ariadne {

SmartThread::SmartThread(String name)
        : _name(name), _function_buffer(1), _got_id_future(_got_id_promise.get_future())
{
    _thread = std::thread([=,this]() {
        _id = std::this_thread::get_id();
        _got_id_promise.set_value();
        while(true) {
            try {
                VoidFunction task = _function_buffer.pull();
                task();
            } catch(BufferStoppedConsumingException& e) { break; }
        }
    });
    _got_id_future.get();
    Logger::instance().register_thread(_id,_name);
}

SmartThread::SmartThread() : SmartThread(String()) {
    _name = to_string(_id);
}

ThreadId SmartThread::id() const {
    return _id;
}

String SmartThread::name() const {
    return _name;
}

SmartThread::~SmartThread() {
    _function_buffer.stop_consuming();
    _thread.join();
    Logger::instance().unregister_thread(_id);
}

} // namespace Ariadne
