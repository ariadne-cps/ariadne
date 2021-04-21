/***************************************************************************
 *            output/loggable_smart_thread.cpp
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
#include "logging.hpp"

namespace Ariadne {

LoggableSmartThread::LoggableSmartThread(String name, std::function<Void(Void)> task, Bool start_immediately)
    : SmartThread(name,task,[this](){ Logger::instance().register_thread(this->id(),this->name()); },
                                         [this](){ Logger::instance().unregister_thread(this->id()); },
                                         start_immediately)
{
}

} // namespace Ariadne
