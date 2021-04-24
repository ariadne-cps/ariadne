/***************************************************************************
 *            concurrency/loggable_thread.hpp
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

/*! \file concurrency/loggable_thread.hpp
 *  \brief A wrapper for smart handling of a thread, also handling registration/deregistration to the Logger
 */

#ifndef ARIADNE_LOGGABLE_THREAD_HPP
#define ARIADNE_LOGGABLE_THREAD_HPP

#include "output/logging.hpp"
#include "thread.hpp"

namespace Ariadne {

//! \brief A thread extension with handling of registration/deregistration to the logger
class LoggableThread : public Thread {
public:
    LoggableThread(VoidFunction task, String name = String()) : Thread(task, name) {
        Logger::instance().register_thread(this->id(),this->name());
    }

    ~LoggableThread() {
        Logger::instance().unregister_thread(this->id());
    }
};

} // namespace Ariadne

#endif // ARIADNE_LOGGABLE_THREAD_HPP
