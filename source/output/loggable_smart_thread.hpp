/***************************************************************************
 *            output/loggable_smart_thread.hpp
 *
 *  Copyright  2007-20  Luca Geretti
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

/*! \file output/loggable_smart_thread.hpp
 *  \brief A smart thread with entry/exit handlers for registration/unregistration
 */

#ifndef ARIADNE_LOGGABLE_SMART_THREAD_HPP
#define ARIADNE_LOGGABLE_SMART_THREAD_HPP

#include "concurrency/smart_thread.hpp"

namespace Ariadne {

class LoggableSmartThread : public SmartThread {
  public:
    LoggableSmartThread(String name, std::function<Void(Void)> task, Bool start_immediately=false);

    virtual ~LoggableSmartThread() = default;
};

} // namespace Ariadne

#endif // ARIADNE_LOGGABLE_SMART_THREAD_HPP
