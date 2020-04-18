/***************************************************************************
 *            utility/stopwatch.hpp
 *
 *  Copyright  2018-20  Luca Geretti
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

/*! \file utility/stopwatch.hpp
 *  \brief Stopwatch class to be used for profiling execution times.
 */

#ifndef ARIADNE_STOPWATCH_HPP
#define ARIADNE_STOPWATCH_HPP

#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

namespace Ariadne {

using TimeUnit = double;

static const clock_t Hz = static_cast<clock_t>(sysconf(_SC_CLK_TCK));

class StopWatch {
public:
    StopWatch() { reset(); }

    TimeUnit elapsed() const { return TimeUnit(_clicked_time.tms_utime - _start_time.tms_utime)/TimeUnit(Hz); }
    void reset() { times(&_start_time); _clicked_time = _start_time; }
    void click() { times(&_clicked_time); }

private:
    tms _start_time;
    tms _clicked_time;
};

} // namespace Ariadne

#endif /* ARIADNE_STOPWATCH_HPP */
