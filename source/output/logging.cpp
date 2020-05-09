/***************************************************************************
 *            output/logging.cpp
 *
 *  Copyright  2007-20  Pieter Collins
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

#include <iostream>

#include "../config.hpp"
#include "../output/logging.hpp"
#include "../utility/writable.hpp"

namespace Ariadne {

std::ofstream log_file_stream;

void redirect_log(const char* filename)
{
    if(log_file_stream.is_open()) {
        log_file_stream.close();
    }
    log_file_stream.open(filename);
    std::clog.rdbuf( log_file_stream.rdbuf() );
}

int global_verbosity=0;

OutputStream& operator<<(OutputStream& os, const WritableInterface& w);

void write_error(OutputStream& os, const char* i, const char* c, const char* t) {
    os << "Error in dynamic_handle_cast: cannot convert object of static type " << c << " and dynamic type " << i << " to type " << t << "\n";
}

void write_error(OutputStream& os, const WritableInterface* w, const char* i, const char* c, const char* t) {
//    os << "Error in dynamic_handle_cast\n";
    os << "Error in dynamic_handle_cast:" << std::flush;
    os << " cannot convert "; assert(w); w->_write(os); os << std::flush;
    os << " of static type " << c << " and dynamic type " << i << " to type " << t << std::endl;
}


} // namespace Ariadne

inline bool startup_function() {
    std::cerr << std::boolalpha;
    std::cout << std::boolalpha;
    std::clog << std::boolalpha;
    return true;
}

static const bool startup = startup_function();

