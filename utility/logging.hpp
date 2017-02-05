/***************************************************************************
 *            logging.hpp
 *
 *  Copyright 2007-17  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file logging.hpp
 *  \brief Support for writing debugging output to a logging stream.
 */

#ifndef ARIADNE_LOGGING_HPP
#define ARIADNE_LOGGING_HPP

#include <iostream>
#include <fstream>

// (Placeholder constant required for compilation)
static const std::string charcode="";

//! Send a message to the global logging stream.
#define ARIADNE_LOG(level,msg) \
    if(verbosity >= level) { \
        std::clog << "[" << charcode << ":" << level << "] "; \
        for(uint _i=0; _i!=level; ++_i) { std::clog<<' '; } \
        std::clog << msg << std::flush; \
    }

namespace Ariadne {

struct Loggable {
  public:
    Loggable() : verbosity(0),charcode("") { }
    mutable int verbosity;
  protected:
    mutable std::string charcode;
};

// Global log output file
extern std::ofstream log_file_stream;

//! \brief Redirect logging output to file \a filename.
void redirect_log(const char* filename);


} // namespace Ariadne

#endif // ARIADNE_LOGGING_HPP
