/***************************************************************************
 *            logging.h
 *
 *  Copyright  2007-8  Alberto Casagrande, Pieter Collins
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

/*! \file logging.h
 *  \brief Support for writing debugging output to a logging stream.
 */

#ifndef ARIADNE_LOGGING_H
#define ARIADNE_LOGGING_H

#include <iostream>
#include <fstream>

//! Send a message to the global logging stream.
#define ARIADNE_LOG(level,msg) \
    if(verbosity >= level) { \
        std::clog << "[" << level << "]"; \
        for(uint i=0; i!=level; ++i) { std::clog<<' '; } \
        std::clog << msg << std::flush; \
    }

namespace Ariadne {

struct Loggable {
  public:
    Loggable() : verbosity(0) { }
    mutable int verbosity;
};

// Global log output file
extern std::ofstream log_file_stream;

//! \brief Redirect logging output to file \a filename.
void redirect_log(const char* filename);

extern int global_verbosity;


} // namespace Ariadne

#endif // ARIADNE_LOGGING_H
