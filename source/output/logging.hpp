/***************************************************************************
 *            output/logging.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/*! \file output/logging.hpp
 *  \brief Support for writing debugging output to a logging stream.
 */

#ifndef ARIADNE_LOGGING_HPP
#define ARIADNE_LOGGING_HPP

#include <iostream>
#include <fstream>
#include <cstring>

// (Placeholder constants required for compilation)
static const std::string charcode="";
static const unsigned int verbosity=0;

//! Send a message to the global logging stream.
#define ARIADNE_LOG(level,msg) \
    if(verbosity >= level) { \
        std::clog << "[" << charcode << ":" << level << "]"; \
        for(uint _i=0; _i!=level; ++_i) { std::clog<<' '; } \
        std::clog << msg << std::flush; \
    }

namespace Ariadne {

inline unsigned int get_verbosity(int argc, const char* argv[]) {
    if(argc>1) {
        if(std::strcmp(argv[1],"-v")==0) {
            if(argc>2) {
                int val = std::atoi(argv[2]);
                if (val < 0) std::cerr << "Verbosity should be a non-negative value.\n";
                return static_cast<unsigned int>(val);
            }
        } else {
            std::cerr << "Unrecognised command-line option \"" << argv[1] << "\"\n";
        }
    }
    return 0;
}

struct Loggable {
  public:
    Loggable() : verbosity(0),charcode("") { }
    mutable unsigned int verbosity;
    void set_verbosity(unsigned int v) const { verbosity=v; }
  protected:
    mutable std::string charcode;
};

// Global log output file
extern std::ofstream log_file_stream;

//! \brief Redirect logging output to file \a filename.
void redirect_log(const char* filename);

} // namespace Ariadne

#endif // ARIADNE_LOGGING_HPP
