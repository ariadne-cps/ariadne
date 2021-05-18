/***************************************************************************
 *            io/command_line_interface.cpp
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

#include "command_line_interface.hpp"

namespace Ariadne {

inline unsigned int get_verbosity(int argc, const char* argv[]) {
    if(argc>1) {
        if(std::strcmp(argv[1],"-v")==0) {
            if(argc>2) {
                int val = std::atoi(argv[2]);
                ARIADNE_ASSERT_MSG(val >= 0,"Verbosity should be a non-negative value.\n");
                return static_cast<unsigned int>(val);
            }
        } else {
            ARIADNE_FAIL_MSG("Unrecognised command-line option \"" << argv[1] << "\"\n");
        }
    }
    return 0;
}

ArgumentStream::ArgumentStream(int argc, const char* argv[]) {
    ARIADNE_PRECONDITION(argc > 0)
    for (int i=1; i<argc; ++i) {
        _args.push(String(argv[i]));
    }
}

String ArgumentStream::peek() const {
    ARIADNE_PRECONDITION(not empty())
    return _args.front();
}

String ArgumentStream::pop() {
    auto val = peek();
    _args.pop();
    return val;
}

Bool ArgumentStream::empty() const {
    return _args.empty();
}

SizeType ArgumentStream::size() const {
    return _args.size();
}

//! \brief Supplies base behavior for argument parsers
class ArgumentParserBase : public ArgumentParserInterface {
  public:

};

CommandLineInterface::CommandLineInterface() : _parsers() { }

Void CommandLineInterface::acquire(int argc, const char* argv[]) const {
    ArgumentStream stream(argc,argv);
}

} // namespace Ariadne