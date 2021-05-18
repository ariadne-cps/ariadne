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

#include "concurrency/task_manager.hpp"
#include "logging.hpp"
#include "command_line_interface.hpp"

namespace Ariadne {

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
  protected:
    ArgumentParserBase(String const& short_id, String const& long_id, String const& instructions)
        : _short_id(short_id), _long_id(long_id), _instructions(instructions) { }
  public:
    String short_id() const { return _short_id; }
    String long_id() const { return _long_id; }
    Bool has_short_id() const { return not _short_id.empty(); }
    String instructions() const { return _instructions; }

    Bool is_consumable(ArgumentStream const& stream) const override;
    String help_description() const override;
  private:
    String _short_id;
    String _long_id;
    String _instructions;
};

String ArgumentParserBase::help_description() const {
    StringStream ss;
    ss << "[";
    if (has_short_id()) ss << "-" << short_id() << " | ";
    ss << "--" << long_id() << "]\t" << instructions();
    return ss.str();
}

Bool ArgumentParserBase::is_consumable(const ArgumentStream &stream) const {
    auto arg = stream.peek();
    String short_argument = "-"+short_id();
    String long_argument = "--"+long_id();
    return (arg == short_argument or arg == long_argument);
}

class VerbosityArgumentParser : public ArgumentParserBase {
  public:
    VerbosityArgumentParser() : ArgumentParserBase(
            "v","verbosity","Choose the logging verbosity as a non-negative integer <value> (default: 0)") { }

    Void consume(ArgumentStream& stream) const override {
        stream.pop(); // Pop out the identifier
        if (stream.empty()) { throw MissingArgumentValueException(); }
        try {
            int val = std::stoi(stream.pop());
            if (val < 0) throw InvalidArgumentValueException("Verbosity should be a non-negative value.");
            Logger::instance().configuration().set_verbosity(static_cast<unsigned int>(val));
        } catch (...) {
            throw InvalidArgumentValueException("Verbosity should be a non-negative integer value.");
        }
    }
};

class ConcurrencyArgumentParser : public ArgumentParserBase {
public:
    ConcurrencyArgumentParser() : ArgumentParserBase(
            "c","concurrency","Choose the concurrency as a non-negative integer <value> or use 'max' to choose the hardware concurrency of this machine (default: 0)") { }

    Void consume(ArgumentStream& stream) const override {
        stream.pop(); // Pop out the identifier
        if (stream.empty()) { throw MissingArgumentValueException(); }
        try {
            String arg = stream.pop();
            if (arg == "max") TaskManager::instance().set_maximum_concurrency();
            else {
                int val = std::stoi(arg);
                if (val < 0) throw InvalidArgumentValueException("Concurrency should be a non-negative value.");
                TaskManager::instance().set_concurrency(static_cast<unsigned int>(val));
            }
        } catch (...) {
            throw InvalidArgumentValueException("Concurrency should be a non-negative integer value.");
        }
    }
};

CommandLineInterface::CommandLineInterface() : _parsers({VerbosityArgumentParser(),ConcurrencyArgumentParser()}) { }

Bool CommandLineInterface::acquire(int argc, const char* argv[]) const {
    ArgumentStream stream(argc,argv);
    while (not stream.empty()) {
        Bool consumed = false;
        for (auto& parser : _parsers) {
            if (parser.is_consumable(stream)) {
                consumed = true;
                try {
                    parser.consume(stream);
                } catch(MissingArgumentValueException& exc) {
                    std::cout << "A value is required by the argument, but it is not supplied." << std::endl;
                    return false;
                } catch(InvalidArgumentValueException& exc) {
                    std::cout << exc.what() << std::endl;
                    return false;
                }
                break;
            }
        }
        if (not consumed) {
            std::cout << "Unrecognised command-line option \"" << stream.peek() << "\"" << std::endl;
            return false;
        }
    }
    return true;
}

} // namespace Ariadne