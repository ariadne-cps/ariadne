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

ArgumentPack::ArgumentPack(String const& id, VoidFunction const& processor) : _id(id), _processor(processor) { }

String const& ArgumentPack::id() const {
    return _id;
}

Void ArgumentPack::process() const {
    _processor();
}

Bool ArgumentPack::operator<(ArgumentPack const& other) const {
    return this->_id < other._id;
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

    ArgumentPack consume(ArgumentStream& stream) const override;

    //! \brief If a value is required for this argument
    virtual Bool requires_value() const = 0;

    Bool is_consumable(ArgumentStream const& stream) const override;
    SizeType help_description_header_size() const override;
    String help_description(SizeType num_chars_to_align_instructions) const override;

    virtual ~ArgumentParserBase() = default;

  protected:
    virtual VoidFunction create_processor(ArgumentStream& stream) const = 0;

  private:
    String _short_id;
    String _long_id;
    String _instructions;
};

SizeType ArgumentParserBase::help_description_header_size() const {
    SizeType result = 4+long_id().size();
    if (has_short_id()) result+= 4+short_id().size();
    if (requires_value()) result += 8;
    return result;
}

String ArgumentParserBase::help_description(SizeType num_chars_to_separate_instructions) const {
    StringStream ss;
    ss << "[";
    if (has_short_id()) { ss << "-" << short_id() << " | "; }
    ss << "--" << long_id() << "]";
    if (requires_value()) { ss << " <value>"; }
    ss << std::string(num_chars_to_separate_instructions,' ') << instructions();
    return ss.str();
}

Bool ArgumentParserBase::is_consumable(const ArgumentStream &stream) const {
    auto arg = stream.peek();
    if (has_short_id()) {
        String short_argument = "-"+short_id();
        if (arg == short_argument) return true;
    }
    String long_argument = "--"+long_id();
    return (arg == long_argument);
}

ArgumentPack ArgumentParserBase::consume(ArgumentStream& stream) const {
    stream.pop(); // Pop out the identifier already recognised as the one for this parser
    return ArgumentPack(_long_id,create_processor(stream));
}

class ValuedArgumentParserBase : public ArgumentParserBase {
  public:
    ValuedArgumentParserBase(String const& s, String const& l, String const& i) : ArgumentParserBase(s,l,i) { }
    Bool requires_value() const override { return true; }
};

class UnvaluedArgumentParserBase : public ArgumentParserBase {
  public:
    UnvaluedArgumentParserBase(String const& s, String const& l, String const& i) : ArgumentParserBase(s,l,i) { }
    Bool requires_value() const override { return false; }
};

class HelpArgumentParser : public UnvaluedArgumentParserBase {
public:
    HelpArgumentParser() : UnvaluedArgumentParserBase(
            "h","help","Show this list of supported arguments.") { }

    VoidFunction create_processor(ArgumentStream& stream) const override {
        return []{};
    }
};

class VerbosityArgumentParser : public ValuedArgumentParserBase {
  public:
    VerbosityArgumentParser() : ValuedArgumentParserBase(
            "v","verbosity","Choose the logging verbosity as a non-negative integer <value> (default: 0)") { }

    VoidFunction create_processor(ArgumentStream& stream) const override {
        VoidFunction f;
        if (stream.empty()) { throw MissingArgumentValueException(); }
        try {
            int val = std::stoi(stream.pop());
            if (val < 0) throw InvalidArgumentValueException("Verbosity should be a non-negative value.");
            f = [val]{ Logger::instance().configuration().set_verbosity(static_cast<unsigned int>(val)); };
        } catch (...) {
            throw InvalidArgumentValueException("Verbosity should be a non-negative integer value.");
        }
        return f;
    }
};

class ConcurrencyArgumentParser : public ValuedArgumentParserBase {
public:
    ConcurrencyArgumentParser() : ValuedArgumentParserBase(
            "c","concurrency","Choose the concurrency as a non-negative integer <value> or use 'max' to choose the hardware concurrency of this machine (default: 0)") { }

    VoidFunction create_processor(ArgumentStream& stream) const override {
        VoidFunction f;
        if (stream.empty()) { throw MissingArgumentValueException(); }
        try {
            String arg = stream.pop();
            if (arg == "max") f = []{ TaskManager::instance().set_maximum_concurrency(); };
            else {
                int val = std::stoi(arg);
                if (val < 0) throw std::exception();
                f = [val]{ TaskManager::instance().set_concurrency(static_cast<unsigned int>(val)); };
            }
        } catch (...) {
            throw InvalidArgumentValueException("Concurrency should be a non-negative integer value.");
        }
        return f;
    }
};

CommandLineInterface::CommandLineInterface() : _parsers({ConcurrencyArgumentParser(),HelpArgumentParser(),VerbosityArgumentParser()}) { }

Bool CommandLineInterface::acquire(int argc, const char* argv[]) const {
    ArgumentStream stream(argc,argv);
    Set<ArgumentPack> packs;
    while (not stream.empty()) {
        Bool consumed = false;
        for (auto& parser : _parsers) {
            if (parser.is_consumable(stream)) {
                consumed = true;
                try {
                    auto p = parser.consume(stream);
                    if (packs.contains(p)) {
                        std::clog << "Argument '" << p.id() << "' specified multiple times.\n\n"
                            << parser.help_description(4) << std::endl;
                    } else if (p.id() == "help") {
                        _print_help();
                        return false;
                    } else if (p.id() == "concurrency") {
                        // Need to process this before the others
                        p.process();
                    } else packs.insert(p);
                } catch(MissingArgumentValueException& exc) {
                    std::clog << "A value is required by the argument, but it is not supplied.\n\n"
                        << parser.help_description(4) << std::endl;
                    return false;
                } catch(InvalidArgumentValueException& exc) {
                    std::clog << exc.what() << "\n\n" << parser.help_description(4) << std::endl;
                    return false;
                }
                break;
            }
        }
        if (not consumed) {
            std::clog << "Unrecognised command-line argument '" << stream.peek() << "'" << std::endl << std::endl;
            _print_help();
            return false;
        }
    }

    for (auto p : packs) p.process();

    return true;
}

void CommandLineInterface::_print_help() const {

    List<SizeType> chars(_parsers.size(),0);
    SizeType max_chars = 0;
    for (SizeType i=0; i<_parsers.size(); ++i)  {
        chars[i] = _parsers[i].help_description_header_size();
        if (chars[i] > max_chars) max_chars = chars[i];
    }

    std::clog << "Supported arguments:" << std::endl;
    for (SizeType i=0; i<_parsers.size(); ++i) {
        std::clog << "    " << _parsers[i].help_description(4+max_chars-chars[i]) << std::endl;
    }
}

} // namespace Ariadne