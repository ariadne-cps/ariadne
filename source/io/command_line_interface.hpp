/***************************************************************************
 *            io/command_line_interface.hpp
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

/*! \file io/command_line_interface.hpp
 *  \brief Support for interacting with the command line.
 */

#ifndef ARIADNE_COMMAND_LINE_INTERFACE_HPP
#define ARIADNE_COMMAND_LINE_INTERFACE_HPP

#include <queue>
#include "utility/macros.hpp"
#include "utility/handle.hpp"
#include "utility/container.hpp"
#include "utility/string.hpp"

namespace Ariadne {

//! \brief Exception for when the argument is not recognised by any parser
class UnrecognisedArgumentException : public std::exception { };
//! \brief Exception for when no value is available, but it should be supplied
class MissingArgumentValueException : public std::exception { };
//! \brief Exception for when a parser recognises the value argument as invalid
class InvalidArgumentValueException : public std::runtime_error {
  public:
    InvalidArgumentValueException(String what) : std::runtime_error(what) { }
};

//! \brief Stream of arguments to be consumed by parsers
class ArgumentStream {
  public:
    //! \brief Construct from CLI arguments
    ArgumentStream(int argc, const char* argv[]);

    //! \brief Peek the head of the stream
    String peek() const;

    //! \brief Extract the head of the stream
    String pop();

    //! \brief If the stream is empty
    Bool empty() const;

    //! \brief Number of the elements in the stream
    SizeType size() const;

  private:
    std::queue<String> _args;
};

//! \brief Interface for a parser for a given CLI argument
class ArgumentParserInterface {
  public:
    //! \brief If the tip of the stream is consumable by this parser
    //! \details Checks only the kind, not the value
    virtual Bool is_consumable(ArgumentStream const& stream) const = 0;

    //! \brief Consume the stream
    //! \throws InvalidArgumentValueException if the value is incorrect
    //! \details Assumes that it has already been checked by is_consumable
    virtual Void consume(ArgumentStream& stream) const = 0;

    //! \brief The description to print for the help
    virtual String help_description() const = 0;
};

//! \brief An argument parser
class ArgumentParser : public Handle<ArgumentParserInterface> {
  public:
    using Handle<ArgumentParserInterface>::Handle;
    Bool is_consumable(ArgumentStream const& stream) const { return this->_ptr->is_consumable(stream); }
    Void consume(ArgumentStream& stream) const { this->_ptr->consume(stream); }
    String help_description() const { return this->_ptr->help_description(); }
};

//! \brief A static class for acquisition of CLI arguments
class CommandLineInterface {
  private:
    CommandLineInterface();
    void _print_help() const;
  public:
    CommandLineInterface(CommandLineInterface const&) = delete;
    Void operator=(CommandLineInterface const&) = delete;

    static CommandLineInterface& instance() {
        static CommandLineInterface instance;
        return instance;
    }

    //! \brief Acquire the CLI arguments
    //! \details Returns whether the acquisition was a success and
    //! any subsequent code should be run
    Bool acquire(int argc, const char* argv[]) const;

  private:
    List<ArgumentParser> const _parsers;
};

} // namespace Ariadne

#endif // ARIADNE_COMMAND_LINE_INTERFACE_HPP
