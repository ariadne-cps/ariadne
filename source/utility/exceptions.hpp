/***************************************************************************
 *            exceptions.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file exceptions.hpp
 *  \brief Exception classes
 */

#ifndef ARIADNE_EXCEPTIONS_HPP
#define ARIADNE_EXCEPTIONS_HPP

#include "../utility/string.hpp"

#include <exception>
#include <stdexcept>

namespace Ariadne {

class Exception {
    String _what;
  public:
    Exception(const String& what);
};

class NotImplemented : public std::logic_error {
  public:
    NotImplemented(const String& str) : std::logic_error(str) { }
};

class DivideByZeroException : public std::runtime_error {
  public:
    DivideByZeroException(const String& str) : std::runtime_error(str) { }
};

class DomainException : public std::runtime_error {
  public:
    DomainException(const String& str) : std::runtime_error(str) { }
};

class IncompatibleSizes : public std::runtime_error {
  public:
    IncompatibleSizes(const String& str) : std::runtime_error(str) { }
};

class InvalidGridPosition : public std::runtime_error {
  public:
    InvalidGridPosition(const String& str) : std::runtime_error(str) { }
};

class InvalidCoordinate : public std::runtime_error {
  public:
    InvalidCoordinate(const String& str) : std::runtime_error(str) { }
};

class InvalidInput : public std::runtime_error {
  public:
    InvalidInput(const String& str) : std::runtime_error(str) { }
};

class OuterChainOverspill : public std::runtime_error {
  public:
    OuterChainOverspill(const String& str) : std::runtime_error(str) { }
};

} // namespace Ariadne

#endif
