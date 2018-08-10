/***************************************************************************
 *            string.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file string.hpp
 *  \brief Wrapper for string class
 */

#ifndef ARIADNE_STRING_HPP
#define ARIADNE_STRING_HPP

#include <string>
#include <sstream>

namespace Ariadne {

class String : public std::string {
  public:
    using std::string::string;
    String(std::string const& str) : std::string(str) { };
    String() = default;
    String(String const&) = default;
};

inline const char* c_str(const String& str) {
    return str.c_str(); }
template<class T> inline String to_string(const T& t) {
    std::stringstream ss; ss << t; return ss.str(); }
template<class T> inline String to_str(T const& t) {
    return to_string(t); }

template<class T> String class_name();
template<> inline String class_name<String>() { return "String"; }

} // namespace Ariadne

#endif /* ARIADNE_STRING_HPP */
