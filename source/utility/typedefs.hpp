/***************************************************************************
 *            utility/typedefs.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file utility/typedefs.hpp
 *  \brief
 */



#ifndef ARIADNE_TYPEDEFS_HPP
#define ARIADNE_TYPEDEFS_HPP

#include <cstdint>
#include <iosfwd>
#include <type_traits>
#include <memory>

typedef unsigned char uchar;
typedef unsigned int uint;

namespace Ariadne {

using OutputStream = std::ostream;
using InputStream = std::istream;
using StringStream = std::stringstream;

class String; // Define as a class for consistency with other value types
using StringType = std::string;

using Void = void;
using Char = char;
using Byte = std::int8_t;
using Bool = bool;
using Nat = uint;
using Int = int;
using Dbl = double;
using Nat32Type = std::uint32_t;
using Int32Type = std::int32_t;
using Nat64Type = std::uint64_t;
using Int64Type = std::int64_t;

using std::declval;

using SizeType = std::size_t;
using PointerDifferenceType = std::ptrdiff_t;
using DegreeType = std::uint16_t;

template<class T> using SharedPointer = std::shared_ptr<T>;
template<class T> using InitializerList = std::initializer_list<T>;
template<class T1, class T2> using Pair = std::pair<T1,T2>;
template<class... TS> using Tuple = std::tuple<TS...>;

template<class T> class Array;
template<class T> class List;
template<class T> class Set;
template<class K, class T> class Map;

} // namespace Ariadne

#endif
