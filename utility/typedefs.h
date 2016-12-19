/***************************************************************************
 *            utility/typedefs.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file utility/typedefs.h
 *  \brief
 */



#ifndef ARIADNE_TYPEDEFS_H
#define ARIADNE_TYPEDEFS_H

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
