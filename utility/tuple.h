/***************************************************************************
 *            tuple.h
 *
 *  Copyright  2007  Pieter Collins
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

/*! \file tuple.h
 *  \brief Pair and Tuple types, and types to be used as lvalues in assignments.
 */

#ifndef ARIADNE_TUPLE_H
#define ARIADNE_TUPLE_H

#include <utility>
#include <tuple>

namespace Ariadne {

template<class T1, class T2> using Pair = std::pair<T1,T2>;
using std::make_pair;
template<class T1, class T2> constexpr Pair<T1&,T2&> make_lpair(T1& t1, T2& t2) { return Pair<T1&,T2&>(t1,t2); }

template<class... TS> using Tuple = std::tuple<TS...>;
using std::make_tuple;
template<class... TS> constexpr Tuple<TS&...> make_ltuple(TS&... ts) { return Tuple<TS&...>(std::tie(ts...)); }


template<class T> inline decltype(auto) get_first(T&& t) { return std::get<0>(std::forward<T>(t)); }
template<class T> inline decltype(auto) get_second(T&& t) { return std::get<1>(std::forward<T>(t)); }
template<class T> inline decltype(auto) get_third(T&& t) { return std::get<2>(std::forward<T>(t)); }
template<class T> inline decltype(auto) get_fourth(T&& t) { return std::get<3>(std::forward<T>(t)); }
template<class T> inline decltype(auto) get_fifth(T&& t) { return std::get<4>(std::forward<T>(t)); }


} // namespace Ariadne

#endif /* ARIADNE_TUPLE_H */
