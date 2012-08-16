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


using std::pair;
using std::make_pair;

// A pair of references, suitable as use as an lvalue for a function returning a pair.
template<class T1, class T2>
struct lpair
{
    inline lpair(T1& t1, T2& t2) : first(t1), second(t2) { }
    inline lpair<T1,T2> operator=(const std::pair<T1,T2>& rv) {
        this->first=rv.first; this->second=rv.second; return *this; }
    T1& first; T2& second;
};

template<class T1,class T2> inline
lpair<T1,T2> make_lpair(T1& t1, T2& t2) {
    return lpair<T1,T2>(t1,t2);
}


/*! \brief A Tuple of possibly different types. */
template<class T1, class T2=void, class T3=void, class T4=void> struct Tuple;

/*! \brief A Tuple of references, suitable as use as an lvalue for a function returning a pair. */
template<class T1, class T2=void, class T3=void, class T4=void> struct LTuple;


template<class T1>
struct Tuple<T1,void,void,void>
{
    inline Tuple(const T1& t1) : first(t1) { }
    inline Tuple(const std::tuple<T1>& t) : first(std::get<0>(t)) { }
    T1 first;
};

template<class T1, class T2>
struct Tuple<T1,T2,void,void>
{
    inline Tuple(const T1& t1, const T2& t2) : first(t1), second(t2) { }
    inline Tuple(const std::tuple<T1,T2>& t) : first(std::get<0>(t)), second(std::get<1>(t)) { }
    T1 first; T2 second;
};

template<class T1, class T2, class T3>
struct Tuple<T1,T2,T3,void>
{
    inline Tuple(const T1& t1, const T2& t2, const T3& t3) : first(t1), second(t2), third(t3) { }
    inline Tuple(const std::tuple<T1,T2,T3>& t) : first(std::get<0>(t)), second(std::get<1>(t)), third(std::get<2>(t)) { }
    T1 first; T2 second; T3 third;
};

template<class T1, class T2, class T3, class T4>
struct Tuple
{
    inline Tuple(const T1& t1, const T2& t2, const T3& t3, const T4& t4)
        : first(t1), second(t2), third(t3), fourth(t4) { }
    inline Tuple(const std::tuple<T1,T2,T3,T4>& t)
        : first(std::get<0>(t)), second(std::get<1>(t)), third(std::get<2>(t)), fourth(std::get<3>(t)) { }
    T1 first; T2 second; T3 third; T4 fourth;
};



template<class T1>
struct LTuple<T1,void,void,void>
{
    inline LTuple(T1& t1) : first(t1) { }
    inline LTuple<T1> operator=(const Tuple<T1>& rv) {
        this->first=rv.first; return *this; }
    inline LTuple<T1> operator=(const T1& rv) {
        this->first=rv; return *this; }
    T1& first;
};

template<class T1, class T2>
struct LTuple<T1,T2,void,void>
{
    inline LTuple(T1& t1, T2& t2) : first(t1), second(t2) { }
    inline LTuple<T1,T2> operator=(const Tuple<T1,T2>& rv) {
        this->first=rv.first; this->second=rv.second; return *this; }
    inline LTuple<T1,T2> operator=(const std::pair<T1,T2>& rv) {
        this->first=rv.first; this->second=rv.second; return *this; }
    T1& first; T2& second;
};

template<class T1, class T2, class T3>
struct LTuple<T1,T2,T3,void>
{
    inline LTuple(T1& t1, T2& t2, T3& t3) : first(t1), second(t2), third(t3) { }
    inline LTuple<T1,T2,T3> operator=(const Tuple<T1,T2,T3>& rv) {
        this->first=rv.first; this->second=rv.second; this->third=rv.third; return *this; }
    T1& first; T2& second; T3& third;
};

template<class T1, class T2, class T3, class T4>
struct LTuple
{
    inline LTuple(T1& t1, T2& t2, T3& t3, T4& t4) : first(t1), second(t2), third(t3), fourth(t4) { }
    inline LTuple<T1,T2,T3,T4> operator=(const Tuple<T1,T2,T3,T4>& rv) {
        this->first=rv.first; this->second=rv.second; this->third=rv.third; this->fourth=rv.fourth; return *this; }
    T1& first; T2& second; T3& third; T4& fourth;
};



template<class T1> inline
Tuple<T1> make_tuple(const T1& t1) {
    return Tuple<T1>(t1);
}

template<class T1,class T2> inline
Tuple<T1,T2> make_tuple(const T1& t1, const T2& t2) {
    return Tuple<T1,T2>(t1,t2);
}

template<class T1,class T2,class T3> inline
Tuple<T1,T2,T3> make_tuple(const T1& t1, const T2& t2, const T3& t3) {
    return Tuple<T1,T2,T3>(t1,t2,t3);
}

template<class T1,class T2,class T3,class T4> inline
Tuple<T1,T2,T3,T4> make_tuple(const T1& t1, const T2& t2, const T3& t3, const T4& t4) {
    return Tuple<T1,T2,T3,T4>(t1,t2,t3,t4);
}


template<class T1> inline
LTuple<T1> make_ltuple(T1& t1) {
    return LTuple<T1>(t1);
}

template<class T1,class T2> inline
LTuple<T1,T2> make_ltuple(T1& t1, T2& t2) {
    return LTuple<T1,T2>(t1,t2);
}

template<class T1,class T2,class T3> inline
LTuple<T1,T2,T3> make_ltuple(T1& t1, T2& t2, T3& t3) {
    return LTuple<T1,T2,T3>(t1,t2,t3);
}

template<class T1,class T2,class T3,class T4> inline
LTuple<T1,T2,T3,T4> make_ltuple(T1& t1, T2& t2, T3& t3, T4& t4) {
    return LTuple<T1,T2,T3,T4>(t1,t2,t3,t4);
}




} // namespace Ariadne

#endif /* ARIADNE_TUPLE_H */
