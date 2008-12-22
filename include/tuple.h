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
 *  \brief Pair and tuple types, and types to be used as lvalues in assignments.
 */

#ifndef ARIADNE_TUPLE_H
#define ARIADNE_TUPLE_H

#include <utility>

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
  

/*! \brief A tuple of possibly different types. */
template<class T1, class T2, class T3=void, class T4=void> struct tuple;
  
/*! \brief A tuple of references, suitable as use as an lvalue for a function returning a pair. */
template<class T1, class T2, class T3=void, class T4=void> struct ltuple;
  
    
template<class T1, class T2> 
struct tuple<T1,T2,void,void>
{
    inline tuple(const T1& t1, const T2& t2) : first(t1), second(t2) { }
    T1 first; T2 second;
};

template<class T1, class T2, class T3> 
struct tuple<T1,T2,T3,void>
{
    inline tuple(const T1& t1, const T2& t2, const T3& t3) : first(t1), second(t2), third(t3) { }
    T1 first; T2 second; T3 third;
};

template<class T1, class T2, class T3, class T4> 
struct tuple
{
    inline tuple(const T1& t1, const T2& t2, const T3& t3, const T4& t4)
        : first(t1), second(t2), third(t3), fourth(t4) { }
    T1 first; T2 second; T3 third; T4 fourth;
};



template<class T1, class T2> 
struct ltuple<T1,T2,void,void>
{
    inline ltuple(T1& t1, T2& t2) : first(t1), second(t2) { }
    inline ltuple<T1,T2> operator=(const tuple<T1,T2>& rv) { 
        this->first=rv.first; this->second=rv.second; return *this; }
    inline ltuple<T1,T2> operator=(const std::pair<T1,T2>& rv) { 
        this->first=rv.first; this->second=rv.second; return *this; }
    T1& first; T2& second;
};

template<class T1, class T2, class T3> 
struct ltuple<T1,T2,T3,void>
{
    inline ltuple(T1& t1, T2& t2, T3& t3) : first(t1), second(t2), third(t3) { }
    inline ltuple<T1,T2,T3> operator=(const tuple<T1,T2,T3>& rv) { 
        this->first=rv.first; this->second=rv.second; this->third=rv.third; return *this; }
    T1& first; T2& second; T3& third;
};

template<class T1, class T2, class T3, class T4> 
struct ltuple
{
    inline ltuple(T1& t1, T2& t2, T3& t3, T4& t4) : first(t1), second(t2), third(t3), fourth(t4) { }
    inline ltuple<T1,T2,T3,T4> operator=(const tuple<T1,T2,T3,T4>& rv) { 
        this->first=rv.first; this->second=rv.second; this->third=rv.third; this->fourth=rv.fourth; return *this; }
    T1& first; T2& second; T3& third; T4& fourth;
};
  


template<class T1,class T2> inline
tuple<T1,T2> make_tuple(const T1& t1, const T2& t2) {
    return tuple<T1,T2>(t1,t2);
}

template<class T1,class T2,class T3> inline
tuple<T1,T2,T3> make_tuple(const T1& t1, const T2& t2, const T3& t3) {
    return tuple<T1,T2,T3>(t1,t2,t3);
}

template<class T1,class T2,class T3,class T4> inline
tuple<T1,T2,T3,T4> make_tuple(const T1& t1, const T2& t2, const T3& t3, const T4& t4) {
    return tuple<T1,T2,T3,T4>(t1,t2,t3,t4);
}


template<class T1,class T2> inline
ltuple<T1,T2> make_ltuple(T1& t1, T2& t2) {
    return ltuple<T1,T2>(t1,t2);
}

template<class T1,class T2,class T3> inline
ltuple<T1,T2,T3> make_ltuple(T1& t1, T2& t2, T3& t3) {
    return ltuple<T1,T2,T3>(t1,t2,t3);
}

template<class T1,class T2,class T3,class T4> inline
ltuple<T1,T2,T3,T4> make_ltuple(T1& t1, T2& t2, T3& t3, T4& t4) {
    return ltuple<T1,T2,T3,T4>(t1,t2,t3,t4);
}



  
} // namespace Ariadne

#endif /* ARIADNE_TUPLE_H */
