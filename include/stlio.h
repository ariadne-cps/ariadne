/***************************************************************************
 *            stlio.h
 *
 *  Copyright  2005-8  Pieter Collins
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

/*! \file stlio.h
 *  \brief Input-output utilities
 */

#ifndef ARIADNE_STLIO_H
#define ARIADNE_STLIO_H

#include <iostream>
#include <stdexcept>

#include <string>
#include <vector>
#include <list>
#include <deque>
#include <valarray>
#include <set>
#include <map>
#include "array.h"
#include "tuple.h"

#include "boost/shared_ptr.hpp"


namespace Ariadne {

template<class InputIterator>
std::ostream&
write_sequence(std::ostream& os, InputIterator first, InputIterator last,
               char opening='[', char closing=']', char separator=',')
{
    os << opening;
    while(first!=last) {
        os << (*first);
        ++first;
        if(first!=last) {
            os << separator;
        }
    }
    os << closing;
    return os;
}


template<class InputIterator>
std::ostream&
write_pointer_sequence(std::ostream& os, InputIterator first, InputIterator last,
                       char opening='[', char closing=']', char separator=',')
{
    os << opening;
    while(first!=last) {
        os << (**first);
        ++first;
        if(first!=last) {
            os << separator;
        }
    }
    os << closing;
    return os;
}


template<class InputIterator>
std::ostream&
write_map_sequence(std::ostream& os, InputIterator first, InputIterator last,
                   char opening='{', char closing='}', char separator=',', char descriptor=':')
{
    os << opening;
    while(first!=last) {
        os << first->first << descriptor << first->second;
        ++first;
        if(first != last) {
            os << separator;
        }
    }
    os << closing;
    return os;
}

template<class InputIterator>
std::ostream&
write_map_pointer_sequence(std::ostream& os, InputIterator first, InputIterator last,
                           char opening='{', char closing='}', char separator=',', char descriptor=':')
{
    os << opening;
    while(first!=last) {
        os << first->first << descriptor << *first->second;
        ++first;
        if(first != last) {
            os << separator;
        }
    }
    os << closing;
    return os;
}

template<class InputIterator>
std::ostream&
write_ariadne_map_sequence(std::ostream& os, InputIterator first, InputIterator last,
                           char opening='{', char closing='}', char separator=',', char descriptor=':')
{
    os << opening;
    while(first!=last) {
        os << first->key() << descriptor << first->data();
        ++first;
        if(first != last) {
            os << separator;
        }
    }
    os << closing;
    return os;
}


template<class Container>
std::istream&
read_sequence(std::istream& is, Container& v,
              char opening='[', char closing=']', char separator=',')
{
    typedef typename Container::value_type T;

    T x;
    char c;

    v.clear();

    try {
        is >> c;
        if(c != opening) {
            throw std::ios_base::failure(std::string("Ariadne::Base::read_vector: Input must begin with ")+opening);
        }

        /* Handle case of empty list */
        is >> c;
        if(c != closing) {
            is.putback(c);
            c=separator;
        }

        while(c != closing) {
            if(is.eof()) {
                throw std::ios_base::failure("Ariadne::Base::read_vector: End-of-file reached");
            }
            if(c!=separator) {
                throw std::ios_base::failure(std::string("Ariadne::Base::read_vector: Items in list must be separated by ")+separator);
            }
            is >> x;
            if(is.fail()) {
                throw std::ios_base::failure("Ariadne::Base::read_vector: Error inputting value in list");
            }
            v.push_back(x);
            is >> c;
        }
    }
    catch(...) {
        // is.seekg(pos);
        throw;
    }

    return is;
}



template<class T> inline
std::ostream& operator<<(std::ostream& os, const Array<T>& ary) {
    return Ariadne::write_sequence(os,ary.begin(),ary.end(),'[',']');
}


template<class T> inline
std::istream&
operator>>(std::istream& is, Array<T>& a) {
    std::vector<T> v;
    read_sequence(is,v);
    a=Array<T>(v.begin(),v.end());
    return is;
}

template<class T1>
inline
std::ostream&
operator<<(std::ostream &os, const Tuple<T1>& t)
{
    return os << '(' << t.first << ',' << ')';
}

template<class T1, class T2>
inline
std::ostream&
operator<<(std::ostream &os, const Tuple<T1,T2>& t)
{
    return os << '(' << t.first << ',' << t.second << ')';
}

template<class T1, class T2, class T3>
inline
std::ostream&
operator<<(std::ostream &os, const Tuple<T1,T2,T3>& t)
{
    return os << '(' << t.first << ',' << t.second << ',' << t.third << ')';
}

template<class T1, class T2, class T3, class T4>
inline
std::ostream&
operator<<(std::ostream &os, const Tuple<T1,T2,T3,T4>& t)
{
    return os << '(' << t.first << ',' << t.second << ',' << t.third << ',' << t.fourth << ')';
}


} // namespace Ariadne


/* FIXME: This is a hack to allow io of STL classes.
   But really we should not modify namespace std.
   Unfortunately, we need to include the code in
   any namespace using operator<<.
*/
namespace std {

template<class S, class T>
inline
std::ostream&
operator<<(std::ostream &os, const std::pair<S,T>& s)
{
    return os << '(' << s.first << ',' << s.second << ')';
}

template<class T>
inline
std::ostream&
operator<< (std::ostream &os, const std::vector<T>& v)
{
    return Ariadne::write_sequence(os,v.begin(),v.end());
}

template<class T>
inline
std::ostream&
operator<< (std::ostream &os, const std::list<T>& l)
{
    return Ariadne::write_sequence(os,l.begin(),l.end());
}

template<class T>
inline
std::ostream&
operator<< (std::ostream &os, const std::deque<T>& d)
{
    return Ariadne::write_sequence(os,d.begin(),d.end());
}

template<class T>
inline
ostream&
operator<< (std::ostream &os, const std::valarray<T>& v) {
    return Ariadne::write_sequence(os,&(v[0]),&(v[v.size()-1]));
}

template<class T, class C>
inline
std::ostream&
operator<<(std::ostream &os, const std::set<T,C>& s)
{
    return Ariadne::write_sequence(os,s.begin(), s.end(), '{', '}');
}

template<class K, class T, class C>
inline
std::ostream&
operator<<(std::ostream &os, const std::map<K,T,C>& m)
{
    return Ariadne::write_map_sequence(os,m.begin(), m.end(), '{', '}');
}

template<class K, class T, class C>
inline
std::ostream&
operator<<(std::ostream &os, const std::map<K,boost::shared_ptr<T>,C>& m)
{
    return Ariadne::write_map_pointer_sequence(os,m.begin(), m.end(), '{', '}');
}


template<class T>
inline
istream&
operator>> (std::istream &is, std::vector<T>& v) {
    return Ariadne::read_sequence(is,v);
}

} // namespace std


#endif /* ARIADNE_STLIO_H */
