/***************************************************************************
 *            utility/container.h
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

/*! \file utility/container.h
 *  \brief
 */



#ifndef ARIADNE_CONTAINER_H
#define ARIADNE_CONTAINER_H

#include "stdlib.h"

#include "metaprogramming.h"
#include "array.h"

namespace Ariadne {

using std::make_tuple;
using std::make_pair;

using SizeType=std::size_t;

template<class T> std::ostream& operator<<(std::ostream& os, const Array<T>& a) {
    bool first=true;
    for(auto x : a) {
        os << (first ? "[" : ",") << x;
        first = false;
    }
    if(first) { os << "["; }
    return os << "]";
}

template<class T> std::ostream& operator<<(std::ostream& os, const SharedArray<T>& a) {
    bool first=true;
    for(auto x : a) {
        os << (first ? "[" : ",") << x;
        first = false;
    }
    if(first) { os << "["; }
    return os << "]";
}

template<class T> class List : public std::vector<T> {
  public:
    typedef typename std::vector<T>::iterator Iterator;
    typedef typename std::vector<T>::const_iterator ConstIterator;
    using std::vector<T>::vector;
    List() : std::vector<T>() { }
    template<class TT, EnableIf<IsConvertible<TT,T>> =dummy>
        List(const List<TT>& l) : std::vector<T>(l.begin(),l.end()) { }
    Void append(const T& t) { this->push_back(t); }
};
template<class T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    bool first=true;
    for(auto x : v) {
        os << (first ? "[" : ",") << x;
        first = false;
    }
    if(first) { os << "["; }
    return os << "]";
}

template<class T> class Set : public std::set<T> {
  public:
    using std::set<T>::set;
    Set() : std::set<T>() { }
    template<class TT, EnableIf<IsConvertible<TT,T>> =dummy>
        Set(const std::initializer_list<TT>& s) : std::set<T>(s.begin(),s.end()) { }
    template<class TT, EnableIf<IsConvertible<TT,T>> =dummy>
        Set(const std::set<TT>& s) : std::set<T>(s.begin(),s.end()) { }
    bool contains(const T& t) const { return this->find(t)!=this->end(); }
};
template<class T> Set<T> join(Set<T> s1, Set<T> const& s2) {
    for(auto x2 : s2) { s1.insert(x2); } return std::move(s1); }

template<class T> std::ostream& operator<<(std::ostream& os, const std::set<T>& v) {
    bool first=true;
    for(auto x : v) {
        os << (first ? "{" : ",") << x;
        first = false;
    }
    if(first) { os << "{"; }
    return os << "}";
}


template<class K, class T> class Map : public std::map<K,T> {
  public:
    using std::map<K,T>::map;
    using std::map<K,T>::insert;
    void insert(const K& k, const T& t) { this->insert(std::make_pair(k,t)); }
    T& operator[](K k) { return this->std::map<K,T>::operator[](k); }
    const T& operator[](K k) const { auto iter=this->find(k); assert(iter!=this->end()); return iter->second; }
    bool has_key(const K& k) const { return this->find(k)!=this->end(); }
};
template<class K, class T> std::ostream& operator<<(std::ostream& os, const std::map<K,T>& m) {
    bool first=true;
    for(auto x : m) {
        os << (first ? "{ " : ", ") << x.first << ":" << x.second;
        first = false;
    }
    if(first) { os << "{"; }
    return os << " }";
}


} // namespace Ariadne

#endif
