/***************************************************************************
 *            container.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file container.h
 *  \brief Extensions to STL containers
 */

#ifndef ARIADNE_CONTAINER_H
#define ARIADNE_CONTAINER_H

#include <cassert>
#include <cstdarg>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include "array.h"
#include "stlio.h"

namespace Ariadne {


typedef std::string String;

inline std::vector<std::string> operator,(const std::string& s1, const char* s2) {
    std::vector<std::string> v; v.push_back(std::string(s1)); v.push_back(std::string(s2)); return v; }
inline std::vector<std::string> operator,(const std::vector<std::string>& v1, const char* s2) {
    std::vector<std::string> r(v1); r.push_back(std::string(s2)); return r; }


template<class T> inline bool contains(const std::set<T>& s, const T& t) {
    return s.find(t)!=s.end(); }
template<class T> inline bool subset(const std::set<T>& s1, const std::set<T>& s2) {
    for(typename std::set<T>::iterator iter=s1.begin(); iter!=s1.end(); ++iter) {
        if(!contains(s2,*iter)) { return false; } } return true; }
template<class T> inline bool disjoint(const std::set<T>& s1, const std::set<T>& s2) {
    for(typename std::set<T>::iterator iter=s1.begin(); iter!=s1.end(); ++iter) {
        if(contains(s2,*iter)) { return false; } } return true; }
template<class T> inline bool disjoint(const std::set<T>& s1, const std::vector<T>& l2) {
    for(typename std::vector<T>::const_iterator iter=l2.begin(); iter!=l2.end(); ++iter) {
        if(contains(s1,*iter)) { return false; } } return true; }
template<class T> inline std::set<T>& insert(std::set<T>& r, const T& t) {
    r.insert(t); return r; }
template<class T> inline std::set<T>& adjoin(std::set<T>& r, const std::set<T>& s) {
    for(typename std::set<T>::iterator iter=s.begin(); iter!=s.end(); ++iter) { r.insert(*iter); } return r; }
template<class T> inline std::set<T>& remove(std::set<T>& r, const std::set<T>& s) {
    typename std::set<T>::iterator iter=r.begin();
    while(iter!=r.end()) { if(contains(s,*iter)) { r.erase(iter++); } else { ++iter; } }
    return r; }
template<class T> inline std::set<T>& restrict(std::set<T>& r, const std::set<T>& s) {
    typename std::set<T>::iterator iter=r.begin();
    while(iter!=r.end()) { if(!contains(s,*iter)) { r.erase(iter++); } else { ++iter; } }
    return r; }
template<class T> inline std::set<T> join(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> r(s1); adjoin(r,s2); return r; }
template<class T> inline std::set<T> intersection(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> r(s1); restrict(r,s2); return r; }
template<class T> inline std::set<T> intersection(const std::set<T>& s1, const std::vector<T>& l2) {
    std::set<T> r; for(typename std::vector<T>::const_iterator iter=l2.begin(); iter!=l2.end(); ++iter) {
        if(contains(s1,*iter)) { r.insert(*iter); } } return r; }
template<class T> inline std::set<T> difference(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> r(s1); remove(r,s2); return r; }
template<class K,class V,class C> inline bool has_key(const std::map<K,V,C>& m, const K& k) {
    return m.find(k)!=m.end(); }


// Temporary home for useful auxiliary classes

template<class T> inline std::string to_string(const T& t) {
    std::stringstream ss; ss<<t; return ss.str(); }
template<> inline std::string to_string(const bool& t) {
    std::stringstream ss; ss<<std::boolalpha<<t; return ss.str(); }
template<class T> inline std::string to_str(const T& t) {
    return to_string(t); }

template<class T1, class T2> class Pair
    : public std::pair<T1,T2>
{
  public:
    Pair() : std::pair<T1,T2>() { }
    template<class X> Pair(const X& x) : std::pair<T1,T2>(x) { }
    template<class X1, class X2> Pair(const X1& x1, const X2& x2) : std::pair<T1,T2>(x1,x2) { }
};


template<class T> class List
    : public std::vector<T>
{
  public:
    List() : std::vector<T>() { }
    List(unsigned int n) : std::vector<T>(n) { }
    List(unsigned int n, const T& t) : std::vector<T>(n,t) { }
    List(const std::vector<T>& l) : std::vector<T>(l) { }
    List(std::initializer_list<T> l) : std::vector<T>(l) { }
    template<class X> List(const List<X>& l) : std::vector<T>(l.begin(),l.end()) { }
    template<class I> List(const I& b, const I& e) : std::vector<T>(b,e) { }
    // Conversion constructor from an value to a single-element list.
    List(const T& t) : std::vector<T>(1u,t) { }
    void append(const T& t) { this->push_back(t); }
    void append(const List<T>& t) { for(unsigned int i=0; i!=t.size(); ++i) { this->push_back(t[i]); } }
    void concatenate(const List<T>& t) { for(unsigned int i=0; i!=t.size(); ++i) { this->push_back(t[i]); } }
};

template<class T> inline List<T> catenate(const List<T>& l1, const List<T>& l2) {
    List<T> r(l1);
    for(typename List<T>::const_iterator iter=l2.begin(); iter!=l2.end(); ++iter) {
        r.append(*iter);
    }
    return r;
}

template<class T> inline List<T> catenate(const List<T>& l1, const T& t2) {
    List<T> r(l1);
    r.append(t2);
    return r;
}


template<class T>
inline List<T> operator,(const T& t1, const T& t2) {
    List<T> v; v.push_back(t1); v.push_back(t2); return v; }
template<class T>
inline List<T> operator,(const std::vector<T>& v, const T& t) {
    List<T> r(v); r.push_back(t); return r; }

template<class T> class LinkedList
    : public std::list<T>
{
  public:
    LinkedList() : std::list<T>() { }
    LinkedList(unsigned int n) : std::list<T>(n) { }
    LinkedList(unsigned int n, const T& t) : std::vector<T>(n,t) { }
    LinkedList(const std::list<T>& l) : std::list<T>(l) { }
    template<class X> LinkedList(const List<X>& l) : std::list<T>(l.begin(),l.end()) { }
    template<class X> LinkedList(const LinkedList<X>& l) : std::list<T>(l.begin(),l.end()) { }
    template<class I> LinkedList(const I& b, const I& e) : std::list<T>(b,e) { }
    // Conversion constructor from an value to a single-element list.
    void append(const T& t) { this->push_back(t); }
    void append(const LinkedList<T>& t) { for(unsigned int i=0; i!=t.size(); ++i) { this->push_back(t[i]); } }
    void concatenate(const LinkedList<T>& t) { for(unsigned int i=0; i!=t.size(); ++i) { this->push_back(t[i]); } }
};

template<class T> class Set
    : public std::set<T>
{
    template<class TT, class KK> class KeyEqual {
        bool operator()(const TT& t, const KK& k) { return t.key()==k; } };
  public:
    Set() : std::set<T>() { }
    Set(const std::set<T>& s) : std::set<T>(s) { }
    template<class TT> explicit Set(const std::set<TT>& s) {
        for(typename std::set<TT>::const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
            this->insert(static_cast<T>(*iter)); } }
    explicit Set(const std::vector<T>& l) : std::set<T>(l.begin(),l.end()) { }
    explicit Set(const T& t) { this->insert(t); }
    template<class I> Set(const I& b, const I& e) : std::set<T>(b,e) { }

    bool contains(const T& t) const {
        return this->find(t)!=this->end(); }
    bool subset(const std::set<T>& s) const {
        for(typename std::set<T>::iterator iter=s.begin(); iter!=s.end(); ++iter) {
            if(!this->contains(*iter)) { return false; } } return true; }
    bool disjoint(const std::set<T>& s) const {
        for(typename std::set<T>::const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
            if(this->contains(*iter)) { return false; } } return true; }
    bool disjoint(const std::vector<T>& lst) const {
        for(typename std::vector<T>::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter) {
            if(this->contains(*iter)) { return false; } } return true; }
    Set<T>& adjoin(const std::set<T>& s) {
        for(typename std::set<T>::const_iterator iter=s.begin(); iter!=s.end(); ++iter) { this->insert(*iter); } return *this; }
    Set<T>& remove(const T& t) {
        typename std::set<T>::iterator iter=this->find(t);
        if(iter!=this->end()) { this->erase(iter); }
        return *this; }
    Set<T>& remove(const std::set<T>& s) {
        typename std::set<T>::iterator iter=this->begin();
        while(iter!=this->end()) { if(s.find(*iter)!=s.end()) { this->erase(iter++); } else { ++iter; } }
        return *this; }
    Set<T>& restrict(const std::set<T>& s) {
        typename std::set<T>::iterator iter=this->begin();
        while(iter!=this->end()) { if(s.find(*iter)!=s.end()) { this->erase(iter++); } else { ++iter; } }
        return *this; }
    template<class TT> Set<T>& remove(const std::vector<TT>& l) {
        for(typename std::vector<TT>::const_iterator iter=l.begin(); iter!=l.end(); ++iter) {
            this->std::set<T>::erase(static_cast<const T&>(*iter)); }
        return *this; }
};

template<class T> inline Set<T> join(const Set<T>& s1, const Set<T>& s2) {
    Set<T> r(s1); adjoin(r,s2); return r; }

template<class T> inline Set<T> make_set(const std::vector<T>& lst) {
    return Set<T>(lst); }

template<class T> inline std::ostream& operator<<(std::ostream& os, const Set<T>& s) {
    return os<<static_cast<const std::set<T>&>(s); }


template<class K, class V> class Map
    : public std::map<K,V>
{
  public:
    Map<K,V>()
        : std::map<K,V>() { }
    Map<K,V>(const std::map<K,V>& m)
        : std::map<K,V>(m) { }
    Map<K,V>(const std::vector< std::pair<K,V> >& l)
        : std::map<K,V>() { for(size_t i=0; i!=l.size(); ++i) { this->insert(l[i]); } }
    template<class W> explicit Map<K,V>(const std::map<K,W>& m) {
        for(typename std::map<K,W>::const_iterator iter=m.begin(); iter!=m.end(); ++iter) {
            this->insert(iter->first,V(iter->second)); } }
    Map(std::initializer_list<std::pair<K,V>> l) : std::map<K,V>(l) { }
    V& operator[](const K& k) {
        return this->std::map<K,V>::operator[](k); }
    const V& operator[](const K& k) const { typename std::map<K,V>::const_iterator p=this->find(k);
        assert(p!=this->end()); return p->second; }
    const V& get(const K& k) const { typename std::map<K,V>::const_iterator p=this->find(k);
        assert(p!=this->end()); return p->second; }
    bool has_key(const K& k) const {
        return this->find(k)!=this->end(); }
    V& value(const K& k) {
        typename std::map<K,V>::iterator iter=this->find(k);
        assert(iter!=this->end()); return iter->second; }
    const V& value(const K& k) const {
        typename std::map<K,V>::const_iterator iter=this->find(k);
        assert(iter!=this->end()); return iter->second; }
    void insert(const std::pair<K,V>& kv) {
        this->std::map<K,V>::insert(kv); }
    void insert(const K& k, const V& v) {
        this->std::map<K,V>::insert(std::make_pair(k,v)); }
    void adjoin(const std::map<K,V>& m) {
        for(typename std::map<K,V>::const_iterator i=m.begin(); i!=m.end(); ++i) { this->insert(*i); } }
    void remove_keys(const Set<K>& s) {
        for(typename Set<K>::const_iterator iter=s.begin(); iter!=s.end(); ++iter) { this->erase(*iter); } }
    Set<K> keys() const {
        Set<K> res; for(typename std::map<K,V>::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
            res.insert(iter->first); } return res; }
    List<V> values() const {
        List<V> res; for(typename std::map<K,V>::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
            res.append(iter->second); } return res; }
    using std::map<K,V>::insert;
};
template<class I, class X, class J> inline X& insert(Map<I,X>& m, const J& k, const X& v) {
    return m.std::template map<I,X>::insert(std::make_pair(k,v)).first->second; }
template<class K,class V> inline Map<K,V> join(const Map<K,V>& m1, const Map<K,V>& m2) {
    Map<K,V> r(m1); for(typename Map<K,V>::const_iterator i=m2.begin(); i!=m2.end(); ++i) { r.insert(*i); } return r; }

template<class K,class V> inline Map<K,V> make_map(const std::vector< std::pair<K,V> >& lst) {
    return Map<K,V>(lst.begin(),lst.end()); }
template<class K,class V> inline Map<K,V> make_map(const std::vector< Pair<K,V> >& lst) {
    return Map<K,V>(lst.begin(),lst.end()); }

template<class K, class V> Map<K,V> restrict_keys(const std::map<K,V>& m, const std::set<K>& k) {
    Map<K,V> result;
    for(typename std::map<K,V>::const_iterator item_iter=m.begin(); item_iter!=m.end(); ++item_iter) {
        if(static_cast<const Set<K>&>(k).contains(item_iter->first)) {
            result.insert(*item_iter);
        }
    }
    return result;
}

template<class T> bool unique_elements(const std::vector<T>& lst) {
    Set<T> found;
    for(typename std::vector<T>::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter) {
        if(found.contains(*iter)) { return false; } else { found.insert(*iter); } }
    return true;
}

template<class T> Set<T> duplicate_elements(const std::vector<T>& lst) {
    Set<T> result; Set<T> found;
    for(typename std::vector<T>::const_iterator iter=lst.begin(); iter!=lst.end(); ++iter) {
        if(found.contains(*iter)) { result.insert(*iter); } else { found.insert(*iter); } }
    return result;
}


template<class T> inline Array<T> make_array(const T& t) { return Array<T>(1u,t); }
template<class T> inline Array<T> make_array(const std::vector<T>& vec) { return Array<T>(vec.begin(),vec.end()); }
template<class T> inline Array<T> make_array(const Array<T>& ary) { return Array<T>(ary); }

template<class T> inline List<T> make_list(const T& t) { return List<T>(1u,t); }
template<class T> inline List<T> make_list(const std::vector<T>& vec) { return List<T>(vec); }
template<class T> inline List<T> make_list(const Array<T>& ary) { return List<T>(ary); }
template<class T> inline List<T> make_list(const Set<T>& set) { return List<T>(set.begin(),set.end()); }


} // namespace Ariadne

#endif /* ARIADNE_CONTAINER_H */
