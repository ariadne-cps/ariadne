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
#include <utility>
#include <vector>
#include <set>
#include <map>
#include "array.h"

namespace Ariadne {


typedef std::string String;

// Functionality for building arrays or std::vector from comma-separated lists
class ArrayBuilderTag {};
static ArrayBuilderTag build_array;

template<class T>
inline std::vector<T> operator,(ArrayBuilderTag, const T& t) {
    std::vector<T> v; v.push_back(t); return v; }

template<class T>
inline std::vector<T> operator,(const T& t1, const T& t2) {
    std::vector<T> v; v.push_back(t1); v.push_back(t2); return v; }
template<class T>
inline std::vector<T> operator,(const std::vector<T>& v, const T& t) {
    std::vector<T> r(v); r.push_back(t); return r; }

inline std::vector<std::string> operator,(const std::string& s1, const char* s2) {
    std::vector<std::string> v; v.push_back(std::string(s1)); v.push_back(std::string(s2)); return v; }
inline std::vector<std::string> operator,(const std::vector<std::string>& v1, const char* s2) {
    std::vector<std::string> r(v1); r.push_back(std::string(s2)); return r; }
inline std::vector<std::string> operator,(ArrayBuilderTag, const char* s) {
    std::vector<std::string> v; v.push_back(std::string(s)); return v; }


template<class T> inline bool contains(const std::set<T>& s, const T& t) {
    return s.find(t)!=s.end(); }
template<class T> inline bool subset(const std::set<T>& s1, const std::set<T>& s2) {
    for(typename std::set<T>::iterator iter=s1.begin(); iter!=s1.end(); ++iter) {
        if(!contains(s2,*iter)) { return false; } } return true; }
template<class T> inline bool disjoint(const std::set<T>& s1, const std::set<T>& s2) {
    for(typename std::set<T>::iterator iter=s1.begin(); iter!=s1.end(); ++iter) {
        if(contains(s2,*iter)) { return false; } } return true; }
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
template<class T> inline std::set<T> difference(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> r(s1); remove(r,s2); return r; }
template<class K,class V,class C> inline bool has_key(const std::map<K,V,C>& m, const K& k) {
    return m.find(k)!=m.end(); }


// Temporary home for useful auxiliary classes

template<class T1, class T2> class Pair
    : public std::pair<T1,T2>
{
  public:
    Pair() : std::pair<T1,T2>() { }
    template<class X> Pair(const X& x) : std::pair<T1,T2>(x) { }
    template<class X1, class X2> Pair(const X1& x1, const X2& x2) : std::pair<T1,T2>(x1,x2) { }
};


template<class T> class Array
    : public Ariadne::array<T>
{
  public:
    Array() : Ariadne::array<T>() { }
    Array(const array<T>& l) : Ariadne::array<T>(l) { }
    template<class I> Array(const I& b, const I& e) : Ariadne::array<T>(b,e) { }
};

template<class T> class List
    : public std::vector<T>
{
  public:
    List() : std::vector<T>() { }
    List(const std::vector<T>& l) : std::vector<T>(l) { }
    template<class I> List(const I& b, const I& e) : std::vector<T>(b,e) { }
    void append(const T& t) { this->push_back(t); }
};


template<class T> class Set
    : public std::set<T>
{
  public:
    Set() : std::set<T>() { }
    Set(const std::set<T>& s) : std::set<T>(s) { }

    bool contains(const T& t) {
        return this->find(t)==this->end(); }
    bool subset(const std::set<T>& s) {
        for(typename std::set<T>::iterator iter=s.begin(); iter!=s.end(); ++iter) {
            if(!this->contains(*iter)) { return false; } } return true; }
    bool disjoint(const std::set<T>& s) {
        for(typename std::set<T>::iterator iter=s.begin(); iter!=s.end(); ++iter) {
            if(this->contains(*iter)) { return false; } } return true; }
    Set<T>& adjoin(const std::set<T>& s) {
        for(typename std::set<T>::iterator iter=s.begin(); iter!=s.end(); ++iter) { this->insert(*iter); } return *this; }
    Set<T>& remove(const std::set<T>& s) {
        typename std::set<T>::iterator iter=this->begin();
        while(iter!=this->end()) { if(s.find(*iter)!=s.end()) { this->erase(iter++); } else { ++iter; } }
        return *this; }
    Set<T>& restrict(const std::set<T>& s) {
        typename std::set<T>::iterator iter=this->begin();
        while(iter!=this->end()) { if(s.find(*iter)!=s.end()) { this->erase(iter++); } else { ++iter; } }
        return *this; }
};

template<class T> inline Set<T> join(const Set<T>& s1, const Set<T>& s2) {
    Set<T> r(s1); adjoin(r,s2); return r; }

template<class T> inline std::ostream& operator<<(std::ostream& os, const Set<T>& s) {
    return os<<static_cast<const std::set<T>&>(s); }

template<class K, class V> class Map
    : public std::map<K,V>
{
  public:
    V& operator[](const K& k) {
        return this->std::map<K,V>::operator[](k); }
    const V& operator[](const K& k) const { typename std::map<K,V>::const_iterator p=this->find(k);
        assert(p!=this->end()); return p->second; }
    bool has_key(const K& k) const {
        return this->find(k)!=this->end(); }
    void insert(const K& k, const V& v) {
        this->std::map<K,V>::insert(std::make_pair(k,v)); }
    using std::map<K,V>::insert;
};


} // namespace Ariadne

#endif /* ARIADNE_CONTAINER_H */
