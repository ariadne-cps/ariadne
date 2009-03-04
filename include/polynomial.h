/***************************************************************************
 *            polynomial.h
 *
 *  Copyright 2008  Pieter Collins
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

/*! \file polynomial.h
 *  \brief Base class for polynomial rings.
 */
#ifndef ARIADNE_POLYNOMIAL_H
#define ARIADNE_POLYNOMIAL_H

#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include "multi_index.h"

namespace Ariadne {

class PolynomialValue;
class PolynomialReference;
class PolynomialConstReference;

class PolynomialValue {
  public:
    PolynomialValue(const MultiIndex& a, const Float& x) : first(a), second(x) { }
    PolynomialValue(const PolynomialReference& dr);
    MultiIndex first;
    Float second;
};

inline bool operator<(const PolynomialValue& dv1, const PolynomialValue& dv2) {
    return dv1.first < dv2.first;
}

inline std::ostream& operator<<(std::ostream& os, const PolynomialValue& dv) {
    return os << dv.first << ":" << dv.second;
}

class PolynomialReference {
  public:
    PolynomialReference(unsigned int n, uchar* p) : first(n,p), second(*(double*)(p+(4*((n+3)/4)))) { }
    PolynomialReference& operator=(const PolynomialReference& dr) {
        first=dr.first; second=dr.second; return *this; }
    PolynomialReference& operator=(const PolynomialValue& dv) {
        first=dv.first; second=dv.second; return *this; }
  public:
    MultiIndexReference first;
    double& second;
};

inline bool operator<(const PolynomialReference& dr1, const PolynomialReference& dr2) {
    return dr1.first < dr2.first;
}

inline PolynomialValue::PolynomialValue(const PolynomialReference& dr)
    : first(dr.first), second(dr.second)
{ }


class PolynomialConstReference {
  public:
    PolynomialConstReference(unsigned int n, const uchar* p) : first(n,p), second(*(const double*)(p+(4*((n+3)/4)))) { }
    PolynomialConstReference(PolynomialReference dr) : first(dr.first), second(dr.second) { }
    PolynomialConstReference(const PolynomialConstReference& dr) : first(dr.first), second(dr.second) { }
    PolynomialConstReference& operator=(const PolynomialConstReference& dr) {
        first=dr.first; const double* second_ptr=&second; second_ptr=(const double*)(&dr.second); return *this; }
    MultiIndexConstReference first;
    const double& second;
};

inline std::ostream& operator<<(std::ostream& os, PolynomialConstReference dr) {
    return os << dr.first << ":" << dr.second;
}

class PolynomialIterator;
class PolynomialConstIterator;

class PolynomialIterator {
    friend class PolynomialConstIterator;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type byte_type;
    typedef MultiIndex::word_type word_type;

    typedef std::random_access_iterator_tag iterator_category;
    typedef int difference_type;
    typedef PolynomialValue value_type;
    typedef PolynomialReference& reference;
    typedef PolynomialReference* pointer;
  public:
    PolynomialIterator(size_type n, void* p) : _n(n), _p(reinterpret_cast<word_type*>(p)), _x((double*)(_p+_offset())) { }
    PolynomialIterator(const PolynomialIterator& i) : _n(i._n), _p(i._p), _x(i._x) { }
    PolynomialIterator& operator=(const PolynomialIterator& i) {
        _n=i._n; _p=i._p; _x=i._x; return *this; }
    int _offset() const { return (_n*sizeof(byte_type))/sizeof(word_type)+1; }
    int _increment() const { return (_n*sizeof(byte_type)+sizeof(double))/sizeof(word_type)+1; }
    PolynomialIterator& operator++() { _p+=_increment(); _x=(double*)(_p+_offset()); return *this; }
    PolynomialIterator& operator--() { _p-=_increment(); _x=(double*)(_p+_offset()); return *this; }
    PolynomialIterator operator++(int) { PolynomialIterator tmp=*this; ++*this; return tmp; }
    PolynomialIterator operator--(int) { PolynomialIterator tmp=*this; --*this; return tmp; }
    PolynomialIterator& operator+=(const difference_type& m) {
        _p+=m*_increment(); _x=(double*)(_p+_offset()); return *this; }
    PolynomialIterator operator+(const difference_type& n) const { PolynomialIterator r(*this); r+=n; return r; }
    PolynomialIterator operator-(const difference_type& n) const { return (*this)+(-n); }
    difference_type operator-(const PolynomialIterator& i) const { return (_p-i._p)/(_increment()); }
    PolynomialReference& operator*() { return reinterpret_cast<PolynomialReference&>(*this); }
    PolynomialReference* operator->() { return reinterpret_cast<PolynomialReference*>(this); }
    PolynomialConstReference const* operator->() const { return reinterpret_cast<const PolynomialConstReference*>(this); }
    bool operator==(const PolynomialIterator& i) const { return _p==i._p; }
    bool operator!=(const PolynomialIterator& i) const { return _p!=i._p; }
    bool operator<(const PolynomialIterator& i) const { return _p<i._p; }
    const void* ptr() const { return _p; }
  private:
    size_type _n;
    word_type* _p;
    double* _x;
};

class PolynomialConstIterator {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::byte_type byte_type;
    typedef MultiIndex::word_type word_type;
  public:
    PolynomialConstIterator(size_type n, const void* p)
        : _n(n), _p(reinterpret_cast<const word_type*>(p)), _x((const double*)(_p+_offset())) { }
    PolynomialConstIterator(const PolynomialIterator& piter) : _n(piter._n), _p(piter._p), _x(piter._x) { }
    int _offset() const { return (_n*sizeof(byte_type))/sizeof(word_type)+1; }
    int _increment() const { return (_n*sizeof(byte_type)+sizeof(double))/sizeof(word_type)+1; }
    PolynomialConstIterator& operator++() { _p+=_increment(); _x=(const double*)(_p+_offset()); return *this; }
    PolynomialConstIterator& operator--() { _p-=_increment(); _x=(const double*)(_p-_offset()); return *this; }
    const PolynomialConstReference& operator*() const {
        return reinterpret_cast<const PolynomialConstReference&>(*this); }
    const PolynomialConstReference* operator->() const {
        return reinterpret_cast<const PolynomialConstReference*>(this); }
    const void* ptr() const { return _p; }
    bool operator==(const PolynomialConstIterator& i) const { return _p==i._p; }
    bool operator!=(const PolynomialConstIterator& i) const { return _p!=i._p; }
  private:
  public:
    size_type _n;
    const word_type* _p;
    const double* _x;
};

inline std::ostream& operator<<(std::ostream& os, const PolynomialConstIterator& pciter) {
    return os << "(n="<<pciter._n<<", p="<<(void*)pciter._p<<", a="<<pciter->first<<", x="<<*pciter._x<<")";
}

class Polynomial
{
  public:
    static Float _zero;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::byte_type byte_type;
    typedef MultiIndex::word_type word_type;
    static const unsigned int sizeof_byte=sizeof(byte_type);
    static const unsigned int sizeof_word=sizeof(word_type);
    static const unsigned int sizeof_double=sizeof(double);
  public:
    typedef PolynomialValue value_type;
    typedef PolynomialReference reference;
    typedef PolynomialConstReference const_reference;
    typedef PolynomialIterator iterator;
    typedef PolynomialConstIterator const_iterator;
  public:
    Polynomial(unsigned int s) : _argument_size(s) { }
    Polynomial(const std::map<MultiIndex,Float>&);
    bool operator==(const Polynomial& p) const { return this->_v == p._v; }
    unsigned int argument_size() const { return _argument_size; }
    unsigned int size() const { return _v.size()/element_size(); }
    void reserve(unsigned int n) { _v.reserve(n*element_size()); }
    void erase(const_iterator iter) { const_cast<double&>(iter->second)=0.0; }
    void insert(const MultiIndex& a, Float x) {
        assert(a.size()==_argument_size); _insert(a,x); }
    void insert(const MultiIndex& a1, const MultiIndex& a2, Float x) {
        assert(a1.size()==_argument_size); assert(a2.size()==_argument_size); _insert(a1,a2,x); }
    void append(const MultiIndex& a, Float x) { _insert(a,x); }
    void append(const MultiIndex& a1, const MultiIndex& a2, Float x) { _insert(a1,a2,x); }
    reference operator[](unsigned int i) { return reference(_argument_size,(byte_type*)(&_v[0]+element_size()*i)); }
    const_reference operator[](unsigned int i) const {
        return const_reference(_argument_size,(byte_type*)(&_v[0]+element_size()*i)); }
    Float& operator[](const MultiIndex& a) { 
        iterator iter=this->find(a); 
        if(iter==this->end()) { this->_insert(a,0.0); iter=this->find(a); }
        return iter->second; }
    const Float& operator[](const MultiIndex& a) const { 
        const_iterator iter=this->find(a); 
        if(iter==this->end()) { return _zero; }
        else { return iter->second; } }
    iterator begin() { return iterator(_argument_size,(byte_type*)_begin()); }
    iterator end() { return iterator(_argument_size,(byte_type*)_end()); }
    const_iterator begin() const { return const_iterator(_argument_size,(byte_type*)_begin()); }
    const_iterator end() const { return const_iterator(_argument_size,(byte_type*)_end()); }
    iterator find(const MultiIndex& a) {
        iterator iter=this->begin(); 
        while(iter!=this->end() && iter->first!=a) { ++iter; } 
        return iter; }
    const_iterator find(const MultiIndex& a) const {
        const_iterator iter=this->begin(); 
        while(iter!=this->end() && iter->first!=a) { ++iter; } 
        return iter; };
    void sort() { std::sort(this->begin(),this->end()); }
    void unique_sort();
    void clear() { _v.clear(); }
    void swap(Polynomial& d) { std::swap(_argument_size,d._argument_size); _v.swap(d._v); }
  public:
    operator std::map<MultiIndex,Float> () const;
  public:
    const std::vector<word_type>& vector() const { return _v; }
    unsigned int vector_size() const { return _v.size(); }
    unsigned int word_size() const { return (_argument_size*+sizeof_word)/sizeof_word; }
    unsigned int element_size() const { return (_argument_size+sizeof_word)/sizeof_word+sizeof_double/sizeof_word; }
  public:
    const word_type* _begin() const { return _v.begin().operator->(); }
    const word_type* _end() const { return _v.end().operator->(); }
    word_type* _begin() { return _v.begin().operator->(); }
    word_type* _end() { return _v.end().operator->(); }
    void _insert(const MultiIndex& a, Float x) {
        _v.resize(_v.size()+element_size());
        word_type* vp=_end()-element_size(); const word_type* ap=a.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap[j]; }
        double* xp=reinterpret_cast<double*>(this->_end())-1; *xp=x; }
    void _insert(const MultiIndex&  a1, const MultiIndex&  a2, Float x) {
        _v.resize(_v.size()+element_size());
        word_type* vp=_end()-element_size(); const word_type* ap1=a1.word_begin(); const word_type* ap2=a2.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap1[j]+ap2[j]; }
        double* xp=reinterpret_cast<double*>(this->_end())-1; *xp=x; }
  private:
    unsigned int _argument_size;
    std::vector<word_type> _v;
};

inline Polynomial::Polynomial(const std::map<MultiIndex,Float>& m)
{
    ARIADNE_ASSERT(!m.empty());
    this->_argument_size=m.begin()->first.size();
    for(std::map<MultiIndex,Float>::const_iterator iter=m.begin(); iter!=m.end(); ++iter) {
        this->append(iter->first,iter->second);
    }
}

inline Polynomial::operator std::map<MultiIndex,Float> () const
{
    std::map<MultiIndex,Float> r;
    for(Polynomial::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        r.insert(std::make_pair(iter->first,iter->second));
    }
    return r;
}



inline std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
    os << "{";
    for(Polynomial::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? ' ' : ',') << iter->first << ":" << iter->second; }
    return os << " }";
}


class MapPolynomial : public std::map<MultiIndex,Float>
{
    typedef MultiIndex::size_type size_type;
  private:
    static double _zero;
  public:
    MapPolynomial(size_type as) : _argument_size(as) { }
    MapPolynomial(const std::map<MultiIndex,Float>& m) : std::map<MultiIndex,Float>(m) { 
        assert(!m.empty()); _argument_size=m.begin()->first.size(); }
    Float& operator[](const MultiIndex& a) { 
        return this->std::map<MultiIndex,Float>::operator[](a); }
    const Float& operator[](const MultiIndex& a) const { 
        const_iterator iter=this->find(a); 
        if(iter==this->end()) { return _zero; }
        else { return iter->second; } }
    size_type argument_size() const { return this->_argument_size; }
    void insert(const MultiIndex& a, Float x) { this->std::map<MultiIndex,Float>::insert(std::make_pair(a,x)); }
    void insert(const MultiIndex& a1, const MultiIndex& a2, Float x) { this->insert(a1+a2,x); }
    void append(const MultiIndex& a, Float x) { this->insert(a,x); }
    void append(const MultiIndex& a1, const MultiIndex& a2, Float x) { this->insert(a1+a2,x); }
    void sort() { }
  private:
    size_type _argument_size;
};


}

#endif /* ARIADNE_POLYNOMIAL_H */
