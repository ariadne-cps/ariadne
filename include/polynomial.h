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
#include <algorithm>
#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include "multi_index.h"

//#define ARIADNE_USE_ARRAY_POLYNOMIAL
#define ARIADNE_USE_MAP_POLYNOMIAL

namespace Ariadne {

#if defined ARIADNE_USE_ARRAY_POLYNOMIAL

template<class X>
class Reference {
  public:
    Reference<X>(X& x) : _x(&x) { }
    void operator=(Reference<X> r) { *_x=*r._x; }
    void operator=(const X& x) { *_x=x; }
    operator X& () { return *_x; }
    operator X const& () const { return *_x; }
  private:
    X* _x;
};

template<class X>
class ConstReference {
  public:
    ConstReference<X>(const X& x) : _x(&x) { }
    operator X const& () const { return *_x; }
  private:
    const X* _x;
};

template<class X> class PolynomialValue;
template<class X> class PolynomialReference;
template<class X> class PolynomialConstReference;

template<class X>
class PolynomialValue {
  public:
    PolynomialValue(const MultiIndex& a, const X& x) : first(a), second(x) { }
    MultiIndex first;
    X second;
};

template<class X>
inline bool operator<(const PolynomialValue<X>& dv1, const PolynomialValue<X>& dv2) {
    return dv1.first < dv2.first;
}

template<class X>
inline std::ostream& operator<<(std::ostream& os, const PolynomialValue<X>& dv) {
    return os << dv.first << ":" << dv.second;
}

template<class X>
class PolynomialReference {
  public:
    PolynomialReference<X>& operator=(const PolynomialReference<X>& dr) {
        first=dr.first; second=dr.second; return *this; }
    PolynomialReference<X>& operator=(const PolynomialValue<X>& dv) {
        first=dv.first; second=dv.second; return *this; }
    operator PolynomialValue<X>() const { 
        return PolynomialValue<X>(this->first,this->second); }
  public:
    MultiIndexReference first;
    Reference<X> second;
};

template<class X>
inline bool operator<(const PolynomialReference<X>& dr1, const PolynomialReference<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline bool operator<(const PolynomialReference<X>& dr1, const PolynomialValue<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline bool operator<(const PolynomialValue<X>& dr1, const PolynomialReference<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline std::ostream& operator<<(std::ostream& os, PolynomialReference<X> dr) {
    return os << dr.first << ":" << dr.second;
}


template<class X>
class PolynomialConstReference {
  public:
    MultiIndexConstReference first;
    ConstReference<X> second;
};

template<class X>
inline bool operator<(const PolynomialConstReference<X>& dr1, const PolynomialConstReference<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline std::ostream& operator<<(std::ostream& os, PolynomialConstReference<X> dr) {
    return os << dr.first << ":" << dr.second;
}

template<class X> class PolynomialIterator;
template<class X> class PolynomialConstIterator;

template<class X>
class PolynomialIterator {
    friend class PolynomialConstIterator<X>;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef X data_type;

    typedef std::random_access_iterator_tag iterator_category;
    typedef int difference_type;
    typedef PolynomialValue<X> value_type;
    typedef PolynomialReference<X>& reference;
    typedef PolynomialReference<X>* pointer;
  public:
    PolynomialIterator(size_type n, void* p) : _n(n), _p(reinterpret_cast<word_type*>(p)), 
        _x(reinterpret_cast<data_type*>(_p+MultiIndex::_word_size(_n))) { }
    PolynomialIterator(const PolynomialIterator<X>& i) : _n(i._n), _p(i._p), _x(i._x) { }
    PolynomialIterator<X>& operator=(const PolynomialIterator<X>& i) {
        _n=i._n; _p=i._p; _x=i._x; return *this; }
    int _offset() const { return MultiIndex::_word_size(_n); }
    int _increment() const { return MultiIndex::_element_size(_n); }
    PolynomialIterator<X>& operator++() { _p+=_increment(); _x=reinterpret_cast<data_type*>(_p+_offset()); return *this; }
    PolynomialIterator<X>& operator--() { _p-=_increment(); _x=reinterpret_cast<data_type*>(_p+_offset()); return *this; }
    PolynomialIterator<X> operator++(int) { PolynomialIterator tmp=*this; ++*this; return tmp; }
    PolynomialIterator<X> operator--(int) { PolynomialIterator tmp=*this; --*this; return tmp; }
    PolynomialIterator<X>& operator+=(const difference_type& m) {
        _p+=m*_increment(); _x=reinterpret_cast<data_type*>(_p+_offset()); return *this; }
    PolynomialIterator<X> operator+(const difference_type& n) const { PolynomialIterator r(*this); r+=n; return r; }
    PolynomialIterator<X> operator-(const difference_type& n) const { return (*this)+(-n); }
    difference_type operator-(const PolynomialIterator<X>& i) const { return (_p-i._p)/(_increment()); }
    PolynomialReference<X>& operator*() { return reinterpret_cast<PolynomialReference<X>&>(*this); }
    PolynomialReference<X>* operator->() { return reinterpret_cast<PolynomialReference<X>*>(this); }
    PolynomialConstReference<X> const* operator->() const { return reinterpret_cast<const PolynomialConstReference<X>*>(this); }
    bool operator==(const PolynomialIterator<X>& i) const { return _p==i._p; }
    bool operator!=(const PolynomialIterator<X>& i) const { return _p!=i._p; }
    bool operator<(const PolynomialIterator<X>& i) const { return _p<i._p; }
    const void* ptr() const { return _p; }
  private:
    size_type _n;
    word_type* _p;
    data_type* _x;
};

template<class X>
class PolynomialConstIterator {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::byte_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef X data_type;
  public:
    PolynomialConstIterator(size_type n, const void* p)
        : _n(n), _p(reinterpret_cast<const word_type*>(p)),
                 _x(reinterpret_cast<const data_type*>(_p+_offset())) { }
    PolynomialConstIterator(const PolynomialIterator<X>& piter)
        : _n(piter._n), _p(piter._p), _x(piter._x) { }
    int _offset() const {
        return MultiIndex::_word_size(_n); }
    int _increment() const {
        return MultiIndex::_element_size(_n); }
    PolynomialConstIterator<X>& operator++() {
        _p+=_increment(); _x=reinterpret_cast<const data_type*>(_p+_offset()); return *this; }
    PolynomialConstIterator<X>& operator--() {
        _p-=_increment(); _x=reinterpret_cast<const data_type*>(_p-_offset()); return *this; }
    const PolynomialConstReference<X>& operator*() const {
        return reinterpret_cast<const PolynomialConstReference<X>&>(*this); }
    const PolynomialConstReference<X>* operator->() const {
        return reinterpret_cast<const PolynomialConstReference<X>*>(this); }
    const void* ptr() const { return _p; }
    bool operator==(const PolynomialConstIterator<X>& i) const { return _p==i._p; }
    bool operator!=(const PolynomialConstIterator<X>& i) const { return _p!=i._p; }
  private:
  public:
    size_type _n;
    const word_type* _p;
    const data_type* _x;
};

template<class X>
inline std::ostream& operator<<(std::ostream& os, const PolynomialConstIterator<X>& pciter) {
    return os << "(n="<<pciter._n<<", p="<<(void*)pciter._p<<", a="<<pciter->first<<", x="<<*pciter._x<<")";
}

template<class X>
class Polynomial
{
  public:
    static X _zero;
  public:
    typedef X Real;
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::byte_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef X data_type;
    static const unsigned int sizeof_byte=sizeof(byte_type);
    static const unsigned int sizeof_word=sizeof(word_type);
    static const unsigned int sizeof_data=sizeof(data_type);
  public:
    typedef PolynomialValue<X> value_type;
    typedef PolynomialReference<X> reference;
    typedef PolynomialConstReference<X> const_reference;
    typedef PolynomialIterator<X> iterator;
    typedef PolynomialConstIterator<X> const_iterator;
  public:
    Polynomial() : _argument_size() { }
    Polynomial(unsigned int s) : _argument_size(s) { }
    Polynomial(unsigned int as, unsigned int deg, double c0, ...);
    template<class XX> Polynomial(const std::map<MultiIndex,XX>&);
    template<class XX> Polynomial(const Polynomial<XX>& p);
    template<class XX> Polynomial<X>& operator=(const Polynomial<XX>& p);
    static Polynomial<X> variable(unsigned int n, unsigned int i);
    static Vector< Polynomial<X> > variables(unsigned int s);
    bool operator==(const Polynomial<X>& p) const { return this->_v == p._v; }
    bool operator!=(const Polynomial<X>& p) const { return this->_v != p._v; }
    unsigned int argument_size() const { return _argument_size; }
    unsigned int size() const { return _v.size()/element_size(); }
    void reserve(unsigned int n) { _v.reserve(n*element_size()); }
    void erase(iterator iter) { iter->second=0.0; }
    void insert(const MultiIndex& a, const Real& x) {
        assert(a.size()==_argument_size); _insert(a,x); }
    void append(const MultiIndex& a, const Real& x) { _append(a,x); }
    void append(const MultiIndex& a1, const MultiIndex& a2, const Real& x) { _append(a1,a2,x); }
    Real& operator[](const MultiIndex& a) {
        iterator iter=this->find(a); 
        if(iter==this->end()) { return this->_insert(a,0.0)->second; }
        return iter->second; }
    const Real& operator[](const MultiIndex& a) const {
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
    void swap(Polynomial<X>& d) { std::swap(_argument_size,d._argument_size); _v.swap(d._v); }
  public:
    operator std::map<MultiIndex,Real> () const;
  public:
    const std::vector<word_type>& vector() const { return _v; }
    unsigned int vector_size() const { return _v.size(); }
    unsigned int word_size() const { return (_argument_size*+sizeof_word)/sizeof_word; }
    unsigned int element_size() const { return (_argument_size+sizeof_word)/sizeof_word+sizeof_data/sizeof_word; }
  public:
    const word_type* _begin() const { return _v.begin().operator->(); }
    const word_type* _end() const { return _v.end().operator->(); }
    word_type* _begin() { return _v.begin().operator->(); }
    word_type* _end() { return _v.end().operator->(); }
    iterator _insert(const MultiIndex& a, const Real& x) {
        this->_append(a,x);
        iterator iter=this->end(); --iter; iterator prev=iter;
        while(iter!=this->begin()) {  
            --prev; if(prev->first<a) { break; } else { *iter=*prev; --iter; } }
        iter->first=a; iter->second=x; return iter; }
    void _append(const MultiIndex& a, const Real& x) {
        _v.resize(_v.size()+element_size());
        word_type* vp=_end()-element_size(); const word_type* ap=a.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end())-1; *xp=x; }
    void _append(const MultiIndex&  a1, const MultiIndex&  a2, const Real& x) {
        _v.resize(_v.size()+element_size());
        word_type* vp=_end()-element_size(); const word_type* ap1=a1.word_begin(); const word_type* ap2=a2.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap1[j]+ap2[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end())-1; *xp=x; }
  private:
    unsigned int _argument_size;
    std::vector<word_type> _v;
};

template<class X> template<class XX>
inline Polynomial<X>::Polynomial(const std::map<MultiIndex,XX>& m)
{
    ARIADNE_ASSERT(!m.empty());
    this->_argument_size=m.begin()->first.size();
    for(typename std::map<MultiIndex,XX>::const_iterator iter=m.begin(); iter!=m.end(); ++iter) {
        this->append(iter->first,X(iter->second));
    }
}

template<class X> template<class XX> inline
Polynomial<X>::Polynomial(const Polynomial<XX>& p)
    : _argument_size(p.argument_size())
{
    for(typename Polynomial<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->first,X(iter->second)); }
}

template<class X>
inline Polynomial<X>::operator std::map<MultiIndex,X> () const
{
    std::map<MultiIndex,X> r;
    for(typename Polynomial<X>::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        r.insert(std::make_pair(iter->first,iter->second));
    }
    return r;
}

template<class X> inline Polynomial<X>::Polynomial(uint as, uint deg, double c0, ...)
{
    MultiIndex a(as); double x;
    va_list args; va_start(args,c0);
    while(a.degree()<=deg) {
        if(a.degree()==0) { x=c0; }
        else { x=va_arg(args,double); }
        if(x!=0) { this->insert(a,x); }
        ++a;
    }
    va_end(args);
    this->sort();
}



#elif defined ARIADNE_USE_MAP_POLYNOMIAL

template<class X>
class Polynomial
    : public std::map<MultiIndex,X>
{
  public:
    typedef MultiIndex::size_type size_type;
    typedef typename std::map<MultiIndex,X>::const_iterator const_iterator;
  private:
    static X _zero;
  public:
    Polynomial() : _argument_size(0) { }
    Polynomial(size_type as) : _argument_size(as) { }
    Polynomial(unsigned int as, unsigned int deg, double c0, ...);
    Polynomial(const std::map<MultiIndex,X>& m) : std::map<MultiIndex,X>(m) {
        assert(!m.empty()); _argument_size=m.begin()->first.size(); }
    template<class XX> Polynomial(const Polynomial<XX>& p);
    static Polynomial<X> variable(unsigned int n, unsigned int i);
    static Vector< Polynomial<X> > variables(unsigned int s);

    X& operator[](const MultiIndex& a) {
        return this->std::map<MultiIndex,X>::operator[](a); }
    const X& operator[](const MultiIndex& a) const {
        const_iterator iter=this->find(a); 
        if(iter==this->end()) { return _zero; }
        else { return iter->second; } }
    size_type argument_size() const { return this->_argument_size; }
    void insert(const MultiIndex& a, const X& x) {
        this->std::map<MultiIndex,X>::insert(std::make_pair(a,x)); }
    void insert(const MultiIndex& a1, const MultiIndex& a2, const X& x) { this->insert(a1+a2,x); }
    void append(const MultiIndex& a, const X& x) { this->insert(a,x); }
    void append(const MultiIndex& a1, const MultiIndex& a2, const X& x) { this->insert(a1+a2,x); }
    void sort() { }
  private:
    size_type _argument_size;
};

template<class X> inline Polynomial<X>::Polynomial(uint as, uint deg, double c0, ...)
    : _argument_size(as)
{
    MultiIndex a(as); double x;
    va_list args; va_start(args,c0);
    while(a.degree()<=deg) {
        if(a.degree()==0) { x=c0; }
        else { x=va_arg(args,double); }
        if(x!=0) { this->insert(a,x); }
        ++a;
    }
    va_end(args);
    this->sort();
}

template<class X> template<class XX> inline
Polynomial<X>::Polynomial(const Polynomial<XX>& p)
    : _argument_size(p.argument_size())
{
    for(typename Polynomial<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->first,X(iter->second)); }
}

#else

#error "No Polynomial class activated."

#endif




template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p) {
    return p; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p) {
    Polynomial<X> r(p.argument_size());  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { r[iter->first]-=iter->second; } return r; }
template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p, const Float& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]+=c; return r; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p, const Float& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]-=c; return r; }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p, const Float& c) {
    if(c==0) { return Polynomial<X>(p.argument_size()); }
    Polynomial<X> r(p); typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->second*=c; } return r; }
template<class X> inline Polynomial<X> operator/(const Polynomial<X>& p, const Float& c) {
    Polynomial<X> r(p); typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->second/=c; } return r; }
template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->first]+=iter->second; } return r; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->first]-=iter->second; } return r; }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    Polynomial<X> r(p1.argument_size()); 
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(Iter iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->first+iter2->first;
            r[a]+=iter1->second*iter2->second; } } return r; }
template<class X> inline Polynomial<X> operator+(const X& c, const Polynomial<X>& p) {
    return p+c; }
template<class X> inline Polynomial<X> operator-(const X& c, const Polynomial<X>& p) {
    return (-p)+c; }
template<class X> inline Polynomial<X> operator*(const X& c, const Polynomial<X>& p) {
    return p*c; }

template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Polynomial<X>& q) {
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->first]+=iter->second; } return p; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Polynomial<X>& q) {
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->first]-=iter->second; } return p; }
template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Float& c) {
    p[MultiIndex(p.argument_size())]+=c; return p; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Float& c) {
    p[MultiIndex(p.argument_size())]-=c; return p; }
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const Float& c) {
    typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->second*=c; } return p; }
template<class X> inline Polynomial<X>& operator/=(Polynomial<X>& p, const Float& c) {
    typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->second/=c; } return p; }

inline Polynomial<Float> midpoint(const Polynomial<Interval>& p) {
    Polynomial<Float> r(p.argument_size());
    for(Polynomial<Interval>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        r.append(iter->first,midpoint(iter->second)); }
    return r;
}

inline Vector< Polynomial<Float> > midpoint(const Vector< Polynomial<Interval> >& p) {
    Vector< Polynomial<Float> > r(p.size());
    for(uint i=0; i!=p.size(); ++i) {
        r[i]=midpoint(p[i]); }
    return r;
}

template<class X, class Y>
Y evaluate(const Polynomial<X>& p, const Vector<Y>& x)
{
    //std::cerr<<ARIADNE_PRETTY_FUNCTION<<std::endl;
    ARIADNE_ASSERT(p.argument_size()==x.size());

    Y zero = x[0]; zero*=0;
    Y one = zero; one+=1;

    Y r=zero;

    for(typename Polynomial<X>::const_iterator iter=p.begin();
        iter!=p.end(); ++iter)
    {
        const MultiIndex& j=iter->first;
        const X& c=iter->second;
        Y t=one;
        for(uint k=0; k!=x.size(); ++k) {
            for(uint l=0; l!=j[k]; ++l) {
                t=t*x[k];
            }
        }
        t*=c;
        r+=t;
    }

    return r;
}


template<class X, class Y>
Vector<Y> evaluate(const Vector< Polynomial<X> >& p, const Vector<Y>& x)
{
    //std::cerr<<ARIADNE_PRETTY_FUNCTION<<std::endl;
    ARIADNE_ASSERT(p.size()>0 && p[0].argument_size()==x.size());

    Y zero = x[0]; zero*=0.0;

    Vector<Y> r(p.size(),zero);

    for(uint i=0; i!=p.size(); ++i) {
        r[i]=evaluate(p[i],x);
    }
    return r;
}


template<class X> inline
Polynomial<X> compose(const Polynomial<X>& p, const Vector< Polynomial<X> >& q)
{
    return evaluate(p,q);
}

template<class X> inline
Vector< Polynomial<X> > compose(const Vector< Polynomial<X> >& p, const Vector< Polynomial<X> >& q)
{
    return evaluate(p,q);
}


template<class X> inline
Polynomial<X> embed(const Polynomial<X>& x, uint new_size, uint start)
{
    ARIADNE_ASSERT(x.argument_size()+start<=new_size);
    Polynomial<X> r(new_size);
    uint old_size=x.argument_size();
    MultiIndex old_index(old_size);
    MultiIndex new_index(new_size);
    for(typename Polynomial<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        old_index=iter->first;
        for(uint j=0; j!=old_size; ++j) {
            uint aj=old_index[j];
            new_index[j+start]=aj;
        }
        r.append(new_index,iter->second);
    }
    return r;
}

template<class X>
std::ostream& operator<<(std::ostream& os, const Polynomial<X>& p) {
    os << "{";
    for(typename Polynomial<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? "" : ",");
        for(unsigned int i=0; i!=iter->first.size(); ++i) {
            os << (i==0?" ":",") << int(iter->first[i]); }
        os << ":" << iter->second; }
    return os << " }";
}


template<class X> Polynomial<X> Polynomial<X>::variable(uint n, uint i) {
    Polynomial<X> p(n); p[MultiIndex::zero(n)]=0.0; p[MultiIndex::unit(n,i)]=X(1);
    return p;
}

template<class X> Vector< Polynomial<X> > Polynomial<X>::variables(uint n) {
    Vector< Polynomial<X> > pv(n,Polynomial<X>(n)); MultiIndex a(n);
    for(uint i=0; i!=n; ++i) { a[i]=1; pv[i][a]=1.0; a[i]=0; }
    return pv; 
}

template<class X> Vector< Polynomial<X> > operator*(const Vector<Float> e, const Polynomial<X>& p) {
    Vector< Polynomial<X> > r(e.size(),Polynomial<X>(p.argument_size()));
    for(uint i=0; i!=r.size(); ++i) { r[i]=X(e[i])*p; }
    return r;
}

}

#endif /* ARIADNE_POLYNOMIAL_H */
