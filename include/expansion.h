/***************************************************************************
 *            expansion.h
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

/*! \file expansion.h
 *  \brief Base class for power series expansions.
 */

// Use vector-based expansion and casting of iterators to get a reference

#ifndef ARIADNE_EXPANSION_H
#define ARIADNE_EXPANSION_H

#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include "vector.h"
#include "multi_index.h"



namespace Ariadne {

template<class X> class Expansion;

template<class X> Expansion<X> embed(unsigned int, const Expansion<X>&, unsigned int);



template<class FwdIter, class Op>
FwdIter unique_key(FwdIter first, FwdIter last, Op op) {
    FwdIter curr=first;
    FwdIter next=curr;
    while(next!=last) {
        if(curr!=next) { *curr=*next; }
        ++next;
        while(next!=last && curr->key()==next->key()) {
            if(curr->key()==next->key()) {
                curr->data()=op(curr->data(),next->data());
                ++next;
            }
        }
        // Removes zero entries; the code below is preferred to the case "curr->data()!=0" for tribool results
        if(curr->data()==0) { }
        else { ++curr; }
    }
    return curr;
}



/* The following code has debugging statements
template<class X>
struct ExpansionValue {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::word_type word_type;

    ~ExpansionValue() {
        //std::cerr<<"destroy "<< *this<<"... "<<std::flush;
        delete[] _p;
        //std::cerr<<" done"<<std::endl;
    }
    ExpansionValue(const MultiIndex& a, const X& x)
        : _n(a.size()), _nw(a.word_size()), _p() {
        //std::cerr<<"create... "<<std::flush;
        _p=new word_type[_nw+_ds]; key()=a; data()=x;
        //std::cerr<<*this<<" done"<<std::endl;
    }
    ExpansionValue(const ExpansionValue<X>& v)
        : _n(v._n), _nw(v._nw), _p() {
        //std::cerr<<"copy construct from "<<v<<"... "<<std::flush;
        _p=new word_type[v._nw+_ds];
        //std::cerr<<" assigning... "<<std::flush;
        _assign(v._p);
        //std::cerr<<*this<<" done"<<std::endl;
    }
    ExpansionValue<X>& operator=(const ExpansionValue<X>& v) {
        //std::cerr<<"operator= "<<*this<<" "<<v<<"... "<<std::flush;
        if(this!=&v) { _resize(v._n,v._nw); _assign(v._p); }
        //std::cerr<<" done"<<std::endl;
        return *this;
    }
    const MultiIndex& key() const { return *reinterpret_cast<const MultiIndex*>(this); }
    const X& data() const { return *reinterpret_cast<const X*>(_p+_nw); }
    MultiIndex& key() { return *reinterpret_cast<MultiIndex*>(this); }
    X& data() { return *reinterpret_cast<X*>(_p+_nw); }
  private:
    void _resize(size_type n, size_type nw) {
        //std::cerr<<"resize"<<*this<<" "<<n<<" "<<nw<<"..."<<std::flush;
        if(_nw!=nw) {
            //std::cerr<<" reallocating..."<<std::flush;
            delete[] _p; _nw=nw; _p=new word_type[_nw+_ds];
        }
        _n=n;
        //std::cerr<<" done"<<std::endl;
    }
    void _assign(const word_type* p) {
        //std::cerr<<"assign"<<*this<<" "<<(void*)_p<<"..."<<std::flush;
        for(size_type j=0; j!=_nw+_ds; ++j) { _p[j]=p[j]; }
        //std::cerr<<" done"<<std::endl;
    }
    //void _assign(const word_type* p) { std::memcpy(_p,p,(_nw+_ds)*sizeof(word_type)); }
  public:
    size_type _n; size_type _nw; word_type* _p;
    static const size_type _ds=sizeof(X)/sizeof(word_type);
};
*/

template<class X>
struct ExpansionValue {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::word_type word_type;

    ~ExpansionValue() { delete[] _p; }
    ExpansionValue(const MultiIndex& a, const X& x)
        : _n(a.size()), _nw(a.word_size()), _p() { _p=new word_type[_nw+_ds]; key()=a; data()=x;  }
    ExpansionValue(const ExpansionValue<X>& v)
        : _n(v._n), _nw(v._nw), _p() { _p=new word_type[v._nw+_ds]; _assign(v._p); }
    ExpansionValue<X>& operator=(const ExpansionValue<X>& v) {
        if(this!=&v) { _resize(v._n,v._nw); _assign(v._p); } return *this; }
    const MultiIndex& key() const { return *reinterpret_cast<const MultiIndex*>(this); }
    const X& data() const { return *reinterpret_cast<const X*>(_p+_nw); }
    MultiIndex& key() { return *reinterpret_cast<MultiIndex*>(this); }
    X& data() { return *reinterpret_cast<X*>(_p+_nw); }
  private:
    void _resize(size_type n, size_type nw) {
        if(_nw!=nw) { delete[] _p; _nw=nw; _p=new word_type[_nw+_ds]; }_n=n; }
    void _assign(const word_type* p) {
        for(size_type j=0; j!=_nw+_ds; ++j) { _p[j]=p[j]; } }
  public:
    size_type _n; size_type _nw; word_type* _p;
    static const size_type _ds=sizeof(X)/sizeof(word_type);
};


struct DataLess {
    template<class X> bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) { return v1.data()<v2.data(); }
    template<class X> bool operator()(const ExpansionValue<X>& v1, const X& d2) { return v1.data()<d2; }
};

struct DataIsZero {
    template<class X> bool operator()(const ExpansionValue<X>& v) { return v.data()==0; }
};

struct KeyEquals {
    template<class X> bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) { return v1.key()==v2.key(); }
    template<class X> bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) { return v1.key()==k2; }
    bool operator()(const MultiIndex& k1, const MultiIndex& k2) { return k2==k2; }
};

struct GradedKeyLess {
    template<class X> bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) {
        return graded_less(v1.key(),v2.key()); }
    template<class X> bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) {
        return graded_less(v1.key(),k2); }
    bool operator()(const MultiIndex& k1, const MultiIndex& k2) {
        return graded_less(k1,k2); }
};

struct LexicographicKeyLess {
    template<class X> bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) {
        return lexicographic_less(v1.key(),v2.key()); }
    template<class X> bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) {
        return lexicographic_less(v1.key(),k2);; }
    bool operator()(const MultiIndex& k1, const MultiIndex& k2) {
        return lexicographic_less(k1,k2); }
};

struct ReverseLexicographicKeyLess {
    template<class X> bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) {
        return reverse_lexicographic_less(v1.key(),v2.key()); }
    template<class X> bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) {
        return reverse_lexicographic_less(v1.key(),k2);; }
    bool operator()(const MultiIndex& k1, const MultiIndex& k2) {
        return reverse_lexicographic_less(k1,k2); }
};



template<class X>
inline std::ostream& operator<<(std::ostream& os, const ExpansionValue<X>& dv) {
    return os << "["<<dv._n<<","<<dv._nw<<","<<(void*)dv._p<<"]"<<dv.key()<<":"<<dv.data();
    //return os << dv.key() << ":" << dv.data();
}


template<class X, class Ref, class Ptr> class ExpansionIterator;
template<class X, class Ref, class Ptr> std::ostream& operator<<(std::ostream&, const ExpansionIterator<X,Ref,Ptr>&);

// Iterator for Expansion<X>
// Has the same data layout as an ExpansionReference, so can be easily cast
// into this class.
// It would be possible to cache the increment in words, but is seems that
// this does not increase performance
template<class X, class Ref, class Ptr>
class ExpansionIterator
     : public boost::iterator_facade<ExpansionIterator<X,Ref,Ptr>,
                                     ExpansionValue<X>,
                                     boost::random_access_traversal_tag,
                                     Ref >

{
    template<class X2, class Ref2, class Ptr2> friend class ExpansionIterator;
    friend std::ostream& operator<< <>(std::ostream&,const ExpansionIterator<X,Ref,Ptr>&);
    typedef ExpansionIterator<X,Ref,Ptr> Iter;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::index_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef X data_type;

    typedef int difference_type;
  public:
    ExpansionIterator()
        : _n(), _nw(), _p() { }
    ExpansionIterator(size_type n, size_type nw, word_type* p)
        : _n(n), _nw(nw), _p(p) { }
    template<class Ref2,class Ptr2> ExpansionIterator(const ExpansionIterator<X,Ref2,Ptr2>& i)
        : _n(i._n), _nw(i._nw), _p(i._p) { }
    template<class Ref2,class Ptr2> bool equal(const ExpansionIterator<X,Ref2,Ptr2>& i) const {
        return _p==i._p; }
    template<class Ref2,class Ptr2> difference_type distance_to(const ExpansionIterator<X,Ref2,Ptr2>& i) const {
        return (i._p-_p)/difference_type(_nw+_ds); }
    Ptr operator->() const {
        return reinterpret_cast<Ptr>(const_cast<Iter*>(this)); }
    Ref operator*() const {
        return reinterpret_cast<Ref>(const_cast<Iter&>(*this)); }
    Iter& increment() {
        return advance(1); }
    Iter& decrement() {
        return advance(-1); }
    Iter& advance(difference_type m) {
        _p+=m*difference_type(_nw+_ds);
        return *this; }
    const word_type* _word_ptr() { return this->_p; }
    const void* _ptr() { return this->_p; }
  private:
    size_type _n;
    size_type _nw;
    word_type* _p;
    static const int _ds=sizeof(data_type)/sizeof(word_type);
};


template<class X, class Ref, class Ptr>
inline std::ostream& operator<<(std::ostream& os, const ExpansionIterator<X,Ref,Ptr>& piter) {
    const MultiIndex::byte_type* ap=reinterpret_cast<const MultiIndex::byte_type*>(piter._p);
    const X* xp=reinterpret_cast<const X*>(piter._p+piter._nw);
    os << "(n="<<piter._n<<", nw="<<piter._nw<<", p="<<(void*)piter._p<<", a=";
    for(MultiIndex::size_type i=0; i!=piter._n; ++i) { os<<(i==0?"(":",")<<int(ap[i]); }
    os << "), x="<<*xp<<")";
    return os;
}

template<class X>
class Expansion
{
  public:
    static X _zero;
  public:
    typedef X RealType;
    typedef unsigned short int smoothness_type;
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::byte_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef int difference_type;
    typedef MultiIndex key_type;
    typedef X data_type;

    static const unsigned int sizeof_byte=sizeof(byte_type);
    static const unsigned int sizeof_word=sizeof(word_type);
    static const unsigned int sizeof_data=sizeof(data_type);
  public:
    typedef ExpansionValue<X> value_type;
    typedef ExpansionValue<X>& reference;
    typedef ExpansionValue<X>const& const_reference;
    typedef ExpansionValue<X>* pointer;
    typedef ExpansionValue<X>const* const_pointer;
    typedef ExpansionIterator<X,ExpansionValue<X>&, ExpansionValue<X>* > iterator;
    typedef ExpansionIterator<X,ExpansionValue<X>const&, ExpansionValue<X>const* > const_iterator;
  public:
    explicit Expansion() : _argument_size() { }
    explicit Expansion(unsigned int as) : _argument_size(as) { }
    template<class XX> explicit Expansion(unsigned int as, unsigned int deg, const XX* c0, ...);
    explicit Expansion(unsigned int as, unsigned int deg, double c0, ...);
    explicit Expansion(unsigned int as, unsigned int nnz, int a00, ...);
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    template<class XX> Expansion(const Expansion<XX>& p);

    Expansion<X>& operator=(const X& x);

    static Expansion<X> variable(unsigned int as, unsigned int i);

    void swap(Expansion<X>& other) {
        std::swap(this->_argument_size,other._argument_size);
        std::swap(this->_coefficients,other._coefficients); }

    bool operator==(const Expansion<X>& other) const {
        if(this->argument_size()!=other.argument_size()) { return false; }
        if(this->number_of_nonzeros()!=other.number_of_nonzeros()) { return false; }
        const_iterator self_iter=this->begin(); const_iterator other_iter=other.begin();
        while(self_iter!=this->end()) {
            if(self_iter->key()!=other_iter->key() || self_iter->data()!=other_iter->data()) { return false; }
            ++self_iter; ++other_iter; }
        return true;
    }
    bool operator!=(const Expansion<X>& other) const { return !this->operator==(other); }

    unsigned int argument_size() const { return this->_argument_size; }
    unsigned int number_of_nonzeros() const { return _coefficients.size()/_element_size(); }
    unsigned int degree() const {
        unsigned char deg=0u; for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
            deg=std::max(deg,iter->key().degree()); } return deg; }
    const std::vector<word_type>& coefficients() const { return this->_coefficients; }

    bool empty() const { return this->_coefficients.empty(); }
    unsigned int size() const { return this->_coefficients.size()/_element_size(); }
    void reserve(size_type nnz) { this->_coefficients.reserve(nnz*_element_size()); }
    void resize(size_type nnz) { this->_coefficients.resize(nnz*_element_size()); }

    template<class CMP> void insert(const MultiIndex& a, const RealType& c, const CMP& cmp) {
        this->_insert(a,c,cmp); }
    template<class CMP> void set(const MultiIndex& a, const RealType& c, const CMP& cmp) {
        this->_at(a,cmp)=c; }
    template<class CMP> RealType& at(const MultiIndex& a, const CMP& cmp) {
        return this->_at(a,cmp); }
    template<class CMP> iterator find(const MultiIndex& a, const CMP& cmp) {
        iterator iter=std::lower_bound(this->begin(),this->end(),a,cmp); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter; }
    template<class CMP> const_iterator find(const MultiIndex& a, const CMP& cmp) const {
        iterator iter=std::lower_bound(this->begin(),this->end(),a,cmp); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter; }

    void append(const MultiIndex& a, const RealType& c) {
        this->_append(a,c); }
    void prepend(const MultiIndex& a, const RealType& c) {
        this->_prepend(a,c); }
    void append(const MultiIndex& a1, const MultiIndex& a2, const RealType& c) {
        this->_append(a1,a2,c); }

    //RealType& operator[](const MultiIndex& a) {
    //    iterator iter=this->find(a);
    //    if(iter==this->end()) { iter=this->_insert(iter,a,static_cast<RealType>(0.0)); }
    //    return iter->data(); }
    const RealType& operator[](const MultiIndex& a) const {
        const_iterator iter=this->find(a);
        if(iter==this->end()) { return _zero; }
        else { return iter->data(); } }

    iterator begin() { return iterator(_argument_size,_index_size(),_begin_ptr()); }
    iterator end() { return iterator(_argument_size,_index_size(),_end_ptr()); }
    iterator find(const MultiIndex& a) {
        iterator iter=this->end(); while(iter!=this->begin()) { --iter; if(iter->key()==a) { return iter; } } return this->end(); }

    const_iterator begin() const { return const_iterator(_argument_size,_index_size(),const_cast<word_type*>(_begin_ptr())); }
    const_iterator end() const { return const_iterator(_argument_size,_index_size(),const_cast<word_type*>(_end_ptr())); }
    const_iterator find(const MultiIndex& a) const {
        const_iterator iter=this->end(); while(iter!=this->begin()) { --iter; if(iter->key()==a) { return iter; } } return this->end(); }

    void erase(iterator iter) { iter->data()=0.0; }
    void clear() { _coefficients.clear(); }

    void remove_zeros() {
        this->resize(std::remove_if(this->begin(),this->end(),DataIsZero())-this->begin()); }
    void combine_terms() {
        iterator curr=this->begin(); const_iterator adv=this->begin(); const_iterator end=this->end();
        while(adv!=end) { curr->key()=adv->key(); curr->data()=adv->data(); ++adv;
            while(adv!=end && curr->key()==adv->key()) { curr->data()+=adv->data(); ++adv; } ++curr; }
        this->resize(curr-this->begin()); }

    template<class CMP> void sort(const CMP& cmp) {
        std::sort(this->begin(),this->end(),cmp);
    }

    void graded_sort() {
        std::sort(this->begin(),this->end(),GradedKeyLess());
    }
    void lexicographic_sort() {
        std::sort(this->begin(),this->end(),LexicographicKeyLess());
    }
    void reverse_lexicographic_sort() {
        std::sort(this->begin(),this->end(),ReverseLexicographicKeyLess());
    }

    void check() const { }

    Expansion<X> embed(unsigned int before_size, const Expansion<X>&, unsigned int after_size) const;
  public:
    size_type _vector_size() const {
        return _coefficients.size(); }
    size_type _index_size() const {
        return (_argument_size+sizeof_word)/sizeof_word; }
    size_type _element_size() const {
        return (_argument_size+sizeof_word)/sizeof_word+sizeof_data/sizeof_word; }
  public:
    const word_type* _begin_ptr() const { return _coefficients.begin().operator->(); }
    const word_type* _end_ptr() const { return _coefficients.end().operator->(); }
    word_type* _begin_ptr() { return _coefficients.begin().operator->(); }
    word_type* _end_ptr() { return _coefficients.end().operator->(); }
    template<class CMP> RealType& _at(const MultiIndex& a,const CMP& cmp) {
        iterator p=std::lower_bound(this->begin(),this->end(),a,cmp);
        if(p!=this->end() && p->key()==a) { return p->data(); }
        else { p=this->_insert(p,a,0.0); return p->data(); }
    }
    template<class CMP> iterator _insert(const MultiIndex& a, const RealType& x,const CMP& cmp) {
        //std::cerr<<"_insert "<<*this<<" "<<a<<" "<<x<<std::endl;
        _coefficients.resize(_coefficients.size()+_element_size());
        iterator p=std::lower_bound(this->begin(),this->end()-1,a,cmp);
        return _allocated_insert(p,a,x); }
    iterator _insert(iterator p, const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_insert "<<*this<<" "<<p._ptr()<<" "<<a<<" "<<x<<std::endl;
        if(_coefficients.size()+_element_size()>_coefficients.capacity()) {
            difference_type i=p-begin();
            _coefficients.resize(_coefficients.size()+_element_size());
            p=begin()+i;
        } else {
            _coefficients.resize(_coefficients.size()+_element_size());
        }
        return _allocated_insert(p,a,x); }
    iterator _allocated_insert(iterator p, const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_allocated_insert "<<*this<<" "<<p<<" "<<p-begin()<<" "<<a<<" "<<x<<std::endl;
        iterator curr=this->end()-1; iterator prev=curr;
        while(curr!=p) { --prev; *curr=*prev; curr=prev; }
        curr->key()=a; curr->data()=x; return p; }
    void _prepend(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_prepend "<<*this<<" "<<a<<" "<<x<<std::endl;
        _coefficients.resize(_coefficients.size()+_element_size());
        _allocated_insert(begin(),a,x); }
    void _append(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_append "<<*this<<" "<<a<<" "<<x<<"... "<<std::flush;
        _coefficients.resize(_coefficients.size()+_element_size());
        word_type* vp=_end_ptr()-_element_size(); const word_type* ap=a.word_begin();
        for(unsigned int j=0; j!=_index_size(); ++j) { vp[j]=ap[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end_ptr())-1; *xp=x;
        //std::cerr<<"done"<<std::endl;
    }
    void _append(const MultiIndex&  a1, const MultiIndex&  a2, const RealType& x) {
        //std::cerr<<"_append "<<*this<<" "<<a1<<" "<<a2<<" "<<x<<std::endl;
        _coefficients.resize(_coefficients.size()+_element_size());
        word_type* vp=_end_ptr()-_element_size();
        const word_type* ap1=a1.word_begin(); const word_type* ap2=a2.word_begin();
        for(unsigned int j=0; j!=_index_size(); ++j) { vp[j]=ap1[j]+ap2[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end_ptr())-1; *xp=x; }
  public:
    std::ostream& write(std::ostream& os, const array<std::string>& variables) const;
  private:
    size_type _argument_size;
    std::vector<word_type> _coefficients;
};

// Disable construction of Expansion<Rational> since above implementation only
// works for "plain old data" types
#if defined HAVE_GMPXX_H and defined ARIADNE_NUMERIC_H
template<> class Expansion<Rational>;
#endif

template<class X> X Expansion<X>::_zero = 0.0;

template<class X> template<class XX>
Expansion<X>::Expansion(unsigned int as, unsigned int deg, const XX* p, ...)
    : _argument_size(as)
{
    MultiIndex a(as);
    while(a.degree()<=deg) {
        if(*p!=0) { this->insert(a,*p); }
        ++a; ++p;
    }
}

template<class X>
Expansion<X>::Expansion(unsigned int as, unsigned int deg, double c0, ...)
    : _argument_size(as)
{
    MultiIndex a(as); double x;
    va_list args; va_start(args,c0);
    while(a.degree()<=deg) {
        if(a.degree()==0) { x=c0; }
        else { x=va_arg(args,double); }
        if(x!=0) { this->append(a,x); }
        ++a;
    }
    va_end(args);
}

template<class X>
Expansion<X>::Expansion(unsigned int as, unsigned int nnz, int a00, ...)
    : _argument_size(as)
{
    MultiIndex a(as);
    double x;
    va_list args;
    va_start(args,a00);
    for(unsigned int i=0; i!=nnz; ++i) {
        for(unsigned int j=0; j!=as; ++j) {
            if(i==0 && j==0) { a[j]=a00; }
            else { a[j]=va_arg(args,int); }
        }
        x=va_arg(args,double);
        if(x!=0) { this->append(a,x); }
    }
    va_end(args);
}

template<class X> template<class XX>
Expansion<X>::Expansion(const std::map<MultiIndex,XX>& m)
{
    ARIADNE_ASSERT(!m.empty());
    this->_argument_size=m.begin()->first.size();
    for(typename std::map<MultiIndex,XX>::const_iterator iter=m.begin(); iter!=m.end(); ++iter) {
        this->append(iter->first,X(iter->second));
    }
}

template<class X> template<class XX> inline
Expansion<X>::Expansion(const Expansion<XX>& p)
    : _argument_size(p.argument_size())
{
    for(typename Expansion<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->key(),X(iter->data())); }
}


template<class X> Expansion<X>& Expansion<X>::operator=(const X& c) {
    this->clear();
    this->append(MultiIndex::zero(this->argument_size()),c);
    return *this;
}

template<class X> Expansion<X> Expansion<X>::variable(unsigned int n, unsigned int i) {
    Expansion<X> p(n); p.append(MultiIndex::unit(n,i),X(1));
    return p;
}




//! \ingroup FunctionModule
//! \brief Convert a power-series expansion into a formula using a version of Horner's rule.
//!
//! For a polynomial in \f$n\f$ variables, Horner's rule is a recursive formula
//! \f[ p(x) = \bigl( \bigl(  x^{d_k-d_{k-1}} q_k(\hat{x})x^{d_0} + \cdots + q_1(\hat{x}) \bigr) x^{d_1-d_0} + q_0(\hat{x}) \bigr) x^{d_0} \f]
//! where \f$\hat{x}=(x_1,\ldots,x_{n-1})\f$ and \f$q_i\f$ is the polynomial of terms in \f$x_n^{d_i}\f$.
//! To evaluate a polynomial using Horner's rule without using recursive function calls, we maintain registers \f$r_k\f$ containing
//! the current evaluation of a polynomial in \f$(x_1,\ldots,x_k)\f$.
//!
//! We list the terms in reverse lexicographic order, defined as \f$\alpha \prec \beta\f$ if \f$\alpha_j>\beta_j\f$,
//! where \f$j=\max\{i\mid \alpha_i\neq\beta_i\}$.\f$c_\alpha\f$.
//! For a given term \f$c_\alpha x^\alpha\f$, let \f$k=\max\{j\mid \alpha_j\neq\beta_j\}\f$, where \f$\beta\f$ is the next multi-index.
//! We update register \f$r_k\f$ by \f[r'_k=(((c_\alpha + r_1) x^{\alpha_1} + r_2 )x^{\alpha_2}+\cdots r_k)x^{\alpha_k-\beta_k}.\f]
//! The result is obtained by updating a fictional register \f$r_{n+1}\f$ at the last step.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
template<class X, class Y> Y horner_evaluate(const Expansion<X>& e, const Vector<Y>& x)
{
    typedef typename Expansion<X>::const_iterator const_iterator;

    const uint n=e.argument_size();
    Y z=x[0]*0; // The zero element of the ring Y
    array< Y > r(e.argument_size(),z); // An array of "registers" containing working p(x[0],...,x[k])
    const_iterator iter=e.begin();
    const_iterator end=e.end();
    uint k=n;   // The current working register
    const uchar* na=iter->key().begin(); // The values of the next multi-index
    uint j=k;   // The lowest register containing a non-zero value
    X c=iter->data();
    const uchar* a=na;
    ++iter;
    while(iter!=end) {
        na=iter->key().begin();
        k=n-1;
        while(a[k]==na[k]) { --k; }
        // Since terms are ordered reverse-lexicographically,
        // previous index must have higher kth value
        assert(a[k]>na[k]);
        // Set r[k]=(((c+r[0])*x[0]^a[0]+r[1])*x[1]^a[1]+...+r[k])*x[k]^(a[k]-na[k])
        // Omit zero terms where possible
        Y t=c;
        for(uint i=0; i!=min(j,k); ++i) {
            for(uint ii=0; ii!=a[i]; ++ii) {
                t=t*x[i];
            }
        }
        for(uint i=min(j,k); i!=k; ++i) {
            t=t+r[i];
            for(uint ii=0; ii!=a[i]; ++ii) {
                t=t*x[i];
            }
            r[i]=z;
        }
        if(j<=k) {
            t=t+r[k];
        }
        for(uint ii=na[k]; ii!=a[k]; ++ii) {
            t=t*x[k];
        }
        r[k]=t;
        //std::cerr<<"a="<<MultiIndex(n,a)<<" c="<<c<<" k="<<k<<" r="<<r<<"\n";
        j=k;
        c=iter->data();
        a=na;
        ++iter;
    }
    // Set r=(((c+r[0])*x[0]^a[0]+r[1])*x[1]^a[1]+...+r[n-1])*x[n-1]^(a[n-1])
    Y t=c;
    for(uint i=0; i!=j; ++i) {
        for(uint ii=0; ii!=a[i]; ++ii) {
            t=t*x[i];
        }
    }
    for(uint i=j; i!=n; ++i) {
        t=t+r[i];
        for(uint ii=0; ii!=a[i]; ++ii) {
            t=t*x[i];
        }
    }
    //std::cerr<<"a="<<MultiIndex(n,a)<<" c="<<c<<" k="<<n<<"\n";
    //std::cerr<<"  r="<<t<<"\n";
    return t;
}

template<class X, class Y>
Y power_evaluate(const Expansion<X>& e, const Vector<Y>& y)
{
    Y zero = y[0]; zero*=0;
    Y one = zero; one+=1;

    Y r=zero;
    Y t=zero;
    for(typename Expansion<X>::const_iterator iter=e.begin();
        iter!=e.end(); ++iter)
    {
        const MultiIndex& j=iter->key();
        const X& c=iter->data();
        t=one;
        for(uint k=0; k!=e.argument_size(); ++k) {
            for(uint l=0; l!=j[k]; ++l) {
                t=t*y[k];
            }
        }
        t*=c;
        r+=t;
    }

    return r;
}


template<class X, class Y>
Y evaluate(const Expansion<X>& e, const Vector<Y>& y)
{
    return power_evaluate(e,y);
}

template<class X, class Y>
Y simple_evaluate(const Expansion<X>& e, const Vector<Y>& y)
{
    return power_evaluate(e,y);
}

template<class X>
Expansion<X> embed(unsigned int before_size, const Expansion<X>& x, unsigned int after_size)
{
    uint old_size=x.argument_size();
    uint new_size=before_size+old_size+after_size;
    Expansion<X> r(new_size);
    MultiIndex old_index(old_size);
    MultiIndex new_index(new_size);
    for(typename Expansion<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        old_index=iter->key();
        for(uint j=0; j!=old_size; ++j) {
            uint aj=old_index[j];
            new_index[j+before_size]=aj;
        }
        r.append(new_index,iter->data());
    }
    return r;
}


inline Expansion<Float> midpoint(const Expansion<Interval>& pse) {
    Expansion<Float> r(pse.argument_size());
    for(Expansion<Interval>::const_iterator iter=pse.begin(); iter!=pse.end(); ++iter) {
        r.append(iter->key(),midpoint(iter->data())); }
    return r;
}


template<class X>
std::ostream& Expansion<X>::write(std::ostream& os, const array<std::string>& variable_names) const
{
    ARIADNE_ASSERT(this->argument_size()==variable_names.size());
    const Expansion<X>& p=*this;
    if(p.size()==0) {
        os << "0";
    } else {
        bool first_term=true;
        for(const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
            MultiIndex a=iter->key();
            X v=iter->data();
            if(v>=0 && !first_term) { os<<"+"; }
            first_term=false;
            bool first_factor=true;
            if(v<0) { os<<"-"; }
            if(abs(v)!=1 || a.degree()==0) { os<<abs(v); first_factor=false; }
            for(uint j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<variable_names[j]; if(a[j]!=1) { os<<"^"<<int(a[j]); }
                }
            }
        }
    }
    return os;
}


template<class X>
std::ostream& operator<<(std::ostream& os, const Expansion<X>& p) {
    array<std::string> variable_names(p.argument_size());
    for(uint j=0; j!=p.argument_size(); ++j) {
        std::stringstream sstr;
        sstr << 'x' << j;
        variable_names[j]=sstr.str();
    }
    return p.write(os,variable_names);
}




template<class X, class Y>
Vector<Y> evaluate(const Vector< Expansion<X> >& x, const Vector<Y>& y)
{
    Vector<Y> r(x.size());
    for(unsigned int i=0; i!=x.size(); ++i) {
        r[i]=evaluate(x[i],y);
    }
    return r;
}


template<class X> Vector< Expansion<X> > operator*(const Expansion<X>& e, const Vector<Float> v) {
    Vector< Expansion<X> > r(v.size(),Expansion<X>(e.argument_size()));
    for(uint i=0; i!=r.size(); ++i) {
        ARIADNE_ASSERT(v[i]==0.0 || v[i]==1.0);
        if(v[i]==1.0) { r[i]=e; }
    }
    return r;
}


inline Vector< Expansion<Float> > midpoint(const Vector< Expansion<Interval> >& pse) {
    Vector< Expansion<Float> > r(pse.size(),Expansion<Float>());
    for(uint i=0; i!=pse.size(); ++i) {
        r[i]=midpoint(pse[i]); }
    return r;
}








}

#endif /* ARIADNE_EXPANSION_H */
