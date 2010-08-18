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



template<class V, class K=typename V::key_type> struct key_equals {
    bool operator()(const V& v, const K& k) const { return v.key()==k; }
};

template<class V, class K=typename V::key_type> struct key_less {
    bool operator()(const V& v, const K& k) const { return v.key()<k; }
};

template<class V> struct data_is_zero {
    bool operator()(const V& v) const { return v.data()==0; }
};

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


template<class X> struct key_less<ExpansionValue<X>,MultiIndex> {
    bool operator()(const ExpansionValue<X>& v, const MultiIndex& k) const { return v.key()<k; }
};


template<class X>
inline bool operator<(const ExpansionValue<X>& dv1, const ExpansionValue<X>& dv2) {
    return dv1.key()<dv2.key();
}

template<class X>
inline bool data_less(const ExpansionValue<X>& dv1, const ExpansionValue<X>& dv2) {
    return dv1.data()<dv2.data();
}


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
    explicit Expansion(unsigned int as, unsigned int deg, double c0, ...);
    explicit Expansion(unsigned int as, unsigned int nnz, int a00, ...);
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    template<class XX> Expansion(const Expansion<XX>& p);

    Expansion<X>& operator=(const X& x);

    static Expansion<X> variable(unsigned int as, unsigned int i);

    void swap(Expansion<X>& p) {
        std::swap(this->_argument_size,p._argument_size);
        std::swap(this->_coefficients,p._coefficients); }

    bool operator==(const Expansion<X>& p) const { return this->_coefficients == p._coefficients; }
    bool operator!=(const Expansion<X>& p) const { return this->_coefficients != p._coefficients; }

    unsigned int argument_size() const { return this->_argument_size; }
    unsigned int number_of_nonzeros() const { return _coefficients.size()/_element_size(); }
    unsigned int degree() const { if(this->empty()) { return 0u; } return (--this->end())->key().degree(); }
    const std::vector<word_type>& coefficients() const { return this->_coefficients; }

    bool empty() const { return this->_coefficients.empty(); }
    unsigned int size() const { return this->_coefficients.size()/_element_size(); }
    void reserve(size_type nnz) { this->_coefficients.reserve(nnz*_element_size()); }
    void resize(size_type nnz) { this->_coefficients.resize(nnz*_element_size()); }

    void insert(const MultiIndex& a, const RealType& c); // Non-inline user function
    void append(const MultiIndex& a, const RealType& c) {
        this->_append(a,c); }
    void prepend(const MultiIndex& a, const RealType& c) {
        this->_prepend(a,c); }
    void append(const MultiIndex& a1, const MultiIndex& a2, const RealType& c) {
        this->_append(a1,a2,c); }

    RealType& operator[](const MultiIndex& a) {
        iterator iter=this->lower_bound(a);
        if(iter==this->end() || iter->key()!=a) { iter=this->_insert(iter,a,static_cast<RealType>(0.0)); }
        return iter->data(); }
    const RealType& operator[](const MultiIndex& a) const {
        const_iterator iter=this->lower_bound(a);
        if(iter==this->end() || iter->key()!=a) { return _zero; }
        else { return iter->data(); } }

    iterator begin() { return iterator(_argument_size,_index_size(),_begin_ptr()); }
    iterator end() { return iterator(_argument_size,_index_size(),_end_ptr()); }
    iterator lower_bound(const MultiIndex& a) {
        return std::lower_bound(this->begin(),this->end(),a,key_less<value_type,key_type>()); }
    iterator find(const MultiIndex& a) {
        iterator iter=this->lower_bound(a); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter; }

    const_iterator begin() const { return const_iterator(_argument_size,_index_size(),const_cast<word_type*>(_begin_ptr())); }
    const_iterator end() const { return const_iterator(_argument_size,_index_size(),const_cast<word_type*>(_end_ptr())); }
    const_iterator lower_bound(const MultiIndex& a) const {
        return std::lower_bound(this->begin(),this->end(),a,key_less<value_type,key_type>()); }
    const_iterator find(const MultiIndex& a) const {
        const_iterator iter=this->lower_bound(a); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter; }

    void erase(iterator iter) { iter->data()=0.0; }
    void clear() { _coefficients.clear(); }

    void sort() {
        //std::cerr<<"sorting... "<<std::flush;
        std::sort(this->begin(),this->end());
    }
    void remove_zeros() {
        //std::cerr<<"remove_zeros... "<<std::flush;
        iterator new_end=std::remove_if(this->begin(),this->end(),data_is_zero<value_type>());
        this->resize(new_end-this->begin());
    }


    void cleanup() { this->sort(); this->remove_zeros(); }

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
    iterator _insert(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_insert "<<*this<<" "<<a<<" "<<x<<std::endl;
        _coefficients.resize(_coefficients.size()+_element_size());
        iterator p=std::lower_bound(this->begin(),this->end()-1,a,key_less<value_type,key_type>());
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

template<class X>
void
Expansion<X>::insert(const MultiIndex& a, const X& x) {
    this->_insert(a,x);
}

// Disable construction of Expansion<Rational> since above implementation only
// works for "plain old data" types
#if defined HAVE_GMPXX_H and defined ARIADNE_NUMERIC_H
template<> class Expansion<Rational>;
#endif

template<class X> X Expansion<X>::_zero = 0.0;

template<class X>
Expansion<X>::Expansion(unsigned int as, unsigned int deg, double c0, ...)
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
    this->cleanup();
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
    this->cleanup();
}

template<class X> template<class XX>
Expansion<X>::Expansion(const std::map<MultiIndex,XX>& m)
{
    ARIADNE_ASSERT(!m.empty());
    this->_argument_size=m.begin()->first.size();
    for(typename std::map<MultiIndex,XX>::const_iterator iter=m.begin(); iter!=m.end(); ++iter) {
        this->append(iter->first,X(iter->second));
    }
    this->remove_zeros();
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
    (*this)[MultiIndex::zero(this->argument_size())]=c; 
    return *this;
}

template<class X> Expansion<X> Expansion<X>::variable(unsigned int n, unsigned int i) {
    Expansion<X> p(n); p[MultiIndex::zero(n)]=0.0; p[MultiIndex::unit(n,i)]=X(1);
    return p;
}




template<class X, class Y>
Y evaluate(const Expansion<X>& x, const Vector<Y>& y)
{
    Y zero = y[0]; zero*=0;
    Y one = zero; one+=1;

    Y r=zero;
    Y t=zero;
    for(typename Expansion<X>::const_iterator iter=x.begin();
        iter!=x.end(); ++iter)
    {
        const MultiIndex& j=iter->key();
        const X& c=iter->data();
        t=one;
        for(uint k=0; k!=x.argument_size(); ++k) {
            for(uint l=0; l!=j[k]; ++l) {
                t=t*y[k];
            }
        }
        t*=c;
        r+=t;
    }

    return r;
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
