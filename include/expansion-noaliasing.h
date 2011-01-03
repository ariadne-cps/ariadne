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
#include <boost/iterator/zip_iterator.hpp>
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




template<class Key, class Data> class KeyDataValue;
template<class Key, class Data> class KeyDataReference;

template<class Key, class Data> class KeyDataValue {
  public:
    template<class K, class D> KeyDataValue(const KeyDataValue<K,D>& m) : _k(m.key()), _d(m.data()) { }
    template<class K, class D> KeyDataValue(const KeyDataReference<K,D>& m) : _k(m.key()), _d(m.data()) { }
    KeyDataValue(Key k, Data d) : _k(k), _d(d) { }
    Key& key() { return _k; }
    const Key& key() const { return _k; }
    Data& data() { return _d; }
    const Data& data() const { return _d; }
  private:
  public:
    Key _k; Data _d;
};

template<class KeyPtr, class DataPtr> class KeyDataReference {
  public:
    typedef typename KeyPtr::value_type Key;
    typedef typename DataPtr::value_type Data;

    KeyDataReference() : _kp(), _dp() { }
    KeyDataReference(KeyPtr ki, DataPtr di) : _kp(ki), _dp(di) { }
    KeyDataReference(const KeyDataReference<KeyPtr,DataPtr>& m) : _kp(m._kp), _dp(m._dp) { }
    template<class KP, class DP> KeyDataReference(const KeyDataReference<KP,DP>& m) : _kp(m._kp), _dp(m._dp) { }
    //template<class K, class D> KeyDataReference<KeyPtr,DataPtr>& operator=(KeyDataReference<K,D> ref) { *_kp=*ref._kp; *_dp=ref._dp; return *this; }
    template<class KP, class DP> KeyDataReference<KeyPtr,DataPtr>& operator=(const KeyDataReference<KP,DP>& ref) { *_kp=*ref._kp; *_dp=*ref._dp; return *this; }
    KeyDataReference<KeyPtr,DataPtr>& operator=(const KeyDataReference<KeyPtr,DataPtr>& ref) {
        *_kp=*ref._kp; *_dp=*ref._dp; return *this; }
    KeyDataReference<KeyPtr,DataPtr>& operator=(const KeyDataValue<Key,Data>& val) { *_kp=val._k; *_dp=val._d; return *this; }
    template<class KP,class DP> bool operator==(const KeyDataReference<KP,DP>& ref) const { return this->key()==ref.key() && this->data()==ref.data(); }
    template<class K,class D> bool operator==(const KeyDataValue<K,D>& val) const { return this->key()==val.key() && this->data()==val.data(); }
    typename KeyPtr::reference key() { return *_kp; }
    typename DataPtr::reference data()  { return *_dp; }
    typename KeyPtr::reference key() const { return *_kp; }
    typename DataPtr::reference data() const { return *_dp; }
    mutable KeyPtr _kp; mutable DataPtr _dp;
};



template<class KeyIter, class DataIter>
class KeyDataIterator
    : public boost::iterator_facade<
        KeyDataIterator<KeyIter,DataIter>,KeyDataValue<typename KeyIter::value_type,typename DataIter::value_type>,
        boost::random_access_traversal_tag,KeyDataReference<KeyIter,DataIter> >
{
  public:
    typedef typename KeyIter::reference key_reference;
    typedef typename DataIter::reference data_reference;
    typedef KeyDataReference<KeyIter,DataIter> reference;
    typedef KeyDataReference<KeyIter,DataIter>* pointer;

    KeyIter key_iter() { return _ref._kp; }
    DataIter data_iter() { return _ref._dp; }

    KeyDataIterator() : _ref() { }
    KeyDataIterator(const KeyIter& kiter, const DataIter& diter) : _ref(kiter,diter) { }
    KeyDataIterator(const KeyDataIterator<KeyIter,DataIter>& other) : _ref(other._ref) { }
    KeyDataIterator& operator=(const KeyDataIterator<KeyIter,DataIter>& other) { _ref._kp=other._ref._kp; _ref._dp=other._ref._dp; return *this; }
    template<class KI, class DI> KeyDataIterator(const KeyDataIterator<KI,DI>& other) : _ref(other._ref) { }
    template<class Iter> bool equal(const Iter& other) const { return this->_ref._kp==other._ref._kp; }
    template<class Iter> bool operator==(const Iter& other) const { return this->_ref._kp==other._ref._kp; }
    template<class Iter> bool operator!=(const Iter& other) const { return this->_ref._kp!=other._ref._kp; }
    template<class Iter> int distance_to(const Iter& other) const { return other._ref._dp-this->_ref._dp; }
    void increment() { advance(1); }
    void decrement() { advance(-1); }
    void advance(int m) { _ref._kp+=m; _ref._dp+=m; }
    reference dereference() const { return _ref; }
    pointer operator->() const { return const_cast<pointer>(&_ref); }
  private:
  public:
    KeyDataReference<KeyIter,DataIter> _ref;
};

template<class K,class D>
inline bool operator<(const KeyDataReference<K,D>& dv1, const KeyDataReference<K,D>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K,class D, class KP, class DP>
inline bool operator<(const KeyDataReference<KP,DP>& dv1, const KeyDataValue<K,D>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K,class D>
inline bool operator<(const KeyDataValue<K,D>& dv1, const KeyDataValue<K,D>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K,class D, class KP, class DP>
inline bool operator<(const KeyDataValue<K,D>& dv1, const KeyDataReference<KP,DP>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K, class D>
std::ostream& operator<<(std::ostream& os, const boost::Tuple<K,D>& tup) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "(" << "??" << ")";
}

template<class K, class D>
std::ostream& operator<<(std::ostream& os, const boost::tuples::cons<K,D>& tup) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "(" << "??" << ")";
}

template<class T>
std::ostream& operator<<(std::ostream& os, boost::zip_iterator< T > iter) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<iter_traversal>";
}

template<class K, class D>
std::ostream& operator<<(std::ostream& os, KeyDataValue<K,D> val) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << val.key() << ":" << val.data();
}

template<class K, class D>
std::ostream& operator<<(std::ostream& os, KeyDataReference<K,D> ref) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << ref.key() << ":" << ref.data();
}

template<class K, class D>
std::ostream& operator<<(std::ostream& os, KeyDataIterator<K,D> iter) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<"<<static_cast<const void*>(iter->key().begin())<<","<<&iter->data()<<">("<<iter->key()<<","<<iter->data()<<")";
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
    typedef MultiIndex::index_type index_type;
    typedef MultiIndex::raw_data_type raw_data_type;
    typedef int difference_type;

    typedef MultiIndex key_type;
    typedef X mapped_type;
    typedef X scalar_type;

    static const unsigned int sizeof_byte=sizeof(byte_type);
    static const unsigned int sizeof_word=sizeof(word_type);
    static const unsigned int sizeof_mapped=sizeof(mapped_type);
  public:
    //typedef ExpansionValue<X> value_type;
    //typedef ExpansionValueReference<X> reference;
    //typedef const ExpansionValueReference<X> const_reference;
    //typedef ExpansionValueReference<X>* pointer;
    //typedef const ExpansionValueReference<X>* const_pointer;
    //typedef ExpansionIterator<X,ExpansionValueReference<X> > iterator;
    //typedef ExpansionIterator<X,ExpansionValueReference<X> > const_iterator;
    typedef MultiIndexList::iterator key_iterator;
    typedef MultiIndexList::const_iterator key_const_iterator;
    typedef typename std::vector<X>::iterator data_iterator;
    typedef typename std::vector<X>::const_iterator data_const_iterator;

    typedef KeyDataIterator<key_iterator,data_iterator> iterator;
    typedef KeyDataIterator<key_const_iterator,data_const_iterator> const_iterator;

    typedef typename iterator::value_type value_type;
    typedef typename iterator::reference reference;
    typedef typename const_iterator::reference const_reference;
    typedef typename iterator::pointer pointer;
    typedef typename const_iterator::pointer const_pointer;
  public:
    Expansion() : _indices(0), _coefficients() { }
    Expansion(unsigned int as) : _indices(as), _coefficients() { }
    Expansion(unsigned int as, unsigned int deg, double c0, ...);
    Expansion(unsigned int as, unsigned int nnz, int a00, ...);
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    template<class XX> Expansion(const Expansion<XX>& p);

    static Expansion<X> variable(unsigned int as, unsigned int i);

    void swap(Expansion<X>& other) {
        std::swap(this->_indices,other._indices);
        std::swap(this->_coefficients,other._coefficients); }

    bool operator==(const Expansion<X>& other) const { return this->_coefficients == other._coefficients && this->_indices==other._indices; }
    bool operator!=(const Expansion<X>& other) const { return !(*this==other); }

    unsigned int argument_size() const { return this->_indices.element_size(); }
    unsigned int number_of_nonzeros() const { return _coefficients.size(); }
    unsigned int degree() const { if(this->empty()) { return 0u; } return this->_indices.back().degree(); }

    bool empty() const { return this->_coefficients.empty(); }
    unsigned int size() const { return this->_coefficients.size(); }
    unsigned int capacity() const { return this->_coefficients.capacity(); }
    void reserve(size_type nnz) { this->_indices.reserve(nnz); this->_coefficients.reserve(nnz); }
    void resize(size_type nnz) { this->_indices.resize(nnz); this->_coefficients.resize(nnz); }

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

    iterator begin() { return iterator(_indices.begin(),_coefficients.begin()); }
    iterator end() { return iterator(_indices.end(),_coefficients.end()); }
    iterator lower_bound(const MultiIndex& a) {
        return std::lower_bound(this->begin(),this->end(),a,key_less<value_type,key_type>()); }
    iterator find(const MultiIndex& a) {
        iterator iter=this->lower_bound(a); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter; }

    const_iterator begin() const { return const_iterator(_indices.begin(),_coefficients.begin()); }
    const_iterator end() const { return const_iterator(_indices.end(),_coefficients.end()); }
    const_iterator lower_bound(const MultiIndex& a) const {
        return std::lower_bound(this->begin(),this->end(),a,key_less<value_type,key_type>()); }
    const_iterator find(const MultiIndex& a) const {
        const_iterator iter=this->lower_bound(a); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter; }

    reference front() { return reference(_indices.begin(),_coefficients.begin()); }
    const_reference front() const { return const_reference(_indices.begin(),_coefficients.begin()); }
    reference back() { return reference(_indices.end()-1,_coefficients.end()-1); }
    const_reference back() const { return const_reference(_indices.end()-1,_coefficients.end()-1); }

    void erase(iterator iter) { iter->data()=0.0; }
    void clear() { _indices.clear(); _coefficients.clear(); }

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
    iterator _insert(iterator p, const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_insert "<<*this<<" "<<p._ptr()<<" "<<a<<" "<<x<<std::endl;
        if(this->size()+1u>this->capacity()) {
            difference_type i=p-begin();
            this->resize(this->size()+1);
            p=begin()+i;
        } else {
            this->resize(this->size()+1);
        }
        return this->_allocated_insert(p,a,x); }
    iterator _insert(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_insert "<<*this<<" "<<a<<" "<<x<<std::endl;
        this->resize(this->size()+1u);
        iterator p=std::lower_bound(this->begin(),this->end()-1,a,key_less<value_type,key_type>());
        return this->_allocated_insert(p,a,x); }
    iterator _allocated_insert(iterator p, const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_allocated_insert "<<*this<<" "<<p<<" "<<p-begin()<<" "<<a<<" "<<x<<std::endl;
        iterator curr=this->end(); --curr; iterator prev=curr;
        while(curr!=p) { --prev; *curr=*prev; curr=prev; }
        curr->key()=a; curr->data()=x; return p; }
    void _prepend(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_prepend "<<*this<<" "<<a<<" "<<x<<std::endl;
        this->resize(this->size()+1u);
        _allocated_insert(this->begin(),a,x); }
    void _append(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_append "<<*this<<" "<<a<<" "<<x<<"... "<<std::flush;
        this->resize(this->size()+1u);
        this->back().key()=a; this->back().data()=x; }
    void _append(const MultiIndex&  a1, const MultiIndex&  a2, const RealType& x) {
        //std::cerr<<"_append "<<*this<<" "<<a1<<" "<<a2<<" "<<x<<std::endl;
        this->resize(this->size()+1u);
        this->back().key()=a1; this->back().key()+=a2; this->back().data()=x; }
  public:
    std::ostream& write(std::ostream& os, const Array<std::string>& variables) const;
  private:
  public:
    MultiIndexList _indices;
    std::vector<scalar_type> _coefficients;
};


/*
template<class X>
std::ostream& operator<<(std::ostream& os, const typename Expansion<X>::const_iterator& citer) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<const_iterator>";
}

template<class X>
std::ostream& operator<<(std::ostream& os, const typename Expansion<X>::reference& ref) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<reference>";
}

std::ostream& operator<<(std::ostream& os, const Expansion<Float>::const_reference& cref) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<const_reference>";
}
*/

template<class X>
void
Expansion<X>::insert(const MultiIndex& a, const X& x) {
    this->_insert(a,x);
}


template<class X> X Expansion<X>::_zero = 0.0;

template<class X>
Expansion<X>::Expansion(unsigned int as, unsigned int deg, double c0, ...)
    : _indices(as), _coefficients()
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
    : _indices(as), _coefficients()
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
    this->_indices=MultiIndexList(m.begin()->first.size());
    for(typename std::map<MultiIndex,XX>::const_iterator iter=m.begin(); iter!=m.end(); ++iter) {
        this->append(iter->first,X(iter->second));
    }
    this->remove_zeros();
}

template<class X> template<class XX> inline
Expansion<X>::Expansion(const Expansion<XX>& p)
    : _indices(p.argument_size()), _coefficients()
{
    for(typename Expansion<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->key(),X(iter->data())); }
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
        //r.append(iter->key(),midpoint(iter->data())); }
        r.append(iter->key(),midpoint(iter->data())); }
    return r;
}


template<class X>
std::ostream& Expansion<X>::write(std::ostream& os, const Array<std::string>& variable_names) const
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
            if(v!=0) {
                if(v>0 && !first_term) { os<<"+"; }
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
    }
    return os;
}


template<class X>
std::ostream& operator<<(std::ostream& os, const Expansion<X>& p) {
    Array<std::string> variable_names(p.argument_size());
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
    Vector< Expansion<Float> > r(pse.size());
    for(uint i=0; i!=pse.size(); ++i) {
        r[i]=midpoint(pse[i]); }
    return r;
}








}

#endif /* ARIADNE_EXPANSION_H */
