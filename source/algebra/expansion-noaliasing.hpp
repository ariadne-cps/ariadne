/***************************************************************************
 *            expansion.hpp
 *
 *  Copyright 2008-17 Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file expansion.hpp
 *  \brief Base class for power series expansions.
 */

// Use vector-based expansion and casting of iterators to get a reference

#ifndef ARIADNE_EXPANSION_HPP
#define ARIADNE_EXPANSION_HPP

#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "../utility/iterator.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/multi_index.hpp"



namespace Ariadne {

template<class X> class Expansion;

template<class X> Expansion<MultiIndex,X> embed(unsigned int, const Expansion<MultiIndex,X>&, unsigned int);



template<class V, class K=typename V::key_type> struct key_equals {
    Bool operator()(const V& v, const K& k) const { return v.key()==k; }
};

template<class V, class K=typename V::key_type> struct key_less {
    Bool operator()(const V& v, const K& k) const { return v.key()<k; }
};

template<class V> struct data_is_zero {
    Bool operator()(const V& v) const { return v.data()==0; }
};

template<class FwdIter, class Op>
FwdIter unique_key(FwdIter first, FwdIter last, Op op) {
    FwdIter curr=first;
    FwdIter next=curr;
    while(next!=last) {
        if(curr!=next) { *curr=*next; }
        ++next;
        while(next!=last && curr->index()==next->index()) {
            if(curr->index()==next->index()) {
                curr->coefficient()=op(curr->coefficient(),next->coefficient());
                ++next;
            }
        }
        // Removes zero entries; the code below is preferred to the case "curr->coefficient()!=0" for ValidatedKleenean results
        if(curr->coefficient()==0) { }
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
    template<class KP,class DP> Bool operator==(const KeyDataReference<KP,DP>& ref) const { return this->index()==ref.key() && this->coefficient()==ref.data(); }
    template<class K,class D> Bool operator==(const KeyDataValue<K,D>& val) const { return this->index()==val.key() && this->coefficient()==val.data(); }
    typename KeyPtr::reference key() { return *_kp; }
    typename DataPtr::reference data()  { return *_dp; }
    typename KeyPtr::reference key() const { return *_kp; }
    typename DataPtr::reference data() const { return *_dp; }
    mutable KeyPtr _kp; mutable DataPtr _dp;
};



template<class KeyIter, class DataIter>
class KeyDataIterator
    : public IteratorFacade<
        KeyDataIterator<KeyIter,DataIter>,KeyDataValue<typename KeyIter::value_type,typename DataIter::value_type>,
        RandomAccessTraversalTag,KeyDataReference<KeyIter,DataIter> >
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
    template<class Iter> Bool equal(const Iter& other) const { return this->_ref._kp==other._ref._kp; }
    template<class Iter> Bool operator==(const Iter& other) const { return this->_ref._kp==other._ref._kp; }
    template<class Iter> Bool operator!=(const Iter& other) const { return this->_ref._kp!=other._ref._kp; }
    template<class Iter> Int distance_to(const Iter& other) const { return other._ref._dp-this->_ref._dp; }
    Void increment() { advance(1); }
    Void decrement() { advance(-1); }
    Void advance(Int m) { _ref._kp+=m; _ref._dp+=m; }
    reference dereference() const { return _ref; }
    pointer operator->() const { return const_cast<pointer>(&_ref); }
  private:
  public:
    KeyDataReference<KeyIter,DataIter> _ref;
};

template<class K,class D>
inline Bool operator<(const KeyDataReference<K,D>& dv1, const KeyDataReference<K,D>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K,class D, class KP, class DP>
inline Bool operator<(const KeyDataReference<KP,DP>& dv1, const KeyDataValue<K,D>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K,class D>
inline Bool operator<(const KeyDataValue<K,D>& dv1, const KeyDataValue<K,D>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K,class D, class KP, class DP>
inline Bool operator<(const KeyDataValue<K,D>& dv1, const KeyDataReference<KP,DP>& dv2) {
    return dv1.key()<dv2.key();
}

template<class K, class D>
OutputStream& operator<<(OutputStream& os, const Tuple<K,D>& tup) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "(" << "??" << ")";
}

template<class K, class D>
OutputStream& operator<<(OutputStream& os, KeyDataValue<K,D> val) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << val.key() << ":" << val.data();
}

template<class K, class D>
OutputStream& operator<<(OutputStream& os, KeyDataReference<K,D> ref) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << ref.key() << ":" << ref.data();
}

template<class K, class D>
OutputStream& operator<<(OutputStream& os, KeyDataIterator<K,D> iter) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<"<<static_cast<const Void*>(iter->index().begin())<<","<<&iter->coefficient()<<">("<<iter->index()<<","<<iter->coefficient()<<")";
}




template<class X>
class Expansion
{
  public:
    static X _zero;
  public:
    typedef X RealType;
    typedef unsigned short Int smoothness_type;
    typedef MultiIndex::SizeType SizeType;
    typedef MultiIndex::ByteType ByteType;
    typedef MultiIndex::WordType WordType;
    typedef MultiIndex::IndexType IndexType;
    typedef MultiIndex::raw_data_type raw_data_type;
    typedef Int difference_type;

    typedef MultiIndex key_type;
    typedef X mapped_type;
    typedef X scalar_type;

    static const unsigned int sizeof_byte=sizeof(ByteType);
    static const unsigned int sizeof_word=sizeof(WordType);
    static const unsigned int sizeof_mapped=sizeof(mapped_type);
  public:
    //typedef ExpansionValue<X> value_type;
    //typedef ExpansionValueReference<X> reference;
    //typedef const ExpansionValueReference<X> const_reference;
    //typedef ExpansionValueReference<X>* pointer;
    //typedef const ExpansionValueReference<X>* const_pointer;
    //typedef ExpansionIterator<X,ExpansionValueReference<X> > Iterator;
    //typedef ExpansionIterator<X,ExpansionValueReference<X> > ConstIterator;
    typedef MultiIndexList::Iterator key_iterator;
    typedef MultiIndexList::ConstIterator key_const_iterator;
    typedef typename std::vector<X>::Iterator data_iterator;
    typedef typename std::vector<X>::ConstIterator data_const_iterator;

    typedef KeyDataIterator<key_iterator,data_iterator> Iterator;
    typedef KeyDataIterator<key_const_iterator,data_const_iterator> ConstIterator;

    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::reference reference;
    typedef typename ConstIterator::reference const_reference;
    typedef typename Iterator::pointer pointer;
    typedef typename ConstIterator::pointer const_pointer;
  public:
    Expansion() : _indices(0), _coefficients() { }
    Expansion(unsigned int as) : _indices(as), _coefficients() { }
    Expansion(unsigned int as, unsigned int deg, double c0, ...);
    Expansion(unsigned int as, unsigned int nnz, Int a00, ...);
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    template<class XX> Expansion(const Expansion<MultiIndex,XX>& p);

    static Expansion<MultiIndex,X> variable(unsigned int as, unsigned int i);

    Void swap(Expansion<MultiIndex,X>& other) {
        std::swap(this->_indices,other._indices);
        std::swap(this->_coefficients,other._coefficients); }

    Bool operator==(const Expansion<MultiIndex,X>& other) const { return this->_coefficients == other._coefficients && this->_indices==other._indices; }
    Bool operator!=(const Expansion<MultiIndex,X>& other) const { return !(*this==other); }

    unsigned int argument_size() const { return this->_indices.element_size(); }
    unsigned int number_of_nonzeros() const { return _coefficients.size(); }
    unsigned int degree() const { if(this->empty()) { return 0u; } return this->_indices.back().degree(); }

    Bool empty() const { return this->_coefficients.empty(); }
    unsigned int size() const { return this->_coefficients.size(); }
    unsigned int capacity() const { return this->_coefficients.capacity(); }
    Void reserve(SizeType nnz) { this->_indices.reserve(nnz); this->_coefficients.reserve(nnz); }
    Void resize(SizeType nnz) { this->_indices.resize(nnz); this->_coefficients.resize(nnz); }

    Void insert(const MultiIndex& a, const RealType& c); // Non-inline user function
    Void append(const MultiIndex& a, const RealType& c) {
        this->_append(a,c); }
    Void prepend(const MultiIndex& a, const RealType& c) {
        this->_prepend(a,c); }
    Void append(const MultiIndex& a1, const MultiIndex& a2, const RealType& c) {
        this->_append(a1,a2,c); }

    RealType& operator[](const MultiIndex& a) {
        Iterator iter=this->lower_bound(a);
        if(iter==this->end() || iter->index()!=a) { iter=this->_insert(iter,a,static_cast<RealType>(0.0)); }
        return iter->coefficient(); }
    const RealType& operator[](const MultiIndex& a) const {
        ConstIterator iter=this->lower_bound(a);
        if(iter==this->end() || iter->index()!=a) { return _zero; }
        else { return iter->coefficient(); } }

    Iterator begin() { return Iterator(_indices.begin(),_coefficients.begin()); }
    Iterator end() { return Iterator(_indices.end(),_coefficients.end()); }
    Iterator lower_bound(const MultiIndex& a) {
        return std::lower_bound(this->begin(),this->end(),a,key_less<value_type,key_type>()); }
    Iterator find(const MultiIndex& a) {
        Iterator iter=this->lower_bound(a); if(iter!=this->end() && iter->index()!=a) { iter=this->end(); } return iter; }

    ConstIterator begin() const { return ConstIterator(_indices.begin(),_coefficients.begin()); }
    ConstIterator end() const { return ConstIterator(_indices.end(),_coefficients.end()); }
    ConstIterator lower_bound(const MultiIndex& a) const {
        return std::lower_bound(this->begin(),this->end(),a,key_less<value_type,key_type>()); }
    ConstIterator find(const MultiIndex& a) const {
        ConstIterator iter=this->lower_bound(a); if(iter!=this->end() && iter->index()!=a) { iter=this->end(); } return iter; }

    reference front() { return reference(_indices.begin(),_coefficients.begin()); }
    const_reference front() const { return const_reference(_indices.begin(),_coefficients.begin()); }
    reference back() { return reference(_indices.end()-1,_coefficients.end()-1); }
    const_reference back() const { return const_reference(_indices.end()-1,_coefficients.end()-1); }

    Void erase(Iterator iter) { iter->coefficient()=0.0; }
    Void clear() { _indices.clear(); _coefficients.clear(); }

    Void sort() {
        //std::cerr<<"sorting... "<<std::flush;
        std::sort(this->begin(),this->end());
    }
    Void remove_zeros() {
        //std::cerr<<"remove_zeros... "<<std::flush;
        Iterator new_end=std::remove_if(this->begin(),this->end(),data_is_zero<value_type>());
        this->resize(new_end-this->begin());
    }


    Void cleanup() { this->sort(); this->remove_zeros(); }

    Void check() const { }

    Expansion<MultiIndex,X> embed(unsigned int before_size, const Expansion<MultiIndex,X>&, unsigned int after_size) const;
  public:
    Iterator _insert(Iterator p, const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_insert "<<*this<<" "<<p._ptr()<<" "<<a<<" "<<x<<std::endl;
        if(this->size()+1u>this->capacity()) {
            difference_type i=p-begin();
            this->resize(this->size()+1);
            p=begin()+i;
        } else {
            this->resize(this->size()+1);
        }
        return this->_allocated_insert(p,a,x); }
    Iterator _insert(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_insert "<<*this<<" "<<a<<" "<<x<<std::endl;
        this->resize(this->size()+1u);
        Iterator p=std::lower_bound(this->begin(),this->end()-1,a,key_less<value_type,key_type>());
        return this->_allocated_insert(p,a,x); }
    Iterator _allocated_insert(Iterator p, const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_allocated_insert "<<*this<<" "<<p<<" "<<p-begin()<<" "<<a<<" "<<x<<std::endl;
        Iterator curr=this->end(); --curr; Iterator prev=curr;
        while(curr!=p) { --prev; *curr=*prev; curr=prev; }
        curr->index()=a; curr->coefficient()=x; return p; }
    Void _prepend(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_prepend "<<*this<<" "<<a<<" "<<x<<std::endl;
        this->resize(this->size()+1u);
        _allocated_insert(this->begin(),a,x); }
    Void _append(const MultiIndex& a, const RealType& x) {
        //std::cerr<<"_append "<<*this<<" "<<a<<" "<<x<<"... "<<std::flush;
        this->resize(this->size()+1u);
        this->back().key()=a; this->back().data()=x; }
    Void _append(const MultiIndex&  a1, const MultiIndex&  a2, const RealType& x) {
        //std::cerr<<"_append "<<*this<<" "<<a1<<" "<<a2<<" "<<x<<std::endl;
        this->resize(this->size()+1u);
        this->back().key()=a1; this->back().key()+=a2; this->back().data()=x; }
  public:
    OutputStream& write(OutputStream& os, const Array<StringType>& variables) const;
  private:
  public:
    MultiIndexList _indices;
    std::vector<scalar_type> _coefficients;
};


/*
template<class X>
OutputStream& operator<<(OutputStream& os, const typename Expansion<MultiIndex,X>::ConstIterator& citer) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<ConstIterator>";
}

template<class X>
OutputStream& operator<<(OutputStream& os, const typename Expansion<MultiIndex,X>::reference& ref) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<reference>";
}

OutputStream& operator<<(OutputStream& os, const Expansion<MultiIndex,FloatDP>::const_reference& cref) {
    //return os << "(" << v.get<0>() << ":" << v.get<1>() << ")";
    return os << "<const_reference>";
}
*/

template<class X>
Void
Expansion<MultiIndex,X>::insert(const MultiIndex& a, const X& x) {
    this->_insert(a,x);
}


template<class X> X Expansion<MultiIndex,X>::_zero = 0.0;

template<class X>
Expansion<MultiIndex,X>::Expansion(unsigned int as, unsigned int deg, double c0, ...)
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
Expansion<MultiIndex,X>::Expansion(unsigned int as, unsigned int nnz, Int a00, ...)
    : _indices(as), _coefficients()
{
    MultiIndex a(as);
    double x;
    va_list args;
    va_start(args,a00);
    for(unsigned int i=0; i!=nnz; ++i) {
        for(unsigned int j=0; j!=as; ++j) {
            if(i==0 && j==0) { a[j]=a00; }
            else { a[j]=va_arg(args,Int); }
        }
        x=va_arg(args,double);
        if(x!=0) { this->append(a,x); }
    }
    va_end(args);
    this->cleanup();
}

template<class X> template<class XX>
Expansion<MultiIndex,X>::Expansion(const std::map<MultiIndex,XX>& m)
{
    ARIADNE_ASSERT(!m.empty());
    this->_indices=MultiIndexList(m.begin()->first.size());
    for(typename std::map<MultiIndex,XX>::ConstIterator iter=m.begin(); iter!=m.end(); ++iter) {
        this->append(iter->first,X(iter->second));
    }
    this->remove_zeros();
}

template<class X> template<class XX> inline
Expansion<MultiIndex,X>::Expansion(const Expansion<MultiIndex,XX>& p)
    : _indices(p.argument_size()), _coefficients()
{
    for(typename Expansion<MultiIndex,XX>::ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->index(),X(iter->coefficient())); }
}


template<class X> Expansion<MultiIndex,X> Expansion<MultiIndex,X>::variable(unsigned int n, unsigned int i) {
    Expansion<MultiIndex,X> p(n); p[MultiIndex::zero(n)]=0.0; p[MultiIndex::unit(n,i)]=X(1);
    return p;
}


template<class X, class Y>
Y evaluate(const Expansion<MultiIndex,X>& x, const Vector<Y>& y)
{
    Y zero = y.zero_element(); zero*=0;
    Y one = zero; one+=1;

    Y r=zero;
    Y t=zero;
    for(typename Expansion<MultiIndex,X>::ConstIterator iter=x.begin();
        iter!=x.end(); ++iter)
    {
        UniformConstReference<MultiIndex> j=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        t=one;
        for(Nat k=0; k!=x.argument_size(); ++k) {
            for(Nat l=0; l!=j[k]; ++l) {
                t=t*y[k];
            }
        }
        t*=c;
        r+=t;
    }

    return r;
}



template<class X>
Expansion<MultiIndex,X> embed(unsigned int before_size, const Expansion<MultiIndex,X>& x, unsigned int after_size)
{
    Nat old_size=x.argument_size();
    Nat new_size=before_size+old_size+after_size;
    Expansion<MultiIndex,X> r(new_size);
    MultiIndex old_index(old_size);
    MultiIndex new_index(new_size);
    for(typename Expansion<MultiIndex,X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        old_index=iter->index();
        for(Nat j=0; j!=old_size; ++j) {
            Nat aj=old_index[j];
            new_index[j+before_size]=aj;
        }
        r.append(new_index,iter->coefficient());
    }
    return r;
}


inline Expansion<MultiIndex,FloatDP> midpoint(const Expansion<MultiIndex,ExactIntervalType>& pse) {
    Expansion<MultiIndex,FloatDP> r(pse.argument_size());
    for(Expansion<MultiIndex,ExactIntervalType>::ConstIterator iter=pse.begin(); iter!=pse.end(); ++iter) {
        //r.append(iter->index(),midpoint(iter->coefficient())); }
        r.append(iter->index(),midpoint(iter->coefficient())); }
    return r;
}


template<class X>
OutputStream& Expansion<MultiIndex,X>::write(OutputStream& os, const Array<StringType>& variable_names) const
{
    ARIADNE_ASSERT(this->argument_size()==variable_names.size());
    const Expansion<MultiIndex,X>& p=*this;
    if(p.size()==0) {
        os << "0";
    } else {
        Bool first_term=true;
        for(ConstIterator iter=p.begin(); iter!=p.end(); ++iter) {
            MultiIndex a=iter->index();
            X v=iter->coefficient();
            if(v!=0) {
                if(v>0 && !first_term) { os<<"+"; }
                first_term=false;
                Bool first_factor=true;
                if(v<0) { os<<"-"; }
                if(abs(v)!=1 || a.degree()==0) { os<<abs(v); first_factor=false; }
                for(Nat j=0; j!=a.size(); ++j) {
                    if(a[j]!=0) {
                        if(first_factor) { first_factor=false; } else { os <<"*"; }
                        os<<variable_names[j]; if(a[j]!=1) { os<<"^"<<Int(a[j]); }
                    }
                }
            }
        }
    }
    return os;
}


template<class X>
OutputStream& operator<<(OutputStream& os, const Expansion<MultiIndex,X>& p) {
    Array<StringType> variable_names(p.argument_size());
    for(Nat j=0; j!=p.argument_size(); ++j) {
        StringStream sstr;
        sstr << 'x' << j;
        variable_names[j]=sstr.str();
    }
    return p.write(os,variable_names);
}




template<class X, class Y>
Vector<Y> evaluate(const Vector< Expansion<MultiIndex,X> >& x, const Vector<Y>& y)
{
    Vector<Y> r(x.size());
    for(unsigned int i=0; i!=x.size(); ++i) {
        r[i]=evaluate(x[i],y);
    }
    return r;
}


template<class X> Vector< Expansion<MultiIndex,X> > operator*(const Expansion<MultiIndex,X>& e, const Vector<FloatDP> v) {
    Vector< Expansion<MultiIndex,X> > r(v.size(),Expansion<MultiIndex,X>(e.argument_size()));
    for(Nat i=0; i!=r.size(); ++i) {
        ARIADNE_ASSERT(v[i]==0.0 || v[i]==1.0);
        if(v[i]==1.0) { r[i]=e; }
    }
    return r;
}


inline Vector< Expansion<MultiIndex,FloatDP> > midpoint(const Vector< Expansion<MultiIndex,ExactIntervalType> >& pse) {
    Vector< Expansion<MultiIndex,FloatDP> > r(pse.size());
    for(Nat i=0; i!=pse.size(); ++i) {
        r[i]=midpoint(pse[i]); }
    return r;
}








}

#endif /* ARIADNE_EXPANSION_HPP */
