/***************************************************************************
 *            expansion.h
 *
 *  Copyright 2008-15  Pieter Collins
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

#include "algebra/vector.h"
#include "algebra/multi_index.h"



namespace Ariadne {

typedef MultiIndex::WordType WordType;

template<class X> class Expansion;
typedef Expansion<Float64> FloatExpansion;



template<class X>
struct ExpansionValue {
    typedef MultiIndex::SizeType SizeType;
    typedef MultiIndex::WordType WordType;

    ~ExpansionValue() { delete[] _p; }
    ExpansionValue(const MultiIndex& a, const X& x)
        : _n(a.size()), _nw(a.word_size()), _p() { _p=new WordType[_nw+_ds]; key()=a; data()=x;  }
    ExpansionValue(const ExpansionValue<X>& v)
        : _n(v._n), _nw(v._nw), _p() { _p=new WordType[v._nw+_ds]; _assign(v._p); }
    ExpansionValue<X>& operator=(const ExpansionValue<X>& v) {
        if(this!=&v) { _resize(v._n,v._nw); _assign(v._p); } return *this; }
    const MultiIndex& key() const { return *reinterpret_cast<const MultiIndex*>(this); }
    const X& data() const { return *reinterpret_cast<const X*>(_p+_nw); }
    MultiIndex& key() { return *reinterpret_cast<MultiIndex*>(this); }
    X& data() { return *reinterpret_cast<X*>(_p+_nw); }
  private:
    Void _resize(SizeType n, SizeType nw) {
        if(_nw!=nw) { delete[] _p; _nw=nw; _p=new WordType[_nw+_ds]; }_n=n; }
    Void _assign(const WordType* p) {
        for(SizeType j=0; j!=_nw+_ds; ++j) { _p[j]=p[j]; } }
  public:
    SizeType _n; SizeType _nw; WordType* _p;
    static const SizeType _ds=sizeof(X)/sizeof(WordType);
};


struct DataLess {
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) { return v1.data()<v2.data(); }
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const X& d2) { return v1.data()<d2; }
};

struct DataIsZero {
    template<class X> Bool operator()(const ExpansionValue<X>& v) { return decide(v.data()==static_cast<X>(0)); }
};

struct KeyEquals {
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) { return v1.key()==v2.key(); }
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) { return v1.key()==k2; }
    Bool operator()(const MultiIndex& k1, const MultiIndex& k2) { return k2==k2; }
};

struct GradedKeyLess {
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) {
        return graded_less(v1.key(),v2.key()); }
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) {
        return graded_less(v1.key(),k2); }
    Bool operator()(const MultiIndex& k1, const MultiIndex& k2) {
        return graded_less(k1,k2); }
};

struct LexicographicKeyLess {
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) {
        return lexicographic_less(v1.key(),v2.key()); }
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) {
        return lexicographic_less(v1.key(),k2);; }
    Bool operator()(const MultiIndex& k1, const MultiIndex& k2) {
        return lexicographic_less(k1,k2); }
};

struct ReverseLexicographicKeyLess {
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const ExpansionValue<X>& v2) {
        return reverse_lexicographic_less(v1.key(),v2.key()); }
    template<class X> Bool operator()(const ExpansionValue<X>& v1, const MultiIndex& k2) {
        return reverse_lexicographic_less(v1.key(),k2);; }
    Bool operator()(const MultiIndex& k1, const MultiIndex& k2) {
        return reverse_lexicographic_less(k1,k2); }
};



template<class X>
inline OutputStream& operator<<(OutputStream& os, const ExpansionValue<X>& dv) {
    return os << "["<<dv._n<<","<<dv._nw<<","<<(Void*)dv._p<<"]"<<dv.key()<<":"<<dv.data();
    //return os << dv.key() << ":" << dv.data();
}


template<class X, class Ref, class Ptr> class ExpansionIterator;
template<class X, class Ref, class Ptr> OutputStream& operator<<(OutputStream&, const ExpansionIterator<X,Ref,Ptr>&);

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
    friend OutputStream& operator<< <>(OutputStream&,const ExpansionIterator<X,Ref,Ptr>&);
    typedef ExpansionIterator<X,Ref,Ptr> Iter;
  public:
    typedef MultiIndex::SizeType SizeType;
    typedef MultiIndex::IndexType ByteType;
    typedef MultiIndex::WordType WordType;
    typedef X DataType;

    typedef int difference_type;
  public:
    ExpansionIterator()
        : _n(), _nw(), _p() { }
    ExpansionIterator(SizeType n, SizeType nw, WordType* p)
        : _n(n), _nw(nw), _p(p) { }
    template<class Ref2,class Ptr2> ExpansionIterator(const ExpansionIterator<X,Ref2,Ptr2>& i)
        : _n(i._n), _nw(i._nw), _p(i._p) { }
    template<class Ref2,class Ptr2> Bool equal(const ExpansionIterator<X,Ref2,Ptr2>& i) const {
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
    const WordType* _word_ptr() { return this->_p; }
    const Void* _ptr() { return this->_p; }
  private:
    SizeType _n;
    SizeType _nw;
    WordType* _p;
    static const SizeType _ds=sizeof(DataType)/sizeof(WordType);
};


template<class X, class Ref, class Ptr>
inline OutputStream& operator<<(OutputStream& os, const ExpansionIterator<X,Ref,Ptr>& piter) {
    const MultiIndex::ByteType* ap=reinterpret_cast<const MultiIndex::ByteType*>(piter._p);
    const X* xp=reinterpret_cast<const X*>(piter._p+piter._nw);
    os << "(n="<<piter._n<<", nw="<<piter._nw<<", p="<<(Void*)piter._p<<", a=";
    for(MultiIndex::SizeType i=0; i!=piter._n; ++i) { os<<(i==0?"(":",")<<Int(ap[i]); }
    os << "), x="<<*xp<<")";
    return os;
}

template<class X>
class Expansion
{
    static const Ariadne::SizeType DEFAULT_CAPACITY=16;
  protected:
    static X _zero;
  public:
    typedef X RealType;
    typedef X CoefficientType;
    typedef unsigned short int SmoothnessType;
    typedef MultiIndex::SizeType SizeType;
    typedef MultiIndex::ByteType ByteType;
    typedef MultiIndex::WordType WordType;
    typedef Int DifferenceType;
    typedef MultiIndex KeyType;
    typedef X DataType;
  public:
    typedef ExpansionValue<X> value_type;
    typedef ExpansionValue<X>& reference;
    typedef ExpansionValue<X>const& const_reference;
    typedef ExpansionValue<X>* pointer;
    typedef ExpansionValue<X>const* const_pointer;
    typedef ExpansionIterator<X,ExpansionValue<X>&, ExpansionValue<X>* > iterator;
    typedef ExpansionIterator<X,ExpansionValue<X>const&, ExpansionValue<X>const* > const_iterator;
  public:
    typedef ExpansionValue<X> ValueType;
    typedef ExpansionValue<X>& Reference;
    typedef ExpansionValue<X>const& ConstReference;
    typedef ExpansionIterator<X,ExpansionValue<X>&, ExpansionValue<X>* > Iterator;
    typedef ExpansionIterator<X,ExpansionValue<X>const&, ExpansionValue<X>const* > ConstIterator;

  public:
    ~Expansion();
    explicit Expansion(); // DEPRECTATED
    explicit Expansion(SizeType as);
    Expansion(SizeType as, DegreeType deg, InitializerList<X> lst); // DEPRECTATED
    Expansion(SizeType as, InitializerList<PairType<InitializerList<Int>,X>> lst);
    Expansion(InitializerList<PairType<InitializerList<Int>,X>> lst);

    template<class XX, EnableIf<IsConvertible<XX,X>> =dummy>
        Expansion(const Expansion<XX>& p);
    template<class XX, EnableIf<IsConstructible<X,XX>> =dummy, DisableIf<IsConvertible<XX,X>> =dummy>
        explicit Expansion(const Expansion<XX>& p);

    // DEPRECTATED Conversions to raw type
    Expansion<RawFloat64>& raw();
    Expansion<RawFloat64>const& raw() const;

    Void swap(Expansion<X>& other);

    template<class XX> friend Bool same(const Expansion<XX>& self, const Expansion<XX>& other);
    Bool operator==(const Expansion<X>& other) const;
    Bool operator!=(const Expansion<X>& other) const;

    SizeType argument_size() const;
    DegreeType degree() const;

    Bool empty() const;
    SizeType number_of_nonzeros() const; // Deprecated
    SizeType number_of_terms() const;
    SizeType size() const; // Number of terms
    SizeType capacity() const;
    Void reserve(SizeType nnz);
    Void resize(SizeType nnz);

    Void append(Pair<MultiIndex,CoefficientType> const& ac);
    Void append(const MultiIndex& a, const CoefficientType& c);
    Void prepend(const MultiIndex& a, const CoefficientType& c);
    Void append_sum(const MultiIndex& a1, const MultiIndex& a2, const CoefficientType& c);

    const CoefficientType& operator[](const MultiIndex& a) const;
//    ExpansionValueReference<X> operator[](const MultiIndex& a);

    Iterator begin();
    ConstIterator begin() const;
    Iterator end();
    ConstIterator end() const;
    Iterator find(const MultiIndex& a);
    ConstIterator find(const MultiIndex& a) const;

    Iterator insert(Iterator pos, const MultiIndex& a, const X& x);
    Void erase(Iterator iter);

    Void clear();
    Void check() const;
  public: // Deprecated
    Void remove_zeros();
    Void combine_terms();
  public:
    Void graded_sort();
    Void lexicographic_sort();
    Void reverse_lexicographic_sort();
  public:
    SizeType _vector_size() const;
    SizeType _index_size() const;
    SizeType _element_size() const;
  protected:
    static const unsigned int sizeof_byte=sizeof(ByteType);
    static const unsigned int sizeof_word=sizeof(WordType);
    static const unsigned int sizeof_data=sizeof(DataType);

    const WordType* _begin_ptr() const;
    const WordType* _end_ptr() const;
    WordType* _begin_ptr();
    WordType* _end_ptr();
    Iterator _insert(Iterator p, const MultiIndex& a, const CoefficientType& x);
    Iterator _allocated_insert(Iterator p, const MultiIndex& a, const CoefficientType& x);
    Void _prepend(const MultiIndex& a, const CoefficientType& x);
    Void _append(const MultiIndex& a, const CoefficientType& x);
    Void _append(const MultiIndex&  a1, const MultiIndex&  a2, const CoefficientType& x);
  protected:
    template<class CMP> Void _sort(const CMP& cmp);
    template<class CMP> CoefficientType& _at(const MultiIndex& a,const CMP& cmp);
    template<class CMP> Iterator _insert(const MultiIndex& a, const CoefficientType& x,const CMP& cmp);
  private:
    template<class T> friend Expansion<T> embed(SizeType, Expansion<T> const&, SizeType);
    Expansion<X> _embed(SizeType before_size, SizeType after_size) const;
    Bool _same(Expansion<X> const& other) const;
  public:
    OutputStream& write(OutputStream& os, const Array<StringType>& variables) const;
    OutputStream& write_polynomial(OutputStream& os) const;
    OutputStream& write_map(OutputStream& os) const;
    OutputStream& write(OutputStream& os) const;
  private:
    SizeType _argument_size;
    std::vector<WordType> _coefficients;
};

template<class X> Expansion<MidpointType<X>> midpoint(const Expansion<X>& pse);
template<class X> Expansion<X> embed(SizeType, Expansion<X> const&, SizeType);
template<class X> OutputStream& operator<<(OutputStream& os, const Expansion<X>& p);

template<class X, class CMP> class SortedExpansionValueReference;

template<class X, class CMP> class SortedExpansion
    : public Expansion<X>
{
  public:
    typedef X CoefficientType;
    typedef typename Expansion<X>::Iterator Iterator;
    typedef typename Expansion<X>::ConstIterator ConstIterator;
  public:
    using Expansion<X>::Expansion;
    SortedExpansion(Expansion<X> e);
    Void sort();

    X const& operator[](const MultiIndex& a) const;
    SortedExpansionValueReference<X,CMP> operator[](const MultiIndex& a);

    Void insert(const MultiIndex& a, const CoefficientType& c);
    Void set(const MultiIndex& a, const CoefficientType& c);
    CoefficientType& at(const MultiIndex& a);
    CoefficientType const& get(const MultiIndex& a) const;
    Iterator find(const MultiIndex& a);
    ConstIterator find(const MultiIndex& a) const;
};

template<class X, class CMP> class SortedExpansionValueReference {
    SortedExpansion<X,CMP>& _e; MultiIndex const& _a;
  public:
    SortedExpansionValueReference(SortedExpansion<X,CMP>& e, const MultiIndex& a) : _e(e), _a(a) { }
    operator const X& () const { return _e.get(_a); }
    //operator X& () { return _e.at(_a); }
    SortedExpansionValueReference<X,CMP>& operator=(const X& x) { _e.set(_a,x); return *this; }
    SortedExpansionValueReference<X,CMP>& operator+=(X const& c) { _e.at(_a)+=c; return *this; }
    SortedExpansionValueReference<X,CMP>& operator-=(X const& c) { _e.at(_a)-=c; return *this; }
    SortedExpansionValueReference<X,CMP>& operator*=(X const& c) { _e.at(_a)*=c; return *this; }
    SortedExpansionValueReference<X,CMP>& operator/=(X const& c) { _e.at(_a)/=c; return *this; }
};

template<class X, class CMP> inline const X& SortedExpansion<X,CMP>::operator[](const MultiIndex& a) const {
    return const_cast<SortedExpansion<X,CMP>&>(*this).at(a); }

template<class X, class CMP> inline SortedExpansionValueReference<X,CMP> SortedExpansion<X,CMP>::operator[](const MultiIndex& a) {
    return SortedExpansionValueReference<X,CMP>(*this,a); }

template<class X> inline Bool same(const Expansion<X>& self, const Expansion<X>& other) {
    return self._same(other);
}

template<class X> inline Expansion<X> embed(SizeType before_size, Expansion<X> const& x, SizeType after_size) {
    return x._embed(before_size,after_size);
}

template<class X> inline OutputStream& operator<<(OutputStream& os, Expansion<X> const& p) {
    return p.write(os);
}

template<class X> template<class XX, EnableIf<IsConvertible<XX,X>>> inline
Expansion<X>::Expansion(const Expansion<XX>& p)
    : _argument_size(p.argument_size())
{
    for(auto iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->key(),iter->data());
    }
}

template<class X> template<class XX, EnableIf<IsConstructible<X,XX>>, DisableIf<IsConvertible<XX,X>>> inline
Expansion<X>::Expansion(const Expansion<XX>& e)
    : _argument_size(e.argument_size())
{
    for(auto iter=e.begin(); iter!=e.end(); ++iter) {
        this->append(iter->key(),X(iter->data()));
    }
}

template<class X> Expansion<MidpointType<X>> midpoint(const Expansion<X>& e) {
    Expansion<MidpointType<X>> r(e.argument_size());
    for(auto iter=e.begin(); iter!=e.end(); ++iter) {
        r.append(iter->key(),midpoint(iter->data())); }
    return r;
}



} // namespace Ariadne

#include "evaluate.h"

#endif /* ARIADNE_EXPANSION_H */
