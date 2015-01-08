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

typedef MultiIndex::word_type WordType;

template<class X> class Expansion;
typedef Expansion<Float> FloatExpansion;


template<class X> Expansion<X> embed(unsigned int, const Expansion<X>&, unsigned int);



template<class FwdIter, class Op> FwdIter unique_key(FwdIter first, FwdIter last, Op op);



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
    template<class X> bool operator()(const ExpansionValue<X>& v) { return v.data()==static_cast<X>(0); }
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
    typedef ExpansionIterator<X,ExpansionValue<X>&, ExpansionValue<X>* > Iterator;
    typedef ExpansionIterator<X,ExpansionValue<X>const&, ExpansionValue<X>const* > ConstIterator;
  public:
    explicit Expansion();
    explicit Expansion(unsigned int as);
    Expansion(unsigned int as, unsigned int deg, std::initializer_list<X> lst);
    Expansion(std::initializer_list< std::pair<std::initializer_list<int>,X> > lst);
    Expansion(unsigned int as, std::initializer_list< std::pair<std::initializer_list<int>,X> > lst);
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    template<class XX, typename std::enable_if<std::is_convertible<XX,X>::value,int>::type=0>
        Expansion(const Expansion<XX>& p);
    template<class XX, typename std::enable_if<std::is_constructible<X,XX>::value and not std::is_convertible<XX,X>::value,int>::type=0>
        explicit Expansion(const Expansion<XX>& p);

    Expansion<RawFloat>& raw();
    Expansion<RawFloat>const& raw() const;

    Expansion<X>& operator=(const X& x);

    static Expansion<X> variable(unsigned int as, unsigned int i);

    void swap(Expansion<X>& other);

    bool operator==(const Expansion<X>& other) const;
    bool operator!=(const Expansion<X>& other) const;

    unsigned int argument_size() const;
    unsigned int number_of_nonzeros() const;
    DegreeType degree() const;
    const std::vector<word_type>& coefficients() const;

    bool empty() const;
    unsigned int size() const;
    void reserve(size_type nnz);
    void resize(size_type nnz);

    template<class CMP> void insert(const MultiIndex& a, const RealType& c, const CMP& cmp);
    template<class CMP> void set(const MultiIndex& a, const RealType& c, const CMP& cmp);
    template<class CMP> RealType& at(const MultiIndex& a, const CMP& cmp);
    template<class CMP> iterator find(const MultiIndex& a, const CMP& cmp);
    template<class CMP> const_iterator find(const MultiIndex& a, const CMP& cmp) const;

    void append(const MultiIndex& a, const RealType& c);
    void prepend(const MultiIndex& a, const RealType& c);
    void append(const MultiIndex& a1, const MultiIndex& a2, const RealType& c);

    const RealType& operator[](const MultiIndex& a) const;

    iterator begin();
    iterator end();
    iterator find(const MultiIndex& a);

    const_iterator begin() const;
    const_iterator end() const;
    const_iterator find(const MultiIndex& a) const;

    reference front();
    reference back();
    const_reference front() const;
    const_reference back() const;

    void erase(iterator iter);
    void clear();

    void remove_zeros();
    void combine_terms();

    template<class CMP> void sort(const CMP& cmp);

    void graded_sort();
    void lexicographic_sort();
    void reverse_lexicographic_sort();

    void check() const;

  public:
    size_type _vector_size() const;
    size_type _index_size() const;
    size_type _element_size() const;
  public:
    const word_type* _begin_ptr() const;
    const word_type* _end_ptr() const;
    word_type* _begin_ptr();
    word_type* _end_ptr();
    template<class CMP> RealType& _at(const MultiIndex& a,const CMP& cmp);
    template<class CMP> iterator _insert(const MultiIndex& a, const RealType& x,const CMP& cmp);
    iterator _insert(iterator p, const MultiIndex& a, const RealType& x);
    iterator _allocated_insert(iterator p, const MultiIndex& a, const RealType& x);
    void _prepend(const MultiIndex& a, const RealType& x);
    void _append(const MultiIndex& a, const RealType& x);
    void _append(const MultiIndex&  a1, const MultiIndex&  a2, const RealType& x);
    Expansion<X> _embed(unsigned int before_size, unsigned int after_size) const;
  public:
    std::ostream& write(std::ostream& os, const Array<std::string>& variables) const;
    std::ostream& write(std::ostream& os) const;
  private:
    size_type _argument_size;
    std::vector<word_type> _coefficients;
};

template<class X, class CMP> class SortedExpansion
    : public Expansion<X>
{
  public:
    typedef typename Expansion<X>::RealType RealType;
    using Expansion<X>::Expansion;

    void insert(const MultiIndex& a, const RealType& c) { this->insert(a,c,CMP()); }
    void set(const MultiIndex& a, const RealType& c) { this->set(a,c,CMP()); }
    RealType& at(const MultiIndex& a) { return this->at(a,CMP()); }
    typename Expansion<X>::Iterator find(const MultiIndex& a) { return this->find(a,CMP()); }
    typename Expansion<X>::ConstIterator find(const MultiIndex& a) const { return this->find(a,CMP()); }
    void sort() { this->sort(CMP()); }
};

// Disable construction of Expansion<Rational> since above implementation only
// works for "plain old data" types
#if defined HAVE_GMPXX_H and defined ARIADNE_NUMERIC_H
template<> class Expansion<Rational>;
#endif



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
//! where \f$j=\max\{i\mid \alpha_i\neq\beta_i\}\f$.
//! For a given term \f$c_\alpha x^\alpha\f$, let \f$k=\max\{j\mid \alpha_j\neq\beta_j\}\f$, where \f$\beta\f$ is the next multi-index.
//! We update register \f$r_k\f$ by \f[r'_k=(((c_\alpha + r_1) x^{\alpha_1} + r_2 )x^{\alpha_2}+\cdots r_k)x^{\alpha_k-\beta_k}.\f]
//! The result is obtained by updating a fictional register \f$r_{n+1}\f$ at the last step.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
template<class X, class Y> Y horner_evaluate(const Expansion<X>& e, const Vector<Y>& x);

template<class X, class Y> Y power_evaluate(const Expansion<X>& e, const Vector<Y>& y);

template<class X, class Y> Y evaluate(const Expansion<X>& e, const Vector<Y>& y);

template<class X, class Y> Y simple_evaluate(const Expansion<X>& e, const Vector<Y>& y);

template<class T> Expansion<MidpointType<T>> midpoint(const Expansion<T>& pse);

template<class X, class Y> Vector<Y> evaluate(const Vector< Expansion<X> >& x, const Vector<Y>& y);

template<class X> Vector< Expansion<X> > operator*(const Expansion<X>& e, const Vector<Float> v);

template<class T> Vector< Expansion<MidpointType<T>> > midpoint(const Vector< Expansion<T> >& pse);



template<class X> inline Expansion<X> embed(unsigned int before_size, const Expansion<X>& x, unsigned int after_size) {
    return x._embed(before_size,after_size);
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Expansion<X>& p) {
    return p.write(os);
}



template<class X> template<class XX, typename std::enable_if<std::is_convertible<XX,X>::value,int>::type> inline
Expansion<X>::Expansion(const Expansion<XX>& p)
    : _argument_size(p.argument_size())
{
    for(typename Expansion<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->key(),X(iter->data())); }
}

template<class X>
template<class XX, typename std::enable_if<std::is_constructible<X,XX>::value and not std::is_convertible<XX,X>::value,int>::type>
inline
Expansion<X>::Expansion(const Expansion<XX>& p)
    : _argument_size(p.argument_size())
{
    for(typename Expansion<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->key(),X(iter->data())); }
}



template<class X> template<class CMP> typename Expansion<X>::RealType& Expansion<X>::_at(const MultiIndex& a,const CMP& cmp) {
    iterator p=std::lower_bound(this->begin(),this->end(),a,cmp);
    if(p!=this->end() && p->key()==a) { return p->data(); }
    else { p=this->_insert(p,a,X(0)); return p->data(); }
}

template<class X> template<class CMP> typename Expansion<X>::iterator Expansion<X>::_insert(const MultiIndex& a, const RealType& x,const CMP& cmp) {
    //std::cerr<<"_insert "<<*this<<" "<<a<<" "<<x<<std::endl;
    _coefficients.resize(_coefficients.size()+_element_size());
    iterator p=std::lower_bound(this->begin(),this->end()-1,a,cmp);
    return _allocated_insert(p,a,x);
}

template<class X> template<class CMP> void Expansion<X>::insert(const MultiIndex& a, const RealType& c, const CMP& cmp) {
    this->_insert(a,c,cmp);
}
template<class X> template<class CMP> void Expansion<X>::set(const MultiIndex& a, const RealType& c, const CMP& cmp) {
    this->_at(a,cmp)=c;
}
template<class X> template<class CMP> X& Expansion<X>::at(const MultiIndex& a, const CMP& cmp) {
    return this->_at(a,cmp);
}
template<class X> template<class CMP> typename Expansion<X>::iterator Expansion<X>::find(const MultiIndex& a, const CMP& cmp) {
    iterator iter=std::lower_bound(this->begin(),this->end(),a,cmp); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter;
}
template<class X> template<class CMP> typename Expansion<X>::const_iterator Expansion<X>::find(const MultiIndex& a, const CMP& cmp) const {
    iterator iter=std::lower_bound(this->begin(),this->end(),a,cmp); if(iter!=this->end() && iter->key()!=a) { iter=this->end(); } return iter;
}


} // namespace Ariadne

#include "evaluate.tcc"

#endif /* ARIADNE_EXPANSION_H */
