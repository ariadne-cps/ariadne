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

#ifndef ARIADNE_EXPANSION_H
#define ARIADNE_EXPANSION_H

#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include "multi_index.h"


#define ARIADNE_USE_ARRAY_EXPANSION

namespace Ariadne {

#if defined DOXYGEN or not defined ARIADNE_USE_ARRAY_EXPANSION

//! \brief A power series expansion with coefficients of some type \a X.
//! Used in Polynomial, Differential and TaylorVariable classes.
template<class X>
class Expansion
{
    typedef X Real;
  public:
    //! \brief The type used to represent sizes.
    typedef unsigned int size_type;
    //! \brief The type used to index the coefficients.
    typedef MultiIndex key_type;
    //! \brief The type of coefficient.
    typedef X data_type;
    //! \brief The type of an (index,coefficient) pair.
    typedef typename std::map<MultiIndex,X>::value_type value_type;
    //! \brief The type of an reference to an (index,coefficient) pair.
    typedef typename std::map<MultiIndex,X>::reference reference;
    //! \brief The type of an reference to an (index,coefficient) pair.
    typedef typename std::map<MultiIndex,X>::const_reference const_reference;
    //! \brief The type of an iterator through (index,coefficient) pairs.
    typedef typename std::map<MultiIndex,X>::iterator iterator;
    //! \brief The type of a constant iterator through (index,coefficient) pairs.
    typedef typename std::map<MultiIndex,X>::const_iterator const_iterator;

    //typedef typename std::map<MultiIndex,X>::pointer pointer;
    //typedef typename std::map<MultiIndex,X>::const_pointer const_pointer;
  private:
    static X _zero;
  public:
    //! \brief A polynomial in no arguments; can only take a constant value.
    Expansion() : _argument_size() { }
    //! \brief The zero polynomial in \a as arguments.
    Expansion(unsigned int as) : _argument_size(as) { }
    //! \brief A dense polynomial with coefficients given by a list of doubles.
    Expansion(unsigned int as, unsigned int deg, double c0, ...);
    //! \brief A sparse polynomial with coefficients given by a list of indices and coefficients.
    Expansion(unsigned int as, unsigned int nnz, int a00, ...);
    //! \brief Construct from a map of (intex,coefficient) pairs.
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    //! \brief Construct from a polynomial with possibly different coefficient types.
    template<class XX> Expansion(const Expansion<XX>& p);

    //! \brief The polynomial \f$x_i\f$ in \a as variables.
    static Expansion<X> variable(unsigned int as, unsigned int i);
    //! \brief A vector of polynomials of the form \f$(x_0,\ldots,x_{n-1})\f$.
    static Vector< Expansion<X> > variables(unsigned int n);

    //! \brief Swap with another polynomial. Included for efficiency.
    void swap(Expansion<X>& p) {
        std::swap(this->_argument_size,p._argument_size);
        std::swap(this->_coefficients,p._coefficients); }

    //! \brief Equality operator.
    bool operator==(const Expansion<X>& p) const { return this->_coefficients == p._coefficients; }
    //! \brief Inequality operator.
    bool operator!=(const Expansion<X>& p) const { return this->_coefficients != p._coefficients; }

    //! \brief The number of argument variables of the polynomial.
    unsigned int argument_size() const { return this->_argument_size; }
    //! \brief The number of nonzero terms.
    unsigned int number_of_nonzeros() const { return this->_coefficients.size(); }
    //! \brief The degree of the polynomial.
    unsigned int degree() const { return (--this->_coefficients.end())->first.degree(); }

    //! \brief Reserve space for \a nnz nonzero elements.
    void reserve(size_type nnz) { }

    //! \brief Insert the term \f$cx^a\f$ into the polynomial.
    //! If polynomial already has a term in \f$x^a\f$, the coefficients are added.
    void insert(const MultiIndex& a, const Real& c) {
        assert(a.size()==this->_argument_size); this->_coefficients[a]+=c; }
    //! \brief Insert the term \f$cx^a\f$ into the polynomial.
    //! The index \a a must be higher than the index of all currently existing terms.
    void append(const MultiIndex& a, const Real& c) {
        this->insert(a,c); }
    //! \brief Insert the term \f$cx^{a_1+a_2}\f$ into the polynomial.
    //! The index \f$a_1+a_2\f$ must be higher than the index of all currently existing terms.
    void append(const MultiIndex& a1, const MultiIndex& a2, const Real& c) {
        this->insert(a1+a2,c); }

    //! \brief A reference to the coefficient of \f$x^a\f$.
    Real& operator[](const MultiIndex& a) {
        return this->_coefficients[a]; }
    //! \brief A constant reference to the coefficient of \f$x^a\f$.
    const Real& operator[](const MultiIndex& a) const {
        const_iterator iter=this->find(a); 
        if(iter==this->end()) { return _zero; }
        else { return iter->second; } }

    //! \brief An iterator to the first nonzero term. 
    iterator begin() { return this->_coefficients.begin(); }
    //! \brief An iterator to the end of the terms.
    iterator end() { return this->_coefficients.end(); }
    //! \brief An iterator to the term in \f$x^a\f$, or the end if there is no term in \f$x^a\f$.
    iterator find(const MultiIndex& a) { return this->_coefficients.find(a); }

    //! \brief A constant iterator to the first nonzero term.
    const_iterator begin() const { return this->_coefficients.begin(); }
    //! \brief A constant iterator to the end of the terms.
    const_iterator end() const { return this->_coefficients.end(); }
    //! \brief A constant iterator to the term in \f$x^a\f$, or the end if there is no term in \f$x^a\f$.
    const_iterator find(const MultiIndex& a) const { return this->_coefficients.find(a); }

    //! \brief Set the coefficient of \f$x^a\f$ to zero.
    void erase(iterator iter) { this->_coefficients.erase(iter); }
    //! \brief Remove all nonzero terms.
    void clear() { this->_coefficients.clear(); }

    //! \brief Clean up the representation of the polynomial, including
    //! sorting the terms and removing zero elements if necessary.
    void cleanup() { }
    //! \brief Check the representation of the polynomial for consistency.
    void check() {
        for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
            ARIADNE_ASSERT(iter->first.size()==this->_argument_size); } }

    //! \brief Evaluate the polynomial on a vector \a x over a ring.
    template<class Y> Y evaluate(const Vector<Y>& x) const;
    //! \brief Embed in a space of higher dimension.
    Expansion<X> embed(unsigned int new_size, unsigned int start) const;
  private:
    size_type _argument_size;
    std::map<MultiIndex,X> _coefficients;
};

#else

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
class Reference<const X> {
  public:
    Reference<const X>(const X& x) : _x(&x) { }
    operator X const& () const { return *_x; }
  private:
    const X* _x;
};

template<>
class Reference<MultiIndex> 
{
    friend class Reference<const MultiIndex>;
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type value_type;
    typedef MultiIndex::word_type word_type;
    typedef MultiIndex::reference reference;
    typedef MultiIndex::const_reference const_reference;
  public:
    Reference(size_type n, value_type* p) : _n(n), _p(p) { }
    Reference(MultiIndex& a) : _n(a.size()), _p(a.begin()) { }
    Reference& operator=(const MultiIndex& a) {
        assert(_n==a.size()); this->_assign(a.begin()); return *this; }
    Reference& operator=(const Reference<MultiIndex>& a) {
        // This is a substitute for the default assignment operator, so must use a const& argument.
        assert(_n==a.size()); this->_assign(a._p); return *this; }
    operator MultiIndex& () { return reinterpret_cast<MultiIndex&>(*this); }
    operator MultiIndex const& () const { return reinterpret_cast<MultiIndex const&>(*this); }
    size_type size() const { return this->_n; }
    value_type degree() const { return this->_p[this->_n]; }
    const_reference operator[](size_type i) const { return this->_p[i]; }
  private:
    void _assign(const value_type* p) {
        for(size_type j=0; j!=MultiIndex::_word_size(this->_n); ++j) {
            reinterpret_cast<word_type*>(this->_p)[j]=reinterpret_cast<const word_type*>(p)[j]; } }
  private:
    size_type _n; value_type* _p;
};

template<>
class Reference<const MultiIndex> {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type value_type;
    typedef MultiIndex::word_type word_type;
    typedef MultiIndex::const_reference const_reference;
  public:
    Reference(size_type n, const value_type* p) : _n(n), _p(const_cast<value_type*>(p)) { }
    Reference(const MultiIndex& a) : _n(a.size()), _p(const_cast<value_type*>(a.begin())) { }
    Reference(Reference<MultiIndex> a) : _n(a._n), _p(a._p) { }
    operator MultiIndex const& () const { return reinterpret_cast<MultiIndex const&>(*this); }
    size_type size() const { return this->_n; }
    value_type degree() const { return this->_p[this->_n]; }
    const_reference operator[](size_type i) const { return this->_p[i]; }
  private:
    size_type _n; value_type* _p;
};



template<class X>
struct ExpansionValue {
    ExpansionValue(const MultiIndex& a, const X& x) : first(a), second(x) { }
    MultiIndex first;
    X second;
};

template<class X>
inline bool operator<(const ExpansionValue<X>& dv1, const ExpansionValue<X>& dv2) {
    return dv1.first < dv2.first;
}

template<class X>
inline std::ostream& operator<<(std::ostream& os, const ExpansionValue<X>& dv) {
    return os << dv.first << ":" << dv.second;
}

template<class X>
struct ExpansionReference {
    ExpansionReference<X>& operator=(const ExpansionReference<X>& dr) {
        first=dr.first; second=dr.second; return *this; }
    ExpansionReference<X>& operator=(const ExpansionValue<X>& dv) {
        first=dv.first; second=dv.second; return *this; }
    operator ExpansionValue<X>() const {
        return ExpansionValue<X>(this->first,this->second); }
    Reference<MultiIndex> first;
    Reference<X> second;
};

template<class X>
struct ExpansionReference<const X> {
    ExpansionReference(const ExpansionReference<X>& r)
        : first(r.first), second(r.second) { }
    operator ExpansionValue<X>() const {
        return ExpansionValue<X>(this->first,this->second); }
    Reference<const MultiIndex> first;
    Reference<const X> second;
};

template<class X>
inline bool operator<(const ExpansionReference<X>& dr1, const ExpansionReference<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline bool operator<(const ExpansionReference<X>& dr1, const ExpansionValue<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline bool operator<(const ExpansionValue<X>& dr1, const ExpansionReference<X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline bool operator<(const ExpansionReference<const X>& dr1, const ExpansionReference<const X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline std::ostream& operator<<(std::ostream& os, ExpansionReference<X> dr) {
    return os << dr.first << ":" << dr.second;
}


template<class X, class Ref>
class ExpansionIterator
     : public boost::iterator_facade<ExpansionIterator<X,Ref>,
                                     ExpansionValue<X>,
                                     boost::random_access_traversal_tag,
                                     Ref >

{
    template<class X2, class Ref2> friend class ExpansionIterator;
    typedef ExpansionIterator<X,Ref> Iter;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef X data_type;

    typedef int difference_type;
  public:
    ExpansionIterator(size_type n, void* p) : _n(n), _p(reinterpret_cast<word_type*>(p)),
        _x(reinterpret_cast<data_type*>(_p+MultiIndex::_word_size(_n))) { }
    template<class Ref2> ExpansionIterator(const ExpansionIterator<X,Ref2>& i)
        : _n(i._n), _p(reinterpret_cast<word_type*>(i._p)),
          _x(reinterpret_cast<data_type*>(_p+MultiIndex::_word_size(_n))) { }
    template<class Ref2> bool equal(const ExpansionIterator<X,Ref2>& i) const {
        return _p==i._p; }
    template<class Ref2> difference_type distance_to(const ExpansionIterator<X,Ref2>& i) const {
        return (i._p-_p)/difference_type(MultiIndex::_element_size(_n)); }
    Ref* operator->() const {
        return reinterpret_cast<Ref*>(const_cast<Iter*>(this)); }
    Ref& operator*() const {
        return reinterpret_cast<Ref&>(const_cast<Iter&>(*this)); }
    Iter& increment() {
        return advance(1); }
    Iter& decrement() {
        return advance(-1); }
    Iter& advance(difference_type m) {
        _p+=m*MultiIndex::_element_size(_n);
        _x=reinterpret_cast<data_type*>(_p+MultiIndex::_word_size(_n));
        return *this; }
  private:
    size_type _n;
    word_type* _p;
    data_type* _x;
};


template<class X, class Ref>
inline std::ostream& operator<<(std::ostream& os, const ExpansionIterator<X,Ref>& piter) {
    return os << "(n="<<piter._n<<", p="<<(void*)piter._p<<", a="<<piter->first<<", x="<<*piter._x<<")";
}

template<class X>
class Expansion
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
    typedef ExpansionValue<X> value_type;
    typedef ExpansionReference<X> reference;
    typedef ExpansionReference<const X> const_reference;
    typedef ExpansionReference<X>* pointer;
    typedef ExpansionReference<const X>* const_pointer;
    typedef ExpansionIterator<X,ExpansionReference<X> > iterator;
    typedef ExpansionIterator<X,ExpansionReference<const X> > const_iterator;
  public:
    Expansion() : _argument_size() { }
    Expansion(unsigned int as) : _argument_size(as) { }
    Expansion(unsigned int as, unsigned int deg, double c0, ...);
    Expansion(unsigned int as, unsigned int nnz, int a00, ...);
    template<class XX> Expansion(const std::map<MultiIndex,XX>& m);
    template<class XX> Expansion(const Expansion<XX>& p);

    static Expansion<X> variable(unsigned int as, unsigned int i);
    static Vector< Expansion<X> > variables(unsigned int n);

    void swap(Expansion<X>& p) {
        std::swap(this->_argument_size,p._argument_size);
        std::swap(this->_coefficients,p._coefficients); }

    bool operator==(const Expansion<X>& p) const { return this->_coefficients == p._coefficients; }
    bool operator!=(const Expansion<X>& p) const { return this->_coefficients != p._coefficients; }

    unsigned int argument_size() const { return this->_argument_size; }
    unsigned int number_of_nonzeros() const { return _coefficients.size()/element_size(); }
    unsigned int degree() const { return (--this->_coefficients.end())->first.degree(); }

    void reserve(size_type nnz) { this->_coefficients.reserve(nnz*element_size()); }

    void insert(const MultiIndex& a, const Real& c) {
        assert(a.size()==_argument_size); this->_insert(a,c); }
    void append(const MultiIndex& a, const Real& c) {
        this->_append(a,c); }
    void append(const MultiIndex& a1, const MultiIndex& a2, const Real& c) {
        this->_append(a1,a2,c); }

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
    iterator find(const MultiIndex& a) {
        iterator iter=this->begin(); 
        while(iter!=this->end() && iter->first!=a) { ++iter; } 
        return iter; }

    const_iterator begin() const { return const_iterator(_argument_size,(byte_type*)_begin()); }
    const_iterator end() const { return const_iterator(_argument_size,(byte_type*)_end()); }
    const_iterator find(const MultiIndex& a) const {
        const_iterator iter=this->begin(); 
        while(iter!=this->end() && iter->first!=a) { ++iter; } 
        return iter; };

    void erase(iterator iter) { iter->second=0.0; }
    void clear() { _coefficients.clear(); }

    void cleanup() {
        if(this->_coefficients.empty()) { return; }
        std::sort(this->begin(),this->end());
        iterator begin=this->begin();
        iterator curr=this->begin();
        iterator end=this->end();
        iterator next=curr; ++next;
        while(next!=end) {
            if(curr->first==next->first) {
                static_cast<X&>(curr->second)+=static_cast<const X&>(next->second); ++next; }
            else {
                ++curr; *curr=*next; ++next; }
        }
        ++curr;
        this->_coefficients.resize((curr-begin)*this->element_size());
    }

    void check() { }

    template<class Y> Y evaluate(const Vector<Y>& x) const;
    Expansion<X> embed(unsigned int new_size, unsigned int start) const;
  public:
    size_type vector_size() const {
        return _coefficients.size(); }
    size_type word_size() const {
        return (_argument_size*+sizeof_word)/sizeof_word; }
    size_type element_size() const {
        return (_argument_size+sizeof_word)/sizeof_word+sizeof_data/sizeof_word; }
  public:
    const word_type* _begin() const { return _coefficients.begin().operator->(); }
    const word_type* _end() const { return _coefficients.end().operator->(); }
    word_type* _begin() { return _coefficients.begin().operator->(); }
    word_type* _end() { return _coefficients.end().operator->(); }
    iterator _insert(const MultiIndex& a, const Real& x) {
        this->_append(a,x);
        iterator iter=this->end(); --iter; iterator prev=iter;
        while(iter!=this->begin()) {  
            --prev; if(prev->first<a) { break; } else { *iter=*prev; --iter; } }
        iter->first=a; iter->second=x; return iter; }
    void _append(const MultiIndex& a, const Real& x) {
        _coefficients.resize(_coefficients.size()+element_size());
        word_type* vp=_end()-element_size(); const word_type* ap=a.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end())-1; *xp=x; }
    void _append(const MultiIndex&  a1, const MultiIndex&  a2, const Real& x) {
        _coefficients.resize(_coefficients.size()+element_size());
        word_type* vp=_end()-element_size();
        const word_type* ap1=a1.word_begin(); const word_type* ap2=a2.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap1[j]+ap2[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end())-1; *xp=x; }
  private:
    size_type _argument_size;
    std::vector<word_type> _coefficients;
};

// Disable construction of Expansion<Rational> since above implementation only
// works for "plain old data" types
#if defined HAVE_GMPXX_H and defined ARIADNE_NUMERIC_H
template<> class Expansion<Rational>;
#endif

#endif


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
}

template<class X> template<class XX> inline
Expansion<X>::Expansion(const Expansion<XX>& p)
    : _argument_size(p.argument_size())
{
    for(typename Expansion<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->first,X(iter->second)); }
}


template<class X> Expansion<X> Expansion<X>::variable(unsigned int n, unsigned int i) {
    Expansion<X> p(n); p[MultiIndex::zero(n)]=0.0; p[MultiIndex::unit(n,i)]=X(1);
    return p;
}

template<class X> Vector< Expansion<X> > Expansion<X>::variables(unsigned int n) {
    Vector< Expansion<X> > pv(n,Expansion<X>(n)); MultiIndex a(n);
    for(uint i=0; i!=n; ++i) { a[i]=1; pv[i][a]=1.0; a[i]=0; }
    return pv; 
}




template<class X> template<class Y>
Y Expansion<X>::evaluate(const Vector<Y>& x) const
{
    ARIADNE_ASSERT(this->argument_size()==x.size());
    Y zero = x[0]; zero*=0;
    Y one = zero; one+=1;

    Y r=zero;

    for(typename Expansion<X>::const_iterator iter=this->begin();
        iter!=this->end(); ++iter)
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
Y evaluate(const Expansion<X>& x, const Vector<Y>& y)
{
    return x.evaluate(y);
}

template<class X, class Y>
Vector<Y> evaluate(const Vector< Expansion<X> >& x, const Vector<Y>& y)
{
    Vector<Y> r(x.size());
    for(unsigned int i=0; i!=x.size(); ++i) {
        r[i]=x[i].evaluate(y);
    }
    return r;
}



template<class X>
Expansion<X> Expansion<X>::embed(unsigned int new_size, unsigned int start) const
{
    ARIADNE_ASSERT(this->argument_size()+start<=new_size);
    Expansion<X> r(new_size);
    uint old_size=this->argument_size();
    MultiIndex old_index(old_size);
    MultiIndex new_index(new_size);
    for(typename Expansion<X>::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
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
std::ostream& operator<<(std::ostream& os, const Expansion<X>& p) {
    if(p.begin()==p.end()) {
        return os << "{"<<MultiIndex::zero(p.argument_size())<<":0.0}"; }

    os << "{";
    for(typename Expansion<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? "" : ",");
        for(unsigned int i=0; i!=iter->first.size(); ++i) {
            os << (i==0?" ":",") << int(iter->first[i]); }
        os << ":" << iter->second; }
    return os << " }";
}


inline Expansion<Float> midpoint(const Expansion<Interval>& pse) {
    Expansion<Float> r(pse.argument_size());
    for(Expansion<Interval>::const_iterator iter=pse.begin(); iter!=pse.end(); ++iter) {
        r.append(iter->first,midpoint(iter->second)); }
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
