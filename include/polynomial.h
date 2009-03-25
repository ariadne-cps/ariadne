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


//#define ARIADNE_USE_MAP_POLYNOMIAL
#define ARIADNE_USE_ARRAY_POLYNOMIAL

namespace Ariadne {

#if defined ARIADNE_USE_MAP_POLYNOMIAL || defined DOXYGEN

//! \brief A polynomial with coefficients of some type \a X.
template<class X>
class Polynomial
{
  public:
    typedef MultiIndex::size_type size_type;
    typedef typename std::map<MultiIndex,X>::value_type value_type;
    typedef typename std::map<MultiIndex,X>::reference reference;
    typedef typename std::map<MultiIndex,X>::pointer pointer;
    typedef typename std::map<MultiIndex,X>::const_reference const_reference;
    typedef typename std::map<MultiIndex,X>::const_pointer const_pointer;
    typedef typename std::map<MultiIndex,X>::key_type key_type;
    typedef X data_type;
    typedef typename std::map<MultiIndex,X>::iterator iterator;
    typedef typename std::map<MultiIndex,X>::const_iterator const_iterator;
  private:
    static X _zero;
  public:
    //! \brief A polynomial in no arguments; can only take a constant value.
    Polynomial() : _argument_size() { }
    //! \brief The zero polynomial in \a as arguments.
    Polynomial(unsigned int as) : _argument_size(as) { }
    //! \brief A dense polynomial with coefficients given by a list of doubles.
    Polynomial(unsigned int as, unsigned int deg, double c0, ...);
    //! \brief Construct from a map of (intex,coefficient) pairs.
    template<class XX> Polynomial(const std::map<MultiIndex,XX>& m);
    //! \brief Construct from a polynomial with possibly different coefficient types.
    template<class XX> Polynomial(const Polynomial<XX>& p);

    //! \brief The polynomial \f$x_i\f$ in \a as variables.
    static Polynomial<X> variable(unsigned int as, unsigned int i);
    //! \brief A vector of polynomials of the form \f$(x_0,\ldots,x_{n-1})\f$.
    static Vector< Polynomial<X> > variables(unsigned int n);

    //! \brief Swap with another polynomial. Included for efficiency.
    void swap(Polynomial<X>& p) {
        std::swap(this->_argument_size,p._argument_size);
        std::swap(this->_coefficients,p._coefficients); }

    //! \brief Equality operator.
    bool operator==(const Polynomial<X>& p) const { return this->_coefficients == p._coefficients; }
    //! \brief Inequality operator.
    bool operator!=(const Polynomial<X>& p) const { return this->_coefficients != p._coefficients; }

    //! \brief The number of argument variables of the polynomial.
    unsigned int argument_size() const { return this->_argument_size; }
    //! \brief The number of nonzero terms.
    unsigned int number_of_nonzeros() const { return this->_coefficients.size(); }
    //! \brief The degree of the polynomial.
    unsigned int degree() const { return (--this->_coefficients.end())->first.degree(); }

    //! \brief Reserve space for \a nnz nonzero elements.
    void reserve(size_type nnz) { }
}
    //! \brief Insert the term \f$cx^a\f$ into the polynomial.
    //! The polynomial must not already have a term in \f$x^a\f$.
    void insert(const MultiIndex& a, const Real& c) {
        assert(a.size()==_argument_size); this->_coefficients.insert(std::make_pair(a,c)); }
    //! \brief Insert the term \f$cx^a\f$ into the polynomial.
    //! The index \a a must be higher than the index of all currently existing terms.
    void append(const MultiIndex& a, const Real& c) {
        this->insert(a,c); }
    //! \brief Insert the term \f$cx^{a_1+a_2}\f$ into the polynomial.
    //! The index \f$a_1+a_2\f$ must be higher than the index of all currently existing terms.
    void append(const MultiIndex& a1, const MultiIndex& a2, const Real& c) {
        this->insert(a1+a2,c)); }

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

    //! \brief Clean up the representation of the polynomial, including
    //! sorting the terms and removing zero elements if necessary.
    void cleanup() { }
    //! \brief Remove all nonzero terms.
    void clear() { this->_coefficients.clear(); }
    //! \brief Check the representation of the polynomial for consistency.
    void check() { for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        ARIADNE_ASSERT(iter->first.size()==this->_argument_size); } }

    //! \brief Evaluate the polynomial on a vector \a x over a ring.
    template<class Y> Y evaluate(const Vector<Y>& x) const;
  private:
    size_type _argument_size;
    std::map<MultiIndex,X> _coefficients;
};

#elif defined ARIADNE_USE_ARRAY_POLYNOMIAL


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
    Reference& operator=(Reference<MultiIndex> a) {
        assert(_n==a._n); this->_assign(a._p); return *this; }
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
struct PolynomialValue {
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
struct PolynomialReference {
    PolynomialReference<X>& operator=(const PolynomialReference<X>& dr) {
        first=dr.first; second=dr.second; return *this; }
    PolynomialReference<X>& operator=(const PolynomialValue<X>& dv) {
        first=dv.first; second=dv.second; return *this; }
    operator PolynomialValue<X>() const { 
        return PolynomialValue<X>(this->first,this->second); }
    Reference<MultiIndex> first;
    Reference<X> second;
};

template<class X>
struct PolynomialReference<const X> {
    PolynomialReference(const PolynomialReference<X>& r)
        : first(r.first), second(r.second) { }
    operator PolynomialValue<X>() const { 
        return PolynomialValue<X>(this->first,this->second); }
    Reference<const MultiIndex> first;
    Reference<const X> second;
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
inline bool operator<(const PolynomialReference<const X>& dr1, const PolynomialReference<const X>& dr2) {
    return dr1.first < dr2.first;
}

template<class X>
inline std::ostream& operator<<(std::ostream& os, PolynomialReference<const X> dr) {
    return os << dr.first << ":" << dr.second;
}

template<class X,class Ref> class PolynomialIterator;

template<class X, class Ref>
class PolynomialIterator
     : public boost::iterator_facade<PolynomialIterator<X,Ref>,
                                     PolynomialValue<X>,
                                     boost::random_access_traversal_tag,
                                     Ref >

{
    template<class X2, class Ref2> friend class PolynomialIterator;
    typedef PolynomialIterator<X,Ref> Iter;
  public:
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type byte_type;
    typedef MultiIndex::word_type word_type;
    typedef X data_type;

    typedef int difference_type;
  public:
    PolynomialIterator(size_type n, void* p) : _n(n), _p(reinterpret_cast<word_type*>(p)), 
        _x(reinterpret_cast<data_type*>(_p+MultiIndex::_word_size(_n))) { }
    template<class Ref2> PolynomialIterator(const PolynomialIterator<X,Ref2>& i)
        : _n(i._n), _p(reinterpret_cast<word_type*>(i._p)),
          _x(reinterpret_cast<data_type*>(_p+MultiIndex::_word_size(_n))) { }
    template<class Ref2> bool equal(const PolynomialIterator<X,Ref2>& i) const {
        return _p==i._p; }
    template<class Ref2> difference_type distance_to(const PolynomialIterator<X,Ref2>& i) const {
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
inline std::ostream& operator<<(std::ostream& os, const PolynomialIterator<X,Ref>& piter) {
    return os << "(n="<<piter._n<<", p="<<(void*)piter._p<<", a="<<piter->first<<", x="<<*piter._x<<")";
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
    typedef PolynomialReference<const X> const_reference;
    typedef PolynomialReference<X>* pointer;
    typedef PolynomialReference<const X>* const_pointer;
    typedef PolynomialIterator<X,PolynomialReference<X> > iterator;
    typedef PolynomialIterator<X,PolynomialReference<const X> > const_iterator;
  public:
    Polynomial() : _argument_size() { }
    Polynomial(unsigned int as) : _argument_size(as) { }
    Polynomial(unsigned int as, unsigned int deg, double c0, ...);
    template<class XX> Polynomial(const std::map<MultiIndex,XX>& m);
    template<class XX> Polynomial(const Polynomial<XX>& p);

    static Polynomial<X> variable(unsigned int as, unsigned int i);
    static Vector< Polynomial<X> > variables(unsigned int n);

    void swap(Polynomial<X>& p) {
        std::swap(this->_argument_size,p._argument_size);
        std::swap(this->_coefficients,p._coefficients); }

    bool operator==(const Polynomial<X>& p) const { return this->_coefficients == p._coefficients; }
    bool operator!=(const Polynomial<X>& p) const { return this->_coefficients != p._coefficients; }

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

    void cleanup();
    void clear() { _coefficients.clear(); }
    void check() { };

    template<class Y> Y evaluate(const Vector<Y>& x) const;
  public:
    unsigned int vector_size() const {
        return _coefficients.size(); }
    unsigned int word_size() const {
        return (_argument_size*+sizeof_word)/sizeof_word; }
    unsigned int element_size() const {
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
        word_type* vp=_end()-element_size(); const word_type* ap1=a1.word_begin(); const word_type* ap2=a2.word_begin();
        for(unsigned int j=0; j!=word_size(); ++j) { vp[j]=ap1[j]+ap2[j]; }
        data_type* xp=reinterpret_cast<data_type*>(this->_end())-1; *xp=x; }
  private:
    unsigned int _argument_size;
    std::vector<word_type> _coefficients;
};


template<class X>
inline void Polynomial<X>::cleanup() {
    if(this->_coefficients.empty()) { return; }
    std::sort(this->begin(),this->end());
    iterator begin=this->begin();
    iterator curr=this->begin();
    iterator end=this->end();
    iterator next=curr; ++next;
    while(next!=end) {
        if(curr->first==next->first) { curr->second+=next->second; ++next; }
        else { ++curr; *curr=*next; ++next; }
    }
    ++curr;
    this->_coefficients.resize((curr-begin)*this->element_size());
}

#else

#error "No Polynomial class activated."

#endif


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
    this->cleanup();
}

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
    std::cerr<<"Polynomial(Polynomial)"<<std::endl;
    for(typename Polynomial<XX>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        this->append(iter->first,X(iter->second)); }
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

template<class X> inline Polynomial<X> operator+(const X& c, const Polynomial<X>& p) {
    return p+c; }
template<class X> inline Polynomial<X> operator-(const X& c, const Polynomial<X>& p) {
    return (-p)+c; }
template<class X> inline Polynomial<X> operator*(const X& c, const Polynomial<X>& p) {
    return p*c; }

template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->first]+=iter->second; } return r; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->first]-=iter->second; } return r; }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1.argument_size()); 
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(Iter iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->first+iter2->first;
            r[a]+=iter1->second*iter2->second; } } return r; }

template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Polynomial<X>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->first]+=iter->second; } return p; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Polynomial<X>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->first]-=iter->second; } return p; }

template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Float& c) {
    p[MultiIndex(p.argument_size())]+=c; return p; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Float& c) {
    p[MultiIndex(p.argument_size())]-=c; return p; }
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const Float& c) {
    typedef typename Polynomial<X>::iterator Iter;
    if(c==0) { p.clear(); }
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
    ARIADNE_ASSERT(p.size()>0 && p[0].argument_size()==x.size());\
    Y zero = x[0]; zero*=0.0;
    Vector<Y> r(p.size(),zero);
    for(uint i=0; i!=p.size(); ++i) {
        r[i]=evaluate(p[i],x);
    }
    return r;
}

template<class X> template<class Y> inline
Y Polynomial<X>::evaluate(const Vector<Y>& x) const
{
    return Ariadne::evaluate(*this,x);
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
    if(p.begin()==p.end()) {
        return os << "{"<<MultiIndex::zero(p.argument_size())<<":0.0}"; }

    os << "{";
    for(typename Polynomial<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        os << (iter==p.begin() ? "" : ",");
        for(unsigned int i=0; i!=iter->first.size(); ++i) {
            os << (i==0?" ":",") << int(iter->first[i]); }
        os << ":" << iter->second; }
    return os << " }";
}


template<class X> Vector< Polynomial<X> > operator*(const Vector<Float> e, const Polynomial<X>& p) {
    Vector< Polynomial<X> > r(e.size(),Polynomial<X>(p.argument_size()));
    for(uint i=0; i!=r.size(); ++i) { r[i]=p; r[i]*=X(e[i]); }
    return r;
}

}

#endif /* ARIADNE_POLYNOMIAL_H */
