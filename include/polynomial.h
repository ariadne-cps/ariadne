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
#include "expansion.h"



namespace Ariadne {

//! \brief A polynomial with coefficients of some type \a X.
template<class X>
class Polynomial
{
    template<class XX> friend class Polynomial;
  public:
    typedef typename Expansion<X>::size_type size_type;
    typedef typename Expansion<X>::value_type value_type;
    typedef typename Expansion<X>::reference reference;
    typedef typename Expansion<X>::const_reference const_reference;
    typedef typename Expansion<X>::iterator iterator;
    typedef typename Expansion<X>::const_iterator const_iterator;
  public:
    //! \brief The zero polynomial in \a as variables.
    Polynomial(unsigned int as=0u) : _expansion(as) { }
    //! \brief Copy/conversion constructor.
    template<class XX> Polynomial(const Polynomial<XX>& p) : _expansion(p._expansion) { }
    //! \brief Copy/conversion constructor.
    template<class XX> explicit Polynomial(const Expansion<XX>& e) : _expansion(e) { }
    //! \brief A dense polynomial with coefficients given by a list of doubles.
    Polynomial(unsigned int as, unsigned int deg, double c0, ...);
    //! \brief A sparse polynomial with coefficients given by a list of indices and coefficients.
    Polynomial(unsigned int as, unsigned int nnz, int a00, ...);

    static Polynomial<X> constant(unsigned int as, const X& c) {
        Polynomial<X> r(as); r[MultiIndex::zero(as)]=c; return r; }
    static Polynomial<X> variable(unsigned int as, unsigned int j) {
        Polynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r; }
    static Vector< Polynomial<X> > variables(unsigned int as) {
        Vector< Polynomial<X> > r(as); for(unsigned int i=0; i!=as; ++i) { r[i]=variable(as,i); } return r; }

    template<class XX> bool operator==(const Polynomial<XX>& p) const {
        const_cast<Polynomial<X>*>(this)->cleanup(); const_cast<Polynomial<XX>&>(p).cleanup();
        return this->_expansion==p._expansion; }
    template<class XX> bool operator!=(const Polynomial<XX>& p) const {
        return !(*this==p); }

    size_type argument_size() const { return this->_expansion.argument_size(); }
    size_type number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
    const X& value() const { return this->_expansion[MultiIndex::zero(this->argument_size())]; }
    X& operator[](const MultiIndex& a) { return this->_expansion[a]; }
    const X& operator[](const MultiIndex& a) const { return this->_expansion[a]; }

    iterator begin() { return this->_expansion.begin(); }
    iterator end() { return this->_expansion.end(); }
    iterator find(const MultiIndex& a) { return this->_expansion.find(a); }
    const_iterator begin() const { return this->_expansion.begin(); }
    const_iterator end() const { return this->_expansion.end(); }
    const_iterator find(const MultiIndex& a) const { return this->_expansion.find(a); }

    void append(const MultiIndex& a1, const MultiIndex& a2, const X& x) { this->_expansion.append(a1,a2,x); }
    void append(const MultiIndex& a, const X& x) { this->_expansion.append(a,x); }
    void insert(const MultiIndex& a, const X& x) { this->_expansion.insert(a,x); }
    void reserve(size_type n) { this->_expansion.reserve(n); }
    void clear() { this->_expansion.clear(); }
    void erase(iterator iter) { this->_expansion.erase(iter); }

    const Expansion<X>& expansion() const { return this->_expansion; }
    Polynomial<X>& truncate(unsigned int d);

    void cleanup();
    void check() const;
  private:
    Expansion<X> _expansion;
};

template<class X>
Polynomial<X>::Polynomial(unsigned int as, unsigned int deg, double c0, ...)
    : _expansion(as)
{
    MultiIndex a(as); double x;
    va_list args; va_start(args,c0);
    while(a.degree()<=deg) {
        if(a.degree()==0) { x=c0; }
        else { x=va_arg(args,double); }
        if(x!=0) { this->_expansion.append(a,x); }
        ++a;
    }
    va_end(args);
    this->cleanup();
}

template<class X>
Polynomial<X>::Polynomial(unsigned int as, unsigned int nnz, int a00, ...)
    : _expansion(as)
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
        if(x!=0) { this->_expansion.append(a,x); }
    }
    va_end(args);
    this->cleanup();
}

template<class X>
void Polynomial<X>::cleanup()
{
    Polynomial<X>* self=const_cast<Polynomial<X>*>(this);
    std::sort(self->_expansion.begin(), self->_expansion.end());
    iterator new_end=unique_key(self->_expansion.begin(), self->_expansion.end(), std::plus<X>());
    self->_expansion.resize(new_end-self->_expansion.begin());
}

template<class X> inline
void Polynomial<X>::check() const
{
    this->_expansion.check();
}

template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p) {
    return p; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p) {
    Polynomial<X> r(p.argument_size());  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { r[iter->key()]-=iter->data(); } return r; }

template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]+=c; return r; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]-=c; return r; }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p, const X& c) {
    if(c==0) { return Polynomial<X>(p.argument_size()); }
    Polynomial<X> r(p); typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->data()*=c; } return r; }
template<class X> inline Polynomial<X> operator/(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->data()/=c; } return r; }

template<class X> inline Polynomial<X> operator+(const X& c, const Polynomial<X>& p) {
    return p+c; }
template<class X> inline Polynomial<X> operator-(const X& c, const Polynomial<X>& p) {
    return (-p)+c; }
template<class X> inline Polynomial<X> operator*(const X& c, const Polynomial<X>& p) {
    return p*c; }

template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->key()]+=iter->data(); } return r; }
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->key()]-=iter->data(); } return r; }
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(Iter iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->key()+iter2->key();
            r[a]+=iter1->data()*iter2->data(); } } return r; }

template<class X, class XX> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Polynomial<XX>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->key()]+=iter->data(); } return p; }
template<class X, class XX> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Polynomial<XX>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->key()]-=iter->data(); } return p; }

template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]+=c; return p; }
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]-=c; return p; }
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const X& c) {
    typedef typename Polynomial<X>::iterator Iter;
    if(c==0) { p.clear(); }
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->data()*=c; } return p; }
template<class X> inline Polynomial<X>& operator/=(Polynomial<X>& p, const X& c) {
    typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->data()/=c; } return p; }

template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const int& c) {
    p+=X(c); return p; }
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const int& c) {
    p*=X(c); return p; }

inline Polynomial<Float> midpoint(const Polynomial<Interval>& p) {
    Polynomial<Float> r(p.argument_size());
    for(Polynomial<Interval>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        r.append(iter->key(),midpoint(iter->data())); }
    return r;
}

inline Vector< Polynomial<Float> > midpoint(const Vector< Polynomial<Interval> >& p) {
    Vector< Polynomial<Float> > r(p.size());
    for(uint i=0; i!=p.size(); ++i) {
        r[i]=midpoint(p[i]); }
    return r;
}


template<class X, class Y> inline
Y evaluate(const Polynomial<X>& p, const Vector<Y>& x)
{
    return evaluate(p.expansion(),x);
}


template<class X, class Y> inline
Vector<Y> evaluate(const Vector< Polynomial<X> >& p, const Vector<Y>& x)
{
    ARIADNE_ASSERT(p.size()>0 && p[0].argument_size()==x.size());
    Y zero = x[0]; zero*=0;
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
Polynomial<X> embed(unsigned int before_size, const Polynomial<X>& x, unsigned int after_size)
{
    return x.expansion().embed(before_size,after_size);
}

template<class X>
Polynomial<X>
derivative(const Polynomial<X>& p, uint j) {
    Polynomial<X> r(p.argument_size());
    MultiIndex ar(p.argument_size());
    for(typename Polynomial<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        const MultiIndex& ap=iter->key();
        if(ap[j]>0) {
            const X& val=iter->data();
            ar=ap;
            ar[j]-=1;
            r.append(ar,val*ap[j]);
        }
    }
    return r;
}

template<class X>
Vector< Polynomial<X> >
derivative(const Vector< Polynomial<X> >& p, uint j) {
    Vector< Polynomial<X> > r(p.size());
    for(uint i=0; i!=p.size(); ++i) {
        r[i]=derivative(p[i],j);
    }
    return r;
}

template<class X>
Polynomial<X>
antiderivative(const Polynomial<X>& p, uint j) {
    Polynomial<X> r(p.argument_size());
    MultiIndex ar(p.argument_size());
    for(typename Polynomial<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        const MultiIndex& ap=iter->key();
        const X& val=iter->data();
        ar=ap;
        ++ar[j];
        r.append(ar,val/ar[j]);
    }
    return r;
}

template<class X>
Vector< Polynomial<X> >
antiderivative(const Vector< Polynomial<X> >& p, uint j) {
    Vector< Polynomial<X> > r(p.size());
    for(uint i=0; i!=p.size(); ++i) {
        r[i]=antiderivative(p[i],j);
    }
    return r;
}

template<class X>
Polynomial<X>
truncate(const Polynomial<X>& p, uint d) {
    Polynomial<X> r(p.argument_size());
    for(typename Polynomial<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        if(iter->key().degree()<=d && iter->data()!=0) {
            r.append(iter->key(),iter->data());
        }
    }
    return r;
}

template<class X>
Vector< Polynomial<X> >
truncate(const Vector< Polynomial<X> >& p, uint d) {
    Vector< Polynomial<X> > r(p.size());
    for(uint i=0; i!=p.size(); ++i) {
        r[i]=truncate(p[i],d);
    }
    return r;
}

template<class X>
Polynomial<X>& Polynomial<X>::truncate(uint d) {
    Polynomial<X> r=Ariadne::truncate(*this);
    this->swap(r);
    return *this;
}

// A polynomial approximation to a flow of degree order
template<class X>
Vector< Polynomial<X> > flow(const Vector< Polynomial<X> >& p, uint order)
{
    uint n=p.size();
    Vector< Polynomial<X> > p0(n,n+1);
    for(uint i=0; i!=n; ++i) { p0[i][MultiIndex::unit(n+1,i)]=1.0; }

    Vector< Polynomial<X> > r=p0;
    for(uint k=0; k!=order; ++k) {
        r=compose(p,r);
        for(uint i=0; i!=n; ++i) {
            r[i].cleanup();
            r[i]=truncate(r[i],k+1);
            r[i]=antiderivative(r[i],n);
            r[i]+=p0[i];
        }
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
        for(unsigned int i=0; i!=iter->key().size(); ++i) {
            os << (i==0?" ":",") << int(iter->key()[i]); }
        os << ":" << iter->data(); }
    return os << " }";
}


template<class X> Vector< Polynomial<X> > operator*(const Polynomial<X>& p, const Vector<Float> e) {
    Vector< Polynomial<X> > r(e.size(),Polynomial<X>(p.argument_size()));
    for(uint i=0; i!=r.size(); ++i) { r[i]=p; r[i]*=X(e[i]); }
    return r;
}



inline Polynomial<Interval> operator+(const Polynomial<Interval>& p, const Float& c) {
    return p+Interval(c); }
inline Polynomial<Interval> operator-(const Polynomial<Interval>& p, const Float& c) {
    return p-Interval(c); }
inline Polynomial<Interval> operator*(const Polynomial<Interval>& p, const Float& c) {
    return p*Interval(c); }
inline Polynomial<Interval> operator/(const Polynomial<Interval>& p, const Float& c) {
    return p/Interval(c); }
inline Polynomial<Interval> operator+(const Float& c, const Polynomial<Interval>& p) {
    return Interval(c)+p; }
inline Polynomial<Interval> operator-(const Float& c, const Polynomial<Interval>& p) {
    return Interval(c)-p; }
inline Polynomial<Interval> operator*(const Float& c, const Polynomial<Interval>& p) {
    return Interval(c)*p; }

inline Polynomial<Interval>& operator+=(Polynomial<Interval>& p, const Float& c) {
    return p+=Interval(c); }
inline Polynomial<Interval>& operator-=(Polynomial<Interval>& p, const Float& c) {
    return p-=Interval(c); }
inline Polynomial<Interval>& operator*=(Polynomial<Interval>& p, const Float& c) {
    return p*=Interval(c); }
inline Polynomial<Interval>& operator/=(Polynomial<Interval>& p, const Float& c) {
    return p/=Interval(c); }
}

#endif /* ARIADNE_POLYNOMIAL_H */
