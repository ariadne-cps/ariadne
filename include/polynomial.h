/***************************************************************************
 *            polynomial.h
 *
 *  Copyright 2008-11  Pieter Collins, Alberto Casagrande
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
#include "taylor_model.h"
#include "differential.h"


namespace Ariadne {

template<class T> class Array;

//! \brief A monomial with coefficients of some type \a X.
template<class X>
class Monomial
    : public ExpansionValue<X>
{
    Monomial(const MultiIndex& a, const X& x) : ExpansionValue<X>(a,x) { }
    Monomial(const ExpansionValue<X>& v) : ExpansionValue<X>(v) { }
};

//! \brief A polynomial with coefficients of some type \a X.
template<class X>
class Polynomial
{
    template<class XX> friend class Polynomial;
  public:
    typedef typename Expansion<X>::size_type size_type;
    typedef typename Expansion<X>::smoothness_type smoothness_type;
    typedef typename Expansion<X>::value_type value_type;
    typedef typename Expansion<X>::reference reference;
    typedef typename Expansion<X>::const_reference const_reference;
    typedef typename Expansion<X>::iterator iterator;
    typedef typename Expansion<X>::const_iterator const_iterator;
  public:
    //@{
    //! \name Constructors

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
    //@}

    //! \brief Create a constant polynomial in \a as variables with value \a c.
    static Polynomial<X> constant(unsigned int as, const X& c) {
        Polynomial<X> r(as); r[MultiIndex::zero(as)]=c; return r; }
    //! \brief Create a polynomial in \a as variables which returns the value of the \a j<sup>th</sup> variable.
    static Polynomial<X> variable(unsigned int as, unsigned int j) {
        ARIADNE_ASSERT(j<as); Polynomial<X> r(as); r[MultiIndex::unit(as,j)]=1; return r; }
    //! \brief Create an array of polynomials in \a as variables,
    //! the i<sup>th</sup> of  which returns the value of the i<sup>th</sup> variable.
    static array< Polynomial<X> > variables(unsigned int as) {
        array< Polynomial<X> > r(as); for(unsigned int i=0; i!=as; ++i) { r[i]=variable(as,i); } return r; }

    //@{
    //! \name Comparisons

    //! \brief Equality operator.
    template<class XX> bool operator==(const Polynomial<XX>& p) const {
        const_cast<Polynomial<X>*>(this)->cleanup(); const_cast<Polynomial<XX>&>(p).cleanup();
        return this->_expansion==p._expansion; }
    //! \brief Inequality operator.
    template<class XX> bool operator!=(const Polynomial<XX>& p) const {
        return !(*this==p); }
    //@}

    //@{
    //! \name Data access

    //! \brief The number of variables in the argument of the polynomial.
    size_type argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The number of structural nonzero terms.
    size_type number_of_nonzeros() const { return this->_expansion.number_of_nonzeros(); }
    //! \brief The order of the highest term.
    size_type degree() const { return this->_expansion.degree(); }
    //! \brief The value of the polynomial at zero.
    const X& value() const { return this->_expansion[MultiIndex::zero(this->argument_size())]; }
    //! \brief A reference to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    X& operator[](const MultiIndex& a) { return this->_expansion[a]; }
    //! \brief A constant referent to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    const X& operator[](const MultiIndex& a) const { return this->_expansion[a]; }
    //! \brief A constant reference to the raw data expansion.
    const Expansion<X>& expansion() const { return this->_expansion; }
    //! \brief A reference to the raw data expansion.
    Expansion<X>& expansion() { return this->_expansion; }
    //@}

    //@{
    //! \name Iterators

    //! \brief An iterator to the beginning of the list of terms.
    iterator begin() { return this->_expansion.begin(); }
    //! \brief An iterator to the end of the list of terms..
    iterator end() { return this->_expansion.end(); }
    //! \brief An iterator to the term in \f$x^a\f$. Returns \c end() if there is no term in \a a.
    iterator find(const MultiIndex& a) { return this->_expansion.find(a); }
    //! \brief A constant iterator to the beginning of the list of terms.
    const_iterator begin() const { return this->_expansion.begin(); }
    //! \brief A constant iterator to the end of the list of terms.
    const_iterator end() const { return this->_expansion.end(); }
    //! \brief An iterator to the term in \f$x^a\f$. Returns \c end() if there is no term in \a a.
    const_iterator find(const MultiIndex& a) const { return this->_expansion.find(a); }
    //@}


    //@{
    //! \name Modifying operations

    //! \brief Append the term \f$c x^{a_1+a_2}\f$ to the list of terms.
    void append(const MultiIndex& a1, const MultiIndex& a2, const X& c) { this->_expansion.append(a1,a2,c); }
    //! \brief Append the term \f$c x^{a_1}\f$ to the list of terms.
    void append(const MultiIndex& a, const X& c) { this->_expansion.append(a,c); }
    //! \brief Insert the term \f$c x^{a_1}\f$ into a sorted list of terms.
    void insert(const MultiIndex& a, const X& c) { this->_expansion.insert(a,c); }
    //! \brief Reserve space for a total of \a n terms.
    void reserve(size_type n) { this->_expansion.reserve(n); }
    //! \brief Remove the term pointed to by \a iter. May be expensive if the term is near the beginning of the list of terms.
    void erase(iterator iter) { this->_expansion.erase(iter); }
    //! \brief Set the polynomial to zero.
    void clear() { this->_expansion.clear(); }
    //! \brief Remove all zero terms from the expansion, and order the expansion lexicographically by term.
    void cleanup();
    //@}

    //@{
    //! \name Modifying operators

    //! \brief Truncate to degree \a d.
    Polynomial<X>& truncate(smoothness_type d);
    //! \brief Differentiate with respect to the \a j<sup>th</sup> variable.
    Polynomial<X>& differentiate(size_type j);
    //! \brief Antidifferentiate (integrate) with respect to the \a j<sup>th</sup> variable.
    Polynomial<X>& antidifferentiate(size_type j);
    //@}

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
Polynomial<X>&
Polynomial<X>::differentiate(size_type j) {
    for(typename Polynomial<X>::iterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex& a=iter->key();
        X& c=iter->data();
        c*=a[j];
        if(a[j]!=0u) { ++a[j]; }
    }
    return *this;
}

template<class X>
Polynomial<X>&
Polynomial<X>::antidifferentiate(size_type j) {
    for(typename Polynomial<X>::iterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex& a=iter->key();
        X& c=iter->data();
        ++a[j];
        c/=a[j];
    }
    return *this;
}

template<class X>
Polynomial<X>& Polynomial<X>::truncate(smoothness_type d) {
    Polynomial<X> r(this->argument_size());
    for(typename Polynomial<X>::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()<=d && iter->data()!=0) {
            r.append(iter->key(),iter->data());
        }
    }
    this->swap(r);
    return *this;
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

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p) {
    return p; }
//! \relates Polynomial<X> \brief .
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p) {
    Polynomial<X> r(p.argument_size());  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { r[iter->key()]-=iter->data(); } return r; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]+=c; return r; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); r[MultiIndex(p.argument_size())]-=c; return r; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p, const X& c) {
    if(c==0) { return Polynomial<X>(p.argument_size()); }
    Polynomial<X> r(p); typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->data()*=c; } return r; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator/(const Polynomial<X>& p, const X& c) {
    Polynomial<X> r(p); typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=r.begin(); iter!=r.end(); ++iter) { iter->data()/=c; } return r; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator+(const X& c, const Polynomial<X>& p) {
    return p+c; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator-(const X& c, const Polynomial<X>& p) {
    return (-p)+c; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator*(const X& c, const Polynomial<X>& p) {
    return p*c; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator+(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->key()]+=iter->data(); } return r; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator-(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1);  typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=p2.begin(); iter!=p2.end(); ++iter) { r[iter->key()]-=iter->data(); } return r; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> operator*(const Polynomial<X>& p1, const Polynomial<X>& p2) {
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    Polynomial<X> r(p1.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter1=p1.begin(); iter1!=p1.end(); ++iter1) {
        for(Iter iter2=p2.begin(); iter2!=p2.end(); ++iter2) {
            MultiIndex a=iter1->key()+iter2->key();
            r[a]+=iter1->data()*iter2->data(); } } return r; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> sqr(const Polynomial<X>& p) {
    return p*p; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X> pow(const Polynomial<X>& p, unsigned int m) {
    Polynomial<X> r=Polynomial<X>::constant(p.argument_size(),1.0); Polynomial<X> q(p);
    while(m) { if(m%2) { r=r*q; } q=q*q; m/=2; } return r; }

//! \relates Polynomial \brief .
template<class X, class XX> inline Polynomial<X>& operator+=(Polynomial<X>& p, const Polynomial<XX>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->key()]+=iter->data(); } return p; }
//! \relates Polynomial \brief .
template<class X, class XX> inline Polynomial<X>& operator-=(Polynomial<X>& p, const Polynomial<XX>& q) {
    ARIADNE_ASSERT(p.argument_size()==q.argument_size());
    typedef typename Polynomial<X>::const_iterator Iter;
    for(Iter iter=q.begin(); iter!=q.end(); ++iter) { p[iter->key()]-=iter->data(); } return p; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X>& operator+=(Polynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]+=c; return p; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X>& operator-=(Polynomial<X>& p, const X& c) {
    p[MultiIndex(p.argument_size())]-=c; return p; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const X& c) {
    typedef typename Polynomial<X>::iterator Iter;
    if(c==0) { p.clear(); }
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->data()*=c; } return p; }
//! \relates Polynomial \brief .
template<class X> inline Polynomial<X>& operator/=(Polynomial<X>& p, const X& c) {
    typedef typename Polynomial<X>::iterator Iter;
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->data()/=c; } return p; }

//! \relates Polynomial \brief .
template<class X> inline Polynomial<X>& operator*=(Polynomial<X>& p, const Monomial<X>& m) {
    typedef typename Polynomial<X>::iterator Iter;
    if(m.data()==0) { p.clear(); }
    for(Iter iter=p.begin(); iter!=p.end(); ++iter) { iter->key()+=m.key(); iter->data()*=m.data(); } return p; }

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


template<class X>
Polynomial<X>
partial_evaluate(const Polynomial<X>& x, uint k, const X& c)
{
    Polynomial<X> r(x.argument_size()-1);
    MultiIndex ra(r.argument_size());
    if(c==0) {
        for(typename Polynomial<X>::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            MultiIndex::value_type xak=xa[k];
            if(xak==0) {
                const X& xv=xiter->data();
                for(uint i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(uint i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
                r.expansion().append(ra,xv);
            }
        }
    } else if(c==1) {
        Polynomial<X> s(x.argument_size()-1);
        array< Polynomial<X> > p(x.degree()+1,Polynomial<X>(x.argument_size()-1));

        for(typename Polynomial<X>::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const X& xv=xiter->data();
            MultiIndex::value_type xak=xa[k];
            for(uint i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(uint i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }

        r=p[0];
        for(uint i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
    } else {
        Polynomial<X> s(x.argument_size()-1);
        array< Polynomial<X> > p(x.degree()+1,Polynomial<X>(x.argument_size()-1));

        array<X> cpowers(x.degree()+1);
        cpowers[0]=static_cast<X>(1); cpowers[1]=c;
        if(x.degree()>=2) { cpowers[2]=sqr(c); }
        for(uint j=3; j<=x.degree(); ++j) {
            cpowers[j]=cpowers[j-2]*cpowers[2];
        }

        for(typename Polynomial<X>::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const X& xv=xiter->data();
            MultiIndex::value_type xak=xa[k];
            for(uint i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(uint i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }
        for(uint i=1; i!=p.size(); ++i) {
            p[i]*=cpowers[i];
        }

        r=p[0];
        for(uint i=1; i!=p.size(); ++i) {
            r+=p[i];
        }
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
    return Polynomial<X>(embed(before_size,x.expansion(),after_size));
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
            r.append(ar,val*static_cast<X>(ap[j]));
        }
    }
    return r;
}

template<class X>
Vector< Polynomial<X> >
derivative(const Polynomial<X>& p) {
    Vector< Polynomial<X> > r(p.argument_size(), Polynomial<X>(p.argument_size()) );
    for(uint j=0; j!=r.size(); ++j) {
        r[j]=derivative(p,j);
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
truncate(const Polynomial<X>& p, unsigned short d) {
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


template<class X> Vector< Polynomial<X> > operator*(const Polynomial<X>& p, const Vector<Float> e) {
    Vector< Polynomial<X> > r(e.size(),Polynomial<X>(p.argument_size()));
    for(uint i=0; i!=r.size(); ++i) { r[i]=p; r[i]*=X(e[i]); }
    return r;
}

template<class X> 
inline Vector< Polynomial<X> > operator*(const Vector<Float> e, const Polynomial<X>& p) {
  return p*e;
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


/*
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
*/

template<class X>
std::ostream& operator<<(std::ostream& os, const Polynomial<X>& p) {
    bool first_term=true;
    bool zero=true;
    for(typename Polynomial<X>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
        MultiIndex a=iter->key();
        X v=iter->data();
        if(v!=0) {
            zero=false;
            if(v>0 && !first_term) { os<<"+"; }
            first_term=false;
            bool first_factor=true;
            if(v<0) { os<<"-"; }
            if(abs(v)!=1 || a.degree()==0) { os<<abs(v); first_factor=false; }
            for(uint j=0; j!=a.size(); ++j) {
                if(a[j]!=0) {
                    if(first_factor) { first_factor=false; } else { os <<"*"; }
                    os<<"x"<<j; if(a[j]!=1) { os<<"^"<<int(a[j]); } }
            }
        }
    }
    if(zero) { os << "0"; }
    return os;
}




} // namespace Ariadne

#endif /* ARIADNE_POLYNOMIAL_H */
