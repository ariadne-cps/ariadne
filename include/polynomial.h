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
    : public Expansion<X>
{
  public:
    //! \brief The zero polynomial in \a as variables.
    Polynomial(unsigned int as=0u) : Expansion<X>(as) { }
    //! \brief Copy/conversion constructor.
    template<class XX> Polynomial(const Expansion<XX>& e) : Expansion<X>(e) { }
    //! \brief A dense polynomial with coefficients given by a list of doubles.
    Polynomial(unsigned int as, unsigned int deg, double c0, ...);
    //! \brief A sparse polynomial with coefficients given by a list of indices and coefficients.
    Polynomial(unsigned int as, unsigned int nnz, int a00, ...);

};

template<class X>
Polynomial<X>::Polynomial(unsigned int as, unsigned int deg, double c0, ...)
    : Expansion<X>(as)
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
Polynomial<X>::Polynomial(unsigned int as, unsigned int nnz, int a00, ...)
    : Expansion<X>(as)
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


template<class X, class Y> inline
Y evaluate(const Polynomial<X>& p, const Vector<Y>& x)
{
    return evaluate(static_cast< Expansion<X> >(p),x);
}


template<class X, class Y> inline
Vector<Y> evaluate(const Vector< Polynomial<X> >& p, const Vector<Y>& x)
{
    ARIADNE_ASSERT(p.size()>0 && p[0].argument_size()==x.size());
    Y zero = x[0]; zero*=0.0;
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


template<class X> Vector< Polynomial<X> > operator*(const Polynomial<X>& p, const Vector<Float> e) {
    Vector< Polynomial<X> > r(e.size(),Polynomial<X>(p.argument_size()));
    for(uint i=0; i!=r.size(); ++i) { r[i]=p; r[i]*=X(e[i]); }
    return r;
}

}

#endif /* ARIADNE_POLYNOMIAL_H */
