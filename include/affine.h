/***************************************************************************
 *            affine.h
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

/*! \file affine.h
 *  \brief Affine scalar and vector functions
 */

#ifndef ARIADNE_AFFINE_H
#define ARIADNE_AFFINE_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "macros.h"
#include "pointer.h"

#include "vector.h"
#include "matrix.h"

namespace Ariadne {

template<class X> class Affine;
template<class X> bool operator==(const Affine<X>&, const Affine<X>&);
template<class X> Affine<X> operator-(const Affine<X>&);
template<class X> Affine<X> operator+(const Affine<X>&, const Affine<X>&);
template<class X> Affine<X> operator-(const Affine<X>&, const Affine<X>&);
template<class X> Affine<X> operator+(const X&, const Affine<X>&);
template<class X> Affine<X> operator+(const Affine<X>&, const X&);
template<class X> Affine<X> operator-(const X&, const Affine<X>&);
template<class X> Affine<X> operator-(const Affine<X>&, const X&);
template<class X> Affine<X> operator*(const X&, const Affine<X>&);
template<class X> Affine<X> operator*(const Affine<X>&, const X&);
template<class X> Affine<X> operator/(const Affine<X>&, const X&);
template<class X> X derivative(const Affine<X>&, uint);

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
template<class X>
class Affine
{
  public:
    explicit Affine() : _c(), _g() { }
    explicit Affine(uint n) : _c(0), _g(n) { }
    explicit Affine(const Vector<X>& g, const X& c) : _c(c), _g(g) { }
    explicit Affine(uint as, double c, double g0, ...) : _c(static_cast<X>(c)), _g(as) {
        ARIADNE_ASSERT(as>=2); _g[0]=static_cast<X>(g0); va_list args; va_start(args,g0);
        for(uint i=1; i!=as; ++i) { _g[i]=static_cast<X>(va_arg(args,double)); } }
    template<class XX> explicit Affine(const Affine<XX>& aff)
        : _c(aff.b()), _g(aff.a()) { }

    Affine<X>& operator=(const X& c) {
        this->_c=c; for(uint i=0; i!=this->_g.size(); ++i) { this->_g[i]=static_cast<X>(0); } return *this; }
    static Affine<X> constant(uint n, X c) {
        return Affine<X>(Vector<X>(n),c); }
    static Affine<X> variable(uint n, uint j) {
        return Affine<X>(Vector<X>::unit(n,j),X(0)); }
    static Vector< Affine<X> > variables(uint n) {
        Vector< Affine<X> > r(n,Affine<X>(n)); for(uint i=0; i!=n; ++i) { r[i]._g[i]=static_cast<X>(1); } return r; }

    const X& operator[](uint i) const { return this->_g[i]; }
    X& operator[](uint i) { return this->_g[i]; }


    const Vector<X>& a() const { return this->_g; }
    const X& b() const { return this->_c; }

    const Vector<X>& gradient() const { return this->_g; }
    const X& gradient(uint i) const { return this->_g[i]; }
    const X& value() const { return this->_c; }

    uint argument_size() const { return this->_g.size(); }

    template<class Y> Y evaluate(const Vector<Y>& x) const {
        Y r=x.zero_element(); for(uint j=0; j!=this->_g.size(); ++j) { r+=this->_g[j]*x[j]; } return r; }

    const X& derivative(uint j) const { return this->_g[j]; }
  private:
    friend bool operator==<>(const Affine<X>&, const Affine<X>&);
    friend Affine<X> operator-<>(const Affine<X>&);
    friend Affine<X> operator+<>(const Affine<X>&, const Affine<X>&);
    friend Affine<X> operator-<>(const Affine<X>&, const Affine<X>&);
    friend Affine<X> operator+<>(const X&, const Affine<X>&);
    friend Affine<X> operator+<>(const Affine<X>&, const X&);
    friend Affine<X> operator-<>(const X&, const Affine<X>&);
    friend Affine<X> operator-<>(const Affine<X>&, const X&);
    friend Affine<X> operator*<>(const X&, const Affine<X>&);
    friend Affine<X> operator*<>(const Affine<X>&, const X&);
    friend Affine<X> operator/<>(const Affine<X>&, const X&);
  private:
    X _c;
    Vector<X> _g;
};

//! \relates Affine
//! \brief Test equality of two affine expressions.
template<class X> inline bool operator==(const Affine<X>& f1, const Affine<X>& f2) {
    return f1._c==f2._c && f1._g == f2._g; }
//! \relates Affine
//! \brief Negation of an affine expression.
template<class X> inline Affine<X> operator-(const Affine<X>& f) {
    return Affine<X>(Vector<X>(-f._g),-f._c); }
//! \relates Affine
//! \brief Addition of two affine expressions.
template<class X> inline Affine<X> operator+(const Affine<X>& f1, const Affine<X>& f2) {
    return Affine<X>(Vector<X>(f1._g+f2._g),f1._c+f2._c); }
//! \relates Affine
//! \brief Subtraction of two affine expressions.
template<class X> inline Affine<X> operator-(const Affine<X>& f1, const Affine<X>& f2) {
    return Affine<X>(Vector<X>(f1._g-f2._g),f1._c-f2._c); }
//! \relates Affine
//! \brief Addition of a constant to an affine expression.
template<class X> inline Affine<X> operator+(const Affine<X>& f1, const X& c2) {
    return Affine<X>(Vector<X>(f1._g),f1._c+c2); }
//! \relates Affine
//! \brief Addition of a constant to an affine expression.
template<class X> inline Affine<X> operator+(const X& c1, const Affine<X>& f2) {
    return Affine<X>(Vector<X>(f2._g),c1+f2._c); }
//! \relates Affine
//! \brief Subtraction of a constant to an affine expression.
template<class X> inline Affine<X> operator-(const Affine<X>& f1, const X& c2) {
    return Affine<X>(Vector<X>(f1._g),f1._c-c2); }
//! \relates Affine
//! \brief Subtraction of an affine expression from a constant.
template<class X> inline Affine<X> operator-(const X& c1, const Affine<X>& f2) {
    return Affine<X>(Vector<X>(-f2._g),c1-f2._c); }
//! \relates Affine
//! \brief Scalar multiplication of an affine expression.
template<class X> inline Affine<X> operator*(const X& c, const Affine<X>& f) {
    return Affine<X>(Vector<X>(c*f._g),c*f._c); }
//! \relates Affine
//! \brief Scalar multiplication of an affine expression.
template<class X> inline Affine<X> operator*(const Affine<X>& f, const X& c) { return c*f; }
//! \relates Affine
//! \brief Scalar division of an affine expression.
template<class X> inline Affine<X> operator/(const Affine<X>& f, const X& c) { return (1/c)*f; }
//! \relates Affine
//! \brief The derivative of an affine expression gives a constant.
template<class X> inline X derivative(const Affine<X>& f, uint k) { return f.derivative(k); }

template<class X> inline Affine<X> operator+(const Affine<X>& f, double c) {
    return f+static_cast<X>(c); }
template<class X> inline Affine<X> operator+(double c, const Affine<X>& f) {
    return static_cast<X>(c)+f; }
template<class X> inline Affine<X> operator-(const Affine<X>& f, double c) {
    return f-static_cast<X>(c); }
template<class X> inline Affine<X> operator-(double c, const Affine<X>& f) {
    return static_cast<X>(c)-f; }
template<class X> inline Affine<X> operator*(const Affine<X>& f, double c) {
    return f*static_cast<X>(c); }
template<class X> inline Affine<X> operator*(double c, const Affine<X>& f) {
    return static_cast<X>(c)*f; }
template<class X> inline Affine<X> operator/(const Affine<X>& f, double c) {
    return f/static_cast<X>(c); }

template<class X> std::ostream& operator<<(std::ostream& os, const Affine<X>& f) {
    bool zero=true;
    if(f.b()!=0) { os<<f.b(); zero=false; }
    for(uint j=0; j!=f.argument_size(); ++j) {
        if(f.a()[j]!=0) {
            if(f.a()[j]>0) { if(!zero) { os<<"+"; } } else { os<<"-"; }
            if(abs(f.a()[j])!=1) { os<<abs(f.a()[j])<<"*"; }
            //ss<<char('x'+j);
            os<<"x"<<j;
            zero=false;
        }
    }
    if(zero) { os << "0"; }
    return os;
}





} // namespace Ariadne

#endif /* ARIADNE_AFFINE_H */
