/***************************************************************************
 *            affine.hpp
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

/*! \file affine.hpp
 *  \brief Affine scalar and vector functions
 */

#ifndef ARIADNE_AFFINE_HPP
#define ARIADNE_AFFINE_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../utility/declarations.hpp"

#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"

#include "../algebra/covector.hpp"
#include "../algebra/matrix.hpp"

namespace Ariadne {

template<class PR> class FloatFactory;

template<class X> class Affine;
typedef Affine<ApproximateNumericType> ApproximateAffine;
typedef Affine<ValidatedNumericType> ValidatedAffine;
typedef Affine<EffectiveNumericType> EffectiveAffine;

/*
template<class X> struct DispatchAlgebraOperations<Affine<X>,X> {
    friend Bool operator==(const Affine<X>& af1, const Affine<X>& af2) { return Operators<Affine<X>>::_eq(af1,af2); }
    friend Affine<X> operator+(const Affine<X>& af) { return Operators<Affine<X>>::_pos(af); }
    friend Affine<X> operator-(const Affine<X>& af) { return Operators<Affine<X>>::_neg(af); }
    friend Affine<X> operator+(const Affine<X>& af1, const Affine<X>& af2) { return Operators<Affine<X>>::_add(af1,af2); }
    friend Affine<X> operator-(const Affine<X>& af1, const Affine<X>& af2) { return Operators<Affine<X>>::_sub(af1,af2); }
    friend Affine<X> operator+(const X& x1, const Affine<X>& af2) { return Operators<Affine<X>>::_add(x1,af2); }
    friend Affine<X> operator+(const Affine<X>& af1, const X& x2) { return Operators<Affine<X>>::_add(af1,x2); }
    friend Affine<X> operator-(const X& x1, const Affine<X>& af2) { return Operators<Affine<X>>::_sub(x1,af2); }
    friend Affine<X> operator-(const Affine<X>& af1, const X& x2) { return Operators<Affine<X>>::_sub(af1,x2); }
    friend Affine<X> operator*(const X& x1, const Affine<X>& af2) { return Operators<Affine<X>>::_mul(x1,af2); }
    friend Affine<X> operator*(const Affine<X>& af1, const X& x2) { return Operators<Affine<X>>::_mul(af1,x2); }
    friend Affine<X> operator/(const Affine<X>& af1, const X& x2) { return Operators<Affine<X>>::_div(af1,x2); }
};
*/

template<class A, class X=typename A::NumericType> struct ProvideAlgebraOperations;

template<class X> struct ProvideAlgebraOperations<Affine<X>,X> {
    //! \relates Affine
    //! \brief Test equality of two affine expressions.
    friend inline Bool operator==(const Affine<X>& f1, const Affine<X>& f2) {
        return f1._c==f2._c && f1._g == f2._g; }
    //! \relates Affine
    //! \brief Negation of an affine expression.
    friend inline Affine<X> operator-(const Affine<X>& f) {
        return Affine<X>(Covector<X>(-f._g),-f._c); }
    //! \relates Affine
    //! \brief Addition of two affine expressions.
    friend inline Affine<X> operator+(const Affine<X>& f1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(f1._g+f2._g),f1._c+f2._c); }
    //! \relates Affine
    //! \brief Subtraction of two affine expressions.
    friend inline Affine<X> operator-(const Affine<X>& f1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(f1._g-f2._g),f1._c-f2._c); }
    //! \relates Affine
    //! \brief Addition of a constant to an affine expression.
    friend inline Affine<X> operator+(const Affine<X>& f1, const X& c2) {
        return Affine<X>(Covector<X>(f1._g),f1._c+c2); }
    //! \relates Affine
    //! \brief Addition of a constant to an affine expression.
    friend inline Affine<X> operator+(const X& c1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(f2._g),c1+f2._c); }
    //! \relates Affine
    //! \brief Subtraction of a constant to an affine expression.
    friend inline Affine<X> operator-(const Affine<X>& f1, const X& c2) {
        return Affine<X>(Covector<X>(f1._g),f1._c-c2); }
    //! \relates Affine
    //! \brief Subtraction of an affine expression from a constant.
    friend inline Affine<X> operator-(const X& c1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(-f2._g),c1-f2._c); }
    //! \relates Affine
    //! \brief Scalar multiplication of an affine expression.
    friend inline Affine<X> operator*(const X& c, const Affine<X>& f) {
        return Affine<X>(Covector<X>(c*f._g),c*f._c); }
    //! \relates Affine
    //! \brief Scalar multiplication of an affine expression.
    friend inline Affine<X> operator*(const Affine<X>& f, const X& c) { return c*f; }
    //! \relates Affine
    //! \brief Scalar division of an affine expression.
    friend inline Affine<X> operator/(const Affine<X>& f, const X& c) { return (1/c)*f; }
    //! \relates Affine
    //! \brief The derivative of an affine expression gives a constant.
    friend inline X derivative(const Affine<X>& f, Nat k) { return f.derivative(k); }

    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator+(Affine<X> const& x, Y const& y) { return x+factory(x).create(y); }
    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator-(Affine<X> const& x, Y const& y) { return x-factory(x).create(y); }
    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator*(Affine<X> const& x, Y const& y) { return x*factory(x).create(y); }
    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator/(Affine<X> const& x, Y const& y) { return x/factory(x).create(y); }
    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator+(Y const& y, Affine<X> const& x) { return factory(x).create(y)+x; }
    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator-(Y const& y, Affine<X> const& x) { return factory(x).create(y)-x; }
    template<class Y, EnableIf<IsGenericNumericType<Y>> =dummy>
        friend decltype(auto) operator*(Y const& y, Affine<X> const& x) { return factory(x).create(y)*x; }
};

template<class X> FloatFactory<PrecisionType<X>> factory(Affine<X> const& a) {
    return FloatFactory<PrecisionType<X>>(a.value().precision()); }

/*
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator+(Affine<X> const& x, Y const& y) { return x+factory(x).create(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator-(Affine<X> const& x, Y const& y) { return x-factory(x).create(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator*(Affine<X> const& x, Y const& y) { return x*factory(x).create(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator/(Affine<X> const& x, Y const& y) { return x/factory(x).create(y); }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator+(Y const& y, Affine<X> const& x) { return factory(x).create(y)+y; }
template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator*(Y const& y, Affine<X> const& x) { return factory(x).create(y)*y; }
*/

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
template<class X>
class Affine
    : public ProvideAlgebraOperations<Affine<X>,X>
{
  public:
    typedef X NumericType;
  public:
    explicit Affine() : _c(), _g() { }
    explicit Affine(Nat n) : _c(0), _g(n) { }
    explicit Affine(const Covector<X>& g, const X& c) : _c(c), _g(g) { }
    explicit Affine(X c, InitializerList<X> g) : _c(c), _g(g) { }
    template<class XX> explicit Affine(const Affine<XX>& aff)
        : _c(aff.b()), _g(aff.a()) { }

    Affine<X>& operator=(const X& c) {
        this->_c=c; for(SizeType i=0; i!=this->_g.size(); ++i) { this->_g[i]=static_cast<X>(0); } return *this; }
    static Affine<X> constant(SizeType n, X c) {
        return Affine<X>(Covector<X>(n),c); }
    static Affine<X> coordinate(SizeType n, SizeType j) {
        return Affine<X>(Covector<X>::unit(n,j),X(0)); }
    static Vector< Affine<X> > coordinates(SizeType n) {
        Vector< Affine<X> > r(n,Affine<X>(n)); for(SizeType i=0; i!=n; ++i) { r[i]._g[i]=1; } return r; }
    static Affine<X> variable(SizeType n, SizeType j) { return coordinate(n,j); }
    static Vector< Affine<X> > variables(SizeType n) { return coordinates(n); }

    const X& operator[](SizeType i) const { return this->_g[i]; }
    X& operator[](Nat i) { return this->_g[i]; }


    const Covector<X>& a() const { return this->_g; }
    const X& b() const { return this->_c; }

    const Covector<X>& gradient() const { return this->_g; }
    const X& gradient(SizeType i) const { return this->_g[i]; }
    const X& value() const { return this->_c; }

    Void resize(SizeType n) { return this->_g.resize(n); }
    SizeType argument_size() const { return this->_g.size(); }

    template<class Y> Y evaluate(const Vector<Y>& x) const;

    const X& derivative(SizeType j) const { return this->_g[j]; }
  private: public:
    X _c;
    Covector<X> _g;
};

template<class X> template<class Y> Y Affine<X>::evaluate(const Vector<Y>& x) const {
    Y r=x.zero_element(); for(SizeType j=0; j!=this->_g.size(); ++j) { r+=this->_g[j]*x[j]; } return r;
}

template<class X> OutputStream& operator<<(OutputStream& os, const Affine<X>& f) {
    os<<f.b();
    for(SizeType j=0; j!=f.argument_size(); ++j) {
        os<<"+" << "(" << f.a()[j] << ")*x" << j;
    }
    return os;
}





} // namespace Ariadne

#endif /* ARIADNE_AFFINE_HPP */
