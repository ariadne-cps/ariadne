/***************************************************************************
 *            univariate_differential.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file univariate_differential.hpp
 *  \brief Differentials with respect to a single variable.
 */

#ifndef ARIADNE_UNIVARIATE_DIFFERENTIAL_HPP
#define ARIADNE_UNIVARIATE_DIFFERENTIAL_HPP

#include "utility/macros.hpp"
#include "utility/array.hpp"
#include "algebra/series.hpp"
#include "numeric/numeric.hpp"

namespace Ariadne {


template<class X> class Series;
template<class X> class UnivariateDifferential;
template<class T, class X> UnivariateDifferential<X> compose(const Series<T>&, const UnivariateDifferential<X>&);

//! \ingroup DifferentiationModule
//! \brief Arbitrary-order derivatives with respect to a single argument.
template<class X> class UnivariateDifferential
//    : public DispatchTranscendentalAlgebraOperations<UnivariateDifferential<X>,X>
{
    Array<X> _ary;
  public:
    typedef X NumericType;
    typedef X ValueType;
    typedef typename X::Paradigm Paradigm;
    typedef UnivariateDifferential<X> SelfType;
    typedef typename Array<X>::Iterator Iterator;
    typedef typename Array<X>::ConstIterator ConstIterator;

    UnivariateDifferential();
    explicit UnivariateDifferential(DegreeType d);
    explicit UnivariateDifferential(DegreeType d, const NumericType& c);
    UnivariateDifferential(DegreeType d, InitializerList<X> e);
    template<class XX> UnivariateDifferential(DegreeType d, XX const* p);
    UnivariateDifferential(DegreeType d, Series<X> const& s); // explicit

    static SelfType constant(DegreeType d, const NumericType& c);
    static SelfType variable(DegreeType d, const NumericType& c);

    SelfType create_zero() const;
    SelfType create_constant(const NumericType& c) const;
    SelfType create_variable(const NumericType& c) const;

    SizeType argument_size() const;
    DegreeType degree() const;
    X zero_coefficient() const;
    const Array<X>& array() const;
    Expansion<X>& array();
    const X& operator[](SizeType k) const;
    X& operator[](SizeType k);

    SelfType& operator=(const NumericType& c);
    template<class T, EnableIf<IsAssignable<X,T>> =dummy>
        SelfType& operator=(const T& c) { X xc=nul(this->value()); xc=c; return (*this)=xc; }

    SelfType& operator+=(const SelfType& x);
    SelfType& operator-=(const SelfType& x);
    SelfType& operator*=(const SelfType& x);
    SelfType& operator+=(const NumericType& c);
    SelfType& operator*=(const NumericType& c);

    SelfType apply(Rec) const;

    const X& value() const;
    const X& gradient() const;
    const X hessian() const;
    const X& half_hessian() const;

    Void clear();
    OutputStream& write(OutputStream& os) const;

    friend UnivariateDifferential<X> compose(Series<X> const& f, UnivariateDifferential<X> const& dx) {
        return UnivariateDifferential<X>::_compose(f,dx); }
    friend UnivariateDifferential<X> derivative(UnivariateDifferential<X> const& dx) {
        return UnivariateDifferential<X>::_derivative(dx); }
    friend UnivariateDifferential<X> antiderivative(UnivariateDifferential<X> const& dx) {
        return UnivariateDifferential<X>::_antiderivative(dx); }
    friend UnivariateDifferential<X> antiderivative(UnivariateDifferential<X> const& dx, X const& c) {
        return UnivariateDifferential<X>::_antiderivative(dx,c); }

    friend OutputStream& operator<<(OutputStream& os, UnivariateDifferential<X> const& dx) {
        return dx.write(os); }
  public:
    static Differential<X> _compose(Series<X> const& f, Differential<X> const& dx);
  private:
    static UnivariateDifferential<X> _compose(Series<X> const& f, UnivariateDifferential<X> const& dx);
    static UnivariateDifferential<X> _derivative(UnivariateDifferential<X> const& dx);
    static UnivariateDifferential<X> _antiderivative(UnivariateDifferential<X> const& dx);
    static UnivariateDifferential<X> _antiderivative(UnivariateDifferential<X> const& dx, X const& c);
};


template<class X> inline const X& UnivariateDifferential<X>::value() const { return _ary[0]; }
template<class X> inline const X& UnivariateDifferential<X>::gradient() const { assert(this->degree()>=1); return _ary[1]; }
template<class X> inline const X& UnivariateDifferential<X>::half_hessian() const { assert(this->degree()>=2); return _ary[2]; }
template<class X> inline const X UnivariateDifferential<X>::hessian() const { assert(this->degree()>=2); return _ary[2]*2; }

template<class X> template<class XX> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, XX const* ptr)
    : _ary(d+1,X(0))
{
    std::copy(ptr,ptr+d+1,_ary.begin());
}



template<class X> UnivariateDifferential<X> create_zero(const UnivariateDifferential<X>& c) {
    return UnivariateDifferential<X>(c.degree(),create_zero(c[0])); }

template<class X> UnivariateDifferential<X>& operator+=(UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    for(Nat i=0; i<=x.degree(); ++i) { x[i]+=y[i]; } return x; }
template<class X> UnivariateDifferential<X>& operator-=(UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    for(Nat i=0; i<=x.degree(); ++i) { x[i]-=y[i]; } return x; }
template<class X> UnivariateDifferential<X>& operator*=(UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    x=x*y; return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>&>
operator+=(UnivariateDifferential<X>& x, const Y& c) {
    x[0]+=c; return x; }
template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>&>
operator-=(UnivariateDifferential<X>& x, const Y& c) {
    x[0]-=c; return x; }
template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>&>
operator*=(UnivariateDifferential<X>& x, const Y& c) {
    for(Nat i=0; i<=x.degree(); ++i) { x[i]*=c; } return x; }
template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>&>
operator/=(UnivariateDifferential<X>& x, const Y& c) {
    for(Nat i=0; i<=x.degree(); ++i) { x[i]/=c; } return x; }

template<class X> UnivariateDifferential<X> operator-(const UnivariateDifferential<X>& x) {
    UnivariateDifferential<X> r(x); for(SizeType i=0; i<=r.degree(); ++i) { r[i]=-r[i]; } return r; }

template<class X> UnivariateDifferential<X> operator+(X c, const UnivariateDifferential<X>& x) {
    UnivariateDifferential<X> r(x); r+=c; return r; }

template<class X> UnivariateDifferential<X> operator-(X c, const UnivariateDifferential<X>& x) {
    UnivariateDifferential<X> r(neg(x)); r+=c; return r; }

template<class X> UnivariateDifferential<X> operator+(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    UnivariateDifferential<X> r(x1); r+=x2; return r; }

template<class X> UnivariateDifferential<X> operator-(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    UnivariateDifferential<X> r(x1); r-=x2; return r; }

template<class X> UnivariateDifferential<X> operator*(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    UnivariateDifferential<X> r(std::min(x1.degree(),x2.degree()),create_zero(x1[0]*x2[0]));
    for(Nat i1=0; i1<=r.degree(); ++i1) {
        for(Nat i2=0; i2<=r.degree()-i1; ++i2) {
            r[i1+i2]+=x1[i1]*x2[i2];
        }
    }
    return r;
}

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator+(UnivariateDifferential<X> x, const Y& c) {
    x[0]+=c; return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator-(UnivariateDifferential<X> x, const Y& c) {
    x[0]-=c; return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator*(UnivariateDifferential<X> x, const Y& c) {
    x*=c; return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator/(UnivariateDifferential<X> x, const Y& c) {
    x*=rec(c); return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator+(const Y& c, UnivariateDifferential<X> x) {
    x+=c; return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator-(const Y& c, UnivariateDifferential<X> x) {
    return c+neg(x); }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator*(const Y& c, UnivariateDifferential<X> x) {
    x*=c; return x; }

template<class X, class Y>
EnableIfNumericType<Y,UnivariateDifferential<X>>
operator/(const Y& c, UnivariateDifferential<X> x) {
    return c*rec(x); }

template<class X> UnivariateDifferential<X> operator/(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    return x1*rec(x2);
}

template<class X> UnivariateDifferential<X> neg(const UnivariateDifferential<X>& dx) {
    UnivariateDifferential<X> r(dx); for(SizeType i=0; i<=r.degree(); ++i) { r[i]=-r[i]; } return r; }

template<class X> UnivariateDifferential<X> rec(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::rec(dx[0]),dx); }

template<class X> UnivariateDifferential<X> sqr(const UnivariateDifferential<X>& dx) {
    return dx*dx; }

template<class X> UnivariateDifferential<X> sqrt(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::sqrt(dx[0]),dx); }

template<class X> UnivariateDifferential<X> pow(const UnivariateDifferential<X>& dx, Int n) {
    return compose(Series<X>::pow(dx[0],n),dx); }

template<class X> UnivariateDifferential<X> exp(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::exp(dx[0]),dx); }

template<class X> UnivariateDifferential<X> log(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::log(dx[0]),dx); }

template<class X> UnivariateDifferential<X> sin(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::sin(dx[0]),dx); }

template<class X> UnivariateDifferential<X> cos(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::cos(dx[0]),dx); }

template<class X> UnivariateDifferential<X> tan(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::tan(dx[0]),dx); }

template<class X> UnivariateDifferential<X> atan(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::atan(dx[0]),dx); }

} // namespace Ariadne


#endif // ARIADNE_UNIVARIATE_DIFFERENTIAL_HPP

