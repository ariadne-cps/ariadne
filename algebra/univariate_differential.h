/***************************************************************************
 *            univariate_differential.h
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

/*! \file univariate_differential.h
 *  \brief Differentials with respect to a single variable.
 */

#ifndef ARIADNE_UNIVARIATE_DIFFERENTIAL_H
#define ARIADNE_UNIVARIATE_DIFFERENTIAL_H

#include "utility/macros.h"
#include "utility/array.h"
#include "algebra/series.h"
#include "numeric/numeric.h"

namespace Ariadne {


template<class X> class Series;
template<class X> class UnivariateDifferential;
template<class T, class X> UnivariateDifferential<X> compose(const Series<T>&, const UnivariateDifferential<X>&);

template<class X>
class UnivariateDifferential
{
    template<class T, class Y> friend UnivariateDifferential<Y> compose(const Series<T>&,const UnivariateDifferential<Y>&);
  public:
    explicit UnivariateDifferential() : _data(1) { }
    explicit UnivariateDifferential(uint d) : _data(d+1) { }
    explicit UnivariateDifferential(uint d, const X& x) : _data(d+1,x) { }
    template<class XX> UnivariateDifferential(uint d, const XX* ptr) : _data(ptr,ptr+(d+1)) { }

    static UnivariateDifferential<X> constant(uint d, const X& c) {
        UnivariateDifferential<X> result(d,create_zero(c)); result._data[0]=c; return result; }
    static UnivariateDifferential<X> variable(uint d, const X& x) {
        UnivariateDifferential<X> result(d,create_zero(x)); result._data[0]=x; result._data[1]=1; return result; }

    UnivariateDifferential<X>& operator=(const X& c) {
        X z=create_zero(_data[0]); _data[0]=c; for(uint i=1; i!=_data.size(); ++i) { _data[i]=z; } return *this; }
    uint degree() const { return this->_data.size()-1; }
    const X& value() const { return this->_data[0]; }
    const X& gradient() const { return this->_data[1]; }
    const X& operator[](uint i) const { return this->_data[i]; }
    X& operator[](uint i) { return this->_data[i]; }
  private:
    Array<X> _data;
};

template<class X> UnivariateDifferential<X> create_zero(const UnivariateDifferential<X>& c) {
    return UnivariateDifferential<X>(c.degree(),create_zero(c[0])); }

template<class X> UnivariateDifferential<X>& operator+=(UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    for(uint i=0; i<=x.degree(); ++i) { x[i]+=y[i]; } return x; }
template<class X> UnivariateDifferential<X>& operator-=(UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    for(uint i=0; i<=x.degree(); ++i) { x[i]-=y[i]; } return x; }
template<class X> UnivariateDifferential<X>& operator*=(UnivariateDifferential<X>& x, const UnivariateDifferential<X>& y) {
    x=x*y; return x; }

template<class X, class Y>
EnableIfNumeric<Y,UnivariateDifferential<X>&>
operator+=(UnivariateDifferential<X>& x, const Y& c) {
    x[0]+=c; return x; }
template<class X, class Y>
EnableIfNumeric<Y,UnivariateDifferential<X>&>
operator*=(UnivariateDifferential<X>& x, const Y& c) {
    for(uint i=0; i<=x.degree(); ++i) { x[i]*=c; } return x; }

template<class X> UnivariateDifferential<X> operator+(X c, const UnivariateDifferential<X>& x) {
    Series<X> r(x); r+=c; return r; }

template<class X> UnivariateDifferential<X> operator-(X c, const UnivariateDifferential<X>& x) {
    UnivariateDifferential<X> r(neg(x)); r+=c; return r; }

template<class X> UnivariateDifferential<X> operator+(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    UnivariateDifferential<X> r(x1); r+=x2; return r; }

template<class X> UnivariateDifferential<X> operator-(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    UnivariateDifferential<X> r(x1); r-=x2; return r; }

template<class X> UnivariateDifferential<X> operator*(const UnivariateDifferential<X>& x1, const UnivariateDifferential<X>& x2) {
    UnivariateDifferential<X> r(min(x1.degree(),x2.degree()),create_zero(x1[0]*x2[0]));
    for(uint i1=0; i1<=r.degree(); ++i1) {
        for(uint i2=0; i2<=r.degree()-i1; ++i2) {
            r[i1+i2]+=x1[i1]*x2[i2];
        }
    }
    return r;
}

template<class X> UnivariateDifferential<X> neg(const UnivariateDifferential<X>& dx) {
    ARIADNE_NOT_IMPLEMENTED; }

template<class X> UnivariateDifferential<X> rec(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::rec(dx.degree(),dx[0]),dx); }

template<class X> UnivariateDifferential<X> sqr(const UnivariateDifferential<X>& dx) {
    return dx*dx; }

template<class X> UnivariateDifferential<X> sqrt(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::sqrt(dx.degree(),dx[0]),dx); }

template<class X> UnivariateDifferential<X> pow(const UnivariateDifferential<X>& dx, int n) {
    return compose(Series<X>::pow(dx.degree(),dx[0],n),dx); }

template<class X> UnivariateDifferential<X> exp(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::exp(dx.degree(),dx[0]),dx); }

template<class X> UnivariateDifferential<X> log(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::log(dx.degree(),dx[0]),dx); }

template<class X> UnivariateDifferential<X> sin(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::sin(dx.degree(),dx[0]),dx); }

template<class X> UnivariateDifferential<X> cos(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::cos(dx.degree(),dx[0]),dx); }

template<class X> UnivariateDifferential<X> tan(const UnivariateDifferential<X>& dx) {
    return compose(Series<X>::tan(dx.degree(),dx[0]),dx); }


template<class X, class Y>
UnivariateDifferential<Y> compose(const Series<X>& x, const UnivariateDifferential<Y>& y)
{
    uint d=std::min(x.degree(),y.degree());

    Y y0 = y[0];
    Y z = create_zero(y0);
    const_cast<Y&>(y[0]) = z;
    UnivariateDifferential<Y> r(d, z);
    r[0]=x[d];
    for(uint n=1; n<=d; ++n) {
        r=r*y;
        r[0]+=x[d-n];
    }
    const_cast<Y&>(y[0])=y0;
    return r;
}

template<class X>
UnivariateDifferential<X>
derivative(const UnivariateDifferential<X>& x)
{
    uint n=x.degree();
    UnivariateDifferential<X> r(min(n,1u)-1u);
    if(n==0) { r[0]=x[0]*0; }
    for(uint i=0; i<n; ++i) {
        r[i]=(i+1)*x[i+1];
    }
    return r;
}

template<class X>
UnivariateDifferential<X>
antiderivative(const UnivariateDifferential<X>& x)
{
    uint n=x.degree();
    UnivariateDifferential<X> r(n+1,create_zero(x[0]));
    for(uint i=0; i<=n; ++i) {
        r[i+1]=x[i]/(i+1);
    }
    return r;
}

template<class X>
UnivariateDifferential<X>
antiderivative(const UnivariateDifferential<X>& x, const X& c)
{
    uint n=x.degree();
    UnivariateDifferential<X> r(n+1,create_zero(x[0]));
    r[0]=c;
    for(uint i=0; i<=n; ++i) {
        r[i+1]=x[i]/(i+1);
    }
    return r;
}



template<class X>
std::ostream& operator<<(std::ostream& os, const UnivariateDifferential<X>& x)
{
    os << "D<"<<x.degree()<<">";
    for(uint i=0; i<=x.degree(); ++i) {
        os << (i==0 ? '[' : ',') << x[i];
    }
    return os << "]";
}

} // namespace Ariadne


#endif // ARIADNE_UNIVARIATE_DIFFERENTIAL_H

