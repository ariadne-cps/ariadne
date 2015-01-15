/***************************************************************************
 *            series.h
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

/*! \file series.h
 *  \brief Taylor series in a single variable.
 */

#ifndef ARIADNE_SERIES_H
#define ARIADNE_SERIES_H

#include "utility/array.h"
#include "numeric/integer.h"
#include "expression/operators.h"

namespace Ariadne {

template<class X>
class PowerSeries
{

  public:
    const X& operator[](DegreeType d) const;

    static PowerSeries<X> rec(const X& x);
    static PowerSeries<X> pow(const X& x, Int n);
    static PowerSeries<X> sqrt(const X& x);
    static PowerSeries<X> exp(const X& x);
    static PowerSeries<X> log(const X& x);

    static PowerSeries<X> sin(const X& x);
    static PowerSeries<X> cos(const X& x);
    static PowerSeries<X> tan(const X& x);
    static PowerSeries<X> asin(const X& x);
    static PowerSeries<X> acos(const X& x);
    static PowerSeries<X> atan(const X& x);
};

template<class X>
class Series
{
  public:
    explicit Series() : _data(1) { }
    explicit Series(Nat d) : _data(d+1) { }
    explicit Series(Nat d, const X& z) : _data(d+1,z) { }
    template<class XX> Series(Nat d, const XX* ptr) : _data(ptr,ptr+(d+1)) { }
    Nat degree() const { return this->_data.size()-1; }
    X& operator[](Nat i) { return this->_data[i]; }
    const X& operator[](Nat i) const { return this->_data[i]; }

    template<class R> Series<X>& operator+=(const Series<R>& y) { _data+=y._data; return *this; }
    template<class R> Series<X>& operator*=(const Series<R>& y) { *this=(*this)*y; return *this; }
    template<class R> Series<X>& operator+=(const R& c) { _data[0]+=c; return *this; }
    template<class R> Series<X>& operator*=(const R& c) { _data*=c; return *this; }

    static Series<X> rec(Nat d, const X& x);
    static Series<X> pow(Nat d, const X& x, Int n);
    static Series<X> sqrt(Nat d, const X& x);
    static Series<X> exp(Nat d, const X& x);
    static Series<X> log(Nat d, const X& x);

    static Series<X> sin(Nat d, const X& x);
    static Series<X> cos(Nat d, const X& x);
    static Series<X> tan(Nat d, const X& x);
    static Series<X> asin(Nat d, const X& x);
    static Series<X> acos(Nat d, const X& x);
    static Series<X> atan(Nat d, const X& x);
  private:
    Array<X> _data;
};

template<class X> inline Series<X> make_series(Rec, Nat d, const X& x) { return Series<X>::rec(d,x); }
template<class X> inline Series<X> make_series(Sqrt, Nat d, const X& x) { return Series<X>::sqrt(d,x); }
template<class X> inline Series<X> make_series(Exp, Nat d, const X& x) { return Series<X>::exp(d,x); }
template<class X> inline Series<X> make_series(Log, Nat d, const X& x) { return Series<X>::log(d,x); }
template<class X> inline Series<X> make_series(Sin, Nat d, const X& x) { return Series<X>::sin(d,x); }
template<class X> inline Series<X> make_series(Cos, Nat d, const X& x) { return Series<X>::cos(d,x); }
template<class X> inline Series<X> make_series(Tan, Nat d, const X& x) { return Series<X>::tan(d,x); }

template<class X> Series<X> operator+(X c, const Series<X>& x) {
    Series<X> r(x); r+=c; return r; }

template<class X> Series<X> operator-(X c, const Series<X>& x) {
    Series<X> r(neg(x)); r+=c; return r; }

template<class X> Series<X> neg(const Series<X>& x) {
    ARIADNE_NOT_IMPLEMENTED; }

template<class X> Series<X> rec(const Series<X>& x) {
    ARIADNE_NOT_IMPLEMENTED; }

template<class X> Series<X> sqr(const Series<X>& x) {
    ARIADNE_NOT_IMPLEMENTED; }

template<class X> Series<X> sqrt(const Series<X>& x) {
    ARIADNE_NOT_IMPLEMENTED; }


template<class X>
Series<X>
Series<X>::rec(Nat d, const X& c)
{
    Series<X> y(d,c*0);
    X mr = (-1)/c;
    for(Nat i=0; i<=y.degree(); ++i) {
        y[i]=-Ariadne::pow(mr,i+1u);
    }
    return y;
}

template<class X>
Series<X>
Series<X>::pow(Nat d, const X& c, Int k)
{

    Nat n=k;
    Series<X> y(d,c*0);
    for(Nat i=0; i<=std::min(Nat(d),n); ++i) {
        Nat j=n-i;
        y[i]=X(bin(n,j))*Ariadne::pow(c,j);
    }
    return y;
}

template<class X>
Series<X>
Series<X>::sqrt(Nat d, const X& c)
{
    Series<X> y(d,c*X(0));
    y[0]=Ariadne::sqrt(c);
    X mhr=(-1/c)/2;
    for(Nat i=1; i<=y.degree(); ++i) {
        // Need to convert Nat to Int to prevent wraparound for 2*1u-3
        y[i]=((2*Int(i)-3)*mhr)/i*y[i-1];
    }
    return y;
}

template<class X>
Series<X>
Series<X>::exp(Nat d, const X& c)
{
    Series<X> y(d,c*0);
    y[0]=Ariadne::exp(c);
    for(Nat i=1; i<=y.degree(); ++i) {
        y[i]=y[i-1]/i;
    }
    return y;
}

template<class X>
Series<X>
Series<X>::log(Nat d, const X& c)
{
    Series<X> y(d,c*0);
    y[0]=Ariadne::log(c);
    X mr=(-1)/c;
    for(Nat i=1; i<=y.degree();++i) {
        y[i]=-Ariadne::pow(mr,i)/i;
    }
    return y;
}

template<class X>
Series<X>
Series<X>::sin(Nat d, const X& c)
{
    Series<X> y(d,c*0);
    y[0]=Ariadne::sin(c);
    if(d>=1) {
        y[1]=Ariadne::cos(c);
        for(Nat i=2; i<=d; ++i) {
            y[i]=-y[i-2]/(i*(i-1));
        }
    }
    return y;
}


template<class X>
Series<X>
Series<X>::cos(Nat d, const X& c)
{
    Series<X> y(d,c*0);
    y[0]=Ariadne::cos(c);
    if(d>=1) {
        y[1]=-Ariadne::sin(c);
        for(Nat i=2; i<=d; ++i) {
            y[i]=-y[i-2]/(i*(i-1));
        }
    }
    return y;
}

template<class X>
Series<X>
Series<X>::tan(Nat d, const X& c)
{
    ARIADNE_NOT_IMPLEMENTED;
}

template<class X>
Series<X>
Series<X>::asin(Nat d, const X& c)
{
    if(d==0) { Series<X> y(d,c*0); y[0]=Ariadne::asin(c); return y; }
    Series<X> y(d-1,c*0); y[0]=c; y[1]=X(1);
    y = Ariadne::rec(Ariadne::sqrt(X(1)-Ariadne::sqr(y)));
    return antiderivative(y,Ariadne::asin(c));
}

template<class X>
Series<X>
Series<X>::acos(Nat d, const X& c)
{
    if(d==0) { Series<X> y(d,c*0); y[0]=Ariadne::acos(c); return y; }
    Series<X> y(d-1,c*0); y[0]=c; y[1]=X(1);
    y = Ariadne::neg(Ariadne::rec(Ariadne::sqrt(X(1)-Ariadne::sqr(y))));
    return antiderivative(y,Ariadne::acos(c));
}

template<class X>
Series<X>
Series<X>::atan(Nat d, const X& c)
{
    if(d==0) { Series<X> y(d,c*0); y[0]=Ariadne::atan(c); return y; }
    Series<X> y(d-1,c*0); y[0]=c; y[1]=X(1);
    y = Ariadne::rec(Ariadne::sqrt(X(1)+Ariadne::sqr(y)));
    return antiderivative(y,Ariadne::atan(c));
}



template<class X>
Series<X>
antiderivative(const Series<X>& x, const X& c)
{
    Nat n=x.degree();
    Series<X> r(n+1,x[0]*0);
    r[0]=c;
    for(Nat i=0; i<=n; ++i) {
        r[i+1]=x[i]/(i+1);
    }
    return r;
}


template<class X>
Series<X>
operator*(const Series<X>& x, const Series<X>& y)
{
    //const Series<X>& x=*this;
    X zero=x[0]*y[0]*0.0;
    Series<X> r(std::min(x.degree(),y.degree()));
    for(Nat k=0; k<=r.degree(); ++k) {
        r[k]=zero;
    }
    for(Nat i=0; i<=r.degree(); ++i) {
        for(Nat j=0; j<=r.degree()-i; ++j) {
            r[i+j]+=bin(i+j,i)*x[i]*y[j];
        }
    }
    return r;
}


template<class X>
OutputStream& operator<<(OutputStream& os, const Series<X>& x)
{
    os << "S";
    for(Nat i=0; i<=x.degree(); ++i) {
        os << (i==0 ? '[' : ',') << x[i];
    }
    return os << "]";
}

} // namespace Ariadne


#endif // ARIADNE_SERIES_H

