/***************************************************************************
 *            fixed_univariate_differential.h
 *
 *  Copyright 2013  Pieter Collins
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

/*! \file fixed_univariate_differential.h
 *  \brief First- and second-order univariate differentials.
 */

#ifndef ARIADNE_FIXED_UNIVARIATE_DIFFERENTIAL_H
#define ARIADNE_FIXED_UNIVARIATE_DIFFERENTIAL_H

#include "utility/macros.h"

namespace Ariadne {

//! \ingroup DifferentiationModule
//! \brief A class representing the value and derivative of a scalar quantity
//! depending on a single argument.
template<class X>
class UnivariateFirstDifferential
{
  public:
    X _value;
    X _gradient;
    static const X _zero;
  public:
    //! \brief The type of used to represent numbers in the coefficient.
    typedef typename X::NumericType NumericType;
    //! \brief The type of used to represent the coefficients.
    typedef X ValueType;

    //! \brief Constructs a differential in \a as variables.
    explicit UnivariateFirstDifferential() : _value(_zero), _gradient(_zero) { }

    //! \brief Constructs a differential with degree zero in \a as variables, initialising to the zero value z.
    explicit UnivariateFirstDifferential(const X& z) : _value(z), _gradient(z) { }

    //! \brief Constructs a differential with degree zero in \a as variables. (Deprecated)
    explicit UnivariateFirstDifferential(const X& v, const X& g) : _value(v), _gradient(g) { }

    //! \brief Conversion constructor from a different numerical type.
    template<class XX> UnivariateFirstDifferential(const UnivariateFirstDifferential<XX>& x)
        : _value(x._value), _gradient(x._gradient) { }

    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    UnivariateFirstDifferential<X>& operator=(const X& c) { _value=c; _gradient=_zero; return *this; }

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static UnivariateFirstDifferential<X> constant(const X& c) {
        UnivariateFirstDifferential<X> r; r._value=c; return r; }
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static UnivariateFirstDifferential<X> variable(const X& v) {
        UnivariateFirstDifferential<X> r; r._value=v; r._gradient=1; return r; }

    //! \brief Equality operator.
    Bool operator==(const UnivariateFirstDifferential<X>& other) const {
        return this->_value==other._value && this->_gradient==other._gradient; }

    //! \brief Inequality operator.
    Bool operator!=(const UnivariateFirstDifferential<X>& other) const {
        return !(*this==other); }

    //! \brief The number of independent variables.
    Nat argument_size() const { return 1u; }
    //! \brief The maximum degree of the stored terms.
    Nat degree() const { return 1u; }
    //! \brief The value of the differential i.e. the coefficient of \f$1\f$.
    const X& value() const { return this->_value; }
    //! \brief The first derivative.
    const X& gradient() const { return this->_gradient; }

    //! \brief Set the coefficient of the constant term to \a c.
    Void set_value(const X& c) { this->_value=c; }

    //! \brief Set all coefficients to zero.
    Void clear() { *this=_zero; }

};

//template<class X>
//const X UnivariateFirstDifferential<X>::_zero=X(0);


template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X>&>
operator+=(UnivariateFirstDifferential<X>& x, const R& c)
{
    x._value+=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X>&>
operator-=(UnivariateFirstDifferential<X>& x, const R& c)
{
    x._value-=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X>&>
operator*=(UnivariateFirstDifferential<X>& x, const R& c)
{
    x._value*=static_cast<X>(c);
    x._gradient*=static_cast<X>(c);
    return x;
}


template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X>&>
operator/=(UnivariateFirstDifferential<X>& x, const R& c)
{
    x._value/=static_cast<X>(c);
    x._gradient/=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator+(const UnivariateFirstDifferential<X>& x, const R& c)
{
    UnivariateFirstDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator+(const R& c, const UnivariateFirstDifferential<X>& x)
{
    UnivariateFirstDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator-(const UnivariateFirstDifferential<X>& x, const R& c)
{
    UnivariateFirstDifferential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator-(const R& c, const UnivariateFirstDifferential<X>& x)
{
    UnivariateFirstDifferential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator*(const UnivariateFirstDifferential<X>& x, const R& c)
{
    UnivariateFirstDifferential<X> r(x); r*=X(c); return r;
}

template<class X>
UnivariateFirstDifferential<X>
operator*(const UnivariateFirstDifferential<X>& x, const Int& c)
{
    UnivariateFirstDifferential<X> r(x); r*=X(0); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator*(const R& c, const UnivariateFirstDifferential<X>& x)
{
    UnivariateFirstDifferential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator/(const UnivariateFirstDifferential<X>& x, const R& c)
{
    UnivariateFirstDifferential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateFirstDifferential<X> >
operator/(const R& c, const UnivariateFirstDifferential<X>& x)
{
    UnivariateFirstDifferential<X> r(rec(x)); r*=X(c); return r;
}


template<class X>
UnivariateFirstDifferential<X>& operator+=(UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    x._value += y._value;
    x._gradient += y._gradient;
    return x;
}

template<class X>
UnivariateFirstDifferential<X>& operator-=(UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    x._value -= y._value;
    x._gradient -= y._gradient;
    return x;
}

template<class X>
UnivariateFirstDifferential<X>& operator*=(UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    x._gradient *= y._value;
    x._gradient += x._value * y._gradient;
    x._value *= y._value;
    return x;
}

template<class X>
UnivariateFirstDifferential<X>& operator/=(UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    x._value /= y._value;
    x._gradient -= x._value * y._gradient;
    x._gradient /= y._value;
    return x;
}


template<class X>
UnivariateFirstDifferential<X> operator+(const UnivariateFirstDifferential<X>& x)
{
    return UnivariateFirstDifferential<X>(+x._value,+x._gradient);
}

template<class X>
UnivariateFirstDifferential<X> operator-(const UnivariateFirstDifferential<X>& x)
{
    return UnivariateFirstDifferential<X>(-x._value,-x._gradient);
}


template<class X>
UnivariateFirstDifferential<X> operator+(const UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    return UnivariateFirstDifferential<X>(x._value+y._value,x._gradient+y._gradient);
}

template<class X>
UnivariateFirstDifferential<X> operator-(const UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    return UnivariateFirstDifferential<X>(x._value-y._value,x._gradient-y._gradient);
}

template<class X>
UnivariateFirstDifferential<X> operator*(const UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    return UnivariateFirstDifferential<X>(x._value*y._value,x._value*y._gradient+y._value*x._gradient);
}

template<class X>
UnivariateFirstDifferential<X> operator/(const UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    return UnivariateFirstDifferential<X>(x._value/y._value,(x._gradient-(x._value/y._value)*y._gradient)/y._value);
}

template<class X> UnivariateFirstDifferential<X> create_zero(const UnivariateFirstDifferential<X>& c) {
    return UnivariateFirstDifferential<X>(create_zero(c.value())); }






template<class X>
UnivariateFirstDifferential<X>
min(const UnivariateFirstDifferential<X>& x1, const UnivariateFirstDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(UnivariateFirstDifferential<X> x1, UnivariateFirstDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}


template<class X>
UnivariateFirstDifferential<X>
max(const UnivariateFirstDifferential<X>& x1,const UnivariateFirstDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"max(UnivariateFirstDifferential<X> x1, UnivariateFirstDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
UnivariateFirstDifferential<X>
abs(const UnivariateFirstDifferential<X>& x)
{
    if(x.value()==0) {
        ARIADNE_THROW(std::runtime_error,"abs(UnivariateFirstDifferential<X> x)","x[0]==0");
    }
    return x.value()>0 ? pos(x) : neg(x);
}


template<class X>
UnivariateFirstDifferential<X>
pos(const UnivariateFirstDifferential<X>& x)
{
    return x;
}

template<class X>
UnivariateFirstDifferential<X>
neg(const UnivariateFirstDifferential<X>& x)
{
    return -x;
}

template<class X>
UnivariateFirstDifferential<X> rec(const UnivariateFirstDifferential<X>& x)
{
    return UnivariateFirstDifferential<X>( rec(x._value), x._gradient * (neg(sqr(rec(x._value)))) );
}

template<class X>
UnivariateFirstDifferential<X> sqr(const UnivariateFirstDifferential<X>& x)
{
    return UnivariateFirstDifferential<X>( sqr(x._value), (2*x._value)*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> pow(const UnivariateFirstDifferential<X>& x, Int n)
{
    return UnivariateFirstDifferential<X>( pow(x._value,n), (n*pow(x._value,n-1))*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> sqrt(const UnivariateFirstDifferential<X>& x)
{
    X sqrt_val = sqrt(x._value);
    return UnivariateFirstDifferential<X>( sqrt_val, rec(2*sqrt_val)*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> exp(const UnivariateFirstDifferential<X>& x)
{
    X exp_val = exp(x._value);
    return UnivariateFirstDifferential<X>( exp_val, exp_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> log(const UnivariateFirstDifferential<X>& x)
{
    X log_val = log(x._value);
    X rec_val = rec(x._value);
    return UnivariateFirstDifferential<X>( log_val, rec_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> sin(const UnivariateFirstDifferential<X>& x)
{
    X sin_val = sin(x._value);
    X cos_val = cos(x._value);
    return UnivariateFirstDifferential<X>( sin_val, cos_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> cos(const UnivariateFirstDifferential<X>& x)
{
    X cos_val = cos(x._value);
    X neg_sin_val = neg(sin(x._value));
    return UnivariateFirstDifferential<X>( cos_val, neg_sin_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> tan(const UnivariateFirstDifferential<X>& x)
{
    X tan_val = tan(x._value);
    X sqr_sec_val = sqr(rec(cos(x._value)));
    return UnivariateFirstDifferential<X>( tan_val, sqr_sec_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> asin(const UnivariateFirstDifferential<X>& x)
{
    X asin_val = asin(x._value);
    X d_asin_val = rec(sqrt(1.0-sqr(x._value)));
    return UnivariateFirstDifferential<X>( asin_val, d_asin_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> acos(const UnivariateFirstDifferential<X>& x)
{
    X acos_val = acos(x._value);
    X d_acos_val = neg(rec(sqrt(1.0-sqr(x._value))));
    return UnivariateFirstDifferential<X>( acos_val, d_acos_val*x._gradient );
}

template<class X>
UnivariateFirstDifferential<X> atan(const UnivariateFirstDifferential<X>& x)
{
    X atan_val = atan(x._value);
    X d_atan_val = rec(1+sqr(x._value));
    return UnivariateFirstDifferential<X>( atan_val, d_atan_val*x._gradient );
}


template<class X>
Bool
operator>=(const UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    return x._value>=y._value;
}


template<class X>
Bool
operator<=(const UnivariateFirstDifferential<X>& x, const UnivariateFirstDifferential<X>& y)
{
    return x._value<=y._value;
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator>=(const UnivariateFirstDifferential<X>& x, const R& c)
{
    return x._value>=static_cast<X>(c);
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator<=(const UnivariateFirstDifferential<X>& x, const R& c)
{
    return x._value<=static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator> (const UnivariateFirstDifferential<X>& x, const R& c)
{
    return x._value> static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator< (const UnivariateFirstDifferential<X>& x, const R& c)
{
    return x._value< static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator>=(const R& c, const UnivariateFirstDifferential<X>& x)
{
    return static_cast<X>(c)>=x._value;
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator<=(const R& c, const UnivariateFirstDifferential<X>& x)
{
    return static_cast<X>(c)<=x._value;
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator> (const R& c, const UnivariateFirstDifferential<X>& x)
{
    return static_cast<X>(c)> x._value;
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator< (const R& c, const UnivariateFirstDifferential<X>& x)
{
    return static_cast<X>(c)< x._value;
}




template<class X>
OutputStream& operator<<(OutputStream& os, const UnivariateFirstDifferential<X>& x)
{
    os << "D<R"<<x.argument_size()<<","<<x.degree()<<">{ ";
    os << x._value << "; " << x._gradient;
    return os << " }";
}








//! \ingroup DifferentiationModule
//! \brief A class representing the value, first and second derivatives of a scalar quantity
//! depending on a single argument.
template<class X>
class UnivariateSecondDifferential
{
  public:
    X _value;
    X _gradient;
    X _hessian;
    static const X _zero;
  public:
    //! \brief The type of used to represent numbers in the coefficient.
    typedef typename X::NumericType NumericType;
    //! \brief The type of used to represent the coefficients.
    typedef X ValueType;

    //! \brief Constructs a differential with degree zero in \a as variables. (Deprecated)
    explicit UnivariateSecondDifferential() : _value(_zero), _gradient(_zero), _hessian(_zero) { }

    //! \brief Constructs a differential with degree zero in \a as variables. (Deprecated)
    explicit UnivariateSecondDifferential(const X& v, const X& g, const X& h) : _value(v), _gradient(g), _hessian(h) { }

    //! \brief Conversion constructor from a different numerical type.
    template<class XX> UnivariateSecondDifferential(const UnivariateSecondDifferential<XX>& x)
        : _value(x._value), _gradient(x._gradient), _hessian(x._hessian) { }

    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    UnivariateSecondDifferential<X>& operator=(const X& c) {
        _value=c; _gradient=_zero; _hessian=_zero; return *this; }

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static UnivariateSecondDifferential<X> constant(const X& c) {
        UnivariateSecondDifferential<X> r; r._value=c; return r; }
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static UnivariateSecondDifferential<X> variable(const X& v) {
        UnivariateSecondDifferential<X> r; r._value=v; r._gradient=1; return r; }

    //! \brief Equality operator.
    Bool operator==(const UnivariateSecondDifferential<X>& other) const {
        return this->_value==other._value && this->_gradient==other._gradient  && this->_hessian==other._hessian; }

    //! \brief Inequality operator.
    Bool operator!=(const UnivariateSecondDifferential<X>& other) const {
        return !(*this==other); }

    //! \brief The number of independent variables.
    Nat argument_size() const { return 1u; }
    //! \brief The maximum degree of the stored terms.
    Nat degree() const { return 2u; }
    //! \brief The value of the differential i.e. the coefficient of \f$1\f$.
    const X& value() const { return this->_value; }
    //! \brief The first derivative.
    const X& gradient() const { return this->_gradient; }
    //! \brief The second derivative.
    X hessian() const { return this->_hessian; }


    //! \brief Set the coefficient of the constant term to \a c.
    Void set_value(const X& c) { this->_value=c; }

    //! \brief Set all coefficients to zero.
    Void clear() { *this=_zero; }

};

template<class X>
const X UnivariateSecondDifferential<X>::_zero=X(0);




template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X>&>
operator+=(UnivariateSecondDifferential<X>& x, const R& c)
{
    x._value+=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X>&>
operator-=(UnivariateSecondDifferential<X>& x, const R& c)
{
    x._value-=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X>&>
operator*=(UnivariateSecondDifferential<X>& x, const R& c)
{
    x._value*=static_cast<X>(c);
    x._gradient*=static_cast<X>(c);
    x._hessian*=static_cast<X>(c);
    return x;
}


template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X>&>
operator/=(UnivariateSecondDifferential<X>& x, const R& c)
{
    x._value/=static_cast<X>(c);
    x._gradient/=static_cast<X>(c);
    x._hessian/=static_cast<X>(c);
    return x;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator+(const UnivariateSecondDifferential<X>& x, const R& c)
{
    UnivariateSecondDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator+(const R& c, const UnivariateSecondDifferential<X>& x)
{
    UnivariateSecondDifferential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator-(const UnivariateSecondDifferential<X>& x, const R& c)
{
    UnivariateSecondDifferential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator-(const R& c, const UnivariateSecondDifferential<X>& x)
{
    UnivariateSecondDifferential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>

EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator*(const UnivariateSecondDifferential<X>& x, const R& c)
{
    UnivariateSecondDifferential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator*(const R& c, const UnivariateSecondDifferential<X>& x)
{
    UnivariateSecondDifferential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator/(const UnivariateSecondDifferential<X>& x, const R& c)
{
    UnivariateSecondDifferential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
EnableIfNumericType<R,UnivariateSecondDifferential<X> >
operator/(const R& c, const UnivariateSecondDifferential<X>& x)
{
    UnivariateSecondDifferential<X> r(rec(x)); r*=X(c); return r;
}



template<class X>
UnivariateSecondDifferential<X>& operator+=(UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    x._value += y._value;
    x._gradient += y._gradient;
    x._hessian += y._hessian;
    return x;
}

template<class X>
UnivariateSecondDifferential<X>& operator-=(UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    x._value -= y._value;
    x._gradient -= y._gradient;
    x._hessian -= y._hessian;
    return x;
}

template<class X>
UnivariateSecondDifferential<X>& operator*=(UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    x._hessian *= y._value;
    x._hessian += 2*x._gradient*y._gradient;
    x._hessian += x._value * y._hessian;
    x._gradient *= y._value;
    x._gradient += x._value * y._gradient;
    x._value *= y._value;
    return x;
}


template<class X>
UnivariateSecondDifferential<X> operator+(const UnivariateSecondDifferential<X>& x)
{
    return UnivariateSecondDifferential<X>(+x._value,+x._gradient,+x._hessian);
}

template<class X>
UnivariateSecondDifferential<X> operator-(const UnivariateSecondDifferential<X>& x)
{
    return UnivariateSecondDifferential<X>(-x._value,-x._gradient,-x._hessian);
}


template<class X>
UnivariateSecondDifferential<X> operator+(const UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    return UnivariateSecondDifferential<X>(x._value+y._value,x._gradient+y._gradient,x._hessian+y._hessian);
}

template<class X>
UnivariateSecondDifferential<X> operator-(const UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    return UnivariateSecondDifferential<X>(x._value-y._value,x._gradient-y._gradient,x._hessian-y._hessian);
}

template<class X>
UnivariateSecondDifferential<X> operator*(const UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    return UnivariateSecondDifferential<X>(x._value*y._value,x._value*y._gradient+y._value*x._gradient,
                                           x._value*y._hessian+2*x._gradient*y._gradient+y._value*x._hessian);
}

template<class X>
UnivariateSecondDifferential<X> operator/(const UnivariateSecondDifferential<X>& x, const UnivariateSecondDifferential<X>& y)
{
    return UnivariateSecondDifferential<X>(x._value/y._value,(x._gradient-(x._value/y._value)*y._gradient)/y._value,
                                           (x._hessian/y._value)-(x._value/y._value)*(y._hessian/y._value)
                                               -2*(x._gradient/y._value)*(y._gradient/y._value)+2*(x._value/y._value)*sqr(y._gradient/y._value));
    //return x*rec(y);
}







template<class X>
UnivariateSecondDifferential<X>
min(const UnivariateSecondDifferential<X>& x1, const UnivariateSecondDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(UnivariateSecondDifferential<X> x1, UnivariateSecondDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}


template<class X>
UnivariateSecondDifferential<X>
max(const UnivariateSecondDifferential<X>& x1,const UnivariateSecondDifferential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"max(UnivariateSecondDifferential<X> x1, UnivariateSecondDifferential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
UnivariateSecondDifferential<X>
abs(const UnivariateSecondDifferential<X>& x)
{
    if(x.value()==0) {
        ARIADNE_THROW(std::runtime_error,"abs(UnivariateSecondDifferential<X> x)","x[0]==0");
    }
    return x.value()>0 ? pos(x) : neg(x);
}

template<class X>
UnivariateSecondDifferential<X>
pos(const UnivariateSecondDifferential<X>& x)
{
    return x;
}

template<class X>
UnivariateSecondDifferential<X>
neg(const UnivariateSecondDifferential<X>& x)
{
    return -x;
}

template<class X>
UnivariateSecondDifferential<X> rec(const UnivariateSecondDifferential<X>& x)
{
    X rec_val = rec(x._value);
    X neg_sqr_rec_val = neg(sqr(rec_val));
    X dbl_cub_rec_val = -2*rec_val*neg_sqr_rec_val;
    return UnivariateSecondDifferential<X>( rec_val, neg_sqr_rec_val*x._gradient, neg_sqr_rec_val*x._hessian + dbl_cub_rec_val * sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> sqr(const UnivariateSecondDifferential<X>& x)
{
    return UnivariateSecondDifferential<X>( sqr(x._value), (2*x._value)*x._gradient, (2*x._value)*x._hessian + sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> pow(const UnivariateSecondDifferential<X>& x, Int n)
{
    return UnivariateSecondDifferential<X>( pow(x._value,n), (n*pow(x._value,n-1))*x._gradient, n*(n-1)*pow(x._value,n-2)*x._hessian );
}

//ddf(y)/dxdx = d/dx ( f'(y) dy/dx) = f''(y) dy/dx dy/dx + f'(y) ddy/dxdx

template<class X>
UnivariateSecondDifferential<X> sqrt(const UnivariateSecondDifferential<X>& x)
{
    X sqrt_val = sqrt(x._value);
    X rec_dbl_sqrt_val = rec(2*sqrt_val);
    X neg_rec_quad_pow_val = neg(rec(4*sqrt_val*x._value));
    return UnivariateSecondDifferential<X>( sqrt_val, rec_dbl_sqrt_val*x._gradient, rec_dbl_sqrt_val*x._hessian + neg_rec_quad_pow_val * sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> exp(const UnivariateSecondDifferential<X>& x)
{
    X exp_val = exp(x._value);
    return UnivariateSecondDifferential<X>( exp_val, exp_val*x._gradient, exp_val*x._hessian+exp_val*sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> log(const UnivariateSecondDifferential<X>& x)
{
    X log_val = log(x._value);
    X rec_val = rec(x._value);
    X neg_sqr_rec_val = neg(sqr(rec_val));
    return UnivariateSecondDifferential<X>( log_val, rec_val*x._gradient, rec_val*x._hessian+neg_sqr_rec_val*sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> sin(const UnivariateSecondDifferential<X>& x)
{
    X sin_val = sin(x._value);
    X cos_val = cos(x._value);
    X neg_sin_val = neg(sin_val);
    return UnivariateSecondDifferential<X>( sin_val, cos_val*x._gradient, cos_val*x._hessian+neg_sin_val*sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> cos(const UnivariateSecondDifferential<X>& x)
{
    X cos_val = cos(x._value);
    X neg_sin_val = neg(sin(x._value));
    X neg_cos_val = neg(cos_val);
    return UnivariateSecondDifferential<X>( cos_val, neg_sin_val*x._gradient, neg_sin_val*x._hessian+neg_cos_val*sqr(x._gradient) );
}

template<class X>
UnivariateSecondDifferential<X> tan(const UnivariateSecondDifferential<X>& x)
{
    X tan_val = tan(x._value);
    X sqr_sec_val = sqr(rec(cos(x._value)));
    X dbl_tan_sqr_sec_val = 2*tan_val*sqr_sec_val;
    return UnivariateSecondDifferential<X>( tan_val, sqr_sec_val*x._gradient, sqr_sec_val*x._hessian+dbl_tan_sqr_sec_val*sqr(x._gradient) );
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator>=(const UnivariateSecondDifferential<X>& x, const R& c)
{
    return x._value>=static_cast<X>(c);
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator<=(const UnivariateSecondDifferential<X>& x, const R& c)
{
    return x._value<=static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator> (const UnivariateSecondDifferential<X>& x, const R& c)
{
    return x._value> static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator< (const UnivariateSecondDifferential<X>& x, const R& c)
{
    return x._value< static_cast<X>(c);
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator>=(const R& c, const UnivariateSecondDifferential<X>& x)
{
    return static_cast<X>(c)>=x._value;
}


template<class X, class R>
EnableIfNumericType<R,Bool>
operator<=(const R& c, const UnivariateSecondDifferential<X>& x)
{
    return static_cast<X>(c)<=x._value;
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator> (const R& c, const UnivariateSecondDifferential<X>& x)
{
    return static_cast<X>(c)> x._value;
}

template<class X, class R>
EnableIfNumericType<R,Bool>
operator< (const R& c, const UnivariateSecondDifferential<X>& x)
{
    return static_cast<X>(c)< x._value;
}




template<class X>
OutputStream& operator<<(OutputStream& os, const UnivariateSecondDifferential<X>& x)
{
     //e.graded_sort();
    os << "D<R"<<x.argument_size()<<","<<x.degree()<<">{";
    os << x._value << ";" << x._gradient << ";" << x._hessian;
    return os << " }";
}





} //namespace Ariadne

#endif // ARIADNE_FIXED_UNIVARIATE_DIFFERENTIAL_H
