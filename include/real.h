
/***************************************************************************
 *            real.h
 *
 *  Copyright 2009  Pieter Collins
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

/*! \file real.h
 *  \brief Real numbers which can be converted to Float or Interval.
 */

#ifndef ARIADNE_REAL_H
#define ARIADNE_REAL_H

#include <iosfwd>
#include <iostream>
#include "numeric.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \brief Computable real numbers.
//!
//! Support over-approximation by an Interval and approximation by a Float.
class Real {
    Interval _ivl;
  public:
    typedef Real ScalarType;
  public:
    //! Default constructor yields the exact value \a 0.
    Real() : _ivl() { }
    //! \brief Construct from a string literal.
    //! This can be used to create real numbers representing decimal values e.g. \c %Real("4.2") .
    explicit Real(const std::string& s);
    Real(unsigned int m) : _ivl(m) { }
    Real(int n) : _ivl(n) { }
    //! \brief Convert from a builtin double-precision floating-point value.
    //! A numeric literal is first processed by the language support, and the resulting %Real may not have
    //! the same value as the mathematical literal. e.g. \c %Real(4.2) has the value 4.2000000000000002 to 16 decimal places.
    Real(double x) : _ivl(x) { }
#ifdef HAVE_GMPXX_H
    //! \brief Construct from a rational number.
    Real(const Rational& q);
#endif
    //! \brief Construct from a floating-point value.
    explicit Real(const Float& x) : _ivl(x) { }
    //! \brief Construct from a interval. The resulting %Real object does not describe a number arbitrarily accurately.
    //! \deprecated This constructor should be avoided in user code.
    explicit Real(const Interval& ivl) : _ivl(ivl) { }
    explicit Real(double l, double u) : _ivl(l,u) { }
    explicit Real(double l, double x, double u) : _ivl(l,u) { }
    Real& operator=(const double& x) { this->_ivl=x; return *this; }
    Real& operator=(const Float& x) { this->_ivl=x; return *this; }
    Real& operator=(const Interval& x) { this->_ivl=x; return *this; }
    // Can't use conversion operators below in g++ since compiler complains
    // about ambiguous conversion to Interval through Interval(Real::operator Float())
    //operator Float() const { return this->_ivl.midpoint(); }
    //operator Interval() const { return this->_ivl; }
  public:
    //! \brief Get an approximation as a builtin double-precision floating-point number.
    double get_d() const { return this->_ivl.get_d(); }
  private:
    friend Float::Float(const Real&);
    friend Interval::Interval(const Real&);
};

inline Float::Float(const Real& x) : v(x._ivl.midpoint().v) { }
inline Interval::Interval(const Real& x) : l(x._ivl.l), u(x._ivl.u) { }
inline Float& Float::operator=(const Real& x) { *this=Float(x); return *this; }
inline Interval& Interval::operator=(const Real& x) { *this=Interval(x); return *this; }

//! \related Real \brief Unary plus operator.
inline Real operator+(const Real& x) { return Real(+static_cast<Interval>(x)); }
//! \related Real \brief Unary negation operator.
inline Real operator-(const Real& x) { return Real(-static_cast<Interval>(x)); }
//! \related Real \brief Binary addition operator.
inline Real operator+(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)+static_cast<Interval>(y)); }
//! \related Real \brief Subtraction operator.
inline Real operator-(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)-static_cast<Interval>(y)); }
//! \related Real \brief Binary multiplication operator.
inline Real operator*(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)*static_cast<Interval>(y)); }
//! \related Real \brief Division operator.
inline Real operator/(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)/static_cast<Interval>(y)); }
//! \related Real \brief Inplace addition operator.
inline Real& operator+=(Real& x, const Real& y) { return x=x+y; }
//! \related Real \brief Inplace subtraction operator.
inline Real& operator-=(Real& x, const Real& y) { return x=x-y; }
//! \related Real \brief Inplace multiplication operator.
inline Real& operator*=(Real& x, const Real& y) { return x=x*y; }
//! \related Real \brief Inplace division operator.
inline Real& operator/=(Real& x, const Real& y) { return x=x/y; }

inline Real operator+(const Real& x, double y) { return x+static_cast<Real>(y); }
inline Real operator-(const Real& x, double y) { return x-static_cast<Real>(y); }
inline Real operator*(const Real& x, double y) { return x*static_cast<Real>(y); }
inline Real operator/(const Real& x, double y) { return x/static_cast<Real>(y); }
inline Real operator+(double x, const Real& y) { return static_cast<Real>(x)+y; }
inline Real operator-(double x, const Real& y) { return static_cast<Real>(x)-y; }
inline Real operator*(double x, const Real& y) { return static_cast<Real>(x)*y; }
inline Real operator/(double x, const Real& y) { return static_cast<Real>(x)/y; }

inline Float operator+(const Real& x, const Float& y) { return static_cast<Float>(x)+y; }
inline Float operator-(const Real& x, const Float& y) { return static_cast<Float>(x)-y; }
inline Float operator*(const Real& x, const Float& y) { return static_cast<Float>(x)*y; }
inline Float operator/(const Real& x, const Float& y) { return static_cast<Float>(x)/y; }
inline Float operator+(const Float& x, const Real& y) { return x+static_cast<Float>(y); }
inline Float operator-(const Float& x, const Real& y) { return x-static_cast<Float>(y); }
inline Float operator*(const Float& x, const Real& y) { return x*static_cast<Float>(y); }
inline Float operator/(const Float& x, const Real& y) { return x/static_cast<Float>(y); }
inline Float& operator+=(Float& x, const Real& y) { return x+=static_cast<Float>(y); }
inline Float& operator-=(Float& x, const Real& y) { return x-=static_cast<Float>(y); }
inline Float& operator*=(Float& x, const Real& y) { return x*=static_cast<Float>(y); }
inline Float& operator/=(Float& x, const Real& y) { return x/=static_cast<Float>(y); }

inline Interval operator+(const Real& x, const Interval& y) { return static_cast<Interval>(x)+y; }
inline Interval operator-(const Real& x, const Interval& y) { return static_cast<Interval>(x)-y; }
inline Interval operator*(const Real& x, const Interval& y) { return static_cast<Interval>(x)*y; }
inline Interval operator/(const Real& x, const Interval& y) { return static_cast<Interval>(x)/y; }
inline Interval operator+(const Interval& x, const Real& y) { return x+static_cast<Interval>(y); }
inline Interval operator-(const Interval& x, const Real& y) { return x-static_cast<Interval>(y); }
inline Interval operator*(const Interval& x, const Real& y) { return x*static_cast<Interval>(y); }
inline Interval operator/(const Interval& x, const Real& y) { return x/static_cast<Interval>(y); }
inline Interval& operator+=(Interval& x, const Real& y) { return x+=static_cast<Interval>(y); }
inline Interval& operator-=(Interval& x, const Real& y) { return x-=static_cast<Interval>(y); }
inline Interval& operator*=(Interval& x, const Real& y) { return x*=static_cast<Interval>(y); }
inline Interval& operator/=(Interval& x, const Real& y) { return x/=static_cast<Interval>(y); }

//! \related Real \brief Equality operator.
//! Returns \c true if the numbers have the same representation or can be evaluated exactly to the same value.
//! Returns \c false if the numbers can be evalated to different values.
//! Returns \c indeterminate if the numbers cannot be shown to be the same or different. Implementation dependent.
inline tribool operator==(const Real& x, const Real& y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
//! \related Real \brief Inequality operator.
//! Returns \c true if the numbers have been proved to be different, such as by evaluation to different values.
//! Returns \c false if the numbers have the same representation or can be evaluated exactly to the same value.
//! Returns \c indeterminate if the numbers cannot be shown to be the same or different. Implementation dependent.
inline tribool operator!=(const Real& x, const Real& y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
//! \related Real \brief Greater-than-or-equal-to comparison operator.
inline tribool operator>=(const Real& x, const Real& y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
//! \related Real \brief Less-than-or-equal-to comparison operator.
inline tribool operator<=(const Real& x, const Real& y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
//! \related Real \brief Strictly-greater-than comparison operator.
inline tribool operator> (const Real& x, const Real& y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
//! \related Real \brief Strictly-less-than comparison operator.
inline tribool operator< (const Real& x, const Real& y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

inline tribool operator==(const Real& x, double y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
inline tribool operator!=(const Real& x, double y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
inline tribool operator>=(const Real& x, double y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
inline tribool operator<=(const Real& x, double y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
inline tribool operator> (const Real& x, double y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
inline tribool operator< (const Real& x, double y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

//! \related Real \brief The absolute value function \c |x|.
inline Real abs(const Real& x) { return Real(abs(static_cast<Interval>(x))); }
//! \related Real \brief The unary plus function \c +x.
inline Real pos(const Real& x) { return Real(pos(static_cast<Interval>(x))); }
//! \related Real \brief The unary negation function \c -x.
inline Real neg(const Real& x) { return Real(neg(static_cast<Interval>(x))); }
//! \related Real \brief The reciprocal function \c 1/x.
inline Real sqr(const Real& x) { return Real(sqr(static_cast<Interval>(x))); }
//! \related Real \brief The reciprocal function \c 1/x.
inline Real rec(const Real& x) { return Real(rec(static_cast<Interval>(x))); }
//! \related Real \brief The binary addition function \c x+y.
inline Real add(const Real& x, const Real& y) { return Real(add(static_cast<Interval>(x),static_cast<Interval>(y))); }
//! \related Real \brief The subtraction function \c x-y.
inline Real sub(const Real& x, const Real& y) { return Real(sub(static_cast<Interval>(x),static_cast<Interval>(y))); }
//! \related Real \brief The binary multiplication function \c x*y.
inline Real mul(const Real& x, const Real& y) { return Real(mul(static_cast<Interval>(x),static_cast<Interval>(y))); }
//! \related Real \brief The division function \c x/y.
inline Real div(const Real& x, const Real& y) { return Real(div(static_cast<Interval>(x),static_cast<Interval>(y))); }
//! \related Real \brief The integer power function \c x^n.
inline Real pow(const Real& x, int n) { return Real(pow(static_cast<Interval>(x),n)); }
//! \related Real \brief The square-root function.
inline Real sqrt(const Real& x) { return Real(sqrt(static_cast<Interval>(x))); }
//! \related Real \brief The exponential function.
inline Real exp(const Real& x) { return Real(exp(static_cast<Interval>(x))); }
//! \related Real \brief The natural logarithm function.
inline Real log(const Real& x) { return Real(log(static_cast<Interval>(x))); }
//! \related Real \brief The constant \a pi.
template<> inline Real pi<>() { return Real(pi<Interval>()); }
//! \related Real \brief The sine function.
inline Real sin(const Real& x) { return Real(sin(static_cast<Interval>(x))); }
//! \related Real \brief The cosine function.
inline Real cos(const Real& x) { return Real(cos(static_cast<Interval>(x))); }
//! \related Real \brief The tangent function.
inline Real tan(const Real& x) { return Real(tan(static_cast<Interval>(x))); }

//! \related Real \brief Write to an output stream.
inline std::ostream& operator<<(std::ostream& os, const Real& x) {
    Interval ivl=static_cast<Interval>(x); return os << "Real(" << ivl.lower() <<',' << ivl.upper() << ")"; }

template<class R, class A> inline R numeric_cast(const A& a);
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Real numeric_cast(const Float& a) { return Real(a); }
template<> inline Real numeric_cast(const Interval& a) { return Real(a); }

} // namespace Ariadne

#endif
