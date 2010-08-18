
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
    explicit Real() : _ivl() { }
    explicit Real(const std::string& s);
    Real(unsigned int m) : _ivl(m) { }
    Real(int n) : _ivl(n) { }
    Real(double x) : _ivl(x) { }
#ifdef HAVE_GMPXX_H
    Real(const Rational& q);
#endif
    explicit Real(const Float& x) : _ivl(x) { }
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
    double get_d() const { return this->_ivl.get_d(); }
  private:
    friend Float::Float(const Real&);
    friend Interval::Interval(const Real&);
};

inline Float::Float(const Real& x) : v(x._ivl.midpoint().v) { }
inline Interval::Interval(const Real& x) : l(x._ivl.l), u(x._ivl.u) { }
inline Float& Float::operator=(const Real& x) { *this=Float(x); return *this; }
inline Interval& Interval::operator=(const Real& x) { *this=Interval(x); return *this; }

inline Real operator+(const Real& x) { return Real(+static_cast<Interval>(x)); }
inline Real operator-(const Real& x) { return Real(-static_cast<Interval>(x)); }
inline Real operator+(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)+static_cast<Interval>(y)); }
inline Real operator-(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)-static_cast<Interval>(y)); }
inline Real operator*(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)*static_cast<Interval>(y)); }
inline Real operator/(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)/static_cast<Interval>(y)); }
inline Real& operator+=(Real& x, const Real& y) { return x=x+y; }
inline Real& operator-=(Real& x, const Real& y) { return x=x-y; }
inline Real& operator*=(Real& x, const Real& y) { return x=x*y; }
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

inline tribool operator==(const Real& x, const Real& y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
inline tribool operator!=(const Real& x, const Real& y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
inline tribool operator>=(const Real& x, const Real& y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
inline tribool operator<=(const Real& x, const Real& y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
inline tribool operator> (const Real& x, const Real& y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
inline tribool operator< (const Real& x, const Real& y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

inline tribool operator==(const Real& x, double y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
inline tribool operator!=(const Real& x, double y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
inline tribool operator>=(const Real& x, double y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
inline tribool operator<=(const Real& x, double y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
inline tribool operator> (const Real& x, double y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
inline tribool operator< (const Real& x, double y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

template<> inline Real pi<>() { return Real(pi<Interval>()); }
inline Real abs(const Real& x) { return Real(abs(static_cast<Interval>(x))); }
inline Real rec(const Real& x) { return Real(rec(static_cast<Interval>(x))); }
inline Real sqr(const Real& x) { return Real(sqr(static_cast<Interval>(x))); }
inline Real sqrt(const Real& x) { return Real(sqrt(static_cast<Interval>(x))); }
inline Real exp(const Real& x) { return Real(exp(static_cast<Interval>(x))); }
inline Real log(const Real& x) { return Real(log(static_cast<Interval>(x))); }
inline Real sin(const Real& x) { return Real(sin(static_cast<Interval>(x))); }
inline Real cos(const Real& x) { return Real(cos(static_cast<Interval>(x))); }
inline Real tan(const Real& x) { return Real(tan(static_cast<Interval>(x))); }

inline std::ostream& operator<<(std::ostream& os, const Real& x) {
    Interval ivl=static_cast<Interval>(x); return os << "Real(" << ivl.lower() <<',' << ivl.upper() << ")"; }

template<class R, class A> inline R numeric_cast(const A& a);
template<> inline double numeric_cast(const Real& a) { return a.get_d(); }
template<> inline Real numeric_cast(const Float& a) { return Real(a); }
template<> inline Real numeric_cast(const Interval& a) { return Real(a); }

} // namespace Ariadne

#endif
