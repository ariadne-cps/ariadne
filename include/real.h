
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

class Real : public Interval {
  public:
    explicit Real() : Interval() { }
    Real(int n) : Interval(n) { }
    Real(double x) : Interval(x) { }
    //explicit Real(Float flt) : Interval(flt) { }
    explicit Real(Interval ivl) : Interval(ivl) { }
    explicit Real(double l, double u) : Interval(l,u) { }
    Real& operator=(const Float& x) { static_cast<Interval&>(*this)=x; return *this; }
    Real& operator=(const Interval& x) { static_cast<Interval&>(*this)=x; return *this; }
    operator Float() const { return this->midpoint(); }
};

inline Real operator+(const Real& x) { return Real(+static_cast<Interval>(x)); }
inline Real operator-(const Real& x) { return Real(-static_cast<Interval>(x)); }
inline Real operator+(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)+static_cast<Interval>(y)); }
inline Real operator-(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)-static_cast<Interval>(y)); }
inline Real operator*(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)*static_cast<Interval>(y)); }
inline Real operator/(const Real& x, const Real& y) { return Real(static_cast<Interval>(x)/static_cast<Interval>(y)); }

//inline Real operator+(const Real& x, const int& y) { return static_cast<Real>(static_cast<const Interval&>(x)+y); }
//inline Real operator-(const Real& x, const int& y) { return static_cast<Real>(static_cast<const Interval&>(x)-y); }
//inline Real operator*(const Real& x, const int& y) { return static_cast<Real>(static_cast<const Interval&>(x)*y); }
//inline Real operator/(const Real& x, const int& y) { return static_cast<Real>(static_cast<const Interval&>(x)/y); }
//inline Real operator/(const int& x, const Real& y) { return static_cast<Real>(x/static_cast<const Interval&>(y)); }
inline Float operator+(const Real& x, const Float& y) { return static_cast<Float>(x)+y; }
inline Float operator-(const Real& x, const Float& y) { return static_cast<Float>(x)-y; }
inline Float operator*(const Real& x, const Float& y) { return static_cast<Float>(x)*y; }
inline Float operator/(const Real& x, const Float& y) { return static_cast<Float>(x)/y; }
inline Float operator+(const Float& x, const Real& y) { return x+static_cast<Float>(y); }
inline Float operator-(const Float& x, const Real& y) { return x-static_cast<Float>(y); }
inline Float operator*(const Float& x, const Real& y) { return x*static_cast<Float>(y); }
inline Float operator/(const Float& x, const Real& y) { return x/static_cast<Float>(y); }

inline tribool operator==(const Real& x, const Real& y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
inline tribool operator!=(const Real& x, const Real& y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
inline tribool operator>=(const Real& x, const Real& y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
inline tribool operator<=(const Real& x, const Real& y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
inline tribool operator> (const Real& x, const Real& y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
inline tribool operator< (const Real& x, const Real& y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

inline tribool operator==(const Real& x, const int& y) { return static_cast<Interval>(x)==static_cast<Interval>(y); }
inline tribool operator!=(const Real& x, const int& y) { return static_cast<Interval>(x)!=static_cast<Interval>(y); }
inline tribool operator>=(const Real& x, const int& y) { return static_cast<Interval>(x)>=static_cast<Interval>(y); }
inline tribool operator<=(const Real& x, const int& y) { return static_cast<Interval>(x)<=static_cast<Interval>(y); }
inline tribool operator> (const Real& x, const int& y) { return static_cast<Interval>(x)> static_cast<Interval>(y); }
inline tribool operator< (const Real& x, const int& y) { return static_cast<Interval>(x)< static_cast<Interval>(y); }

template<> inline Real pi<>() { return Real(pi<Interval>()); }
inline Real rec(const Real& x) { return Real(rec(static_cast<Interval>(x))); }
inline Real sqr(const Real& x) { return Real(sqr(static_cast<Interval>(x))); }
inline Real sqrt(const Real& x) { return Real(sqrt(static_cast<Interval>(x))); }
inline Real exp(const Real& x) { return Real(exp(static_cast<Interval>(x))); }
inline Real log(const Real& x) { return Real(log(static_cast<Interval>(x))); }
inline Real sin(const Real& x) { return Real(sin(static_cast<Interval>(x))); }
inline Real cos(const Real& x) { return Real(cos(static_cast<Interval>(x))); }
inline Real tan(const Real& x) { return Real(tan(static_cast<Interval>(x))); }

inline std::ostream& operator<<(std::ostream& os, const Real& x) {
    return os << "Real(" << x.lower() <<',' << x.upper() << ")"; }


} // namespace Ariadne

#endif
