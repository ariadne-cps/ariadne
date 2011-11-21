/***************************************************************************
 *            real.cc
 *
 *  Copyright 2008-10  Pieter Collins
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

#include <iostream>
#include <iomanip>
#include <cassert>
#include <limits>

#include "config.h"
#include "real.h"


namespace Ariadne {

namespace{
static const double _pi_up=3.1415926535897936;
static const double _pi_approx=3.1415926535897931;
static const double _pi_down=3.1415926535897931;
static const double _infinity=std::numeric_limits<double>::infinity();
}

Real::~Real() { }
Real::Real() : _ivl() { }
Real::Real(double l, double u) : _ivl(l,u) { }
Real::Real(double l, double x, double u) : _ivl(l,u) { }

Real::Real(unsigned int m) : _ivl(m) { }
Real::Real(int n) : _ivl(n) { }
Real::Real(double x) : _ivl(x) { }
Real::Real(const ExactFloat& x) : _ivl(x.value()) { }
Real::Real(const Real& x) : _ivl(x._ivl) { }
Real& Real::operator=(const double& x) { this->_ivl=x; return *this; }
Real& Real::operator=(const ExactFloat& x) { this->_ivl=x.value(); return *this; }
Real& Real::operator=(const Real& x) { this->_ivl=x._ivl; return *this; }
double Real::get_d() const { return this->_ivl.get_d(); }

Float::Float(const Real& x) : v(x._ivl.midpoint().v) { }
Interval::Interval(const Real& x) : l(x._ivl.l), u(x._ivl.u) { }
Float& Float::operator=(const Real& x) { *this=Float(x); return *this; }
Interval& Interval::operator=(const Real& x) { *this=Interval(x); return *this; }

Real _make_real(const Interval& ivl) { return Real(ivl.lower().get_d(),ivl.upper().get_d()); }

Real operator+(const Real& x) { return _make_real(+static_cast<Interval>(x)); }
Real operator-(const Real& x) { return _make_real(-static_cast<Interval>(x)); }
Real operator+(const Real& x, const Real& y) { return _make_real(static_cast<Interval>(x)+static_cast<Interval>(y)); }
Real operator-(const Real& x, const Real& y) { return _make_real(static_cast<Interval>(x)-static_cast<Interval>(y)); }
Real operator*(const Real& x, const Real& y) { return _make_real(static_cast<Interval>(x)*static_cast<Interval>(y)); }
Real operator/(const Real& x, const Real& y) { return _make_real(static_cast<Interval>(x)/static_cast<Interval>(y)); }

Float mag(const Real& x) { return mag(static_cast<Interval>(x)); }

const Real pi=Real(_pi_down,_pi_approx,_pi_up);
const Real infinity=Real(_infinity,_infinity,_infinity);

Real abs(const Real& x) { return _make_real(abs(static_cast<Interval>(x))); }
Real pos(const Real& x) { return _make_real(pos(static_cast<Interval>(x))); }
Real neg(const Real& x) { return _make_real(neg(static_cast<Interval>(x))); }
Real sqr(const Real& x) { return _make_real(sqr(static_cast<Interval>(x))); }
Real rec(const Real& x) { return _make_real(rec(static_cast<Interval>(x))); }
Real add(const Real& x, const Real& y) { return _make_real(add(static_cast<Interval>(x),static_cast<Interval>(y))); }
Real sub(const Real& x, const Real& y) { return _make_real(sub(static_cast<Interval>(x),static_cast<Interval>(y))); }
Real mul(const Real& x, const Real& y) { return _make_real(mul(static_cast<Interval>(x),static_cast<Interval>(y))); }
Real div(const Real& x, const Real& y) { return _make_real(div(static_cast<Interval>(x),static_cast<Interval>(y))); }
Real pow(const Real& x, uint m) { return _make_real(pow(static_cast<Interval>(x),m)); }
Real pow(const Real& x, int n) { return _make_real(pow(static_cast<Interval>(x),n)); }
Real sqrt(const Real& x) { return _make_real(sqrt(static_cast<Interval>(x))); }
Real exp(const Real& x) { return _make_real(exp(static_cast<Interval>(x))); }
Real log(const Real& x) { return _make_real(log(static_cast<Interval>(x))); }
Real sin(const Real& x) { return _make_real(sin(static_cast<Interval>(x))); }
Real cos(const Real& x) { return _make_real(cos(static_cast<Interval>(x))); }
Real tan(const Real& x) { return _make_real(tan(static_cast<Interval>(x))); }
Real asin(const Real& x) { return _make_real(asin(static_cast<Interval>(x))); }
Real acos(const Real& x) { return _make_real(acos(static_cast<Interval>(x))); }
Real atan(const Real& x) { return _make_real(atan(static_cast<Interval>(x))); }



std::ostream& operator<<(std::ostream& os, const Real& x)
{
    Interval ivl=static_cast<Interval>(x);
    return os << ivl.midpoint();
    if(ivl.lower()==ivl.upper()) { return os << ivl.lower(); }
    else { return os << ivl; }
    return os << "Real(" << ivl.lower() <<',' << ivl.upper() << ")";

}



#ifdef HAVE_GMPXX_H

Real::Real(const std::string& str)
{
    Rational q;
    bool decimal_point=false;
    uint decimal_places=0;
    int sign=1;
    const char* c_ptr=str.c_str();
    while(*c_ptr != 0) {
        const char& c=*c_ptr;
        if(c=='.') {
            if(decimal_point) {
                ARIADNE_THROW(std::runtime_error,"Real(String)","real literal \""<<str<<"\" has more than one decimal point.");
            }
            else {
                decimal_point=true;
            }
        } else if(c>='0' && c<='9') {
            q=q*10+(c-'0');
            if(decimal_point) {
                ++decimal_places;
            }
        } else if(c=='-' && c_ptr==str.c_str()) {
            sign=-1;
        } else {
            ARIADNE_THROW(std::runtime_error,"Real(String)","invalid symbol '"<<c<<"' in string literal \""<<str<<"\"");
        }
        ++c_ptr;
    }
    for(uint i=0; i!=decimal_places; ++i) {
        q=q/10;
    }
    *this=Real(sign*q);
}

Real::Real(const Rational& q)
{
    rounding_mode_t rnd=get_rounding_mode();
    double x=q.get_d();
    volatile double ml=-x;
    volatile double u=x;
    set_rounding_upward();
    while(-ml>static_cast<const mpq_class&>(q)) {
        ml+=std::numeric_limits<double>::min();
    }
    while(u<static_cast<const mpq_class&>(q)) {
        u+=std::numeric_limits<double>::min();
    }
    *this=Real(-ml,x,u);
    set_rounding_mode(rnd);
}

#else

Real::Real(const std::string& str)
{
    ARIADNE_THROW(std::runtime_error,"Need GMP library to convert string literal to Real.");
}

#endif



} // namespace Ariadne

