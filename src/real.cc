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


std::ostream& operator<<(std::ostream& os, const Real& x)
{
    Interval ivl=static_cast<Interval>(x);
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
        } else {
            ARIADNE_THROW(std::runtime_error,"Real(String)","invalid symbol '"<<c<<"' in string literal \""<<str<<"\"");
        }
        ++c_ptr;
    }
    for(uint i=0; i!=decimal_places; ++i) {
        q=q/10;
    }
    *this=Real(q);
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

