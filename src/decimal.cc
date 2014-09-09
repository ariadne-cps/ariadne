/***************************************************************************
 *            decimal.cc
 *
 *  Copyright 2014  Pieter Collins
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

#include "config.h"

#include "macros.h"
#include "decimal.h"

#include "rational.h"
#include "interval.h"
#include "float.h"

namespace Ariadne {

Decimal operator"" _dec(long double x)
{
    return Decimal(static_cast<double>(x));
}

Decimal::Decimal(double x)
{
    std::stringstream ss;
    int s=+1;
    if(x<0) { s=-1; ss<<'-'; }
    double y=s*x;
    static const uint sf=9;
    static const double acc=1e-9;
    static const double tol=1e-15;
    int exp=0; double pow=1;
    while(y<1.0) { y*=10; exp-=1; pow/=10; }
    while(y>=1.0) { y/=10; exp+=1; pow*=10; }
    long int n=std::floor(y/acc+0.5);
    double re=y-n*acc;
    if(std::fabs(re)>=tol) {
        ARIADNE_THROW(std::runtime_error,"Decimal(double)","double-precision floating-point number must have a relative error of "<<tol<<" with respect to its approximation to "<<sf<<" significant figures; number "<<std::setprecision(17)<<x<<" has a relative error of "<<re<<"");
    }
    long int m=1e9; int c=0;
    if(exp<=0) { ss << '0'; }
    if(exp<0) { ss << '.'; for(int i=0; i!=-exp; ++i) { ss << "0"; } }
    while(n!=0) {
        if(exp==0) { ss << "."; }
        m/=10;
        c=n/m;
        n=n%m;
        ss << c;
        --exp;
    }
    while(exp>0) { --exp; ss << '0'; }
    if(exp==0) { ss << '.'; }
    _str=ss.str();
}

Decimal::Decimal(std::string str)
    : _str(str)
{
    // Parse string to ensure correctness
    bool found_decimal_point=false;
    const char* c_ptr=str.c_str();
    const char& c=*c_ptr;
    if(c=='-' || c=='+') {
        ++c_ptr;
    }
    while(*c_ptr != 0) {
        const char& c=*c_ptr;
        if(c=='.') {
            if(found_decimal_point) {
                ARIADNE_THROW(std::runtime_error,"Decimal(String)","real literal \""<<str<<"\" has more than one decimal point.");
            }
            else {
                found_decimal_point=true;
            }
        } else if(c<'0' || c>'9') {
            ARIADNE_THROW(std::runtime_error,"Decimal(String)","invalid symbol '"<<c<<"' in string literal \""<<str<<"\"");
        }
        ++c_ptr;
    }
    if(!found_decimal_point) {
        _str+='.';
    }
}

std::ostream& operator<<(std::ostream& os, Decimal const& d) {
    return os << d._str;
}

#ifdef HAVE_GMPXX_H

Decimal::operator Rational() const {
    Rational q;
    bool decimal_point=false;
    uint decimal_places=0;
    int sign=1;
    const char* c_ptr=this->_str.c_str();
    const char& c=*c_ptr;
    if(c=='-') {
        sign=-1;
        ++c_ptr;
    } else if(c=='+') {
        ++c_ptr;
    }
    while(*c_ptr != 0) {
        const char& c=*c_ptr;
        if(c=='.') {
            if(decimal_point) {
                ARIADNE_THROW(std::runtime_error,"Decimal::operator Rational()","Invalid decimal literal: "<<this->_str<<"\" has more than one decimal point.");
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
            ARIADNE_THROW(std::runtime_error,"Decimal::operator Rational()","invalid symbol '"<<c<<"' in decimal literal \""<<this->_str<<"\"");
        }
        ++c_ptr;
    }
    for(uint i=0; i!=decimal_places; ++i) {
        q=q/10;
    }
    return sign*q;
}

#endif // HAVE_GMPXX_H


} // namespace Ariadne

