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

#include "utility/standard.h"
#include "config.h"

#include "utility/macros.h"
#include "numeric/decimal.h"
#include "numeric/rational.h"
#include "numeric/float.h"

namespace Ariadne {

Decimal operator"" _decimal(long double x)
{
    return Decimal(static_cast<double>(x));
}

Decimal operator"" _dec(long double x)
{
    return operator"" _decimal(x);
}

Decimal operator-(Decimal const& d)
{
    if(d._str[0]=='-') { return Decimal(std::string(d._str.c_str()+1)); }
    else { return Decimal(std::string("-")+d._str); }
}

Decimal::Decimal(double x)
{
    StringStream ss;
    Int s=+1;
    if(x<0) { s=-1; ss<<'-'; }
    double y=s*x;
    static const Nat sf=9;
    static const double acc=1e-9;
    static const double tol=1e-15;
    Int exp=0; double pow=1;
    while(y<1.0) { y*=10; exp-=1; pow/=10; }
    while(y>=1.0) { y/=10; exp+=1; pow*=10; }
    long int n=std::floor(y/acc+0.5);
    double re=y-n*acc;
    if(std::fabs(re)>=tol) {
        ARIADNE_THROW(std::runtime_error,"Decimal(double)","double-precision floating-point number must have a relative error of "<<tol<<" with respect to its approximation to "<<sf<<" significant figures; number "<<std::setprecision(17)<<x<<" has a relative error of "<<re<<"");
    }
    long int m=1e9; Int c=0;
    if(exp<=0) { ss << '0'; }
    if(exp<0) { ss << '.'; for(Int i=0; i!=-exp; ++i) { ss << "0"; } }
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

Decimal::Decimal(StringType str)
    : _str(str)
{
    // Parse string to ensure correctness
    Bool found_decimal_point=false;
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

OutputStream& operator<<(OutputStream& os, Decimal const& d) {
    return os << d._str;
}

Decimal::operator Rational() const {
    Rational q;
    Bool decimal_point=false;
    Nat decimal_places=0;
    Int sign=1;
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
    for(Nat i=0; i!=decimal_places; ++i) {
        q=q/10;
    }
    return sign*q;
}

Decimal::operator ExactNumber() const {
    return Rational(*this).operator ExactNumber();
}

} // namespace Ariadne

