/***************************************************************************
 *            taylor_variable.cc
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
 
#include "numeric.h"
#include "sparse_differential.h"
#include "taylor_variable.h"

namespace Ariadne {

inline void acc(Interval& e, Float& x, const Float& y) {
    Float z=x;
    x+=y;
    e+=(Interval(y)+Interval(z)-Interval(x));
}

inline void acc(Interval& e, Float& x, const Interval& y) {
    Float z=x;
    x+=midpoint(y);
    e+=(y+Interval(z)-Interval(x));
}

inline void acc(Interval& e, Float& x, const Float& y, const Float& z) {
    Float xold=x;
    x+=y*z;
    e+=Interval(xold)+Interval(y)*Interval(z)-Interval(x);
}


TaylorVariable&
TaylorVariable::operator+=(const TaylorVariable& x)
{
    for(const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        const Float& xref=iter->second;
        Float& tref=(*this)[iter->first];
        acc(this->_error,tref,xref);
    }
    return *this;
}

TaylorVariable&
TaylorVariable::operator-=(const TaylorVariable& x)
{
    for(const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        Float xneg=-iter->second;
        Float& tref=(*this)[iter->first];
        acc(this->_error,tref,xneg);
    }
    return *this;
}


TaylorVariable&
TaylorVariable::operator+=(const Float& x)
{
    return (*this)+=Interval(x);
}

TaylorVariable&
TaylorVariable::operator-=(const Float& x)
{
    return (*this)-=Interval(x);
}


TaylorVariable&
TaylorVariable::operator+=(const Interval& x)
{
    acc(this->_error,const_cast<Float&>(this->_expansion.value()),x);
    return *this;
}

TaylorVariable&
TaylorVariable::operator-=(const Interval& x)
{
    acc(this->_error,const_cast<Float&>(this->_expansion.value()),-x);
    return *this;
}

TaylorVariable&
TaylorVariable::operator*=(const Interval& x)
{
    if(x==0) {
        this->_expansion*=0;
        this->_error=0;
    } else {
        this->_error*=mag(x);
        for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
            Float& tref=const_cast<Float&>(iter->second);
            Interval s=tref*x;
            tref=midpoint(s);
            this->_error+=(s-midpoint(s));
        }
    }
    return *this;
}

TaylorVariable&
TaylorVariable::operator/=(const Interval& x)
{
    return (*this) *= Ariadne::rec(x);
}


TaylorVariable 
operator+(const TaylorVariable& x) {
    return x; 
}

TaylorVariable 
operator-(const TaylorVariable& x) {
    return neg(x); 
}

TaylorVariable 
operator+(const TaylorVariable& x, const TaylorVariable& y) {
    return add(x,y); 
}

TaylorVariable 
operator-(const TaylorVariable& x, const TaylorVariable& y) {
    return add(x,neg(y)); 
}

TaylorVariable 
operator*(const TaylorVariable& x, const TaylorVariable& y) {
    return mul(x,y);
}

TaylorVariable 
operator/(const TaylorVariable& x, const TaylorVariable& y) {
    return mul(x,rec(y));
}



TaylorVariable 
operator+(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); r+=Interval(c); return r;
}

TaylorVariable 
operator-(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); r-=Interval(c); return r;
}

TaylorVariable 
operator*(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); r*=Interval(c); return r;
}

TaylorVariable 
operator/(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); r/=Interval(c); return r;
}



TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable abs(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y) {
    TaylorVariable r=x;
    r+=y;
    return r;
}

Interval sum(const TaylorVariable& x) {
    typedef TaylorVariable::const_iterator const_iterator;
    Interval r=0;
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r+=xiter->second;
    }
    return r;
}

TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    typedef TaylorVariable::const_iterator const_iterator;
    TaylorVariable r(x.argument_size());
    Interval& e=r.error();
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            acc(e,r[xiter->first+yiter->first],xiter->second,yiter->second);
        }
    }
    e += x.error() * sum(y) + sum(x) * y.error() + x.error() * y.error();
    return r;
}

TaylorVariable neg(const TaylorVariable& x) {
    return TaylorVariable(-x.expansion(),-x.error());
}

struct TaylorSeries {
    TaylorSeries(uint d) : expansion(d), error(0) { }
    Float& operator[](uint i) { return expansion[i]; }
    array<Float> expansion;
    Interval error;
};

Interval 
TaylorVariable::range() const {
    Interval r=this->error();
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=iter->second*Interval(-1,1);
    }
    return r;
}
 

TaylorVariable 
compose(const TaylorSeries& ts, const TaylorVariable& tv)
{
    Float& vref=const_cast<Float&>(tv.expansion().value());
    Float vtmp=vref; 
    vref=0.0;
    TaylorVariable r(tv);
    r+=ts.expansion[ts.expansion.size()-1];
    for(uint i=1; i!=ts.expansion.size(); ++i) {
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
    }
    vref=vtmp;
    return r;
}

TaylorVariable rec(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable pow(const TaylorVariable& x, int n) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable sqrt(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable exp(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable log(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable sin(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable cos(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorVariable tan(const TaylorVariable& x) {
    ARIADNE_NOT_IMPLEMENTED;
}



std::ostream& 
operator<<(std::ostream& os, const TaylorVariable& tv) {
    return os << "TaylorVariable( expansion=" << tv.expansion() << ", error=" << tv.error() << " )";
}

 
} //namespace Ariadne


