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

inline void acc(Interval& e, const Interval& d) {
    e+=d;
}

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


void
TaylorVariable::sweep(const Float& m)
{
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        Float a=abs(iter->second);
        if(abs(iter->second)<=m) {
            this->_error+=Interval(-a,a);
            this->_expansion.data().erase(iter++);
        } else {
            ++iter;
        }
    }
}

void
TaylorVariable::clean()
{
    this->sweep(0.0);
}



TaylorVariable&
TaylorVariable::operator+=(const TaylorVariable& x)
{
    for(const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
        const Float& xref=iter->second;
        Float& tref=(*this)[iter->first];
        acc(this->_error,tref,xref);
    }
    acc(this->_error,x._error);
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
TaylorVariable::operator*=(const Float& x)
{
    if(x==0) {
        this->_expansion*=0;
        this->_error=0;
    } else {
        this->_error*=x;
        for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
            Float& tref=const_cast<Float&>(iter->second);
            Interval s=tref*Interval(x);
            tref=midpoint(s);
            this->_error+=(s-midpoint(s));
        }
    }
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
TaylorVariable::operator/=(const Float& x)
{
    this->_error/=x;
    Interval r=Interval(1)/x;
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Float& tref=const_cast<Float&>(iter->second);
        Interval s=tref*r;
        tref=midpoint(s);
        this->_error+=(s-midpoint(s));
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

TaylorVariable 
operator+(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); r+=Interval(c); return r;
}

TaylorVariable 
operator-(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); r-=Interval(c); return neg(r);
}

TaylorVariable 
operator*(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); r*=Interval(c); return r;
}

TaylorVariable 
operator/(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); r/=Interval(c); return rec(r);
}



TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    Interval xr=x.range();
    Interval yr=y.range();
    if(xr.lower()>=yr.upper()) {
        return x;
    } else if(yr.lower()>=xr.upper()) {
        return y;
    } else {
        TaylorVariable z(x.argument_size());
        z.set_value(max(x.value(),y.value()));
        z.set_error(max(xr,yr)-max(x.value(),y.value()));
        return z;
    }
}

TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y) {
    return -max(-x,-y);
}

TaylorVariable abs(const TaylorVariable& x) {
    Interval xr=x.range();    
    if(xr.lower()>=0.0) {
        return x;
    } else if(xr.upper()<=0.0) {
        return -x;
    } else {
        TaylorVariable z(x.argument_size());
        Float xv=x.value();
        z.set_value(abs(xv));
        z.set_error(abs(xr)-abs(xv));
        return z;
    }

}

Interval _sum(const TaylorVariable& x) {
    typedef TaylorVariable::const_iterator const_iterator;
    Interval r=0;
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r+=xiter->second;
    }
    return r;
}

TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y) {
    TaylorVariable r=x;
    r+=y;
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
    e += x.error() * _sum(y) + _sum(x) * y.error() + x.error() * y.error();
    return r;
}

TaylorVariable neg(const TaylorVariable& x) {
    return TaylorVariable(-x.expansion(),-x.error());
}

Vector<Interval> 
TaylorVariable::domain() const
{
    return Vector<Interval>(this->argument_size(),Interval(-1,1));
}

Interval 
TaylorVariable::range() const {
    Interval r=this->error();
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter==this->begin()) {
            r+=iter->second;
        } else {
            r+=iter->second*Interval(-1,1);
        }
    }
    return r;
}
 
Interval 
TaylorVariable::evaluate(const Vector<Interval>& v) const
{
    ARIADNE_ASSERT(subset(v,this->domain()));
    Interval r=this->error();
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Interval t=iter->second;
        for(uint j=0; j!=iter->first.size(); ++j) {
            t*=pow(v[j],iter->first[j]);
        }
        r+=t;
    }
    return r;
}

template<class X> class Series;  

struct TaylorSeries {
    typedef Series<Interval>(*series_function_pointer)(uint,const Interval&); 
    TaylorSeries(uint d) : expansion(d+1), error(0) { }
    TaylorSeries(uint degree, series_function_pointer function, 
                 const Float& centre, const Interval& domain);
    uint degree() const { return expansion.size()-1; }
    Float& operator[](uint i) { return expansion[i]; }
    array<Float> expansion;
    Interval error;
    void sweep(Float e) { 
        for(uint i=0; i<=degree(); ++i) {
            if(abs(expansion[i])<=e) { 
                error+=expansion[i]*Interval(-1,1); 
                expansion[i]=0; } } }
};

TaylorSeries::TaylorSeries(uint d, series_function_pointer fn, 
                           const Float& c, const Interval& r)
    : expansion(d+1), error(0) 
{
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);
    Interval p=1;
    for(uint i=0; i!=d; ++i) {
        this->expansion[i]=midpoint(centre_series[i]);
        this->error+=(centre_series[i]-this->expansion[i])*p;
        p*=r;
    }
    this->expansion[d]=midpoint(centre_series[d]);
    this->error+=(range_series[d]-this->expansion[d])*p;
}


std::ostream& 
operator<<(std::ostream& os, const TaylorSeries& ts) {
    return os<<"TS("<<ts.expansion<<","<<ts.error<<")";
}


TaylorVariable 
_compose(const TaylorSeries& ts, const TaylorVariable& tv, Float eps)
{
    std::cerr<<"\ncompose\n";
    std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    Float& vref=const_cast<Float&>(tv.expansion().value());
    Float vtmp=vref; 
    vref=0.0;
    TaylorVariable r(tv.argument_size());
    r+=ts.expansion[ts.expansion.size()-1];
    for(uint i=1; i!=ts.expansion.size(); ++i) {
        std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
        r.sweep(eps);
    }
    std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

TaylorVariable 
compose(const TaylorSeries& ts, const TaylorVariable& tv)
{
    return _compose(ts,tv,0.0);
}


TaylorVariable rec(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::rec,
                                x.value(),x.range()),x);
}

TaylorVariable sqr(const TaylorVariable& x) {
    TaylorVariable r=x*x;
    return r;
}

TaylorVariable pow(const TaylorVariable& x, int n) {
    TaylorVariable r(x.argument_size()); r+=1;
    TaylorVariable p(x);
    while(n) {
        if(n%2) { r=r*p; } 
        p=sqr(p);
        n/=2;
    }
    return r;
}

TaylorVariable sqrt(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::sqrt,
                                x.value(),x.range()),x);
}

TaylorVariable exp(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::exp,
                                x.value(),x.range()),x);
}

TaylorVariable log(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::log,
                                x.value(),x.range()),x);
}

TaylorVariable sin(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::sin,
                                x.value(),x.range()),x);
}

TaylorVariable cos(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::cos,
                                x.value(),x.range()),x);
}

TaylorVariable tan(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::tan,
                                x.value(),x.range()),x);
}

TaylorVariable asin(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::asin,
                                x.value(),x.range()),x);
}

TaylorVariable acos(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::acos,
                                x.value(),x.range()),x);
}

TaylorVariable atan(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::atan,
                                x.value(),x.range()),x);
}


inline int pow2(uint k) { return 1<<k; }
inline int powm1(uint k) { return (k%2) ? -1 : +1; }


pair<TaylorVariable,TaylorVariable>
split(const TaylorVariable& tv, uint j) 
{
    uint as=tv.argument_size();
    uint deg=tv.expansion().degree();
    const SparseDifferential<Float>& expansion=tv.expansion();

    MultiIndex index(as);
    Float value;
    
    SparseDifferential<Float> expansion1(as,deg);
	for(SparseDifferential<Float>::const_iterator iter=expansion.begin();
        iter!=expansion.end(); ++iter) 
    {
        index=iter->first;
        value=iter->second;
        uint k=index[j];
        for(uint l=0; l<=k; ++l) {
            index.set(j,l);
            expansion1[index] += bin(k,l) * value / pow2(k);
        }
    }

    SparseDifferential<Float> expansion2(as,deg);
	for(SparseDifferential<Float>::const_iterator iter=expansion.begin();
        iter!=expansion.end(); ++iter) 
    {
        index=iter->first;
        value=iter->second;
        uint k=index[j];
        for(uint l=0; l<=k; ++l) {
            index.set(j,l);
            // Need brackets in expression below to avoid converting negative signed
            // integer to unsigned
            expansion2[index] += powm1(k-l) * (bin(k,l) * value / pow2(k));
        }
    }
    // FIXME: Add roundoff errors when computing new expansions
    
    Interval error1=tv.error();
    Interval error2=tv.error();
        
    return make_pair(TaylorVariable(expansion1,error1),TaylorVariable(expansion2,error2));
}


std::string 
TaylorVariable::str() const
{
    std::stringstream ss;
    for(SparseDifferential<Float>::const_iterator iter=this->_expansion.begin();
        iter!=this->_expansion.end(); ++iter) 
    {
        const MultiIndex& j=iter->first;
        const Float& c=iter->second;
        ss << j << ": " << c << "\n";
    }
    ss<<"error: "<<this->_error<<"\n";
    return ss.str();
}


std::ostream& 
operator<<(std::ostream& os, const TaylorVariable& tv) {
    return os << "TaylorVariable( expansion=" << tv.expansion() << ", error=" << tv.error() << " )";
}

 
} //namespace Ariadne


