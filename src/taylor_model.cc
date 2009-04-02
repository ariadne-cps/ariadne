/***************************************************************************
 *            taylor_model.cc
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include <iomanip>

#include "config.h"
#include "rounding.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "expansion.h"
#include "series.h"
#include "differential.h"
#include "taylor_model.h"
#include "exceptions.h"

namespace Ariadne {


const double em=2.2204460492503131e-16;
const double ec=em/2;

double TaylorModel::_default_sweep_threshold=1e-18;
uint TaylorModel::_default_maximum_degree=16;



TaylorModel::TaylorModel()
    : _expansion(), _error(),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(0,_default_maximum_degree)
{ }

TaylorModel::TaylorModel(uint as)
    : _expansion(as), _error(0),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(as,_default_maximum_degree)
{
}

TaylorModel::TaylorModel(const std::map<MultiIndex,Float>& d, const Float& e)
    : _expansion(d), _error(e),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(d.begin()->first.size(),_default_maximum_degree)
{
    ARIADNE_ASSERT(!d.empty());
    if(e<0) { std::cerr<<e<<std::endl; } ARIADNE_ASSERT(this->_error>=0);
}

TaylorModel::TaylorModel(const Expansion<Float>& f, const Float& e)
    : _expansion(f), _error(e),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(f.begin()->first.size(),_default_maximum_degree)
{
}


void
TaylorModel::swap(TaylorModel& tm)
{
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
}



TaylorModel&
TaylorModel::operator=(const Float& c)
{
    this->_expansion.clear();
    if(c!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),c);
    }
    this->_error=0;
    return *this;
}

TaylorModel&
TaylorModel::operator=(const Interval& c)
{
    this->_expansion.clear();
    Float m=c.midpoint();
    if(m!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    this->_error=c.radius();
    return *this;
}


namespace { // Internal code for arithmetic

void scal(TaylorModel& r, const Float& c)
{
    //std::cerr<<"TaylorModel::scal(Float c) c="<<c<<std::endl;
    Float& re=r.error();

    // Shortcuts for special cases
    if(c==0.0) {
        r.expansion().clear();
        re=0;
        return;
    }
    if(c==1.0) {
        return;
    }
    if(c==0.5 || c==2.0 || c==-2.0 || c==-1.0 || c==-0.5) {
        // Operation can be performed exactly
        for(TaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
            iter->second*=c;
        }
        re*=abs(c);
        return;
    }

    // General case with error analysis
    set_rounding_upward();
    volatile Float te=0; // Twice the maximum accumulated error
    for(TaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        volatile Float u=riter->second*c;
        volatile Float t=-riter->second;
        volatile Float ml=t*c;
        te+=(u+ml);
    }
    re*=abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    for(TaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->second*=c;
    }
    return;
}



void scal(TaylorModel& r, const Interval& c)
{
    //std::cerr<<"TaylorModel::scal(Interval c) c="<<c<<std::endl;
    Float& re=r.error();
    set_rounding_upward();
    volatile Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    for(TaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->second;
        volatile Float mrv=-rv;
        if(rv>=0) {
            u=rv*c.u;
            ml=mrv*c.l;
        } else {
            u=rv*c.l;
            ml=mrv*c.u;
        }
        te+=(u+ml);
    }
    re*=mag(c);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=(c.u+c.l)/2;
    for(TaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->second*=m;
    }
    return;
}



void acc(TaylorModel& r, const Float& c)
{
    // Compute self+=c
    if(c==0) { return; }
    Float& rv=r.value();
    Float& re=r.error();
    set_rounding_upward();
    volatile Float rvu=rv+c;
    volatile Float mrvl=(-rv)-c;
    //std::cerr<<"re="<<re.u<<" ";
    re+=(rvu+mrvl)/2;
    //std::cerr<<"nre="<<re.u<<"\n";
    set_rounding_to_nearest();
    rv+=c;
    return;
}



void acc(TaylorModel& r, const Interval& c)
{
    // Compute self+=c
    if(c.l==-c.u) {
        set_rounding_upward();
        r.error()+=c.upper();
        set_rounding_to_nearest();
        return;
    }

    Float& rv=r.value();
    Float& re=r.error();
    set_rounding_upward();
    volatile Float rvu=rv+c.u;
    volatile Float mrvl=(-rv)-c.l;
    re+=(rvu+mrvl)/2;
    set_rounding_to_nearest();
    volatile Float m=(c.u+c.l)/2;
    rv+=m;
    return;
}


void acc1(TaylorModel& x, const TaylorModel& y)
{
    // Compute self+=y
    Float& xe=x.error();

    set_rounding_upward();
    Float te=0;
    for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        const Float& yv=yiter->second;;
        Float& xv=x[yiter->first];
        if(xv!=0) {
            volatile Float u=xv+yv;
            volatile Float t=-xv;
            volatile Float ml=t-yv;
            te+=(u+ml);
            //std::cerr<<" xv="<<xv<<" yv="<<yv<<" u="<<u<<" ml="<<ml<<" d="<<(u+ml)<<" te="<<te<<"\n";
        }
    }
    xe+=te/2;
    xe+=y.error();

    set_rounding_to_nearest();
    for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        x[yiter->first]+=yiter->second;
    }
    x.expansion().cleanup();
    return;
}


void acc2(TaylorModel& x, const TaylorModel& y)
{
    // Compute x+=y
    TaylorModel r(x.argument_size());

    set_rounding_upward();
    Float te=0.0;
    TaylorModel::const_iterator xiter=x.begin();
    TaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->first<yiter->first) {
            ++xiter;
        } else if(yiter->first<xiter->first) {
            ++yiter;
        } else {
            const Float& xv=xiter->second;
            const Float& yv=yiter->second;
            volatile Float u=xv+yv;
            volatile Float t=-xv;
            volatile Float ml=t-yv;
            te+=(u+ml);
            ++xiter; ++yiter;
        }
    }
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);

    set_rounding_to_nearest();
    xiter=x.begin();
    yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->first<yiter->first) {
            r.expansion().append(xiter->first,xiter->second);
            ++xiter;
        } else if(yiter->first<xiter->first) {
            r.expansion().append(yiter->first,yiter->second);
            ++yiter;
        } else {
            Float c=xiter->second+yiter->second;
            if(c!=0) { r.expansion().append(xiter->first,c); }
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x.end()) {
        r.expansion().append(xiter->first,xiter->second);
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().append(yiter->first,yiter->second);
        ++yiter;
    }

    x.expansion().swap(r.expansion());
    x.error()=r.error();

    return;
}


inline void acc(TaylorModel& r, const TaylorModel& x)
{
    acc2(r,x);
}


struct Ivl { double u; double ml; };

void mulacc1(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    // Compute self+=x*y
    typedef TaylorModel::const_iterator const_iterator;
    typedef std::map<MultiIndex,Ivl>::const_iterator ivl_const_iterator;
    Float& re=r.error();
    std::map<MultiIndex,Ivl> z;

    set_rounding_upward();
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->second;;
            const Float& yv=yiter->second;;
            Ivl& zv=z[xiter->first+yiter->first];
            zv.u+=xv*yv;
            volatile double t=-xv;
            zv.ml+=t*yv;
        }
    }

    for(const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        Ivl& zv=z[riter->first];
        const Float& rv=riter->second;
        zv.u+=rv; zv.ml-=rv;
    }

    volatile Float te=0;
    for(ivl_const_iterator ziter=z.begin(); ziter!=z.end(); ++ziter) {
        const Ivl& zv=ziter->second;
        te+=(zv.u+zv.ml);
    }
    te/=2;

    Float xs=0;
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->second);
    }

    Float ys=0;
    for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->second);
    }

    const Float& xe=x.error();
    const Float& ye=y.error();

    re+=xs*ye+ys*xe+te+xe*ye;

    set_rounding_to_nearest();
    for(ivl_const_iterator ziter=z.begin(); ziter!=z.end(); ++ziter) {
        const Ivl& zv=ziter->second;
        r[ziter->first]=(zv.u-zv.ml)/2;
    }

    return;
}


void mulacc2(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    ARIADNE_ASSERT(r.begin()->first.degree()==0);
    TaylorModel t(x.argument_size());
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        set_rounding_upward();
        double te=0.0;
        for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->second;
            const Float& yv=yiter->second;
            volatile Float u=xv*yv;
            volatile Float ml=-xv; ml=ml*yv;
            te+=(u+ml);
        }
        t.error()=te/2;
        set_rounding_to_nearest();
        for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            t.expansion().append(xiter->first,yiter->first,xiter->second*yiter->second);
        }
        r+=t;
        //std::cerr<<"  t="<<t<<"\n  r="<<r<<std::endl;
        t.expansion().clear();
        t.error()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->second);
    }

    Float ys=0;
    for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->second);
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    r.clean();

    set_rounding_to_nearest();
    return;
}


inline void mulacc(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    mulacc2(r,x,y);
}

} // namespace


///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


TaylorModel&
TaylorModel::clean()
{
    return this->clean(this->_maximum_degree,this->_sweep_threshold);
}

TaylorModel&
TaylorModel::clean(uint deg, double eps)
{
    TaylorModel r(this->argument_size());
    r.error()=this->error();

    set_rounding_upward();
    for(TaylorModel::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->first.degree()==0 || (iter->first.degree() <= deg && abs(iter->second)>=eps)) {
            r.expansion().append(iter->first,iter->second);
        } else {
            r.error()+=abs(iter->second);
        }
    }
    set_rounding_to_nearest();

    this->expansion().swap(r.expansion());
    this->error()=r.error();

    return *this;
}

TaylorModel&
TaylorModel::sweep()
{
    return this->sweep(this->_sweep_threshold);
}

TaylorModel&
TaylorModel::sweep(double m)
{
    this->clean(255,m);
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        Float a=abs(iter->second);
        if(abs(iter->second)<=m) {
            this->_error+=a;
           _expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    return *this;
}

TaylorModel&
TaylorModel::truncate()
{
    return this->truncate(this->_maximum_degree);
}

TaylorModel&
TaylorModel::truncate(uint d)
{
    this->clean(d,0.0);
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(iter->first.degree()>d) {
            e+=abs(iter->second);
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error+=e;
    set_rounding_to_nearest();
    return *this;
}

TaylorModel&
TaylorModel::truncate(const MultiIndex& a)
{
    ARIADNE_ASSERT(a.size()==this->argument_size());
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        bool erase=false;
        const MultiIndex& b=iter->first;
        for(uint i=0; i!=a.size(); ++i) {
            if(a[i] < b[i]) {
                erase=true;
                break;
            }
        }
        if(erase) {
            e+=abs(iter->second);
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error+=e;
    set_rounding_to_nearest();
    return *this;
}

TaylorModel&
TaylorModel::truncate(const MultiIndexBound& b)
{
    ARIADNE_ASSERT(b.size()==this->argument_size());
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(!(iter->first<=b)) {
            e+=abs(iter->second);
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error+=e;
    set_rounding_to_nearest();
    return *this;
}

TaylorModel&
TaylorModel::clobber()
{
    this->_error=0;
    return *this;
}

TaylorModel&
TaylorModel::clobber(uint o)
{
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(iter->first.degree()>o) {
           _expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error=0;
    return *this;
}

TaylorModel&
TaylorModel::clobber(uint so, uint to)
{
    uint n=this->argument_size()-1;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        const MultiIndex& a=iter->first;
        if(a[n]>to || a.degree()>so+a[n]) {
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error=0;
    return *this;
}


///////////////////////////////////////////////////////////////////////////////

// Arithmetic operators

TaylorModel&
operator+=(TaylorModel& x, const TaylorModel& y)
{
    if(x.argument_size()==0) { x=TaylorModel::zero(y.argument_size()); }
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    acc(x,y); return x;

}

TaylorModel&
operator-=(TaylorModel& x, const TaylorModel& y)
{
    if(x.argument_size()==0) { x=TaylorModel::zero(y.argument_size()); }
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    acc(x,neg(y)); return x;
}


TaylorModel&
operator+=(TaylorModel& x, const Float& c)
{
    acc(x,c); return x;
}

TaylorModel&
operator-=(TaylorModel& x, const Float& c)
{
    acc(x,-c); return x;
}


TaylorModel&
operator+=(TaylorModel& x, const Interval& c)
{
    acc(x,c); return x;
}

TaylorModel&
operator-=(TaylorModel& x, const Interval& c)
{
    acc(x,-c); return x;
}

TaylorModel&
operator*=(TaylorModel& x, const Float& c)
{
    scal(x,c); return x;
}

TaylorModel&
operator*=(TaylorModel& x, const Interval& c)
{
    scal(x,c); return x;
}


TaylorModel&
operator/=(TaylorModel& x, const Float& c)
{
    if(c==0) {
        ARIADNE_THROW(DivideByZeroException,"operator/=(TaylorModel x,Float c)","x="<<x<<" c="<<c);
    }
    scal(x,Interval(1.0)/c); return x;
}


TaylorModel&
operator/=(TaylorModel& x, const Interval& c)
{
    if(c.upper()>=0 && c.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"operator/=(TaylorModel x,Interval c)","x="<<x<<" c="<<c);
    }
    scal(x,1.0/c); return x;
}




TaylorModel
operator+(const TaylorModel& x) {
    return x;
}

TaylorModel
operator-(const TaylorModel& x) {
    return neg(x);
}


TaylorModel
operator+(const TaylorModel& x, const TaylorModel& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    TaylorModel r(x); acc(r,y); return r;
}

TaylorModel
operator-(const TaylorModel& x, const TaylorModel& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    TaylorModel r=neg(y); acc(r,x); return r;
}

TaylorModel
operator*(const TaylorModel& x, const TaylorModel& y) {
    if(x.argument_size()!=y.argument_size()) { std::cerr<<"operator*(TaylorModel x, TaylorModel y)\n  x="<<x<<" y="<<y<<"\n"; }
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    TaylorModel r(x.argument_size()); r.set_value(0.0); mulacc(r,x,y); return r;
}

TaylorModel
operator/(const TaylorModel& x, const TaylorModel& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    TaylorModel r(x.argument_size()); r.set_value(0.0); mulacc(r,x,rec(y)); return r;
}



TaylorModel
operator+(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); acc(r,c); return r;
}

TaylorModel
operator-(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); acc(r,-c); return r;
}

TaylorModel
operator*(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); scal(r,c); return r;
}

TaylorModel
operator/(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); scal(r,Interval(1)/c); return r;
}

TaylorModel
operator+(const Float& c, const TaylorModel& x) {
    TaylorModel r(x); acc(r,c); return r;
}

TaylorModel
operator-(const Float& c, const TaylorModel& x) {
    TaylorModel r=neg(x); acc(r,c); return r;
}

TaylorModel
operator*(const Float& c, const TaylorModel& x) {
    TaylorModel r(x); scal(r,c); return r;
}

TaylorModel
operator/(const Float& c, const TaylorModel& x) {
    TaylorModel r(x); scal(r,Interval(1)/c); return r;
}



TaylorModel
operator+(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); acc(r,c); return r;
}

TaylorModel
operator-(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); acc(r,-c); return r;
}

TaylorModel
operator*(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); scal(r,c); return r;
}

TaylorModel
operator/(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); scal(r,1/c); return r;
}

TaylorModel
operator+(const Interval& c, const TaylorModel& x) {
    TaylorModel r(x); acc(r,c); return r;
}

TaylorModel
operator-(const Interval& c, const TaylorModel& x) {
    TaylorModel r=neg(x); acc(r,c); return r;
}

TaylorModel
operator*(const Interval& c, const TaylorModel& x) {
    TaylorModel r(x); scal(r,c); return r;
}

TaylorModel
operator/(const Interval& c, const TaylorModel& x) {
    TaylorModel r=rec(x); scal(r,c); return r;
}



//////////////////////////////////////////////////////////////////////////////

// Exact functions (max, min, abs, neg) and arithmetical functions (sqr, pow)


TaylorModel max(const TaylorModel& x, const TaylorModel& y) {
    Interval xr=x.range();
    Interval yr=y.range();
    if(xr.lower()>=yr.upper()) {
        return x;
    } else if(yr.lower()>=xr.upper()) {
        return y;
    } else {
        TaylorModel z(x.argument_size());
        z.value()=max(x.value(),y.value());
        z.error()=(max(xr.u,yr.u)-max(x.value(),y.value()));
        return z;
    }
}


TaylorModel min(const TaylorModel& x, const TaylorModel& y) {
    return -max(-x,-y);
}

TaylorModel abs(const TaylorModel& x) {
    Interval xr=x.range();
    if(xr.lower()>=0.0) {
        return x;
    } else if(xr.upper()<=0.0) {
        return -x;
    } else {
        TaylorModel z(x.argument_size());
        Float xv=x.value();
        z.value()=(abs(xv));
        z.error()=(abs(xr.u)-abs(xv));
        return z;
    }

}

TaylorModel neg(const TaylorModel& x) {
    TaylorModel r(x.argument_size());
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r.expansion().append(xiter->first,-xiter->second);
    }
    r.error()=x.error();
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

TaylorModel sqr(const TaylorModel& x) {
    TaylorModel r=x*x;
    return r;
}

TaylorModel pow(const TaylorModel& x, int n) {
    TaylorModel r(x.argument_size()); r+=1;
    TaylorModel p(x);
    while(n) {
        if(n%2) { r=r*p; }
        p=sqr(p);
        n/=2;
    }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

Vector<Interval>
TaylorModel::domain() const
{
    return Vector<Interval>(this->argument_size(),Interval(-1,1));
}

Interval
TaylorModel::range() const {
    Interval r(-this->error(),+this->error());
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->first.degree()==0) {
            r+=iter->second;
        } else {
            r+=abs(iter->second)*Interval(-1,1);
        }
    }
    return r;

    /* FIXME: The following code does not work with optimisation turned on using gcc */
    set_rounding_mode(upward);
    volatile Float t=this->error();
    volatile Float v=0.0;
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->first.degree()==0) {
            v=iter->second;
        } else {
            t+=abs(iter->second);
        }
    }
    set_rounding_mode(to_nearest);
    return v+Interval(-t,t);
}


Interval
TaylorModel::evaluate(const Vector<Float>& v) const
{
    return this->evaluate(Vector<Interval>(v));
}

Interval
TaylorModel::evaluate(const Vector<Interval>& v) const
{
    ARIADNE_ASSERT(this->argument_size()==v.size());
    ARIADNE_ASSERT(subset(v,this->domain()));
    Interval r=this->error()*Interval(-1,+1);
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Interval t(iter->second);
        for(uint j=0; j!=iter->first.size(); ++j) {
            t*=pow(v[j],iter->first[j]);
        }
        r+=t;
    }
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;

namespace {

typedef Series<Interval>(*series_function_pointer)(uint,const Interval&);

class TaylorSeries {
    typedef Series<Interval>(*series_function_pointer)(uint,const Interval&);
  public:
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
    Interval e=r-c;
    //std::cerr<<"\nc="<<c<<" r="<<r<<" e="<<e<<"\n";
    //std::cerr<<"centre_series="<<centre_series<<"\nrange_series="<<range_series<<"\n";
    for(uint i=0; i!=d; ++i) {
        this->expansion[i]=midpoint(centre_series[i]);
        this->error+=(centre_series[i]-this->expansion[i])*p;
        p*=e;
    }
    //this->expansion[d]=midpoint(centre_series[d]);
    this->expansion[d]=midpoint(range_series[d]);
    this->error+=(range_series[d]-this->expansion[d])*p;
    //std::cerr<<"expansion="<<this->expansion<<"\nerror="<<this->error<<"\n";
}


std::ostream&
operator<<(std::ostream& os, const TaylorSeries& ts) {
    return os<<"TS("<<ts.expansion<<","<<ts.error<<")";
}


TaylorModel
_compose(const TaylorSeries& ts, const TaylorModel& tv, Float eps)
{
    //std::cerr<<"_compose(TaylorSeries,TaylorModel,Error)\n";
    //std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    Float& vref=const_cast<Float&>(tv.value());
    Float vtmp=vref;
    vref=0.0;
    TaylorModel r(tv.argument_size());
    r+=ts.expansion[ts.expansion.size()-1];
    for(uint i=1; i!=ts.expansion.size(); ++i) {
        //std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
        r.sweep(eps);
    }
    //std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    //std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

TaylorModel
compose(const TaylorSeries& ts, const TaylorModel& tv)
{
    return _compose(ts,tv,ec);
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
TaylorModel
_compose1(const series_function_pointer& fn, const TaylorModel& tv, Float eps)
{
    static const uint DEGREE=18;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tv.value();
    Interval r=tv.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    Float truncation_error_estimate=mag(range_series[d])*pow(mag(r-c),d);
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error_estimate
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorModel x=tv-c;
    TaylorModel res(tv.argument_size());
    res+=range_series[d];
    for(uint i=0; i!=d; ++i) {
        //std::cerr<<"i="<<i<<" r="<<res<<"\n";
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    //std::cerr<<"i="<<d<<" r="<<res<<"\n";
    return res;
}

// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. The radius of convergence of this method is still quite low,
// typically only half of the radius of convergence of the power series itself
TaylorModel
_compose2(const series_function_pointer& fn, const TaylorModel& tv, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tv.value();
    Interval r=tv.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    Float truncation_error=mag(range_series[d]-centre_series[d])*pow(mag(r-c),d);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorModel x=tv-c;
    TaylorModel res(tv.argument_size());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
TaylorModel
_compose3(const series_function_pointer& fn, const TaylorModel& tv, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tv.value();
    Interval r=tv.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    Interval se=range_series[d]-centre_series[d];
    Interval e=r-c;
    Interval p=pow(e,d-1); p.l*=-e.l; p.u*=e.u;
    //std::cerr<<"se="<<se<<" e="<<e<<" p="<<p<<std::endl;
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float truncation_error=max(se.l*p.l,se.u*p.u);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorModel x=tv;
    TaylorModel res(tv.argument_size());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


TaylorModel
_compose(const series_function_pointer& fn, const TaylorModel& tv, Float eps) {
    //std::cerr<<"_compose(SeriesFunction,TaylorModel,Error)\n";
    return _compose3(fn,tv,eps);
}


} // namespace



///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   bounded domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)

TaylorModel sqrt(const TaylorModel& x) {
    //std::cerr<<"rec(TaylorModel)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    assert(r.l>0);
    Float a=(r.l+r.u)/2;
    set_rounding_upward();
    Float eps=(r.u-r.l)/(r.u+r.l);
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=uint(log((1-eps)*x._sweep_threshold)/log(eps)+1);
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    TaylorModel y=(x/a)-1.0;
    //std::cerr<<"y="<<y<<std::endl;
    TaylorModel z(x.argument_size());
    Series<Interval> sqrt_series=Series<Interval>::sqrt(d,Interval(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(uint i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        z.sweep(); z.truncate();
        //std::cerr<<"z="<<z<<std::endl;
    }
    Float trunc_err=pow(eps,d)/(1-eps)*mag(sqrt_series[d]);
    //std::cerr<<"te="<<trunc_err<<" te*[-1,+1]="<<trunc_err*Interval(-1,1)<<std::endl;
    z.error()+=trunc_err;
    //std::cerr<<"z="<<z<<std::endl;
    Interval sqrta=sqrt(Interval(a));
    //std::cerr<<"sqrt(a)="<<sqrta<<std::endl;
    z*=sqrt(Interval(a));
    //std::cerr<<"z="<<z<<std::endl;
    return z;
}

TaylorModel rec(const TaylorModel& x) {
    //std::cerr<<"rec(TaylorModel)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    if(r.upper()>=0 && r.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"rec(TaylorModel x)","x="<<x);
    }
    Float a=(r.l+r.u)/2;
    set_rounding_upward();
    Float eps=abs((r.u-r.l)/(r.u+r.l));
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=uint(log((1-eps)*x._sweep_threshold)/log(eps))+1;
    //std::cerr<<"  x="<<x<<"\n";
    //std::cerr<<"  r="<<r<<"\n";
    //std::cerr<<"  a="<<a<<"\n";
    TaylorModel y=1-(x/a);
    //std::cerr<<"  y="<<y<<"\n";
    TaylorModel z(x.argument_size());
    z+=Float(d%2?-1:+1);
    for(uint i=0; i!=d; ++i) {
        z=1.0 + z * y;
        z.sweep(); z.truncate();
        //std::cerr<<"z.sweep_threshold="<<z.sweep_threshold()<<"\n";
    }
    //std::cerr<<"  z="<<z<<"\n";
    Float te=pow(eps,d)/(1-eps);
    //std::cerr<<"  te="<<te<<"\n";
    set_rounding_upward();
    Float nze=te+z.error();
    //std::cerr<<"  nze="<<nze<<"\n";

    set_rounding_to_nearest();
    z.set_error(nze);
    //std::cerr<<"  z="<<z<<"\n";
    z/=a;
    //std::cerr<<"  z="<<z<<"\n";
    return z;
}

TaylorModel log(const TaylorModel& x) {
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    assert(r.l>0);
    Float a=(r.l+r.u)/2;
    set_rounding_upward();
    Float eps=(r.u-r.l)/(r.u+r.l);
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=uint(log((1-eps)*x._sweep_threshold)/log(eps)+1);
    TaylorModel y=x/a-1;
    TaylorModel z(x.argument_size());
    z+=Float(d%2?-1:+1)/d;
    for(uint i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
        z.sweep(); z.truncate();
    }
    z=z*y;
    z.sweep();
    Float trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(Interval(a))+trunc_err*Interval(-1,1);
}

TaylorModel exp(const TaylorModel& x) {
    static const uint DEG=18;
    return _compose(&Series<Interval>::exp,x,x._sweep_threshold);
    return compose(TaylorSeries(DEG,&Series<Interval>::exp,
                                x.value(),x.range()),x);
}

TaylorModel sin(const TaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::sin,
                                x.value(),x.range()),x);
}

TaylorModel cos(const TaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::cos,
                                x.value(),x.range()),x);
}

TaylorModel tan(const TaylorModel& x) {
    return sin(x)*rec(cos(x));
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::tan,
                                x.value(),x.range()),x);
}

TaylorModel asin(const TaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::asin,
                                x.value(),x.range()),x);
}

TaylorModel acos(const TaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::acos,
                                x.value(),x.range()),x);
}

TaylorModel atan(const TaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::atan,
                                x.value(),x.range()),x);
}


namespace {
inline int pow2(uint k) { return 1<<k; }
inline int powm1(uint k) { return (k%2) ? -1 : +1; }
}


///////////////////////////////////////////////////////////////////////////////

// Inplace function operators (rescale, restrict, antidifferentiate)


TaylorModel& TaylorModel::rescale(const Interval& ocd, const Interval& ncd)
{
    TaylorModel& x=*this;

    const double a=ocd.l; const double b=ocd.u;
    const double c=ncd.l; const double d=ncd.u;

    // Scale the interval [a,b] onto [c,d]
    // The function is given by x:-> alpha*x+beta where
    //alpha=(d-c)/(b-a) and beta=(cb-ad)/(b-a)
    Interval tmp=1.0/sub_ivl(b,a);
    Interval alpha=sub_ivl(d,c)*tmp;
    Interval beta=(mul_ivl(c,b)-mul_ivl(a,d))*tmp;

    x*=alpha;
    x+=beta;

    return x;
}

TaylorModel& TaylorModel::restrict(const Vector<Interval>& nd)
{
    TaylorModel& x=*this;
    ARIADNE_ASSERT(x.argument_size()==nd.size());
    const uint as=x.argument_size();
    const uint d=x.expansion().degree();
    if(d==0) { return x; }

    array< array<Interval> > sf(as,array<Interval>(d+1u));
    for(uint j=0; j!=as; ++j) {
        sf[j][0]=1.0;
        sf[j][1]=rad_ivl(nd[j]);
        for(uint k=2; k<=d; ++k) {
            sf[j][k]=sf[j][k-1]*sf[j][1];
        }
    }

    // TODO: separate into roundoff computation and value computation
    for(iterator iter=x.begin(); iter!=x.end(); ++iter) {
        for(uint j=0; j!=as; ++j) {
            Float& c=iter->second;
            Interval ci=c;
            for(uint j=0; j!=as; ++j) {
                ci*=sf[j][iter->first[j]];
            }
            iter->second=ci.midpoint();
            x._error=add_up(x._error,mag(ci-iter->second));
        }
    }

    return x;
}


TaylorModel& TaylorModel::antidifferentiate(uint k)
{
    TaylorModel& x=*this;
    ARIADNE_ASSERT(k<x.argument_size());

    Float& xe=x.error();
    volatile double ml,u;
    set_rounding_upward();
    Float te=0; // Twice the maximum accumulated error
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const uint c=xiter->first[k]+1;
        const double& xv=xiter->second;
        volatile double mxv=-xv;
        if(xv>=0) {
            u=xv/c;
            ml=mxv/c;
        } else {
            u=xv/c;
            ml=mxv/c;
        }
        te+=(u+ml);
    }
    xe+=te/2;

    set_rounding_to_nearest();
    for(TaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const MultiIndex& aconst=xiter->first;
        MultiIndex& a=const_cast<MultiIndex&>(aconst);
        ++a[k];
        const uint c=xiter->first[k];
        xiter->second/=c;
    }

    return x;
}


///////////////////////////////////////////////////////////////////////////////

// Scalar function operators (evaluate, split, unscale, embed, restrict)
// and predicates (refines)

Interval
evaluate(const TaylorModel& tv, const Vector<Interval>& x)
{
    return tv.evaluate(x);
}

Vector<Interval>
evaluate(const Vector<TaylorModel>& tv, const Vector<Interval>& x)
{
    Vector<Interval> r(tv.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}


pair<TaylorModel,TaylorModel>
split(const TaylorModel& tv, uint j)
{
    // TODO: improve efficiency of implementation
    uint as=tv.argument_size();
    Vector<TaylorModel> s=TaylorModel::variables(as);
    s[j]=TaylorModel::scaling(as,j,Interval(-1,0));
    TaylorModel r1=compose(tv,s);
    s[j]=TaylorModel::scaling(as,j,Interval(0,+1));
    TaylorModel r2=compose(tv,s);
    return make_pair(r1,r2);
}

TaylorModel
split(const TaylorModel& tv, uint j, bool b)
{
    uint as=tv.argument_size();
    Interval domj=(b?Interval(0,+1):Interval(-1,0));
    Vector<TaylorModel> s=TaylorModel::variables(as);
    s[j]=TaylorModel::scaling(as,j,domj);
    TaylorModel r=compose(tv,s);
    return r;
}


TaylorModel
unscale(const TaylorModel& tv, const Interval& ivl)
{
    // Scale tv so that the interval ivl maps into [-1,1]
    // The result is given by  (tv-c)*s where c is the centre
    // and s the reciprocal of the radius of ivl
    const Float& l=ivl.l;
    const Float& u=ivl.u;

    TaylorModel r=tv;
    Interval c=Interval(l/2)+Interval(u/2);
    Interval s=2/(Interval(u)-Interval(l));
    r-=c;
    r*=s;

    return r;
}


TaylorModel
rescale(const TaylorModel& tv, const Interval& ivl)
{
    // Scale tv so that the interval [-1,1] maps into ivl
    // The result is given by  (tv*s)+c where c is the centre
    // and r the radius of ivl
    const Float& l=ivl.l;
    const Float& u=ivl.u;

    TaylorModel r=tv;
    Interval c=add_ivl(l/2,u/2);
    Interval s=sub_ivl(u/2,l/2);
    r*=s;
    r+=c;

    return r;
}



Float norm(const TaylorModel& h) {
    set_rounding_mode(upward);
    Float r=h.error();
    for(TaylorModel::const_iterator iter=h.begin(); iter!=h.end(); ++iter) {
        r+=abs(iter->second);
    }
    set_rounding_mode(to_nearest);
    return r;
}

bool
refines(const TaylorModel& tv1, const TaylorModel& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    TaylorModel d=tv2;
    d.error()=0.0;
    d-=tv1;
    return norm(d) <= tv2.error();
}


bool
disjoint(const TaylorModel& tv1, const TaylorModel& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    Float tv1e=tv1.error();
    Float tv2e=tv2.error();
    const_cast<TaylorModel&>(tv1).error()=0.0;
    const_cast<TaylorModel&>(tv2).error()=0.0;
    TaylorModel d=tv2-tv1;
    const_cast<TaylorModel&>(tv1).error()=tv1e;
    const_cast<TaylorModel&>(tv2).error()=tv2e;
    //std::cerr<<"\ntv1="<<tv1<<"\ntv2="<<tv2<<"\nd="<<d
    //         <<"\n|d|="<<norm(d)<<" |e|="<<add_up(tv1.error(),tv2.error())<<"\n\n";
    return norm(d)>add_up(tv1.error(),tv2.error());
}


TaylorModel embed(const TaylorModel& x, uint as)
{
    return TaylorModel(embed(0u,x.expansion(),as),x.error());
}

TaylorModel embed(uint as, const TaylorModel& x)
{
    return TaylorModel(embed(as,x.expansion(),0u),x.error());
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

std::ostream&
operator<<(std::ostream& os, const TaylorModel& tv) {
    //os << "TaylorModel";
    return os << "(" << tv.expansion() << "," << tv.error() << ")";
}


///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


Vector<TaylorModel> TaylorModel::zeros(uint rs, uint as)
{
    Vector<TaylorModel> result(rs);
    for(uint i=0; i!=rs; ++i) {
        result[i]=TaylorModel::zero(as);
    }
    return result;
}

Vector<TaylorModel> TaylorModel::constants(uint as, const Vector<Float>& c)
{
    Vector<TaylorModel> result(c.size());
    for(uint i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel::constant(as,c[i]);
    }
    return result;
}

Vector<TaylorModel> TaylorModel::constants(uint as, const Vector<Interval>& c)
{
    Vector<TaylorModel> result(c.size());
    for(uint i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel::constant(as,c[i]);
    }
    return result;
}

Vector<TaylorModel> TaylorModel::variables(uint as)
{
    Vector<TaylorModel> result(as);
    for(uint i=0; i!=as; ++i) { result[i]=TaylorModel::variable(as,i); }
    return result;
}

Vector<TaylorModel> TaylorModel::scalings(const Vector<Interval>& d)
{
    Vector<TaylorModel> result(d.size());
    for(uint i=0; i!=d.size(); ++i) {
        result[i]=TaylorModel::scaling(d.size(),i,d[i]);
    }
    return result;
}

Vector<TaylorModel> TaylorModel::unscalings(const Vector<Interval>& cd)
{
    Vector<TaylorModel> result(cd.size());
    for(uint i=0; i!=cd.size(); ++i) {
        result[i]=TaylorModel::unscaling(cd.size(),i,cd[i]);
    }
    return result;
}

Vector<TaylorModel> TaylorModel::rescalings(const Vector<Interval>& cd,const Vector<Interval>& d)
{
    ARIADNE_ASSERT(cd.size()==d.size());
    Vector<TaylorModel> result(d.size());
    for(uint i=0; i!=d.size(); ++i) {
        result[i]=TaylorModel::rescaling(d.size(),i,cd[i],d[i]);
    }
    return result;
}


///////////////////////////////////////////////////////////////////////////////

// Vector-valued versions of scalar operators




Vector<TaylorModel> embed(const Vector<TaylorModel>& x, uint as)
{
    Vector<TaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=embed(x[i],as);
    }
    return r;
}

Vector<TaylorModel> embed(uint as, const Vector<TaylorModel>& x)
{
    Vector<TaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=embed(as,x[i]);
    }
    return r;
}

Vector<TaylorModel> combine(const Vector<TaylorModel>& x1, const Vector<TaylorModel>& x2) {
    return join(embed(x1,x2[0].argument_size()),embed(x1[0].argument_size(),x2));
}

Vector<TaylorModel> combine(const Vector<TaylorModel>& x1, const TaylorModel& x2) {
    return join(embed(x1,x2.argument_size()),embed(x1[0].argument_size(),x2));
}

bool
refines(const Vector<TaylorModel>& tv1, const Vector<TaylorModel>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(uint i=0; i!=tv1.size(); ++i) {
        if(!refines(tv1[i],tv2[i])) { return false; }
    }
    return true;
}

bool
disjoint(const Vector<TaylorModel>& tv1, const Vector<TaylorModel>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(uint i=0; i!=tv1.size(); ++i) {
        if(disjoint(tv1[i],tv2[i])) { return true; }
    }
    return false;
}

pair< Vector<TaylorModel>, Vector<TaylorModel> >
split(const Vector<TaylorModel>& tv, uint j)
{
    Vector<TaylorModel> r1(tv.size());
    Vector<TaylorModel> r2(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        make_lpair(r1[i],r2[i])=split(tv[i],j);
    }
    return make_pair(r1,r2);
}

Vector<TaylorModel>
split(const Vector<TaylorModel>& tv, uint j, bool h)
{
    Vector<TaylorModel> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=split(tv[i],j,h);
    }
    return r;
}

Vector<TaylorModel>
unscale(const Vector<TaylorModel>& tvs, const Vector<Interval>& ivls)
{
    Vector<TaylorModel> r(tvs.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=unscale(tvs[i],ivls[i]);
    }
    return r;
}

TaylorModel
rescale(const TaylorModel& tv, const Interval& old_codom, const Interval& new_codom)
{
    TaylorModel res(tv); res.rescale(old_codom,new_codom); return res;
}

TaylorModel&
rescale(TaylorModel& tv, const Interval& old_codom, const Interval& new_codom)
{
    tv.rescale(old_codom,new_codom); return tv;
}

Vector<TaylorModel>
rescale(const Vector<TaylorModel>& tvs, const Vector<Interval>& old_codom, const Vector<Interval>& new_codom)
{
    Vector<TaylorModel> r(tvs.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=rescale(tvs[i],old_codom[i],new_codom[i]);
    }
    return r;
}

Vector<TaylorModel>&
rescale(Vector<TaylorModel>& tvs, const Vector<Interval>& old_codom, const Vector<Interval>& new_codom)
{
    for(uint i=0; i!=tvs.size(); ++i) {
        rescale(tvs[i],old_codom[i],new_codom[i]);
    }
    return tvs;
}

Vector<TaylorModel>
scale(const Vector<TaylorModel>& tvs, const Vector<Interval>& new_codom)
{
    Vector<TaylorModel> r(tvs);
    for(uint i=0; i!=tvs.size(); ++i) {
        r[i].rescale(Interval(-1,+1),new_codom[i]);
    }
    return r;
}

TaylorModel antiderivative(const TaylorModel& x, uint k) {
    TaylorModel r(x);
    r.antidifferentiate(k);
    return r;
}

Vector<TaylorModel> antiderivative(const Vector<TaylorModel>& x, uint k) {
    Vector<TaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=antiderivative(x[i],k);
    }
    return r;
}

TaylorModel antiderivative(const TaylorModel& x, const Interval& dk, uint k) {
    Interval dkr=rad_ivl(dk.l,dk.u);
    TaylorModel r(x);
    r.antidifferentiate(k);
    r*=dkr;
    return r;
}

Vector<TaylorModel> antiderivative(const Vector<TaylorModel>& x, const Interval& dk, uint k) {
    Vector<TaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=antiderivative(x[i],dk,k);
    }
    return r;
}


Vector<Interval> ranges(const Vector<TaylorModel>& f) {
    Vector<Interval> r(f.size()); for(uint i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

inline Vector<TaylorModel>& clobber(Vector<TaylorModel>& h) {
    for(uint i=0; i!=h.size(); ++i) { h[i].set_error(0.0); } return h; }

inline Vector<Float> errors(const Vector<TaylorModel>& h) {
    Vector<Float> e(h.size()); for(uint i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

inline Vector<Float> norms(const Vector<TaylorModel>& h) {
    Vector<Float> r(h.size()); for(uint i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return r; }

Float norm(const Vector<TaylorModel>& h) {
    return norm(norms(h));
}


///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
Matrix<Interval>
jacobian(const Vector<TaylorModel>& f, const Vector<Interval>& x)
{
    Vector< Differential<Interval> > dx=Differential<Interval>::variables(1u,x);
    Vector< Differential<Interval> > df(f.size());
    for(uint i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    return jacobian(df);
}

// Compute the Jacobian over an arbitrary domain
Matrix<Interval>
jacobian2(const Vector<TaylorModel>& f, const Vector<Interval>& x)
{
    Vector< Differential<Interval> > dx(x.size());
    for(uint i=0; i!=x.size()-f.size(); ++i) {
        dx[i]=Differential<Interval>::constant(f.size(),1u,x[i]); }
    for(uint i=0; i!=f.size(); ++i) {
        uint j=i+(x.size()-f.size());
        dx[j]=Differential<Interval>::variable(f.size(),1u,x[j],i); }
    Vector< Differential<Interval> > df(f.size());
    for(uint i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<Interval> J=jacobian(df);
    return J;
}

/*
// Compute the Jacobian over an arbitrary domain
Matrix<Interval>
jacobian(const Vector<TaylorModel>& f, const Vector<Interval>& d)
{
    uint rs=f.size();
    uint as=f[0].argument_size();
    Matrix<Interval> J(rs,as);
    for(uint i=0; i!=rs; ++i) {
        for(TaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            const MultiIndex& a=iter->first;
            const double& x=iter->second;
            for(uint k=0; k!=as; ++k) {
                const uint c=a[k];
                if(c>0) {
                    if(iter->first.degree()==1) {
                        J[i][k]+=x;
                    }
                    else {
                        Interval p(-x,+x);
                        p*=Float(c);
                        for(uint l=0; l!=as; ++l) {
                            if(l==k) { if(a[l]>1) { p*=pow(d[l],a[l]-1); } }
                            else { if(a[k]>0) { p*=pow(d[k],a[k]); } }
                        }
                        J[i][k]+=p;
                    }
                }
            }
        }
    }
    return J;
}
*/

// Compute the Jacobian over the unit domain
Matrix<Float>
jacobian_value(const Vector<TaylorModel>& f)
{
    uint rs=f.size();
    uint as=f[0].argument_size();
    Matrix<Float> J(rs,as);
    MultiIndex a(as);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=as; ++j) {
            a[j]=1; const double x=f[i][a]; J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<Float>
jacobian2_value(const Vector<TaylorModel>& f)
{
    const uint rs=f.size();
    const uint fas=f[0].argument_size();
    const uint has=fas-rs;
    Matrix<Float> J(rs,rs);
    MultiIndex a(fas);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=rs; ++j) {
            a[has+j]=1; const double x=f[i][a]; J[i][j]=x; a[has+j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
Matrix<Interval>
jacobian_range(const Vector<TaylorModel>& f)
{
    uint rs=f.size();
    uint as=f[0].argument_size();
    Matrix<Interval> J(rs,as);
    for(uint i=0; i!=rs; ++i) {
        for(TaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=as; ++k) {
                const uint c=iter->first[k];
                if(c>0) {
                    const double& x=iter->second;
                    if(iter->first.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->first<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<Interval>
jacobian2_range(const Vector<TaylorModel>& f)
{
    uint rs=f.size();
    uint fas=f[0].argument_size();
    uint has=fas-rs;
    Matrix<Interval> J(rs,has);
    for(uint i=0; i!=rs; ++i) {
        for(TaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=rs; ++k) {
                const uint c=iter->first[has+k];
                if(c>0) {
                    const double& x=iter->second;
                    if(iter->first.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->first<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}

///////////////////////////////////////////////////////////////////////////////

// Vector operators (compose, solve, implicit, flow


TaylorModel
compose(const TaylorModel& x,
        const Vector<Interval>& d,
        const Vector<TaylorModel>& y)
{
    return compose(x,scale(y,d));
}

Vector<TaylorModel>
compose(const Vector<TaylorModel>& x,
        const Vector<Interval>& d,
        const Vector<TaylorModel>& y)
{
    return compose(x,scale(y,d));
}

TaylorModel
compose(const TaylorModel& x,
        const Vector<TaylorModel>& ys)
{
    return compose(Vector<TaylorModel>(1u,x),ys)[0];
}

Vector<TaylorModel>
compose(const Vector<TaylorModel>& x,
        const Vector<TaylorModel>& ys)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(ys.size()==x[0].argument_size());
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }
    for(uint i=1; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys[0].argument_size(),"ys="<<ys); }

    uint as=ys[0].argument_size();

    Vector<TaylorModel> r(x.size(),TaylorModel(as));
    TaylorModel t(as);
    for(uint i=0; i!=x.size(); ++i) {
        r[i].set_error(x[i].error());
        for(TaylorModel::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            t=iter->second;
            for(uint j=0; j!=iter->first.size(); ++j) {
                TaylorModel p=pow(ys[j],iter->first[j]);
                t=t*p;
            }
            r[i]+=t;
        }
    }

    return r;
}


Vector<Interval>
_solve1(const Vector<TaylorModel>& f, uint n)
{
    //std::cerr<<"solve("<<f<<")";
    Vector<Interval> x(f.size(),Interval(-1,+1));

    bool contracting=false;
    Vector<Float> midpt_x; Vector<Interval> new_x, delta_x; Matrix<Interval> Df, Dfinv;
    for(uint i=0; i!=n; ++i) {
        midpt_x=midpoint(x);
        Df=jacobian(f,x);
        Dfinv=inverse(Df);
        delta_x=prod(Dfinv,evaluate(f,Vector<Interval>(midpt_x)));
        new_x=midpt_x-delta_x;

        if(!contracting) {
            if(disjoint(x,new_x)) {
                ARIADNE_THROW(std::runtime_error,"solve(Vector<TaylorModel>)","No solution found");
            } else if(subset(new_x,x)) {
                contracting=true;
            }
        }

        //std::cerr<<" x="<<x<<" f(x)="<<evaluate(f,x)<<" f(m(x))="<<evaluate(f,midpt_x)
        //         <<" Df(x)="<<Df<<" inv(Df(x))="<<Dfinv<<" delta_x="<<delta_x<<" new_x="<<new_x<<"\n";

        x=intersection(x,new_x);
    }

    if(!contracting) {
        ARIADNE_THROW(std::runtime_error,"solve(Vector<TaylorModel>)","Could not validate solution");
    }

    return x;
}

Vector<Interval>
solve(const Vector<TaylorModel>& f)
{
    for(uint i=0; i!=f.size(); ++i) { ARIADNE_ASSERT(f[i].argument_size()==f.size()); }
    const int number_of_steps=6;
    return _solve1(f,number_of_steps);
}


TaylorModel
implicit(const TaylorModel& f) {
    return implicit(Vector<TaylorModel>(1u,f))[0];
}



// Basic non-rigorous implicit function iteration; set y = midpoint(inverse(D2f)*compose(f,join(id,h)))
Vector<TaylorModel> _implicit1(const Vector<TaylorModel>& f, uint n)
{
    //std::cerr<<__FUNCTION__<<std::endl;
    uint rs=f.size(); uint fas=f[0].argument_size(); uint has=fas-rs;

    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=TaylorModel::zeros(has,rs);
    Vector<TaylorModel> fidh,dh;

    Matrix<Float> D2finv=inverse(jacobian2_value(f));

    //std::cerr<<"\n  f="<<f<<std::endl;
    //std::cerr<<"  h="<<h<<"\n\n";
    for(uint k=0; k!=n; ++k) {
        fidh=compose(f,join(id,h));
        //std::cerr<<"    fidh="<<fidh<<"\n";
        dh=prod(D2finv,fidh);
        //std::cerr<<"    dh="<<dh<<"\n";
        h-=dh;
        clobber(h);
        //std::cerr<<"\n  h="<<h<<"\n\n";
    }
    return h;
}

// Basic implicit function iteration; set y = inverse(D2f)*compose(f,join(id,midpoint(h)))
Vector<TaylorModel> _implicit2(const Vector<TaylorModel>& f, uint n)
{
    //std::cerr<<__FUNCTION__<<std::endl;
    uint rs=f.size(); uint fas=f[0].argument_size(); uint has=fas-rs;

    Vector<Interval> domain_h(rs,Interval(-1,+1));
    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=TaylorModel::constants(has,domain_h);
    Vector<TaylorModel> fidh,dh;

    //std::cerr<<"\n  f="<<f<<std::endl;
    //std::cerr<<"  h="<<h<<"\n\n";
    for(uint k=0; k!=n; ++k) {
        Matrix<Interval> D2finv=inverse(jacobian2(f,join(domain_h,ranges(h))));
        //std::cerr<<"    D2finv="<<D2finv<<"\n";
        clobber(h);
        fidh=compose(f,join(id,h));
        //std::cerr<<"    fidh="<<fidh<<"\n";
        dh=prod(D2finv,fidh);
        //std::cerr<<"    dh="<<dh<<"\n";
        h-=dh;
        //std::cerr<<"\n  h="<<h<<"\n\n";
    }
    return h;}

// Basic implicit function iteration without updating the matrix D2f
Vector<TaylorModel> _implicit3(const Vector<TaylorModel>& f, uint n)
{
    //std::cerr<<__FUNCTION__<<std::endl;
    uint rs=f.size(); uint fas=f[0].argument_size(); uint has=fas-rs;

    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=TaylorModel::zeros(has,rs);
    Vector<TaylorModel> fidh,dh;

    Matrix<Interval> D2finv=inverse(jacobian2_range(f));

    //std::cerr<<"\n  f="<<f<<std::endl;
    //std::cerr<<"  h="<<h<<"\n\n";
    for(uint k=0; k!=n; ++k) {
        clobber(h);
        fidh=compose(f,join(id,h));
        //std::cerr<<"    fidh="<<fidh<<"\n";
        dh=prod(D2finv,fidh);
        //std::cerr<<"    dh="<<dh<<"\n";
        h-=dh;
        //std::cerr<<"\n  h="<<h<<"\n\n";
    }
    return h;
}

// In this method, we multiply the inverse derivative with f before evaluating.
// The works very poorly due to inaccuracies in the interval computations
Vector<TaylorModel> _implicit4(const Vector<TaylorModel>& f, uint n)
{
    //std::cerr<<__FUNCTION__<<std::endl;
    uint rs=f.size(); uint fas=f[0].argument_size(); uint has=fas-rs;

    Vector<Interval> domain_h(rs,Interval(-1,+1));
    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=TaylorModel::constants(has,domain_h);
    Vector<TaylorModel> fidh,dh;

    //std::cerr<<"\n  f="<<f<<std::endl;
    //std::cerr<<"  h="<<h<<"\n\n";
    for(uint k=0; k!=n; ++k) {
        Matrix<Interval> D2finv=inverse(jacobian2(f,join(domain_h,ranges(h))));
        //std::cerr<<"    D2finv="<<D2finv<<"\n";
        Vector<Interval> h_err=errors(h);
        clobber(h);
        Vector<TaylorModel> D2finvf=prod(D2finv,f);
        dh=compose(D2finvf,join(id,h));
        //std::cerr<<"    dh="<<dh<<"\n";
        //std::cerr<<"    dh_err="<<norms(dh)<<" h_err="<<h_err<<"\n";
        h-=dh;
        //std::cerr<<"\n  h="<<h<<"\n\n";
    }
    return h;
}

//Compute the implicit function by preconditioning f by the inverse
// of the Jacobian value matrix and using a Gauss-Seidel iteration scheme
// on the system y=g(y)
Vector<TaylorModel> _implicit5(const Vector<TaylorModel>& f, uint n)
{
    //std::cerr<<__FUNCTION__<<std::endl;
    uint rs=f.size(); uint fas=f[0].argument_size(); uint has=fas-rs;

    Vector<Interval> domain_h(rs,Interval(-1,+1));
    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=TaylorModel::constants(has,domain_h);
    Vector<TaylorModel> idh=join(id,h);

    // Compute the Jacobian of f with respect to the second arguments at the centre of the domain
    Matrix<Float> D2f=jacobian2_value(f);

    // Compute g=-D2finv*(f-D2f*y) = -D2finv*f-y
    Vector<TaylorModel> g=-f;
    for(uint i=0; i!=rs; ++i) {
        for(uint j=has; j!=fas; ++j) {
            g[i][MultiIndex::unit(fas,j)]=0;
        }
    }
    Matrix<Float> J=inverse(D2f);
    g=prod(J,g);
    for(uint i=0; i!=rs; ++i) {
        g[i].clean();
    }

    // Iterate h'=h(g(x,h(x)))
    Vector<TaylorModel> h_new;
    bool contracting=false;
    for(uint k=0; k!=n; ++k) {
        h_new=compose(g,join(id,h));
        if(!contracting) {
            if(refines(h_new,h)) {
                contracting=true;
            }
            if(disjoint(h_new,h)) {
                ARIADNE_THROW(ImplicitFunctionException,"implicit(Vector<TaylorModel> f)",
                              "(with f="<<f<<"): Application of Newton solver to "<<h<<" yields "<<h_new<<
                              " which is disjoint. No solution");
            }
        }
        h=h_new;
    }

    if(contracting) {
        return h;
    } else {
        ARIADNE_THROW(ImplicitFunctionException,"implicit(Vector<TaylorModel>)",
                      "Could not verify solution of implicit function to "<<h<<"\n");
    }
}



Vector<TaylorModel>
implicit(const Vector<TaylorModel>& f)
{
    // Check that the arguments are suitable
    ARIADNE_ASSERT(f.size()>0);
    for(uint i=1; i!=f.size(); ++i) { ARIADNE_ASSERT(f[i].argument_size()==f[0].argument_size()); }

    // Set some useful size constants
    const uint rs=f.size();
    const uint fas=f[0].argument_size();
    const uint has=fas-rs;

    // Check to see if a solution exists
    Matrix<Interval> D2f=jacobian2_range(f);
    Matrix<Interval> D2finv;
    try {
        D2finv=inverse(D2f);
    }
    catch(...) {
        ARIADNE_THROW(ImplicitFunctionException,
                      "implicit(Vector<TaylorModel>)",
                      "Jacobian "<<D2f<<" is not invertible");
    }

    /* The following code is unneeded; instead we use intersection in the
     * interval-Newton style contractor
    // Check to see if the implicit function iteration is a contraction in the first step
    Vector<Interval> dom(fas);
    Vector<Interval> h_codom(rs,Interval(-1,+1));
    Vector<Interval> h_dom(has,Interval(-1,+1));
    for(uint i=0; i!=has; ++i) { dom[i]=Interval(-1,+1); }
    for(uint i=0; i!=rs; ++i) { dom[i+has]=Interval(0); }
    Vector<Interval> bound = prod(D2finv,evaluate(f,dom));
    if(!subset(bound,h_codom)) {
        ARIADNE_THROW(ImplicitFunctionException,
                      "implicit(Vector<TaylorModel>)",
                      "Single-step iteration bound "<<bound<<" is not a contraction");
    }

    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=TaylorModel::constants(has,h_dom);
    Vector<TaylorModel> fidh,deltah;
    Vector<Float> h_err;
    */

    uint number_of_steps=6;
    Vector<TaylorModel> id=TaylorModel::variables(has);
    Vector<TaylorModel> h=_implicit5(f,number_of_steps);
    /*
    std::cerr<<"\n  f="<<f<<std::endl;
    std::cerr<<"  h="<<h<<"\n\n";
    for(uint k=0; k!=6; ++k) {
        Matrix<Interval> D2finv=inverse(jacobian2(f,join(h_dom,ranges(h))));
        std::cerr<<"    D2finv="<<D2finv<<"\n";
        h_err=errors(h);
        clobber(h);
        fidh=compose(f,join(id,h));
        std::cerr<<"    fidh="<<fidh<<"\n";
        deltah=prod(D2finv,fidh);
        std::cerr<<"    deltah="<<deltah<<"\n";
        std::cerr<<"      h_err="<<h_err<<" dh_nrm="<<norms(deltah)<<" dh_err="<<errors(deltah)<<"\n";
        h-=deltah;
        std::cerr<<"\n  h="<<h<<"\n\n";
    }
    */

/*    // Perform a final rigorous step to check inclusion
    Vector<TaylorModel> old_h,idh;
    old_h=h;
    id=TaylorModel::variables(has);
    idh=join(id,h);
    D2finv=inverse(jacobian2(f,ranges(idh)));
    clobber(h);
    idh=join(id,h);
    fidh=compose(f,idh);
    deltah=prod(D2finv,fidh);
    h-=deltah;

    if(!refines(h,old_h)) { std::cerr<<"Warning: h="<<h<<" does not refine "<<old_h<<"\n"; }
*/
    // Check that the result has the correct sizes.
    ARIADNE_ASSERT(h.size()==f.size());
    for(uint i=0; i!=h.size(); ++i) {
        ARIADNE_ASSERT(h[0].argument_size()+f.size()==f[i].argument_size());
    }

    return h;
}


inline Vector<TaylorModel>
_flow_step(const Vector<TaylorModel>& vf, const Vector<TaylorModel>& yz, const Interval& h_rad,
           const Vector<TaylorModel>& phi)
{
    const uint n=vf.size();
    Vector<TaylorModel> new_phi=compose(vf,phi);
    for(uint i=0; i!=n; ++i) {
        new_phi[i].antidifferentiate(n);
        new_phi[i]*=h_rad;
        new_phi[i]+=yz[i];
    }
    return new_phi;
}

Vector<TaylorModel>
flow(const Vector<TaylorModel>& vf, const Vector<Interval>& d, const Interval& h, uint order)
{
    ARIADNE_ASSERT(vf.size()>0);
    ARIADNE_ASSERT(vf.size()==d.size());
    for(uint i=0; i!=vf.size(); ++i) { ARIADNE_ASSERT(vf.size()==vf[i].argument_size()); }
    ARIADNE_ASSERT(h.l<=0.0 && h.u>=0);

    ARIADNE_ASSERT(h.l==0.0 || h.l==-h.u);

    Vector<Interval> vf_domain(vf.size(),Interval(-1,+1));
    ARIADNE_ASSERT(subset(d,vf_domain));

    // Don't check for flow bounds since this should have been done in calling function
/*
    Vector<Interval> vf_range(vf.size());
    for(uint i=0; i!=vf.size(); ++i) { vf_range[i]=vf[i].range(); }
    Vector<Interval> flow_range=d+vf_range*h;
    if(!subset(flow_range,vf_domain)) {
        ARIADNE_THROW(FlowBoundsException,"flow(Vector<TaylorModel>,Vector<Interval>,Interval,Nat)",
                      "range "<<flow_range<<" of scaled flow "<<vf<<
                      " on domain "<<d<<" for time interval "<<h<<
                      " is not a subset of the unit box");
    }
*/
    uint n=vf.size();
    Vector<TaylorModel> yz(n,TaylorModel(n+1));
    for(uint i=0; i!=yz.size(); ++i) {
        yz[i]=TaylorModel::scaling(n+1,i,d[i]); }
    Vector<TaylorModel> y(n);
    for(uint i=0; i!=y.size(); ++i) {
        y[i]=TaylorModel::constant(n+1,Interval(-1,+1)); }

    Vector<TaylorModel> old_y(n,TaylorModel(n+1));

    Float h_rad=h.upper();

    //std::cerr << "\ny[0]=" << y << std::endl << std::endl;
    //std::cerr<<"  y="<<y<<"\n";
    for(uint j=0; j!=order+18; ++j) {
        for(uint i=0; i!=n; ++i) {
            y[i].swap(old_y[i]);
        }

        y=compose(vf,y);
        for(uint i=0; i!=n; ++i) {
            y[i].antidifferentiate(n);
            y[i]*=h_rad;
            y[i]+=yz[i];
        }
    }

    //for(uint i=0; i!=n; ++i) { y[i].clobber(so,to); }
    for(uint i=0; i!=n; ++i) { y[i].sweep(0.0); }

    if(h.l==0) {
        //Vector<TaylorModel> s=TaylorModel::variables(n+1);
        //s[n]=TaylorModel::scaling(n+1,n,Interval(0,1));
        //std::cerr<<"\nsplit="<<split(y,n)<<"\nscale="<<compose(y,s)<<"\n\n";
        //y=split(y,n).second;
        //y=compose(y,s);
        y=split(y,n,true);
    }

    return y;

}





} //namespace Ariadne


