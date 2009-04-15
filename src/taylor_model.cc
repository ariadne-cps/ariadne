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
#include "taylor_series.h"
#include "exceptions.h"

namespace Ariadne {


const double em=2.2204460492503131e-16;
const double ec=em/2;

double TaylorModel::_default_sweep_threshold=1e-18;
uint TaylorModel::_default_maximum_degree=16;

TaylorModel::Accuracy::Accuracy() : _sweep_threshold(_default_sweep_threshold), _maximum_degree(_default_maximum_degree) { }
TaylorModel::Accuracy::Accuracy(double st, uint md) : _sweep_threshold(st), _maximum_degree(md) { }

TaylorModel::Accuracy
max(const TaylorModel::Accuracy& acc1, const TaylorModel::Accuracy& acc2) {
    return TaylorModel::Accuracy(std::min(acc1._sweep_threshold,acc1._sweep_threshold),
        std::max(acc1._maximum_degree,acc2._maximum_degree));
}

TaylorModel::Accuracy
min(const TaylorModel::Accuracy& acc1, const TaylorModel::Accuracy& acc2) {
    return TaylorModel::Accuracy(std::max(acc1._sweep_threshold,acc1._sweep_threshold),
        std::min(acc1._maximum_degree,acc2._maximum_degree));
}

inline bool TaylorModel::Accuracy::discard(const Float& x) const { return abs(x)<this->_sweep_threshold; }
inline bool TaylorModel::Accuracy::discard(const MultiIndex& a) const { return a.degree()>this->_maximum_degree; }
inline bool TaylorModel::Accuracy::discard(const MultiIndex& a, const Float& x) const { return this->discard(x) || this->discard(a); }

std::ostream& operator<<(std::ostream& os, const TaylorModel::Accuracy& acc) {
    return os<<"( sweep_threshold="<<acc._sweep_threshold<<", maximum_degree="<<acc._maximum_degree<<" )";
}


Vector<Interval> unscale(const Vector<Interval>& x, const Vector<Interval>& d) {
    Vector<Interval> r(x);
    for(uint i=0; i!=r.size(); ++i) {
        (r[i]-=med_ivl(d[i]))/=rad_ivl(d[i]); }
    return r;
}

TaylorModel::TaylorModel()
    : _expansion(0), _error(), _accuracy_ptr(new Accuracy())
{ }

TaylorModel::TaylorModel(uint as)
    : _expansion(as), _error(0), _accuracy_ptr(new Accuracy())
{
}

TaylorModel::TaylorModel(uint as, shared_ptr<Accuracy> acc)
    : _expansion(as), _error(0), _accuracy_ptr(acc)
{
}

TaylorModel::TaylorModel(const std::map<MultiIndex,Float>& d, const Float& e)
    : _expansion(d), _error(e), _accuracy_ptr(new Accuracy())
{
    ARIADNE_ASSERT(!d.empty());
    ARIADNE_ASSERT_MSG(this->_error>=0,"e="<<e);
}

TaylorModel::TaylorModel(const Expansion<Float>& f, const Float& e)
    : _expansion(f), _error(e), _accuracy_ptr(new Accuracy())
{
}

TaylorModel::TaylorModel(const Expansion<Float>& f, const Float& e, shared_ptr<Accuracy> a)
    : _expansion(f), _error(e), _accuracy_ptr(a)
{
}


void
TaylorModel::swap(TaylorModel& tm)
{
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
}

void
TaylorModel::clear()
{
    this->_expansion.clear();
    this->_error=0.0;
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

// Inplace negation
void _neg(TaylorModel& r)
{
    for(TaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data()=-iter->data();
    }
}

inline void _scal_exact(TaylorModel& r, const Float& c)
{
    // Operation can be performed exactly
    register Float pc=c;
    for(TaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data()*=pc;
    }
    r.error()*=abs(pc);
    return;
}


inline void _scal_approx(TaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error();
    set_rounding_upward();
    volatile Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    register Float pc=c;
    register Float mc=-c;
    for(TaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(TaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}

void _scal(TaylorModel& r, const Float& c)
{
    if(c==0) { r.expansion().clear(); r.error()=0; }
    else if(c==1) { }
    else if(c==2 || c==0.5) { _scal_exact(r,c); }
    else { _scal_approx(r,c); }
}


void _scal(TaylorModel& r, const Interval& c)
{
    //std::cerr<<"TaylorModel::scal(Interval c) c="<<c<<std::endl;
    Float& re=r.error();
    set_rounding_upward();
    volatile Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    register Float cu=c.u;
    register Float mcl=-c.l;
    for(TaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data();
        if(rv>=0) {
            u=rv*cu;
            ml=rv*mcl;
        } else {
            Float mrv=-rv;
            u=mrv*mcl;
            ml=mrv*cu;
        }
        te+=(u+ml);
    }
    re*=mag(c);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=(c.u+c.l)/2;
    for(TaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}

void _scal2(TaylorModel& r, const Interval& c)
{
    //std::cerr<<"TaylorModel::scal(Interval c) c="<<c<<std::endl;
    Float& re=r.error();
    volatile Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    volatile Float cu=c.u;
    volatile Float mcl=-c.l;
    set_rounding_to_nearest();
    Float cm=(c.l+c.u)/2;
    for(TaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        Float& rv=riter->data();
        set_rounding_upward();
        if(rv>=0) {
            u=rv*cu;
            ml=rv*mcl;
        } else {
            Float mrv=-rv;
            u=mrv*mcl;
            ml=mrv*cu;
        }
        te+=(u+ml);
        set_rounding_to_nearest();
        rv*=cm;
    }
    set_rounding_upward();
    re*=mag(c);
    re+=te/2;
    set_rounding_to_nearest();

    return;
}


struct UnitMultiIndex { uint argument_size; uint unit_index; };


inline void _incr(TaylorModel& r, const MultiIndex& a)
{
    for(TaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        static_cast<MultiIndex&>(iter->key())+=a;
    }
}

inline void _incr(TaylorModel& r, uint j)
{
    for(TaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++static_cast<MultiIndex&>(iter->key())[j];
    }
}


inline
void _acc(TaylorModel& r, const Float& c)
{
    // Compute self+=c
    if(c==0) { return; }
    if(r.expansion().empty()) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else if(r.begin()->key().degree()>0) {
        r.expansion().prepend(MultiIndex(r.argument_size()),c);
    } else {
        Float& rv=r.begin()->data();
        Float& re=r.error();
        set_rounding_upward();
        volatile Float rvu=rv+c;
        volatile Float mrvl=(-rv)-c;
        //std::cerr<<"re="<<re.u<<" ";
        re+=(rvu+mrvl)/2;
        //std::cerr<<"nre="<<re.u<<"\n";
        set_rounding_to_nearest();
        rv+=c;
    }
    return;
}



inline
void _acc(TaylorModel& r, const Interval& c)
{
    // Compute self+=c
    if(c.l==-c.u) { // The midpoint of the interval is zero, so no need to change constant term
        set_rounding_upward();
        r.error()+=c.upper();
        set_rounding_to_nearest();
        return;
    }
    if(r.expansion().empty()) { // Append a constant term zero (slightly faster than prepend)
        r.expansion().append(MultiIndex(r.argument_size()),0.0);
    } else if(r.begin()->key().degree()>0) { // Prepend a constant term zero
        r.expansion().prepend(MultiIndex(r.argument_size()),0.0);
    }

    Float& rv=r.begin()->data();
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



// Compute r=x+y, assuming r is empty
// First compute the errors and then at the end compute the expansion
// This saves rounding mode changes, but requires an extra loop
inline void _add1(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    TaylorModel::const_iterator xiter=x.begin();
    TaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            ++yiter;
        } else {
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
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
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().append(yiter->key(),yiter->data());
            ++yiter;
        } else {
            Float c=xiter->data()+yiter->data();
            if(c!=0) { r.expansion().append(xiter->key(),c); }
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x.end()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().append(yiter->key(),yiter->data());
        ++yiter;
    }

}

// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
inline void _add2(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    set_rounding_upward();
    Float te=0.0;
    TaylorModel::const_iterator xiter=x.begin();
    TaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().append(yiter->key(),yiter->data());
            ++yiter;
        } else {
            volatile Float& xv=const_cast<Float&>(static_cast<const Float&>(xiter->data()));
            volatile Float& yv=const_cast<Float&>(static_cast<const Float&>(yiter->data()));
            set_rounding_upward();
            volatile Float u=xv+yv;
            volatile Float ml=-xv; ml-=yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            Float c=xiter->data()+yiter->data();
            if(c!=0) { r.expansion().append(xiter->key(),xiter->data()+yiter->data()); }
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().append(yiter->key(),yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);
    set_rounding_to_nearest();

}

inline void _add(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    _add2(r,x,y);
}

inline void _acc(TaylorModel& r, const TaylorModel& x)
{
    TaylorModel s(r.argument_size()); _add(s,r,x); s.swap(r);
}

inline void _sub(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    TaylorModel::const_iterator xiter=x.begin();
    TaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().append(yiter->key(),yiter->data());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            volatile Float u=xv-yv;
            volatile Float t=-xv;
            volatile Float ml=t+yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().append(xiter->key(),xiter->data()+yiter->data());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().append(yiter->key(),yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);
    set_rounding_to_nearest();

}



struct Ivl { double u; double ml; };

inline void _mul1(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    // Compute r+=x*y
    typedef TaylorModel::const_iterator const_iterator;
    typedef std::map<MultiIndex,Ivl>::const_iterator ivl_const_iterator;
    Float& re=r.error();
    std::map<MultiIndex,Ivl> z;

    set_rounding_upward();
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->data();;
            const Float& yv=yiter->data();;
            Ivl& zv=z[xiter->key()+yiter->key()];
            zv.u+=xv*yv;
            volatile double t=-xv;
            zv.ml+=t*yv;
        }
    }

    for(const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        Ivl& zv=z[riter->key()];
        const Float& rv=riter->data();
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
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
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


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
void _mul2(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    TaylorModel t(x.argument_size());
    TaylorModel s(x.argument_size());
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        set_rounding_upward();
        double te=0.0;
        for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            volatile Float u=xv*yv;
            volatile Float ml=-xv; ml=ml*yv;
            te+=(u+ml);
        }
        t.error()=te/2;
        set_rounding_to_nearest();
        for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            t.expansion().append(xiter->key(),yiter->key(),xiter->data()*yiter->data());
        }
        _add2(s,r,t);
        r.expansion().swap(s.expansion());
        r.error()=s.error();
        s.expansion().clear();
        s.error()=0.0;
        t.expansion().clear();
        t.error()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    r.clean();

    set_rounding_to_nearest();
    return;
}


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Change the rounding mode to avoid iterating using opposite rounding
inline void _mul3(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    TaylorModel t(x.argument_size());
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        volatile Float pxv=xiter->data();
        volatile Float nxv=-pxv;
        volatile Float te=0.0;
        for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            set_rounding_upward();
            const Float& yv=yiter->data();
            te+=(pxv*yv)+(nxv*yv);
            set_rounding_to_nearest();
            t.expansion().append(xiter->key(),yiter->key(),xiter->data()*yiter->data());
        }
        t.error()=te/2;
        r+=t;
        //std::cerr<<"  t="<<t<<"\n  r="<<r<<std::endl;
        t.clear();
    }

    set_rounding_upward();
    Float xs=0;
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    r.clean();

    set_rounding_to_nearest();
    return;
}

inline void _add(MultiIndex& r, const MultiIndex& a1, const MultiIndex& a2) {
    for(MultiIndex::size_type j=0; j!=r.word_size(); ++j) {
        r.word_at(j)=a1.word_at(j)+a2.word_at(j);
    }
}


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
void _mul4(TaylorModel& r, const TaylorModel& x, const TaylorModel& y, const TaylorModel::Accuracy& accuracy)
{
    const uint as=r.argument_size();
    TaylorModel t(as);
    TaylorModel s(as);
    MultiIndex ta(as);
    Float tv;
    volatile Float u;
    volatile Float ml;
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        double tte=0.0; // trucation error
        double tre=0.0; // roundoff error
        const MultiIndex& xa=xiter->key();
        volatile Float& xv=const_cast<Float&>(static_cast<const Float&>(xiter->data()));
        volatile Float mxv=-xv;
        volatile Float axv=abs(xv);
        for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const MultiIndex& ya=yiter->key();
            volatile Float& yv=const_cast<Float&>(static_cast<const Float&>(yiter->data()));
            set_rounding_to_nearest();
            ta=xa+ya;
            tv=xv*yv;
            set_rounding_upward();
            if(accuracy.discard(ta,tv)) {
                tte+=axv*abs(yv);
            } else {
                t.expansion().append(ta,tv);
                u=xv*yv;
                ml=mxv*yv;
                tre+=(u+ml);
            }
        }
        t.error()=tte+(tte/2);

        _add2(s,r,t);
        r.expansion().swap(s.expansion());
        r.error()=s.error();
        s.expansion().clear();
        s.error()=0.0;
        t.expansion().clear();
        t.error()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(TaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    set_rounding_to_nearest();
    return;
}

inline void _mul(TaylorModel& r, const TaylorModel& x, const TaylorModel& y) {
    //_mul2(r,x,y);
    _mul4(r,x,y,r.accuracy());
}



} // namespace

void _mul_full(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    _mul2(r,x,y);
}

void _mul_clear(TaylorModel& r, const TaylorModel& x, const TaylorModel& y)
{
    _mul4(r,x,y,r.accuracy());
}


///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


TaylorModel&
TaylorModel::clean()
{
    return this->clean(*this->_accuracy_ptr);
}


// Clean in-place by shifting elements up to new position
inline TaylorModel& _clean1(TaylorModel& tm, const TaylorModel::Accuracy& accuracy)
{
    set_rounding_upward();
    TaylorModel::iterator curr=tm.begin();
    TaylorModel::const_iterator adv=curr;
    while(adv!=tm.end()) {
        //std::cerr<<"discard("<<adv->key()<<","<<adv->data()<<")="<<accuracy.discard(adv->key(),adv->data())<<"\n";
        if(accuracy.discard(adv->key(),adv->data())) {
            tm.error()+=abs(adv->data());
        } else {
            if(curr!=adv) { curr->key()=adv->key(); curr->data()=adv->data(); }
            ++curr;
        }
        ++adv;
    }
    set_rounding_to_nearest();
    tm.expansion().resize(curr-tm.begin());
    return tm;
}

// Clean by making a copy. Appears to be significantly slower than _clean1
inline TaylorModel& _clean2(TaylorModel& tm, const TaylorModel::Accuracy& accuracy)
{
    // Alternative code takes significantly longer
    TaylorModel r(tm.argument_size());
    r.expansion().reserve(tm.expansion().size());
    r.error()=tm.error();

    set_rounding_upward();
    for(TaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(accuracy.discard(iter->key(),iter->data())) {
            r.error()+=abs(iter->data());
        } else {
            r.expansion().append(iter->key(),iter->data());
        }
    }
    set_rounding_to_nearest();

    tm.expansion().swap(r.expansion());
    tm.error()=r.error();
    return tm;
}


TaylorModel&
TaylorModel::clean(const Accuracy& accuracy)
{
    return _clean1(*this,accuracy);
}


TaylorModel&
TaylorModel::sweep()
{
    return this->sweep(this->sweep_threshold());
}

/*
TaylorModel&
TaylorModel::sweep(double m)
{
    TaylorModel r(this->argument_size(),this->_accuracy_ptr);
    r.expansion().reserve(this->expansion().size());
    for(TaylorModel::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Float av=abs(iter->data());
        if(av>=m) {
            r.expansion().append(iter->key(),iter->data());
        } else {
            r.error()+=av;
        }
    }
    this->expansion().swap(r.expansion());
    this->error()=r.error();
    return *this;
}
*/

TaylorModel&
TaylorModel::sweep(double m)
{
    TaylorModel::const_iterator end=this->end();
    TaylorModel::const_iterator adv=this->begin();
    TaylorModel::iterator curr=this->begin();
    set_rounding_upward();
    while(adv!=end) {
        if(abs(adv->data())>=m) {
            *curr=*adv; ++curr;
        } else {
            this->error()+=adv->data();
        }
        ++adv;
    }
    this->expansion().resize(curr-this->begin());
    return *this;
}

TaylorModel&
TaylorModel::truncate()
{
    return this->truncate(this->maximum_degree());
}

/*
TaylorModel&
TaylorModel::truncate(uint d)
{
    //std::cerr<<"truncate(tm="<<*this<<", d="<<d<<")"<<std::endl;
    volatile uchar bd=d;
    TaylorModel r(this->argument_size(),this->_accuracy_ptr);
    r.expansion().reserve(this->expansion().size());
    for(TaylorModel::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        uchar adeg=iter->key().degree();
        if(adeg<=d) {
        //if(adeg<=byte_d) {
        //if(a.degree()<=byte_d) {
            
            r.expansion().append(iter->key(),iter->data());
        } else {
            r.error()+=abs(iter->data());
        }
    }
    this->expansion().swap(r.expansion());
    this->error()=r.error();
    //std::cerr<<"  res="<<*this<<std::endl;

    return *this;
}
*/

// This code works in place, and seems to word even when the copying code
// succumbs to compiler optimisation errors
TaylorModel&
TaylorModel::truncate(uint d)
{
    //std::cerr<<"truncate(tm="<<*this<<", d="<<d<<")"<<std::endl;
    TaylorModel::const_iterator end=this->end();
    TaylorModel::const_iterator adv=this->begin();
    TaylorModel::iterator curr=this->begin();
    set_rounding_upward();
    while(adv!=end) {
        if(adv->key().degree()<=d) {
            *curr=*adv; ++curr;
        } else {
            this->error()+=adv->data();
        }
        ++adv;
    }
    this->expansion().resize(curr-this->begin());
    return *this;
}

TaylorModel&
TaylorModel::truncate(const MultiIndex& b)
{
    //std::cerr<<"truncate(tm="<<*this<<", d="<<d<<")"<<std::endl;
    TaylorModel::const_iterator end=this->end();
    TaylorModel::const_iterator adv=this->begin();
    TaylorModel::iterator curr=this->begin();
    set_rounding_upward();
    while(adv!=end) {
        if(adv->key()<=b) {
            *curr=*adv; ++curr;
        } else {
            this->error()+=adv->data();
        }
        ++adv;
    }
    this->expansion().resize(curr-this->begin());
    return *this;
}



TaylorModel&
TaylorModel::truncate(const MultiIndexBound& b)
{
    ARIADNE_NOT_IMPLEMENTED;
    ARIADNE_ASSERT(b.size()==this->argument_size());
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(!(iter->key()<=b)) {
            e+=abs(iter->data());
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
        if(iter->key().degree()>o) {
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
        const MultiIndex& a=iter->key();
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

// Accuracy control

void TaylorModel::set_maximum_degree(uint d) {
    this->_accuracy_ptr->_maximum_degree=d;
}

void TaylorModel::set_sweep_threshold(double e) {
    this->_accuracy_ptr->_sweep_threshold=e;
}

uint TaylorModel::maximum_degree() const {
    return this->_accuracy_ptr->_maximum_degree;
}

double TaylorModel::sweep_threshold() const {
    return this->_accuracy_ptr->_sweep_threshold;
}



///////////////////////////////////////////////////////////////////////////////

// Arithmetic operators

TaylorModel&
operator+=(TaylorModel& x, const TaylorModel& y)
{
    if(x.argument_size()==0) { x=TaylorModel::zero(y.argument_size()); }
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    _acc(x,y); return x;

}

TaylorModel&
operator-=(TaylorModel& x, const TaylorModel& y)
{
    if(x.argument_size()==0) { x=TaylorModel::zero(y.argument_size()); }
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    _acc(x,neg(y)); return x;
}


TaylorModel&
operator+=(TaylorModel& x, const Float& c)
{
    _acc(x,c); return x;
}

TaylorModel&
operator-=(TaylorModel& x, const Float& c)
{
    _acc(x,-c); return x;
}


TaylorModel&
operator+=(TaylorModel& x, const Interval& c)
{
    _acc(x,c); return x;
}

TaylorModel&
operator-=(TaylorModel& x, const Interval& c)
{
    _acc(x,-c); return x;
}

TaylorModel&
operator*=(TaylorModel& x, const Float& c)
{
    _scal(x,c); return x;
}

TaylorModel&
operator*=(TaylorModel& x, const Interval& c)
{
    _scal(x,c); return x;
}


TaylorModel&
operator/=(TaylorModel& x, const Float& c)
{
    if(c==0) {
        ARIADNE_THROW(DivideByZeroException,"operator/=(TaylorModel x,Float c)","x="<<x<<" c="<<c);
    }
    _scal(x,Interval(1.0)/c); return x;
}


TaylorModel&
operator/=(TaylorModel& x, const Interval& c)
{
    if(c.upper()>=0 && c.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"operator/=(TaylorModel x,Interval c)","x="<<x<<" c="<<c);
    }
    _scal(x,1.0/c); return x;
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
    TaylorModel r(x.argument_size()); _add(r,x,y); return r;
}

TaylorModel
operator-(const TaylorModel& x, const TaylorModel& y) {
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x=("<<x.argument_size()<<")"<<x<<"y="<<y<<"\n");
    TaylorModel r=neg(y); _acc(r,x); return r;
}

TaylorModel
operator*(const TaylorModel& x, const TaylorModel& y) {
    if(x.argument_size()!=y.argument_size()) { std::cerr<<"operator*(TaylorModel x, TaylorModel y)\n  x="<<x<<" y="<<y<<"\n"; }
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    TaylorModel r(y.argument_size(),y.accuracy_ptr()); _mul(r,x,y); return r;
}

TaylorModel
operator/(const TaylorModel& x, const TaylorModel& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    TaylorModel r(x.argument_size()); _mul(r,x,rec(y)); return r;
}



TaylorModel
operator+(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); _acc(r,c); return r;
}

TaylorModel
operator-(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); _acc(r,-c); return r;
}

TaylorModel
operator*(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); _scal(r,c); return r;
}

TaylorModel
operator/(const TaylorModel& x, const Float& c) {
    TaylorModel r(x); _scal(r,Interval(1)/c); return r;
}

TaylorModel
operator+(const Float& c, const TaylorModel& x) {
    TaylorModel r(x); _acc(r,c); return r;
}

TaylorModel
operator-(const Float& c, const TaylorModel& x) {
    TaylorModel r=neg(x); _acc(r,c); return r;
}

TaylorModel
operator*(const Float& c, const TaylorModel& x) {
    TaylorModel r(x); _scal(r,c); return r;
}

TaylorModel
operator/(const Float& c, const TaylorModel& x) {
    TaylorModel r(rec(x)); _scal(r,Interval(1)/c); return r;
}



TaylorModel
operator+(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); _acc(r,c); return r;
}

TaylorModel
operator-(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); _acc(r,-c); return r;
}

TaylorModel
operator*(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); _scal(r,c); return r;
}

TaylorModel
operator/(const TaylorModel& x, const Interval& c) {
    TaylorModel r(x); _scal(r,1/c); return r;
}

TaylorModel
operator+(const Interval& c, const TaylorModel& x) {
    TaylorModel r(x); _acc(r,c); return r;
}

TaylorModel
operator-(const Interval& c, const TaylorModel& x) {
    TaylorModel r=neg(x); _acc(r,c); return r;
}

TaylorModel
operator*(const Interval& c, const TaylorModel& x) {
    TaylorModel r(x); _scal(r,c); return r;
}

TaylorModel
operator/(const Interval& c, const TaylorModel& x) {
    TaylorModel r=rec(x); _scal(r,c); return r;
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
        r.expansion().append(xiter->key(),-xiter->data());
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
        if(iter->key().degree()==0) {
            r+=iter->data();
        } else {
            r+=abs(iter->data())*Interval(-1,1);
        }
    }
    return r;

    /* FIXME: The following code does not work with optimisation turned on using gcc */
    set_rounding_mode(upward);
    volatile Float t=this->error();
    volatile Float v=0.0;
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()==0) {
            v=iter->data();
        } else {
            t+=abs(iter->data());
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
    // Don't check for subset of domain,since this is an internal function, and
    // subset may fail due to intermediate roundoff errors
    //ARIADNE_ASSERT(subset(v,this->domain()));
    Interval r=this->error()*Interval(-1,+1);
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Interval t(iter->data());
        for(uint j=0; j!=iter->key().size(); ++j) {
            t*=pow(v[j],iter->key()[j]);
        }
        r+=t;
    }
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
class TaylorSeries;


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
compose(const TaylorSeries& ts, const TaylorModel& tm)
{
    return _compose(ts,tm,ec);
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
TaylorModel
_compose1(const series_function_pointer& fn, const TaylorModel& tm, Float eps)
{
    static const uint DEGREE=18;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tm.value();
    Interval r=tm.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    Float truncation_error_estimate=mag(range_series[d])*pow(mag(r-c),d);
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error_estimate
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorModel x=tm-c;
    TaylorModel res(tm.argument_size(),tm.accuracy_ptr());
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
_compose2(const series_function_pointer& fn, const TaylorModel& tm, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tm.value();
    Interval r=tm.range();
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

    TaylorModel x=tm-c;
    TaylorModel res(tm.argument_size(),tm.accuracy_ptr());
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
_compose3(const series_function_pointer& fn, const TaylorModel& tm, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tm.value();
    Interval r=tm.range();
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

    TaylorModel x=tm;
    TaylorModel res(tm.argument_size(),tm.accuracy_ptr());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        //res.sweep(eps);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


TaylorModel
_compose(const series_function_pointer& fn, const TaylorModel& tm, Float eps) {
    //std::cerr<<"_compose(SeriesFunction,TaylorModel,Error)\n";
    return _compose3(fn,tm,eps);
}




///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   bounded domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)

namespace {
inline int pow2(uint k) { return 1<<k; }
inline int powm1(uint k) { return (k%2) ? -1 : +1; }
double rec_fac_up(uint n) { set_rounding_upward(); double r=1.0; for(uint i=1; i<=n; ++i) { r/=i; } return r; }
}

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
    uint d=uint(log((1-eps)*x.sweep_threshold())/log(eps)+1);
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    TaylorModel y=(x/a)-1.0;
    //std::cerr<<"y="<<y<<std::endl;
    TaylorModel z(x.argument_size(),x.accuracy_ptr());
    Series<Interval> sqrt_series=Series<Interval>::sqrt(d,Interval(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(uint i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        z.clean();
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

    uint d=uint(log((1-eps)*x.sweep_threshold())/log(eps))+1;

    TaylorModel y=1-(x/a);
    TaylorModel z(x.argument_size(),x.accuracy_ptr());
    z+=Float(d%2?-1:+1);
    //std::cerr<<"  y="<<y<<"\n";
    //std::cerr<<"  z="<<z<<"\n";
    for(uint i=0; i!=d; ++i) {
        z=1.0 + z * y;
        //std::cerr<<"  z="<<z<<"\n";
    }

    // Compute the truncation error
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
    uint d=uint(log((1-eps)*x.sweep_threshold())/log(eps)+1);
    TaylorModel y=x/a-1;
    TaylorModel z(x.argument_size(),x.accuracy_ptr());
    z+=Float(d%2?-1:+1)/d;
    for(uint i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
        z.clean();
    }
    z=z*y;
    z.sweep();
    Float trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(Interval(a))+trunc_err*Interval(-1,1);
}

// Use special code to utilise exp(ax+b)=exp(x)^a*exp(b)
TaylorModel exp(const TaylorModel& x) {
    // FIXME: Truncation error may be incorrect


    // Scale to unit interval
    TaylorModel y=x;
    Float xval=x.value();
    y-=xval;
    Float xrad=mag(y.range());
    uint sfp=0; // A number such that 2^sfp>rad(x.range())
    while(Float(1<<sfp)<xrad) { ++sfp; }
    double sf=1.0/(1<<sfp);
    _scal_exact(y,sf);
    Float yrad=xrad*sf;

    uint degree=x.maximum_degree();

    // Since x is in unit domain, truncation error is no worse than maximum omitted term, i.e. xr/fac(d+1)
    TaylorModel res(x.argument_size(),x.accuracy_ptr());
    res.set_error(pow_up(yrad,degree+1));
    for(uint i=0; i!=degree; ++i) {
        res/=(degree-i);
        res=y*res+1.0;
    }

    // Square r a total of sfp times
    TaylorModel square(x.argument_size(),x.accuracy_ptr());
    for(uint i=0; i!=sfp; ++i) {
        _mul(square,res,res);
        res.swap(square);
        square.clear();
    }

    // Multiply by exp(xv)
    res*=Ariadne::exp(Interval(xval));

    return res;
    //return _compose(&Series<Interval>::exp,x,x.sweep_threshold());
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form x*f(x^2)
TaylorModel sin(const TaylorModel& x) {
    // FIXME: Truncation error may be incorrect
    TaylorModel y(x.argument_size(),x.accuracy_ptr());
    TaylorModel s(x.argument_size(),x.accuracy_ptr());
    TaylorModel r(x.argument_size(),x.accuracy_ptr());
    TaylorModel t(x.argument_size(),x.accuracy_ptr());

    double two_pi=2*pi<double>();
    int n=floor(x.value()/two_pi + 0.5);
    y=x-(n*2*Ariadne::pi<Interval>());

    if(y.error()>two_pi/2 || mag(y.range())*rec_fac_up(x.maximum_degree())>1) {
        r.error()=1.0;
    } else {
        _mul(s,y,y); 

        int d=(x.maximum_degree()+3)/2;
        Float srad=mag(s.range());
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute x(1-y/6+y^2/120-y^3/5040+... = x(1-y/6*(1-y/20*(1-y/42*...)
        r=1.0;
        for(int i=0; i!=d; ++i) {
            r/=double(-2*(d-i)*(2*(d-i)+1));
            _mul(t,s,r); r.swap(t); t.clear();
            r+=1.0;
        }
        _mul(t,y,r); r.swap(t);

        r.error()+=truncation_error;
    }

    return r;
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form f(x^2)
TaylorModel cos(const TaylorModel& x) {
    // FIXME: Truncation error may be incorrect
    TaylorModel y(x.argument_size(),x.accuracy_ptr());
    TaylorModel s(x.argument_size(),x.accuracy_ptr());
    TaylorModel r(x.argument_size(),x.accuracy_ptr());
    TaylorModel t(x.argument_size(),x.accuracy_ptr());

    double two_pi=2*pi<double>();
    int n=floor(x.value()/two_pi + 0.5);

    y=x-2*n*pi<Interval>();

    if(y.error()>two_pi/2 || mag(y.range())*rec_fac_up(x.maximum_degree())>1) {
        r.error()=1.0;
    } else {
        _mul(s,y,y);

        int d=(x.maximum_degree()+3)/2;
        Float srad=mag(s.range());
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute 1-y/2+y^2/24-y^3/720+... = (1-y/2*(1-y/12*(1-y/30*...)
        r=1.0;
        for(int i=0; i!=d; ++i) {
            r/=double(-2*(d-i)*(2*(d-i)-1));
            _mul(t,s,r); r.swap(t); t.clear();
            r+=1.0;
        }

        r.error()+=truncation_error;
    }

    return r;
}

TaylorModel tan(const TaylorModel& x) {
    return sin(x)*rec(cos(x));
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

    ARIADNE_ASSERT_MSG(ocd.radius()>0,"Illegal scaling from interval "<<ocd<<" with zero radius to interval "<<ncd);
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
            Float& c=iter->data();
            Interval ci=c;
            for(uint j=0; j!=as; ++j) {
                ci*=sf[j][iter->key()[j]];
            }
            iter->data()=ci.midpoint();
            x._error=add_up(x._error,mag(ci-iter->data()));
        }
    }

    return x;
}

inline double div_rnd(double x, unsigned int m) { return (volatile double&)x/(volatile unsigned int&)m; }
inline double div_rnd(double x, int n) { return (volatile double&)x/(volatile int&)n; }

// Compute antiderivative by first computing the error and then the individual terms
// Seems to have problems with setting the rounding mode; with certain compiler flags
// sometimes the truncation error is not set
void _antidifferentiate1(TaylorModel& x, uint k)
{
    ARIADNE_ASSERT(k<x.argument_size());

    //std::cerr<<"xe="<<xe<<"\n";
    set_rounding_mode(upward);
    volatile Float tre=0; // Twice the maximum accumulated roundoff error
    for(TaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const uint c=xiter->key()[k]+1;
        register double xv=xiter->data();
        register double mxv=-xv;
        assert(c>0);
        //volatile double ml,u;
        if(xv>=0) {
            //volatile double u=xv/c;
            //volatile double ml=mxv/c;
            tre+=add_rnd(div_rnd(xv,c),div_rnd(mxv,c));
        } else {
            tre+=add_rnd(div_rnd(mxv,c),div_rnd(xv,c));
            //volatile double u=mxv/c;
            //volatile double ml=xv/c;
            //te+=(u+ml);
        }
        //std::cerr<<"  te="<<te;;
    }
    //std::cerr<<"  te="<<tre;;
    x.error()+=(tre/2);
    //xe+=te/2;
    //std::cerr<<"xe="<<xe<<"\n";

    set_rounding_to_nearest();
    for(TaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex& a=xiter->key();
        Float& v=xiter->data();

        // Use the following code as opposed to ++a[k]; c=a[k];
        // as this seems to be safer against compiler optimisations through MultiIndexElementReference
        a.increment(k);
        const int& c=a.at(k);
        assert(c>0);
        v/=c;
    }
}

// Compute antiderivative by computing term-by-term, switching the rounding mode
// May be slow if no assembly rounding mode switching is available.
void _antidifferentiate2(TaylorModel& x, uint k)
{
    ARIADNE_ASSERT(k<x.argument_size());

    //std::cerr<<"xe="<<xe<<"\n";
    set_rounding_mode(upward);
    volatile Float tre=0; // Twice the maximum accumulated roundoff error
    volatile double u,ml;
    uint c;
    for(TaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex& xa=xiter->key();
        Float& xv=xiter->data();
        Float mxv=-xv;
        c=xa.at(k)+1;

        set_rounding_upward();
        if(xv>=0) {
            u=xv/c;
            ml=mxv/c;
        } else {
            u=mxv/c;
            ml=xv/c;
        }
        tre+=(u+ml);

        set_rounding_to_nearest();
        xa.increment(k);
        xv/=c;
    }

    set_rounding_upward();
    x.error()+=(tre/2);
    set_rounding_to_nearest();
}


TaylorModel& TaylorModel::antidifferentiate(uint k)
{
    _antidifferentiate2(*this,k); return *this;
}

// Compute derivative by computing term-by-term, switching the rounding mode
// May be slow if no assembly rounding mode switching is available.


TaylorModel derivative(const TaylorModel& x, uint k)
{
    ARIADNE_ASSERT(k<x.argument_size());
    TaylorModel r(x.argument_size(),x.accuracy_ptr());

    //std::cerr<<"xe="<<xe<<"\n";
    volatile Float tre=0; // Twice the maximum accumulated roundoff error
    double xv;
    volatile double u,ml,mxv;
    MultiIndex a(x.argument_size());
    uint c;
    for(TaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        c=xiter->key().at(k);
        if(c!=0) {
            a=xiter->key();
            xv=xiter->data();
            mxv=-xv;

            if(c!=1 && c!=2 && c!=4) {
                set_rounding_upward();
                u=xv*c;
                ml=mxv*c;
                tre+=(u+ml);
                set_rounding_to_nearest();
            }
            --a[k];
            xv*=c;
            r.expansion().append(a,xv);
        }
    }
    r.error()=(tre/2);

    return r;
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
    r1.set_accuracy(tv.accuracy_ptr());
    s[j]=TaylorModel::scaling(as,j,Interval(0,+1));
    TaylorModel r2=compose(tv,s);
    r2.set_accuracy(tv.accuracy_ptr());
    ARIADNE_ASSERT(r1.accuracy_ptr()==tv.accuracy_ptr());
    ARIADNE_ASSERT(r2.accuracy_ptr()==tv.accuracy_ptr());
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
        r+=abs(iter->data());
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
    return TaylorModel(embed(0u,x.expansion(),as),x.error(),x.accuracy_ptr());
}

TaylorModel embed(uint as, const TaylorModel& x)
{
    return TaylorModel(embed(as,x.expansion(),0u),x.error(),x.accuracy_ptr());
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
    Vector<TaylorModel> r(x);
    for(uint i=0; i!=x.size(); ++i) {
        r[i].antidifferentiate(k);
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

TaylorModel intersection(const TaylorModel& x, const TaylorModel& y) {
    TaylorModel r(x.argument_size(),x.accuracy_ptr());
    Float twice_max_error=0.0;

    const Float& xe=x.error();
    const Float& ye=y.error();
    volatile Float rv,xv,yv,xu,yu,mxl,myl,u,ml;
    //const MultiIndex* aptr;
    MultiIndex a;

    TaylorModel::const_iterator xiter=x.begin();
    TaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() || yiter!=y.end()) {
        // Can't use const MultiIndex& here since references change as the iterators change
        // We would need to use a smart reference
        //const MultiIndex& xa=xiter->key();
        //const MultiIndex& ya=yiter->key();
        if(xiter==x.end()) {
            a=yiter->key();
            yv=yiter->data();
            xv=0.0;
            ++yiter;
        } else if(yiter==y.end()) {
            a=xiter->key();
            xv=xiter->data();
            yv=0.0;
            ++xiter;
        } else if(xiter->key()==yiter->key()) {
            a=xiter->key();
            xv=xiter->data();
            yv=yiter->data();
            ++xiter;
            ++yiter;
        } else if(xiter->key()<yiter->key()) {
            a=xiter->key();
            xv=xiter->data();
            yv=0.0;
            ++xiter;
        } else { // xa>ya 
            a=yiter->key();
            yv=yiter->data();
            xv=0.0;
            ++yiter;
        }

        set_rounding_upward();
        xu=xv+xe;
        yu=yv+ye;
        mxl=xe-xv;
        myl=ye-yv;
        u=min(xu,yu);
        ml=min(mxl,myl);
        if(u+ml<0) {
            ARIADNE_THROW(IntersectionException,"intersection(TaylorModel,TaylorModel)",x<<" and "<<y<<" are disjoint.");
        }
        twice_max_error=max(twice_max_error,(u+ml));

        set_rounding_to_nearest();
        xu=xv+xe;
        yu=yv+ye;
        mxl=xe-xv;
        myl=ye-yv;
        u=min(xu,yu);
        ml=min(mxl,myl);
        rv=(u-ml)/2;
        if(rv!=0.0) { r.expansion().append(a,(u-ml)/2); }
    }

    r.error()=twice_max_error/2;
    
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
            const MultiIndex& a=iter->key();
            const double& x=iter->data();
            for(uint k=0; k!=as; ++k) {
                const uint c=a[k];
                if(c>0) {
                    if(iter->key().degree()==1) {
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
                const uint c=iter->key()[k];
                if(c>0) {
                    const double& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
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
    Matrix<Interval> J(rs,rs);
    for(uint i=0; i!=rs; ++i) {
        for(TaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=rs; ++k) {
                const uint c=iter->key()[has+k];
                if(c>0) {
                    const double& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}

///////////////////////////////////////////////////////////////////////////////

// Vector operators (compose, solve, implicit, flow




// Compose by computing each and every term individually without caching
// Easy to implement, but far too slow
inline
Vector<TaylorModel>
_compose1(const Vector<TaylorModel>& x,
          const Vector<TaylorModel>& ys)
{
    //std::cerr<<"compose1"<<std::endl;
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
            t=iter->data();
            for(uint j=0; j!=iter->key().size(); ++j) {
                TaylorModel p=pow(ys[j],iter->key()[j]);
                t=t*p;
            }
            r[i]+=t;
        }
    }

    return r;
}


inline
Vector<TaylorModel>
_compose2(const Vector<TaylorModel>& x,
          const Vector<TaylorModel>& ys)
{
    //std::cerr<<"compose2"<<std::endl;
    uint yrs=ys.size();
    uint xas=ys.size();
    uint as=ys[0].argument_size();
    shared_ptr<TaylorModel::Accuracy> accuracy_ptr=ys[0].accuracy_ptr();
    
    array<uchar> max_power(ys.size());
    for(uint j=0; j!=ys.size(); ++j) { max_power[j]=1; }

    for(uint i=0; i!=x.size(); ++i) {
        for(TaylorModel::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            assert(xas==iter->key().size());
            for(uint j=0; j!=iter->key().size(); ++j) {
                max_power[j]=max(max_power[j],iter->key()[j]);
            }
        }
    }

    array< array< TaylorModel > > powers(yrs);
    for(uint j=0; j!=yrs; ++j) {
        powers[j].resize(max_power[j]+1);
        powers[j][0]=ys[j]*0;
        powers[j][1]=ys[j];
        for(uint k=2; k!=powers[j].size(); ++k) {
            powers[j][k]=powers[j][k/2]*powers[j][(k+1)/2];
        }
    }

    Vector<TaylorModel> r(x.size(),TaylorModel(as,accuracy_ptr));
    TaylorModel t(as,accuracy_ptr);
    MultiIndex a;
    Float c;
    for(uint i=0; i!=x.size(); ++i) {
        r[i].set_error(x[i].error());
        for(TaylorModel::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            a=iter->key();
            c=iter->data();
            t=c;
            for(uint j=0; j!=a.size(); ++j) {
                if(a[j]>0) {
                    t=t*powers[j][a[j]];
                }
            }
            r[i]+=t;
        }
    }

    ARIADNE_ASSERT(r[0].accuracy_ptr()==ys[0].accuracy_ptr());

    return r;
}

Vector<TaylorModel>
_compose(const Vector<TaylorModel>& x,
         const Vector<TaylorModel>& ys)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(ys.size()==x[0].argument_size());
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }
    for(uint i=1; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys[0].argument_size(),"ys="<<ys); }

    return _compose2(x,ys);
}


Vector<TaylorModel>
unchecked_compose(const Vector<TaylorModel>& x,
                  const Vector<TaylorModel>& ys)
{
    return _compose(x,ys);
}

TaylorModel
unchecked_compose(const TaylorModel& x,
                  const Vector<TaylorModel>& ys)
{
    return _compose(Vector<TaylorModel>(1u,x),ys)[0];
}

Vector<TaylorModel>
compose(const Vector<TaylorModel>& x,
        const Vector<TaylorModel>& ys)
{
    return _compose(x,ys);
}

TaylorModel
compose(const TaylorModel& x,
        const Vector<TaylorModel>& ys)
{
    return _compose(Vector<TaylorModel>(1u,x),ys)[0];
}

TaylorModel
compose(const TaylorModel& x,
        const Vector<Interval>& d,
        const Vector<TaylorModel>& y)
{
    return _compose(Vector<TaylorModel>(1u,x),unscale(y,d))[0];
}

Vector<TaylorModel>
compose(const Vector<TaylorModel>& x,
        const Vector<Interval>& d,
        const Vector<TaylorModel>& y)
{
    return _compose(x,unscale(y,d));
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
    uint number_non_contracting=rs;
    array<bool> contracting(rs,false);
    for(uint k=0; k!=n; ++k) {
        h_new=compose(g,join(id,h));
        for(uint i=0; i!=rs; ++i) {
            if(!contracting[i]) {
                if(refines(h_new,h)) {
                    contracting[i]=true;
                    --number_non_contracting;
                }
                if(disjoint(h_new,h)) {
                    ARIADNE_THROW(ImplicitFunctionException,"implicit(Vector<TaylorModel> f)",
                                "(with f="<<f<<"): Application of Newton solver to "<<h<<" yields "<<h_new<<
                                " which is disjoint. No solution");
                }
            }
        }
        for(uint i=0; i!=rs; ++i) {
            h[i]=intersection(h[i],h_new[i]);
        }
    }

    if(number_non_contracting==0) {
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
_flow_step(const Vector<TaylorModel>& vf, const Vector<TaylorModel>& yz,
           const Vector<TaylorModel>& phi)
{
    const uint n=vf.size();
    Vector<TaylorModel> new_phi=compose(vf,phi);
    for(uint i=0; i!=n; ++i) {
        new_phi[i].antidifferentiate(n);
        new_phi[i]+=yz[i];
    }
    return new_phi;
}

Vector<TaylorModel>
unchecked_flow(const Vector<TaylorModel>& vf, const Vector<TaylorModel>& y0, uint order)
{
    uint n=vf.size();
    assert(y0.size()==n);
    assert(y0[0].argument_size()==n || y0[0].argument_size()==n+1);

    // Set inital set at time zero; embed in space including time variable if necessary
    Vector<TaylorModel> yz;
    if(y0[0].argument_size()==n) { yz=embed(y0,1); } else { yz=y0; }

    // Set initial conditions
    // The Perron-Frobenius operator should act as a contraction on this set
    Vector<TaylorModel> y(n);
    for(uint i=0; i!=y.size(); ++i) {
        y[i]=TaylorModel::constant(n+1,yz[i].range()+vf[i].range()*Interval(-1,+1));
    }

    Vector<TaylorModel> new_y(n,TaylorModel(n+1));


    // Initialise to handle initial conditions and initial time step
/*
    new_y=compose(vf,y);
    for(uint i=0; i!=n; ++i) {
        new_y[i].antidifferentiate(n);
        new_y[i]*=h_rad;
        new_y[i]+=yz[i];
    }
    for(uint i=0; i!=n; ++i) {
        new_y[i].swap(y[i]);
    }
*/

    //std::cerr << "\ny[0]=" << y << std::endl << std::endl;
    //std::cerr<<"  y="<<y<<"\n";
    for(uint j=0; j!=order; ++j) {

        new_y=compose(vf,y);
        for(uint i=0; i!=n; ++i) {
            new_y[i].antidifferentiate(n);
            new_y[i]+=yz[i];

            //new_y[i]=intersection(y[i],new_y[i]);

            // The following code ignores disjoint set errors
            // and simply uses the next iterate as the new approximation
            // I'm not sure if this is entirely safe
            try {
                new_y[i]=intersection(y[i],new_y[i]);
            }
            catch(const IntersectionException& e) {
                std::cerr<<"Warning: "<<e.what()<<"\n";
            }
        }

        for(uint i=0; i!=n; ++i) {
            new_y[i].swap(y[i]);
        }
    }

    //for(uint i=0; i!=n; ++i) { y[i].clobber(so,to); }
    for(uint i=0; i!=n; ++i) { y[i].sweep(0.0); }

    return y;
}



Vector<TaylorModel>
flow(const Vector<TaylorModel>& vf, const Vector<TaylorModel>& yz, uint o)
{
    return unchecked_flow(vf,yz,o);
}



} //namespace Ariadne


