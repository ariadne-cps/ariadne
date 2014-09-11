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
#include <limits>

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
#include "function.h"
#include "exceptions.h"

#include "algebra_mixin.tcc"

#define VOLATILE ;
#include <include/multi_index-noaliasing.h>
#include <include/function_mixin.h>
#include <include/vector.h>

namespace Ariadne {

namespace {

bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }

} // namespace


const double em=2.2204460492503131e-16;
const double ec=em/2;

/*

typedef Float ErrorFloat;
typedef Float ApproxFloat;

void axpy(Float& te, ApproxFloat& r, const Float& a, const Float& x, const Float& y) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=a*xv+y;
    VOLATILE Float ml=a*mxv-y;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=a*xv+y;
}

void add_op(Float& te, ApproxFloat& r, const Float& x, const Float& y) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=xv+y;
    VOLATILE Float ml=mxv-y;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv+y;
}

void sub_op(Float& te, ApproxFloat& r, const Float& x, const Float& y) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=xv-y;
    VOLATILE Float ml=mxv+y;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv-y;
}

void mul_op(Float& te, ApproxFloat& r, const Float& s, const Float& x) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u=xv*s;
    VOLATILE Float ml=mxv*s;
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv*s;
}

void mul_op(Float& te, ApproxFloat& r, const Float& sl, const Float& sm, const Float& su, const Float& x) {
    VOLATILE Float& xv=const_cast<VOLATILE Float&>(x);
    VOLATILE Float mxv=-x;
    set_rounding_upward();
    VOLATILE Float u,ml;
    if(x>=0) {
        u=xv*su;
        ml=mxv*sl;
    } else {
        u=xv*sl;
        ml=mxv*su;
    }
    te+=(u+ml);
    set_rounding_to_nearest();
    r=xv*sm;
}

*/



Vector<ValidatedNumberType> unscale(const Vector<ValidatedNumberType>& x, const Vector<Interval>& d) {
    Vector<ValidatedNumberType> r(x);
    for(uint i=0; i!=r.size(); ++i) {
        if(d[i].lower()==d[i].upper()) {
            if(x==d) {
                r[i]=ValidatedNumberType(0.0,0.0);
            } else {
                r[i]=ValidatedNumberType(-inf,+inf);
            }
        } else {
            r[i]=(2*r[i]-add_ivl(d[i].lower(),d[i].upper()))/sub_ivl(d[i].upper(),d[i].lower());
        }
    }
    return r;
}



/* */

ValidatedTaylorModel::TaylorModel()
    : _expansion(0), _error(0), _sweeper()
{ }


ValidatedTaylorModel::TaylorModel(uint as, Sweeper swp)
    : _expansion(as), _error(0), _sweeper(swp)
{
}

ValidatedTaylorModel::TaylorModel(const Expansion<Float>& f, const Float& e, Sweeper swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->unique_sort();
    this->sweep();
}


ValidatedTaylorModel
ValidatedTaylorModel::create() const
{
    return ValidatedTaylorModel(this->argument_size(),this->_sweeper);
}

ValidatedTaylorModel
ValidatedTaylorModel::create_ball(Float e) const
{
    ARIADNE_PRECONDITION(e>=0);
    ValidatedTaylorModel r(this->argument_size(),this->_sweeper);
    r.set_error(e);
    return r;
}

void
ValidatedTaylorModel::swap(ValidatedTaylorModel& tm)
{
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
    std::swap(this->_sweeper,tm._sweeper);
}

void
ValidatedTaylorModel::clear()
{
    this->_expansion.clear();
    this->_error=0.0;
}

uint
ValidatedTaylorModel::degree() const
{
    uchar deg=0u;
    for(ValidatedTaylorModel::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        deg=std::max(deg,iter->key().degree());
    }
    return deg;
}

ValidatedTaylorModel&
ValidatedTaylorModel::operator=(double c)
{
    return this->operator=(Float(c));
}

ValidatedTaylorModel&
ValidatedTaylorModel::operator=(const Float& c)
{
    this->_expansion.clear();
    if(c!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),c);
    }
    this->_error=0;
    return *this;
}

ValidatedTaylorModel&
ValidatedTaylorModel::operator=(const Real& c)
{
    return this->operator=(ValidatedNumberType(c));
}

ValidatedTaylorModel&
ValidatedTaylorModel::operator=(const ValidatedNumberType& c)
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
void _neg(ValidatedTaylorModel& r)
{
    for(ValidatedTaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data()=-iter->data();
    }
}

inline void _scal_exact(ValidatedTaylorModel& r, const Float& c)
{
    // Operation can be performed exactly
    register Float pc=c;
    for(ValidatedTaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data()*=pc;
    }
    r.error()*=abs(pc);
    return;
}




inline void _scal_approx(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error();
    set_rounding_upward();
    Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    Float pc=c;
    Float mc=-c;
    for(ValidatedTaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}



inline void _scal_approx4(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error();
    set_rounding_upward();
    Float u,ml;
    //Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    Float pc=c;
    Float mc=-pc;
    for(ValidatedTaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data();
        //volatile Float& rv=const_cast<volatile Float&>(riter->data());
        u=rv*pc;
        ml=rv*mc;
        te+=u+ml;
    }
    re*=abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}

inline void _scal_approx3(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    double& re=internal_cast<double&>(r.error());
    set_rounding_upward();
    volatile double u,ml;
    double te=0; // Twice the maximum accumulated error
    register double pc=internal_cast<const double&>(c);
    register double mc=-pc;
    for(ValidatedTaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const double& rv=internal_cast<const double&>(riter->data());
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=std::abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}


// Version using only Float class directly
inline void _scal_approx2(ValidatedTaylorModel& r, const Float& c)
{
    // General case with error analysis
    Float& re=r.error();
    set_rounding_upward();
    Float u,ml;
    register Float te=0; // Twice the maximum accumulated error
    register Float pc=c;
    register Float mc=-c;
    for(ValidatedTaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->data();
        u=mul_rnd(rv,pc);
        ml=mul_rnd(rv,mc);
        te+=(u+ml);
    }
    re*=abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=c;
    for(ValidatedTaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}

// Modified version to try to obtain optimal performance
inline void _scal_approx1(ValidatedTaylorModel& rr, const Float& cc)
{
    // General case with error analysis
    Expansion<double>& r=reinterpret_cast<Expansion<double>&>(rr.expansion());
    double& re=reinterpret_cast<double&>(rr.error());
    const double& c=reinterpret_cast<const double&>(cc);

    set_rounding_upward();
    volatile double u,ml;
    double te=0; // Twice the maximum accumulated error
    register double pc=c;
    register double mc=-c;
    for(Expansion<double>::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const double& rv=riter->data();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=std::abs(pc);
    re+=te/2;

    set_rounding_to_nearest();
    double m=c;
    for(Expansion<double>::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}


// Original version except with double
inline void _scal_approx0(ValidatedTaylorModel& rr, const Float& cc)
{
    Expansion<double>& r=reinterpret_cast<Expansion<double>&>(rr.expansion());
    volatile double& re=static_cast<volatile double&>(rr.error());
    volatile double c=static_cast<const double&>(cc);
    //const double& c=cc.v;

    set_rounding_upward();
    volatile double u,ml;
    double te=0; // Twice the maximum accumulated error
    register double pc=c;
    register double mc=-c;
    for(Expansion<double>::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const double& rv=riter->data();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=std::abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    double m=c;
    for(Expansion<double>::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }
    return;
}

void _scal(ValidatedTaylorModel& r, const Float& c)
{
    // No measurable speedup in general case by avoiding checks
    _scal_approx4(r,c); return;
    if(c==0.0) { r.expansion().clear(); r.error()=0; }
    else if(c==1.0) { }
    else if(c==2.0 || c==0.5) { _scal_exact(r,c); }
    else { _scal_approx(r,c); }
}


void _scal(ValidatedTaylorModel& r, const ValidatedNumberType& c)
{
    ARIADNE_ASSERT_MSG(c.lower()<=c.upper(),c);
    ARIADNE_ASSERT_MSG(r.error()>=0,"r="<<r);

    if(r.error()==inf) { r.expansion().clear(); return; }

    //std::cerr<<"ValidatedTaylorModel::scal(ValidatedNumberType c) c="<<c<<std::endl;
    Float& re=r.error();
    set_rounding_upward();
    VOLATILE Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    register Float cu=c.upper();
    register Float mcl=-c.lower();
    for(ValidatedTaylorModel::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
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
    Float m=(c.upper()+c.lower())/2;
    for(ValidatedTaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=m;
    }

    ARIADNE_ASSERT(r.error()>=0);
    return;
}

void _scal2(ValidatedTaylorModel& r, const ValidatedNumberType& c)
{
    //std::cerr<<"ValidatedTaylorModel::scal(ValidatedNumberType c) c="<<c<<std::endl;
    Float& re=r.error();
    VOLATILE Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    VOLATILE Float cu=c.upper();
    VOLATILE Float mcl=-c.lower();
    set_rounding_to_nearest();
    Float cm=(c.lower()+c.upper())/2;
    for(ValidatedTaylorModel::iterator riter=r.begin(); riter!=r.end(); ++riter) {
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

    ARIADNE_ASSERT(r.error()>=0);
    return;
}


struct UnitMultiIndex { uint argument_size; uint unit_index; };


inline void _incr(ValidatedTaylorModel& r, const MultiIndex& a)
{
    for(ValidatedTaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        static_cast<MultiIndex&>(iter->key())+=a;
    }
}

inline void _incr(ValidatedTaylorModel& r, uint j)
{
    for(ValidatedTaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++static_cast<MultiIndex&>(iter->key())[j];
    }
}


inline
void _acc(ValidatedTaylorModel& r, const Float& c)
{
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT(r.error()>=0);
    if(c==0) { return; }
    if(r.expansion().empty()) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else if((r.end()-1)->key().degree()>0) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else {
        Float& rv=(r.end()-1)->data();
        Float& re=r.error();
        set_rounding_upward();
        VOLATILE Float rvu=rv+c;
        VOLATILE Float mrvl=(-rv)-c;
        //std::cerr<<"re="<<re.u<<" ";
        re+=(rvu+mrvl)/2;
        //std::cerr<<"nre="<<re.u<<"\n";
        set_rounding_to_nearest();
        rv+=c;
    }
    ARIADNE_ASSERT_MSG(r.error()>=0,"c="<<c<<" r="<<r);
    return;
}



inline
void _acc(ValidatedTaylorModel& r, const ValidatedNumberType& c)
{
    // Compute self+=c
    ARIADNE_ASSERT_MSG(r.error()>=0,r);

    if(c.lower()==-inf || c.upper()==+inf) {
        r.clear();
        r.set_error(+inf);
        return;
    }

    if(c.lower()==-c.upper()) { // The midpoint of the interval is zero, so no need to change constant term
        set_rounding_upward();
        r.error()+=c.upper();
        set_rounding_to_nearest();
        return;
    }
    if(r.expansion().empty()) { // Append a constant term zero
        r.expansion().append(MultiIndex(r.argument_size()),0.0);
    } else if((r.end()-1)->key().degree()>0) { // Append a constant term zero
        r.expansion().append(MultiIndex(r.argument_size()),0.0);
    }

    Float& rv=(r.end()-1)->data();
    Float& re=r.error();
    set_rounding_upward();
    VOLATILE Float rvu=rv+c.upper();
    VOLATILE Float mrvl=(-rv)-c.lower();
    re+=(rvu+mrvl)/2;
    set_rounding_to_nearest();
    VOLATILE Float m=(c.upper()+c.lower())/2;
    rv+=m;
    ARIADNE_DEBUG_ASSERT(r.error()>=0);
}



// Compute r=x+y, assuming r is empty
// First compute the errors and then at the end compute the expansion
// This saves rounding mode changes, but requires an extra loop
inline void _add1(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    ValidatedTaylorModel::const_iterator xiter=x.begin();
    ValidatedTaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            ++yiter;
        } else {
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            VOLATILE Float u=xv+yv;
            VOLATILE Float t=-xv;
            VOLATILE Float ml=t-yv;
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
    ARIADNE_ASSERT(r.error()>=0);

}

// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
inline void _add2(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    set_rounding_upward();
    Float te=0.0;
    ValidatedTaylorModel::const_iterator xiter=x.begin();
    ValidatedTaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().append(yiter->key(),yiter->data());
            ++yiter;
        } else {
            assert(xiter->key()==yiter->key());
            VOLATILE Float& xv=const_cast<Float&>(static_cast<const Float&>(xiter->data()));
            VOLATILE Float& yv=const_cast<Float&>(static_cast<const Float&>(yiter->data()));
            set_rounding_upward();
            VOLATILE Float u=xv+yv;
            VOLATILE Float ml=-xv; ml-=yv;
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
    ARIADNE_ASSERT(r.error()>=0);

}

inline void _add(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    _add2(r,x,y);
    ARIADNE_ASSERT(r.error()>=0);
}

inline void _acc(ValidatedTaylorModel& r, const ValidatedTaylorModel& x)
{
    ValidatedTaylorModel s(r.argument_size(),r.sweeper()); _add(s,r,x); s.swap(r);
    ARIADNE_ASSERT(r.error()>=0);
}

inline void _sub(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    ValidatedTaylorModel::const_iterator xiter=x.begin();
    ValidatedTaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().append(yiter->key(),-yiter->data());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            VOLATILE Float u=xv-yv;
            VOLATILE Float t=-xv;
            VOLATILE Float ml=t+yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().append(xiter->key(),xiter->data()-yiter->data());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().append(yiter->key(),-yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error()>=0);
}

inline void _sma(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const Float& c, const ValidatedTaylorModel& y)
{
    if(c==1.0) { _add(r,x,y); return; }
    if(c==-1.0) { _sub(r,x,y); return; }

    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    VOLATILE Float u,ml;
    register Float mc=-c;

    ValidatedTaylorModel::const_iterator xiter=x.begin();
    ValidatedTaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            set_rounding_upward();
            const Float& yv=yiter->data();
            u=c*yv;
            ml=mc*yv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().append(yiter->key(),c*yiter->data());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            u=xv+c*yv;
            ml=mc*yv-xv;
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().append(xiter->key(),xiter->data()+c*yiter->data());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        set_rounding_upward();
        const Float& yv=yiter->data();
        VOLATILE Float u=c*yv;
        VOLATILE Float ml=mc*yv;
        te+=(u+ml);
        set_rounding_to_nearest();
        r.expansion().append(yiter->key(),c*yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error()>=0);
}

inline void _sma(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedNumberType& c, const ValidatedTaylorModel& y)
{
    if(c.lower()==c.upper()) {
        _sma(r,x,c.lower(),y); return;
    }

    ARIADNE_ASSERT_MSG(c.lower()<=c.upper(),c);
    ARIADNE_ASSERT_MSG(x.error()>=0,"x="<<x);
    ARIADNE_ASSERT_MSG(y.error()>=0,"y="<<y);

    //std::cerr<<"ValidatedTaylorModel::scal(ValidatedNumberType c) c="<<c<<std::endl;
    set_rounding_upward();
    VOLATILE Float u,ml,myv;
    Float te=0; // Twice the maximum accumulated error
    register Float cu=c.upper();
    register Float mcl=-c.lower();
    register Float cm=c.midpoint();

    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    ValidatedTaylorModel::const_iterator xiter=x.begin();
    ValidatedTaylorModel::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            set_rounding_upward();
            const Float& yv=yiter->data();
            if(yv>=0.0) {
                u=cu*yv;
                ml=mcl*yv;
            } else {
                myv = -yv;
                u = mcl*myv;
                ml = cu*myv;
            }
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().append(yiter->key(),cm*yiter->data());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            if(yv>0.0) {
                u=cu*yv+xv;
                ml=mcl*yv-xv;
            } else {
                myv=-yv;
                u=mcl*myv+xv;
                ml=cu*myv-xv;
            }
            te+=(u+ml);
            set_rounding_to_nearest();
            r.expansion().append(xiter->key(),xiter->data()+cm*yiter->data());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        set_rounding_upward();
        const Float& yv=yiter->data();
        if(yv>0.0) {
            u=cu*yv;
            ml=mcl*yv;
        } else {
            myv = -yv;
            u=mcl*myv;
            ml=cu*myv;
        }
        te+=(u+ml);
        set_rounding_to_nearest();
        r.expansion().append(yiter->key(),cm*yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_ASSERT(r.error()>=0);
}

struct Ivl { Float u; Float ml; };

/*
inline void _mul1(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    // Compute r+=x*y
    typedef ValidatedTaylorModel::const_iterator const_iterator;
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
            VOLATILE Float t=-xv;
            zv.ml+=t*yv;
        }
    }

    for(const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        Ivl& zv=z[riter->key()];
        const Float& rv=riter->data();
        zv.u+=rv; zv.ml-=rv;
    }

    VOLATILE Float te=0;
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
*/


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
void _mul2(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    ValidatedTaylorModel t(x.argument_size(),r.sweeper());
    ValidatedTaylorModel s(x.argument_size(),r.sweeper());
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        set_rounding_upward();
        Float te=0.0;
        for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->data();
            const Float& yv=yiter->data();
            VOLATILE Float u=xv*yv;
            VOLATILE Float ml=-xv; ml=ml*yv;
            te+=(u+ml);
        }
        t.error()=te/2;
        set_rounding_to_nearest();
        for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
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
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    r.sweep();

    set_rounding_to_nearest();
    return;
}


// Compute r+=x*y
// Compute monomial-by-monomial in y
// Change the rounding mode to avoid iterating using opposite rounding
inline void _mul3(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    ValidatedTaylorModel t(x.argument_size(),r.sweeper());
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        VOLATILE Float pxv=xiter->data();
        VOLATILE Float nxv=-pxv;
        VOLATILE Float te=0.0;
        for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
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
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    r.sweep();

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
void _mul4(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y, const Sweeper& accuracy)
{
    const uint as=r.argument_size();
    ValidatedTaylorModel t(as,r.sweeper());
    ValidatedTaylorModel s(as,r.sweeper());
    MultiIndex ta(as);
    Float tv;
    VOLATILE Float u;
    VOLATILE Float ml;
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        Float tte=0.0; // trucation error
        Float tre=0.0; // roundoff error
        const MultiIndex& xa=xiter->key();
        VOLATILE Float& xv=const_cast<Float&>(static_cast<const Float&>(xiter->data()));
        VOLATILE Float mxv=-xv;
        VOLATILE Float axv=abs(xv);
        for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const MultiIndex& ya=yiter->key();
            VOLATILE Float& yv=const_cast<Float&>(static_cast<const Float&>(yiter->data()));
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
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data());
    }

    Float ys=0;
    for(ValidatedTaylorModel::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data());
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    set_rounding_to_nearest();
    return;
}

inline void _mul(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    //_mul2(r,x,y);
    _mul4(r,x,y,r.sweeper());
    ARIADNE_ASSERT(r.error()>=0);
}



} // namespace

void _mul_full(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    _mul2(r,x,y);
}

void _mul_clear(ValidatedTaylorModel& r, const ValidatedTaylorModel& x, const ValidatedTaylorModel& y)
{
    _mul4(r,x,y,r.sweeper());
}


///////////////////////////////////////////////////////////////////////////////

// Inplace arithmetical operations for Algebra concept

void ValidatedTaylorModel::iadd(const ValidatedNumberType& c)
{
    _acc(*this,c);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error()>=0,*this);
}

void ValidatedTaylorModel::imul(const ValidatedNumberType& c)
{
    _scal(*this,c);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error()>=0,*this);
}

void ValidatedTaylorModel::isma(const ValidatedNumberType& c, const ValidatedTaylorModel& y)
{
    ValidatedTaylorModel& x=*this;
    ValidatedTaylorModel r=this->create();
    _sma(r,x,c,y);
    this->swap(r);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error()>=0,*this);
}

void ValidatedTaylorModel::ifma(const ValidatedTaylorModel& x1, const ValidatedTaylorModel& x2)
{
    _mul(*this,x1,x2);
    this->sweep();
    ARIADNE_ASSERT_MSG(this->error()>=0,*this);
}



///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


ValidatedTaylorModel&
ValidatedTaylorModel::unique_sort()
{
    this->_expansion.reverse_lexicographic_sort();

    ValidatedTaylorModel::const_iterator advanced =this->begin();
    ValidatedTaylorModel::const_iterator end =this->end();
    ValidatedTaylorModel::iterator current=this->begin();
    Float te=0.0;
    while(advanced!=end) {
        current->key()=advanced->key();
        VOLATILE Float u=advanced->data();
        VOLATILE Float ml=-advanced->data();
        VOLATILE Float v=advanced->data();
        ++advanced;
        while(advanced!=end && advanced->key()==current->key()) {
            const Float& xv=advanced->data();
            set_rounding_upward();
            u+=xv;
            ml-=xv;
            set_rounding_to_nearest();
            v+=xv;
            ++advanced;
        }
        current->data()=v;
        te+=(u+ml);
        ++current;
    }
    set_rounding_upward();
    this->error()+=te;
    set_rounding_to_nearest();
    this->_expansion.resize(current-this->begin());

    return *this;
}


ValidatedTaylorModel&
ValidatedTaylorModel::sweep()
{
    this->_sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

ValidatedTaylorModel&
ValidatedTaylorModel::sweep(const SweeperInterface& sweeper)
{
    sweeper.sweep(this->_expansion,this->_error);
    return *this;
}



ValidatedTaylorModel&
ValidatedTaylorModel::clobber()
{
    this->_error=0;
    return *this;
}



///////////////////////////////////////////////////////////////////////////////

// Accuracy control

Float
ValidatedTaylorModel::tolerance() const
{
    const ThresholdSweeper* ptr=dynamic_cast<const ThresholdSweeper*>(&static_cast<const SweeperInterface&>(this->_sweeper));
    if(ptr) {
        return ptr->sweep_threshold();
    } else {
        return std::numeric_limits<double>::epsilon();
    }
}



//////////////////////////////////////////////////////////////////////////////

// Exact functions (max, min, abs, neg) and arithmetical functions (sqr, pow)


ValidatedTaylorModel max(const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    Interval xr=x.range();
    Interval yr=y.range();
    if(xr.lower()>=yr.upper()) {
        return x;
    } else if(yr.lower()>=xr.upper()) {
        return y;
    } else {
        return ((x+y)+abs(x-y))/Float(2);
    }
}


ValidatedTaylorModel min(const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    return -max(-x,-y);
}

ValidatedTaylorModel abs(const ValidatedTaylorModel& x) {
    Interval xr=x.range();
    if(xr.lower()>=0.0) {
        return x;
    } else if(xr.upper()<=0.0) {
        return -x;
    } else {
        // Use power series expansion $abs(x)=\sum_{i=0}^{7} p_i x^{2i} \pm e$ for $x\in[-1,+1]$ with
        // p=[0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667] and e=0.035
        // TODO: Find more accurate and stable formula
        static const uint n=7u;
        static const double p[n]={0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667};
        static const double err=0.035;
        ValidatedTaylorModel r(x.argument_size(),x.sweeper());
        Float xmag=mag(xr);
        ValidatedTaylorModel s=x/xmag;
        s=sqr(s);
        r=p[n-1];
        for(uint i=0; i!=(n-1); ++i) {
            uint j=(n-2)-i;
            r=s*r+p[j];
        }
        r+=Interval(-err,+err);
        return r*xmag;
    }
}

ValidatedTaylorModel neg(const ValidatedTaylorModel& x) {
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r.expansion().append(xiter->key(),-xiter->data());
    }
    r.error()=x.error();
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

ValidatedTaylorModel sqr(const ValidatedTaylorModel& x) {
    ValidatedTaylorModel r=x*x;
    return r;
}

ValidatedTaylorModel pow(const ValidatedTaylorModel& x, int n) {
    ValidatedTaylorModel r(x.argument_size(),x.sweeper()); r+=1;
    ValidatedTaylorModel p(x);
    while(n) {
        if(n%2) { r=r*p; }
        p=sqr(p);
        n/=2;
    }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

/* FIXME: The following code does not work with optimisation turned on using gcc */
Interval _range1(const ValidatedTaylorModel& tm) {
    set_rounding_mode(upward);
    VOLATILE Float t=tm.error();
    VOLATILE Float v=0.0;
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            v=iter->data();
        } else {
            t+=abs(iter->data());
        }
    }
    set_rounding_mode(to_nearest);
    return v+Interval(-t,t);
}


Interval _range2(const ValidatedTaylorModel& tm) {
    Interval r(-tm.error(),+tm.error());
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            r+=iter->data();
        } else {
            r+=abs(iter->data())*Interval(-1,1);
        }
    }
    return r;
}

// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
Interval _range3(const ValidatedTaylorModel& tm) {
    const uint as=tm.argument_size();
    Array<Float> linear_terms(as,0.0);
    Array<Float> quadratic_terms(as,0.0);
    Interval r(-tm.error(),+tm.error());
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            r+=iter->data();
        } else if(iter->key().degree()==1) {
            for(uint j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==1) { linear_terms[j]=iter->data(); break; }
            }
        } else if(iter->key().degree()==2) {
            for(uint j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==2) { quadratic_terms[j]=iter->data(); break; }
                if(iter->key()[j]==1) { r+=abs(iter->data())*Interval(-1,1); break; }
            }
        } else {
            r+=abs(iter->data())*Interval(-1,1);
        }
    }
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(uint j=0; j!=as; ++j) {
        const Float& a=quadratic_terms[j];
        const Float& b=linear_terms[j];
        Interval ql=abs(a)*Interval(-1,1) + abs(b)*Interval(-1,+1);
        Interval qf=a*(sqr(Interval(-1,+1)+Interval(b)/(2*a)))-sqr(Interval(b))/(4*a);
        r += intersection(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
    }
    return r;
}



Interval
ValidatedTaylorModel::range() const {
    return Ariadne::_range3(*this);
}


Vector<Interval>
ValidatedTaylorModel::domain() const
{
    return Vector<Interval>(this->argument_size(),Interval(-1,1));
}



//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
class TaylorSeries;


ValidatedTaylorModel
_compose(const TaylorSeries& ts, const ValidatedTaylorModel& tv, double eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    //std::cerr<<"_compose(TaylorSeries,ValidatedTaylorModel,Error)\n";
    //std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    Float& vref=const_cast<Float&>(tv.value());
    Float vtmp=vref;
    vref=0.0;
    ValidatedTaylorModel r(tv.argument_size(),tv.sweeper());
    r+=ts.expansion[ts.expansion.size()-1];
    for(uint i=1; i!=ts.expansion.size(); ++i) {
        //std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
        r.sweep(threshold_sweeper);
    }
    //std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    //std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

ValidatedTaylorModel
compose(const TaylorSeries& ts, const ValidatedTaylorModel& tm)
{
    return _compose(ts,tm,ec);
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
ValidatedTaylorModel
_compose1(const series_function_pointer& fn, const ValidatedTaylorModel& tm, double eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    static const uint DEGREE=18;
    static const double TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tm.value();
    Interval r=tm.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    Float truncation_error_estimate=mag(range_series[d])*pow(mag(r-c),d);
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error_estimate
                     <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    ValidatedTaylorModel x=tm-c;
    ValidatedTaylorModel res(tm.argument_size(),tm.sweeper());
    res+=range_series[d];
    for(uint i=0; i!=d; ++i) {
        //std::cerr<<"i="<<i<<" r="<<res<<"\n";
        res=centre_series[d-i-1]+x*res;
        res.sweep(threshold_sweeper);
    }
    //std::cerr<<"i="<<d<<" r="<<res<<"\n";
    return res;
}

// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. The radius of convergence of this method is still quite low,
// typically only half of the radius of convergence of the power series itself
ValidatedTaylorModel
_compose2(const series_function_pointer& fn, const ValidatedTaylorModel& tm, double eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
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
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    ValidatedTaylorModel x=tm-c;
    ValidatedTaylorModel res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(threshold_sweeper);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
ValidatedTaylorModel
_compose3(const series_function_pointer& fn, const ValidatedTaylorModel& tm, Float eps)
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
    Interval p=pow(e,d-1);
    p=Interval(-p.lower()*e.lower(),p.upper()*e.upper());
    //std::cerr<<"se="<<se<<" e="<<e<<" p="<<p<<std::endl;
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float truncation_error=max(se.lower()*p.lower(),se.upper()*p.upper());
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    ValidatedTaylorModel x=tm;
    ValidatedTaylorModel res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        //res.sweep(eps);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


ValidatedTaylorModel
_compose(const series_function_pointer& fn, const ValidatedTaylorModel& tm, Float eps) {
    //std::cerr<<"_compose(SeriesFunction,ValidatedTaylorModel,Error)\n";
    return _compose3(fn,tm,eps);
}




///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   bounded domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)


ValidatedTaylorModel sqrt(const ValidatedTaylorModel& x) {
    //std::cerr<<"rec(ValidatedTaylorModel)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();

    if(r.lower()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x.range());
    }

    assert(r.lower()>0);
    Float a=(r.lower()+r.upper())/2;
    set_rounding_upward();
    Float eps=(r.upper()-r.lower())/(r.upper()+r.lower());
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=integer_cast<int>((log((1-eps)*x.tolerance())/log(eps)+1));
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    ValidatedTaylorModel y=(x/a)-1.0;
    //std::cerr<<"y="<<y<<std::endl;
    ValidatedTaylorModel z(x.argument_size(),x.sweeper());
    Series<Interval> sqrt_series=Series<Interval>::sqrt(d,Interval(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(uint i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        z.sweep();
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

ValidatedTaylorModel rec(const ValidatedTaylorModel& x) {
    //std::cerr<<"rec(ValidatedTaylorModel)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    if(r.upper()>=0 && r.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"rec(ValidatedTaylorModel x)","x="<<x<<", x.range()="<<x.range());
    }
    Float a=(r.lower()+r.upper())/2;
    set_rounding_upward();
    Float eps=abs((r.upper()-r.lower())/(r.upper()+r.lower()));
    set_rounding_to_nearest();
    assert(eps<1);

    uint d=integer_cast<uint>((log((1-eps)*x.tolerance())/log(eps))+1);

    ValidatedTaylorModel y=1-(x/a);
    ValidatedTaylorModel z(x.argument_size(),x.sweeper());
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

ValidatedTaylorModel log(const ValidatedTaylorModel& x) {
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    if(r.lower()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x.range());
    }
    Float a=(r.lower()+r.upper())/2;
    set_rounding_upward();
    Float eps=(r.upper()-r.lower())/(r.upper()+r.lower());
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=integer_cast<uint>((log((1-eps)*x.tolerance())/log(eps)+1));
    ValidatedTaylorModel y=x/a-1;
    ValidatedTaylorModel z(x.argument_size(),x.sweeper());
    z+=Float(d%2?-1:+1)/d;
    for(uint i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
        z.sweep();
    }
    z=z*y;
    z.sweep();
    Float trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(Interval(a))+trunc_err*Interval(-1,1);
}

static const uint MAXIMUM_DEGREE = 18;

// Use special code to utilise exp(ax+b)=exp(x)^a*exp(b)
ValidatedTaylorModel exp(const ValidatedTaylorModel& x) {
    // FIXME: Truncation error may be incorrect

    // Scale to unit interval
    ValidatedTaylorModel y=x;
    Float xval=x.value();
    y-=xval;
    Float xrad=mag(y.range());
    uint sfp=0; // A number such that 2^sfp>rad(x.range())
    while(Float(1<<sfp)<xrad) { ++sfp; }
    Float sf=1.0/(1<<sfp);
    _scal_exact(y,sf);
    Float yrad=xrad*sf;

    uint degree=MAXIMUM_DEGREE;

    // Since x is in unit domain, truncation error is no worse than maximum omitted term, i.e. xr/fac(d+1)
    ValidatedTaylorModel res(x.argument_size(),x.sweeper());
    res.set_error(pow_up(yrad,degree+1));
    for(uint i=0; i!=degree; ++i) {
        res/=(degree-i);
        res=y*res+1.0;
    }

    // Square r a total of sfp times
    ValidatedTaylorModel square(x.argument_size(),x.sweeper());
    for(uint i=0; i!=sfp; ++i) {
        _mul(square,res,res);
        res.swap(square);
        square.clear();
    }

    // Multiply by exp(xv)
    res*=Ariadne::exp(Interval(xval));

    return res;
    //return _compose(&Series<Interval>::exp,x,x.tolerance());
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form x*f(x^2)
ValidatedTaylorModel sin(const ValidatedTaylorModel& x) {
    // FIXME: Truncation error may be incorrect
    ValidatedTaylorModel y(x.argument_size(),x.sweeper());
    ValidatedTaylorModel s(x.argument_size(),x.sweeper());
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    ValidatedTaylorModel t(x.argument_size(),x.sweeper());

    Float two_pi_approx=2*pi_approx;
    int n=integer_cast<int>(floor(x.value()/two_pi_approx + 0.5));
    y=x-(n*2*pi_ivl);

    if(y.error()>two_pi_approx/2) {
        r.error()=1.0;
    } else {
        _mul(s,y,y);

        int d=(MAXIMUM_DEGREE+3)/2;
        Float srad=mag(s.range());
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute x(1-y/6+y^2/120-y^3/5040+... = x(1-y/6*(1-y/20*(1-y/42*...)
        r=1.0;
        for(int i=0; i!=d; ++i) {
            r/=Float(-2*(d-i)*(2*(d-i)+1));
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
ValidatedTaylorModel cos(const ValidatedTaylorModel& x) {
    // FIXME: Truncation error may be incorrect
    ValidatedTaylorModel y(x.argument_size(),x.sweeper());
    ValidatedTaylorModel s(x.argument_size(),x.sweeper());
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    ValidatedTaylorModel t(x.argument_size(),x.sweeper());

    Float two_pi=2*pi_approx;
    int n=integer_cast<int>(floor(x.value()/two_pi + 0.5));

    y=x-2*n*pi_ivl;

    if(y.error()>two_pi/2) {
        r.error()=1.0;
    } else {
        _mul(s,y,y);

        int d=(MAXIMUM_DEGREE+3)/2;
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

ValidatedTaylorModel tan(const ValidatedTaylorModel& x) {
    return sin(x)*rec(cos(x));
}

ValidatedTaylorModel asin(const ValidatedTaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<ValidatedNumberType>::asin,
                                x.value(),x.range()),x);
}

ValidatedTaylorModel acos(const ValidatedTaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<ValidatedNumberType>::acos,
                                x.value(),x.range()),x);
}

ValidatedTaylorModel atan(const ValidatedTaylorModel& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<ValidatedNumberType>::atan,
                                x.value(),x.range()),x);
}



///////////////////////////////////////////////////////////////////////////////

// Inplace function operators (rescale, restrict, antidifferentiate)


ValidatedTaylorModel& ValidatedTaylorModel::rescale(const Interval& ocd, const Interval& ncd)
{
    ValidatedTaylorModel& x=*this;

    const Float a=ocd.lower(); const Float b=ocd.upper();
    const Float c=ncd.lower(); const Float d=ncd.upper();

    // Scale the interval [a,b] onto [c,d]
    // The function is given by x:-> alpha*x+beta where
    // alpha=(d-c)/(b-a) and beta=(cb-ad)/(b-a)
    // If a==b, return a ValidatedTaylorModel with unbounded error

    ARIADNE_ASSERT_MSG(ocd.radius()>=0,"Illegal scaling from interval "<<ocd<<" with zero radius to interval "<<ncd);
    if(ocd.lower()==ocd.upper()) {
        x.clear();
        x.set_error(+inf);
    } else {
        ValidatedNumberType tmp=1.0/sub_ivl(b,a);
        ValidatedNumberType alpha=sub_ivl(d,c)*tmp;
        ValidatedNumberType beta=(mul_ivl(c,b)-mul_ivl(a,d))*tmp;
        x*=alpha;
        x+=beta;
    }

    return x;
}



// Replace the kth variable x[k] by a*x[k]+b.
ValidatedTaylorModel preaffine(const ValidatedTaylorModel& tm, uint k, const ValidatedNumberType& a, const ValidatedNumberType& b) {
    uint d=tm.degree();
    uint as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    ValidatedTaylorModel r(as,swp);
    r.set_error(tm.error());

    // Create a temporary TaylorModels containing just terms x[k]^i
    Array<ValidatedTaylorModel> atm(d+1,ValidatedTaylorModel(as,swp));
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        MultiIndex a=iter->key();
        const Float& c=iter->data();
        uint ak=a[k];
        a[k]=0;
        atm[ak].expansion().append(a,c);
    }

    ValidatedTaylorModel xk=ValidatedTaylorModel::variable(as,k,swp);

    for(uint i=0; i<=d; ++i) {
        for(uint j=i; j<=d; ++j) {
            ValidatedNumberType c=(Ariadne::bin(j,i)*Ariadne::pow(a,i)*Ariadne::pow(b,j-i));
            r+=c*atm[j];
            atm[j]*=xk;
        }
    }
    return r;
}


ValidatedTaylorModel discard(const ValidatedTaylorModel& tm, Array<uint>& discarded_variables) {
    return recondition(tm,discarded_variables,0u,0u);
}

ValidatedTaylorModel recondition(const ValidatedTaylorModel& tm, Array<uint>& discarded_variables, uint number_of_error_variables) {
    return recondition(tm,discarded_variables,number_of_error_variables,number_of_error_variables);
}

ValidatedTaylorModel recondition(const ValidatedTaylorModel& tm, Array<uint>& discarded_variables, uint number_of_error_variables, uint index_of_error)
{
    for(uint i=0; i!=discarded_variables.size()-1; ++i) {
        ARIADNE_PRECONDITION(discarded_variables[i]<discarded_variables[i+1]);
    }
    ARIADNE_PRECONDITION(discarded_variables[discarded_variables.size()-1]<tm.argument_size());
    ARIADNE_PRECONDITION(index_of_error<=number_of_error_variables);

    const uint number_of_discarded_variables = discarded_variables.size();
    const uint number_of_kept_variables = tm.argument_size() - discarded_variables.size();

    // Make an Array of the variables to be kept
    Array<uint> kept_variables(number_of_kept_variables);
    uint kd=0; uint kk=0;
    for(uint j=0; j!=tm.argument_size(); ++j) {
        if(kd==number_of_discarded_variables || j!=discarded_variables[kd]) {
            kept_variables[kk]=j; ++kk;
        } else {
            ++kd;
        }
    }

    // Construct result and reserve memory
    ValidatedTaylorModel r(number_of_kept_variables+number_of_error_variables,tm.sweeper());
    r.expansion().reserve(tm.number_of_nonzeros()+1u);
    MultiIndex a(number_of_kept_variables+number_of_error_variables);

    // Set the uniform error of the original model
    // If index_of_error == number_of_error_variables, then the error is kept as a uniform error bound
    Float* error_ptr;
    if(number_of_error_variables==index_of_error) {
        error_ptr = &r.error();
    } else {
        a[number_of_kept_variables+index_of_error]=1;
        r.expansion().append(a,tm.error());
        a[number_of_kept_variables+index_of_error]=0;
        error_ptr = &r.begin()->data();
    }

    set_rounding_upward();
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        bool keep=true;
        for(uint k=0; k!=number_of_discarded_variables; ++k) {
            if(iter->key()[discarded_variables[k]]!=0) {
                *error_ptr = add_rnd(*error_ptr,abs(iter->data()));
                keep=false;
                break;
            }
        }
        if(keep) {
            for(uint k=0; k!=number_of_kept_variables; ++k) {
                a[k]=iter->key()[kept_variables[k]];
            }
            r.expansion().append(a,iter->data());
        }
    }
    set_rounding_to_nearest();

    return r;
}


ValidatedTaylorModel restrict(const ValidatedTaylorModel& tm, uint k, const Interval& nd) {
    ARIADNE_ASSERT(k<tm.argument_size());
    ARIADNE_ASSERT(nd.lower()>=-1 && nd.upper()<=+1);
    if(nd.lower()==-1 && nd.upper()==1) {
        return tm;
    } else if(nd.lower()==-1 && nd.upper()==0) {
        return split(tm,k,false);
    } else if(nd.lower()==0 && nd.upper()==1) {
        return split(tm,k,true);
    } else if(nd.lower()==-0.5 && nd.upper()==0.5) {
        return split(tm,k,indeterminate);
    } else {
        ValidatedNumberType a=sub_ivl(nd.upper()/2,nd.lower()/2);
        ValidatedNumberType b=add_ivl(nd.upper()/2,nd.lower()/2);
        return preaffine(tm,k,a,b);
    }
}


ValidatedTaylorModel& ValidatedTaylorModel::restrict(const Vector<Interval>& nd)
{
    ValidatedTaylorModel& x=*this;
    ARIADNE_ASSERT(x.argument_size()==nd.size());
    const uint as=x.argument_size();
    const uint d=x.expansion().degree();
    if(d==0) { return x; }

    Array< Array<ValidatedNumberType> > sf(as,Array<ValidatedNumberType>(d+1u));
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
            Interval ci=Interval(c);
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
void _antidifferentiate1(ValidatedTaylorModel& x, uint k)
{
    ARIADNE_ASSERT(k<x.argument_size());

    //std::cerr<<"xe="<<xe<<"\n";
    set_rounding_mode(upward);
    VOLATILE Float tre=0; // Twice the maximum accumulated roundoff error
    for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const uint c=xiter->key()[k]+1;
        register Float xv=xiter->data();
        register Float mxv=-xv;
        assert(c>0);
        //VOLATILE Float ml,u;
        if(xv>=0) {
            //VOLATILE Float u=xv/c;
            //VOLATILE Float ml=mxv/c;
            tre+=add_rnd(div_rnd(xv,c),div_rnd(mxv,c));
        } else {
            tre+=add_rnd(div_rnd(mxv,c),div_rnd(xv,c));
            //VOLATILE Float u=mxv/c;
            //VOLATILE Float ml=xv/c;
            //te+=(u+ml);
        }
        //std::cerr<<"  te="<<te;;
    }
    //std::cerr<<"  te="<<tre;;
    x.error()+=(tre/2);
    //xe+=te/2;
    //std::cerr<<"xe="<<xe<<"\n";

    set_rounding_to_nearest();
    for(ValidatedTaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
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
void _antidifferentiate2(ValidatedTaylorModel& x, uint k)
{
    ARIADNE_ASSERT(k<x.argument_size());

    //std::cerr<<"xe="<<xe<<"\n";
    set_rounding_mode(upward);
    VOLATILE Float tre=0; // Twice the maximum accumulated roundoff error
    VOLATILE Float u,ml;
    uint c;
    for(ValidatedTaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
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


ValidatedTaylorModel& ValidatedTaylorModel::antidifferentiate(uint k)
{
    _antidifferentiate2(*this,k); return *this;
}

// Compute derivative by computing term-by-term, switching the rounding mode
// May be slow if no assembly rounding mode switching is available.


ValidatedTaylorModel derivative(const ValidatedTaylorModel& x, uint k)
{
    ARIADNE_ASSERT(k<x.argument_size());
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());

    //std::cerr<<"xe="<<xe<<"\n";
    VOLATILE Float tre=0; // Twice the maximum accumulated roundoff error
    Float xv;
    VOLATILE Float u,ml,mxv;
    MultiIndex a(x.argument_size());
    uint c;
    for(ValidatedTaylorModel::iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
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

ValidatedNumberType
evaluate(const ValidatedTaylorModel& tm, const Vector<ValidatedNumberType>& x)
{
    return horner_evaluate(tm.expansion(),x)+ValidatedNumberType(-tm.error(),+tm.error());
}

Vector<ValidatedNumberType>
evaluate(const Vector<ValidatedTaylorModel>& tv, const Vector<ValidatedNumberType>& x)
{
    Vector<ValidatedNumberType> r(tv.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}


template<class T> class Powers {
  public:
    explicit Powers(const T& t) { _values.push_back(t*0+1); _values.push_back(t); }
    explicit Powers(const T& z, const T& t) { _values.push_back(z); _values.push_back(t); }
    const T& operator[](uint i) const { while(_values.size()<=i) { _values.push_back(_values[1]*_values.back()); } return _values[i]; }
  private:
    mutable std::vector<T> _values;
};

template<> class Powers<ValidatedNumberType> {
  public:
    explicit Powers(const ValidatedNumberType& t) { _values.push_back(ValidatedNumberType(1)); _values.push_back(t); }
    explicit Powers(const ValidatedNumberType& z, const ValidatedNumberType& t) { _values.push_back(z); _values.push_back(t); }
    const ValidatedNumberType& operator[](uint i) const {
        while(_values.size()<=i) {
            if(_values.size()%2==0) { _values.push_back(sqr(_values[_values.size()/2])); }
            else { _values.push_back(_values[1]*_values.back()); } }
        return _values[i]; }
  private:
    mutable std::vector<ValidatedNumberType> _values;
};


ValidatedTaylorModel substitute(const ValidatedTaylorModel& x, uint k, const ValidatedTaylorModel& s) {
    ARIADNE_ASSERT(x.argument_size()==s.argument_size()+1u);
    const uint n=s.argument_size();
    Sweeper swp=s.sweeper();
    Vector<ValidatedTaylorModel> y(n+1,ValidatedTaylorModel::zero(n,swp));
    for(uint i=0; i!=n; ++i) {
        y[i]=ValidatedTaylorModel::variable(n,i,swp);
    }
    y[n]=s;
    return compose(x,y);
}

Vector<ValidatedTaylorModel> substitute(const Vector<ValidatedTaylorModel>& x, uint k, const ValidatedTaylorModel& s) {
    const uint n=s.argument_size();
    Sweeper swp=s.sweeper();
    Vector<ValidatedTaylorModel> y(n+1,ValidatedTaylorModel::zero(n,swp));
    for(uint i=0; i!=n; ++i) {
        y[i]=ValidatedTaylorModel::variable(n,i,swp);
    }
    y[n]=s;
    return compose(x,y);
}


ValidatedTaylorModel
partial_evaluate(const ValidatedTaylorModel& x, uint k, ExactNumberType c)
{
    return partial_evaluate(x,k,ValidatedNumberType(c));
}


ValidatedTaylorModel
partial_evaluate(const ValidatedTaylorModel& x, uint k, ValidatedNumberType c)
{
    ValidatedTaylorModel r(x.argument_size()-1,x.sweeper());
    MultiIndex ra(r.argument_size());
    if(c==0) {
        for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            MultiIndex::index_type xak=xa[k];
            if(xak==0) {
                const Float& xv=xiter->data();
                for(uint i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(uint i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
                r.expansion().append(ra,xv);
            }
        }
        r.set_error(x.error());
    } else if(c==1) {
        ValidatedTaylorModel s(x.argument_size()-1,x.sweeper());
        Array<ValidatedTaylorModel> p(x.degree()+1,ValidatedTaylorModel(x.argument_size()-1,x.sweeper()));

        for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const Float& xv=xiter->data();
            MultiIndex::index_type xak=xa[k];
            for(uint i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(uint i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }

        r=p[0];
        r.set_error(x.error());
        for(uint i=1; i!=p.size(); ++i) {
            _add(s,r,p[i]);
            r.swap(s);
            s.clear();
        }
    } else {
        ValidatedTaylorModel s(x.argument_size()-1,x.sweeper());
        Array<ValidatedTaylorModel> p(x.degree()+1,ValidatedTaylorModel(x.argument_size()-1,x.sweeper()));

        Powers<ValidatedNumberType> cpowers(c);

        for(ValidatedTaylorModel::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const Float& xv=xiter->data();
            MultiIndex::index_type xak=xa[k];
            for(uint i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(uint i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().append(ra,xv);
        }
        for(uint i=1; i!=p.size(); ++i) {
            p[i]*=cpowers[i];
        }

        r=p[0];
        r.set_error(x.error());
        for(uint i=1; i!=p.size(); ++i) {
            _add(s,r,p[i]);
            r.swap(s);
            s.clear();
        }
    }

    return r;
}

Vector<ValidatedTaylorModel>
partial_evaluate(const Vector<ValidatedTaylorModel>& tv, uint k, ExactNumberType c)
{
    // FIXME: Fails if tv.size()==0
    Vector<ValidatedTaylorModel> r(tv.size(),ValidatedTaylorModel::zero(tv.zero_element().argument_size(),tv.zero_element().sweeper()));
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=partial_evaluate(tv[i],k,c);
    }
    return r;
}

Vector<ValidatedTaylorModel>
partial_evaluate(const Vector<ValidatedTaylorModel>& tv, uint k, ValidatedNumberType c)
{
    Vector<ValidatedTaylorModel> r(tv.size(),ValidatedTaylorModel::zero(tv.zero_element().argument_size(),tv.zero_element().sweeper()));
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=partial_evaluate(tv[i],k,c);
    }
    return r;
}


pair<ValidatedTaylorModel,ValidatedTaylorModel>
split(const ValidatedTaylorModel& tv, uint j)
{
    // TODO: improve efficiency of implementation
    uint as=tv.argument_size();
    Sweeper swp=tv.sweeper();
    Vector<ValidatedTaylorModel> s=ValidatedTaylorModel::variables(as,swp);
    s[j]=ValidatedTaylorModel::scaling(as,j,Interval(-1,0),swp);
    ValidatedTaylorModel r1=compose(tv,s);
    r1.set_sweeper(tv.sweeper());
    s[j]=ValidatedTaylorModel::scaling(as,j,Interval(0,+1),swp);
    ValidatedTaylorModel r2=compose(tv,s);
    r2.set_sweeper(tv.sweeper());
    return make_pair(r1,r2);
}




ValidatedTaylorModel
_split1(const ValidatedTaylorModel& tm, uint k, tribool b)
{
    const uint deg=tm.degree();
    const uint as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    ValidatedTaylorModel r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(ValidatedTaylorModel::iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const uchar ak=iter->key()[k];
        Float& c=iter->data();
        c/=(1<<ak);
    }

    int tr=( indeterminate(b) ? 0 : definitely(b) ? +1 : -1 );

    if(tr==0) { return r; }

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<ValidatedTaylorModel> ary(deg+1,ValidatedTaylorModel(as,swp));
    for(ValidatedTaylorModel::const_iterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->key();
        const Float& c=iter->data();
        uchar ak=a[k];
        a[k]=0u;
        ary[ak].expansion().append(a,c);
    }

    Float re=r.error();
    r.clear();
    r.set_error(re);

    for(uint i=0; i<=deg; ++i) {
        for(uint j=i; j<=deg; ++j) {
            int sf=bin(j,i);
            if(tr==-1 && (j-i)%2==1) { sf=-sf; }
            r+=ary[j]*sf;
            for(ValidatedTaylorModel::iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->key()[k];
            }
         }
    }

    return r;
}


ValidatedTaylorModel
split(const ValidatedTaylorModel& tm, uint j, tribool b)
{
    return _split1(tm,j,b);
}

ValidatedTaylorModel
unscale(const ValidatedTaylorModel& tv, const Interval& ivl)
{
    // Scale tv so that the interval ivl maps into [-1,1]
    // The result is given by  (tv-c)*s where c is the centre
    // and s the reciprocal of the radius of ivl

    // If the radius of ivl is less than a constant, then
    // there are three possible policies. We can either map to zero or to
    // everything, or map a constant to zero and other models to everything.
    // The motivation for mapping to zero is that the domain is
    // restricted to the point and this is
    // The motivation for mapping to everything is that any function on the
    // resulting interval should be independent of the unneeded component

    const ExactNumberType& l=ivl.lower();
    const ExactNumberType& u=ivl.upper();
    ARIADNE_ASSERT_MSG(l<=u,"Cannot unscale ValidatedTaylorModel "<<tv<<" from empty interval "<<ivl);

    if(l==u) {
        ValidatedTaylorModel r=tv.create();
        r+=ivl;
        // Uncomment out line below to make unscaling to a singleton interval undefined
        //r.set_error(+inf);
        return r;
    } else {
        ValidatedTaylorModel r=tv;
        ValidatedNumberType c=ValidatedNumberType(l/2)+ValidatedNumberType(u/2);
        ValidatedNumberType s=2/(ValidatedNumberType(u)-ValidatedNumberType(l));
        r-=c;
        r*=s;

        return r;
    }
}


ValidatedTaylorModel
rescale(const ValidatedTaylorModel& tv, const Interval& ivl)
{
    // Scale tv so that the interval [-1,1] maps into ivl
    // The result is given by  (tv*s)+c where c is the centre
    // and r the radius of ivl
    const Float& l=ivl.lower();
    const Float& u=ivl.upper();

    ValidatedTaylorModel r=tv;
    ValidatedNumberType c=add_ivl(l/2,u/2);
    ValidatedNumberType s=sub_ivl(u/2,l/2);
    r*=s;
    r+=c;

    return r;
}



Float ValidatedTaylorModel::radius() const {
    set_rounding_mode(upward);
    Float r=this->error();
    for(ValidatedTaylorModel::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=abs(iter->data());
        }
    }
    set_rounding_mode(to_nearest);
    return r;
}

Float ValidatedTaylorModel::norm() const {
    set_rounding_mode(upward);
    Float r=this->error();
    for(ValidatedTaylorModel::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=abs(iter->data());
    }
    set_rounding_mode(to_nearest);
    return r;
}

Float norm(const ValidatedTaylorModel& tv) {
    return tv.norm();
}

bool
refines(const ValidatedTaylorModel& tv1, const ValidatedTaylorModel& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    ValidatedTaylorModel d=tv2;
    d.error()=0.0;
    d-=tv1;
    return norm(d) <= tv2.error();
}


bool
disjoint(const ValidatedTaylorModel& tv1, const ValidatedTaylorModel& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    Float tv1e=tv1.error();
    Float tv2e=tv2.error();
    const_cast<ValidatedTaylorModel&>(tv1).error()=0.0;
    const_cast<ValidatedTaylorModel&>(tv2).error()=0.0;
    ValidatedTaylorModel d=tv2-tv1;
    const_cast<ValidatedTaylorModel&>(tv1).error()=tv1e;
    const_cast<ValidatedTaylorModel&>(tv2).error()=tv2e;
    //std::cerr<<"\ntv1="<<tv1<<"\ntv2="<<tv2<<"\nd="<<d
    //         <<"\n|d|="<<norm(d)<<" |e|="<<add_up(tv1.error(),tv2.error())<<"\n\n";
    return norm(d)>add_up(tv1.error(),tv2.error());
}


ValidatedTaylorModel embed(const ValidatedTaylorModel& x, uint as)
{
    return ValidatedTaylorModel(embed(0u,x.expansion(),as),x.error(),x.sweeper());
}

ValidatedTaylorModel embed(uint as, const ValidatedTaylorModel& x)
{
    return ValidatedTaylorModel(embed(as,x.expansion(),0u),x.error(),x.sweeper());
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

std::ostream&
operator<<(std::ostream& os, const ValidatedTaylorModel& tm) {
    // Set the variable names to be 'parameter' s0,s1,..
    Array<std::string> variable_names(tm.argument_size());
    for(uint j=0; j!=tm.argument_size(); ++j) {
        std::stringstream sstr;
        sstr << 's' << j;
        variable_names[j]=sstr.str();
    }

    //os << "ValidatedTaylorModel";
    os << "TM["<<tm.argument_size()<<"](";
    Expansion<Float> e=tm.expansion();
    e.graded_sort();
    e.write(os,variable_names);
    return os << "+/-" << tm.error() << ")";
}


///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


Vector<ValidatedTaylorModel> ValidatedTaylorModel::zeros(uint rs, uint as, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(rs,ValidatedTaylorModel::zero(as,swp));
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::constants(uint as, const Vector<ExactNumberType>& c, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(c.size(),ValidatedTaylorModel::zero(as,swp));
    for(uint i=0; i!=c.size(); ++i) {
        result[i]=ValidatedTaylorModel::constant(as,c[i],swp);
    }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::constants(uint as, const Vector<ValidatedNumberType>& c, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(c.size(),ValidatedTaylorModel::zero(as,swp));
    for(uint i=0; i!=c.size(); ++i) {
        result[i]=ValidatedTaylorModel::constant(as,c[i],swp);
    }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::variables(uint as, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(as,ValidatedTaylorModel::zero(as,swp));
    for(uint i=0; i!=as; ++i) { result[i]=ValidatedTaylorModel::variable(as,i,swp); }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::scalings(const Vector<Interval>& d, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(d.size(),ValidatedTaylorModel::zero(d.size(),swp));
    for(uint i=0; i!=d.size(); ++i) {
        result[i]=ValidatedTaylorModel::scaling(d.size(),i,d[i],swp);
    }
    return result;
}

Vector<ValidatedTaylorModel> ValidatedTaylorModel::unscalings(const Vector<Interval>& cd, Sweeper swp)
{
    Vector<ValidatedTaylorModel> result(cd.size(),ValidatedTaylorModel::zero(cd.size(),swp));
    for(uint i=0; i!=cd.size(); ++i) {
        result[i]=ValidatedTaylorModel::unscaling(cd.size(),i,cd[i],swp);
    }
    return result;
}



///////////////////////////////////////////////////////////////////////////////

// Vector-valued versions of scalar operators




Vector<ValidatedTaylorModel> embed(const Vector<ValidatedTaylorModel>& x, uint as)
{
    Vector<ValidatedTaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=embed(x[i],as);
    }
    return r;
}

Vector<ValidatedTaylorModel> embed(uint as, const Vector<ValidatedTaylorModel>& x)
{
    Vector<ValidatedTaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=embed(as,x[i]);
    }
    return r;
}

Vector<ValidatedTaylorModel> combine(const Vector<ValidatedTaylorModel>& x1, const Vector<ValidatedTaylorModel>& x2) {
    return join(embed(x1,x2.zero_element().argument_size()),embed(x1.zero_element().argument_size(),x2));
}

Vector<ValidatedTaylorModel> combine(const Vector<ValidatedTaylorModel>& x1, const ValidatedTaylorModel& x2) {
    return join(embed(x1,x2.argument_size()),embed(x1.zero_element().argument_size(),x2));
}

Vector<ValidatedTaylorModel> combine(const ValidatedTaylorModel& x1, const Vector<ValidatedTaylorModel>& x2) {
    return join(embed(x1,x2.zero_element().argument_size()),embed(x1.argument_size(),x2));
}

Vector<ValidatedTaylorModel> combine(const ValidatedTaylorModel& x1, const ValidatedTaylorModel& x2) {
    return join(embed(x1,x2.argument_size()),embed(x1.argument_size(),x2));
}

bool
refines(const Vector<ValidatedTaylorModel>& tv1, const Vector<ValidatedTaylorModel>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(uint i=0; i!=tv1.size(); ++i) {
        if(!refines(tv1[i],tv2[i])) { return false; }
    }
    return true;
}

bool
disjoint(const Vector<ValidatedTaylorModel>& tv1, const Vector<ValidatedTaylorModel>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(uint i=0; i!=tv1.size(); ++i) {
        if(disjoint(tv1[i],tv2[i])) { return true; }
    }
    return false;
}

pair< Vector<ValidatedTaylorModel>, Vector<ValidatedTaylorModel> >
split(const Vector<ValidatedTaylorModel>& tv, uint j)
{
    Vector<ValidatedTaylorModel> r1(tv.size());
    Vector<ValidatedTaylorModel> r2(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        make_lpair(r1[i],r2[i])=split(tv[i],j);
    }
    return make_pair(r1,r2);
}

Vector<ValidatedTaylorModel>
split(const Vector<ValidatedTaylorModel>& tv, uint j, bool h)
{
    Vector<ValidatedTaylorModel> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=split(tv[i],j,h);
    }
    return r;
}

Vector<ValidatedTaylorModel>
unscale(const Vector<ValidatedTaylorModel>& tvs, const Vector<Interval>& ivls)
{
    Vector<ValidatedTaylorModel> r(tvs.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=unscale(tvs[i],ivls[i]);
    }
    return r;
}

ValidatedTaylorModel
rescale(const ValidatedTaylorModel& tv, const Interval& old_codom, const Interval& new_codom)
{
    ValidatedTaylorModel res(tv); res.rescale(old_codom,new_codom); return res;
}

ValidatedTaylorModel&
rescale(ValidatedTaylorModel& tv, const Interval& old_codom, const Interval& new_codom)
{
    tv.rescale(old_codom,new_codom); return tv;
}

Vector<ValidatedTaylorModel>
rescale(const Vector<ValidatedTaylorModel>& tvs, const Vector<Interval>& old_codom, const Vector<Interval>& new_codom)
{
    Vector<ValidatedTaylorModel> r(tvs.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=rescale(tvs[i],old_codom[i],new_codom[i]);
    }
    return r;
}

Vector<ValidatedTaylorModel>&
rescale(Vector<ValidatedTaylorModel>& tvs, const Vector<Interval>& old_codom, const Vector<Interval>& new_codom)
{
    for(uint i=0; i!=tvs.size(); ++i) {
        rescale(tvs[i],old_codom[i],new_codom[i]);
    }
    return tvs;
}

Vector<ValidatedTaylorModel>
scale(const Vector<ValidatedTaylorModel>& tvs, const Vector<Interval>& new_codom)
{
    Vector<ValidatedTaylorModel> r(tvs);
    for(uint i=0; i!=tvs.size(); ++i) {
        r[i].rescale(Interval(-1,+1),new_codom[i]);
    }
    return r;
}

ValidatedTaylorModel antiderivative(const ValidatedTaylorModel& x, uint k) {
    ValidatedTaylorModel r(x);
    r.antidifferentiate(k);
    return r;
}

Vector<ValidatedTaylorModel> antiderivative(const Vector<ValidatedTaylorModel>& x, uint k) {
    Vector<ValidatedTaylorModel> r(x);
    for(uint i=0; i!=x.size(); ++i) {
        r[i].antidifferentiate(k);
    }
    return r;
}

ValidatedTaylorModel antiderivative(const ValidatedTaylorModel& x, const Interval& dk, uint k) {
    ValidatedNumberType dkr=rad_ivl(dk.lower(),dk.upper());
    ValidatedTaylorModel r(x);
    r.antidifferentiate(k);
    r*=dkr;
    return r;
}

Vector<ValidatedTaylorModel> antiderivative(const Vector<ValidatedTaylorModel>& x, const Interval& dk, uint k) {
    Vector<ValidatedTaylorModel> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=antiderivative(x[i],dk,k);
    }
    return r;
}

ValidatedTaylorModel intersection(const ValidatedTaylorModel& x, const ValidatedTaylorModel& y) {
    ValidatedTaylorModel r(x.argument_size(),x.sweeper());
    Float twice_max_error=0.0;

    const Float& xe=x.error();
    const Float& ye=y.error();
    VOLATILE Float rv,xv,yv,xu,yu,mxl,myl,u,ml;
    //const MultiIndex* aptr;
    MultiIndex a;

    ValidatedTaylorModel::const_iterator xiter=x.begin();
    ValidatedTaylorModel::const_iterator yiter=y.begin();
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
            ARIADNE_THROW(IntersectionException,"intersection(ValidatedTaylorModel,ValidatedTaylorModel)",x<<" and "<<y<<" are disjoint.");
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

Vector<Interval> ranges(const Vector<ValidatedTaylorModel>& f) {
    Vector<Interval> r(f.size()); for(uint i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

inline Vector<ValidatedTaylorModel>& clobber(Vector<ValidatedTaylorModel>& h) {
    for(uint i=0; i!=h.size(); ++i) { h[i].set_error(0.0); } return h; }

inline Vector<Float> errors(const Vector<ValidatedTaylorModel>& h) {
    Vector<Float> e(h.size()); for(uint i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

inline Vector<Float> norms(const Vector<ValidatedTaylorModel>& h) {
    Vector<Float> r(h.size()); for(uint i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return r; }

Float norm(const Vector<ValidatedTaylorModel>& h) {
    return norm(norms(h));
}

///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumberType>
jacobian(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumberType>& x)
{
    Vector< Differential<ValidatedNumberType> > dx=Differential<ValidatedNumberType>::variables(1u,x);
    Vector< Differential<ValidatedNumberType> > df(f.size(),x.size(),1u);
    for(uint i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    return jacobian(df);
}

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumberType>
jacobian2(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumberType>& x)
{
    Vector< Differential<ValidatedNumberType> > dx(x.size());
    for(uint i=0; i!=x.size()-f.size(); ++i) {
        dx[i]=Differential<ValidatedNumberType>::constant(f.size(),1u,x[i]); }
    for(uint i=0; i!=f.size(); ++i) {
        uint j=i+(x.size()-f.size());
        dx[j]=Differential<ValidatedNumberType>::variable(f.size(),1u,x[j],i); }
    Vector< Differential<ValidatedNumberType> > df(f.size());
    for(uint i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumberType> J=jacobian(df);
    return J;
}

/*
// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumberType>
jacobian(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumberType>& d)
{
    uint rs=f.size();
    uint as=f.zero_element().argument_size();
    Matrix<ValidatedNumberType> J(rs,as);
    for(uint i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            const MultiIndex& a=iter->key();
            const Float& x=iter->data();
            for(uint k=0; k!=as; ++k) {
                const uint c=a[k];
                if(c>0) {
                    if(iter->key().degree()==1) {
                        J[i][k]+=x;
                    }
                    else {
                        ValidatedNumberType p(-x,+x);
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
jacobian_value(const Vector<ValidatedTaylorModel>& f)
{
    uint rs=f.size();
    uint as=f.zero_element().argument_size();
    Matrix<Float> J(rs,as);
    MultiIndex a(as);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=as; ++j) {
            a[j]=1; const Float x=f[i][a]; J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<Float>
jacobian2_value(const Vector<ValidatedTaylorModel>& f)
{
    const uint rs=f.size();
    const uint fas=f.zero_element().argument_size();
    const uint has=fas-rs;
    Matrix<Float> J(rs,rs);
    MultiIndex a(fas);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=rs; ++j) {
            a[has+j]=1; const Float x=f[i][a]; J[i][j]=x; a[has+j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
Matrix<Interval>
jacobian_range(const Vector<ValidatedTaylorModel>& f)
{
    uint rs=f.size();
    uint as=f.zero_element().argument_size();
    Matrix<Interval> J(rs,as);
    for(uint i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=as; ++k) {
                const uint c=iter->key()[k];
                if(c>0) {
                    const Float& x=iter->data();
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
jacobian2_range(const Vector<ValidatedTaylorModel>& f)
{
    uint rs=f.size();
    uint fas=f.zero_element().argument_size();
    uint has=fas-rs;
    Matrix<Interval> J(rs,rs);
    for(uint i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=rs; ++k) {
                const uint c=iter->key()[has+k];
                if(c>0) {
                    const Float& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=Interval(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}


ValidatedNumberType
jacobian2(const ValidatedTaylorModel& f, const Vector<ValidatedNumberType>& x)
{
    return jacobian2(Vector<ValidatedTaylorModel>(1u,f),x)[0][0];
}

Float
jacobian2_value(const ValidatedTaylorModel& f)
{
    return jacobian2_value(Vector<ValidatedTaylorModel>(1u,f))[0][0];
}

Interval jacobian2_range(const ValidatedTaylorModel& f) {
    return jacobian2_range(Vector<ValidatedTaylorModel>(1u,f))[0][0];
}

///////////////////////////////////////////////////////////////////////////////

// Vector operators (evaluate, compose




// Compose by computing each and every term individually without caching
// Easy to implement, but far too slow
inline
Vector<ValidatedTaylorModel>
_compose1(const Vector<ValidatedTaylorModel>& x,
          const Vector<ValidatedTaylorModel>& ys)
{
    //std::cerr<<"compose1"<<std::endl;
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(ys.size()==x.zero_element().argument_size());
    for(uint i=0; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x.zero_element().argument_size()); }
    for(uint i=0; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys.zero_element().argument_size(),"ys="<<ys); }

    uint as=ys.zero_element().argument_size();
    Sweeper swp=ys.zero_element().sweeper();

    Vector<ValidatedTaylorModel> r(x.size(),ValidatedTaylorModel(as,swp));
    ValidatedTaylorModel t(as,swp);
    for(uint i=0; i!=x.size(); ++i) {
        r[i].set_error(x[i].error());
        for(ValidatedTaylorModel::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            t=iter->data();
            for(uint j=0; j!=iter->key().size(); ++j) {
                ValidatedTaylorModel p=pow(ys[j],iter->key()[j]);
                t=t*p;
            }
            r[i]+=t;
        }
    }

    return r;
}


inline
Vector<ValidatedTaylorModel>
_compose2(const Vector<ValidatedTaylorModel>& x,
          const Vector<ValidatedTaylorModel>& ys)
{
    //std::cerr<<"compose2"<<std::endl;
    uint yrs=ys.size();
    uint xas=ys.size();
    uint as=ys.zero_element().argument_size();
    Sweeper sweeper=ys.zero_element().sweeper();

    Array<uchar> max_power(ys.size());
    for(uint j=0; j!=ys.size(); ++j) { max_power[j]=1; }

    for(uint i=0; i!=x.size(); ++i) {
        for(ValidatedTaylorModel::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            assert(xas==iter->key().size());
            for(uint j=0; j!=iter->key().size(); ++j) {
                max_power[j]=max(max_power[j],iter->key()[j]);
            }
        }
    }

    uchar max_max_power = 1;
    for(uint j=0; j!=ys.size(); ++j) {
        max_max_power = max(max_max_power,max_power[j]);
    }

    Array< Array< ValidatedTaylorModel > > powers(yrs, Array<ValidatedTaylorModel>(max_max_power+1,ValidatedTaylorModel::zero(as,sweeper)));
    for(uint j=0; j!=yrs; ++j) {
        powers[j][0]=ys[j]*0;
        powers[j][1]=ys[j];
        for(uint k=2; k!=powers[j].size(); ++k) {
            powers[j][k]=powers[j][k/2]*powers[j][(k+1)/2];
        }
    }

    Vector<ValidatedTaylorModel> r(x.size(),ValidatedTaylorModel(as,sweeper));
    ValidatedTaylorModel t(as,sweeper);
    MultiIndex a;
    Float c;
    for(uint i=0; i!=x.size(); ++i) {
        r[i].set_error(x[i].error());
        for(ValidatedTaylorModel::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
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

    return r;
}

Vector<ValidatedTaylorModel>
_compose(const Vector<ValidatedTaylorModel>& x,
         const Vector<ValidatedTaylorModel>& ys)
{
    ARIADNE_ASSERT_MSG(x.size()>0,"x="<<x<<", ys="<<ys);
    ARIADNE_ASSERT_MSG(ys.size()==x[0].argument_size(),"x="<<x<<", ys="<<ys);
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }
    for(uint i=1; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys[0].argument_size(),"ys="<<ys); }

    return _compose2(x,ys);
}


Vector<ValidatedTaylorModel>
unchecked_compose(const Vector<ValidatedTaylorModel>& x,
                  const Vector<ValidatedTaylorModel>& ys)
{
    return _compose(x,ys);
}

ValidatedTaylorModel
unchecked_compose(const ValidatedTaylorModel& x,
                  const Vector<ValidatedTaylorModel>& ys)
{
    return _compose(Vector<ValidatedTaylorModel>(1u,x),ys)[0];
}

Vector<ValidatedTaylorModel>
compose(const Vector<ValidatedTaylorModel>& x,
        const Vector<ValidatedTaylorModel>& ys)
{
    return _compose(x,ys);
}

ValidatedTaylorModel
compose(const ValidatedTaylorModel& x,
        const Vector<ValidatedTaylorModel>& ys)
{
    return Ariadne::power_evaluate(x.expansion(),ys)+ValidatedNumberType(-x.error(),+x.error());
    return _compose(Vector<ValidatedTaylorModel>(1u,x),ys)[0];
}

ValidatedTaylorModel
compose(const ValidatedTaylorModel& x,
        const Vector<Interval>& d,
        const Vector<ValidatedTaylorModel>& y)
{
    return _compose(Vector<ValidatedTaylorModel>(1u,x),unscale(y,d))[0];
}

Vector<ValidatedTaylorModel>
compose(const Vector<ValidatedTaylorModel>& x,
        const Vector<Interval>& d,
        const Vector<ValidatedTaylorModel>& y)
{
    return _compose(x,unscale(y,d));
}

Vector<ValidatedTaylorModel>
unchecked_compose(const Vector<ValidatedTaylorModel>& x,
        const Vector<Interval>& d,
        const Vector<ValidatedTaylorModel>& y)
{
    return _compose(x,unscale(y,d));
}



template TaylorModel<ValidatedTag> rec(const TaylorModel<ValidatedTag>&);
template TaylorModel<ValidatedTag> sqrt(const TaylorModel<ValidatedTag>&);
template TaylorModel<ValidatedTag> exp(const TaylorModel<ValidatedTag>&);
template TaylorModel<ValidatedTag> log(const TaylorModel<ValidatedTag>&);
template TaylorModel<ValidatedTag> sin(const TaylorModel<ValidatedTag>&);
template TaylorModel<ValidatedTag> cos(const TaylorModel<ValidatedTag>&);
template TaylorModel<ValidatedTag> tan(const TaylorModel<ValidatedTag>&);







ApproximateTaylorModel::TaylorModel(uint as)
    : _expansion(as), _sweeper()
{
}

ApproximateTaylorModel::TaylorModel(uint as, Sweeper swp)
    : _expansion(as), _sweeper(swp)
{
}

void ApproximateTaylorModel::iadd(const Float& c)
{
    // Compute self+=c
    ApproximateTaylorModel& x=*this;
    if(c==0) { return; }
    if(x._expansion.empty()) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else if((x._expansion.end()-1)->key().degree()>0) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else {
        Float& rv=(x._expansion.end()-1)->data();
        rv+=c;
    }
}

void ApproximateTaylorModel::imul(const Float& c)
{
    // Compute self*=c
    if(c==0) { this->clear(); return; }
    if(c==1) { return; }
    for(ExpansionType::iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        iter->data() *= c;
    }
}



void ApproximateTaylorModel::isma(const Float& c, const ApproximateTaylorModel& y)
{
    const ApproximateTaylorModel& x=*this;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<", y="<<y);
    ApproximateTaylorModel r(x.argument_size());

    const_iterator xiter=x._expansion.begin();
    const_iterator yiter=y._expansion.begin();
    while(xiter!=x._expansion.end() && yiter!=y._expansion.end()) {
        if(xiter->key()<yiter->key()) {
            r._expansion.append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r._expansion.append(yiter->key(),c*yiter->data());
            ++yiter;
        } else {
            r._expansion.append(xiter->key(),xiter->data()+c*yiter->data());
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x._expansion.end()) {
        r._expansion.append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y._expansion.end()) {
        r._expansion.append(yiter->key(),c*yiter->data());
        ++yiter;
    }

    this->swap(r);
}


void ApproximateTaylorModel::ifma(const ApproximateTaylorModel& x, const ApproximateTaylorModel& y)
{
    ApproximateTaylorModel& r(*this);
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<",y="<<y);

    ApproximateTaylorModel t(x.argument_size());
    MultiIndex sa;
    for(ExpansionType::const_iterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        ExpansionType::const_iterator yiter=y._expansion.begin();
        ExpansionType::const_iterator riter=r._expansion.begin();
        while(riter!=r._expansion.end() && yiter!=y._expansion.end()) {
            const MultiIndex& ra=riter->key();
            sa=xiter->key()+yiter->key();
            if(sa==ra) {
                t._expansion.append(sa,riter->data()+xiter->data()*yiter->data());
                ++riter; ++yiter;
            } else if(sa<ra) {
                t._expansion.append(sa,xiter->data()*yiter->data());
                ++yiter;
            } else {
                t._expansion.append(ra,riter->data());
                ++riter;
            }
        }
        while(riter!=r._expansion.end()) {
            t._expansion.append(riter->key(),riter->data());
            ++riter;
        }
        while(yiter!=y._expansion.end()) {
            t._expansion.append(xiter->key(),yiter->key(),xiter->data()*yiter->data());
            ++yiter;

        }

        r._expansion.swap(t._expansion);
        t._expansion.clear();
    }
}




Float ApproximateTaylorModel::norm() const {
    return (*this)[MultiIndex(this->argument_size())];
}

Float ApproximateTaylorModel::average() const {
    return (*this)[MultiIndex(this->argument_size())];
}

Float ApproximateTaylorModel::radius() const {
    Float r=0.0;
    for(Expansion<Float>::const_iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=abs(iter->data());
        }
    }
    return r;
}

Interval ApproximateTaylorModel::range() const {
    Float rad=this->radius(); return this->average()+Interval(-rad,+rad);
}

Float ApproximateTaylorModel::tolerance() const {
    return dynamic_cast<const ThresholdSweeper&>(static_cast<const SweeperInterface&>(this->_sweeper)).sweep_threshold();
}


std::ostream& ApproximateTaylorModel::write(std::ostream& os) const {
    return os << "TM["<<this->argument_size()<<"](" << this->_expansion << ")";
}


template TaylorModel<ApproximateTag> neg(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> rec(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> sqr(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> pow(const TaylorModel<ApproximateTag>&, int);
template TaylorModel<ApproximateTag> sqrt(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> exp(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> log(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> sin(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> cos(const TaylorModel<ApproximateTag>&);
template TaylorModel<ApproximateTag> tan(const TaylorModel<ApproximateTag>&);


std::ostream& ApproximateTaylorModel::str(std::ostream& os) const {
    return os << this->_expansion;
}
} //namespace Ariadne


