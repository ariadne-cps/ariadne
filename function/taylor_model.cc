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

#include "numeric/numeric.h"
#include "config.h"

#include <iomanip>
#include <limits>

#include "numeric/rounding.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/expansion.h"
#include "algebra/series.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"
#include "function/taylor_series.h"
#include "function/function.h"
#include "utility/exceptions.h"

#include "algebra/algebra_mixin.tcc"

#define VOLATILE ;
#include "algebra/multi_index-noaliasing.h"
#include "function/function_mixin.h"
#include "algebra/vector.h"

namespace Ariadne {

template class Series<ValidatedFloat>;

namespace {

static const double HALF_MACHINE_EPSILON = 2.2204460492503131e-16/2;

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }

} // namespace




Vector<ValidatedNumber> unscale(const Vector<ValidatedNumber>& x, const Vector<ExactInterval>& d) {
    Vector<ValidatedNumber> r(x);
    for(SizeType i=0; i!=r.size(); ++i) {
        if(d[i].lower().raw()==d[i].upper().raw()) {
            if(x[i]==d[i].midpoint()) {
                r[i]=ValidatedNumber(0.0,0.0);
            } else {
                r[i]=ValidatedNumber(-inf,+inf);
            }
        } else {
            //r[i]=(2*r[i]-add_ivl(d[i].lower().raw(),d[i].upper().raw()))/sub_ivl(d[i].upper().raw(),d[i].lower().raw());
            r[i]=(2*r[i]-(d[i].lower()+d[i].upper()))/(d[i].upper()-d[i].lower());
        }
    }
    return r;
}



TaylorModel<ValidatedFloat>::TaylorModel()
    : _expansion(0), _error(0), _sweeper()
{ }


TaylorModel<ValidatedFloat>::TaylorModel(SizeType as, Sweeper swp)
    : _expansion(as), _error(0), _sweeper(swp)
{
}

TaylorModel<ValidatedFloat>::TaylorModel(const Expansion<CoefficientType>& f, const ErrorType& e, Sweeper swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->cleanup();
}

TaylorModel<ValidatedFloat>::TaylorModel(const Expansion<Float>& f, const Float& e, Sweeper swp)
    : TaylorModel(reinterpret_cast<Expansion<CoefficientType>const&>(f),reinterpret_cast<ErrorType const&>(e),swp)
{
}

TaylorModel<ValidatedFloat> TaylorModel<ValidatedFloat>::scaling(SizeType as, SizeType j, const ExactInterval& codom, Sweeper swp) {
    TaylorModel<ValidatedFloat> r(as,swp);
    r.set_gradient(j,1);
    r*=codom.radius();
    r+=codom.midpoint();
    return r;
}

TaylorModel<ValidatedFloat> TaylorModel<ValidatedFloat>::create() const {
    return TaylorModel<ValidatedFloat>(this->argument_size(),this->_sweeper);
}

TaylorModel<ValidatedFloat> TaylorModel<ValidatedFloat>::create_zero() const {
    return TaylorModel<ValidatedFloat>(this->argument_size(),this->_sweeper);
}

TaylorModel<ValidatedFloat> TaylorModel<ValidatedFloat>::create_constant(NumericType c) const {
    return TaylorModel<ValidatedFloat>::constant(this->argument_size(),c,this->_sweeper);
}

TaylorModel<ValidatedFloat> TaylorModel<ValidatedFloat>::create_coordinate(SizeType j) const {
    ARIADNE_PRECONDITION(j<this->argument_size());
    TaylorModel<ValidatedFloat> r(this->argument_size(),this->_sweeper);
    r._expansion.append(MultiIndex::unit(this->argument_size(),j),1);
    return r;
}

TaylorModel<ValidatedFloat> TaylorModel<ValidatedFloat>::create_ball(ErrorType e) const {
    ARIADNE_PRECONDITION(e.raw()>=0);
    TaylorModel<ValidatedFloat> r(this->argument_size(),this->_sweeper);
    r._error=e;
    return r;
}

Void TaylorModel<ValidatedFloat>::swap(TaylorModel<ValidatedFloat>& tm) {
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
    std::swap(this->_sweeper,tm._sweeper);
}

Void TaylorModel<ValidatedFloat>::clear() {
    this->_expansion.clear();
    this->_error=0u;
}

DegreeType TaylorModel<ValidatedFloat>::degree() const {
    DegreeType deg=0u;
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        deg=max(deg,iter->key().degree());
    }
    return deg;
}


TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::operator=(const ValidatedNumber& c) {
    this->_expansion.clear();
    ExactFloat m=c.midpoint();
    if(m!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    this->_error=c.radius();
    return *this;
}


namespace { // Internal code for arithmetic

ExactFloat add_err(ExactFloat const& x1, ExactFloat const& x2, ErrorFloat& e) {
    ExactFloat mx1=-x1;
    set_rounding_to_nearest();
    ExactFloat r(x1.raw() + x2.raw());
    set_rounding_upward();
    Float u=x1.raw()+x2.raw();
    Float ml=mx1.raw()-x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

ExactFloat add_err(ExactFloat const& x, ValidatedFloat const& c, ErrorFloat& e) {
    Float const& xv=x.raw();
    Float const cl=c.lower().raw();
    Float const cu=c.upper().raw();
    Float& re=e.raw();
    set_rounding_to_nearest();
    Float cm=(cl+cu)/2;
    Float rv=xv+cm;
     set_rounding_upward();
    Float u=xv+cu;
    Float ml=(-xv)-cl;
    re += (u+ml)/2;
    return ExactFloat(rv);
}

ExactFloat sub_err(ExactFloat const& x1, ExactFloat const& x2, ErrorFloat& e) {
    ExactFloat mx1=-x1;
    set_rounding_to_nearest();
    ExactFloat r(x1.raw() - x2.raw());
    set_rounding_upward();
    Float u=x1.raw()-x2.raw();
    Float ml=mx1.raw()+x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

ExactFloat mul_no_err(ExactFloat const& x1, ExactFloat const& x2) {
    set_rounding_to_nearest();
    ExactFloat r(x1.raw() * x2.raw());
    set_rounding_upward();
    return r;
}

ExactFloat mul_err(ExactFloat const& x1, ExactFloat const& x2, ErrorFloat& e) {
    ExactFloat mx1=-x1;
    set_rounding_to_nearest();
    ExactFloat r(x1.raw() * x2.raw());
    set_rounding_upward();
    Float u=x1.raw()*x2.raw();
    Float ml=mx1.raw()*x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

ExactFloat mul_err(ExactFloat const& x, ValidatedFloat const& vc, ApproximateFloat ca, ErrorFloat& e) {
    Float const& xv=x.raw();
    Float const cu=vc.upper().raw();
    Float const cl=vc.lower().raw();
    Float& re=e.raw();
    set_rounding_to_nearest();
    Float rv=xv*ca.raw();
    set_rounding_upward();
    Float u,ml;
    if(xv>=0) {
        Float mcl=-cl;
        u=xv*cu;
        ml=xv*mcl;
    } else {
        Float mcu=-cu;
        u=xv*cl;
        ml=xv*mcu;
    }
    re+=(u+ml)/2;
    return ExactFloat(rv);
}

ExactFloat div_err(ExactFloat const& x1, ExactFloat const& x2, ErrorFloat& e) {
    ExactFloat mx1=-x1;
    set_rounding_to_nearest();
    ExactFloat r(x1.raw() / x2.raw());
    set_rounding_upward();
    Float u=x1.raw()/x2.raw();
    Float ml=mx1.raw()/x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

// Inplace negation
Void _neg(TaylorModel<ValidatedFloat>& r)
{
    for(auto iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data().raw()=-iter->data().raw();
    }
}

inline Void _scal_exact(TaylorModel<ValidatedFloat>& r, const RawFloat& c)
{
    // Operation can be performed exactly
    for(TaylorModel<ValidatedFloat>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data().raw()*=c;
    }
    r.error().raw()*=abs(c);
    return;
}

#ifndef ARIADNE_AVOID_ROUNDING_MODE_CHANGES

inline Void _scal_approx(TaylorModel<ValidatedFloat>& r, const RawFloat& c) {
    // General case with error analysis
    // Perform one loop, with rounding mode changes each time
    RawFloat& re=r.error().raw();
    RawFloat u,ml;
    RawFloat te=0; // Twice the maximum accumulated error
    RawFloat pc=c;
    RawFloat mc=-c;    set_rounding_upward(); // Should be unnecessary; default rounding mode is upward
    for(TaylorModel<ValidatedFloat>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        RawFloat& rv=riter->data().raw();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
        set_rounding_to_nearest();
        rv*=c;
        set_rounding_upward();
    }

    re*=abs(c);
    re+=te/2;

    return;
}

#else // ARIADNE_AVOID_ROUNDING_MODE_CHANGES

inline Void _scal_approx(TaylorModel<ValidatedFloat>& r, const RawFloat& c) {
    // General case with error analysis
    // Perform two loops to avoid unnecessary rounding mode changes.
    RawFloat& re=r.error().raw();
    RawFloat u,ml;
    RawFloat te=0; // Twice the maximum accumulated error
    RawFloat pc=c;
    RawFloat mc=-c;
    set_rounding_upward();
    for(TaylorModel<ValidatedFloat>::ConstIterator riter=r.begin(); riter!=r.end(); ++riter) {
        const RawFloat& rv=riter->data().raw();
        u=rv*pc;
        ml=rv*mc;
        te+=(u+ml);
    }
    re*=abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    RawFloat m=c;
    for(TaylorModel<ValidatedFloat>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data().raw()*=m;
    }
    return;
}
#endif // ARIADNE_AVOID_ROUNDING_MODE_CHANGES



Void _scal(TaylorModel<ValidatedFloat>& r, const RawFloat& c)
{
    // No measurable speedup in general case by avoiding checks
    if(c==0.0) { r.expansion().clear(); r.error().raw()=0; }
    else if(c==1.0) { }
    else if(c==2.0 || c==0.5) { _scal_exact(r,c); }
    else { _scal_approx(r,c); }
}


inline Void _scal(TaylorModel<ValidatedFloat>& r, const ExactFloat& c) {
    ErrorFloat e=0u; // The maximum accumulated error
    for(TaylorModel<ValidatedFloat>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data() = mul_err(riter->data(),c,e);
    }

    ErrorFloat& re=r.error();
    re*=abs(c);
    re+=e;
}

Void _scal(TaylorModel<ValidatedFloat>& r, const ValidatedFloat& c)
{
    //std::cerr<<"TaylorModel<ValidatedFloat>::scal(ValidatedFloat c) c="<<c<<std::endl;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);

    if(r.error().raw()==inf) {
        r.expansion().clear(); return;
    }

    ErrorFloat e=0u;
    ApproximateFloat ca=c.value();
    for(TaylorModel<ValidatedFloat>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        ExactFloat& rv=riter->data();
        rv=mul_err(rv,c,ca,e);
    }
    r.error()*=mag(c);
    r.error()+=e;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


struct UnitMultiIndex { SizeType argument_size; SizeType unit_index; };


inline Void _incr(TaylorModel<ValidatedFloat>& r, const MultiIndex& a) {
    for(TaylorModel<ValidatedFloat>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        static_cast<MultiIndex&>(iter->key())+=a;
    }
}

inline Void _incr(TaylorModel<ValidatedFloat>& r, SizeType j) {
    for(TaylorModel<ValidatedFloat>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++static_cast<MultiIndex&>(iter->key())[j];
    }
}


inline Void _acc(TaylorModel<ValidatedFloat>& r, const ExactFloat& c) {
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    if(c==0) { return; }
    if(r.expansion().empty()) {
        r.expansion().raw().append(MultiIndex(r.argument_size()),c.raw());
    } else if((r.end()-1)->key().degree()>0) {
        r.expansion().raw().append(MultiIndex(r.argument_size()),c.raw());
    } else {
        ExactFloat& rv=(r.end()-1)->data();
        ErrorFloat& re=r.error();
        rv=add_err(rv,c,re);
    }
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    return;
}



inline Void _acc(TaylorModel<ValidatedFloat>& r, const ValidatedFloat& c)
{
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);

    if(c.lower().raw()==-inf || c.upper().raw()==+inf) {
        r.clear();
        r.set_error(+infty);
        return;
    }

    if(c.lower().raw()==-c.upper().raw()) { // The midpoint of the interval is zero, so no need to change constant term
        r.error()+=mag(c);
        return;
    }

    if(r.expansion().empty()) { // Append a constant term zero
        r._append(MultiIndex(r.argument_size()),0);
    } else if((r.end()-1)->key().degree()>0) { // Append a constant term zero
        r._append(MultiIndex(r.argument_size()),0);
    }

    ExactFloat& rv=(r.end()-1)->data();
    ErrorFloat& re=r.error();
    rv=add_err(rv,c,re);

    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}


#ifndef ARIADNE_AVOID_ROUNDING_MODE_CHANGES
// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
inline Void _add(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y)
{
    ARIADNE_PRECONDITION(r.number_of_nonzeros()==0);
    ErrorFloat e=0u;
    TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin();
    TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r._append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r._append(yiter->key(),yiter->data());
            ++yiter;
        } else {
            ARIADNE_DEBUG_ASSERT(xiter->key()==yiter->key());
            r._append(xiter->key(),add_err(xiter->data(),yiter->data(),e));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r._append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r._append(yiter->key(),yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=(x.error()+y.error())+e;
    set_rounding_to_nearest();

    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}

#else // ARIADNE_AVOID_ROUNDING_MODE_CHANGES

// Compute r=x+y, assuming r is empty
// First compute the errors and then at the end compute the expansion
// This saves rounding mode changes, but requires an extra loop
inline Void _add(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y)
{
    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    Float te=0.0;
    TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin();
    TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            ++yiter;
        } else {
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            VOLATILE Float u=xv+yv;
            VOLATILE Float t=-xv;
            VOLATILE Float ml=t-yv;
            te+=(u+ml);
            ++xiter; ++yiter;
        }
    }
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);

    set_rounding_to_nearest();
    xiter=x.begin();
    yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r.expansion().raw().append(yiter->key(),yiter->data().raw());
            ++yiter;
        } else {
            RawFloat c=xiter->data().raw()+yiter->data().raw();
            if(c!=0) { r.expansion().raw().append(xiter->key(),c); }
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().raw().append(yiter->key(),yiter->data().raw());
        ++yiter;
    }
    ARIADNE_DEBUG_ASSERT_MSG(r.error()>=0,r);

}

#endif // ARIADNE_AVOID_ROUNDING_MODE_CHANGES

inline Void _acc(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x)
{
    TaylorModel<ValidatedFloat> s(r.argument_size(),r.sweeper()); _add(s,r,x); s.swap(r);
    ARIADNE_DEBUG_ASSERT_MSG(r.error()>=0,r);
}

inline Void _sub(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y)
{
    ARIADNE_PRECONDITION(r.number_of_nonzeros()==0);
    ErrorFloat e=0u;
    TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin();
    TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r._append(xiter->key(),xiter->data());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            r._append(yiter->key(),-yiter->data());
            ++yiter;
        } else {
            ARIADNE_DEBUG_ASSERT(xiter->key()==yiter->key());
            r._append(xiter->key(),sub_err(xiter->data(),yiter->data(),e));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r._append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r._append(yiter->key(),-yiter->data());
        ++yiter;
    }

    set_rounding_upward();
    r.error()=(x.error()+y.error())+e;
    set_rounding_to_nearest();

    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


inline Void _sma(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x, const ValidatedFloat& c, const TaylorModel<ValidatedFloat>& y)
{
    ARIADNE_ASSERT_MSG(c.lower().raw()<=c.upper().raw(),c);
    ARIADNE_ASSERT_MSG(x.error().raw()>=0,"x="<<x);
    ARIADNE_ASSERT_MSG(y.error().raw()>=0,"y="<<y);

    //std::cerr<<"TaylorModel<ValidatedFloat>::scal(ValidatedNumber c) c="<<c<<std::endl;
    set_rounding_upward();
    VOLATILE Float u,ml,myv;
    Float te=0; // Twice the maximum accumulated error
    Float cu=c.upper().raw();
    Float mcl=-c.lower().raw();
    Float cm=c.midpoint().raw();

    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin();
    TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            r.expansion().raw().append(xiter->key(),xiter->data().raw());
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            set_rounding_upward();
            const Float& yv=yiter->data().raw();
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
            r.expansion().raw().append(yiter->key(),cm*yiter->data().raw());
            ++yiter;
        } else {
            set_rounding_upward();
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
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
            r.expansion().raw().append(xiter->key(),xiter->data().raw()+cm*yiter->data().raw());
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r.expansion().raw().append(xiter->key(),xiter->data().raw());
        ++xiter;
    }
    while(yiter!=y.end()) {
        set_rounding_upward();
        const Float& yv=yiter->data().raw();
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
        r.expansion().raw().append(yiter->key(),cm*yiter->data().raw());
        ++yiter;
    }

    set_rounding_upward();
    r.error().raw()=x.error().raw();
    r.error().raw()+=y.error().raw();
    r.error().raw()+=(te/2);
    set_rounding_to_nearest();

    ARIADNE_DEBUG_ASSERT_MSG(r.error()>=0,r);
}



#ifndef ARIADNE_AVOID_ROUNDING_MODE_CHANGES

// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
Void _mul(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y)
{
    const SizeType as=r.argument_size();
    TaylorModel<ValidatedFloat> t(as,r.sweeper());
    TaylorModel<ValidatedFloat> s(as,r.sweeper());
    MultiIndex ta(as);
    for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        ErrorFloat te=0u; // trucation error
        ErrorFloat re=0u; // roundoff error
        const MultiIndex& xa=xiter->key();
        const ExactFloat& xv=xiter->data();
        ErrorFloat mag_xv=mag(xv);
        for(TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const MultiIndex& ya=yiter->key();
            const ExactFloat& yv=yiter->data();
            ta=xa+ya;
            ExactFloat tv=mul_no_err(xv,yv);
            if(r.sweeper().discard(ta,tv.raw())) {
                te+=mag_xv*mag(yv);
            } else {
                t._append(ta,tv);
                re+=(xv*yv).error();
            }
        }
        t.error()=te+re;

        _add(s,r,t);
        r.expansion().swap(s.expansion());
        r.error()=s.error();
        s.expansion().clear();
        s.error()=0u;
        t.expansion().clear();
        t.error()=0u;
    }

    ErrorFloat xs=0u;
    for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=mag(xiter->data());
    }

    ErrorFloat ys=0u;
    for(TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=mag(yiter->data());
    }

    ErrorFloat& re=r.error();
    const ErrorFloat& xe=x.error();
    const ErrorFloat& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    return;
}

#else // ARIADNE_AVOID_ROUNDING_MODE_CHANGES

// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
inline
Void _mul(TaylorModel<ValidatedFloat>& r, const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y)
{
    TaylorModel<ValidatedFloat> t(x.argument_size(),r.sweeper());
    TaylorModel<ValidatedFloat> s(x.argument_size(),r.sweeper());
    for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        set_rounding_upward();
        Float te=0.0;
        for(TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->data().raw();
            const Float& yv=yiter->data().raw();
            VOLATILE Float u=xv*yv;
            VOLATILE Float ml=-xv; ml=ml*yv;
            te+=(u+ml);
        }
        t.error().raw()=te/2;
        set_rounding_to_nearest();
        for(TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            t.expansion().raw().append_sum(xiter->key(),yiter->key(),xiter->data().raw()*yiter->data().raw());
        }
        _add(s,r,t);
        r.expansion().swap(s.expansion());
        r.error().raw()=s.error().raw();
        s.expansion().clear();
        s.error().raw()=0.0;
        t.expansion().clear();
        t.error().raw()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->data().raw());
    }

    Float ys=0;
    for(TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->data().raw());
    }

    Float& re=r.error().raw();
    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();
    re+=xs*ye+ys*xe+xe*ye;

    r.sweep();

    set_rounding_to_nearest();
    return;
}

#endif

} // namespace


///////////////////////////////////////////////////////////////////////////////

// Inplace arithmetical operations for Algebra concept

Void TaylorModel<ValidatedFloat>::iadd(const ValidatedNumber& c)
{
    _acc(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}

Void TaylorModel<ValidatedFloat>::imul(const ValidatedNumber& c)
{
    _scal(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}

Void TaylorModel<ValidatedFloat>::isma(const ValidatedNumber& c, const TaylorModel<ValidatedFloat>& y)
{
    TaylorModel<ValidatedFloat>& x=*this;
    TaylorModel<ValidatedFloat> r=this->create();
    _sma(r,x,c,y);
    this->swap(r);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}

Void TaylorModel<ValidatedFloat>::ifma(const TaylorModel<ValidatedFloat>& x1, const TaylorModel<ValidatedFloat>& x2)
{
    _mul(*this,x1,x2);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}



///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::sort() {
    this->_expansion.reverse_lexicographic_sort();
};

TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::unique()
{
    TaylorModel<ValidatedFloat>::ConstIterator advanced =this->begin();
    TaylorModel<ValidatedFloat>::ConstIterator end =this->end();
    TaylorModel<ValidatedFloat>::Iterator current=this->begin();
    ErrorFloat e=0u;
    while(advanced!=end) {
        current->key()=advanced->key();
        ExactFloat rv=advanced->data();
        ++advanced;
        while(advanced!=end && advanced->key()==current->key()) {
            const ExactFloat& xv=advanced->data();
            rv=add_err(rv,xv,e);
        }
        current->data()=rv;
        ++current;
    }
    this->error()+=e;
    this->_expansion.resize(current-this->begin());

    return *this;
}

TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::sweep() {
    this->_sweeper.sweep(this->_expansion.raw(),this->_error.raw());
    return *this;
}

TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::sweep(const Sweeper& sweeper) {
    sweeper.sweep(this->_expansion.raw(),this->_error.raw());
    return *this;
}

TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::cleanup() {
    this->sort();
    this->unique();
    this->sweep();
    return *this;
}

TaylorModel<ValidatedFloat>& TaylorModel<ValidatedFloat>::clobber() {
    this->_error=0u;
    return *this;
}



///////////////////////////////////////////////////////////////////////////////

// Accuracy control

Float TaylorModel<ValidatedFloat>::tolerance() const {
    const ThresholdSweeper* ptr=dynamic_cast<const ThresholdSweeper*>(&static_cast<const SweeperInterface&>(this->_sweeper));
    if(ptr) {
        return ptr->sweep_threshold();
    } else {
        return std::numeric_limits<double>::epsilon();
    }
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

Box<UnitInterval>
TaylorModel<ValidatedFloat>::domain() const
{
    return Box<UnitInterval>(this->argument_size(),UnitInterval());
}

ExactInterval
TaylorModel<ValidatedFloat>::codomain() const
{
    return make_exact_interval(this->range());
}


UpperInterval _range1(const TaylorModel<ValidatedFloat>& tm) {
    set_rounding_mode(upward);
    ExactFloat v=0;
    ErrorFloat e=0u;
    for(auto iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            v=iter->data();
        } else {
            e+=abs(iter->data());
        }
    }
    UpperInterval r(v-e,v+e);
    set_rounding_mode(to_nearest);
    return r;
}


// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
UpperInterval _range2(const TaylorModel<ValidatedFloat>& tm) {
    const SizeType as=tm.argument_size();
    Array<Float> linear_terms(as,0.0);
    Array<Float> quadratic_terms(as,0.0);
    UpperInterval r(-tm.error().raw(),+tm.error().raw());
    for(auto iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            r+=iter->data().raw();
        } else if(iter->key().degree()==1) {
            for(SizeType j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==1) { linear_terms[j]=iter->data().raw(); break; }
            }
        } else if(iter->key().degree()==2) {
            for(SizeType j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==2) { quadratic_terms[j]=iter->data().raw(); break; }
                if(iter->key()[j]==1) { r+=abs(iter->data().raw())*UpperInterval(-1,1); break; }
            }
        } else {
            r+=abs(iter->data().raw())*UpperInterval(-1,1);
        }
    }
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(SizeType j=0; j!=as; ++j) {
        const Float& a=quadratic_terms[j];
        const Float& b=linear_terms[j];
        UpperInterval ql=abs(a)*UpperInterval(-1,1) + abs(b)*UpperInterval(-1,+1);
        UpperInterval qf=a*(sqr(UpperInterval(-1,+1)+UpperInterval(b)/(2*a)))-sqr(UpperInterval(b))/(4*a);
        r += intersection(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
    }
    return r;
}



UpperInterval
TaylorModel<ValidatedFloat>::range() const {
    return Ariadne::_range2(*this);
}

UpperInterval TaylorModel<ValidatedFloat>::gradient_range(SizeType j) const {
    SizeType as=this->argument_size();
    UpperInterval g(0,0);
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex const& a=iter->key();
        const Nat c=a[j];
        if(c>0) {
            const ExactFloat& x=iter->data();
            if(a.degree()==1) { g+=UpperInterval(x); }
            else { g+=UpperInterval(-1,1)*x*c; }
        }
    }
    return g;
}



//////////////////////////////////////////////////////////////////////////////

// Exact functions (max, min, abs, neg) and arithmetical functions (sqr, pow)


TaylorModel<ValidatedFloat> max(const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y) {
    UpperInterval xr=x.range();
    UpperInterval yr=y.range();
    if(definitely(xr.lower()>=yr.upper())) {
        return x;
    } else if(definitely(yr.lower()>=xr.upper())) {
        return y;
    } else {
        return ((x+y)+abs(x-y))/ExactFloat(2);
    }
}


TaylorModel<ValidatedFloat> min(const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y) {
    UpperInterval xr=x.range();
    UpperInterval yr=y.range();
    if(definitely(xr.upper()<=yr.lower())) {
        return x;
    } else if(definitely(yr.upper()<=xr.lower())) {
        return y;
    } else {
        return ((x+y)-abs(x-y))/ExactFloat(2);
    }
}

TaylorModel<ValidatedFloat> abs(const TaylorModel<ValidatedFloat>& x) {
    UpperInterval xr=x.range();
    if(definitely(xr.lower()>=0)) {
        return x;
    } else if(definitely(xr.upper()<=0)) {
        return -x;
    } else {
        // Use power series expansion $abs(x)=\sum_{i=0}^{7} p_i x^{2i} \pm e$ for $x\in[-1,+1]$ with
        // p=[0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667] and e=0.035
        // TODO: Find more accurate and stable formula
        static const Nat n=7u;
        static const Dbl p[n]={0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667};
        static const Dbl err=0.035;
        TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());
        ExactFloat xmag=make_exact(mag(xr));
        TaylorModel<ValidatedFloat> s=x/xmag;
        s=sqr(s);
        r=static_cast<ExactFloat>(p[n-1]);
        for(Nat i=0; i!=(n-1); ++i) {
            Nat j=(n-2)-i;
            r=s*r+p[j];
        }
        r+=ValidatedFloat(-err,+err);
        return r*xmag;
    }
}

TaylorModel<ValidatedFloat> neg(const TaylorModel<ValidatedFloat>& x) {
    TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());
    for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r._append(xiter->key(),-xiter->data());
    }
    r.error()=x.error();
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

TaylorModel<ValidatedFloat> sqr(const TaylorModel<ValidatedFloat>& x) {
    return x*x;
}

TaylorModel<ValidatedFloat> pow(const TaylorModel<ValidatedFloat>& x, Int n) {
    TaylorModel<ValidatedFloat> r=x.create_constant(1);
    TaylorModel<ValidatedFloat> p(x);
    while(n) { if(n%2) { r=r*p; } p=sqr(p); n/=2; }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
template<class X> class TaylorSeries;


TaylorModel<ValidatedFloat>
compose(const TaylorSeries<ValidatedFloat>& ts, const TaylorModel<ValidatedFloat>& tv)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(HALF_MACHINE_EPSILON));
    ExactFloat& vref=const_cast<ExactFloat&>(tv.value());
    ExactFloat vtmp=vref;
    vref=0;
    TaylorModel<ValidatedFloat> r(tv.argument_size(),tv.sweeper());
    r+=ts[ts.degree()];
    for(SizeType i=1; i<=ts.degree(); ++i) {
        r=r*tv;
        r+=ts[ts.degree()-i];
        r.sweep(threshold_sweeper);
    }
    r.error()+=ts.error();
    vref=vtmp;
    r.sweep();
    return r;
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
TaylorModel<ValidatedFloat>
_compose1(const AnalyticFunction& fn, const TaylorModel<ValidatedFloat>& tm, RawFloat eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    static const DegreeType DEGREE=18;
    static const RawFloat TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    ExactFloat c=tm.value();
    UpperInterval r=tm.range();
    Series<ValidatedFloat> centre_series=fn.series(c);
    Series<ValidatedFloat> range_series=fn.series(static_cast<ValidatedFloat>(r));

    ErrorFloat truncation_error_estimate=(mag(range_series[d])*pow(mag(r-c),d));
    if(truncation_error_estimate.raw()>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error_estimate
                     <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    TaylorModel<ValidatedFloat> x=tm-c;
    TaylorModel<ValidatedFloat> res(tm.argument_size(),tm.sweeper());
    res+=range_series[d];
    for(DegreeType i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(threshold_sweeper);
    }
    return res;
}

// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. The radius of convergence of this method is still quite low,
// typically only half of the radius of convergence of the power series itself
TaylorModel<ValidatedFloat>
_compose2(const AnalyticFunction& fn, const TaylorModel<ValidatedFloat>& tm, RawFloat eps)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(eps));
    static const DegreeType DEGREE=20;
    static const RawFloat TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    ExactFloat c=tm.value();
    UpperInterval r=tm.range();
    Series<ValidatedFloat> centre_series=fn.series(c);
    Series<ValidatedFloat> range_series=fn.series(static_cast<ValidatedFloat>(r));

    ErrorFloat truncation_error=mag(range_series[d]-centre_series[d])*pow(mag(r-c),d);
    if(truncation_error.raw()>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    TaylorModel<ValidatedFloat> x=tm-c;
    TaylorModel<ValidatedFloat> res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(threshold_sweeper);
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
TaylorModel<ValidatedFloat>
_compose3(const AnalyticFunction& fn, const TaylorModel<ValidatedFloat>& tm, RawFloat eps)
{
    static const DegreeType DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    Nat d=DEGREE;
    ExactFloat c=tm.value();
    UpperInterval r=tm.range();
    Series<ValidatedFloat> centre_series=fn.series(c);
    Series<ValidatedFloat> range_series=fn.series(static_cast<ValidatedFloat>(r));

    ValidatedFloat se=range_series[d]-centre_series[d];
    ValidatedFloat e=static_cast<ValidatedFloat>(r)-c;
    ValidatedFloat p=pow(e,d-1);
    p=ValidatedFloat(-p.lower().raw()*e.lower().raw(),p.upper().raw()*e.upper().raw());
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float truncation_error=max(se.lower().raw()*p.lower().raw(),se.upper().raw()*p.upper().raw());
    if(truncation_error>TRUNCATION_ERROR) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n");
    }

    TaylorModel<ValidatedFloat> x=tm;
    TaylorModel<ValidatedFloat> res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        // Don't sweep here...
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}


TaylorModel<ValidatedFloat>
compose(const AnalyticFunction& fn, const TaylorModel<ValidatedFloat>& tm, RawFloat eps) {
    return _compose3(fn,tm,eps);
}

TaylorModel<ValidatedFloat>
compose(const AnalyticFunction& fn, const TaylorModel<ValidatedFloat>& tm) {
    return compose(fn,tm,HALF_MACHINE_EPSILON);
}




///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions
//   bounded domain (rec,sqrt,log,tan)
//   unbounded domain (exp,sin,cos)

static const DegreeType MAXIMUM_DEGREE = 18;


TaylorModel<ValidatedFloat> sqrt(const TaylorModel<ValidatedFloat>& x) {
    //std::cerr<<"rec(TaylorModel<ValidatedFloat>)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    UpperInterval r=x.range();

    if(r.lower().raw()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x.range());
    }

    assert(r.lower().raw()>0);
    Float a=(r.lower().raw()+r.upper().raw())/2;
    set_rounding_upward();
    Float eps=(r.upper().raw()-r.lower().raw())/(r.upper().raw()+r.lower().raw());
    set_rounding_to_nearest();
    assert(eps<1);
    Nat d=integer_cast<Int>((log((1-eps)*x.tolerance())/log(eps)+1));
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    TaylorModel<ValidatedFloat> y=(x/a)-1.0;
    //std::cerr<<"y="<<y<<std::endl;
    TaylorModel<ValidatedFloat> z(x.argument_size(),x.sweeper());
    Series<ValidatedFloat> sqrt_series=Series<ValidatedFloat>::sqrt(ValidatedFloat(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(Nat i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        z.sweep();
        //std::cerr<<"z="<<z<<std::endl;
    }
    Float trunc_err=(pow(eps,d)/(1-eps))*mag(sqrt_series[d]).raw();
    //std::cerr<<"te="<<trunc_err<<" te*[-1,+1]="<<trunc_err*ValidatedFloat(-1,1)<<std::endl;
    z.error().raw()+=trunc_err;
    //std::cerr<<"z="<<z<<std::endl;
    ValidatedFloat sqrta=sqrt(ValidatedFloat(a));
    //std::cerr<<"sqrt(a)="<<sqrta<<std::endl;
    z*=sqrt(ValidatedFloat(a));
    //std::cerr<<"z="<<z<<std::endl;
    return z;
}

TaylorModel<ValidatedFloat> rec(const TaylorModel<ValidatedFloat>& x) {
    //std::cerr<<"rec(TaylorModel<ValidatedFloat>)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    ValidatedFloat r=static_cast<ValidatedFloat>(x.range());
    if(r.upper().raw()>=0 && r.lower().raw()<=0) {
        ARIADNE_THROW(DivideByZeroException,"rec(TaylorModel<ValidatedFloat> x)","x="<<x<<", x.range()="<<x.range());
    }
    Float a=(r.lower().raw()+r.upper().raw())/2;
    set_rounding_upward();
    Float eps=abs((r.upper().raw()-r.lower().raw())/(r.upper().raw()+r.lower().raw()));
    set_rounding_to_nearest();
    assert(eps<1);

    Nat d=integer_cast<Nat>((log((1-eps)*x.tolerance())/log(eps))+1);

    TaylorModel<ValidatedFloat> y=1-(x/a);
    TaylorModel<ValidatedFloat> z(x.argument_size(),x.sweeper());
    z+=Float(d%2?-1:+1);
    //std::cerr<<"  y="<<y<<"\n";
    //std::cerr<<"  z="<<z<<"\n";
    for(Nat i=0; i!=d; ++i) {
        z=1.0 + z * y;
        //std::cerr<<"  z="<<z<<"\n";
    }

    // Compute the truncation error
    Float te=pow(eps,d)/(1-eps);
    //std::cerr<<"  te="<<te<<"\n";
    set_rounding_upward();
    Float nze=te+z.error().raw();
    //std::cerr<<"  nze="<<nze<<"\n";
    set_rounding_to_nearest();
    z._set_error(nze);
    //std::cerr<<"  z="<<z<<"\n";
    z/=a;
    //std::cerr<<"  z="<<z<<"\n";
    return z;
}

TaylorModel<ValidatedFloat> log(const TaylorModel<ValidatedFloat>& x) {
    // Use a special routine to minimise errors
    // Given range [rl,ru], rescale by constant a such that rl/a=1-d; ru/a=1+d
    UpperInterval r=x.range();
    if(r.lower().raw()<=0) {
        ARIADNE_THROW(DomainException,"sqrt",x.range());
    }
    Float a=(r.lower().raw()+r.upper().raw())/2;
    set_rounding_upward();
    Float eps=(r.upper().raw()-r.lower().raw())/(r.upper().raw()+r.lower().raw());
    set_rounding_to_nearest();
    assert(eps<1);
    Nat d=integer_cast<Nat>((log((1-eps)*x.tolerance())/log(eps)+1));
    TaylorModel<ValidatedFloat> y=x/a-1;
    TaylorModel<ValidatedFloat> z(x.argument_size(),x.sweeper());
    z+=Float(d%2?-1:+1)/d;
    for(Nat i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
        z.sweep();
    }
    z=z*y;
    z.sweep();
    Float trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(ValidatedFloat(a))+ValidatedFloat(-trunc_err,+trunc_err);
}

// Use special code to utilise exp(ax+b)=exp(x)^a*exp(b)
TaylorModel<ValidatedFloat> exp(const TaylorModel<ValidatedFloat>& x) {
    // FIXME: Truncation error may be incorrect

    // Scale to unit interval
    TaylorModel<ValidatedFloat> y=x;
    Float xval=x.value().raw();
    y-=xval;
    Float xrad=mag(y.range()).raw();
    Nat sfp=0; // A number such that 2^sfp>rad(x.range())
    while(Float(1<<sfp)<xrad) { ++sfp; }
    Float sf=1.0/(1<<sfp);
    _scal_exact(y,sf);
    Float yrad=xrad*sf;

    const DegreeType degree=MAXIMUM_DEGREE;

    // Since x is in unit domain, truncation error is no worse than maximum omitted term, i.e. xr/fac(d+1)
    TaylorModel<ValidatedFloat> res(x.argument_size(),x.sweeper());
    res._set_error(pow_up(yrad,degree+1));
    for(Nat i=0; i!=degree; ++i) {
        res/=(degree-i);
        res=y*res+1.0;
    }

    // Square r a total of sfp times
    TaylorModel<ValidatedFloat> square(x.argument_size(),x.sweeper());
    for(Nat i=0; i!=sfp; ++i) {
        _mul(square,res,res);
        res.swap(square);
        square.clear();
    }

    // Multiply by exp(xv)
    res*=Ariadne::exp(ValidatedFloat(xval));

    return res;
    //return _compose(&Series<ValidatedFloat>::exp,x,x.tolerance());
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form x*f(x^2)
TaylorModel<ValidatedFloat> sin(const TaylorModel<ValidatedFloat>& x) {
    // FIXME: Truncation error may be incorrect
    TaylorModel<ValidatedFloat> y(x.argument_size(),x.sweeper());
    TaylorModel<ValidatedFloat> s(x.argument_size(),x.sweeper());
    TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());
    TaylorModel<ValidatedFloat> t(x.argument_size(),x.sweeper());

    ExactNumber one(1.0);

    Float two_pi_approx=2*pi_approx;
    Int n=integer_cast<Int>(floor(x.value().raw()/two_pi_approx + 0.5));
    y=x-(n*2*pi_val);

    if(y.error().raw()>two_pi_approx/2) {
        r.error().raw()=one.raw();
    } else {
        _mul(s,y,y);

        Int d=(MAXIMUM_DEGREE+3)/2;
        Float srad=mag(s.range()).raw();
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute x(1-y/6+y^2/120-y^3/5040+... = x(1-y/6*(1-y/20*(1-y/42*...)
        r=one;
        for(Int i=0; i!=d; ++i) {
            r/=Float(-2*(d-i)*(2*(d-i)+1));
            _mul(t,s,r); r.swap(t); t.clear();
            r+=one;
        }
        _mul(t,y,r); r.swap(t);

        r.error().raw()+=truncation_error;
    }

    return r;
}

// Use special code to utilise sin(x+2pi)=sin(x)
// and that the power series is of the form f(x^2)
TaylorModel<ValidatedFloat> cos(const TaylorModel<ValidatedFloat>& x) {
    // FIXME: Truncation error may be incorrect
    TaylorModel<ValidatedFloat> y(x.argument_size(),x.sweeper());
    TaylorModel<ValidatedFloat> s(x.argument_size(),x.sweeper());
    TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());
    TaylorModel<ValidatedFloat> t(x.argument_size(),x.sweeper());

    Float two_pi=2*pi_approx;
    Int n=integer_cast<Int>(floor(x.value().raw()/two_pi + 0.5));

    y=x-2*n*pi_val;
    ExactNumber one(1.0);

    if(y.error().raw()>two_pi/2) {
        r.error().raw()=one.raw();
    } else {
        _mul(s,y,y);

        Int d=(MAXIMUM_DEGREE+3)/2;
        Float srad=mag(s.range()).raw();
        Float truncation_error=pow_up(srad,d+1)*rec_fac_up((d+1)*2);

        // Compute 1-y/2+y^2/24-y^3/720+... = (1-y/2*(1-y/12*(1-y/30*...)
        r=one;
        for(Int i=0; i!=d; ++i) {
            r/=double(-2*(d-i)*(2*(d-i)-1));
            _mul(t,s,r); r.swap(t); t.clear();
            r+=one;
        }

        r.error().raw()+=truncation_error;
    }

    return r;
}

TaylorModel<ValidatedFloat> tan(const TaylorModel<ValidatedFloat>& x) {
    return sin(x)*rec(cos(x));
}

TaylorModel<ValidatedFloat> atan(const TaylorModel<ValidatedFloat>& x) {
    static const DegreeType degree=MAXIMUM_DEGREE;
    TaylorSeries<ValidatedFloat> series(Atan(),x.codomain(),x.value(),degree);
    return compose(series,x);
}




///////////////////////////////////////////////////////////////////////////////

// Inplace function operators (rescale, restrict, antidifferentiate)





// Replace the kth variable x[k] by a*x[k]+b.
TaylorModel<ValidatedFloat> preaffine(const TaylorModel<ValidatedFloat>& tm, SizeType k, const ValidatedNumber& a, const ValidatedNumber& b) {
    DegreeType d=tm.degree();
    SizeType as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    TaylorModel<ValidatedFloat> r(as,swp);
    r._set_error(tm.error().raw());

    // Create a temporary TaylorModels containing just terms x[k]^i
    Array<TaylorModel<ValidatedFloat>> atm(d+1,TaylorModel<ValidatedFloat>(as,swp));
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        MultiIndex a=iter->key();
        const Float& c=iter->data().raw();
        DegreeType ak=a[k];
        a[k]=0;
        atm[ak].expansion().raw().append(a,c);
    }

    TaylorModel<ValidatedFloat> xk=TaylorModel<ValidatedFloat>::coordinate(as,k,swp);

    for(DegreeType i=0; i<=d; ++i) {
        for(DegreeType j=i; j<=d; ++j) {
            ValidatedNumber c=(Ariadne::bin(j,i)*Ariadne::pow(a,i)*Ariadne::pow(b,j-i));
            r+=c*atm[j];
            atm[j]*=xk;
        }
    }
    return r;
}

// Given an array of ordered indices below some maximum, make an array of the indices not in the array
Array<SizeType> complement(SizeType nmax, Array<SizeType> vars) {
    Array<SizeType> cmpl(nmax-vars.size());
    SizeType kr=0; SizeType kv=0;
    for(SizeType j=0; j!=nmax; ++j) {
        if(kv==vars.size() || j!=vars[kv]) {
            cmpl[kr]=j; ++kr;
        } else {
            ++kv;
        }
    }
    return std::move(cmpl);
}


TaylorModel<ValidatedFloat> embed_error(const TaylorModel<ValidatedFloat>& tm) {
    const SizeType as=tm.argument_size();
    TaylorModel<ValidatedFloat> rtm(as+1u,tm.sweeper());
    MultiIndex ra(as+1u);

    // The new error term is first in reverse lexicographic order.
    ExactFloat err_coef=make_exact(tm.error());
    ra[as]=1;
    rtm._append(ra,err_coef);
    ra[as]=0;

    // Copy new terms
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
        MultiIndex const& xa=iter->key();
        ExactFloat const& xv=iter->data();
        for(SizeType j=0; j!=as; ++j) { ra[j]=xa[j]; }
        rtm._append(ra,xv);
    }
    return std::move(rtm);
}

TaylorModel<ValidatedFloat> discard_variables(const TaylorModel<ValidatedFloat>& tm, Array<SizeType> const& discarded_variables) {
    for(SizeType i=0; i!=discarded_variables.size()-1; ++i) {
        ARIADNE_PRECONDITION(discarded_variables[i]<discarded_variables[i+1]);
    }
    ARIADNE_PRECONDITION(discarded_variables[discarded_variables.size()-1]<tm.argument_size());

    const SizeType number_of_variables = tm.argument_size();
    const SizeType number_of_discarded_variables = discarded_variables.size();
    const SizeType number_of_kept_variables = number_of_variables - number_of_discarded_variables;

    // Make an Array of the variables to be kept
    Array<SizeType> kept_variables=complement(number_of_variables,discarded_variables);

    // Construct result and reserve memory
    TaylorModel<ValidatedFloat> rtm(number_of_kept_variables,tm.sweeper());
    rtm.expansion().reserve(tm.number_of_nonzeros()+1u);

    // Set the uniform error of the original model
    // If index_of_error == number_of_error_variables, then the error is kept as a uniform error bound
    MultiIndex ra(number_of_kept_variables);
    ErrorFloat derr=0u; // Magnitude of discarded terms
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        MultiIndex const& xa=iter->key();
        ExactFloat const& xv=iter->data();
        Bool keep=true;
        for(SizeType k=0; k!=number_of_discarded_variables; ++k) {
            if(xa[discarded_variables[k]]!=0) {
                derr += mag(iter->data());
                keep=false;
                break;
            }
        }
        if(keep) {
            for(SizeType k=0; k!=number_of_kept_variables; ++k) {
                ra[k]=xa[kept_variables[k]];
            }
            rtm._append(ra,xv);
        }
    }
    rtm.error()+=derr;

    return rtm;
}



TaylorModel<ValidatedFloat> restrict(const TaylorModel<ValidatedFloat>& tm, SizeType k, const ExactInterval& nd) {
    ARIADNE_ASSERT(k<tm.argument_size());
    ARIADNE_ASSERT(nd.lower().raw()>=-1 && nd.upper().raw()<=+1);
    if(nd.lower().raw()==-1 && nd.upper().raw()==1) {
        return tm;
    } else if(nd.lower().raw()==-1 && nd.upper().raw()==0) {
        return split(tm,k,SplitPart::LOWER);
    } else if(nd.lower().raw()==0 && nd.upper().raw()==1) {
        return split(tm,k,SplitPart::UPPER);
    } else if(nd.lower().raw()==-0.5 && nd.upper().raw()==0.5) {
        return split(tm,k,SplitPart::MIDDLE);
    } else {
        ValidatedNumber a=(nd.upper()-nd.lower())/2;
        ValidatedNumber b=(nd.lower()+nd.upper())/2;
        return preaffine(tm,k,a,b);
    }
}



Void TaylorModel<ValidatedFloat>::antidifferentiate(SizeType k) {
    TaylorModel<ValidatedFloat>& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());

    ErrorFloat e=0;
    for(TaylorModel<ValidatedFloat>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex& xa=xiter->key();
        ExactFloat& xv=xiter->data();
        xa[k]+=1;
        Nat c=xa[k];
        xv=div_err(xv,c,e);
    }

    x.error()+=e;
}

TaylorModel<ValidatedFloat> antiderivative(const TaylorModel<ValidatedFloat>& x, SizeType k) {
    TaylorModel<ValidatedFloat> r(x);
    r.antidifferentiate(k);
    return r;
}


// Compute derivative inplace by computing term-by-term, switching the rounding mode
// Note that since some terms may be eliminated, requiring two iterators.
Void TaylorModel<ValidatedFloat>::differentiate(SizeType k) {
    TaylorModel<ValidatedFloat> const& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());

    // ARIADNE_PRECONDITION_MSG(x.error().raw()==0,x);
    this->clobber();

    TaylorModel<ValidatedFloat>& r=*this;
    ErrorFloat& re=r.error();
    TaylorModel<ValidatedFloat>::Iterator riter=r.begin();
    for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex const& xa=xiter->key();
        ExactFloat const& xv=xiter->data();
        Nat c=xa[k];
        if(c!=0) {
            MultiIndex& ra=riter->key();
            ExactFloat& rv=riter->data();
            ra=xa; ra[k]-=1;
            rv=mul_err(xv,c,re);
            ++riter;
        }
    }

    r.expansion().resize(riter - r.begin());
}




TaylorModel<ValidatedFloat> derivative(const TaylorModel<ValidatedFloat>& x, SizeType k) {
    TaylorModel<ValidatedFloat> rx=x; rx.differentiate(k); return rx;

    ARIADNE_ASSERT(k<x.argument_size());

    MultiIndex ra(x.argument_size()); ExactFloat rv; Nat c;

    TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());
    ErrorFloat& re=r.error();
    for(TaylorModel<ValidatedFloat>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex const& xa=xiter->key();
        ExactFloat const& xv=xiter->data();
        c=xa[k];
        if(c!=0) {
            ra=xa;
            ra[k]-=1;
            rv=mul_err(xv,c,re);
            r._append(ra,rv);
        }
    }

    return r;
}



///////////////////////////////////////////////////////////////////////////////

// Scalar function operators (evaluate, split, unscale, embed, restrict)
// and predicates (refines)

ValidatedNumber evaluate(const TaylorModel<ValidatedFloat>& tm, const Vector<ValidatedNumber>& x) {
    return horner_evaluate(tm.expansion(),x)+ValidatedNumber(-tm.error(),+tm.error());
}

ApproximateNumber evaluate(const TaylorModel<ValidatedFloat>& tm, const Vector<ApproximateNumber>& x) {
    return horner_evaluate(tm.expansion(),x);
}


Covector<ValidatedNumber> gradient(const TaylorModel<ValidatedFloat>& f, const Vector<ValidatedNumber>& x)
{
    Vector< Differential<ValidatedNumber> > dx=Differential<ValidatedNumber>::variables(1u,x);
    Differential<ValidatedNumber> df=horner_evaluate(f.expansion(),dx);
    return gradient(df);
}



TaylorModel<ValidatedFloat> compose(TaylorModel<ValidatedFloat> const& x, Vector<TaylorModel<ValidatedFloat>> const& y) {
    return horner_evaluate(x.expansion(),y)+ValidatedFloat(-x.error(),+x.error());
}



template<class T> class Powers {
  public:
    explicit Powers(const T& t) { _values.push_back(t*0+1); _values.push_back(t); }
    explicit Powers(const T& z, const T& t) { _values.push_back(z); _values.push_back(t); }
    const T& operator[](SizeType i) const { while(_values.size()<=i) { _values.push_back(_values[1]*_values.back()); } return _values[i]; }
  private:
    mutable std::vector<T> _values;
};

template<> class Powers<ValidatedNumber> {
  public:
    explicit Powers(const ValidatedNumber& t) { _values.push_back(ValidatedNumber(1)); _values.push_back(t); }
    explicit Powers(const ValidatedNumber& z, const ValidatedNumber& t) { _values.push_back(z); _values.push_back(t); }
    const ValidatedNumber& operator[](SizeType i) const {
        while(_values.size()<=i) {
            if(_values.size()%2==0) { _values.push_back(sqr(_values[_values.size()/2])); }
            else { _values.push_back(_values[1]*_values.back()); } }
        return _values[i]; }
  private:
    mutable std::vector<ValidatedNumber> _values;
};



TaylorModel<ValidatedFloat>
partial_evaluate(const TaylorModel<ValidatedFloat>& x, SizeType k, ValidatedNumber c)
{
    TaylorModel<ValidatedFloat> r(x.argument_size()-1,x.sweeper());
    MultiIndex ra(r.argument_size());
    if(c==ValidatedNumber(0)) {
        for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            MultiIndex::IndexType xak=xa[k];
            if(xak==0) {
                const Float& xv=xiter->data().raw();
                for(SizeType i=0; i!=k; ++i) { ra[i]=xa[i]; }
                for(SizeType i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
                r.expansion().raw().append(ra,xv);
            }
        }
        r._set_error(x.error().raw());
    } else if(c==ValidatedNumber(1)) {
        TaylorModel<ValidatedFloat> s(x.argument_size()-1,x.sweeper());
        Array<TaylorModel<ValidatedFloat>> p(x.degree()+1,TaylorModel<ValidatedFloat>(x.argument_size()-1,x.sweeper()));

        for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const Float& xv=xiter->data().raw();
            MultiIndex::IndexType xak=xa[k];
            for(SizeType i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(SizeType i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().raw().append(ra,xv);
        }

        r=p[0];
        r._set_error(x.error().raw());
        for(SizeType i=1; i!=p.size(); ++i) {
            _add(s,r,p[i]);
            r.swap(s);
            s.clear();
        }
    } else {
        TaylorModel<ValidatedFloat> s(x.argument_size()-1,x.sweeper());
        Array<TaylorModel<ValidatedFloat>> p(x.degree()+1,TaylorModel<ValidatedFloat>(x.argument_size()-1,x.sweeper()));

        Powers<ValidatedNumber> cpowers(c);

        for(TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
            const MultiIndex& xa=xiter->key();
            const Float& xv=xiter->data().raw();
            MultiIndex::IndexType xak=xa[k];
            for(SizeType i=0; i!=k; ++i) { ra[i]=xa[i]; }
            for(SizeType i=k; i!=ra.size(); ++i) { ra[i]=xa[i+1]; }
            assert(ra.degree()+xak==xa.degree());
            p[xak].expansion().raw().append(ra,xv);
        }
        for(SizeType i=1; i!=p.size(); ++i) {
            p[i]*=cpowers[i];
        }

        r=p[0];
        r._set_error(x.error().raw());
        for(SizeType i=1; i!=p.size(); ++i) {
            _add(s,r,p[i]);
            r.swap(s);
            s.clear();
        }
    }

    return r;
}



TaylorModel<ValidatedFloat> embed(SizeType as1, const TaylorModel<ValidatedFloat>& tm2, SizeType as3) {
    return TaylorModel<ValidatedFloat>(embed(as1,tm2.expansion(),as3),tm2.error(),tm2.sweeper());
}


TaylorModel<ValidatedFloat> split(const TaylorModel<ValidatedFloat>& tm, SizeType k, SplitPart h) {
    const DegreeType deg=tm.degree();
    const SizeType as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    TaylorModel<ValidatedFloat> r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(TaylorModel<ValidatedFloat>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const uchar ak=iter->key()[k];
        Float& c=iter->data().raw();
        c/=(1<<ak);
    }


    if(h==SplitPart::MIDDLE) { return r; }
    Int tr=( h==SplitPart::UPPER ? +1 : -1 );

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<TaylorModel<ValidatedFloat>> ary(deg+1,TaylorModel<ValidatedFloat>(as,swp));
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->key();
        const Float& c=iter->data().raw();
        uchar ak=a[k];
        a[k]=0u;
        ary[ak].expansion().raw().append(a,c);
    }

    Float re=r.error().raw();
    r.clear();
    r._set_error(re);

    for(DegreeType i=0; i<=deg; ++i) {
        for(DegreeType j=i; j<=deg; ++j) {
            Int sf=bin(j,i);
            if(tr==-1 && (j-i)%2==1) { sf=-sf; }
            r+=ary[j]*sf;
            for(TaylorModel<ValidatedFloat>::Iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->key()[k];
            }
         }
    }

    return r;
}

TaylorModel<ValidatedFloat>
unscale(const TaylorModel<ValidatedFloat>& tv, const ExactInterval& ivl) {
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

    auto& l=ivl.lower().raw();
    auto& u=ivl.upper().raw();
    ARIADNE_ASSERT_MSG(l<=u,"Cannot unscale TaylorModel<ValidatedFloat> "<<tv<<" from empty interval "<<ivl);

    if(l==u) {
        TaylorModel<ValidatedFloat> r=tv.create();
        r+=static_cast<ValidatedFloat>(ivl);
        // Uncomment out line below to make unscaling to a singleton interval undefined
        //r.set_error(+inf);
        return r;
    } else {
        TaylorModel<ValidatedFloat> r=tv;
        ValidatedNumber c=ValidatedNumber(l/2)+ValidatedNumber(u/2);
        ValidatedNumber s=2/(ValidatedNumber(u)-ValidatedNumber(l));
        r-=c;
        r*=s;

        return r;
    }
}

TaylorModel<ValidatedFloat> compose(const Unscaling& u, const TaylorModel<ValidatedFloat>& y) {
    return u(y);
}




///////////////////////////////////////////////////////////////////////////////

// Banach algebra operations


TaylorModel<ValidatedFloat>::CoefficientType TaylorModel<ValidatedFloat>::average() const {
    return (*this)[MultiIndex::zero(this->argument_size())];
}

TaylorModel<ValidatedFloat>::NormType TaylorModel<ValidatedFloat>::radius() const {
    ErrorFloat r=0u;
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=mag(iter->data());
        }
    }
    r+=this->error();
    return r;
}

TaylorModel<ValidatedFloat>::NormType TaylorModel<ValidatedFloat>::norm() const {
    ErrorFloat r=0u;
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=mag(iter->data());
    }
    r+=this->error();
    return r;
}

NormType norm(const TaylorModel<ValidatedFloat>& tv) {
    return tv.norm();
}


Bool refines(const TaylorModel<ValidatedFloat>& tv1, const TaylorModel<ValidatedFloat>& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    TaylorModel<ValidatedFloat> d=tv2;
    d.error().raw()=0.0;
    d-=tv1;
    return norm(d).raw() <= tv2.error().raw();
}


Bool inconsistent(const TaylorModel<ValidatedFloat>& tm1, const TaylorModel<ValidatedFloat>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (norm(tm1-tm2) > (tm1.error()+tm2.error())*2u);
}

TaylorModel<ValidatedFloat> refinement(const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y) {
    TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());
    Float twice_max_error=0.0;

    const Float& xe=x.error().raw();
    const Float& ye=y.error().raw();
    VOLATILE Float rv,xv,yv,xu,yu,mxl,myl,u,ml;
    //const MultiIndex* aptr;
    MultiIndex a;

    TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin();
    TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() || yiter!=y.end()) {
        // Can't use const MultiIndex& here since references change as the iterators change
        // We would need to use a smart reference
        //const MultiIndex& xa=xiter->key();
        //const MultiIndex& ya=yiter->key();
        if(xiter==x.end()) {
            a=yiter->key();
            yv=yiter->data().raw();
            xv=0.0;
            ++yiter;
        } else if(yiter==y.end()) {
            a=xiter->key();
            xv=xiter->data().raw();
            yv=0.0;
            ++xiter;
        } else if(xiter->key()==yiter->key()) {
            a=xiter->key();
            xv=xiter->data().raw();
            yv=yiter->data().raw();
            ++xiter;
            ++yiter;
        } else if(xiter->key()<yiter->key()) {
            a=xiter->key();
            xv=xiter->data().raw();
            yv=0.0;
            ++xiter;
        } else { // xa>ya
            a=yiter->key();
            yv=yiter->data().raw();
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
            ARIADNE_THROW(IntersectionException,"intersection(TaylorModel<ValidatedFloat>,TaylorModel<ValidatedFloat>)",x<<" and "<<y<<" are disjoint.");
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
        if(rv!=0.0) { r.expansion().raw().append(a,(u-ml)/2); }
    }

    r.error().raw()=twice_max_error/2;

    return r;
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

OutputStream&
operator<<(OutputStream& os, const TaylorModel<ValidatedFloat>& tm) {
    // Set the variable names to be 'parameter' s0,s1,..
    Array<StringType> variable_names(tm.argument_size());
    for(SizeType j=0; j!=tm.argument_size(); ++j) {
        StringStream sstr;
        sstr << 's' << j;
        variable_names[j]=sstr.str();
    }

    //os << "TaylorModel<ValidatedFloat>";
    os << "TM["<<tm.argument_size()<<"](";
    Expansion<RawFloat> e=tm.expansion().raw();
    e.graded_sort();
    e.write(os,variable_names);
    return os << "+/-" << tm.error().raw() << ")";
}


///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


Vector<TaylorModel<ValidatedFloat>> TaylorModel<ValidatedFloat>::zeros(SizeType rs, SizeType as, Sweeper swp)
{
    Vector<TaylorModel<ValidatedFloat>> result(rs,TaylorModel<ValidatedFloat>::zero(as,swp));
    return result;
}


Vector<TaylorModel<ValidatedFloat>> TaylorModel<ValidatedFloat>::constants(SizeType as, const Vector<ValidatedNumber>& c, Sweeper swp)
{
    Vector<TaylorModel<ValidatedFloat>> result(c.size(),TaylorModel<ValidatedFloat>::zero(as,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel<ValidatedFloat>::constant(as,c[i],swp);
    }
    return result;
}

Vector<TaylorModel<ValidatedFloat>> TaylorModel<ValidatedFloat>::coordinates(SizeType as, Sweeper swp)
{
    Vector<TaylorModel<ValidatedFloat>> result(as,TaylorModel<ValidatedFloat>::zero(as,swp));
    for(SizeType i=0; i!=as; ++i) { result[i]=TaylorModel<ValidatedFloat>::coordinate(as,i,swp); }
    return result;
}

Vector<TaylorModel<ValidatedFloat>> TaylorModel<ValidatedFloat>::scalings(const Vector<ExactInterval>& d, Sweeper swp)
{
    Vector<TaylorModel<ValidatedFloat>> result(d.size(),TaylorModel<ValidatedFloat>::zero(d.size(),swp));
    for(SizeType i=0; i!=d.size(); ++i) {
        result[i]=TaylorModel<ValidatedFloat>::scaling(d.size(),i,d[i],swp);
    }
    return result;
}




///////////////////////////////////////////////////////////////////////////////

// Vector-valued versions of scalar operators


Vector<ValidatedNumber> evaluate(const Vector<TaylorModel<ValidatedFloat>>& tv, const Vector<ValidatedNumber>& x)
{
    Vector<ValidatedNumber> r(tv.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=evaluate(tv[i],x); }
    return r;
}

Vector<TaylorModel<ValidatedFloat>> partial_evaluate(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType k, const ValidatedNumber& c)
{
    Vector<TaylorModel<ValidatedFloat>> r(x.size(),ValidatedTaylorModel::zero(x.zero_element().argument_size()-1,x.zero_element().sweeper()));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=partial_evaluate(x[i],k,c); }
    return r;
}

Vector<TaylorModel<ValidatedFloat>> embed(SizeType as1, const Vector<TaylorModel<ValidatedFloat>>& x2, SizeType as3) {
    Vector<TaylorModel<ValidatedFloat>> r(x2.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=embed(as1,x2[i],as3); }
    return r;
}

Vector<TaylorModel<ValidatedFloat>> combine(const Vector<TaylorModel<ValidatedFloat>>& x1, const Vector<TaylorModel<ValidatedFloat>>& x2) {
    return join(embed(0u,x1,x2.zero_element().argument_size()),embed(x1.zero_element().argument_size(),x2,0u));
}

Vector<TaylorModel<ValidatedFloat>> combine(const Vector<TaylorModel<ValidatedFloat>>& x1, const TaylorModel<ValidatedFloat>& x2) {
    return join(embed(0u,x1,x2.argument_size()),embed(x1.zero_element().argument_size(),x2,0u));
}

Vector<TaylorModel<ValidatedFloat>> combine(const TaylorModel<ValidatedFloat>& x1, const Vector<TaylorModel<ValidatedFloat>>& x2) {
    return join(embed(0u,x1,x2.zero_element().argument_size()),embed(x1.argument_size(),x2,0u));
}

Vector<TaylorModel<ValidatedFloat>> combine(const TaylorModel<ValidatedFloat>& x1, const TaylorModel<ValidatedFloat>& x2) {
    return {embed(0u,x1,x2.argument_size()),embed(x1.argument_size(),x2,0u)};
}

Bool refines(const Vector<TaylorModel<ValidatedFloat>>& tv1, const Vector<TaylorModel<ValidatedFloat>>& tv2) {
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(SizeType i=0; i!=tv1.size(); ++i) { if(not refines(tv1[i],tv2[i])) { return false; } }
    return true;
}

Bool inconsistent(const Vector<TaylorModel<ValidatedFloat>>& tv1, const Vector<TaylorModel<ValidatedFloat>>& tv2) {
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(SizeType i=0; i!=tv1.size(); ++i) { if(decide(inconsistent(tv1[i],tv2[i]))) { return true; } }
    return false;
}

Vector<TaylorModel<ValidatedFloat>>
split(const Vector<TaylorModel<ValidatedFloat>>& tv, SizeType j, SplitPart h) {
    Vector<TaylorModel<ValidatedFloat>> r(tv.size());
    for(SizeType i=0; i!=tv.size(); ++i) { r[i]=split(tv[i],j,h); }
    return r;
}

Vector<TaylorModel<ValidatedFloat>>
unscale(const Vector<TaylorModel<ValidatedFloat>>& tvs, const Vector<ExactInterval>& ivls) {
    Vector<TaylorModel<ValidatedFloat>> r(tvs.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=unscale(tvs[i],ivls[i]); }
    return r;
}


Vector<TaylorModel<ValidatedFloat>> antiderivative(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType k) {
    Vector<TaylorModel<ValidatedFloat>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return r;
}

Vector<TaylorModel<ValidatedFloat>> derivative(const Vector<TaylorModel<ValidatedFloat>>& x, SizeType k) {
    Vector<TaylorModel<ValidatedFloat>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=derivative(x[i],k); }
    return r;
}

Vector<UpperInterval> ranges(const Vector<TaylorModel<ValidatedFloat>>& f) {
    Vector<UpperInterval> r(f.size()); for(SizeType i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

Vector<ErrorFloat> errors(const Vector<TaylorModel<ValidatedFloat>>& h) {
    Vector<ErrorFloat> e(h.size()); for(SizeType i=0; i!=h.size(); ++i) { e[i]=h[i].error(); } return e; }

Vector<ErrorFloat> norms(const Vector<TaylorModel<ValidatedFloat>>& h) {
    Vector<ErrorFloat> r(h.size()); for(SizeType i=0; i!=h.size(); ++i) { r[i]=norm(h[i]); } return r; }

ErrorFloat norm(const Vector<TaylorModel<ValidatedFloat>>& h) {
    ErrorFloat r=0u; for(SizeType i=0; i!=h.size(); ++i) { r=max(r,norm(h[i])); } return r;
}

///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumber>
jacobian(const Vector<TaylorModel<ValidatedFloat>>& f, const Vector<ValidatedNumber>& x)
{
    Vector< Differential<ValidatedNumber> > dx=Differential<ValidatedNumber>::variables(1u,x);
    Vector< Differential<ValidatedNumber> > df(f.size(),x.size(),1u);
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumber> J=jacobian(df);
    return J;
}

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumber>
jacobian(const Vector<TaylorModel<ValidatedFloat>>& f, const Vector<ValidatedNumber>& x, const Array<SizeType>& p)
{
    Vector<Differential<ValidatedNumber>> dx(x.size());
    for(SizeType j=0; j!=x.size(); ++j) {
        dx[j]=Differential<ValidatedNumber>::constant(p.size(),1u,x[j]); }
    for(SizeType k=0; k!=p.size(); ++k) {
        SizeType j=p[k];
        dx[j]=Differential<ValidatedNumber>::variable(p.size(),1u,x[j],k); }
    Vector< Differential<ValidatedNumber> > df(f.size());
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumber> J=jacobian(df);
    return J;
}

// Compute the Jacobian at the origin
Matrix<ExactFloat>
jacobian_value(const Vector<TaylorModel<ValidatedFloat>>& f)
{
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<ExactFloat> J(rs,as);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType j=0; j!=as; ++j) {
            a[j]=1; const ExactFloat x=f[i][a]; J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian at the origin with respect to the variables args.
Matrix<ExactFloat>
jacobian_value(const Vector<TaylorModel<ValidatedFloat>>& f, const Array<SizeType>& p)
{
    const SizeType rs=f.size();
    const SizeType as=f.zero_element().argument_size();
    const SizeType ps=p.size();
    Matrix<ExactFloat> J(rs,ps);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType k=0; k!=ps; ++k) {
            SizeType j=p[k]; a[j]=1; const ExactFloat x=f[i][a]; J[i][k]=x; a[j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
Matrix<UpperInterval>
jacobian_range(const Vector<TaylorModel<ValidatedFloat>>& f)
{
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<UpperInterval> J(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        for(TaylorModel<ValidatedFloat>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            MultiIndex const& a=iter->key();
            for(SizeType k=0; k!=as; ++k) {
                const Nat c=a[k];
                if(c>0) {
                    const ExactFloat& x=iter->data();
                    if(a.degree()==1) { J[i][k]+=UpperInterval(x); }
                    else { J[i][k]+=UpperInterval(-1,1)*x*c; }
                }
            }
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain, with respect to the variables p.
Matrix<UpperInterval>
jacobian_range(const Vector<TaylorModel<ValidatedFloat>>& f, const Array<SizeType>& p)
{
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    SizeType ps=p.size();
    Matrix<UpperInterval> J(rs,ps);
    for(SizeType i=0; i!=rs; ++i) {
        for(TaylorModel<ValidatedFloat>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            MultiIndex const& a=iter->key();
            for(SizeType k=0; k!=ps; ++k) {
                SizeType j=p[k];
                const Nat c=a[j];
                if(c>0) {
                    const ExactFloat& x=iter->data();
                    if(a.degree()==1) { J[i][k]+=UpperInterval(x); }
                    else { J[i][k]+=UpperInterval(-1,1)*x*c; }
                }
            }
        }
    }
    return J;
}


ValidatedNumber gradient(const TaylorModel<ValidatedFloat>& f, const Vector<ValidatedNumber>& x, SizeType k) {
    return jacobian({f},x,{k})[0][0];
}

ExactFloat gradient_value(const TaylorModel<ValidatedFloat>& f, SizeType k) {
    return jacobian_value({f},{k})[0][0];
}

UpperInterval gradient_range(const TaylorModel<ValidatedFloat>& f, SizeType k) {
    return jacobian_range({f},{k})[0][0];
}

///////////////////////////////////////////////////////////////////////////////

// Vector operators (evaluate, compose)


// Compose by computing each and every term individually without caching
// Easy to implement, but far too slow
inline
Vector<TaylorModel<ValidatedFloat>>
_compose1(const Vector<TaylorModel<ValidatedFloat>>& x,
          const Vector<TaylorModel<ValidatedFloat>>& ys)
{
    //std::cerr<<"compose1"<<std::endl;
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(ys.size()==x.zero_element().argument_size());
    for(SizeType i=0; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x.zero_element().argument_size()); }
    for(SizeType i=0; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys.zero_element().argument_size(),"ys="<<ys); }

    SizeType as=ys.zero_element().argument_size();
    Sweeper swp=ys.zero_element().sweeper();

    Vector<TaylorModel<ValidatedFloat>> r(x.size(),TaylorModel<ValidatedFloat>(as,swp));
    TaylorModel<ValidatedFloat> t(as,swp);
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]._set_error(x[i].error().raw());
        for(TaylorModel<ValidatedFloat>::ConstIterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            t=static_cast<ExactNumber>(iter->data().raw());
            for(SizeType j=0; j!=iter->key().size(); ++j) {
                TaylorModel<ValidatedFloat> p=pow(ys[j],iter->key()[j]);
                t=t*p;
            }
            r[i]+=t;
        }
    }

    return r;
}


inline
Vector<TaylorModel<ValidatedFloat>>
_compose2(const Vector<TaylorModel<ValidatedFloat>>& x,
          const Vector<TaylorModel<ValidatedFloat>>& ys)
{
    //std::cerr<<"compose2"<<std::endl;
    SizeType yrs=ys.size();
    SizeType xas=ys.size();
    SizeType as=ys.zero_element().argument_size();
    Sweeper sweeper=ys.zero_element().sweeper();

    Array<uchar> max_power(ys.size());
    for(SizeType j=0; j!=ys.size(); ++j) { max_power[j]=1; }

    for(SizeType i=0; i!=x.size(); ++i) {
        for(TaylorModel<ValidatedFloat>::ConstIterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            assert(xas==iter->key().size());
            for(SizeType j=0; j!=iter->key().size(); ++j) {
                max_power[j]=max(max_power[j],iter->key()[j]);
            }
        }
    }

    uchar max_max_power = 1;
    for(SizeType j=0; j!=ys.size(); ++j) {
        max_max_power = max(max_max_power,max_power[j]);
    }

    Array< Array< TaylorModel<ValidatedFloat> > > powers(yrs, Array<TaylorModel<ValidatedFloat>>(max_max_power+1,TaylorModel<ValidatedFloat>::zero(as,sweeper)));
    for(SizeType j=0; j!=yrs; ++j) {
        powers[j][0]=ys[j]*0;
        powers[j][1]=ys[j];
        for(SizeType k=2; k!=powers[j].size(); ++k) {
            powers[j][k]=powers[j][k/2]*powers[j][(k+1)/2];
        }
    }

    Vector<TaylorModel<ValidatedFloat>> r(x.size(),TaylorModel<ValidatedFloat>(as,sweeper));
    TaylorModel<ValidatedFloat> t(as,sweeper);
    MultiIndex a;
    Float c;
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]._set_error(x[i].error().raw());
        for(TaylorModel<ValidatedFloat>::ConstIterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            a=iter->key();
            c=iter->data().raw();
            t=static_cast<ExactNumber>(c);
            for(SizeType j=0; j!=a.size(); ++j) {
                if(a[j]>0) {
                    t=t*powers[j][a[j]];
                }
            }
            r[i]+=t;
        }
    }

    return r;
}

Vector<TaylorModel<ValidatedFloat>>
_compose(const Vector<TaylorModel<ValidatedFloat>>& x,
         const Vector<TaylorModel<ValidatedFloat>>& ys)
{
    ARIADNE_ASSERT_MSG(x.size()>0,"x="<<x<<", ys="<<ys);
    ARIADNE_ASSERT_MSG(ys.size()==x[0].argument_size(),"x="<<x<<", ys="<<ys);
    for(SizeType i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }
    for(SizeType i=1; i!=ys.size(); ++i) { ARIADNE_ASSERT_MSG(ys[i].argument_size()==ys[0].argument_size(),"ys="<<ys); }

    return _compose2(x,ys);
}


Vector<TaylorModel<ValidatedFloat>>
compose(const VectorUnscaling& u,
        const Vector<TaylorModel<ValidatedFloat>>& y)
{
    ARIADNE_ASSERT_MSG(u.size()==y.size(),"u="<<u.domain()<<", y="<<y);
    Vector<TaylorModel<ValidatedFloat>> r(u.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(u[i],y[i]); }
    return std::move(r);
}



Vector<TaylorModel<ValidatedFloat>>
compose(const Vector<TaylorModel<ValidatedFloat>>& x,
        const Vector<TaylorModel<ValidatedFloat>>& y)
{
    Vector<TaylorModel<ValidatedFloat>> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=compose(x[i],y); }
    return std::move(r);
}



TaylorModel<ValidatedFloat>
compose(const TaylorModel<ValidatedFloat>& x,
        const VectorUnscaling& u,
        const Vector<TaylorModel<ValidatedFloat>>& y)
{
    return compose(x,compose(u,y));
}

Vector<TaylorModel<ValidatedFloat>>
compose(const Vector<TaylorModel<ValidatedFloat>>& x,
        const VectorUnscaling& u,
        const Vector<TaylorModel<ValidatedFloat>>& y)
{
    return compose(x,compose(u,y));
}




template TaylorModel<ValidatedFloat> rec(const TaylorModel<ValidatedFloat>&);
template TaylorModel<ValidatedFloat> sqrt(const TaylorModel<ValidatedFloat>&);
template TaylorModel<ValidatedFloat> exp(const TaylorModel<ValidatedFloat>&);
template TaylorModel<ValidatedFloat> log(const TaylorModel<ValidatedFloat>&);
template TaylorModel<ValidatedFloat> sin(const TaylorModel<ValidatedFloat>&);
template TaylorModel<ValidatedFloat> cos(const TaylorModel<ValidatedFloat>&);
template TaylorModel<ValidatedFloat> tan(const TaylorModel<ValidatedFloat>&);







TaylorModel<ApproximateFloat>::TaylorModel(SizeType as)
    : _expansion(as), _sweeper()
{
}

TaylorModel<ApproximateFloat>::TaylorModel(SizeType as, Sweeper swp)
    : _expansion(as), _sweeper(swp)
{
}

TaylorModel<ApproximateFloat> TaylorModel<ApproximateFloat>::create_constant(ApproximateFloat c) const {
    TaylorModel<ApproximateFloat> r(this->argument_size(),this->_sweeper);
    r._expansion.append(MultiIndex::zero(this->argument_size()),c);
}

TaylorModel<ApproximateFloat> TaylorModel<ApproximateFloat>::create_ball(ErrorType) const {
    return TaylorModel<ApproximateFloat>(this->argument_size(),this->_sweeper);
}


Void TaylorModel<ApproximateFloat>::iadd(const ApproximateFloat& c)
{
    // Compute self+=c
    TaylorModel<ApproximateFloat>& x=*this;
    if(c==0) { return; }
    if(x._expansion.empty()) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else if((x._expansion.end()-1)->key().degree()>0) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else {
        ApproximateFloat& rv=(x._expansion.end()-1)->data();
        rv+=c;
    }
}

Void TaylorModel<ApproximateFloat>::imul(const ApproximateFloat& c)
{
    // Compute self*=c
    if(c==0) { this->clear(); return; }
    if(c==1) { return; }
    for(ExpansionType::Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        iter->data() *= c;
    }
}



Void TaylorModel<ApproximateFloat>::isma(const ApproximateFloat& c, const TaylorModel<ApproximateFloat>& y)
{
    const TaylorModel<ApproximateFloat>& x=*this;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<", y="<<y);
    TaylorModel<ApproximateFloat> r(x.argument_size());

    ConstIterator xiter=x._expansion.begin();
    ConstIterator yiter=y._expansion.begin();
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


Void TaylorModel<ApproximateFloat>::ifma(const TaylorModel<ApproximateFloat>& x, const TaylorModel<ApproximateFloat>& y)
{
    TaylorModel<ApproximateFloat>& r(*this);
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<",y="<<y);

    TaylorModel<ApproximateFloat> t(x.argument_size());
    MultiIndex sa;
    for(ExpansionType::ConstIterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        ExpansionType::ConstIterator yiter=y._expansion.begin();
        ExpansionType::ConstIterator riter=r._expansion.begin();
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
            t._expansion.append_sum(xiter->key(),yiter->key(),xiter->data()*yiter->data());
            ++yiter;

        }

        r._expansion.swap(t._expansion);
        t._expansion.clear();
    }
}




TaylorModel<ApproximateFloat>::NormType TaylorModel<ApproximateFloat>::norm() const {
    return ExactFloat(((*this)[MultiIndex(this->argument_size())]).raw());
}

TaylorModel<ApproximateFloat>::CoefficientType TaylorModel<ApproximateFloat>::average() const {
    return ExactFloat(((*this)[MultiIndex(this->argument_size())]).raw());
}

TaylorModel<ApproximateFloat>::NormType TaylorModel<ApproximateFloat>::radius() const {
    ApproximateFloat r=0.0;
    for(Expansion<ApproximateFloat>::ConstIterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=abs(iter->data());
        }
    }
    return ErrorType(r.raw());
}

ApproximateInterval TaylorModel<ApproximateFloat>::range() const {
    ApproximateFloat av=this->average();
    ApproximateFloat rad=this->radius();
    return ApproximateInterval(av-rad,av+rad);
}

Float TaylorModel<ApproximateFloat>::tolerance() const {
    return dynamic_cast<const ThresholdSweeper&>(static_cast<const SweeperInterface&>(this->_sweeper)).sweep_threshold();
}


OutputStream& TaylorModel<ApproximateFloat>::write(OutputStream& os) const {
    return os << "TM["<<this->argument_size()<<"](" << this->_expansion << ")";
}


template TaylorModel<ApproximateFloat> neg(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> rec(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> sqr(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> pow(const TaylorModel<ApproximateFloat>&, Int);
template TaylorModel<ApproximateFloat> sqrt(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> exp(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> log(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> sin(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> cos(const TaylorModel<ApproximateFloat>&);
template TaylorModel<ApproximateFloat> tan(const TaylorModel<ApproximateFloat>&);


OutputStream& TaylorModel<ApproximateFloat>::str(OutputStream& os) const {
    return os << this->_expansion;
}



} //namespace Ariadne


