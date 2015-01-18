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

static const double MACHINE_EPSILON = 2.2204460492503131e-16;

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }

} // namespace




Vector<ValidatedNumber> unscale(const Vector<ValidatedNumber>& x, const Vector<ExactInterval>& d) {
    Vector<ValidatedNumber> r(x);
    for(SizeType i=0; i!=r.size(); ++i) {
        if(d[i].lower()==d[i].upper()) {
            if(x[i]==d[i].midpoint()) {
                r[i]=ValidatedNumber(0.0,0.0);
            } else {
                r[i]=ValidatedNumber(-inf,+inf);
            }
        } else {
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
    ARIADNE_DEBUG_PRECONDITION(e.raw()>=0);
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

struct ValidatedApproximateFloat {
    ValidatedFloat _v; ApproximateFloat _a;
    ValidatedApproximateFloat(ValidatedFloat x) : _v(x), _a(x) { }
    LowerFloat lower() const { return _v.lower(); }
    ApproximateFloat middle() const { return _a; }
    UpperFloat upper() const { return _v.upper(); }
    Float const& lower_raw() const { return _v.lower_raw(); }
    Float const& middle_raw() const { return _a.raw(); }
    Float const& upper_raw() const { return _v.upper_raw(); }
};

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

ExactFloat add_err(ExactFloat const& x, ValidatedApproximateFloat const& c, ErrorFloat& e) {
    Float const& xv=x.raw();
    Float const& cl=c.lower_raw();
    Float const& cm=c.middle_raw();
    Float const& cu=c.upper_raw();
    Float& re=e.raw();
    set_rounding_to_nearest();
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

ExactFloat mul_err(ExactFloat const& x, ValidatedApproximateFloat const& c, ErrorFloat& e) {
    Float const& xv=x.raw();
    Float const& cu=c.upper_raw();
    Float const& cm=c.middle_raw();
    Float const& cl=c.lower_raw();
    Float& re=e.raw();
    set_rounding_to_nearest();
    Float rv=xv*cm;
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

ExactFloat fma_err(ExactFloat const& x, ValidatedApproximateFloat const& c, ExactFloat y, ErrorFloat& e) {
    Float const& xv=x.raw();
    Float const& cu=c.upper_raw();
    Float const& cm=c.middle_raw();
    Float const& cl=c.lower_raw();
    Float const& yv=y.raw();
    Float& re=e.raw();
    set_rounding_to_nearest();
    Float rv=xv+cm*yv;
    set_rounding_upward();
    Float u,ml;
    if(yv>=0) {
        Float mcl=-cl;
        u=cu*yv+xv;
        ml=mcl*yv-xv;
    } else {
        Float mcu=-cu;
        u=cl*yv+xv;
        ml=mcu*yv-xv;
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
        iter->data()=-iter->data();
    }
}


Void _scal(TaylorModel<ValidatedFloat>& r, const TwoExp& c)
{
    if(ExactFloat(c)==1) { return; }
    for(TaylorModel<ValidatedFloat>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=c;
    }
    r.error()*=ErrorFloat(c);
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
    ValidatedApproximateFloat clmu=c;
    for(TaylorModel<ValidatedFloat>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        ExactFloat& rv=riter->data();
        rv=mul_err(rv,clmu,e);
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
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else if((r.end()-1)->key().degree()>0) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
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

    VOLATILE Float u,ml,myv;
    ErrorFloat te=0; // Twice the maximum accumulated error
    ErrorFloat err=0u; // Twice the maximum accumulated error
    ValidatedApproximateFloat clmu=c;

    // Compute r=x+y, assuming r is empty
    set_rounding_upward();
    TaylorModel<ValidatedFloat>::ConstIterator xiter=x.begin();
    TaylorModel<ValidatedFloat>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            const ExactFloat& xv=xiter->data();
            r.expansion().append(xiter->key(),xv);
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            const ExactFloat& yv=yiter->data();
            r.expansion().append(yiter->key(),mul_err(yv,c,err));
            ++yiter;
        } else {
            const ExactFloat& xv=xiter->data();
            const ExactFloat& yv=yiter->data();
            r.expansion().append(xiter->key(),fma_err(xv,clmu,yv,err));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        const ExactFloat& xv=xiter->data();
        r.expansion().append(xiter->key(),xv);
        ++xiter;
    }
    while(yiter!=y.end()) {
        const ExactFloat& yv=yiter->data();
        r.expansion().append(yiter->key(),mul_err(yv,clmu,err));
        ++yiter;
    }

    r.error()=x.error()+y.error();
    r.error()+=err;

    ARIADNE_DEBUG_ASSERT_MSG(r.error()>=0,r);
}



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

// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
UpperInterval TaylorModel<ValidatedFloat>::range() const {
    const TaylorModel<ValidatedFloat>& tm=*this;
    const SizeType as=tm.argument_size();
    ExactFloat constant_term=0;
    Array<ExactFloat> linear_terms(as,0);
    Array<ExactFloat> quadratic_terms(as,0);
    ErrorFloat err=0u;
    for(auto iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->key().degree()==0) {
            constant_term=iter->data();
        } else if(iter->key().degree()==1) {
            for(SizeType j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==1) { linear_terms[j]=iter->data(); break; }
            }
        } else if(iter->key().degree()==2) {
            for(SizeType j=0; j!=tm.argument_size(); ++j) {
                if(iter->key()[j]==2) { quadratic_terms[j]=iter->data(); break; }
                if(iter->key()[j]==1) { err+=mag(iter->data()); break; }
            }
        } else {
            err+=mag(iter->data());
        }
    }
    err=err+tm.error();
    ValidatedFloat r(constant_term-err,constant_term+err);
    const ValidatedFloat unit_ivl(-1,+1);
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(SizeType j=0; j!=as; ++j) {
        const ExactFloat& a=quadratic_terms[j];
        const ExactFloat& b=linear_terms[j];
        ValidatedFloat ql=abs(a)*unit_ivl + abs(b)*unit_ivl;
        ValidatedFloat qf=a*(sqr(unit_ivl+b/a/2))-sqr(b)/a/4;
        r += refinement(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
    }
    return r;
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
    Sweeper threshold_sweeper(new ThresholdSweeper(MACHINE_EPSILON));
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




// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error.
TaylorModel<ValidatedFloat>
compose(const AnalyticFunction& fn, const TaylorModel<ValidatedFloat>& tm) {
    std::cerr<<"\n\nfn="<<fn<<"\ntm="<<tm<<"\n";

    static const DegreeType MAX_DEGREE=20;
    static const Float MAX_TRUNCATION_ERROR=MACHINE_EPSILON;
    SweeperInterface const& sweeper=tm.sweeper();

    Float max_truncation_error=MAX_TRUNCATION_ERROR;
    ThresholdSweeper const* threshold_sweeper_ptr = dynamic_cast<ThresholdSweeper const*>(&sweeper);
    if(threshold_sweeper_ptr) { max_truncation_error=threshold_sweeper_ptr->sweep_threshold(); }

    DegreeType max_degree=MAX_DEGREE;
    GradedSweeper const* graded_sweeper_ptr = dynamic_cast<GradedSweeper const*>(&sweeper);
    if(graded_sweeper_ptr) { max_degree=graded_sweeper_ptr->degree(); }

    std::cerr<<"max_truncation_error="<<max_truncation_error<<"\nmax_degree="<<(uint)max_degree<<"\n";

    Nat d=max_degree;
    ExactFloat c=tm.value();
    ValidatedFloat r=make_singleton(tm.range());
    Series<ValidatedFloat> centre_series=fn.series(c);
    Series<ValidatedFloat> range_series=fn.series(r);
    std::cerr<<"c="<<c<<"\nr="<<r<<"\n";
    std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";


    ErrorFloat se=mag(range_series[d]-centre_series[d]);
    ErrorFloat e=mag(r-c);
    ErrorFloat p=pow(e,d);
    std::cerr<<"se="<<se<<"\ne="<<e<<"\np="<<p<<"\n";
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    ErrorFloat truncation_error=se*p;
    std::cerr<<"truncation_error="<<truncation_error<<"\n\n";
    if(truncation_error.raw()>max_truncation_error) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<max_truncation_error<<"\n");
    }

    TaylorModel<ValidatedFloat> x=tm-c;
    TaylorModel<ValidatedFloat> res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        // Don't sweep here...
    }
    res+=ValidatedFloat(-truncation_error,+truncation_error);
    return res;
}



///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions are implemented using generic Algebra code

///////////////////////////////////////////////////////////////////////////////

// Inplace operators manipulating the error term

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





///////////////////////////////////////////////////////////////////////////////

// Differentiation operators


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

// Scalar function operators (evaluate, split, unscale, embed)
// and predicates (refines)

ValidatedNumber evaluate(const TaylorModel<ValidatedFloat>& tm, const Vector<ValidatedNumber>& x) {
    return horner_evaluate(tm.expansion(),x)+ValidatedNumber(-tm.error(),+tm.error());
}

ApproximateNumber evaluate(const TaylorModel<ValidatedFloat>& tm, const Vector<ApproximateNumber>& x) {
    return horner_evaluate(tm.expansion(),x);
}


Covector<ValidatedNumber> gradient(const TaylorModel<ValidatedFloat>& tm, const Vector<ValidatedNumber>& x)
{
    Vector< Differential<ValidatedNumber> > dx=Differential<ValidatedNumber>::variables(1u,x);
    Differential<ValidatedNumber> df=horner_evaluate(tm.expansion(),dx)+ValidatedNumber(-tm.error(),+tm.error());
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
    const SizeType as=x.argument_size();
    Vector<TaylorModel<ValidatedFloat>> y(as,TaylorModel<ValidatedFloat>(as-1,x.sweeper()));
    for(SizeType i=0; i!=k; ++i) { y[i]=TaylorModel<ValidatedFloat>::coordinate(as-1,i,x.sweeper()); }
    y[k]=TaylorModel<ValidatedFloat>::constant(as-1,c,x.sweeper());
    for(SizeType i=k+1; i!=as; ++i) { y[i]=TaylorModel<ValidatedFloat>::coordinate(as-1,i-1,x.sweeper()); }
    return compose(x,y);
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
        ExactFloat& c=iter->data();
        c/=two_exp(ak);
    }

    if(h==SplitPart::MIDDLE) { return r; }
    Int tr=( h==SplitPart::UPPER ? +1 : -1 );

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<TaylorModel<ValidatedFloat>> ary(deg+1,TaylorModel<ValidatedFloat>(as,swp));
    for(TaylorModel<ValidatedFloat>::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->key();
        const ExactFloat& c=iter->data();
        uchar ak=a[k];
        a[k]=0u;
        ary[ak].expansion().append(a,c);
    }

    ErrorFloat re=r.error();
    r.clear();
    r.set_error(re);

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
unscale(const TaylorModel<ValidatedFloat>& tm, const ExactInterval& ivl) {
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

    ARIADNE_ASSERT_MSG(ivl.lower()<=ivl.upper(),"Cannot unscale TaylorModel<ValidatedFloat> "<<tm<<" from empty interval "<<ivl);

    if(ivl.lower()==ivl.upper()) {
        TaylorModel<ValidatedFloat> r=tm.create();
        r+=ivl.midpoint();
        // Uncomment out line below to make unscaling to a singleton interval undefined
        //r.set_error(+inf);
        return r;
    } else {
        TaylorModel<ValidatedFloat> r=tm;
        ValidatedNumber c=ivl.centre();
        ValidatedNumber s=rec(ivl.radius());
        r-=c;
        r*=s;

        return r;
    }
}

TaylorModel<ValidatedFloat> compose(const Unscaling& u, const TaylorModel<ValidatedFloat>& y) {
    unscale(y,u.domain());
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

NormType norm(const TaylorModel<ValidatedFloat>& tm) {
    return tm.norm();
}


Bool refines(const TaylorModel<ValidatedFloat>& tm1, const TaylorModel<ValidatedFloat>& tm2)
{
    ARIADNE_ASSERT(tm1.argument_size()==tm2.argument_size());
    TaylorModel<ValidatedFloat> d=tm2;
    d.error()=0u;
    d-=tm1;
    return norm(d).raw() <= tm2.error().raw();
}


Bool inconsistent(const TaylorModel<ValidatedFloat>& tm1, const TaylorModel<ValidatedFloat>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (norm(tm1-tm2) > (tm1.error()+tm2.error())*2u);
}

TaylorModel<ValidatedFloat> refinement(const TaylorModel<ValidatedFloat>& x, const TaylorModel<ValidatedFloat>& y) {
    TaylorModel<ValidatedFloat> r(x.argument_size(),x.sweeper());

    ErrorFloat max_error=0u;

    const ErrorFloat& xe=x.error();
    const ErrorFloat& ye=y.error();
    ExactFloat rv,xv,yv;
    Float xu,yu,mxl,myl,u,ml;
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
            yv=yiter->data();
            xv=0;
            ++yiter;
        } else if(yiter==y.end()) {
            a=xiter->key();
            xv=xiter->data();
            yv=0;
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
            yv=0;
            ++xiter;
        } else { // xa>ya
            a=yiter->key();
            yv=yiter->data();
            xv=0;
            ++yiter;
        }

        ValidatedFloat rve=refinement( xv.pm(xe), yv.pm(ye) );
        if(rve.error().raw()<0.0) {
            ARIADNE_THROW(IntersectionException,"refinement(TaylorModel<ValidatedFloat>,TaylorModel<ValidatedFloat>)",x<<" and "<<y<<" are inconsistent.");
        }

        if(rve.value()!=0) { r.expansion().append(a,rve.value()); }
        max_error=max(max_error,rve.error());
    }

    r.error()=max_error;

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
    Expansion<ExactFloat> e=tm.expansion();
    e.graded_sort();
    e.write(os,variable_names);
    return os << "+/-" << tm.error() << ")";
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
    ApproximateFloat r=0.0;
    for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        r+=mag(iter->data());
    }
    return NormType(r);
}

TaylorModel<ApproximateFloat>::CoefficientType TaylorModel<ApproximateFloat>::average() const {
    return (*this)[MultiIndex(this->argument_size())];
}

TaylorModel<ApproximateFloat>::NormType TaylorModel<ApproximateFloat>::radius() const {
    ApproximateFloat r=0.0;
    for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=mag(iter->data());
        }
    }
    return NormType(r);
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


