/***************************************************************************
 *            taylor_model.tcc
 *
 *  Copyright 2008-15  Pieter Collins
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
#include "algebra/covector.h"
#include "algebra/matrix.h"
#include "algebra/expansion.h"
#include "algebra/series.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"
#include "function/taylor_series.h"
#include "function/function.h"
#include "utility/exceptions.h"

#include "algebra/algebra_operations.tcc"

#define VOLATILE ;
#include "algebra/multi_index-noaliasing.h"
#include "function/function_mixin.h"
#include "algebra/vector.h"

namespace Ariadne {

namespace {

static const double MACHINE_EPSILON = 2.2204460492503131e-16;

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }

} // namespace





template<class F> TaylorModel<ValidatedTag,F>::TaylorModel()
    : _expansion(0), _error(0u), _sweeper()
{
}


template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(SizeType as, Sweeper swp)
    : _expansion(as), _error(0u), _sweeper(swp)
{
}

template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(const Expansion<CoefficientType>& f, const ErrorType& e, Sweeper swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->cleanup();
}

template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(const Expansion<Float64>& f, const Float64& e, Sweeper swp)
    : TaylorModel(reinterpret_cast<Expansion<CoefficientType>const&>(f),reinterpret_cast<ErrorType const&>(e),swp)
{
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::scaling(SizeType as, SizeType j, const ExactIntervalType& codom, Sweeper swp) {
    TaylorModel<ValidatedTag,F> r(as,swp);
    r.set_gradient(j,1);
    r*=codom.radius();
    r+=codom.midpoint();
    return r;
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create() const {
    return TaylorModel<ValidatedTag,F>(this->argument_size(),this->_sweeper);
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create_zero() const {
    return TaylorModel<ValidatedTag,F>(this->argument_size(),this->_sweeper);
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create_constant(NumericType c) const {
    return TaylorModel<ValidatedTag,F>::constant(this->argument_size(),c,this->_sweeper);
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create_coordinate(SizeType j) const {
    ARIADNE_PRECONDITION(j<this->argument_size());
    TaylorModel<ValidatedTag,F> r(this->argument_size(),this->_sweeper);
    r._expansion.append(MultiIndex::unit(this->argument_size(),j),1);
    return r;
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create_ball(ErrorType e) const {
    ARIADNE_DEBUG_PRECONDITION(e.raw()>=0);
    TaylorModel<ValidatedTag,F> r(this->argument_size(),this->_sweeper);
    r._error=e;
    return r;
}

template<class F> Void TaylorModel<ValidatedTag,F>::swap(TaylorModel<ValidatedTag,F>& tm) {
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
    std::swap(this->_sweeper,tm._sweeper);
}

template<class F> Void TaylorModel<ValidatedTag,F>::clear() {
    this->_expansion.clear();
    this->_error=0u;
}

template<class F> DegreeType TaylorModel<ValidatedTag,F>::degree() const {
    DegreeType deg=0u;
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        deg=std::max(deg,iter->key().degree());
    }
    return deg;
}


template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::operator=(const ValidatedNumericType& c) {
    this->_expansion.clear();
    ExactFloat64 m=c.value();
    if(m!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    this->_error=c.error();
    return *this;
}


namespace { // Internal code for arithmetic

struct ValidatedApproximateFloat {
    BoundedFloat64 _v; ApproximateFloat64 _a;
    ValidatedApproximateFloat(BoundedFloat64 x) : _v(x), _a(x) { }
    LowerFloat64 lower() const { return _v.lower(); }
    ApproximateFloat64 middle() const { return _a; }
    UpperFloat64 upper() const { return _v.upper(); }
    Float64 const& lower_raw() const { return _v.lower_raw(); }
    Float64 const& middle_raw() const { return _a.raw(); }
    Float64 const& upper_raw() const { return _v.upper_raw(); }
};

ExactFloat64 add_err(ExactFloat64 const& x1, ExactFloat64 const& x2, ErrorFloat64& e) {
    ExactFloat64 mx1=-x1;
    Float64::set_rounding_to_nearest();
    ExactFloat64 r(x1.raw() + x2.raw());
    Float64::set_rounding_upward();
    Float64 u=x1.raw()+x2.raw();
    Float64 ml=mx1.raw()-x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

ExactFloat64 add_err(ExactFloat64 const& x, ValidatedApproximateFloat const& c, ErrorFloat64& e) {
    Float64 const& xv=x.raw();
    Float64 const& cl=c.lower_raw();
    Float64 const& cm=c.middle_raw();
    Float64 const& cu=c.upper_raw();
    Float64& re=e.raw();
    Float64::set_rounding_to_nearest();
    Float64 rv=xv+cm;
    Float64::set_rounding_upward();
    Float64 u=xv+cu;
    Float64 ml=(-xv)-cl;
    re += (u+ml)/2;
    return ExactFloat64(rv);
}

ExactFloat64 sub_err(ExactFloat64 const& x1, ExactFloat64 const& x2, ErrorFloat64& e) {
    ExactFloat64 mx1=-x1;
    Float64::set_rounding_to_nearest();
    ExactFloat64 r(x1.raw() - x2.raw());
    Float64::set_rounding_upward();
    Float64 u=x1.raw()-x2.raw();
    Float64 ml=mx1.raw()+x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

ExactFloat64 mul_no_err(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    Float64::set_rounding_to_nearest();
    ExactFloat64 r(x1.raw() * x2.raw());
    Float64::set_rounding_upward();
    return r;
}

ExactFloat64 mul_err(ExactFloat64 const& x1, ExactFloat64 const& x2, ErrorFloat64& e) {
    ExactFloat64 mx1=-x1;
    Float64::set_rounding_to_nearest();
    ExactFloat64 r(x1.raw() * x2.raw());
    Float64::set_rounding_upward();
    Float64 u=x1.raw()*x2.raw();
    Float64 ml=mx1.raw()*x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

ExactFloat64 mul_err(ExactFloat64 const& x, ValidatedApproximateFloat const& c, ErrorFloat64& e) {
    Float64 const& xv=x.raw();
    Float64 const& cu=c.upper_raw();
    Float64 const& cm=c.middle_raw();
    Float64 const& cl=c.lower_raw();
    Float64& re=e.raw();
    Float64::set_rounding_to_nearest();
    Float64 rv=xv*cm;
    Float64::set_rounding_upward();
    Float64 u,ml;
    if(xv>=0) {
        Float64 mcl=-cl;
        u=xv*cu;
        ml=xv*mcl;
    } else {
        Float64 mcu=-cu;
        u=xv*cl;
        ml=xv*mcu;
    }
    re+=(u+ml)/2;
    return ExactFloat64(rv);
}

ExactFloat64 fma_err(ExactFloat64 const& x, ValidatedApproximateFloat const& c, ExactFloat64 y, ErrorFloat64& e) {
    Float64 const& xv=x.raw();
    Float64 const& cu=c.upper_raw();
    Float64 const& cm=c.middle_raw();
    Float64 const& cl=c.lower_raw();
    Float64 const& yv=y.raw();
    Float64& re=e.raw();
    Float64::set_rounding_to_nearest();
    Float64 rv=xv+cm*yv;
    Float64::set_rounding_upward();
    Float64 u,ml;
    if(yv>=0) {
        Float64 mcl=-cl;
        u=cu*yv+xv;
        ml=mcl*yv-xv;
    } else {
        Float64 mcu=-cu;
        u=cl*yv+xv;
        ml=mcu*yv-xv;
    }
    re+=(u+ml)/2;
    return ExactFloat64(rv);
}

ExactFloat64 div_err(ExactFloat64 const& x1, ExactFloat64 const& x2, ErrorFloat64& e) {
    ExactFloat64 mx1=-x1;
    Float64::set_rounding_to_nearest();
    ExactFloat64 r(x1.raw() / x2.raw());
    Float64::set_rounding_upward();
    Float64 u=x1.raw()/x2.raw();
    Float64 ml=mx1.raw()/x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

// Inplace negation
template<class F> Void _neg(TaylorModel<ValidatedTag,F>& r)
{
    for(auto iter=r.begin(); iter!=r.end(); ++iter) {
        iter->data()=-iter->data();
    }
}


template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const TwoExp& c)
{
    if(ExactFloat64(c)==1) { return; }
    for(typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data()*=c;
    }
    r.error()*=ErrorFloat64(c);
 }



template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const ExactFloat64& c) {
    ErrorFloat64 e=0u; // The maximum accumulated error
    for(typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->data() = mul_err(riter->data(),c,e);
    }
    ErrorFloat64& re=r.error();
    re*=abs(c);
    re+=e;
}

template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const BoundedFloat64& c)
{
    //std::cerr<<"TaylorModel<ValidatedTag,F>::scal(BoundedFloat64 c) c="<<c<<std::endl;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);

    if(r.error().raw()==inf) {
        r.expansion().clear(); return;
    }

    ErrorFloat64 e=0u;
    ValidatedApproximateFloat clmu=c;
    for(typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        ExactFloat64& rv=riter->data();
        rv=mul_err(rv,clmu,e);
    }
    r.error()*=mag(c);
    r.error()+=e;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


struct UnitMultiIndex { SizeType argument_size; SizeType unit_index; };


template<class F> inline Void _incr(TaylorModel<ValidatedTag,F>& r, const MultiIndex& a) {
    for(typename TaylorModel<ValidatedTag,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        static_cast<MultiIndex&>(iter->key())+=a;
    }
}

template<class F> inline Void _incr(TaylorModel<ValidatedTag,F>& r, SizeType j) {
    for(typename TaylorModel<ValidatedTag,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++static_cast<MultiIndex&>(iter->key())[j];
    }
}


template<class F> inline Void _acc(TaylorModel<ValidatedTag,F>& r, const ExactFloat64& c) {
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    if(c==0) { return; }
    if(r.expansion().empty()) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else if((r.end()-1)->key().degree()>0) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else {
        ExactFloat64& rv=(r.end()-1)->data();
        ErrorFloat64& re=r.error();
        rv=add_err(rv,c,re);
    }
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    return;
}



template<class F> inline Void _acc(TaylorModel<ValidatedTag,F>& r, const BoundedFloat64& c)
{
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);

    if(c.lower().raw()==-inf || c.upper().raw()==+inf) {
        r.clear();
        r.set_error(mag(infty));
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

    ExactFloat64& rv=(r.end()-1)->data();
    ErrorFloat64& re=r.error();
    rv=add_err(rv,c,re);

    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}


// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
template<class F> inline Void _add(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y)
{
    ARIADNE_PRECONDITION(r.number_of_nonzeros()==0);
    ErrorFloat64 e=0u;
    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
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

    Float64::set_rounding_upward();
    r.error()=(x.error()+y.error())+e;
    Float64::set_rounding_to_nearest();

    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


template<class F> inline Void _acc(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x)
{
    TaylorModel<ValidatedTag,F> s(r.argument_size(),r.sweeper()); _add(s,r,x); s.swap(r);
    ARIADNE_DEBUG_ASSERT_MSG(r.error()>=0,r);
}

template<class F> inline Void _sub(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y)
{
    ARIADNE_PRECONDITION(r.number_of_nonzeros()==0);
    ErrorFloat64 e=0u;
    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
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

    Float64::set_rounding_upward();
    r.error()=(x.error()+y.error())+e;
    Float64::set_rounding_to_nearest();

    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


template<class F> inline Void _sma(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x, const BoundedFloat64& c, const TaylorModel<ValidatedTag,F>& y)
{
    ARIADNE_ASSERT_MSG(c.lower().raw()<=c.upper().raw(),c);
    ARIADNE_ASSERT_MSG(x.error().raw()>=0,"x="<<x);
    ARIADNE_ASSERT_MSG(y.error().raw()>=0,"y="<<y);

    VOLATILE Float64 u,ml,myv;
    ErrorFloat64 te=0u; // Twice the maximum accumulated error
    ErrorFloat64 err=0u; // Twice the maximum accumulated error
    ValidatedApproximateFloat clmu=c;

    // Compute r=x+y, assuming r is empty
    Float64::set_rounding_upward();
    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()<yiter->key()) {
            const ExactFloat64& xv=xiter->data();
            r.expansion().append(xiter->key(),xv);
            ++xiter;
        } else if(yiter->key()<xiter->key()) {
            const ExactFloat64& yv=yiter->data();
            r.expansion().append(yiter->key(),mul_err(yv,c,err));
            ++yiter;
        } else {
            const ExactFloat64& xv=xiter->data();
            const ExactFloat64& yv=yiter->data();
            r.expansion().append(xiter->key(),fma_err(xv,clmu,yv,err));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        const ExactFloat64& xv=xiter->data();
        r.expansion().append(xiter->key(),xv);
        ++xiter;
    }
    while(yiter!=y.end()) {
        const ExactFloat64& yv=yiter->data();
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
template<class F> inline Void _mul(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y)
{
    const SizeType as=r.argument_size();
    TaylorModel<ValidatedTag,F> t(as,r.sweeper());
    TaylorModel<ValidatedTag,F> s(as,r.sweeper());
    MultiIndex ta(as);
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        ErrorFloat64 te=0u; // trucation error
        ErrorFloat64 re=0u; // roundoff error
        const MultiIndex& xa=xiter->key();
        const ExactFloat64& xv=xiter->data();
        ErrorFloat64 mag_xv=mag(xv);
        for(typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const MultiIndex& ya=yiter->key();
            const ExactFloat64& yv=yiter->data();
            ta=xa+ya;
            ExactFloat64 tv=mul_no_err(xv,yv);
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

    ErrorFloat64 xs=0u;
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=mag(xiter->data());
    }

    ErrorFloat64 ys=0u;
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=mag(yiter->data());
    }

    ErrorFloat64& re=r.error();
    const ErrorFloat64& xe=x.error();
    const ErrorFloat64& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    return;
}


} // namespace


///////////////////////////////////////////////////////////////////////////////

// Inplace arithmetical operations for Algebra concept

template<class F> Void TaylorModel<ValidatedTag,F>::iadd(const ValidatedNumericType& c)
{
    _acc(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}

template<class F> Void TaylorModel<ValidatedTag,F>::imul(const ValidatedNumericType& c)
{
    _scal(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}

template<class F> Void TaylorModel<ValidatedTag,F>::isma(const ValidatedNumericType& c, const TaylorModel<ValidatedTag,F>& y)
{
    TaylorModel<ValidatedTag,F>& x=*this;
    TaylorModel<ValidatedTag,F> r=this->create();
    _sma(r,x,c,y);
    this->swap(r);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}

template<class F> Void TaylorModel<ValidatedTag,F>::ifma(const TaylorModel<ValidatedTag,F>& x1, const TaylorModel<ValidatedTag,F>& x2)
{
    _mul(*this,x1,x2);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error()>=0,*this);
}



///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::sort() {
    this->_expansion.reverse_lexicographic_sort();
};

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::unique()
{
    typename TaylorModel<ValidatedTag,F>::ConstIterator advanced =this->begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator end =this->end();
    typename TaylorModel<ValidatedTag,F>::Iterator current=this->begin();
    ErrorFloat64 e=0u;
    while(advanced!=end) {
        current->key()=advanced->key();
        ExactFloat64 rv=advanced->data();
        ++advanced;
        while(advanced!=end && advanced->key()==current->key()) {
            const ExactFloat64& xv=advanced->data();
            rv=add_err(rv,xv,e);
            ++advanced;
        }
        current->data()=rv;
        ++current;
    }
    this->error()+=e;
    this->_expansion.resize(current-this->begin());

    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::sweep() {
    this->_sweeper.sweep(this->_expansion.raw(),this->_error.raw());
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::sweep(const Sweeper& sweeper) {
    sweeper.sweep(this->_expansion.raw(),this->_error.raw());
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::cleanup() {
    this->sort();
    this->unique();
    this->sweep();
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::clobber() {
    this->_error=0u;
    return *this;
}



///////////////////////////////////////////////////////////////////////////////

// Accuracy control

template<class F> Float64 TaylorModel<ValidatedTag,F>::tolerance() const {
    const ThresholdSweeper* ptr=dynamic_cast<const ThresholdSweeper*>(&static_cast<const SweeperInterface&>(this->_sweeper));
    if(ptr) {
        return ptr->sweep_threshold();
    } else {
        return std::numeric_limits<double>::epsilon();
    }
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

template<class F> Box<UnitIntervalType> TaylorModel<ValidatedTag,F>::domain() const
{
    return Box<UnitIntervalType>(this->argument_size(),UnitIntervalType());
}

template<class F> ExactIntervalType TaylorModel<ValidatedTag,F>::codomain() const
{
    UpperIntervalType rng=this->range();
    return ExactIntervalType(rng.lower().raw(),rng.upper().raw());
}

// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
template<class F> UpperIntervalType TaylorModel<ValidatedTag,F>::range() const {
    const TaylorModel<ValidatedTag,F>& tm=*this;
    const SizeType as=tm.argument_size();
    ExactFloat64 constant_term=0;
    Array<ExactFloat64> linear_terms(as,0);
    Array<ExactFloat64> quadratic_terms(as,0);
    ErrorFloat64 err=0u;
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
    BoundedFloat64 r(constant_term-err,constant_term+err);
    const BoundedFloat64 unit_ivl(-1,+1);
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(SizeType j=0; j!=as; ++j) {
        const ExactFloat64& a=quadratic_terms[j];
        const ExactFloat64& b=linear_terms[j];
        BoundedFloat64 ql=abs(a)*unit_ivl + abs(b)*unit_ivl;
        if(a!=0) { // Explicitly test for zero
            BoundedFloat64 qf=a*(sqr(unit_ivl+b/a/2))-sqr(b)/a/4;
            r += refinement(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
        } else {
            r += ql;
        }
    }
    return UpperIntervalType(r);
}


template<class F> UpperIntervalType TaylorModel<ValidatedTag,F>::gradient_range(SizeType j) const {
    SizeType as=this->argument_size();
    BoundedFloat64 g(0,0);
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex const& a=iter->key();
        const Nat c=a[j];
        if(c>0) {
            const ExactFloat64& x=iter->data();
            if(a.degree()==1) { g+=x; }
            else { g+=BoundedFloat64(-1,1)*x*c; }
        }
    }
    return UpperIntervalType(g);
}

template<class F> Covector<UpperIntervalType> TaylorModel<ValidatedTag,F>::gradient_range() const {
    SizeType as=this->argument_size();
    Covector<BoundedFloat64> g(this->argument_size(),BoundedFloat64(0,0));
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        MultiIndex const& a=iter->key();
        const ExactFloat64& x=iter->data();
        for(SizeType j=0; j!=this->argument_size(); ++j) {
            const Nat c=a[j];
            if(c>0) {
                if(a.degree()==1) { g[j]+=x; }
                else { g[j]+=BoundedFloat64(-1,1)*x*c; }
            }
        }
    }
    return Covector<UpperIntervalType>(g);
}





//////////////////////////////////////////////////////////////////////////////

// ExactTag functions (max, min, abs, neg) and arithmetical functions (sqr, pow)


template<class F> TaylorModel<ValidatedTag,F> max(const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y) {
    UpperIntervalType xr=x.range();
    UpperIntervalType yr=y.range();
    if(definitely(xr.lower()>=yr.upper())) {
        return x;
    } else if(definitely(yr.lower()>=xr.upper())) {
        return y;
    } else {
        return ((x+y)+abs(x-y))/ExactFloat64(2);
    }
}


template<class F> TaylorModel<ValidatedTag,F> min(const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y) {
    UpperIntervalType xr=x.range();
    UpperIntervalType yr=y.range();
    if(definitely(xr.upper()<=yr.lower())) {
        return x;
    } else if(definitely(yr.upper()<=xr.lower())) {
        return y;
    } else {
        return ((x+y)-abs(x-y))/ExactFloat64(2);
    }
}

template<class F> TaylorModel<ValidatedTag,F> abs(const TaylorModel<ValidatedTag,F>& x) {
    UpperIntervalType xr=x.range();
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
        TaylorModel<ValidatedTag,F> r(x.argument_size(),x.sweeper());
        ExactFloat64 xmag=cast_exact(mag(xr));
        TaylorModel<ValidatedTag,F> s=x/xmag;
        s=sqr(s);
        r=static_cast<ExactFloat64>(p[n-1]);
        for(Nat i=0; i!=(n-1); ++i) {
            Nat j=(n-2)-i;
            r=s*r+static_cast<ExactFloat64>(p[j]);
        }
        r+=BoundedFloat64(-err,+err);
        return r*xmag;
    }
}

template<class F> TaylorModel<ValidatedTag,F> neg(const TaylorModel<ValidatedTag,F>& x) {
    TaylorModel<ValidatedTag,F> r(x.argument_size(),x.sweeper());
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r._append(xiter->key(),-xiter->data());
    }
    r.error()=x.error();
    return r;
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

template<class F> TaylorModel<ValidatedTag,F> sqr(const TaylorModel<ValidatedTag,F>& x) {
    return x*x;
}

template<class F> TaylorModel<ValidatedTag,F> pow(const TaylorModel<ValidatedTag,F>& x, Int n) {
    TaylorModel<ValidatedTag,F> r=x.create_constant(1);
    TaylorModel<ValidatedTag,F> p(x);
    while(n) { if(n%2) { r=r*p; } p=sqr(p); n/=2; }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
template<class X> class TaylorSeries;


template<class F> TaylorModel<ValidatedTag,F>
compose(const TaylorSeries<BoundedFloat64>& ts, const TaylorModel<ValidatedTag,F>& tv)
{
    Sweeper threshold_sweeper(new ThresholdSweeper(MACHINE_EPSILON));
    ExactFloat64& vref=const_cast<ExactFloat64&>(tv.value());
    ExactFloat64 vtmp=vref;
    vref=0;
    TaylorModel<ValidatedTag,F> r(tv.argument_size(),tv.sweeper());
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
template<class F> TaylorModel<ValidatedTag,F>
compose(const AnalyticFunction& fn, const TaylorModel<ValidatedTag,F>& tm) {
    std::cerr<<"\n\nfn="<<fn<<"\ntm="<<tm<<"\n";

    static const DegreeType MAX_DEGREE=20;
    static const Float64 MAX_TRUNCATION_ERROR=MACHINE_EPSILON;
    SweeperInterface const& sweeper=tm.sweeper();

    Float64 max_truncation_error=MAX_TRUNCATION_ERROR;
    ThresholdSweeper const* threshold_sweeper_ptr = dynamic_cast<ThresholdSweeper const*>(&sweeper);
    if(threshold_sweeper_ptr) { max_truncation_error=threshold_sweeper_ptr->sweep_threshold(); }

    DegreeType max_degree=MAX_DEGREE;
    GradedSweeper const* graded_sweeper_ptr = dynamic_cast<GradedSweeper const*>(&sweeper);
    if(graded_sweeper_ptr) { max_degree=graded_sweeper_ptr->degree(); }

    std::cerr<<"max_truncation_error="<<max_truncation_error<<"\nmax_degree="<<(uint)max_degree<<"\n";

    Nat d=max_degree;
    ExactFloat64 c=tm.value();
    BoundedFloat64 r=cast_singleton(tm.range());
    Series<BoundedFloat64> centre_series=fn.series(c);
    Series<BoundedFloat64> range_series=fn.series(r);
    std::cerr<<"c="<<c<<"\nr="<<r<<"\n";
    std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";


    ErrorFloat64 se=mag(range_series[d]-centre_series[d]);
    ErrorFloat64 e=mag(r-c);
    ErrorFloat64 p=pow(e,d);
    std::cerr<<"se="<<se<<"\ne="<<e<<"\np="<<p<<"\n";
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    ErrorFloat64 truncation_error=se*p;
    std::cerr<<"truncation_error="<<truncation_error<<"\n\n";
    if(truncation_error.raw()>max_truncation_error) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<max_truncation_error<<"\n");
    }

    TaylorModel<ValidatedTag,F> x=tm-c;
    TaylorModel<ValidatedTag,F> res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        // Don't sweep here...
    }
    res+=BoundedFloat64(-truncation_error,+truncation_error);
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


template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_embed_error(const TaylorModel<ValidatedTag,F>& tm) {
    const SizeType as=tm.argument_size();
    TaylorModel<ValidatedTag,F> rtm(as+1u,tm.sweeper());
    MultiIndex ra(as+1u);

    // The new error term is first in reverse lexicographic order.
    ExactFloat64 err_coef=cast_exact(tm.error());
    ra[as]=1;
    rtm._append(ra,err_coef);
    ra[as]=0;

    // Copy new terms
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
        MultiIndex const& xa=iter->key();
        ExactFloat64 const& xv=iter->data();
        for(SizeType j=0; j!=as; ++j) { ra[j]=xa[j]; }
        rtm._append(ra,xv);
    }
    return std::move(rtm);
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_discard_variables(const TaylorModel<ValidatedTag,F>& tm, Array<SizeType> const& discarded_variables) {
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
    TaylorModel<ValidatedTag,F> rtm(number_of_kept_variables,tm.sweeper());
    rtm.expansion().reserve(tm.number_of_nonzeros()+1u);

    // Set the uniform error of the original model
    // If index_of_error == number_of_error_variables, then the error is kept as a uniform error bound
    MultiIndex ra(number_of_kept_variables);
    ErrorFloat64 derr=0u; // Magnitude of discarded terms
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        MultiIndex const& xa=iter->key();
        ExactFloat64 const& xv=iter->data();
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


template<class F> Void TaylorModel<ValidatedTag,F>::antidifferentiate(SizeType k) {
    TaylorModel<ValidatedTag,F>& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());

    ErrorFloat64 e=0u;
    for(typename TaylorModel<ValidatedTag,F>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex& xa=xiter->key();
        ExactFloat64& xv=xiter->data();
        xa[k]+=1;
        Nat c=xa[k];
        xv=div_err(xv,c,e);
    }
    x.error()+=e;
}

template<class F> TaylorModel<ValidatedTag,F> antiderivative(const TaylorModel<ValidatedTag,F>& x, SizeType k) {
    TaylorModel<ValidatedTag,F> r(x);
    r.antidifferentiate(k);
    return r;
}


// Compute derivative inplace by computing term-by-term, switching the rounding mode
// Note that since some terms may be eliminated, requiring two iterators.
template<class F> Void TaylorModel<ValidatedTag,F>::differentiate(SizeType k) {
    TaylorModel<ValidatedTag,F> const& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());
    // ARIADNE_PRECONDITION_MSG(x.error().raw()==0,x);
    this->clobber();

    TaylorModel<ValidatedTag,F>& r=*this;
    ErrorFloat64& re=r.error();
    typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin();
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex const& xa=xiter->key();
        ExactFloat64 const& xv=xiter->data();
        Nat c=xa[k];
        if(c!=0) {
            MultiIndex& ra=riter->key();
            ExactFloat64& rv=riter->data();
            ra=xa; ra[k]-=1;
            rv=mul_err(xv,c,re);
            ++riter;
        }
    }

    r.expansion().resize(riter - r.begin());
}




template<class F> TaylorModel<ValidatedTag,F> derivative(const TaylorModel<ValidatedTag,F>& x, SizeType k) {
    TaylorModel<ValidatedTag,F> rx=x; rx.differentiate(k); return rx;

    ARIADNE_ASSERT(k<x.argument_size());

    MultiIndex ra(x.argument_size()); ExactFloat64 rv; Nat c;

    TaylorModel<ValidatedTag,F> r(x.argument_size(),x.sweeper());
    ErrorFloat64& re=r.error();
    for(typename TaylorModel<ValidatedTag,F>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex const& xa=xiter->key();
        ExactFloat64 const& xv=xiter->data();
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

template<class F> ValidatedNumericType TaylorModel<ValidatedTag,F>::_evaluate(const TaylorModel<ValidatedTag,F>& tm, const Vector<ValidatedNumericType>& x) {
    return horner_evaluate(tm.expansion(),x)+ValidatedNumericType(-tm.error(),+tm.error());
}

template<class F> ApproximateNumericType TaylorModel<ValidatedTag,F>::_evaluate(const TaylorModel<ValidatedTag,F>& tm, const Vector<ApproximateNumericType>& x) {
    return horner_evaluate(tm.expansion(),x);
}


template<class F> Covector<ValidatedNumericType> TaylorModel<ValidatedTag,F>::_gradient(const TaylorModel<ValidatedTag,F>& tm, const Vector<ValidatedNumericType>& x) {
    Vector< Differential<ValidatedNumericType> > dx=Differential<ValidatedNumericType>::variables(1u,x);
    Differential<ValidatedNumericType> df=horner_evaluate(tm.expansion(),dx)+ValidatedNumericType(-tm.error(),+tm.error());
    return gradient(df);
}



template<class F> TaylorModel<ValidatedTag,F>
TaylorModel<ValidatedTag,F>::_compose(TaylorModel<ValidatedTag,F> const& x, Vector<TaylorModel<ValidatedTag,F>> const& y) {
    return horner_evaluate(x.expansion(),y)+BoundedFloat64(-x.error(),+x.error());
}

template<class F> Void TaylorModel<ValidatedTag,F>::unscale(ExactIntervalType const& ivl) {
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
    TaylorModel<ValidatedTag,F>& tm=*this;
    ARIADNE_ASSERT_MSG(ivl.lower()<=ivl.upper(),"Cannot unscale TaylorModel<ValidatedTag,F> "<<tm<<" from empty interval "<<ivl);

    if(ivl.lower()==ivl.upper()) {
        tm=ivl.midpoint();
        // Uncomment out line below to make unscaling to a singleton interval undefined
        //tm.clear(); tm.set_error(+inf);
    } else {
        ValidatedNumericType c=ivl.centre();
        ValidatedNumericType s=rec(ivl.radius());
        tm-=c;
        tm*=s;
    }
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_compose(const Unscaling& u, const TaylorModel<ValidatedTag,F>& y) {
    TaylorModel<ValidatedTag,F> r=y; r.unscale(u.domain()); return std::move(r);
}

template<class F> TaylorModel<ValidatedTag,F>
TaylorModel<ValidatedTag,F>::_compose(const TaylorModel<ValidatedTag,F>& x, const VectorUnscaling& u, const Vector<TaylorModel<ValidatedTag,F>>& y) {
    return compose(x,compose(u,y));
}



template<class T> class Powers {
  public:
    explicit Powers(const T& t) { _values.push_back(t*0+1); _values.push_back(t); }
    explicit Powers(const T& z, const T& t) { _values.push_back(z); _values.push_back(t); }
    const T& operator[](SizeType i) const { while(_values.size()<=i) { _values.push_back(_values[1]*_values.back()); } return _values[i]; }
  private:
    mutable std::vector<T> _values;
};

template<> class Powers<ValidatedNumericType> {
  public:
    explicit Powers(const ValidatedNumericType& t) { _values.push_back(ValidatedNumericType(1)); _values.push_back(t); }
    explicit Powers(const ValidatedNumericType& z, const ValidatedNumericType& t) { _values.push_back(z); _values.push_back(t); }
    const ValidatedNumericType& operator[](SizeType i) const {
        while(_values.size()<=i) {
            if(_values.size()%2==0) { _values.push_back(sqr(_values[_values.size()/2])); }
            else { _values.push_back(_values[1]*_values.back()); } }
        return _values[i]; }
  private:
    mutable std::vector<ValidatedNumericType> _values;
};



template<class F> TaylorModel<ValidatedTag,F>
TaylorModel<ValidatedTag,F>::_partial_evaluate(const TaylorModel<ValidatedTag,F>& x, SizeType k, ValidatedNumericType c)
{
    const SizeType as=x.argument_size();
    Vector<TaylorModel<ValidatedTag,F>> y(as,TaylorModel<ValidatedTag,F>(as-1,x.sweeper()));
    for(SizeType i=0; i!=k; ++i) { y[i]=TaylorModel<ValidatedTag,F>::coordinate(as-1,i,x.sweeper()); }
    y[k]=TaylorModel<ValidatedTag,F>::constant(as-1,c,x.sweeper());
    for(SizeType i=k+1; i!=as; ++i) { y[i]=TaylorModel<ValidatedTag,F>::coordinate(as-1,i-1,x.sweeper()); }
    return compose(x,y);
}



template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_embed(SizeType as1, const TaylorModel<ValidatedTag,F>& tm2, SizeType as3) {
    return TaylorModel<ValidatedTag,F>(embed(as1,tm2.expansion(),as3),tm2.error(),tm2.sweeper());
}


template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_split(const TaylorModel<ValidatedTag,F>& tm, SizeType k, SplitPart h) {
    const DegreeType deg=tm.degree();
    const SizeType as=tm.argument_size();
    Sweeper swp=tm.sweeper();

    TaylorModel<ValidatedTag,F> r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(typename TaylorModel<ValidatedTag,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const uchar ak=iter->key()[k];
        ExactFloat64& c=iter->data();
        c/=two_exp(ak);
    }

    if(h==SplitPart::MIDDLE) { return r; }
    Int tr=( h==SplitPart::UPPER ? +1 : -1 );

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<TaylorModel<ValidatedTag,F>> ary(deg+1,TaylorModel<ValidatedTag,F>(as,swp));
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->key();
        const ExactFloat64& c=iter->data();
        uchar ak=a[k];
        a[k]=0u;
        ary[ak].expansion().append(a,c);
    }

    ErrorFloat64 re=r.error();
    r.clear();
    r.set_error(re);

    for(DegreeType i=0; i<=deg; ++i) {
        for(DegreeType j=i; j<=deg; ++j) {
            Int sf=bin(j,i);
            if(tr==-1 && (j-i)%2==1) { sf=-sf; }
            r+=ary[j]*sf;
            for(typename TaylorModel<ValidatedTag,F>::Iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->key()[k];
            }
         }
    }

    return r;
}



///////////////////////////////////////////////////////////////////////////////

// Banach algebra operations


template<class F> typename TaylorModel<ValidatedTag,F>::CoefficientType TaylorModel<ValidatedTag,F>::average() const {
    return (*this)[MultiIndex::zero(this->argument_size())];
}

template<class F> typename TaylorModel<ValidatedTag,F>::NormType TaylorModel<ValidatedTag,F>::radius() const {
    ErrorFloat64 r=0u;
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=mag(iter->data());
        }
    }
    r+=this->error();
    return r;
}

template<class F> typename TaylorModel<ValidatedTag,F>::NormType TaylorModel<ValidatedTag,F>::norm() const {
    ErrorFloat64 r=0u;
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=mag(iter->data());
    }
    r+=this->error();
    return r;
}

template<class F> typename TaylorModel<ValidatedTag,F>::NormType norm(const TaylorModel<ValidatedTag,F>& tm) {
    return tm.norm();
}


template<class F> Bool TaylorModel<ValidatedTag,F>::_refines(const TaylorModel<ValidatedTag,F>& tm1, const TaylorModel<ValidatedTag,F>& tm2)
{
    ARIADNE_ASSERT(tm1.argument_size()==tm2.argument_size());
    TaylorModel<ValidatedTag,F> d=tm2;
    d.error()=0u;
    d-=tm1;
    return d.norm().raw() <= tm2.error().raw();
}


template<class F> Bool TaylorModel<ValidatedTag,F>::_consistent(const TaylorModel<ValidatedTag,F>& tm1, const TaylorModel<ValidatedTag,F>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (Ariadne::norm(tm1-tm2).raw() <= (tm1.error()+tm2.error()).raw()*2u);
}

template<class F> Bool TaylorModel<ValidatedTag,F>::_inconsistent(const TaylorModel<ValidatedTag,F>& tm1, const TaylorModel<ValidatedTag,F>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (Ariadne::mag(tm1.value()-tm2.value()).raw() > (tm1.error()+tm2.error()).raw());
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_refinement(const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y) {
    TaylorModel<ValidatedTag,F> r(x.argument_size(),x.sweeper());

    ErrorFloat64 max_error=0u;

    const ErrorFloat64& xe=x.error();
    const ErrorFloat64& ye=y.error();
    ExactFloat64 rv,xv,yv;
    Float64 xu,yu,mxl,myl,u,ml;
    MultiIndex a;

    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
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

        BoundedFloat64 rve=refinement( xv.pm(xe), yv.pm(ye) );
        if(rve.error().raw()<0.0) {
            ARIADNE_THROW(IntersectionException,"refinement(TaylorModel<ValidatedTag,F>,TaylorModel<ValidatedTag,F>)",x<<" and "<<y<<" are inconsistent.");
        }

        if(rve.value()!=0) { r.expansion().append(a,rve.value()); }
        max_error=max(max_error,rve.error());
    }

    r.error()=max_error;

    return r;
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

template<class F> OutputStream& TaylorModel<ValidatedTag,F>::str(OutputStream& os) const {
    TaylorModel<ValidatedTag,F> const& tm=*this;

    // Set the variable names to be 'parameter' s0,s1,..
    Array<StringType> variable_names(tm.argument_size());
    for(SizeType j=0; j!=tm.argument_size(); ++j) {
        StringStream sstr;
        sstr << 's' << j;
        variable_names[j]=sstr.str();
    }

    //os << "TaylorModel<ValidatedTag,F>";
    os << "TM["<<tm.argument_size()<<"](";
    Expansion<ExactFloat64> e=tm.expansion();
    e.graded_sort();
    e.write(os,variable_names);
    return os << "+/-" << tm.error() << ")";
}


///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::zeros(SizeType rs, SizeType as, Sweeper swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(rs,TaylorModel<ValidatedTag,F>::zero(as,swp));
    return result;
}


template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::constants(SizeType as, const Vector<ValidatedNumericType>& c, Sweeper swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(c.size(),TaylorModel<ValidatedTag,F>::zero(as,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel<ValidatedTag,F>::constant(as,c[i],swp);
    }
    return result;
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::coordinates(SizeType as, Sweeper swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(as,TaylorModel<ValidatedTag,F>::zero(as,swp));
    for(SizeType i=0; i!=as; ++i) { result[i]=TaylorModel<ValidatedTag,F>::coordinate(as,i,swp); }
    return result;
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::scalings(const Vector<ExactIntervalType>& d, Sweeper swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(d.size(),TaylorModel<ValidatedTag,F>::zero(d.size(),swp));
    for(SizeType i=0; i!=d.size(); ++i) {
        result[i]=TaylorModel<ValidatedTag,F>::scaling(d.size(),i,d[i],swp);
    }
    return result;
}



///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
template<class F> Matrix<ValidatedNumericType>
jacobian(const Vector<TaylorModel<ValidatedTag,F>>& f, const Vector<ValidatedNumericType>& x) {
    Vector< Differential<ValidatedNumericType> > dx=Differential<ValidatedNumericType>::variables(1u,x);
    Vector< Differential<ValidatedNumericType> > df(f.size(),x.size(),1u);
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumericType> J=jacobian(df);
    return J;
}

// Compute the Jacobian over an arbitrary domain
template<class F> Matrix<ValidatedNumericType>
jacobian(const Vector<TaylorModel<ValidatedTag,F>>& f, const Vector<ValidatedNumericType>& x, const Array<SizeType>& p) {
    Vector<Differential<ValidatedNumericType>> dx(x.size());
    for(SizeType j=0; j!=x.size(); ++j) {
        dx[j]=Differential<ValidatedNumericType>::constant(p.size(),1u,x[j]); }
    for(SizeType k=0; k!=p.size(); ++k) {
        SizeType j=p[k];
        dx[j]=Differential<ValidatedNumericType>::variable(p.size(),1u,x[j],k); }
    Vector< Differential<ValidatedNumericType> > df(f.size());
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumericType> J=jacobian(df);
    return J;
}

// Compute the Jacobian at the origin
template<class F> Matrix<ExactFloat64>
jacobian_value(const Vector<TaylorModel<ValidatedTag,F>>& f) {
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<ExactFloat64> J(rs,as);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType j=0; j!=as; ++j) {
            a[j]=1; const ExactFloat64 x=f[i][a]; J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian at the origin with respect to the variables args.
template<class F> Matrix<ExactFloat64>
jacobian_value(const Vector<TaylorModel<ValidatedTag,F>>& f, const Array<SizeType>& p) {
    const SizeType rs=f.size();
    const SizeType as=f.zero_element().argument_size();
    const SizeType ps=p.size();
    Matrix<ExactFloat64> J(rs,ps);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType k=0; k!=ps; ++k) {
            SizeType j=p[k]; a[j]=1; const ExactFloat64 x=f[i][a]; J[i][k]=x; a[j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
template<class F> Matrix<UpperIntervalType>
jacobian_range(const Vector<TaylorModel<ValidatedTag,F>>& f) {
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<BoundedFloat64> J(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            MultiIndex const& a=iter->key();
            for(SizeType k=0; k!=as; ++k) {
                const Nat c=a[k];
                if(c>0) {
                    const ExactFloat64& x=iter->data();
                    if(a.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=BoundedFloat64(-1,1)*x*c; }
                }
            }
        }
    }
    return Matrix<UpperIntervalType>(J);
}

// Compute the Jacobian over the unit domain, with respect to the variables p.
template<class F> Matrix<UpperIntervalType>
jacobian_range(const Vector<TaylorModel<ValidatedTag,F>>& f, const Array<SizeType>& p) {
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    SizeType ps=p.size();
    Matrix<BoundedFloat64> J(rs,ps);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            MultiIndex const& a=iter->key();
            for(SizeType k=0; k!=ps; ++k) {
                SizeType j=p[k];
                const Nat c=a[j];
                if(c>0) {
                    const ExactFloat64& x=iter->data();
                    if(a.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=BoundedFloat64(-1,1)*x*c; }
                }
            }
        }
    }
    return Matrix<UpperIntervalType>(J);
}











template<class F> TaylorModel<ApproximateTag,F>::TaylorModel(SizeType as)
    : _expansion(as), _sweeper()
{
}

template<class F> TaylorModel<ApproximateTag,F>::TaylorModel(SizeType as, Sweeper swp)
    : _expansion(as), _sweeper(swp)
{
}

template<class F> TaylorModel<ApproximateTag,F> TaylorModel<ApproximateTag,F>::create_constant(ApproximateFloat64 c) const {
    TaylorModel<ApproximateTag,F> r(this->argument_size(),this->_sweeper);
    r._expansion.append(MultiIndex::zero(this->argument_size()),c);
}

template<class F> TaylorModel<ApproximateTag,F> TaylorModel<ApproximateTag,F>::create_ball(ErrorType) const {
    return TaylorModel<ApproximateTag,F>(this->argument_size(),this->_sweeper);
}


template<class F> Void TaylorModel<ApproximateTag,F>::iadd(const ApproximateFloat64& c)
{
    // Compute self+=c
    TaylorModel<ApproximateTag,F>& x=*this;
    if(decide(c==0)) { return; }
    if(x._expansion.empty()) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else if((x._expansion.end()-1)->key().degree()>0) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else {
        ApproximateFloat64& rv=(x._expansion.end()-1)->data();
        rv+=c;
    }
}

template<class F> Void TaylorModel<ApproximateTag,F>::imul(const ApproximateFloat64& c)
{
    // Compute self*=c
    if(decide(c==0)) { this->clear(); return; }
    if(decide(c==1)) { return; }
    for(ExpansionType::Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        iter->data() *= c;
    }
}



template<class F> Void TaylorModel<ApproximateTag,F>::isma(const ApproximateFloat64& c, const TaylorModel<ApproximateTag,F>& y)
{
    const TaylorModel<ApproximateTag,F>& x=*this;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<", y="<<y);
    TaylorModel<ApproximateTag,F> r(x.argument_size());

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


template<class F> Void TaylorModel<ApproximateTag,F>::ifma(const TaylorModel<ApproximateTag,F>& x, const TaylorModel<ApproximateTag,F>& y)
{
    TaylorModel<ApproximateTag,F>& r(*this);
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<",y="<<y);

    TaylorModel<ApproximateTag,F> t(x.argument_size());
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


template<class F> Void TaylorModel<ApproximateTag,F>::unscale(ExactIntervalType const& dom) {
    TaylorModel<ApproximateTag,F>& x=*this;
    x-=dom.midpoint();
    x/=dom.radius();
}



template<class F> typename TaylorModel<ApproximateTag,F>::NormType TaylorModel<ApproximateTag,F>::norm() const {
    ApproximateFloat64 r=0.0;
    for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        r+=mag(iter->data());
    }
    return NormType(r);
}

template<class F> typename TaylorModel<ApproximateTag,F>::CoefficientType TaylorModel<ApproximateTag,F>::average() const {
    return (*this)[MultiIndex(this->argument_size())];
}

template<class F> typename TaylorModel<ApproximateTag,F>::NormType TaylorModel<ApproximateTag,F>::radius() const {
    ApproximateFloat64 r=0.0;
    for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        if(iter->key().degree()!=0) {
            r+=mag(iter->data());
        }
    }
    return NormType(r);
}

template<class F> UnitBoxType TaylorModel<ApproximateTag,F>::domain() const {
    return UnitBoxType(this->argument_size(),UnitIntervalType());
}

template<class F> ApproximateIntervalType TaylorModel<ApproximateTag,F>::range() const {
    ApproximateFloat64 av=this->average();
    ApproximateFloat64 rad=this->radius();
    return ApproximateIntervalType(av-rad,av+rad);
}

template<class F> TaylorModel<ApproximateTag,F>& TaylorModel<ApproximateTag,F>::sweep() {
    this->_sweeper.sweep(this->_expansion.raw());
    return *this;
}

template<class F> Float64 TaylorModel<ApproximateTag,F>::tolerance() const {
    return dynamic_cast<const ThresholdSweeper&>(static_cast<const SweeperInterface&>(this->_sweeper)).sweep_threshold();
}


template<class F> OutputStream& TaylorModel<ApproximateTag,F>::write(OutputStream& os) const {
    return os << "TM["<<this->argument_size()<<"](" << this->_expansion << ")";
}

template<class F> OutputStream& TaylorModel<ApproximateTag,F>::str(OutputStream& os) const {
    return os << this->_expansion;
}

} //namespace Ariadne


