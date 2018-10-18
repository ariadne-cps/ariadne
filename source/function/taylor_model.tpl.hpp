/***************************************************************************
 *            taylor_model.tpl.hpp
 *
 *  Copyright 2008-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../numeric/numeric.hpp"

#include <iomanip>
#include <limits>

#include "../numeric/rounding.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/covector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/expansion.hpp"
#include "../algebra/series.hpp"
#include "../algebra/differential.hpp"
#include "../function/taylor_model.hpp"
#include "../function/taylor_series.hpp"
#include "../function/function.hpp"
#include "../utility/exceptions.hpp"

#include "../algebra/expansion.inl.hpp"
#include "../algebra/evaluate.tpl.hpp"
#include "../algebra/algebra_operations.tpl.hpp"

#define VOLATILE ;
#include "../algebra/multi_index-noaliasing.hpp"
#include "../function/function_mixin.hpp"
#include "../algebra/vector.hpp"

namespace Ariadne {

FloatDP operator+(FloatDP x1, FloatDP x2);
FloatDP operator-(FloatDP x1, FloatDP x2);
FloatDP operator*(FloatDP x1, FloatDP x2);
FloatDP operator/(FloatDP x1, FloatDP x2);
FloatDP& operator+=(FloatDP& x1, FloatDP x2);
FloatDP& operator-=(FloatDP& x1, FloatDP x2);
FloatDP& operator*=(FloatDP& x1, FloatDP x2);
FloatDP& operator/=(FloatDP& x1, FloatDP x2);
FloatMP operator+(FloatMP const& x1, FloatMP const& x2);
FloatMP operator-(FloatMP const& x1, FloatMP const& x2);
FloatMP operator*(FloatMP const& x1, FloatMP const& x2);
FloatMP operator/(FloatMP const& x1, FloatMP const& x2);
FloatMP& operator+=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator-=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator*=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator/=(FloatMP& x1, FloatMP const& x2);

namespace {

static const double MACHINE_EPSILON = 2.2204460492503131e-16;

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }


Interval<FloatDPValue> const& convert_interval(IntervalDomainType const& ivl, DoublePrecision) { return ivl; }
Interval<FloatMPValue> convert_interval(IntervalDomainType const& ivl, MultiplePrecision pr) {
    return Interval<FloatMPValue>(FloatMP(ivl.lower().get_d(),pr),FloatMP(ivl.upper().get_d(),pr)); }

Interval<FloatDPValue> convert_exact_interval(Interval<FloatDPUpperBound> const& ivl, DoublePrecision pr) {
    return cast_exact(ivl); }
Interval<FloatDPValue> convert_exact_interval(Interval<FloatMPUpperBound> const& ivl, DoublePrecision pr) {
    FloatDP l(Dyadic(ivl.lower().raw()),downward,pr); FloatDP u(Dyadic(ivl.lower().raw()),upward,pr);
    return Interval<FloatDPValue>(FloatDPValue(l),FloatDPValue(u)); }

inline Box<Interval<FloatDPValue>> const& convert_box(BoxDomainType const& bx, DoublePrecision) { return bx; }

} // namespace



template<class F>
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,FloatValue<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatValue<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatValue<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatValue<PR>>::Iterator curr=p.begin();

    // FIXME: Not needed, but added to pair with rounding mode change below
    F::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(adv->index(),adv->coefficient().raw())) {
            //te+=abs(adv->coefficient());
            te+=cast_positive(abs(adv->coefficient()));
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    // FIXME: Removing reset of rounding mode causes error in TestNonlinear programming
    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class F>
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,FloatApproximation<PR>>& p) const
{
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::Iterator curr=p.begin();
    while(adv!=end) {
        if(this->_discard(adv->index(),adv->coefficient().raw())) {
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    p.resize(static_cast<SizeType>(curr-p.begin()));
}


namespace {

template<class I, class X> decltype(auto) radius(Expansion<I,X> const& p) {
    typedef decltype(mag(declval<X>())) R;
    R r=mag(p.zero_coefficient());
    for (auto term : p) { if (term.index().degree() != 0) { r += mag(term.coefficient()); } }
    return r;
}


} //namespace


template<class F>
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,FloatValue<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatValue<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatValue<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatValue<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=radius(p)+e;

    // FIXME: Not needed, but added to pair with rounding mode change below
    F::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(adv->coefficient().raw(),nrm.raw())) {
            //te+=abs(adv->coefficient());
            te+=cast_positive(abs(adv->coefficient()));
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;

    // FIXME: Removing reset of rounding mode causes error in TestNonlinear programming
    F::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}


template<class F>
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,FloatApproximation<PR>>& p) const
{
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatApproximation<PR>>::Iterator curr=p.begin();

    FloatApproximation<PR> nrm=radius(p);

    while(adv!=end) {
        if(this->_discard(adv->coefficient().raw(),nrm.raw())) {
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    p.resize(static_cast<SizeType>(curr-p.begin()));
}



template<class F> TaylorModel<ValidatedTag,F>::TaylorModel()
    : _expansion(0), _error(0u), _sweeper()
{
}


template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(SizeType as, SweeperType swp)
    : _expansion(as), _error(0u), _sweeper(swp)
{
}

template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(const Expansion<MultiIndex,CoefficientType>& f, const ErrorType& e, SweeperType swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->cleanup();
}

template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(const Expansion<MultiIndex,F>& f, const F& e, SweeperType swp)
//    : TaylorModel(reinterpret_cast<Expansion<MultiIndex,CoefficientType>const&>(f),reinterpret_cast<ErrorType const&>(e),swp)
    : TaylorModel(static_cast<Expansion<MultiIndex,CoefficientType>>(f),static_cast<ErrorType>(e),swp)
{
}

template<class F> TaylorModel<ValidatedTag,F>::TaylorModel(const Expansion<MultiIndex,double>& f, const double& e, SweeperType swp)
    : TaylorModel(Expansion<MultiIndex,CoefficientType>(Expansion<MultiIndex,Dyadic>(f),swp.precision()),ErrorType(Dyadic(e),swp.precision()),swp)
{
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::scaling(SizeType as, SizeType j, const IntervalDomainType& codom, SweeperType swp) {
    TaylorModel<ValidatedTag,F> r(as,swp);
    auto ivl=convert_interval(codom,r.precision());
    r.set_gradient(j,1);
    r*=ivl.radius();
    r+=ivl.midpoint();
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

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create_constant(GenericNumericType c) const {
    return this->create_constant(NumericType(c,this->precision()));
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::create_coordinate(SizeType j) const {
    ARIADNE_PRECONDITION(j<this->argument_size());
    TaylorModel<ValidatedTag,F> r(this->argument_size(),this->_sweeper);
    CoefficientType one(1,this->precision());
    r._expansion.append(MultiIndex::unit(this->argument_size(),j),one);
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
        deg=std::max(deg,iter->index().degree());
    }
    return deg;
}


template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::operator=(const NumericType& c) {
    this->_expansion.clear();
    FloatValue<PR> m=c.value();
    if(m!=0) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    this->_error=c.error();
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::operator=(const GenericNumericType& c) {
    return *this = NumericType(c,this->precision());
}

namespace { // Internal code for arithmetic

template<class F> struct ValidatedApproximation {
    Bounds<F> _v; Approximation<F> _a;
    ValidatedApproximation(Bounds<F>const& x) : _v(x), _a(x) { }
    LowerBound<F> lower() const { return _v.lower(); }
    Approximation<F> middle() const { return _a; }
    UpperBound<F> upper() const { return _v.upper(); }
    F const& lower_raw() const { return _v.lower_raw(); }
    F const& middle_raw() const { return _a.raw(); }
    F const& upper_raw() const { return _v.upper_raw(); }
};

template<class F> Value<F> add_err(Value<F> const& x1, Value<F> const& x2, Error<F>& e) {
    Value<F> mx1=-x1;
    F::set_rounding_to_nearest();
    Value<F> r(x1.raw() + x2.raw());
    F::set_rounding_upward();
    F u=x1.raw()+x2.raw();
    F ml=mx1.raw()-x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

template<class F> Value<F> add_err(Value<F> const& x, ValidatedApproximation<F> const& c, Error<F>& e) {
    F const& xv=x.raw();
    F const& cl=c.lower_raw();
    F const& cm=c.middle_raw();
    F const& cu=c.upper_raw();
    F& re=e.raw();
    F::set_rounding_to_nearest();
    F rv=xv+cm;
    F::set_rounding_upward();
    F u=xv+cu;
    F ml=(-xv)-cl;
    re += (u+ml)/2;
    return Value<F>(rv);
}

template<class F> Value<F> add_err(Value<F> const& x, Bounds<F> const& c, Error<F>& e) {
    return add_err(x,ValidatedApproximation<F>(c),e);
}

template<class F> Value<F> sub_err(Value<F> const& x1, Value<F> const& x2, Error<F>& e) {
    Value<F> mx1=-x1;
    F::set_rounding_to_nearest();
    Value<F> r(x1.raw() - x2.raw());
    F::set_rounding_upward();
    F u=x1.raw()-x2.raw();
    F ml=mx1.raw()+x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

template<class F> Value<F> mul_no_err(Value<F> const& x1, Value<F> const& x2) {
    F::set_rounding_to_nearest();
    Value<F> r(x1.raw() * x2.raw());
    F::set_rounding_upward();
    return r;
}

template<class F> Value<F> mul_err(Value<F> const& x1, Value<F> const& x2, Error<F>& e) {
    Value<F> mx1=-x1;
    F::set_rounding_to_nearest();
    Value<F> r(x1.raw() * x2.raw());
    F::set_rounding_upward();
    F u=x1.raw()*x2.raw();
    F ml=mx1.raw()*x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

template<class F> Value<F> mul_err(Value<F> const& x, ValidatedApproximation<F> const& c, Error<F>& e) {
    F const& xv=x.raw();
    F const& cu=c.upper_raw();
    F const& cm=c.middle_raw();
    F const& cl=c.lower_raw();
    F& re=e.raw();
    F::set_rounding_to_nearest();
    F rv=xv*cm;
    F::set_rounding_upward();
    F u,ml;
    if(xv>=0) {
        F mcl=-cl;
        u=xv*cu;
        ml=xv*mcl;
    } else {
        F mcu=-cu;
        u=xv*cl;
        ml=xv*mcu;
    }
    re+=(u+ml)/2;
    return Value<F>(rv);
}

template<class F> Value<F> mul_err(Value<F> const& x, Bounds<F> const& c, Error<F>& e) {
    return mul_err(x,ValidatedApproximation<F>(c),e);
}

template<class F> Value<F> fma_err(Value<F> const& x, ValidatedApproximation<F> const& c, Value<F> y, Error<F>& e) {
    F const& xv=x.raw();
    F const& cu=c.upper_raw();
    F const& cm=c.middle_raw();
    F const& cl=c.lower_raw();
    F const& yv=y.raw();
    F& re=e.raw();
    F::set_rounding_to_nearest();
    F rv=xv+cm*yv;
    F::set_rounding_upward();
    F u,ml;
    if(yv>=0) {
        F mcl=-cl;
        u=cu*yv+xv;
        ml=mcl*yv-xv;
    } else {
        F mcu=-cu;
        u=cl*yv+xv;
        ml=mcu*yv-xv;
    }
    re+=(u+ml)/2;
    return Value<F>(rv);
}

template<class F> Value<F> mul_err(Value<F> const& x1, Nat n2, Error<F>& e) {
    return mul_err(x1,Value<F>(n2,x1.precision()),e);
}

template<class F> Value<F> div_err(Value<F> const& x1, Value<F> const& x2, Error<F>& e) {
    Value<F> mx1=-x1;
    F::set_rounding_to_nearest();
    Value<F> r(x1.raw() / x2.raw());
    F::set_rounding_upward();
    F u=x1.raw()/x2.raw();
    F ml=mx1.raw()/x2.raw();
    e.raw() += (u+ml)/2;
    return r;
}

template<class F> Value<F> div_err(Value<F> const& x1, Nat n2, Error<F>& e) {
    return div_err(x1,Value<F>(n2,x1.precision()),e);
}

// Inplace negation
template<class F> Void _neg(TaylorModel<ValidatedTag,F>& r)
{
    for(auto iter=r.begin(); iter!=r.end(); ++iter) {
        iter->coefficient()=-iter->coefficient();
    }
}


template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const TwoExp& c)
{
    typedef typename F::PrecisionType PR;
    if(FloatValue<PR>(c)==1) { return; }
    for(typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->coefficient()*=c;
    }
    r.error()*=FloatError<PR>(c);
 }

template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const Value<F>& c) {
    typedef typename F::PrecisionType PR;
    FloatError<PR> e=0u; // The maximum accumulated error
    for(typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->coefficient() = mul_err(riter->coefficient(),c,e);
    }
    FloatError<PR>& re=r.error();
    re*=abs(c);
    re+=e;
}

template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const Bounds<F>& c)
{
    typedef typename F::PrecisionType PR;
    //std::cerr<<"TaylorModel<ValidatedTag,F>::scal(Float64Bounds c) c="<<c<<std::endl;
    ARIADNE_ASSERT_MSG(is_finite(c.lower().raw()) && is_finite(c.upper().raw()),"scal(tm,c): tm="<<r<<", c="<<c);
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);

    const F inf = F::inf(r.precision());
    if(r.error().raw()==inf) {
        r.expansion().clear(); return;
    }

    FloatError<PR> e=nul(r.error());
    ValidatedApproximation<F> clmu=c;
    for(typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        UniformReference<FloatValue<PR>> rv=riter->coefficient();
        rv=mul_err(rv,clmu,e);
    }
    r.error()*=mag(c);
    r.error()+=e;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


struct UnitMultiIndex { SizeType argument_size; SizeType unit_index; };


template<class F> inline Void _incr(TaylorModel<ValidatedTag,F>& r, const MultiIndex& a) {
    for(typename TaylorModel<ValidatedTag,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->index()+=a;
    }
}

template<class F> inline Void _incr(TaylorModel<ValidatedTag,F>& r, SizeType j) {
    for(typename TaylorModel<ValidatedTag,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++iter->index()[j];
    }
}


template<class F> inline Void _acc(TaylorModel<ValidatedTag,F>& r, const Value<F>& c) {
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    typedef typename F::PrecisionType PR;
    if(c==0) { return; }
    if(r.expansion().empty() || r.expansion().back().index().degree()>0) {
        r.expansion().append(MultiIndex(r.argument_size()),c);
    } else {
        UniformReference<FloatValue<PR>> rv=(r.end()-1)->coefficient();
        FloatError<PR>& re=r.error();
        rv=add_err(rv,c,re);
    }
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    return;
}



template<class F> inline Void _acc(TaylorModel<ValidatedTag,F>& r, const Bounds<F>& c)
{
    // Compute self+=c
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
    typedef typename F::PrecisionType PR;

    const F inf = F::inf(r.precision());
    if(c.lower().raw()==-inf || c.upper().raw()==+inf) {
        r.clear();
        r.set_error(mag(FloatValue<PR>(inf)));
        return;
    }

    if(c.lower().raw()==-c.upper().raw()) { // The midpoint of the interval is zero, so no need to change constant term
        r.error()+=mag(c);
        return;
    }

    if(r.expansion().empty() || r.expansion().back().index().degree()>0) { // Append a constant term zero
        r._append(MultiIndex(r.argument_size()),FloatValue<PR>(0,r.precision()));
    }

    UniformReference<FloatValue<PR>> rv=(r.end()-1)->coefficient();
    FloatError<PR>& re=r.error();
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
    typedef typename F::PrecisionType PR;
    FloatError<PR> e=nul(r.error());
    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->index()<yiter->index()) {
            r._append(xiter->index(),xiter->coefficient());
            ++xiter;
        } else if(yiter->index()<xiter->index()) {
            r._append(yiter->index(),yiter->coefficient());
            ++yiter;
        } else {
            ARIADNE_DEBUG_ASSERT(xiter->index()==yiter->index());
            r._append(xiter->index(),add_err(xiter->coefficient(),yiter->coefficient(),e));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r._append(xiter->index(),xiter->coefficient());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r._append(yiter->index(),yiter->coefficient());
        ++yiter;
    }

    FloatDP::set_rounding_upward();
    r.error()=(x.error()+y.error())+e;
    FloatDP::set_rounding_to_nearest();

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
    typedef typename F::PrecisionType PR;
    FloatError<PR> e=0u;
    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->index()<yiter->index()) {
            r._append(xiter->index(),xiter->coefficient());
            ++xiter;
        } else if(yiter->index()<xiter->index()) {
            r._append(yiter->index(),-yiter->coefficient());
            ++yiter;
        } else {
            ARIADNE_DEBUG_ASSERT(xiter->index()==yiter->index());
            r._append(xiter->index(),sub_err(xiter->coefficient(),yiter->coefficient(),e));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        r._append(xiter->index(),xiter->coefficient());
        ++xiter;
    }
    while(yiter!=y.end()) {
        r._append(yiter->index(),-yiter->coefficient());
        ++yiter;
    }

    FloatDP::set_rounding_upward();
    r.error()=(x.error()+y.error())+e;
    FloatDP::set_rounding_to_nearest();

    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


template<class F> inline Void _sma(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x, const FloatBounds<typename F::PrecisionType>& c, const TaylorModel<ValidatedTag,F>& y)
{
    typedef typename F::PrecisionType PR;
    ARIADNE_ASSERT_MSG(c.lower().raw()<=c.upper().raw(),c);
    ARIADNE_ASSERT_MSG(x.error().raw()>=0,"x="<<x);
    ARIADNE_ASSERT_MSG(y.error().raw()>=0,"y="<<y);

    VOLATILE FloatDP u,ml,myv;
    FloatError<PR> err=nul(r.error()); // Twice the maximum accumulated error
    ValidatedApproximation<F> clmu=c;

    // Compute r=x+y, assuming r is empty
    RawFloat<PR>::set_rounding_upward();
    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->index()<yiter->index()) {
            UniformConstReference<FloatValue<PR>> xv=xiter->coefficient();
            r.expansion().append(xiter->index(),xv);
            ++xiter;
        } else if(yiter->index()<xiter->index()) {
            UniformConstReference<FloatValue<PR>> yv=yiter->coefficient();
            r.expansion().append(yiter->index(),mul_err(yv,c,err));
            ++yiter;
        } else {
            UniformConstReference<FloatValue<PR>> xv=xiter->coefficient();
            UniformConstReference<FloatValue<PR>> yv=yiter->coefficient();
            r.expansion().append(xiter->index(),fma_err(xv,clmu,yv,err));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        UniformConstReference<FloatValue<PR>> xv=xiter->coefficient();
        r.expansion().append(xiter->index(),xv);
        ++xiter;
    }
    while(yiter!=y.end()) {
        UniformConstReference<FloatValue<PR>> yv=yiter->coefficient();
        r.expansion().append(yiter->index(),mul_err(yv,clmu,err));
        ++yiter;
    }

    r.error()=x.error()+y.error();
    r.error()+=err;

    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}



// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
template<class F> inline Void _mul(TaylorModel<ValidatedTag,F>& r, const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y)
{
    typedef typename F::PrecisionType PR;
    const SizeType as=r.argument_size();
    TaylorModel<ValidatedTag,F> t(as,r.sweeper());
    TaylorModel<ValidatedTag,F> s(as,r.sweeper());
    MultiIndex ta(as);
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        FloatError<PR> te=nul(r.error()); // trucation error
        FloatError<PR> re=nul(r.error()); // roundoff error
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<FloatValue<PR>> xv=xiter->coefficient();
        for(typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            UniformConstReference<MultiIndex> ya=yiter->index();
            UniformConstReference<FloatValue<PR>> yv=yiter->coefficient();
            ta=xa+ya;
            FloatValue<PR> tv=mul_no_err(xv,yv);
            // NOTE: Previously, we allowed to discard terms immediately since Sweeper() had a discard methd
            // if(r.sweeper().discard(ta,tv)) { te+=mag(xv)*mag(yv); }
            t._append(ta,tv);
            re+=(xv*yv).error();
        }
        t.error()=te+re;

        t.sweep();

        _add(s,r,t);
        r.expansion().swap(s.expansion());
        r.error()=s.error();
        s.expansion().clear();
        s.error()=0u;
        t.expansion().clear();
        t.error()=0u;
    }

    FloatError<PR> xs=nul(r.error());
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=mag(xiter->coefficient());
    }

    FloatError<PR> ys=nul(r.error());
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=mag(yiter->coefficient());
    }

    FloatError<PR>& re=r.error();
    const FloatError<PR>& xe=x.error();
    const FloatError<PR>& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    return;
}


} // namespace


///////////////////////////////////////////////////////////////////////////////

// Inplace arithmetical operations for Algebra concept

template<class F> Void TaylorModel<ValidatedTag,F>::iadd(const NumericType& c)
{
    _acc(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class F> Void TaylorModel<ValidatedTag,F>::imul(const NumericType& c)
{
    _scal(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class F> Void TaylorModel<ValidatedTag,F>::isma(const NumericType& c, const TaylorModel<ValidatedTag,F>& y)
{
    TaylorModel<ValidatedTag,F>& x=*this;
    TaylorModel<ValidatedTag,F> r=this->create();
    _sma(r,x,c,y);
    this->swap(r);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class F> Void TaylorModel<ValidatedTag,F>::ifma(const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y)
{
    _mul(*this,x,y);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class F> struct AlgebraOperations<TaylorModel<ValidatedTag,F>>
    : NormedAlgebraOperations<TaylorModel<ValidatedTag,F>>
{
    typedef TaylorModel<ValidatedTag,F> ModelType;
    typedef typename TaylorModel<ValidatedTag,F>::NumericType NumericType;
    typedef typename ModelType::RangeType RangeType;
  public:
    using NormedAlgebraOperations<TaylorModel<ValidatedTag,F>>::apply;
    static ModelType apply(Nul, ModelType const& x) {
        return ModelType(x.argument_size(),x.sweeper()); }
    static ModelType apply(Pos, ModelType x) {
        return std::move(x); }
    static ModelType apply(Neg, ModelType x) {
        x.imul(NumericType(-1)); return std::move(x); }
    static ModelType apply(Add, ModelType const& x, ModelType const& y) {
        auto r=x; r.isma(NumericType(+1),y); return std::move(r); }
    static ModelType apply(Sub, ModelType const& x, ModelType const& y) {
        auto r=x; r.isma(NumericType(-1),y); return std::move(r); }
    static ModelType apply(Mul, ModelType const& x, ModelType const& y) {
        auto r=nul(x); r.ifma(x,y); return std::move(r); }
    static ModelType apply(Add, ModelType x, NumericType const& c) {
        auto& r=x; r.iadd(c); return std::move(r); }
    static ModelType apply(Mul, ModelType x, NumericType const& c) {
        auto& r=x; r.imul(c); return std::move(r); }
    static ModelType apply(Max, ModelType const& x, ModelType const& y);
    static ModelType apply(Min, ModelType const& x, ModelType const& y);
    static ModelType apply(Abs, ModelType const& x);
};


///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::sort() {
    this->_expansion.sort();
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::unique()
{
    typename TaylorModel<ValidatedTag,F>::ConstIterator advanced =this->begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator end =this->end();
    typename TaylorModel<ValidatedTag,F>::Iterator current=this->begin();
    FloatError<PR> e=nul(this->error());
    while(advanced!=end) {
        current->index()=advanced->index();
        FloatValue<PR> rv=advanced->coefficient();
        ++advanced;
        while(advanced!=end && advanced->index()==current->index()) {
            UniformConstReference<FloatValue<PR>> xv=advanced->coefficient();
            rv=add_err(rv,xv,e);
            ++advanced;
        }
        current->coefficient()=rv;
        ++current;
    }
    this->error()+=e;
    this->_expansion.resize(static_cast<SizeType>(current-this->begin()));

    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::sweep() {
//    this->_sweeper.sweep(reinterpret_cast<Expansion<MultiIndex,F>&>(this->_expansion),reinterpret_cast<F&>(this->_error));
    this->_sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::sweep(const SweeperType& sweeper) {
    sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::simplify() {
    return this->sweep();
}

template<class F> TaylorModel<ValidatedTag,F>& TaylorModel<ValidatedTag,F>::simplify(const PropertiesType& properties) {
    return this->sweep(properties);
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

template<class F> F TaylorModel<ValidatedTag,F>::tolerance() const {
    const ThresholdSweeper<F>* ptr=dynamic_cast<const ThresholdSweeper<F>*>(&static_cast<const SweeperInterface<F>&>(this->_sweeper));
    return (ptr) ? ptr->sweep_threshold() : F(std::numeric_limits<double>::epsilon(),this->precision());
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

template<class F> UnitBox TaylorModel<ValidatedTag,F>::domain() const
{
    return UnitBox(this->argument_size(),UnitInterval());
}

template<class F> auto TaylorModel<ValidatedTag,F>::codomain() const -> CodomainType
{
    RangeType rng=this->range();
    return convert_exact_interval(rng,dp);
}


// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
template<class F> auto TaylorModel<ValidatedTag,F>::range() const -> RangeType {
    const TaylorModel<ValidatedTag,F>& tm=*this;
    const SizeType as=tm.argument_size();
    const PrecisionType prec = tm.precision();
    const FloatValue<PR> zero(prec);
    FloatValue<PR> constant_term(zero);
    Array<FloatValue<PR>> linear_terms(as,zero);
    Array<FloatValue<PR>> quadratic_terms(as,zero);
    FloatError<PR> err(prec);
    for(auto iter=tm.begin(); iter!=tm.end(); ++iter) {
        if(iter->index().degree()==0) {
            constant_term=iter->coefficient();
        } else if(iter->index().degree()==1) {
            for(SizeType j=0; j!=tm.argument_size(); ++j) {
                if(iter->index()[j]==1) { linear_terms[j]=iter->coefficient(); break; }
            }
        } else if(iter->index().degree()==2) {
            for(SizeType j=0; j!=tm.argument_size(); ++j) {
                if(iter->index()[j]==2) { quadratic_terms[j]=iter->coefficient(); break; }
                if(iter->index()[j]==1) { err+=mag(iter->coefficient()); break; }
            }
        } else {
            err+=mag(iter->coefficient());
        }
    }
    err=err+tm.error();
    FloatBounds<PR> r(constant_term-err,constant_term+err);
    const FloatBounds<PR> unit_ivl(-1,+1);
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(SizeType j=0; j!=as; ++j) {
        const FloatValue<PR>& a=quadratic_terms[j];
        const FloatValue<PR>& b=linear_terms[j];
        FloatBounds<PR> ql=abs(a)*unit_ivl + abs(b)*unit_ivl;
        if(a!=0) { // Explicitly test for zero
            FloatBounds<PR> qf=a*(sqr(unit_ivl+b/a/2))-sqr(b)/a/4;
            r += refinement(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
        } else {
            r += ql;
        }
    }
    return RangeType(r);
}


template<class F> auto TaylorModel<ValidatedTag,F>::gradient_range(SizeType j) const -> RangeType {
    FloatBounds<PR> g(0,0);
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        UniformConstReference<MultiIndex> a=iter->index();
        const Nat c=a[j];
        if(c>0) {
            UniformConstReference<FloatValue<PR>> x=iter->coefficient();
            if(a.degree()==1) { g+=x; }
            else { g+=FloatBounds<PR>(-1,1)*x*c; }
        }
    }
    return RangeType(g);
}

template<class F> auto TaylorModel<ValidatedTag,F>::gradient_range() const -> Covector<RangeType> {
    Covector<FloatBounds<PR>> g(this->argument_size(),FloatBounds<PR>(0,0));
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<FloatValue<PR>> x=iter->coefficient();
        for(SizeType j=0; j!=this->argument_size(); ++j) {
            const Nat c=a[j];
            if(c>0) {
                if(a.degree()==1) { g[j]+=x; }
                else { g[j]+=FloatBounds<PR>(-1,1)*x*c; }
            }
        }
    }
    return Covector<RangeType>(g);
}





//////////////////////////////////////////////////////////////////////////////

// ExactTag functions (max, min, abs, neg) and arithmetical functions (sqr, pow)


template<class F> TaylorModel<ValidatedTag,F> AlgebraOperations<TaylorModel<ValidatedTag,F>>::apply(Max, const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y) {
    typedef typename F::PrecisionType PR;
    RangeType xr=x.range();
    RangeType yr=y.range();
    if(definitely(xr.lower()>=yr.upper())) {
        return x;
    } else if(definitely(yr.lower()>=xr.upper())) {
        return y;
    } else {
        return ((x+y)+abs(x-y))/FloatValue<PR>(2);
    }
}


template<class F> TaylorModel<ValidatedTag,F> AlgebraOperations<TaylorModel<ValidatedTag,F>>::apply(Min, const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y) {
    typedef typename F::PrecisionType PR;
    RangeType xr=x.range();
    RangeType yr=y.range();
    if(definitely(xr.upper()<=yr.lower())) {
        return x;
    } else if(definitely(yr.upper()<=xr.lower())) {
        return y;
    } else {
        return ((x+y)-abs(x-y))/FloatValue<PR>(2);
    }
}

template<class F> TaylorModel<ValidatedTag,F> AlgebraOperations<TaylorModel<ValidatedTag,F>>::apply(Abs, const TaylorModel<ValidatedTag,F>& x) {
    typedef typename F::PrecisionType PR;
    RangeType xr=x.range();
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
        FloatValue<PR> xmag=cast_exact(mag(xr));
        TaylorModel<ValidatedTag,F> s=x/xmag;
        s=sqr(s);
        r=static_cast<FloatValue<PR>>(p[n-1u]);
        for(Nat i=0; i!=(n-1u); ++i) {
            Nat j=(n-2)-i;
            r=s*r+static_cast<FloatValue<PR>>(p[j]);
        }
        r+=FloatBounds<PR>(-err,+err);
        return r*xmag;
    }
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
compose(const TaylorSeries<typename F::PrecisionType>& ts, const TaylorModel<ValidatedTag,F>& tv)
{
    typedef typename F::PrecisionType PR;
    Sweeper<F> threshold_sweeper(new ThresholdSweeper<F>(tv.precision(),MACHINE_EPSILON));
    FloatValue<PR>& vref=const_cast<FloatValue<PR>&>(tv.value());
    FloatValue<PR> vtmp=vref;
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
    typedef typename F::PrecisionType PR;

    static const DegreeType MAX_DEGREE=20;
    static const FloatDP MAX_TRUNCATION_ERROR=MACHINE_EPSILON;
    SweeperInterface<F> const& sweeper=tm.sweeper();

    F max_truncation_error=MAX_TRUNCATION_ERROR;
    ThresholdSweeper<F> const* threshold_sweeper_ptr = dynamic_cast<ThresholdSweeper<F> const*>(&sweeper);
    if(threshold_sweeper_ptr) { max_truncation_error=threshold_sweeper_ptr->sweep_threshold(); }

    DegreeType max_degree=MAX_DEGREE;
    GradedSweeper<F> const* graded_sweeper_ptr = dynamic_cast<GradedSweeper<F> const*>(&sweeper);
    if(graded_sweeper_ptr) { max_degree=graded_sweeper_ptr->degree(); }

    //std::cerr<<"max_truncation_error="<<max_truncation_error<<"\nmax_degree="<<(uint)max_degree<<"\n";

    Nat d=max_degree;
    FloatValue<PR> c=tm.value();
    FloatDPBounds r=cast_singleton(tm.range());
    Series<FloatDPBounds> centre_series=fn.series(c);
    Series<FloatDPBounds> range_series=fn.series(r);
    //std::cerr<<"c="<<c<<"\nr="<<r<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";


    FloatError<PR> se=mag(range_series[d]-centre_series[d]);
    FloatError<PR> e=mag(r-c);
    FloatError<PR> p=pow(e,d);
    //std::cerr<<"se="<<se<<"\ne="<<e<<"\np="<<p<<"\n";
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    FloatError<PR> truncation_error=se*p;
    //std::cerr<<"truncation_error="<<truncation_error<<"\n\n";
    if(truncation_error.raw()>max_truncation_error) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<max_truncation_error<<"\n");
    }

    TaylorModel<ValidatedTag,F> x=tm-c;
    TaylorModel<ValidatedTag,F> res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(Nat i=0; i!=d; ++i) {
        res=centre_series[d-i-1u]+x*res;
        // Don't sweep here...
    }
    res+=FloatDPBounds(-truncation_error,+truncation_error);
    return res;
}



///////////////////////////////////////////////////////////////////////////////

// Algebraic and trancendental functions are implemented using generic Algebra code

///////////////////////////////////////////////////////////////////////////////

// Inplace operators manipulating the error term

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_embed_error(const TaylorModel<ValidatedTag,F>& tm) {
    const SizeType as=tm.argument_size();
    TaylorModel<ValidatedTag,F> rtm(as+1u,tm.sweeper());
    MultiIndex ra(as+1u);

    // The new error term is first in reverse lexicographic order.
    FloatValue<PR> err_coef=cast_exact(tm.error());
    ra[as]=1;
    rtm._append(ra,err_coef);
    ra[as]=0;

    // Copy new terms
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
        UniformConstReference<MultiIndex> xa=iter->index();
        UniformConstReference<FloatValue<PR>> xv=iter->coefficient();
        for(SizeType j=0; j!=as; ++j) { ra[j]=xa[j]; }
        rtm._append(ra,xv);
    }
    return std::move(rtm);
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_discard_variables(const TaylorModel<ValidatedTag,F>& tm, Array<SizeType> const& discarded_variables) {
    for(SizeType i=0; i!=discarded_variables.size()-1u; ++i) {
        ARIADNE_PRECONDITION(discarded_variables[i]<discarded_variables[i+1u]);
    }
    ARIADNE_PRECONDITION(discarded_variables[discarded_variables.size()-1u]<tm.argument_size());

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
    FloatError<PR> derr=mag(tm.error()); // Magnitude of discarded terms
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        UniformConstReference<MultiIndex> xa=iter->index();
        UniformConstReference<FloatValue<PR>> xv=iter->coefficient();
        Bool keep=true;
        for(SizeType k=0; k!=number_of_discarded_variables; ++k) {
            if(xa[discarded_variables[k]]!=0) {
                derr += mag(iter->coefficient());
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

    FloatError<PR> e=nul(this->error());
    for(typename TaylorModel<ValidatedTag,F>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformReference<FloatValue<PR>> xv=xiter->coefficient();
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
    FloatError<PR>& re=r.error();
    typename TaylorModel<ValidatedTag,F>::Iterator riter=r.begin();
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<FloatValue<PR>> xv=xiter->coefficient();
        Nat c=xa[k];
        if(c!=0) {
            UniformReference<MultiIndex> ra=riter->index();
            UniformReference<FloatValue<PR>> rv=riter->coefficient();
            ra=xa; ra[k]-=1;
            rv=mul_err(xv,c,re);
            ++riter;
        }
    }

    r.expansion().resize(static_cast<SizeType>(riter - r.begin()));
}




template<class F> TaylorModel<ValidatedTag,F> derivative(const TaylorModel<ValidatedTag,F>& x, SizeType k) {
    TaylorModel<ValidatedTag,F> rx=x; rx.differentiate(k); return rx;

    ARIADNE_ASSERT(k<x.argument_size());
    typedef typename F::PrecisionType PR;

    MultiIndex ra(x.argument_size()); FloatValue<PR> rv; Nat c;

    TaylorModel<ValidatedTag,F> r(x.argument_size(),x.sweeper());
    FloatError<PR>& re=r.error();
    for(typename TaylorModel<ValidatedTag,F>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<FloatValue<PR>> xv=xiter->coefficient();
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

template<class F> auto TaylorModel<ValidatedTag,F>::_evaluate(const TaylorModel<ValidatedTag,F>& tm, const Vector<FloatBounds<PR>>& x) -> FloatBounds<PR> {
    return horner_evaluate(tm.expansion(),x)+FloatBounds<typename F::PrecisionType>(-tm.error(),+tm.error());
}

template<class F> auto TaylorModel<ValidatedTag,F>::_evaluate(const TaylorModel<ValidatedTag,F>& tm, const Vector<FloatApproximation<PR>>& x) -> FloatApproximation<PR> {
    return horner_evaluate(tm.expansion(),x);
}


template<class F> auto TaylorModel<ValidatedTag,F>::_gradient(const TaylorModel<ValidatedTag,F>& tm, const Vector<FloatBounds<PR>>& x) -> Covector<FloatBounds<PR>> {
    Vector< Differential<FloatBounds<PR>> > dx=Differential<FloatBounds<PR>>::variables(1u,x);
    Differential<FloatBounds<PR>> df=horner_evaluate(tm.expansion(),dx)+FloatBounds<PR>(-tm.error(),+tm.error());
    return gradient(df);
}



template<class F> TaylorModel<ValidatedTag,F>
TaylorModel<ValidatedTag,F>::_compose(TaylorModel<ValidatedTag,F> const& x, Vector<TaylorModel<ValidatedTag,F>> const& y) {
    return horner_evaluate(x.expansion(),y)+NumericType(-x.error(),+x.error());
}

template<class F> Void TaylorModel<ValidatedTag,F>::unscale(IntervalDomainType const& codom) {
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
    Interval<FloatValue<PR>> ivl=convert_interval(codom,this->precision());
    ARIADNE_ASSERT_MSG(ivl.lower()<=ivl.upper(),"Cannot unscale TaylorModel<ValidatedTag,F> "<<tm<<" from empty interval "<<ivl);

    if(ivl.lower()==ivl.upper()) {
        tm=0;
        // Uncomment out line below to make unscaling to a singleton interval undefined
        //tm.clear(); tm.set_error(+inf);
    } else {
        auto c=ivl.centre();
        auto s=rec(ivl.radius());
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
TaylorModel<ValidatedTag,F>::_partial_evaluate(const TaylorModel<ValidatedTag,F>& x, SizeType k, NumericType c)
{
    const SizeType as=x.argument_size();
    Vector<TaylorModel<ValidatedTag,F>> y(as,TaylorModel<ValidatedTag,F>(as-1u,x.sweeper()));
    for(SizeType i=0; i!=k; ++i) { y[i]=TaylorModel<ValidatedTag,F>::coordinate(as-1u,i,x.sweeper()); }
    y[k]=TaylorModel<ValidatedTag,F>::constant(as-1u,c,x.sweeper());
    for(SizeType i=k+1; i!=as; ++i) { y[i]=TaylorModel<ValidatedTag,F>::coordinate(as-1u,i-1u,x.sweeper()); }
    return compose(x,y);
}



template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_embed(SizeType as1, const TaylorModel<ValidatedTag,F>& tm2, SizeType as3) {
    return TaylorModel<ValidatedTag,F>(embed(as1,tm2.expansion(),as3),tm2.error(),tm2.sweeper());
}


template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_split(const TaylorModel<ValidatedTag,F>& tm, SizeType k, SplitPart h) {
    const DegreeType deg=tm.degree();
    const SizeType as=tm.argument_size();
    SweeperType swp=tm.sweeper();

    TaylorModel<ValidatedTag,F> r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(typename TaylorModel<ValidatedTag,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const DegreeType ak=iter->index()[k];
        UniformReference<FloatValue<PR>> c=iter->coefficient();
        c/=pow(two,ak);
    }

    if(h==SplitPart::MIDDLE) { return r; }
    Int tr=( h==SplitPart::UPPER ? +1 : -1 );

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<TaylorModel<ValidatedTag,F>> ary(deg+1u,TaylorModel<ValidatedTag,F>(as,swp));
    for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->index();
        UniformConstReference<FloatValue<PR>> c=iter->coefficient();
        DegreeType ak=a[k];
        a[k]=0u;
        ary[ak].expansion().append(a,c);
    }

    FloatError<PR> re=r.error();
    r.clear();
    r.set_error(re);

    for(DegreeType i=0; i<=deg; ++i) {
        for(DegreeType j=i; j<=deg; ++j) {
            Int sf=bin(j,i);
            if(tr==-1 && (j-i)%2==1) { sf=-sf; }
            r+=ary[j]*sf;
            for(typename TaylorModel<ValidatedTag,F>::Iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->index()[k];
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
    FloatError<PR> r(0u,this->precision());
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->index().degree()!=0) {
            r+=mag(iter->coefficient());
        }
    }
    r+=this->error();
    return r;
}

template<class F> typename TaylorModel<ValidatedTag,F>::NormType TaylorModel<ValidatedTag,F>::norm() const {
    FloatError<PR> r(0u,this->precision());
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=mag(iter->coefficient());
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
    return (mag(tm1.value()-tm2.value()).raw() > (tm1.error()+tm2.error()).raw());
}

template<class F> TaylorModel<ValidatedTag,F> TaylorModel<ValidatedTag,F>::_refinement(const TaylorModel<ValidatedTag,F>& x, const TaylorModel<ValidatedTag,F>& y) {
    TaylorModel<ValidatedTag,F> r(x.argument_size(),x.sweeper());

    FloatError<PR> max_error=nul(r.error());

    const FloatError<PR>& xe=x.error();
    const FloatError<PR>& ye=y.error();
    FloatValue<PR> rv,xv,yv;
    FloatDP xu,yu,mxl,myl,u,ml;
    MultiIndex a;

    typename TaylorModel<ValidatedTag,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<ValidatedTag,F>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() || yiter!=y.end()) {
        // Can't use const MultiIndex& here since references change as the iterators change
        // We would need to use a smart reference
        //UniformConstReference<MultiIndex> xa=xiter->index();
        //UniformConstReference<MultiIndex> ya=yiter->index();
        if(xiter==x.end()) {
            a=yiter->index();
            yv=yiter->coefficient();
            xv=0;
            ++yiter;
        } else if(yiter==y.end()) {
            a=xiter->index();
            xv=xiter->coefficient();
            yv=0;
            ++xiter;
        } else if(xiter->index()==yiter->index()) {
            a=xiter->index();
            xv=xiter->coefficient();
            yv=yiter->coefficient();
            ++xiter;
            ++yiter;
        } else if(xiter->index()<yiter->index()) {
            a=xiter->index();
            xv=xiter->coefficient();
            yv=0;
            ++xiter;
        } else { // xa>ya
            a=yiter->index();
            yv=yiter->coefficient();
            xv=0;
            ++yiter;
        }

        FloatBounds<PR> rve=refinement( xv.pm(xe), yv.pm(ye) );
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
    Expansion<MultiIndex,FloatValue<PR>> e=tm.expansion();
    e.sort(GradedIndexLess());
    e.write(os,variable_names);
    return os << "+/-" << tm.error() << ")";
}

template<class F> OutputStream& TaylorModel<ValidatedTag,F>::repr(OutputStream& os) const {
    return this->str(os);
}

///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::zeros(SizeType rs, SizeType as, SweeperType swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(rs,TaylorModel<ValidatedTag,F>::zero(as,swp));
    return result;
}


template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::constants(SizeType as, const Vector<ValidatedNumericType>& c, SweeperType swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(c.size(),TaylorModel<ValidatedTag,F>::zero(as,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel<ValidatedTag,F>::constant(as,c[i],swp);
    }
    return result;
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::coordinates(SizeType as, SweeperType swp)
{
    Vector<TaylorModel<ValidatedTag,F>> result(as,TaylorModel<ValidatedTag,F>::zero(as,swp));
    for(SizeType i=0; i!=as; ++i) { result[i]=TaylorModel<ValidatedTag,F>::coordinate(as,i,swp); }
    return result;
}

template<class F> Vector<TaylorModel<ValidatedTag,F>> TaylorModel<ValidatedTag,F>::scalings(const Vector<IntervalDomainType>& d, SweeperType swp)
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
    Vector<Differential<ValidatedNumericType>> dx(x.size(),x.size(),1u);
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
template<class F> Matrix<Value<F>>
jacobian_value(const Vector<TaylorModel<ValidatedTag,F>>& f) {
    typedef typename F::PrecisionType PR;
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<FloatValue<PR>> J(rs,as);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType j=0; j!=as; ++j) {
            a[j]=1; const FloatValue<PR> x=f[i][a]; J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian at the origin with respect to the variables args.
template<class F> Matrix<Value<F>>
jacobian_value(const Vector<TaylorModel<ValidatedTag,F>>& f, const Array<SizeType>& p) {
    typedef typename F::PrecisionType PR;
    const SizeType rs=f.size();
    const SizeType as=f.zero_element().argument_size();
    const SizeType ps=p.size();
    Matrix<FloatValue<PR>> J(rs,ps);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType k=0; k!=ps; ++k) {
            SizeType j=p[k]; a[j]=1; const FloatValue<PR> x=f[i][a]; J[i][k]=x; a[j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
template<class F> Matrix<UpperIntervalType>
jacobian_range(const Vector<TaylorModel<ValidatedTag,F>>& f) {
    typedef typename F::PrecisionType PR;
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<FloatBounds<PR>> J(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            for(SizeType k=0; k!=as; ++k) {
                const Nat c=a[k];
                if(c>0) {
                    UniformConstReference<FloatValue<PR>> x=iter->coefficient();
                    if(a.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=FloatBounds<PR>(-1,1)*x*c; }
                }
            }
        }
    }
    return Matrix<UpperIntervalType>(J);
}

// Compute the Jacobian over the unit domain, with respect to the variables p.
template<class F> Matrix<UpperIntervalType>
jacobian_range(const Vector<TaylorModel<ValidatedTag,F>>& f, const Array<SizeType>& p) {
    typedef typename F::PrecisionType PR;
    SizeType rs=f.size();
    SizeType ps=p.size();
    Matrix<FloatBounds<PR>> J(rs,ps);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<ValidatedTag,F>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            for(SizeType k=0; k!=ps; ++k) {
                SizeType j=p[k];
                const Nat c=a[j];
                if(c>0) {
                    UniformConstReference<FloatValue<PR>> x=iter->coefficient();
                    if(a.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=FloatBounds<PR>(-1,1)*x*c; }
                }
            }
        }
    }
    return Matrix<UpperIntervalType>(J);
}











template<class F> TaylorModel<ApproximateTag,F>::TaylorModel()
    : _expansion(0), _sweeper()
{
}

template<class F> TaylorModel<ApproximateTag,F>::TaylorModel(SizeType as)
    : _expansion(as), _sweeper()
{
}

template<class F> TaylorModel<ApproximateTag,F>::TaylorModel(SizeType as, SweeperType swp)
    : _expansion(as), _sweeper(swp)
{
}

template<class F> TaylorModel<ApproximateTag,F> TaylorModel<ApproximateTag,F>::create_constant(FloatApproximation<PR> c) const {
    TaylorModel<ApproximateTag,F> r(this->argument_size(),this->_sweeper);
    r._expansion.append(MultiIndex::zero(this->argument_size()),c);
    return std::move(r);
}

template<class F> TaylorModel<ApproximateTag,F> TaylorModel<ApproximateTag,F>::create_constant(ApproximateNumber c) const {
    return this->create_constant(NumericType(c,this->precision()));
}

template<class F> TaylorModel<ApproximateTag,F> TaylorModel<ApproximateTag,F>::create_ball(ErrorType) const {
    return TaylorModel<ApproximateTag,F>(this->argument_size(),this->_sweeper);
}


template<class F> Void TaylorModel<ApproximateTag,F>::iadd(const FloatApproximation<PR>& c)
{
    // Compute self+=c
    TaylorModel<ApproximateTag,F>& x=*this;
    if(decide(c==0)) { return; }
    if(x._expansion.empty()) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else if((x._expansion.end()-1)->index().degree()>0) {
        x._expansion.append(MultiIndex(x.argument_size()),c);
    } else {
        UniformReference<FloatApproximation<PR>> rv=(x._expansion.end()-1)->coefficient();
        rv+=c;
    }
}

template<class F> Void TaylorModel<ApproximateTag,F>::imul(const FloatApproximation<PR>& c)
{
    // Compute self*=c
    if(decide(c==0)) { this->clear(); return; }
    if(decide(c==1)) { return; }
    for(Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        iter->coefficient() *= c;
    }
}



template<class F> Void TaylorModel<ApproximateTag,F>::isma(const FloatApproximation<PR>& c, const TaylorModel<ApproximateTag,F>& y)
{
    const TaylorModel<ApproximateTag,F>& x=*this;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<", y="<<y);
    TaylorModel<ApproximateTag,F> r(x.argument_size());

    ConstIterator xiter=x._expansion.begin();
    ConstIterator yiter=y._expansion.begin();
    while(xiter!=x._expansion.end() && yiter!=y._expansion.end()) {
        if(xiter->index()<yiter->index()) {
            r._expansion.append(xiter->index(),xiter->coefficient());
            ++xiter;
        } else if(yiter->index()<xiter->index()) {
            r._expansion.append(yiter->index(),c*yiter->coefficient());
            ++yiter;
        } else {
            r._expansion.append(xiter->index(),xiter->coefficient()+c*yiter->coefficient());
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x._expansion.end()) {
        r._expansion.append(xiter->index(),xiter->coefficient());
        ++xiter;
    }
    while(yiter!=y._expansion.end()) {
        r._expansion.append(yiter->index(),c*yiter->coefficient());
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
    for(ConstIterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        ConstIterator yiter=y._expansion.begin();
        ConstIterator riter=r._expansion.begin();
        while(riter!=r._expansion.end() && yiter!=y._expansion.end()) {
            UniformConstReference<MultiIndex> ra=riter->index();
            sa=xiter->index()+yiter->index();
            if(sa==ra) {
                t._expansion.append(sa,riter->coefficient()+xiter->coefficient()*yiter->coefficient());
                ++riter; ++yiter;
            } else if(sa<ra) {
                t._expansion.append(sa,xiter->coefficient()*yiter->coefficient());
                ++yiter;
            } else {
                t._expansion.append(ra,riter->coefficient());
                ++riter;
            }
        }
        while(riter!=r._expansion.end()) {
            t._expansion.append(riter->index(),riter->coefficient());
            ++riter;
        }
        while(yiter!=y._expansion.end()) {
            t._expansion.append_sum(xiter->index(),yiter->index(),xiter->coefficient()*yiter->coefficient());
            ++yiter;

        }

        r._expansion.swap(t._expansion);
        t._expansion.clear();
    }
}


template<class F> struct AlgebraOperations<TaylorModel<ApproximateTag,F>>
    : NormedAlgebraOperations<TaylorModel<ApproximateTag,F>>
{
    typedef typename F::PrecisionType PR;
    typedef CanonicalNumericType<ApproximateTag,PR> ApproximateNumericType;
  public:
    using NormedAlgebraOperations<TaylorModel<ApproximateTag,F>>::apply;
    static TaylorModel<ApproximateTag,F> apply(Nul, TaylorModel<ApproximateTag,F> const& x) {
        return TaylorModel<ApproximateTag,F>(x.argument_size(),x.sweeper()); }
    static TaylorModel<ApproximateTag,F> apply(Pos, TaylorModel<ApproximateTag,F> x) {
        return std::move(x); }
    static TaylorModel<ApproximateTag,F> apply(Neg, TaylorModel<ApproximateTag,F> x) {
        x.imul(ApproximateNumericType(-1,x.precision())); return std::move(x); }
    static TaylorModel<ApproximateTag,F> apply(Add, TaylorModel<ApproximateTag,F> const& x, TaylorModel<ApproximateTag,F> const& y) {
        auto r=x; r.isma(ApproximateNumericType(-1,x.precision()),y); return std::move(r); }
    static TaylorModel<ApproximateTag,F> apply(Sub, TaylorModel<ApproximateTag,F> const& x, TaylorModel<ApproximateTag,F> const& y) {
        auto r=x; r.isma(ApproximateNumericType(-1,x.precision()),y); return std::move(r); }
    static TaylorModel<ApproximateTag,F> apply(Mul, TaylorModel<ApproximateTag,F> const& x, TaylorModel<ApproximateTag,F> const& y) {
        auto r=nul(x); r.ifma(x,y); return std::move(r); }
    static TaylorModel<ApproximateTag,F> apply(Add, TaylorModel<ApproximateTag,F> x, ApproximateNumericType const& c) {
        auto& r=x; r.iadd(c); return std::move(r); }
    static TaylorModel<ApproximateTag,F> apply(Mul, TaylorModel<ApproximateTag,F> x, ApproximateNumericType const& c) {
        auto& r=x; r.imul(c); return std::move(r); }
    static TaylorModel<ApproximateTag,F> apply(Max, TaylorModel<ApproximateTag,F> const& x, TaylorModel<ApproximateTag,F> const& y) {
        ARIADNE_NOT_IMPLEMENTED; }
    static TaylorModel<ApproximateTag,F> apply(Min, TaylorModel<ApproximateTag,F> const& x, TaylorModel<ApproximateTag,F> const& y) {
        ARIADNE_NOT_IMPLEMENTED; }
    static TaylorModel<ApproximateTag,F> apply(Abs, TaylorModel<ApproximateTag,F> const& x) {
        ARIADNE_NOT_IMPLEMENTED; }
};

template<class F> Void TaylorModel<ApproximateTag,F>::unscale(IntervalDomainType const& dom) {
    auto ivl=convert_interval(dom,this->precision());
    TaylorModel<ApproximateTag,F>& x=*this;
    x-=ivl.midpoint();
    x/=ivl.radius();
}



template<class F> typename TaylorModel<ApproximateTag,F>::NormType TaylorModel<ApproximateTag,F>::norm() const {
    FloatApproximation<PR> r(0u,this->precision());
    for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        r+=mag(iter->coefficient());
    }
    return NormType(r);
}

template<class F> typename TaylorModel<ApproximateTag,F>::CoefficientType TaylorModel<ApproximateTag,F>::average() const {
    return (*this)[MultiIndex(this->argument_size())];
}

template<class F> typename TaylorModel<ApproximateTag,F>::NormType TaylorModel<ApproximateTag,F>::radius() const {
    FloatApproximation<PR> r(0u,this->precision());
    for(auto iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        if(iter->index().degree()!=0) {
            r+=mag(iter->coefficient());
        }
    }
    return NormType(r);
}

template<class F> UnitBox TaylorModel<ApproximateTag,F>::domain() const {
    return UnitBox(this->argument_size(),UnitInterval());
}

template<class F> auto TaylorModel<ApproximateTag,F>::range() const -> RangeType {
    FloatApproximation<PR> av=this->average();
    FloatApproximation<PR> rad=this->radius();
    return Interval<FloatApproximation<PR>>(av-rad,av+rad);
}

template<class F> TaylorModel<ApproximateTag,F>& TaylorModel<ApproximateTag,F>::sweep() {
    this->_sweeper.sweep(this->_expansion);
    return *this;
}

template<class F> F TaylorModel<ApproximateTag,F>::tolerance() const {
    const ThresholdSweeper<F>* ptr=dynamic_cast<const ThresholdSweeper<F>*>(&static_cast<const SweeperInterface<F>&>(this->_sweeper));
    return (ptr) ? ptr->sweep_threshold() : F(std::numeric_limits<double>::epsilon(),this->precision());
}


template<class F> OutputStream& TaylorModel<ApproximateTag,F>::write(OutputStream& os) const {
    return os << "TM["<<this->argument_size()<<"](" << this->_expansion << ")";
}

template<class F> OutputStream& TaylorModel<ApproximateTag,F>::str(OutputStream& os) const {
    return os << this->_expansion;
}

template<class F> OutputStream& TaylorModel<ApproximateTag,F>::repr(OutputStream& os) const {
    return this->str(os);
}

} //namespace Ariadne


