/***************************************************************************
 *            taylor_model.tpl.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

template<class F> F UnknownError<F>::raw() const {
    return F(0u,this->precision()); }
template<class F> typename F::PrecisionType UnknownError<F>::precision() const {
    if constexpr (IsSame<F,FloatMP>::value) { return MultiplePrecision(64); }
    else { return typename F::PrecisionType(); } }
template<class F> UnknownError<F>::operator PositiveApproximation<F> () const {
    return PositiveApproximation<F>(0,this->precision()); }

template<class F> UnknownError<F> nul(UnknownError<F>) { return UnknownError<F>(); }
template<class F> UnknownError<F> mag(UnknownError<F>) { return UnknownError<F>(); }
template<class F> UnknownError<F> operator+(UnknownError<F>, UnknownError<F>) { return UnknownError<F>(); }
template<class F> UnknownError<F> operator*(UnknownError<F>, UnknownError<F>) { return UnknownError<F>(); }
template<class F> UnknownError<F>& operator+=(UnknownError<F>& e, UnknownError<F> const&) { return e; }
template<class F> UnknownError<F>& operator*=(UnknownError<F>& e, UnknownError<F> const&) { return e; }

template<class F> UnknownError<F> operator+(UnknownError<F>, PositiveApproximation<F> const&) { return UnknownError<F>(); }
template<class F> UnknownError<F>& operator+=(UnknownError<F>& e, PositiveApproximation<F> const&) { return e; }
template<class F> UnknownError<F>& operator*=(UnknownError<F>& e, PositiveApproximation<F> const&) { return e; }
template<class F> OutputStream& operator<<(OutputStream& os, UnknownError<F> const& e) { return os << "???"; }


template<class F> Bounds<F> fma(Bounds<F> const& x1, Bounds<F> const& x2, Bounds<F> y) { return x1*x2+y; }
template<class F> Bounds<F>& operator/=(Bounds<F>& x, TwoExp y) { return x*=Dyadic(rec(y)); }

template<class F> Approximation<F>& operator/=(Approximation<F>& x, TwoExp y) { return x/=Dyadic(y); }
template<class F> Approximation<F> fma(Approximation<F>const& x, Approximation<F> const& y, Approximation<F> z) {
    z._a = fma(near,x._a,y._a,z._a);; return std::move(z); }
template<class F> Approximation<F> pm(UnknownError<F> e) { return Approximation<F>(e.precision()); }



template<class F> PositiveUpperBound<F> mag(Interval<UpperBound<F>> const& ivl) {
    return cast_positive(max(-ivl.lower(),ivl.upper()));
}
template<class F> PositiveApproximation<F> mag(Interval<Approximation<F>> const& ivl) {
    return cast_positive(max(-ivl.lower(),ivl.upper()));
}


namespace {

static const double MACHINE_EPSILON = 2.2204460492503131e-16;

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }


Interval<FloatDPValue> convert_interval(Interval<FloatDPValue> const& ivl, DoublePrecision pr) {
    return ivl; }
Interval<FloatMPValue> convert_interval(Interval<FloatDPValue> const& ivl, MultiplePrecision pr) {
    return Interval<FloatMPValue>(FloatMP(ivl.lower().raw(),pr),FloatMP(ivl.upper().raw(),pr)); }

Interval<FloatDPUpperBound> const& convert_interval(Interval<FloatDPUpperBound> const& ivl, DoublePrecision) {
    return ivl; }
Interval<FloatDPUpperBound> convert_interval(Interval<FloatMPUpperBound> const& ivl, DoublePrecision pr) {
    return Interval<FloatDPUpperBound>(FloatDP(ivl.lower().raw(),down,pr),FloatDP(ivl.upper().raw(),up,pr)); }

Interval<FloatDPApproximation> const& convert_interval(Interval<FloatDPApproximation> const& ivl, DoublePrecision) {
    return ivl; }
Interval<FloatDPApproximation> convert_interval(Interval<FloatMPApproximation> const& ivl, DoublePrecision pr) {
    return Interval<FloatDPApproximation>(FloatDP(ivl.lower().raw(),near,pr),FloatDP(ivl.upper().raw(),near,pr)); }

template<class FLT> Bool is_same_as_zero(Approximation<FLT> const& xa) { return xa.raw()==0; }
template<class FLT> Bool is_same_as_zero(Bounds<FLT> const& xb) { return xb.lower_raw()==0 && xb.upper_raw()==0; }
template<class FLT> Bool is_same_as_zero(Value<FLT> const& xv) { return xv.raw()==0; }

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
Void SweeperBase<F>::_sweep(Expansion<MultiIndex,FloatBounds<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatBounds<PR>>::Iterator curr=p.begin();

    // FIXME: Not needed, but added to pair with rounding mode change below
    F::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(adv->index(),mag(adv->coefficient()).raw())) {
            //te+=abs(adv->coefficient());
            te+mag(adv->coefficient());
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
Void RelativeSweeperBase<F>::_sweep(Expansion<MultiIndex,FloatBounds<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatBounds<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=radius(p)+e;

    // FIXME: Not needed, but added to pair with rounding mode change below
    F::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(mag(adv->coefficient()).raw(),nrm.raw())) {
            //te+=abs(adv->coefficient());
            te+=mag(adv->coefficient());
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



template<class P, class F> TaylorModel<P,F>::TaylorModel()
    : _expansion(0), _error(0u), _sweeper()
{
}


template<class P, class F> TaylorModel<P,F>::TaylorModel(SizeType as, SweeperType swp)
    : _expansion(as), _error(0u), _sweeper(swp)
{
}

template<class P, class F> TaylorModel<P,F>::TaylorModel(const Expansion<MultiIndex,CoefficientType>& f, const ErrorType& e, SweeperType swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->cleanup();
}

template<class P, class F> TaylorModel<P,F>::TaylorModel(const Expansion<MultiIndex,RawFloatType>& f, const RawFloatType& e, SweeperType swp)
//    : TaylorModel(reinterpret_cast<Expansion<MultiIndex,CoefficientType>const&>(f),reinterpret_cast<ErrorType const&>(e),swp)
    : TaylorModel(static_cast<Expansion<MultiIndex,CoefficientType>>(f),static_cast<ErrorType>(e),swp)
{
}

template<class P, class F> TaylorModel<P,F>::TaylorModel(const Expansion<MultiIndex,double>& f, const double& e, SweeperType swp)
    : TaylorModel(Expansion<MultiIndex,CoefficientType>(Expansion<MultiIndex,Dyadic>(f),swp.precision()),ErrorType(Dyadic(e),swp.precision()),swp)
{
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::scaling(SizeType as, SizeType j, const IntervalDomainType& codom, SweeperType swp) {
    TaylorModel<P,F> r(as,swp);
    auto ivl=convert_interval(codom,r.precision());
    r.set_gradient(j,1);
    r*=ivl.radius();
    r+=ivl.midpoint();
    return r;
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::create() const {
    return TaylorModel<P,F>(this->argument_size(),this->_sweeper);
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::create_zero() const {
    return TaylorModel<P,F>(this->argument_size(),this->_sweeper);
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::create_constant(NumericType c) const {
    return TaylorModel<P,F>::constant(this->argument_size(),c,this->_sweeper);
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::create_constant(GenericNumericType c) const {
    return this->create_constant(NumericType(c,this->precision()));
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::create_coordinate(SizeType j) const {
    ARIADNE_PRECONDITION(j<this->argument_size());
    TaylorModel<P,F> r(this->argument_size(),this->_sweeper);
    CoefficientType one(1,this->precision());
    r._expansion.append(MultiIndex::unit(this->argument_size(),j),one);
    return r;
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::create_ball(ErrorType e) const {
    ARIADNE_DEBUG_PRECONDITION(e.raw()>=0);
    TaylorModel<P,F> r(this->argument_size(),this->_sweeper);
    r._error=e;
    return r;
}

template<class P, class F> Void TaylorModel<P,F>::swap(TaylorModel<P,F>& tm) {
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
    std::swap(this->_sweeper,tm._sweeper);
}

template<class P, class F> Void TaylorModel<P,F>::clear() {
    this->_expansion.clear();
    this->_error=0u;
}

template<class P, class F> DegreeType TaylorModel<P,F>::degree() const {
    DegreeType deg=0u;
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        deg=std::max(deg,iter->index().degree());
    }
    return deg;
}


template<class F> Value<F> set_err(Bounds<F> const& x, Error<F>& e) {
    e+=x.error();
    return x.value();
}

template<class F> Approximation<F> const& set_err(Approximation<F> const& x, UnknownError<F>& e) {
    return x;
}


template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::operator=(const NumericType& c) {
    this->clear();
    CoefficientType m=set_err(c,this->_error);
    if(not is_same_as_zero(m)) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    return *this;
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::operator=(const GenericNumericType& c) {
    return *this = NumericType(c,this->precision());
}

namespace { // Internal code for arithmetic

template<class F> struct ValidatedApproximation {
    Bounds<F> _v; Approximation<F> _a;
    ValidatedApproximation(Bounds<F>const& x) : _v(x), _a(x) { }
    operator Bounds<F> const& () const { return _v; }
    LowerBound<F> lower() const { return _v.lower(); }
    Approximation<F> middle() const { return _a; }
    UpperBound<F> upper() const { return _v.upper(); }
    F const& lower_raw() const { return _v.lower_raw(); }
    F const& middle_raw() const { return _a.raw(); }
    F const& upper_raw() const { return _v.upper_raw(); }
};
template<class F> Approximation<F> const& make_validated_approximation(Approximation<F> const& x) { return x; }
template<class F> ValidatedApproximation<F> make_validated_approximation(Bounds<F> const& x) { return ValidatedApproximation<F>(x); }


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

template<class F> Bounds<F> add_err(Bounds<F> const& x1, Bounds<F> const& x2, Error<F>& e) {
    return add(x1,x2);
}

template<class F> Approximation<F> add_err(Approximation<F> const& x1, Approximation<F> const& x2, UnknownError<F>& e) {
    return add(x1,x2);
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

template<class F> Approximation<F> sub_err(Approximation<F> const& x1, Approximation<F> const& x2, UnknownError<F>& e) {
    return sub(x1,x2);
}

template<class F> Value<F> mul_no_err(Value<F> const& x1, Value<F> const& x2) {
    F::set_rounding_to_nearest();
    Value<F> r(x1.raw() * x2.raw());
    F::set_rounding_upward();
    return r;
}

template<class F> Approximation<F> mul_no_err(Approximation<F> const& x1, Approximation<F> const& x2) {
    return x1*x2;
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

template<class F> Value<F> mul_err(Value<F> const& x1, Nat n2, Error<F>& e) {
    return mul_err(x1,Value<F>(n2,x1.precision()),e);
}

template<class F> Bounds<F> mul_err(Bounds<F> const& x1, Bounds<F> const& x2, Error<F>& e) {
    return mul(x1,x2);
}

template<class F> Bounds<F> mul_err(Bounds<F> const& x1, ValidatedApproximation<F> const& x2, Error<F>& e) {
    return mul(x1,static_cast<Bounds<F>>(x2));
}

template<class F> Bounds<F> mul_err(Bounds<F> const& x1, Nat const& n2, Error<F>& e) {
    return mul(x1,n2);
}

template<class F> Approximation<F> mul_err(Approximation<F> const& x1, Approximation<F> const& x2, UnknownError<F>& e) {
    return mul(x1,x2);
}

template<class F> Approximation<F> mul_err(Approximation<F> const& x1, Nat n2, UnknownError<F>& e) {
    return mul_err(x1,Approximation<F>(n2,x1.precision()),e);
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

template<class F> Bounds<F> div_err(Bounds<F> const& x1, Bounds<F> const& x2, Error<F>& e) {
    return div(x1,x2);
}

template<class F> Bounds<F> div_err(Bounds<F> const& x1, Nat n2, Error<F>& e) {
    return div(x1,n2);
}

template<class F> Approximation<F> div_err(Approximation<F> const& x1, Approximation<F> const& x2, UnknownError<F>& e) {
    return div(x1,x2);
}

template<class F> Approximation<F> div_err(Approximation<F> const& x1, Nat n2, UnknownError<F>& e) {
    return div_err(x1,Approximation<F>(n2,x1.precision()),e);
}





template<class F> Value<F> fma_err(Value<F> const& x, ValidatedApproximation<F> const& c, Value<F> y, Error<F>& e) {
    F const& xv=x.raw();
    F const& cu=c.upper_raw();
    F const& cm=c.middle_raw();
    F const& cl=c.lower_raw();
    F const& yv=y.raw();
    F& re=e.raw();
    F::set_rounding_to_nearest();
    F rv=xv*cm+yv;
    F::set_rounding_upward();
    F u,ml;
    if(xv>=0) {
        F mcl=-cl;
        u=xv*cu+yv;
        ml=xv*mcl-yv;
    } else {
        F mcu=-cu;
        u=xv*cl+yv;
        ml=xv*mcu-yv;
    }
    re+=(u+ml)/2;
    return Value<F>(rv);
}

template<class F> Value<F> fma_err(Value<F> const& x, Bounds<F> const& c, Value<F> y, Error<F>& e) {
    return fma_err(x,ValidatedApproximation<F>(c),y,e);
}

template<class F> Bounds<F> fma_err(Bounds<F> const& x, Bounds<F> const& y, Bounds<F> z, Error<F>& e) {
    return fma(x,y,z);
}

template<class F> Bounds<F> fma_err(Bounds<F> const& x, ValidatedApproximation<F> const& y, Bounds<F> z, Error<F>& e) {
    return fma(x,static_cast<Bounds<F>>(y),z);
}

template<class F> Approximation<F> fma_err(Approximation<F> const& x, Approximation<F> const& y, Approximation<F> z, UnknownError<F>& e) {
    return fma(x,y,z);
}

// Inplace negation
template<class P, class F> Void _neg(TaylorModel<P,F>& r)
{
    for(auto iter=r.begin(); iter!=r.end(); ++iter) {
        iter->coefficient()=-iter->coefficient();
    }
}


template<class P, class F> Void _scal(TaylorModel<P,F>& r, const TwoExp& c)
{
    if (c.exponent()==0) { return; }
    for(typename TaylorModel<P,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->coefficient()*=c;
    }
    r.error()*=Error<F>(c);
 }

template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const Value<F>& c) {
    typedef ValidatedTag P;
    using ErrorType = typename TaylorModel<P,F>::ErrorType;
    ErrorType e(0u,r.error_precision()); // The maximum accumulated error
    for(typename TaylorModel<P,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->coefficient() = mul_err(riter->coefficient(),c,e);
    }
    auto& re=r.error();
    re*=abs(c);
    re+=e;
}

template<class F> Void _scal(TaylorModel<ValidatedTag,F>& r, const Bounds<F>& c)
{
    typedef ValidatedTag P;
    typedef typename F::PrecisionType PR;
    //std::cerr<<"TaylorModel<P,F>::scal(Float64Bounds c) c="<<c<<std::endl;
    ARIADNE_ASSERT_MSG(is_finite(c.lower().raw()) && is_finite(c.upper().raw()),"scal(tm,c): tm="<<r<<", c="<<c);
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);

    const F inf = F::inf(r.precision());
    if(r.error().raw()==inf) {
        r.expansion().clear(); return;
    }

    FloatError<PR> e=nul(r.error());
    ValidatedApproximation<F> clmu=c;
    for(typename TaylorModel<P,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        UniformReference<FloatValue<PR>> rv=riter->coefficient();
        rv=mul_err(rv,clmu,e);
    }
    r.error()*=mag(c);
    r.error()+=e;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}

template<class F> Void _scal(TaylorModel<ValidatedTag,Bounds<F>>& r, const Bounds<F>& c)
{
    typedef ValidatedTag P;
    using CoefficientType = typename TaylorModel<P,Bounds<F>>::CoefficientType;
    using ErrorType = typename TaylorModel<P,Bounds<F>>::ErrorType;
    //std::cerr<<"TaylorModel<P,F>::scal(Float64Bounds c) c="<<c<<std::endl;
    ARIADNE_ASSERT_MSG(is_finite(c.lower().raw()) && is_finite(c.upper().raw()),"scal(tm,c): tm="<<r<<", c="<<c);

    const F inf = F::inf(r.precision());
    if(r.error().raw()==inf) {
        r.expansion().clear(); return;
    }

    ErrorType e=nul(r.error());
    for(typename TaylorModel<P,Bounds<F>>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        UniformReference<CoefficientType> rv=riter->coefficient();
        rv=mul_err(rv,c,e);
    }
    r.error()*=mag(c);
    r.error()+=e;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}

template<class F> Void _scal(TaylorModel<ApproximateTag,F>& r, const Approximation<F>& c)
{
    typedef ApproximateTag P;
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using ErrorType = typename TaylorModel<P,F>::ErrorType;
    //std::cerr<<"TaylorModel<P,F>::scal(Float64Bounds c) c="<<c<<std::endl;
    ARIADNE_ASSERT_MSG(is_finite(c.raw()),"scal(tm,c): tm="<<r<<", c="<<c);

    const F inf = F::inf(r.precision());
    if(r.error().raw()==inf) {
        r.expansion().clear(); return;
    }

    ErrorType e=nul(r.error());
    for(typename TaylorModel<P,F>::Iterator riter=r.begin(); riter!=r.end(); ++riter) {
        UniformReference<CoefficientType> rv=riter->coefficient();
        rv=mul_err(rv,c,e);
    }
    r.error()*=mag(c);
    r.error()+=e;
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
}


struct UnitMultiIndex { SizeType argument_size; SizeType unit_index; };


template<class P, class F> inline Void _incr(TaylorModel<P,F>& r, const MultiIndex& a) {
    for(typename TaylorModel<P,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->index()+=a;
    }
}

template<class P, class F> inline Void _incr(TaylorModel<P,F>& r, SizeType j) {
    for(typename TaylorModel<P,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++iter->index()[j];
    }
}


template<class P, class F> inline Void _acc(TaylorModel<P,F>& r, const Value<F>& c) {
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



template<class P, class F> inline Void _acc(TaylorModel<P,F>& r, const Bounds<F>& c)
{
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using ErrorType = typename TaylorModel<P,F>::ErrorType;

    // Compute self+=c
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);

    const F inf = F::inf(r.precision());
    if(c.lower().raw()==-inf || c.upper().raw()==+inf) {
        r.clear();
        r.set_error(mag(CoefficientType(inf)));
        return;
    }

    if(c.lower().raw()==-c.upper().raw()) { // The midpoint of the interval is zero, so no need to change constant term
        r.error()+=mag(c);
        return;
    }

    if(r.expansion().empty() || r.expansion().back().index().degree()>0) { // Append a constant term zero
        r._append(MultiIndex(r.argument_size()),CoefficientType(0,r.precision()));
    }

    UniformReference<CoefficientType> rv=(r.end()-1)->coefficient();
    ErrorType& re=r.error();
    rv=add_err(rv,c,re);

    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}

template<class P, class F> inline Void _acc(TaylorModel<P,Bounds<F>>& r, const Bounds<F>& c)
{
    using CoefficientType = typename TaylorModel<P,Bounds<F>>::CoefficientType;
    using ErrorType = typename TaylorModel<P,Bounds<F>>::ErrorType;

    // Compute self+=c
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);

    const F inf = F::inf(r.precision());
    if(c.lower().raw()==-inf || c.upper().raw()==+inf) {
        r.clear();
        r.set_error(mag(CoefficientType(inf)));
        return;
    }

    if(r.expansion().empty() || r.expansion().back().index().degree()>0) { // Append a constant term zero
        r._append(MultiIndex(r.argument_size()),CoefficientType(0,r.precision()));
    }

    UniformReference<CoefficientType> rv=(r.end()-1)->coefficient();
    ErrorType& re=r.error();
    rv=add_err(rv,c,re);

    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}

template<class P, class F> inline Void _acc(TaylorModel<P,F>& r, const Approximation<F>& c)
{
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using ErrorType = typename TaylorModel<P,F>::ErrorType;

    // Compute self+=c
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);

    const F inf = F::inf(r.precision());
    if(c.raw()==-inf || c.raw()==+inf) {
        r.clear();
        return;
    }

    if(r.expansion().empty() || r.expansion().back().index().degree()>0) { // Append a constant term zero
        r._append(MultiIndex(r.argument_size()),CoefficientType(0,r.precision()));
    }

    UniformReference<CoefficientType> rv=(r.end()-1)->coefficient();
    ErrorType& re=r.error();
    rv=add_err(rv,c,re);
}


// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
template<class P, class F> inline Void _add(TaylorModel<P,F>& r, const TaylorModel<P,F>& x, const TaylorModel<P,F>& y)
{
    ARIADNE_PRECONDITION(r.number_of_nonzeros()==0);
    using ErrorType = typename TaylorModel<P,F>::ErrorType;
    ErrorType e=nul(r.error());
    typename TaylorModel<P,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<P,F>::ConstIterator yiter=y.begin();
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


template<class P, class F> inline Void _acc(TaylorModel<P,F>& r, const TaylorModel<P,F>& x)
{
    TaylorModel<P,F> s(r.argument_size(),r.sweeper()); _add(s,r,x); s.swap(r);
    ARIADNE_DEBUG_ASSERT_MSG(r.error()>=0,r);
}


template<class P, class F> inline Void _sub(TaylorModel<P,F>& r, const TaylorModel<P,F>& x, const TaylorModel<P,F>& y)
{
    ARIADNE_PRECONDITION(r.number_of_nonzeros()==0);
    typename TaylorModel<P,F>::ErrorType e(r.error().precision());
    typename TaylorModel<P,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<P,F>::ConstIterator yiter=y.begin();
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


template<class P, class F> inline Void _sma(TaylorModel<P,F>& r, const TaylorModel<P,F>& x, const typename TaylorModel<P,F>::NumericType& c, const TaylorModel<P,F>& y)
{
    typedef typename F::PrecisionType PR;
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using ErrorType = typename TaylorModel<P,F>::ErrorType;

    if constexpr (IsSame<decltype(c),FloatBounds<PR>>::value) {
        ARIADNE_DEBUG_ASSERT_MSG(c.lower().raw()<=c.upper().raw(),c);
        ARIADNE_DEBUG_ASSERT_MSG(x.error().raw()>=0,"x="<<x);
        ARIADNE_DEBUG_ASSERT_MSG(y.error().raw()>=0,"y="<<y);
    }

    ErrorType err=nul(r.error()); // Twice the maximum accumulated error
    auto clmu=make_validated_approximation(c);

    // Compute r=x+y, assuming r is empty
    RawFloat<PR>::set_rounding_upward();
    typename TaylorModel<P,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<P,F>::ConstIterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->index()<yiter->index()) {
            UniformConstReference<CoefficientType> xv=xiter->coefficient();
            r.expansion().append(xiter->index(),mul_err(xv,c,err));
            ++xiter;
        } else if(yiter->index()<xiter->index()) {
            UniformConstReference<CoefficientType> yv=yiter->coefficient();
            r.expansion().append(yiter->index(),yv);
            ++yiter;
        } else {
            UniformConstReference<CoefficientType> xv=xiter->coefficient();
            UniformConstReference<CoefficientType> yv=yiter->coefficient();
            r.expansion().append(xiter->index(),fma_err(xv,clmu,yv,err));
            ++xiter; ++yiter;
        }
    }

    while(xiter!=x.end()) {
        UniformConstReference<CoefficientType> xv=xiter->coefficient();
        r.expansion().append(xiter->index(),mul_err(xv,clmu,err));
        ++xiter;
    }
    while(yiter!=y.end()) {
        UniformConstReference<CoefficientType> yv=yiter->coefficient();
        r.expansion().append(yiter->index(),yv);
        ++yiter;
    }

    r.error()=x.error()*mag(c)+y.error();
    r.error()+=err;

    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}



// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
template<class P, class F> inline Void _mul(TaylorModel<P,F>& r, const TaylorModel<P,F>& x, const TaylorModel<P,F>& y)
{
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using ErrorType = typename TaylorModel<P,F>::ErrorType;

    const SizeType as=r.argument_size();
    TaylorModel<P,F> t(as,r.sweeper());
    TaylorModel<P,F> s(as,r.sweeper());
    MultiIndex ta(as);
    for(typename TaylorModel<P,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        ErrorType te=nul(r.error()); // trucation error
        ErrorType re=nul(r.error()); // roundoff error
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<CoefficientType> xv=xiter->coefficient();
        for(typename TaylorModel<P,F>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            UniformConstReference<MultiIndex> ya=yiter->index();
            UniformConstReference<CoefficientType> yv=yiter->coefficient();
            ta=xa+ya;
            CoefficientType tv=mul_err(xv,yv,re);
            // NOTE: Previously, we allowed to discard terms immediately since Sweeper() had a discard methd
            // if(r.sweeper().discard(ta,tv)) { te+=mag(xv)*mag(yv); }
            t._append(ta,tv);
//            re+=(xv*yv).error();
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

    ErrorType xs=nul(r.error());
    for(typename TaylorModel<P,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=mag(xiter->coefficient());
    }

    ErrorType ys=nul(r.error());
    for(typename TaylorModel<P,F>::ConstIterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=mag(yiter->coefficient());
    }

    ErrorType& re=r.error();
    const ErrorType& xe=x.error();
    const ErrorType& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    return;
}


} // namespace


///////////////////////////////////////////////////////////////////////////////

// Inplace arithmetical operations for Algebra concept

template<class P, class F> Void TaylorModel<P,F>::iadd(const NumericType& c)
{
    _acc(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class P, class F> Void TaylorModel<P,F>::imul(const NumericType& c)
{
    _scal(*this,c);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class P, class F> Void TaylorModel<P,F>::isma(const NumericType& c, const TaylorModel<P,F>& y)
{
    TaylorModel<P,F>& x=*this;
    TaylorModel<P,F> r=this->create();
    _sma(r,y,c,x);
    this->swap(r);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

template<class P, class F> Void TaylorModel<P,F>::ifma(const TaylorModel<P,F>& x, const TaylorModel<P,F>& y)
{
    _mul(*this,x,y);
    this->sweep();
    ARIADNE_DEBUG_ASSERT_MSG(this->error().raw()>=0,*this);
}

/*
template<class P, class F> struct AlgebraOperations<TaylorModel<P,F>>
    : NormedAlgebraOperations<TaylorModel<P,F>>
{
    typedef TaylorModel<P,F> ModelType;
    typedef typename TaylorModel<P,F>::NumericType NumericType;
    typedef typename ModelType::RangeType RangeType;
  public:
    using NormedAlgebraOperations<TaylorModel<P,F>>::apply;
    static ModelType apply(Nul, ModelType const& x) {
        return ModelType(x.argument_size(),x.sweeper()); }
    static ModelType apply(Pos, ModelType x) {
        return x; }
    static ModelType apply(Neg, ModelType x) {
        x.imul(NumericType(-1)); return x; }
    static ModelType apply(Add, ModelType const& x, ModelType const& y) {
        auto r=x; r.isma(NumericType(+1),y); return r; }
    static ModelType apply(Sub, ModelType const& x, ModelType const& y) {
        auto r=x; r.isma(NumericType(-1),y); return r; }
    static ModelType apply(Mul, ModelType const& x, ModelType const& y) {
        auto r=nul(x); r.ifma(x,y); return r; }
    static ModelType apply(Add, ModelType x, NumericType const& c) {
        x.iadd(c); return x; }
    static ModelType apply(Mul, ModelType x, NumericType const& c) {
        x.imul(c); return x; }
    static ModelType apply(Max, ModelType const& x, ModelType const& y);
    static ModelType apply(Min, ModelType const& x, ModelType const& y);
// TODO: Should be able to automatically generate these operations
    static ModelType apply(Max, ModelType const& x, NumericType const& c);
    static ModelType apply(Min, ModelType const& x, NumericType const& c);
    static ModelType apply(Max, NumericType const& c, ModelType const& x);
    static ModelType apply(Min, NumericType const& c, ModelType const& x);
    static ModelType apply(Abs, ModelType const& x);
};
*/

template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Nul, ModelType const& x) -> ModelType {
    return ModelType(x.argument_size(),x.sweeper()); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Pos, ModelType x) -> ModelType {
    return std::move(x); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Neg, ModelType x) -> ModelType {
    x.imul(NumericType(-1)); return std::move(x); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Add, ModelType const& x, ModelType const& y) -> ModelType {
    auto r=x; r.isma(NumericType(+1),y); return std::move(r); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Sub, ModelType const& x, ModelType const& y) -> ModelType {
    auto r=x; r.isma(NumericType(-1),y); return std::move(r); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Mul, ModelType const& x, ModelType const& y) -> ModelType {
    auto r=nul(x); r.ifma(x,y); return std::move(r); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Add, ModelType x, NumericType const& c) -> ModelType {
    auto& r=x; r.iadd(c); return std::move(r); }
template<class P, class F> auto AlgebraOperations<TaylorModel<P,F>>::apply(Mul, ModelType x, NumericType const& c) -> ModelType {
    auto& r=x; r.imul(c); return std::move(r); }

// TODO: Should be able to automatically generate these operations
/*
template<class P, class F> auto AlgebraOperations<P,F>::apply(Max, ModelType const& x, ModelType const& y) -> ModelType;
template<class P, class F> auto AlgebraOperations<P,F>::apply(Min, ModelType const& x, ModelType const& y) -> ModelType;

template<class P, class F> auto AlgebraOperations<P,F>::apply(Max, ModelType const& x, NumericType const& c) -> ModelType;
template<class P, class F> auto AlgebraOperations<P,F>::apply(Min, ModelType const& x, NumericType const& c) -> ModelType;
template<class P, class F> auto AlgebraOperations<P,F>::apply(Max, NumericType const& c, ModelType const& x) -> ModelType;
template<class P, class F> auto AlgebraOperations<P,F>::apply(Min, NumericType const& c, ModelType const& x) -> ModelType;
template<class P, class F> auto AlgebraOperations<P,F>::apply(Abs, ModelType const& x) -> ModelType;
*/

///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::sort() {
    this->_expansion.sort();
    return *this;
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::unique()
{
    typename TaylorModel<P,F>::ConstIterator advanced =this->begin();
    typename TaylorModel<P,F>::ConstIterator end =this->end();
    typename TaylorModel<P,F>::Iterator current=this->begin();
    ErrorType e=nul(this->error());
    while(advanced!=end) {
        current->index()=advanced->index();
        CoefficientType rv=advanced->coefficient();
        ++advanced;
        while(advanced!=end && advanced->index()==current->index()) {
            UniformConstReference<CoefficientType> xv=advanced->coefficient();
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

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::sweep() {
//    this->_sweeper.sweep(reinterpret_cast<Expansion<MultiIndex,F>&>(this->_expansion),reinterpret_cast<F&>(this->_error));
    this->_sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::sweep(const SweeperType& sweeper) {
    sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::simplify() {
    return this->sweep();
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::simplify(const PropertiesType& properties) {
    return this->sweep(properties);
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::cleanup() {
    this->sort();
    this->unique();
    this->sweep();
    return *this;
}

template<class P, class F> TaylorModel<P,F>& TaylorModel<P,F>::clobber() {
    this->_error=0u;
    return *this;
}



///////////////////////////////////////////////////////////////////////////////

// Accuracy control

template<class P, class F> auto TaylorModel<P,F>::tolerance() const -> RawFloatType {
    typedef RawFloatType FLT;
    const ThresholdSweeper<FLT>* ptr=dynamic_cast<const ThresholdSweeper<FLT>*>(&static_cast<const SweeperInterface<FLT>&>(this->_sweeper));
    return (ptr) ? ptr->sweep_threshold() : FLT(std::numeric_limits<double>::epsilon(),this->precision());
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

template<class P, class F> UnitBox TaylorModel<P,F>::domain() const
{
    return UnitBox(this->argument_size(),UnitInterval());
}

template<class P, class F> auto TaylorModel<P,F>::codomain() const -> CodomainType
{
    RangeType rng=this->range();
    return cast_exact_interval(convert_interval(rng,dp));
}


// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
template<class P, class F> auto TaylorModel<P,F>::range() const -> RangeType {
    const TaylorModel<P,F>& tm=*this;
    const SizeType as=tm.argument_size();
    const PrecisionType prec = tm.precision();
    const CoefficientType zero(prec);
    CoefficientType constant_term(zero);
    Array<CoefficientType> linear_terms(as,zero);
    Array<CoefficientType> quadratic_terms(as,zero);
    ErrorType err(prec);
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
if constexpr(IsSame<P,ValidatedTag>::value) {
    FloatBounds<PR> r(constant_term-err,constant_term+err);
    const FloatBounds<PR> unit_ivl(-1,+1);
    // If the ratio b/a is very large, then roundoff error can cause a significant
    // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
    for(SizeType j=0; j!=as; ++j) {
        const CoefficientType& a=quadratic_terms[j];
        const CoefficientType& b=linear_terms[j];
        FloatBounds<PR> ql=abs(a)*unit_ivl + abs(b)*unit_ivl;
        if(not is_same_as_zero(a)) { // Explicitly test for zero
            FloatBounds<PR> qf=a*(sqr(unit_ivl+b/a/2))-sqr(b)/a/4;
            r += refinement(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
        } else {
            r += ql;
        }
    }
    return RangeType(r);
} else { assert(false); }
}




//////////////////////////////////////////////////////////////////////////////

// ExactTag functions (max, min, abs, neg) and arithmetical functions (sqr, pow)

template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Max, const TaylorModel<P,F>& x, const TaylorModel<P,F>& y) {
    typedef typename TaylorModel<P,F>::RangeType RangeType;
    RangeType xr=x.range();
    RangeType yr=y.range();
    if(definitely(xr.lower()>=yr.upper())) {
        return x;
    } else if(definitely(yr.lower()>=xr.upper())) {
        return y;
    } else {
        return hlf((x+y)+abs(x-y));;
    }
}


template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Min, const TaylorModel<P,F>& x, const TaylorModel<P,F>& y) {
    typedef typename TaylorModel<P,F>::RangeType RangeType;
    RangeType xr=x.range();
    RangeType yr=y.range();
    if(definitely(xr.upper()<=yr.lower())) {
        return x;
    } else if(definitely(yr.upper()<=xr.lower())) {
        return y;
    } else {
        return hlf((x+y)-abs(x-y));
    }
}

template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Abs, const TaylorModel<P,F>& x) {
    typedef typename F::PrecisionType PR;
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    typedef typename TaylorModel<P,F>::RangeType RangeType;
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
        TaylorModel<P,F> r(x.argument_size(),x.sweeper());
        CoefficientType xmag=cast_exact(mag(xr));
        TaylorModel<P,F> s=x/xmag;
        s=sqr(s);
        r=static_cast<CoefficientType>(p[n-1u]);
        for(Nat i=0; i!=(n-1u); ++i) {
            Nat j=(n-2)-i;
            r=s*r+static_cast<CoefficientType>(p[j]);
        }
        r+=FloatBounds<PR>(-err,+err);
        return r*xmag;
    }
}

template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Max op, const TaylorModel<P,F>& x, const NumericType& c) {
    return apply(op, x, x.create_constant(c));
}
template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Min op, const TaylorModel<P,F>& x, const NumericType& c) {
    return apply(op, x, x.create_constant(c));
}
template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Max op, const NumericType& c, const TaylorModel<P,F>& x) {
    return apply(op, x.create_constant(c), x);
}
template<class P, class F> TaylorModel<P,F> AlgebraOperations<TaylorModel<P,F>>::apply(Min op, const NumericType& c, const TaylorModel<P,F>& x) {
    return apply(op, x.create_constant(c), x);
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

template<class P, class F> TaylorModel<P,F> sqr(const TaylorModel<P,F>& x) {
    return x*x;
}

template<class P, class F> TaylorModel<P,F> pow(const TaylorModel<P,F>& x, Int n) {
    TaylorModel<P,F> r=x.create_constant(1);
    TaylorModel<P,F> p(x);
    while(n) { if(n%2) { r=r*p; } p=sqr(p); n/=2; }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
template<class X> class TaylorSeries;


template<class P, class F> TaylorModel<P,F>
compose(const TaylorSeries<typename F::PrecisionType>& ts, const TaylorModel<P,F>& tv)
{
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    Sweeper<F> threshold_sweeper(new ThresholdSweeper<F>(tv.precision(),MACHINE_EPSILON));
    CoefficientType& vref=const_cast<CoefficientType&>(tv.value());
    CoefficientType vtmp=vref;
    vref=0;
    TaylorModel<P,F> r(tv.argument_size(),tv.sweeper());
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
template<class P, class F> TaylorModel<P,F>
compose(const AnalyticFunction& fn, const TaylorModel<P,F>& tm) {
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;

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
    CoefficientType c=tm.value();
    FloatDPBounds r=cast_singleton(tm.range());
    Series<FloatDPBounds> centre_series=fn.series(c);
    Series<FloatDPBounds> range_series=fn.series(r);
    //std::cerr<<"c="<<c<<"\nr="<<r<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";


    ErrorType se=mag(range_series[d]-centre_series[d]);
    ErrorType e=mag(r-c);
    ErrorType p=pow(e,d);
    //std::cerr<<"se="<<se<<"\ne="<<e<<"\np="<<p<<"\n";
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    ErrorType truncation_error=se*p;
    //std::cerr<<"truncation_error="<<truncation_error<<"\n\n";
    if(truncation_error.raw()>max_truncation_error) {
        ARIADNE_WARN("Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<max_truncation_error<<"\n");
    }

    TaylorModel<P,F> x=tm-c;
    TaylorModel<P,F> res(tm.argument_size(),tm.sweeper());
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

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::_embed_error(const TaylorModel<P,F>& tm) {
    if constexpr (IsSame<ErrorType,UnknownError<F>>::value) { return tm; }// do nothing
    else {
        const SizeType as=tm.argument_size();
        TaylorModel<P,F> rtm(as+1u,tm.sweeper());
        MultiIndex ra(as+1u);

        // The new error term is first in reverse lexicographic order.
        CoefficientType err_coef=cast_exact(tm.error());
        ra[as]=1;
        rtm._append(ra,err_coef);
        ra[as]=0;

        // Copy new terms
        for(typename TaylorModel<P,F>::ConstIterator iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
            UniformConstReference<MultiIndex> xa=iter->index();
            UniformConstReference<CoefficientType> xv=iter->coefficient();
            for(SizeType j=0; j!=as; ++j) { ra[j]=xa[j]; }
            rtm._append(ra,xv);
        }
        return std::move(rtm);
    }
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::_discard_variables(const TaylorModel<P,F>& tm, Array<SizeType> const& discarded_variables) {
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
    TaylorModel<P,F> rtm(number_of_kept_variables,tm.sweeper());
    rtm.expansion().reserve(tm.number_of_nonzeros()+1u);

    // Set the uniform error of the original model
    // If index_of_error == number_of_error_variables, then the error is kept as a uniform error bound
    MultiIndex ra(number_of_kept_variables);
    ErrorType derr=mag(tm.error()); // Magnitude of discarded terms
    for(typename TaylorModel<P,F>::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        UniformConstReference<MultiIndex> xa=iter->index();
        UniformConstReference<CoefficientType> xv=iter->coefficient();
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


template<class P, class F> Void TaylorModel<P,F>::antidifferentiate(SizeType k) {
    TaylorModel<P,F>& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());

    ErrorType e=nul(this->error());
    for(typename TaylorModel<P,F>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformReference<CoefficientType> xv=xiter->coefficient();
        xa[k]+=1;
        Nat c=xa[k];
        xv=div_err(xv,c,e);
    }
    x.error()+=e;
}

template<class P, class F> TaylorModel<P,F> antiderivative(const TaylorModel<P,F>& x, SizeType k) {
    TaylorModel<P,F> r(x);
    r.antidifferentiate(k);
    return r;
}


// Compute derivative inplace by computing term-by-term, switching the rounding mode
// Note that since some terms may be eliminated, requiring two iterators.
template<class P, class F> Void TaylorModel<P,F>::differentiate(SizeType k) {
    TaylorModel<P,F> const& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());
    // ARIADNE_PRECONDITION_MSG(x.error().raw()==0,x);
    this->clobber();

    TaylorModel<P,F>& r=*this;
    ErrorType& re=r.error();
    typename TaylorModel<P,F>::Iterator riter=r.begin();
    for(typename TaylorModel<P,F>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<CoefficientType> xv=xiter->coefficient();
        Nat c=xa[k];
        if(c!=0) {
            UniformReference<MultiIndex> ra=riter->index();
            UniformReference<CoefficientType> rv=riter->coefficient();
            ra=xa; ra[k]-=1;
            rv=mul_err(xv,c,re);
            ++riter;
        }
    }

    r.expansion().resize(static_cast<SizeType>(riter - r.begin()));
}




template<class P, class F> TaylorModel<P,F> derivative(const TaylorModel<P,F>& x, SizeType k) {
    TaylorModel<P,F> rx=x; rx.differentiate(k); return rx;

    ARIADNE_ASSERT(k<x.argument_size());
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;

    MultiIndex ra(x.argument_size()); CoefficientType rv; Nat c;

    TaylorModel<P,F> r(x.argument_size(),x.sweeper());
    typename TaylorModel<P,F>::ErrorType& re=r.error();
    for(typename TaylorModel<P,F>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<CoefficientType> xv=xiter->coefficient();
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

template<class P, class F> auto TaylorModel<P,F>::_evaluate(const TaylorModel<P,F>& tm, const Vector<ValidatedNumericType>& x) -> NumericType {
    if constexpr (IsSame<NumericType,ValidatedNumericType>::value) {
        return horner_evaluate(tm.expansion(),x)+pm(tm.error());
    } else {
        return horner_evaluate(tm.expansion(),Vector<NumericType>(x))+pm(tm.error());
    }
}

template<class P, class F> auto TaylorModel<P,F>::_evaluate(const TaylorModel<P,F>& tm, const Vector<ApproximateNumericType>& x) -> ApproximateNumericType {
    return horner_evaluate(tm.expansion(),x);
}


template<class P, class F> auto TaylorModel<P,F>::_gradient(const TaylorModel<P,F>& tm, const Vector<NumericType>& x) -> Covector<NumericType> {
    Vector<Differential<NumericType>> dx=Differential<NumericType>::variables(1u,x);
    Differential<NumericType> df=horner_evaluate(tm.expansion(),dx)+pm(tm.error());
    return gradient(df);
}



template<class P, class F> TaylorModel<P,F>
TaylorModel<P,F>::_compose(TaylorModel<P,F> const& x, Vector<TaylorModel<P,F>> const& y) {
    return horner_evaluate(x.expansion(),y)+pm(x.error());
}

template<class P, class F> Void TaylorModel<P,F>::unscale(IntervalDomainType const& codom) {
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
    TaylorModel<P,F>& tm=*this;
    Interval<CoefficientType> ivl=convert_interval(codom,this->precision());
    ARIADNE_ASSERT_MSG(decide(ivl.lower()<=ivl.upper()),"Cannot unscale TaylorModel<P,F> "<<tm<<" from empty interval "<<ivl);

    if(codom.lower()==codom.upper()) {
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

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::_compose(const Unscaling& u, const TaylorModel<P,F>& y) {
    TaylorModel<P,F> r=y; r.unscale(u.domain()); return r;
}

template<class P, class F> TaylorModel<P,F>
TaylorModel<P,F>::_compose(const TaylorModel<P,F>& x, const VectorUnscaling& u, const Vector<TaylorModel<P,F>>& y) {
    return compose(x,compose(u,y));
}



template<class T> class Powers {
  public:
    explicit Powers(const T& t) {
        _values.push_back(t*0+1); _values.push_back(t); }
    explicit Powers(const T& z, const T& t) {
        _values.push_back(z); _values.push_back(t); }
    const T& operator[](SizeType i) const {
        while(_values.size()<=i) { _values.push_back(_values[1]*_values.back()); } return _values[i]; }
  private:
    mutable std::vector<T> _values;
};

template<class F> class Powers<Bounds<F>> {
    typedef Bounds<F> T;
  public:
    explicit Powers(const T& t) {
        _values.push_back(T(1,t.precision())); _values.push_back(t); }
    explicit Powers(const T& z, const T& t) { _values.push_back(z); _values.push_back(t); }
    const T& operator[](SizeType i) const {
        while(_values.size()<=i) {
            if(_values.size()%2==0) { _values.push_back(sqr(_values[_values.size()/2])); }
            else { _values.push_back(_values[1]*_values.back()); } }
        return _values[i]; }
  private:
    mutable std::vector<T> _values;
};



template<class P, class F> TaylorModel<P,F>
TaylorModel<P,F>::_partial_evaluate(const TaylorModel<P,F>& x, SizeType k, NumericType c)
{
    const SizeType as=x.argument_size();
    Vector<TaylorModel<P,F>> y(as,TaylorModel<P,F>(as-1u,x.sweeper()));
    for(SizeType i=0; i!=k; ++i) { y[i]=TaylorModel<P,F>::coordinate(as-1u,i,x.sweeper()); }
    y[k]=TaylorModel<P,F>::constant(as-1u,c,x.sweeper());
    for(SizeType i=k+1; i!=as; ++i) { y[i]=TaylorModel<P,F>::coordinate(as-1u,i-1u,x.sweeper()); }
    return compose(x,y);
}



template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::_embed(SizeType as1, const TaylorModel<P,F>& tm2, SizeType as3) {
    return TaylorModel<P,F>(embed(as1,tm2.expansion(),as3),tm2.error(),tm2.sweeper());
}


template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::_split(const TaylorModel<P,F>& tm, SizeType k, SplitPart h) {
    const DegreeType deg=tm.degree();
    const SizeType as=tm.argument_size();
    SweeperType swp=tm.sweeper();

    TaylorModel<P,F> r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(typename TaylorModel<P,F>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const DegreeType ak=iter->index()[k];
        UniformReference<CoefficientType> c=iter->coefficient();
        c/=pow(two,ak);
    }

    if(h==SplitPart::MIDDLE) { return r; }
    Int tr=( h==SplitPart::UPPER ? +1 : -1 );

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<TaylorModel<P,F>> ary(deg+1u,TaylorModel<P,F>(as,swp));
    for(typename TaylorModel<P,F>::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
        MultiIndex a=iter->index();
        UniformConstReference<CoefficientType> c=iter->coefficient();
        DegreeType ak=a[k];
        a[k]=0u;
        ary[ak].expansion().append(a,c);
    }

    ErrorType re=r.error();
    r.clear();
    r.set_error(re);

    for(DegreeType i=0; i<=deg; ++i) {
        for(DegreeType j=i; j<=deg; ++j) {
            Int sf=bin(j,i);
            if(tr==-1 && (j-i)%2==1) { sf=-sf; }
            r+=ary[j]*sf;
            for(typename TaylorModel<P,F>::Iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->index()[k];
            }
         }
    }

    return r;
}



///////////////////////////////////////////////////////////////////////////////

// Banach algebra operations


template<class P, class F> auto TaylorModel<P,F>::average() const -> CoefficientType {
    return (*this)[MultiIndex::zero(this->argument_size())];
}

template<class P, class F> auto TaylorModel<P,F>::radius() const -> NormType {
    typename TaylorModel<P,F>::NormType r(0u,this->precision());
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->index().degree()!=0) {
            r+=mag(iter->coefficient());
        }
    }
    r+=this->error();
    return r;
}

template<class P, class F> auto TaylorModel<P,F>::norm() const -> NormType {
    NormType r(0u,this->precision());
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=mag(iter->coefficient());
    }
    r+=this->error();
    return r;
}

template<class P, class F> typename TaylorModel<P,F>::NormType norm(const TaylorModel<P,F>& tm) {
    return tm.norm();
}


template<class P, class F> Bool TaylorModel<P,F>::_refines(const TaylorModel<P,F>& tm1, const TaylorModel<P,F>& tm2)
{
    ARIADNE_ASSERT(tm1.argument_size()==tm2.argument_size());
    TaylorModel<P,F> d=tm2;
    d.error()=0u;
    d-=tm1;
    return d.norm().raw() <= tm2.error().raw();
}


template<class P, class F> Bool TaylorModel<P,F>::_consistent(const TaylorModel<P,F>& tm1, const TaylorModel<P,F>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (Ariadne::norm(tm1-tm2).raw() <= (tm1.error()+tm2.error()).raw()*2u);
}

template<class P, class F> Bool TaylorModel<P,F>::_inconsistent(const TaylorModel<P,F>& tm1, const TaylorModel<P,F>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (mag(tm1.value()-tm2.value()).raw() > (tm1.error()+tm2.error()).raw());
}

template<class P, class F> TaylorModel<P,F> TaylorModel<P,F>::_refinement(const TaylorModel<P,F>& x, const TaylorModel<P,F>& y) {
if constexpr (IsSame<P,ValidatedTag>::value) {
    TaylorModel<P,F> r(x.argument_size(),x.sweeper());

    ErrorType max_error=nul(r.error());

    const ErrorType& xe=x.error();
    const ErrorType& ye=y.error();
    CoefficientType rv,xv,yv;
    FloatDP xu,yu,mxl,myl,u,ml;
    MultiIndex a;

    typename TaylorModel<P,F>::ConstIterator xiter=x.begin();
    typename TaylorModel<P,F>::ConstIterator yiter=y.begin();
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
            ARIADNE_THROW(IntersectionException,"refinement(TaylorModel<P,F>,TaylorModel<P,F>)",x<<" and "<<y<<" are inconsistent.");
        }

        if(rve.value()!=0) { r.expansion().append(a,rve.value()); }
        max_error=max(max_error,rve.error());
    }

    r.error()=max_error;

    return r;
} else { ARIADNE_NOT_IMPLEMENTED; }
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

template<class P, class F> OutputStream& TaylorModel<P,F>::str(OutputStream& os) const {
    TaylorModel<P,F> const& tm=*this;

    // Set the variable names to be 'parameter' s0,s1,..
    Array<StringType> variable_names(tm.argument_size());
    for(SizeType j=0; j!=tm.argument_size(); ++j) {
        StringStream sstr;
        sstr << 's' << j;
        variable_names[j]=sstr.str();
    }

    //os << "TaylorModel<P,F>";
    os << "TM["<<tm.argument_size()<<"](";
    Expansion<MultiIndex,CoefficientType> e=tm.expansion();
    e.sort(GradedIndexLess());
    e.write(os,variable_names);
    return os << "+/-" << tm.error() << ")";
}

template<class P, class F> OutputStream& TaylorModel<P,F>::repr(OutputStream& os) const {
    return this->str(os);
}

///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


template<class P, class F> Vector<TaylorModel<P,F>> TaylorModel<P,F>::zeros(SizeType rs, SizeType as, SweeperType swp)
{
    Vector<TaylorModel<P,F>> result(rs,TaylorModel<P,F>::zero(as,swp));
    return result;
}


template<class P, class F> Vector<TaylorModel<P,F>> TaylorModel<P,F>::constants(SizeType as, const Vector<ValidatedNumericType>& c, SweeperType swp)
{
    Vector<TaylorModel<P,F>> result(c.size(),TaylorModel<P,F>::zero(as,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel<P,F>::constant(as,c[i],swp);
    }
    return result;
}

template<class P, class F> Vector<TaylorModel<P,F>> TaylorModel<P,F>::coordinates(SizeType as, SweeperType swp)
{
    Vector<TaylorModel<P,F>> result(as,TaylorModel<P,F>::zero(as,swp));
    for(SizeType i=0; i!=as; ++i) { result[i]=TaylorModel<P,F>::coordinate(as,i,swp); }
    return result;
}

template<class P, class F> Vector<TaylorModel<P,F>> TaylorModel<P,F>::scalings(const BoxDomainType& d, SweeperType swp)
{
    Vector<TaylorModel<P,F>> result(d.size(),TaylorModel<P,F>::zero(d.size(),swp));
    for(SizeType i=0; i!=d.size(); ++i) {
        result[i]=TaylorModel<P,F>::scaling(d.size(),i,d[i],swp);
    }
    return result;
}



///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
template<class F> Matrix<Bounds<F>>
jacobian(const Vector<TaylorModel<ValidatedTag,F>>& f, const Vector<Bounds<F>>& x) {
    Vector< Differential<Bounds<F>> > dx=Differential<Bounds<F>>::variables(1u,x);
    Vector< Differential<Bounds<F>> > df(f.size(),x.size(),1u);
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<Bounds<F>> J=jacobian(df);
    return J;
}

// Compute the Jacobian over an arbitrary domain
template<class F> Matrix<Bounds<F>>
jacobian(const Vector<TaylorModel<ValidatedTag,F>>& f, const Vector<Bounds<F>>& x, const Array<SizeType>& p) {
    Vector<Differential<Bounds<F>>> dx(x.size(),x.size(),1u);
    for(SizeType j=0; j!=x.size(); ++j) {
        dx[j]=Differential<Bounds<F>>::constant(p.size(),1u,x[j]); }
    for(SizeType k=0; k!=p.size(); ++k) {
        SizeType j=p[k];
        dx[j]=Differential<Bounds<F>>::variable(p.size(),1u,x[j],k); }
    Vector< Differential<Bounds<F>> > df(f.size());
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<Bounds<F>> J=jacobian(df);
    return J;
}

// Compute the Jacobian at the origin
template<class P, class F> Matrix<Value<F>>
jacobian_value(const Vector<TaylorModel<P,F>>& f) {
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<CoefficientType> J(rs,as);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType j=0; j!=as; ++j) {
            a[j]=1; const CoefficientType x=f[i][a]; J[i][j]=x; a[j]=0;
        }
    }
    return J;
}

// Compute the Jacobian at the origin with respect to the variables args.
template<class P, class F> Matrix<Value<F>>
jacobian_value(const Vector<TaylorModel<P,F>>& f, const Array<SizeType>& p) {
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    const SizeType rs=f.size();
    const SizeType as=f.zero_element().argument_size();
    const SizeType ps=p.size();
    Matrix<CoefficientType> J(rs,ps);
    MultiIndex a(as);
    for(SizeType i=0; i!=rs; ++i) {
        for(SizeType k=0; k!=ps; ++k) {
            SizeType j=p[k]; a[j]=1; const CoefficientType x=f[i][a]; J[i][k]=x; a[j]=0;
        }
    }
    return J;
}



// Compute the Jacobian over the unit domain
template<class P, class F> Matrix<typename TaylorModel<P,F>::RangeType>
jacobian_range(const Vector<TaylorModel<P,F>>& f) {
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using RangeType = typename TaylorModel<P,F>::RangeType;
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<RangeType> J(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<P,F>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            for(SizeType k=0; k!=as; ++k) {
                const Nat c=a[k];
                if(c>0) {
                    UniformConstReference<CoefficientType> x=iter->coefficient();
                    if(a.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=RangeType(-1,1)*x*c; }
                }
            }
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain, with respect to the variables p.
template<class P, class F> Matrix<typename TaylorModel<P,F>::RangeType>
jacobian_range(const Vector<TaylorModel<P,F>>& f, const Array<SizeType>& p) {
    using CoefficientType = typename TaylorModel<P,F>::CoefficientType;
    using RangeType = typename TaylorModel<P,F>::RangeType;
    SizeType rs=f.size();
    SizeType ps=p.size();
    Matrix<RangeType> J(rs,ps);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<P,F>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            for(SizeType k=0; k!=ps; ++k) {
                SizeType j=p[k];
                const Nat c=a[j];
                if(c>0) {
                    UniformConstReference<CoefficientType> x=iter->coefficient();
                    if(a.degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=RangeType(-1,1)*x*c; }
                }
            }
        }
    }
    return J;
}


template<class F> TaylorModel<ValidatedTag,F> value_coefficients(TaylorModel<ValidatedTag,Bounds<F>> const& tm) {
    TaylorModel<ValidatedTag,F> r(tm.argument_size(),tm.sweeper());
    auto& e=r.error();

    for (auto term : tm.expansion()) {
        Value<F> c=set_err(term.coefficient(),e);
        r.expansion().append(term.index(), c);
    }
    e+=tm.error();
    return r;
}

template<class F> TaylorModel<ValidatedTag,Bounds<F>> exact_coefficients(TaylorModel<ValidatedTag,Bounds<F>> const& tm) {
    TaylorModel<ValidatedTag,Bounds<F>> r(tm.argument_size(),tm.sweeper());
    auto& e=r.error();

    for (auto term : tm.expansion()) {
        Value<F> c=set_err(term.coefficient(),e);
        r.expansion().append(term.index(), c);
    }
    e+=tm.error();
    return r;
}




} //namespace Ariadne


