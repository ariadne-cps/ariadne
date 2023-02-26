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

#include "numeric/numeric.hpp"

#include <iomanip>
#include <limits>

#include "numeric/rounding.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/expansion.hpp"
#include "algebra/series.hpp"
#include "algebra/differential.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_series.hpp"
#include "function/function.hpp"
#include "utility/exceptions.hpp"

#include "algebra/expansion.inl.hpp"
#include "algebra/evaluate.tpl.hpp"
#include "algebra/algebra_operations.tpl.hpp"

#define VOLATILE ;
#include "algebra/multi_index-noaliasing.hpp"
#include "function/function_mixin.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

typedef IntegralConstant<Int,0> Zero;
typedef IntegralConstant<Int,1> One;
static const Zero zero = Zero();
static const One one = One();


template<class FLT> FLT UnknownError<FLT>::raw() const {
    return FLT(0u,this->precision()); }
template<class FLT> typename FLT::PrecisionType UnknownError<FLT>::precision() const {
    if constexpr (Same<FLT,FloatMP>) { return MultiplePrecision(64_bits); }
    else { return typename FLT::PrecisionType(); } }
template<class FLT> UnknownError<FLT>::operator PositiveApproximation<FLT> () const {
    return PositiveApproximation<FLT>(0u,this->precision()); }

template<class FLT> UnknownError<FLT> nul(UnknownError<FLT>) { return UnknownError<FLT>(); }
template<class FLT> UnknownError<FLT> mag(UnknownError<FLT>) { return UnknownError<FLT>(); }
template<class FLT> UnknownError<FLT> operator+(UnknownError<FLT>, UnknownError<FLT>) { return UnknownError<FLT>(); }
template<class FLT> UnknownError<FLT> operator*(UnknownError<FLT>, UnknownError<FLT>) { return UnknownError<FLT>(); }
template<class FLT> UnknownError<FLT>& operator+=(UnknownError<FLT>& e, UnknownError<FLT> const&) { return e; }
template<class FLT> UnknownError<FLT>& operator*=(UnknownError<FLT>& e, UnknownError<FLT> const&) { return e; }

template<class FLT> UnknownError<FLT> operator+(UnknownError<FLT>, PositiveApproximation<FLT> const&) { return UnknownError<FLT>(); }
template<class FLT> UnknownError<FLT>& operator+=(UnknownError<FLT>& e, PositiveApproximation<FLT> const&) { return e; }
template<class FLT> UnknownError<FLT>& operator*=(UnknownError<FLT>& e, PositiveApproximation<FLT> const&) { return e; }
template<class FLT> OutputStream& operator<<(OutputStream& os, UnknownError<FLT> const& e) { return os << "???"; }
template<class FLT> UnknownError<FLT>& operator+=(UnknownError<FLT>& e, ExactDouble const&) { return e; }


template<class FLT> ZeroError<FLT> const& nul(ZeroError<FLT> const& e) { return e; }
template<class FLT> PositiveUpperBound<FLT> const& operator+(PositiveUpperBound<FLT> const& r, ZeroError<FLT> const&) { return r; }
template<class FLT> PositiveUpperBound<FLT>& operator+=(PositiveUpperBound<FLT>&, ZeroError<FLT> const&) { }
template<class FLT> ZeroError<FLT>& operator+=(ZeroError<FLT>&, ZeroError<FLT> const&) { }


template<class FLT> Bounds<FLT> fma(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Bounds<FLT> y) { return x1*x2+y; }
template<class FLT> Bounds<FLT>& operator/=(Bounds<FLT>& x, TwoExp y) { return x*=Dyadic(rec(y)); }

template<class FLT> Approximation<FLT>& operator/=(Approximation<FLT>& x, TwoExp y) { return x/=Dyadic(y); }
template<class FLT> Approximation<FLT> fma(Approximation<FLT>const& x, Approximation<FLT> const& y, Approximation<FLT> z) {
    z._a = fma(near,x._a,y._a,z._a);; return std::move(z); }
template<class FLT> Approximation<FLT> pm(UnknownError<FLT> e) { return Approximation<FLT>(e.precision()); }


template<class FLT> UpperInterval<FLT>& operator/=(UpperInterval<FLT>& ivl1, ValidatedNumber const& y2) {
    return ivl1=ivl1/y2; }
template<class FLT> UpperInterval<FLT>& operator/=(UpperInterval<FLT>& ivl1, Dyadic const& w2) {
    return ivl1=ivl1/w2; }

//namespace {
// Operations for Interval<UpperBound<FLT>> coefficients
template<class X> X const& cast_singleton(X const& x) { return x; }
template<class FLT> Bounds<FLT> const& cast_singleton(Interval<UpperBound<FLT>> const& ivl) { return reinterpret_cast<Bounds<FLT>const&>(ivl); }

template<class FLT> ApproximateInterval<FLT> operator+(ApproximateInterval<FLT> const& x1, ApproximateInterval<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1.lower_bound()+x2.lower_bound(),x1.upper_bound()+x2.upper_bound()); }
template<class FLT> ApproximateInterval<FLT> operator-(ApproximateInterval<FLT> const& x1, ApproximateInterval<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1.lower_bound()-x2.upper_bound(),x1.upper_bound()-x2.lower_bound()); }
template<class FLT> ApproximateInterval<FLT> operator*(ApproximateInterval<FLT> const& x1, ApproximateInterval<FLT> const& x2) {
    return make_interval(reinterpret_cast<Bounds<FLT>const&>(x1)*reinterpret_cast<Bounds<FLT>const&>(x2)); }

template<class FLT> ApproximateInterval<FLT> operator+(ApproximateInterval<FLT> x1, Approximation<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1.lower_bound()+x2,x1.upper_bound()+x2); }
template<class FLT> ApproximateInterval<FLT> operator-(ApproximateInterval<FLT> x1, Approximation<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1.lower_bound()-x2,x1.upper_bound()-x2); }
template<class FLT> ApproximateInterval<FLT> operator*(ApproximateInterval<FLT> x1, Approximation<FLT> const& x2) {
    if (x2.raw()>=0) { return ApproximateInterval<FLT>(x1.lower_bound()*x2,x1.upper_bound()*x2); }
    else { return ApproximateInterval<FLT>(x1.upper_bound()*x2,x1.lower_bound()*x2); }
}
template<class FLT> ApproximateInterval<FLT> operator+(Approximation<FLT> const& x1, ApproximateInterval<FLT> x2) {
    return ApproximateInterval<FLT>(x1+x2.lower_bound(),x1+x2.upper_bound()); }
template<class FLT> ApproximateInterval<FLT> operator-(Approximation<FLT> const& x1, ApproximateInterval<FLT> x2) {
    return ApproximateInterval<FLT>(x1-x2.upper_bound(),x1-x2.lower_bound()); }
template<class FLT> ApproximateInterval<FLT> operator*(Approximation<FLT> const& x1, ApproximateInterval<FLT> x2) {
    return x2*x1; }
/*
template<class FLT> ApproximateInterval<FLT> operator+(UpperInterval<FLT> x1, Approximation<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1.lower_bound()+x2,x1.upper_bound()+x2); }
template<class FLT> ApproximateInterval<FLT> operator-(UpperInterval<FLT> x1, Approximation<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1.lower_bound()-x2,x1.upper_bound()-x2); }
*/
template<class FLT> ApproximateInterval<FLT> operator*(UpperInterval<FLT> const& x1, Approximation<FLT> const& x2) {
    if (x2.raw()>=0) { return ApproximateInterval<FLT>(x1.lower_bound()*x2,x1.upper_bound()*x2); }
    else { return ApproximateInterval<FLT>(x1.upper_bound()*x2,x1.lower_bound()*x2); } }
/*
template<class FLT> ApproximateInterval<FLT> operator+(Approximation<FLT> const& x1, UpperInterval<FLT> x2) {
    return ApproximateInterval<FLT>(x1+x2.lower_bound(),x1+x2.upper_bound()); }
template<class FLT> ApproximateInterval<FLT> operator-(Approximation<FLT> const& x1, UpperInterval<FLT> x2) {
    return ApproximateInterval<FLT>(x1-x2.upper_bound(),x1-x2.lower_bound()); }
*/
template<class FLT> ApproximateInterval<FLT> operator*(Approximation<FLT> const& x1, UpperInterval<FLT> const& x2) {
    return x2*x1; }

template<class FLT> ApproximateInterval<FLT> operator*(UpperInterval<FLT> const& x1, ApproximateInterval<FLT> const& x2) {
    return ApproximateInterval<FLT>(x1)*x2; }
template<class FLT> ApproximateInterval<FLT> operator*(ApproximateInterval<FLT> const& x1, UpperInterval<FLT> const& x2) {
    return x1*ApproximateInterval<FLT>(x2); }


//} // namespace


template<class FLT> PositiveUpperBound<FLT> mag(Interval<UpperBound<FLT>> const& ivl) {
    return cast_positive(max(-ivl.lower_bound(),ivl.upper_bound()));
}
template<class FLT> PositiveApproximation<FLT> mag(Interval<Approximation<FLT>> const& ivl) {
    return cast_positive(max(-ivl.lower_bound(),ivl.upper_bound()));
}


template<class FLT> ValidatedIntervalTaylorModel<FLT> max(ValidatedIntervalTaylorModel<FLT> const& x1, Bounds<FLT> const& x2);
template<class FLT> ValidatedIntervalTaylorModel<FLT> max(Bounds<FLT> const& x1, ValidatedIntervalTaylorModel<FLT> const& x2);
template<class FLT> ValidatedIntervalTaylorModel<FLT> min(ValidatedIntervalTaylorModel<FLT> const& x1, Bounds<FLT> const& x2);
template<class FLT> ValidatedIntervalTaylorModel<FLT> max(Bounds<FLT> const& x1, ValidatedIntervalTaylorModel<FLT> const& x2);

struct UnitMultiIndex {
    SizeType n; DegreeType j;
    friend MultiIndex& operator+=(MultiIndex& a, UnitMultiIndex const& aj) { a[aj.j]+=1u; return a; }
};
struct SumMultiIndex {
    MultiIndex const& a1; MultiIndex const& a2;
    SizeType size() const { return a1.size(); }
    DegreeType operator[] (SizeType j) const { return a1[j]+a2[j]; }
    DegreeType degree() const { return a1.degree()+a2.degree(); }
};


namespace {

static const ExactDouble MACHINE_EPSILON(2.2204460492503131e-16);

Bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    return reverse_lexicographic_less(a1,a2); }


inline Interval<FloatDP> convert_interval(Interval<FloatDP> const& ivl, DoublePrecision pr) {
    return ivl; }
inline Interval<FloatMP> convert_interval(Interval<FloatDP> const& ivl, MultiplePrecision pr) {
    return Interval<FloatMP>(FloatMP(ivl.lower_bound().raw(),pr),FloatMP(ivl.upper_bound().raw(),pr)); }

inline Interval<FloatDPUpperBound> const& convert_interval(Interval<FloatDPUpperBound> const& ivl, DoublePrecision) {
    return ivl; }
inline Interval<FloatDPUpperBound> convert_interval(Interval<FloatMPUpperBound> const& ivl, DoublePrecision pr) {
    return Interval<FloatDPUpperBound>(FloatDP(ivl.lower_bound().raw(),down,pr),FloatDP(ivl.upper_bound().raw(),up,pr)); }

inline Interval<FloatDPApproximation> const& convert_interval(Interval<FloatDPApproximation> const& ivl, DoublePrecision) {
    return ivl; }
inline Interval<FloatDPApproximation> convert_interval(Interval<FloatMPApproximation> const& ivl, DoublePrecision pr) {
    return Interval<FloatDPApproximation>(FloatDP(ivl.lower_bound().raw(),near,pr),FloatDP(ivl.upper_bound().raw(),near,pr)); }

template<class FLT> Bool is_same_as_zero(Approximation<FLT> const& xa) { return xa.raw()==0; }
template<class FLT> Bool is_same_as_zero(Bounds<FLT> const& xb) { return xb.lower_raw()==0 && xb.upper_raw()==0; }
template<class FLT> Bool is_same_as_zero(FLT const& xv) { return xv==0; }
template<class FLT> Bool is_same_as_zero(UpperInterval<FLT> const& xv) { return xv.lower_bound().raw()==0 && xv.upper_bound().raw()==0; }

} // namespace



template<class FLT>
Void SweeperBase<FLT>::_sweep(Expansion<MultiIndex,Float<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,Float<PR>>::Iterator curr=p.begin();

    // FIXME: Not needed, but added to pair with rounding mode change below
    FLT::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(adv->index(),adv->coefficient())) {
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
    FLT::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class FLT>
Void SweeperBase<FLT>::_sweep(Expansion<MultiIndex,FloatBounds<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatBounds<PR>>::Iterator curr=p.begin();

    // FIXME: Not needed, but added to pair with rounding mode change below
    FLT::set_rounding_upward();
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
    FLT::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class FLT>
Void SweeperBase<FLT>::_sweep(Expansion<MultiIndex,FloatApproximation<PR>>& p) const
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


template<class FLT>
Void SweeperBase<FLT>::_sweep(Expansion<MultiIndex,FloatUpperInterval<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::Iterator curr=p.begin();

    // FIXME: Not needed, but added to pair with rounding mode change below
    FLT::set_rounding_upward();
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
    FLT::set_rounding_to_nearest();
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


template<class FLT>
Void RelativeSweeperBase<FLT>::_sweep(Expansion<MultiIndex,Float<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,Float<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,Float<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=radius(p)+e;

    // FIXME: Not needed, but added to pair with rounding mode change below
    FLT::set_rounding_upward();
    FloatError<PR> te(e.precision());
    while(adv!=end) {
        if(this->_discard(adv->coefficient(),nrm.raw())) {
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
    FLT::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class FLT>
Void RelativeSweeperBase<FLT>::_sweep(Expansion<MultiIndex,FloatBounds<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatBounds<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatBounds<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=radius(p)+e;

    // FIXME: Not needed, but added to pair with rounding mode change below
    FLT::set_rounding_upward();
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
    FLT::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}


template<class FLT>
Void RelativeSweeperBase<FLT>::_sweep(Expansion<MultiIndex,FloatApproximation<PR>>& p) const
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


template<class FLT>
Void RelativeSweeperBase<FLT>::_sweep(Expansion<MultiIndex,FloatUpperInterval<PR>>& p, FloatError<PR>& e) const
{
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator end=p.end();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::ConstIterator adv=p.begin();
    typename Expansion<MultiIndex,FloatUpperInterval<PR>>::Iterator curr=p.begin();

    FloatError<PR> nrm=radius(p)+e;

    // FIXME: Not needed, but added to pair with rounding mode change below
    FLT::set_rounding_upward();
    FloatError<PR> te(0u,e.precision());
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
    FLT::set_rounding_to_nearest();
    p.resize(static_cast<SizeType>(curr-p.begin()));
}

template<class FLT> FLT create_default() {
    if constexpr (Same<FLT,FloatDP>) { return FloatDP(dp); }
    else if constexpr (Same<FLT,FloatDP>) { return FloatDP(dp); }
    else if constexpr (Same<FLT,FloatDPError>) { return FloatDPError(dp); }
    else if constexpr (Same<FLT,FloatMP>) { return FloatMP(FloatMP::get_default_precision()); }
    else if constexpr (Same<FLT,FloatMP>) { return FloatMP(FloatMP::get_default_precision()); }
    else if constexpr (Same<FLT,FloatMPError>) { return FloatMPError(FloatMP::get_default_precision()); }
    else if constexpr (Same<FLT,Interval<FloatDPUpperBound>>) { return FloatDPUpperInterval(0,0); }
    else { abort(); }
}


template<class P, class FLT> TaylorModel<P,FLT>::TaylorModel()
    : _expansion(0,create_default<CoefficientType>()), _error(create_default<ErrorType>()), _sweeper()
{
}

template<class P, class FLT> TaylorModel<P,FLT>::TaylorModel(SizeType as, SweeperType swp)
    : _expansion(as,CoefficientType(0,swp.precision())), _error(swp.precision()), _sweeper(swp)
{
}

template<class P, class FLT> TaylorModel<P,FLT>::TaylorModel(const Expansion<MultiIndex,CoefficientType>& f, const ErrorType& e, SweeperType swp)
    : _expansion(f), _error(e), _sweeper(swp)
{
    this->cleanup();
}

template<class P, class FLT> TaylorModel<P,FLT>::TaylorModel(const Expansion<MultiIndex,RawFloatType>& f, const RawFloatType& e, SweeperType swp)
//    : TaylorModel(reinterpret_cast<Expansion<MultiIndex,CoefficientType>const&>(f),reinterpret_cast<ErrorType const&>(e),swp)
    : TaylorModel(static_cast<Expansion<MultiIndex,CoefficientType>>(f),static_cast<ErrorType>(e),swp)
{
}

template<class P, class FLT> TaylorModel<P,FLT>::TaylorModel(const Expansion<MultiIndex,ExactDouble>& f, const ExactDouble& e, SweeperType swp)
    : TaylorModel(Expansion<MultiIndex,CoefficientType>(Expansion<MultiIndex,Dyadic>(f),swp.precision()),ErrorType(Dyadic(e),swp.precision()),swp)
{
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::scaling(SizeType as, SizeType j, const IntervalDomainType& codom, SweeperType swp) {
    TaylorModel<P,FLT> r(as,swp);
    auto ivl=convert_interval(codom,r.precision());
    r.set_gradient(j,1);
    r*=ivl.radius();
    r+=static_cast<NumericType>(ivl.midpoint());
    return r;
}

template<class P, class FLT> auto
TaylorModel<P,FLT>::constant(SizeType as, const NumericType& c, SweeperType swp) -> TaylorModel<P,FLT> {
    TaylorModel<P,FLT> r(as,swp); r.set_value(CoefficientType(1,r.precision())); r*=c; return r;
}
template<class P, class FLT> auto
TaylorModel<P,FLT>::constant(SizeType as, const GenericNumericType& c, SweeperType swp) -> TaylorModel<P,FLT> {
    TaylorModel<P,FLT> r(as,swp); r.set_value(CoefficientType(1,r.precision())); r*=c; return r;
}



template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::create() const {
    return TaylorModel<P,FLT>(this->argument_size(),this->_sweeper);
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::create_zero() const {
    return TaylorModel<P,FLT>(this->argument_size(),this->_sweeper);
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::create_constant(NumericType c) const {
    return TaylorModel<P,FLT>::constant(this->argument_size(),c,this->_sweeper);
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::create_constant(GenericNumericType c) const {
    return this->create_constant(NumericType(c,this->precision()));
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::create_coordinate(SizeType j) const {
    ARIADNE_PRECONDITION(j<this->argument_size());
    TaylorModel<P,FLT> r(this->argument_size(),this->_sweeper);
    CoefficientType o(1,this->precision());
    r._expansion.append(MultiIndex::unit(this->argument_size(),j),o);
    return r;
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::create_ball(ErrorType e) const {
    ARIADNE_DEBUG_PRECONDITION(e.raw()>=0);
    TaylorModel<P,FLT> r(this->argument_size(),this->_sweeper);
    r._error=e;
    return r;
}

template<class P, class FLT> Void TaylorModel<P,FLT>::swap(TaylorModel<P,FLT>& tm) {
    this->_expansion.swap(tm._expansion);
    std::swap(this->_error,tm._error);
    std::swap(this->_sweeper,tm._sweeper);
}

template<class P, class FLT> Void TaylorModel<P,FLT>::clear() {
    this->_expansion.clear();
    this->_error=0u;
}

template<class P, class FLT> DegreeType TaylorModel<P,FLT>::degree() const {
    DegreeType deg=0u;
    for(auto iter=this->begin(); iter!=this->end(); ++iter) {
        deg=std::max(deg,iter->index().degree());
    }
    return deg;
}


template<class FLT> UpperInterval<FLT> set_err(UpperInterval<FLT> const& x, Error<FLT>& e) {
    return x;
}

template<class FLT> FLT set_err(Bounds<FLT> const& x, Error<FLT>& e) {
     e+=x.error(); return x.value();
}

template<class FLT> Approximation<FLT> const& set_err(Approximation<FLT> const& x, UnknownError<FLT>& e) {
    return x;
}


template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::operator=(const NumericType& c) {
    this->clear();
    CoefficientType m=set_err(c,this->_error);
    if(not is_same_as_zero(m)) {
        this->_expansion.append(MultiIndex::zero(this->argument_size()),m);
    }
    return *this;
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::operator=(const GenericNumericType& c) {
    return *this = NumericType(c,this->precision());
}

namespace { // Internal code for arithmetic

template<class FLT> struct ValidatedApproximation {
    Bounds<FLT> _v; Approximation<FLT> _a;
    ValidatedApproximation(Bounds<FLT>const& x) : _v(x), _a(x) { }
    operator Bounds<FLT> const& () const { return _v; }
    LowerBound<FLT> lower_bound() const { return _v.lower_bound(); }
    Approximation<FLT> middle() const { return _a; }
    UpperBound<FLT> upper_bound() const { return _v.upper_bound(); }
    FLT const& lower_raw() const { return _v.lower_raw(); }
    FLT const& middle_raw() const { return _a.raw(); }
    FLT const& upper_raw() const { return _v.upper_raw(); }
    friend OutputStream& operator<<(OutputStream& os, ValidatedApproximation<FLT> const& x) {
        return os << "{"<<x._v.lower_bound()<<":"<<x._a<<":"<<x._v.upper_bound()<<"}"; }
};
template<class FLT> Approximation<FLT> const& make_validated_approximation(Approximation<FLT> const& x) { return x; }
template<class FLT> ValidatedApproximation<FLT> make_validated_approximation(Bounds<FLT> const& x) { return ValidatedApproximation<FLT>(x); }



template<ARawFloat FLT> Rounded<FLT> const& cast_rounded(FLT const& x) { return reinterpret_cast<Rounded<FLT>const&>(x); }
template<ARawFloat FLT> Rounded<FLT>& cast_rounded(FLT& x) { return reinterpret_cast<Rounded<FLT>&>(x); }
template<ARawFloat FLT> Rounded<FLT>& cast_rounded(Error<FLT>& x) { return reinterpret_cast<Rounded<FLT>&>(x); }

template<class FLT, class PRE> Ball<FLT,RawFloatType<PRE>> add(FLT const& x1, FLT const& x2, PRE pre) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(x1.raw() + x2.raw());
    FLT::set_rounding_upward();
    FLT u=x1.raw()+x2.raw();
    FLT ml=mx1.raw()-x2.raw();
    Error<RawFloatType<PRE>> e(pre);
    e += max(u-r,ml+r);
    return Ball(r,e);
}
template<class FLT, class PRE> Ball<FLT,RawFloatType<PRE>> mul(FLT const& x1, FLT const& x2, PRE pre) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(x1.raw() * x2.raw());
    FLT::set_rounding_upward();
    FLT u=x1.raw()*x2.raw();
    FLT ml=mx1.raw()*x2.raw();
    Error<RawFloatType<PRE>> e(pre);
    e += max(u-r,ml+r);
    return Ball(r,e);
}



template<class FLT> FLT add_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    Rounded<FLT> mx1=-x1;
    FLT::set_rounding_to_nearest();
    Rounded<FLT> r(cast_rounded(x1) + cast_rounded(x2));
    FLT::set_rounding_upward();
    Rounded<FLT> u=cast_rounded(x1)+cast_rounded(x2);
    Rounded<FLT> ml=mx1-cast_rounded(x2);
    cast_rounded(e) += hlf(u+ml);
    return cast_exact(r);
}

template<class FLT> FLT add_err(FLT const& x, ValidatedApproximation<FLT> const& c, Error<FLT>& e) {
    Rounded<FLT> const& xv=x.raw();
    Rounded<FLT> const& cl=c.lower_raw();
    Rounded<FLT> const& cm=c.middle_raw();
    Rounded<FLT> const& cu=c.upper_raw();
    Rounded<FLT>& re=cast_rounded(e);
    FLT::set_rounding_to_nearest();
    Rounded<FLT> rv=xv+cm;
    FLT::set_rounding_upward();
    Rounded<FLT> u=xv+cu;
    Rounded<FLT> ml=(-xv)-cl;
    re += hlf(u+ml);
    return FLT(rv);
}

template<class FLT> FLT add_err(FLT const& x, Bounds<FLT> const& c, Error<FLT>& e) {
    return add_err(x,ValidatedApproximation<FLT>(c),e);
}

template<class FLT> Bounds<FLT> add_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return add(x1,x2);
}

template<class FLT> UpperInterval<FLT> add_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return add(x1,x2);
}

template<class FLT> UpperInterval<FLT> add_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return add(x1,n2);
}

template<class FLT> Approximation<FLT> add_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return add(x1,x2);
}

template<class FLT> FLT sub_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    Rounded<FLT> mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(cast_rounded(x1) - cast_rounded(x2));
    FLT::set_rounding_upward();
    Rounded<FLT> u=cast_rounded(x1)-cast_rounded(x2);
    Rounded<FLT> ml=mx1+cast_rounded(x2);
    cast_rounded(e) += hlf(u+ml);
    return r;
}

template<class FLT> Bounds<FLT> sub_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return sub(x1,x2);
}

template<class FLT> UpperInterval<FLT> sub_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return sub(x1,x2);
}

template<class FLT> UpperInterval<FLT> sub_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return sub(x1,n2);
}

template<class FLT> Approximation<FLT> sub_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return sub(x1,x2);
}

template<class FLT> FLT mul_no_err(FLT const& x1, FLT const& x2) {
    FLT::set_rounding_to_nearest();
    FLT r(x1.raw() * x2.raw());
    FLT::set_rounding_upward();
    return r;
}

template<class FLT> Approximation<FLT> mul_no_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2) {
    return x1*x2;
}

template<class FLT> FLT mul_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    Rounded<FLT> mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(cast_rounded(x1) * cast_rounded(x2));
    FLT::set_rounding_upward();
    Rounded<FLT> u=cast_rounded(x1) * cast_rounded(x2);
    Rounded<FLT> ml=mx1*cast_rounded(x2);
    cast_rounded(e) += (u+ml)/2;
    return r;
}

template<class FLT> FLT mul_err(FLT const& x, ValidatedApproximation<FLT> const& c, Error<FLT>& e) {
    Rounded<FLT> const& xv=x.raw();
    Rounded<FLT> const& cu=c.upper_raw();
    Rounded<FLT> const& cm=c.middle_raw();
    Rounded<FLT> const& cl=c.lower_raw();
    Rounded<FLT>& re=cast_rounded(e);
    FLT::set_rounding_to_nearest();
    Rounded<FLT> rv=xv*cm;
    FLT::set_rounding_upward();
    if(cast_exact(xv)>=0) {
        Rounded<FLT> mcl=-cl;
        Rounded<FLT> u=xv*cu;
        Rounded<FLT> ml=xv*mcl;
        re+=hlf(u+ml);
    } else {
        Rounded<FLT> mcu=-cu;
        Rounded<FLT> u=xv*cl;
        Rounded<FLT> ml=xv*mcu;
        re+=hlf(u+ml);
    }
    return FLT(rv);
}

template<class FLT> FLT mul_err(FLT const& x, Bounds<FLT> const& c, Error<FLT>& e) {
    return mul_err(x,ValidatedApproximation<FLT>(c),e);
}

template<class FLT> FLT mul_err(FLT const& x1, Nat n2, Error<FLT>& e) {
    return mul_err(x1,FLT(n2,x1.precision()),e);
}

template<class FLT> Bounds<FLT> mul_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return mul(x1,x2);
}

template<class FLT> Bounds<FLT> mul_err(Bounds<FLT> const& x1, ValidatedApproximation<FLT> const& x2, Error<FLT>& e) {
    return mul(x1,static_cast<Bounds<FLT>>(x2));
}

template<class FLT> Bounds<FLT> mul_err(Bounds<FLT> const& x1, Nat const& n2, Error<FLT>& e) {
    return mul(x1,n2);
}

template<class FLT> UpperInterval<FLT> mul_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return mul(x1,x2);
}

template<class FLT> UpperInterval<FLT> mul_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return mul(x1,n2);
}

template<class FLT> Approximation<FLT> mul_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return mul(x1,x2);
}

template<class FLT> Approximation<FLT> mul_err(Approximation<FLT> const& x1, Nat n2, UnknownError<FLT>& e) {
    return mul_err(x1,Approximation<FLT>(n2,x1.precision()),e);
}

template<class FLT> FLT div_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    Rounded<FLT> mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(cast_rounded(x1) / cast_rounded(x2));
    FLT::set_rounding_upward();
    Rounded<FLT> u=cast_rounded(x1)/cast_rounded(x2);
    Rounded<FLT> ml=mx1/x2;
    cast_rounded(e) += (u+ml)/2;
    return r;
}

template<class FLT> FLT div_err(FLT const& x1, Nat n2, Error<FLT>& e) {
    return div_err(x1,FLT(n2,x1.precision()),e);
}

template<class FLT> Bounds<FLT> div_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return div(x1,x2);
}

template<class FLT> Bounds<FLT> div_err(Bounds<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return div(x1,n2);
}

template<class FLT> UpperInterval<FLT> div_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return div(x1,x2);
}

template<class FLT> UpperInterval<FLT> div_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return div(x1,n2);
}

template<class FLT> Approximation<FLT> div_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return div(x1,x2);
}

template<class FLT> Approximation<FLT> div_err(Approximation<FLT> const& x1, Nat n2, UnknownError<FLT>& e) {
    return div_err(x1,Approximation<FLT>(n2,x1.precision()),e);
}





template<class FLT> FLT fma_err(FLT const& x, FLT const& y, FLT z, Error<FLT>& e) {
    Rounded<FLT> const& xv=x.raw();
    Rounded<FLT> const& yv=y.raw();
    Rounded<FLT>const& zv=z.raw();
    Rounded<FLT>& re=cast_rounded(e);
    FLT::set_rounding_to_nearest();
    Rounded<FLT> rv=xv*yv+zv;
    FLT::set_rounding_upward();
    Rounded<FLT> myv=-yv;
    Rounded<FLT> u=xv*yv+zv;
    Rounded<FLT> ml=xv*myv-zv;
    re+=(u+ml)/2;
    return FLT(rv);
}

template<class FLT> FLT fma_err(ValidatedApproximation<FLT> const& c, FLT const& x, FLT y, Error<FLT>& e) {
    Rounded<FLT> const& xv=x.raw();
    Rounded<FLT> const& cu=c.upper_raw();
    Rounded<FLT> const& cm=c.middle_raw();
    Rounded<FLT> const& cl=c.lower_raw();
    Rounded<FLT> const& yv=y.raw();
    Rounded<FLT>& re=cast_rounded(e);
    FLT::set_rounding_to_nearest();
    Rounded<FLT> rv=xv*cm+yv;
    FLT::set_rounding_upward();
    Rounded<FLT> u,ml;
    if(xv>=0) {
        Rounded<FLT> mcl=-cl;
        u=xv*cu+yv;
        ml=xv*mcl-yv;
    } else {
        Rounded<FLT> mcu=-cu;
        u=xv*cl+yv;
        ml=xv*mcu-yv;
    }
    re+=(u+ml)/2;
    return FLT(rv);
}

template<class FLT> FLT fma_err(Bounds<FLT> const& c, FLT const& x, FLT y, Error<FLT>& e) {
    return fma_err(ValidatedApproximation<FLT>(c),x,y,e);
}

template<class FLT> Bounds<FLT> fma_err(Bounds<FLT> const& x, Bounds<FLT> const& y, Bounds<FLT> z, Error<FLT>& e) {
    return fma(x,y,z);
}

template<class FLT> UpperInterval<FLT> fma_err(UpperInterval<FLT> const& x, UpperInterval<FLT> const& y, UpperInterval<FLT> z, Error<FLT>& e) {
    return fma(x,y,z);
}

template<class FLT> Bounds<FLT> fma_err(ValidatedApproximation<FLT> const& x, Bounds<FLT> const& y, Bounds<FLT> z, Error<FLT>& e) {
    return fma(static_cast<Bounds<FLT>>(x),y,z);
}

template<class FLT> Approximation<FLT> fma_err(Approximation<FLT> const& x, Approximation<FLT> const& y, Approximation<FLT> z, UnknownError<FLT>& e) {
    return fma(x,y,z);
}



template<class V, class OP> inline decltype(auto) _sparse_apply(OP const& op_err, V& r) {
    auto e=r.error();
    e=nul(e);
    for(auto iter=r.begin(); iter!=r.end(); ++iter) {
        iter->coefficient()=op_err(iter->coefficient(),e);
    }
    r.error()+=e;
}

template<class V, class OP> inline Void _sparse_apply(OP const& op_err, V& r, const V& x, const V& y) {
    auto e=r.error();
    e=nul(e); // The maximum accumulated error
    auto xiter=x.begin();
    auto yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->index()<yiter->index()) {
            auto xv=xiter->coefficient();
            r._append(xiter->index(),op_err(xv,zero,e));
            ++xiter;
        } else if(yiter->index()<xiter->index()) {
            auto yv=yiter->coefficient();
            r._append(yiter->index(),op_err(zero,yv,e));
            ++yiter;
        } else {
            auto xv=xiter->coefficient();
            auto yv=yiter->coefficient();
            r._append(xiter->index(),op_err(xv,yv,e));
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x.end()) {
        auto xv=xiter->coefficient();
        r._append(xiter->index(),op_err(xv,zero,e));
        ++xiter;
    }
    while(yiter!=y.end()) {
        auto yv=yiter->coefficient();
        r._append(yiter->index(),op_err(zero,yv,e));
        ++yiter;
    }
    r.error()+=e;
}

template<class X> X const& add(X const& x, Zero) { return x; }
template<class X> X const& add(Zero, X const& x) { return x; }
template<class X> X const& sub(X const& x, Zero) { return x; }
template<class X> X sub(Zero, X const& x) { return neg(x); }
template<class X> Zero mul(X const& x, Zero) { return zero; }
template<class X> Zero mul(Zero, X const& x) { return zero; }


struct NegErr {
    template<class V, class E> decltype(auto) operator()(V&& v, E& e) const { return neg(std::forward<V>(v)); }
};
struct AddErr {
    template<class V, class E> inline V operator() (V const& xv, V const& yv, E& e) const { return add_err(xv,yv,e); }
    template<class V, class E> inline V const& operator() (V const& xv, Zero, E& e) const { return xv; }
    template<class V, class E> inline V const& operator() (Zero, V const& yv, E& e) const { return yv; }
};
struct SubErr {
    template<class V, class E> inline V operator() (V const& xv, V const& yv, E& e) const { return sub_err(xv,yv,e); }
    template<class V, class E> inline V const& operator() (V const& xv, Zero, E& e) const { return xv; }
    template<class V, class E> inline V operator() (Zero, V const& yv, E& e) const { return neg(yv); }
};
struct MulErr {
    template<class VX, class VY, class E> inline decltype(auto) operator() (VX const& xv, VY const& yv, E& e) const { return mul_err(xv,yv,e); }
    template<class V, class E> inline Zero operator() (V const& xv, Zero, E& e) const { return zero; }
    template<class V, class E> inline Zero operator() (Zero, V const& yv, E& e) const { return zero; }
};
struct FmaErr {
    template<class C, class VX, class VY, class E> inline decltype(auto) operator() (C const& c, VX const& xv, VY const& yv, E& e) const { return fma_err(c,xv,yv,e); }
    template<class C, class V, class E> inline V operator() (C const& c, V const& xv, Zero, E& e) const { return mul_err(c,xv,e); }
    template<class C, class V, class E> inline V const& operator() (C const&, Zero, V const& yv, E& e) const { return yv; }
};

// Inplace negation
template<class P, class FLT> Void _neg(TaylorModel<P,FLT>& r)
{
    _sparse_apply(NegErr(),r);
}

template<class P, class FLT> Void _scal(TaylorModel<P,FLT>& r, const TwoExp& c) {
    if (c.exponent()==0) { return; }
    r.error()*=Error<FLT>(c);
    _sparse_apply(std::bind(MulErr(),std::placeholders::_1,c,std::placeholders::_2),r);
}

template<class P, class FLT, class C> Void _scal(TaylorModel<P,FLT>& r, const C& c) {
    r.error()*=mag(c);
    _sparse_apply(std::bind(MulErr(),std::placeholders::_1,c,std::placeholders::_2),r);
}

/*
template<class FLT> Void _scal(TaylorModel<ValidatedTag,FLT>& r, const Bounds<FLT>& c) {
    ARIADNE_ASSERT_MSG(is_finite(c.lower_bound().raw()) && is_finite(c.upper_bound().raw()),"scal(tm,c): tm="<<r<<", c="<<c);
    ValidatedApproximation<FLT> clmu=c;
    r.error()*=mag(c);
    _sparse_apply(SMulErr{clmu},r);
}
*/

struct UnitMultiIndex { SizeType argument_size; SizeType unit_index; };


template<class P, class FLT> inline Void _incr(TaylorModel<P,FLT>& r, const MultiIndex& a) {
    for(typename TaylorModel<P,FLT>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        iter->index()+=a;
    }
}

template<class P, class FLT> inline Void _incr(TaylorModel<P,FLT>& r, SizeType j) {
    for(typename TaylorModel<P,FLT>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        ++iter->index()[j];
    }
}


template<class P, class FLT, class C> inline Void _acc(TaylorModel<P,FLT>& r, const C& c) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    using ErrorType = typename TaylorModel<P,FLT>::ErrorType;

    if (is_same_as_zero(c)) { return; }
    if(r.expansion().empty() || r.expansion().back().index().degree()>0) {
        CoefficientType rv=set_err(c,r.error());
        r.expansion().append(MultiIndex(r.argument_size()),rv);
    } else {
        UniformReference<CoefficientType> rv=(r.end()-1)->coefficient();
        ErrorType& re=r.error();
        rv=add_err(rv,c,re);
    }
    return;
}




// Compute r=x+y, assuming r is empty.
// Use a rounding mode change every iteration, as this appears to be faster
//   than using two loops
// Use opposite rounding to compute difference of upward and downward roundings,
//   as this seems to be marginally faster than changing the rounding mode
template<class P, class FLT> inline TaylorModel<P,FLT> _add(const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y)
{
    //ARIADNE_PRECONDITION(x.sweeper()==y.sweeper());
    TaylorModel<P,FLT> r(x.argument_size(),x.sweeper());
    _sparse_apply(AddErr(),r,x,y);
    r.error()+=(x.error()+y.error());
    r.sweep();
    return r;
}


template<class P, class FLT> inline TaylorModel<P,FLT> _sub(const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y)
{
    TaylorModel<P,FLT> r(x.argument_size(),x.sweeper());
    _sparse_apply(SubErr(),r,x,y);
    r.error()+=(x.error()+y.error());
    r.sweep();
    ARIADNE_DEBUG_ASSERT(r.error().raw()>=0);
    return r;
}

template<class P, class FLT> inline Void _sma(TaylorModel<P,FLT>& r, const TaylorModel<P,FLT>& x, const typename TaylorModel<P,FLT>::NumericType& c, const TaylorModel<P,FLT>& y)
{
    typedef typename FLT::PrecisionType PR;
    using namespace std::placeholders;

    if constexpr (Same<decltype(c),FloatBounds<PR>>) {
        ARIADNE_DEBUG_ASSERT_MSG(c.lower_bound().raw()<=c.upper_bound().raw(),c);
        ARIADNE_DEBUG_ASSERT_MSG(x.error().raw()>=0,"x="<<x);
        ARIADNE_DEBUG_ASSERT_MSG(y.error().raw()>=0,"y="<<y);
    }

    auto clmu=make_validated_approximation(c);
    _sparse_apply(std::bind(FmaErr(), clmu,_1,_2,_3));
//    _sparse_apply(SFmaErr<decltype(clmu)>{clmu},r,x,y);
    r.error()+=(x.error()*mag(c)+y.error());
    ARIADNE_DEBUG_ASSERT_MSG(r.error().raw()>=0,r);
}



// Compute r+=x*y
// Compute monomial-by-monomial in y
// Avoid changing rounding mode
template<class P, class FLT> inline Void _ifma(TaylorModel<P,FLT>& r, const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y)
{
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    using ErrorType = typename TaylorModel<P,FLT>::ErrorType;

    const SizeType as=r.argument_size();
    TaylorModel<P,FLT> t(as,r.sweeper());
    MultiIndex ta(as);
    CoefficientType tv(t.precision());
    ErrorType te=t.error();

    for(auto xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformConstReference<CoefficientType> xv=xiter->coefficient();

        auto riter = r.begin(); auto yiter=y.begin();
        while (riter!=r.end() && yiter!=y.end()) {
            auto ra=riter->index();
            auto rv=riter->coefficient();
            auto ya=yiter->index();
            auto yv=yiter->coefficient();
            ta = xa + ya;
            if (ra == ta) {
                tv=fma_err(xv,yv,rv,te);
                t._append(ta,tv);
                ++riter; ++yiter;
            } else if (ra < ta) {
                t._append(ra,rv);
                ++riter;
            } else { // ta<ra
                tv=mul_err(xv,yv,te);
                t._append(ta,tv);
                ++yiter;
            }
        }
        while (riter!=r.end()) {
            auto ra=riter->index();
            auto rv=riter->coefficient();
            t._append(ra,rv);
            ++riter;
        }
        while (yiter!=y.end()) {
            auto ya=yiter->index();
            auto yv=yiter->coefficient();
            ta = xa + ya;
            tv=mul_err(xv,yv,te);
            t._append(ta,tv);
            ++yiter;
        }

        t.error()=r.error()+te;
        te = 0u;

        t.sweep();
        r.expansion().swap(t.expansion());
        r.error()=t.error();
        t.clear();
    }

    ErrorType xs=nul(r.error());
    for(auto xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=mag(xiter->coefficient());
    }

    ErrorType ys=nul(r.error());
    for(auto yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=mag(yiter->coefficient());
    }

    ErrorType& re=r.error();
    const ErrorType& xe=x.error();
    const ErrorType& ye=y.error();
    re+=xe*ye;
    re+=xs*ye+ys*xe;

}

template<class P, class FLT> inline TaylorModel<P,FLT> _fma(const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y, TaylorModel<P,FLT> z) {
    ARIADNE_PRECONDITION(x.argument_size()==y.argument_size());
    //ARIADNE_PRECONDITION(x.sweeper()==y.sweeper());
    _ifma(z,x,y);
    return z;
}

template<class P, class FLT> inline TaylorModel<P,FLT> _mul(const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y) {
    ARIADNE_PRECONDITION(x.argument_size()==y.argument_size());
    //ARIADNE_PRECONDITION(x.sweeper()==y.sweeper());
    TaylorModel<P,FLT> r(x.argument_size(),x.sweeper());
    _ifma(r,x,y);
    return r;
}



} // namespace


///////////////////////////////////////////////////////////////////////////////

/*
template<class P, class FLT> struct AlgebraOperations<TaylorModel<P,FLT>>
    : NormedAlgebraOperations<TaylorModel<P,FLT>>
{
    typedef TaylorModel<P,FLT> ModelType;
    typedef typename TaylorModel<P,FLT>::NumericType NumericType;
    typedef typename ModelType::RangeType RangeType;
  public:
    using NormedAlgebraOperations<TaylorModel<P,FLT>>::apply;
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

template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Nul, ModelType const& x) -> ModelType {
    return ModelType(x.argument_size(),x.sweeper()); }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Pos, ModelType x) -> ModelType {
    return x; }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Neg, ModelType x) -> ModelType {
    _neg(x); return x; }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Add, ModelType const& x, ModelType const& y) -> ModelType {
    return _add(x,y); }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Sub, ModelType const& x, ModelType const& y) -> ModelType {
    return _sub(x,y); }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Mul, ModelType const& x, ModelType const& y) -> ModelType {
    return _mul(x,y); }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Add, ModelType x, NumericType const& c) -> ModelType {
    _acc(x,c); return x; }
template<class P, class FLT> auto AlgebraOperations<TaylorModel<P,FLT>>::apply(Mul, ModelType x, NumericType const& c) -> ModelType {
    _scal(x,c); return x; }

// TODO: Should be able to automatically generate these operations
/*
template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Max, ModelType const& x, ModelType const& y) -> ModelType;
template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Min, ModelType const& x, ModelType const& y) -> ModelType;

template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Max, ModelType const& x, NumericType const& c) -> ModelType;
template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Min, ModelType const& x, NumericType const& c) -> ModelType;
template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Max, NumericType const& c, ModelType const& x) -> ModelType;
template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Min, NumericType const& c, ModelType const& x) -> ModelType;
template<class P, class FLT> auto AlgebraOperations<P,FLT>::apply(Abs, ModelType const& x) -> ModelType;
*/

///////////////////////////////////////////////////////////////////////////////

// Truncation and error control


template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::sort() {
    this->_expansion.sort();
    return *this;
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::unique()
{
    typename TaylorModel<P,FLT>::ConstIterator advanced =this->begin();
    typename TaylorModel<P,FLT>::ConstIterator end =this->end();
    typename TaylorModel<P,FLT>::Iterator current=this->begin();
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

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::sweep() {
//    this->_sweeper.sweep(reinterpret_cast<Expansion<MultiIndex,FLT>&>(this->_expansion),reinterpret_cast<FLT&>(this->_error));
    this->_sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::sweep(const SweeperType& sweeper) {
    sweeper.sweep(this->_expansion,this->_error);
    return *this;
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::simplify() {
    return this->sweep();
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::simplify(const PropertiesType& properties) {
    return this->sweep(properties);
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::cleanup() {
    this->sort();
    this->unique();
    this->sweep();
    return *this;
}

template<class P, class FLT> TaylorModel<P,FLT>& TaylorModel<P,FLT>::clobber() {
    this->_error=0u;
    return *this;
}



///////////////////////////////////////////////////////////////////////////////

// Accuracy control

template<class P, class FLT> auto TaylorModel<P,FLT>::tolerance() const -> RawFloatType {
    using F = typename ModelNumericTraits<P,FLT>::RawFloatType;
    const ThresholdSweeper<F>* ptr=dynamic_cast<const ThresholdSweeper<F>*>(&static_cast<const SweeperInterface<F>&>(this->_sweeper));
    return (ptr) ? ptr->sweep_threshold() : F(cast_exact(std::numeric_limits<double>::epsilon()),this->precision());
}



//////////////////////////////////////////////////////////////////////////////

// Basic function operators (domain, range, evaluate)

template<class P, class FLT> UnitBox TaylorModel<P,FLT>::domain() const
{
    return UnitBox(this->argument_size(),UnitInterval());
}

template<class P, class FLT> auto TaylorModel<P,FLT>::codomain() const -> CodomainType
{
    RangeType rng=this->range();
    return cast_exact_interval(convert_interval(rng,dp));
}


// Compute the range by grouping all quadratic terms x[i]^2 with linear terms x[i]
// The range of ax^2+bx+c is a([-1,1]+b/2a)^2+(c-b^2/4a)
template<class P, class FLT> auto TaylorModel<P,FLT>::range() const -> RangeType {
    const TaylorModel<P,FLT>& tm=*this;
    const SizeType as=tm.argument_size();
    const PrecisionType prec = tm.precision();
    const CoefficientType zero(prec);
    CoefficientType constant_term(zero);
    Array<CoefficientType> linear_terms(as,zero);
    Array<CoefficientType> quadratic_terms(as,zero);
    NormType err(prec);
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
    RangeType r(-err,+err);
    if constexpr(Same<P,ValidatedTag>) {
        r=r+constant_term;
        const RangeType unit_ivl(-1,+1,this->precision());
        // If the ratio b/a is very large, then roundoff error can cause a significant
        // additional error. We compute both |a|+|b| and a([-1,+1]+b/2a)-b^2/4a and take best bound
        for(SizeType j=0; j!=as; ++j) {
            const CoefficientType& a=quadratic_terms[j];
            const CoefficientType& b=linear_terms[j];
            RangeType ql=abs(a)*unit_ivl + abs(b)*unit_ivl;
            if(not is_same_as_zero(a)) { // Explicitly test for zero
                RangeType qf=a*(sqr(unit_ivl+b/a/2))-sqr(b)/a/4;
                r += refinement(ql,qf); // NOTE: ql must be the first term in case of NaN in qf
            } else {
                r += ql;
            }
        }
    } else {
        ARIADNE_ASSERT_MSG(false,"Range only available for a Validated TaylorModel.");
    }
    return r;
}




//////////////////////////////////////////////////////////////////////////////

// ExactTag functions (max, min, abs, neg) and arithmetical functions (sqr, pow)

template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Max, const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y) {
    typedef typename TaylorModel<P,FLT>::RangeType RangeType;
    RangeType xr=x.range();
    RangeType yr=y.range();
    if(definitely(xr.lower_bound()>=yr.upper_bound())) {
        return x;
    } else if(definitely(yr.lower_bound()>=xr.upper_bound())) {
        return y;
    } else {
        return hlf((x+y)+abs(x-y));;
    }
}


template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Min, const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y) {
    typedef typename TaylorModel<P,FLT>::RangeType RangeType;
    RangeType xr=x.range();
    RangeType yr=y.range();
    if(definitely(xr.upper_bound()<=yr.lower_bound())) {
        return x;
    } else if(definitely(yr.upper_bound()<=xr.lower_bound())) {
        return y;
    } else {
        return hlf((x+y)-abs(x-y));
    }
}

template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Abs, const TaylorModel<P,FLT>& x) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    typedef typename TaylorModel<P,FLT>::RangeType RangeType;
    RangeType xr=x.range();
    if(definitely(xr.lower_bound()>=0)) {
        return x;
    } else if(definitely(xr.upper_bound()<=0)) {
        return -x;
    } else {
        // Use power series expansion $abs(x)=\sum_{i=0}^{7} p_i x^{2i} \pm e$ for $x\in[-1,+1]$ with
        // p=[0.0112167620474, 5.6963263292747541, -31.744583789655049, 100.43002481377681, -162.01366698662306, 127.45243493284417, -38.829743345344667] and e=0.035
        // TODO: Find more accurate and stable formula
        static const Nat n=7u;
        static const ExactDouble p[n]={0.0112167620474_pr, 5.6963263292747541_pr, -31.744583789655049_pr, 100.43002481377681_pr, -162.01366698662306_pr, 127.45243493284417_pr, -38.829743345344667_pr};
        static const ExactDouble err=0.035_pr;
        TaylorModel<P,FLT> r(x.argument_size(),x.sweeper());
        CoefficientType xmag=static_cast<CoefficientType>(cast_exact(mag(xr)));
        TaylorModel<P,FLT> s=x/xmag;
        s=sqr(s);
        r=p[n-1u];
        for(Nat i=0; i!=(n-1u); ++i) {
            Nat j=(n-2)-i;
            r=s*r+p[j];
        }
        r.error()+=err;
        return r*xmag;
    }
}

template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Max op, const TaylorModel<P,FLT>& x, const NumericType& c) {
    return apply(op, x, x.create_constant(c));
}
template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Min op, const TaylorModel<P,FLT>& x, const NumericType& c) {
    return apply(op, x, x.create_constant(c));
}
template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Max op, const NumericType& c, const TaylorModel<P,FLT>& x) {
    return apply(op, x.create_constant(c), x);
}
template<class P, class FLT> TaylorModel<P,FLT> AlgebraOperations<TaylorModel<P,FLT>>::apply(Min op, const NumericType& c, const TaylorModel<P,FLT>& x) {
    return apply(op, x.create_constant(c), x);
}

//////////////////////////////////////////////////////////////////////////////

// Arithmetical functions (sqr, pow)

template<class P, class FLT> TaylorModel<P,FLT> sqr(const TaylorModel<P,FLT>& x) {
    return x*x;
}

template<class P, class FLT> TaylorModel<P,FLT> pow(const TaylorModel<P,FLT>& x, Int n) {
    TaylorModel<P,FLT> r=x.create_constant(1);
    TaylorModel<P,FLT> p(x);
    while(n) { if(n%2) { r=r*p; } p=sqr(p); n/=2; }
    return r;
}



//////////////////////////////////////////////////////////////////////////////

// Composition with power series

template<class X> class Series;
template<class X> class TaylorSeries;


template<class P, class FLT> TaylorModel<P,FLT>
compose(const TaylorSeries<typename FLT::PrecisionType>& ts, const TaylorModel<P,FLT>& tv)
{
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    Sweeper<FLT> threshold_sweeper(new ThresholdSweeper<FLT>(tv.precision(),MACHINE_EPSILON));
    CoefficientType& vref=const_cast<CoefficientType&>(tv.value());
    CoefficientType vtmp=vref;
    vref=0;
    TaylorModel<P,FLT> r(tv.argument_size(),tv.sweeper());
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
template<class P, class FLT> TaylorModel<P,FLT>
compose(const AnalyticFunction& fn, const TaylorModel<P,FLT>& tm) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    using ErrorType = typename TaylorModel<P,FLT>::ErrorType;

    static const DegreeType MAX_DEGREE=20;
    static const ExactDouble MAX_TRUNCATION_ERROR=MACHINE_EPSILON;
    SweeperInterface<FLT> const& sweeper=tm.sweeper();

    FLT max_truncation_error(MAX_TRUNCATION_ERROR,sweeper.precision());
    ThresholdSweeper<FLT> const* threshold_sweeper_ptr = dynamic_cast<ThresholdSweeper<FLT> const*>(&sweeper);
    if(threshold_sweeper_ptr) { max_truncation_error=threshold_sweeper_ptr->sweep_threshold(); }

    DegreeType max_degree=MAX_DEGREE;
    GradedSweeper<FLT> const* graded_sweeper_ptr = dynamic_cast<GradedSweeper<FLT> const*>(&sweeper);
    if(graded_sweeper_ptr) { max_degree=graded_sweeper_ptr->degree(); }

    //std::cerr<<"max_truncation_error="<<max_truncation_error<<"\nmax_degree="<<(uint)max_degree<<"\n";

    DegreeType d=max_degree;
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
                 <<" is greater than maximum allowable truncation error "<<max_truncation_error);
    }

    TaylorModel<P,FLT> x=tm-c;
    TaylorModel<P,FLT> res(tm.argument_size(),tm.sweeper());
    res+=centre_series[d];
    for(DegreeType i=0; i!=d; ++i) {
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

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::_embed_error(const TaylorModel<P,FLT>& tm) {
    if constexpr (Same<ErrorType,UnknownError<FLT>>) { return tm; }// do nothing
    else {
        const SizeType as=tm.argument_size();
        TaylorModel<P,FLT> rtm(as+1u,tm.sweeper());
        MultiIndex ra(as+1u);

        // The new error term is first in reverse lexicographic order.
        CoefficientType err_coef=static_cast<CoefficientType>(cast_exact(tm.error()));
        ra[as]=1;
        rtm._append(ra,err_coef);
        ra[as]=0;

        // Copy new terms
        for(typename TaylorModel<P,FLT>::ConstIterator iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
            UniformConstReference<MultiIndex> xa=iter->index();
            UniformConstReference<CoefficientType> xv=iter->coefficient();
            for(SizeType j=0; j!=as; ++j) { ra[j]=xa[j]; }
            rtm._append(ra,xv);
        }
        return rtm;
    }
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::_discard_variables(const TaylorModel<P,FLT>& tm, Array<SizeType> const& discarded_variables) {
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
    TaylorModel<P,FLT> rtm(number_of_kept_variables,tm.sweeper());
    rtm.expansion().reserve(tm.number_of_nonzeros()+1u);

    // Set the uniform error of the original model
    // If index_of_error == number_of_error_variables, then the error is kept as a uniform error bound
    MultiIndex ra(number_of_kept_variables);
    ErrorType derr=mag(tm.error()); // Magnitude of discarded terms
    for(typename TaylorModel<P,FLT>::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
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


template<class P, class FLT> Void TaylorModel<P,FLT>::antidifferentiate(SizeType k) {
    TaylorModel<P,FLT>& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());

    ErrorType e=nul(this->error());
    for(typename TaylorModel<P,FLT>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        UniformConstReference<MultiIndex> xa=xiter->index();
        UniformReference<CoefficientType> xv=xiter->coefficient();
        xa[k]+=1;
        Nat c=xa[k];
        xv=div_err(xv,c,e);
    }
    x.error()+=e;
}

template<class P, class FLT> TaylorModel<P,FLT> antiderivative(const TaylorModel<P,FLT>& x, SizeType k) {
    TaylorModel<P,FLT> r(x);
    r.antidifferentiate(k);
    return r;
}


// Compute derivative inplace by computing term-by-term, switching the rounding mode
// Note that since some terms may be eliminated, requiring two iterators.
template<class P, class FLT> Void TaylorModel<P,FLT>::differentiate(SizeType k) {
    TaylorModel<P,FLT> const& x=*this;
    ARIADNE_PRECONDITION(k<x.argument_size());
    // ARIADNE_PRECONDITION_MSG(x.error().raw()==0,x);
    this->clobber();

    TaylorModel<P,FLT>& r=*this;
    ErrorType& re=r.error();
    typename TaylorModel<P,FLT>::Iterator riter=r.begin();
    for(typename TaylorModel<P,FLT>::ConstIterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
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




template<class P, class FLT> TaylorModel<P,FLT> derivative(const TaylorModel<P,FLT>& x, SizeType k) {
    TaylorModel<P,FLT> rx=x; rx.differentiate(k); return rx;

    ARIADNE_ASSERT(k<x.argument_size());
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;

    MultiIndex ra(x.argument_size()); CoefficientType rv; Nat c;

    TaylorModel<P,FLT> r(x.argument_size(),x.sweeper());
    typename TaylorModel<P,FLT>::ErrorType& re=r.error();
    for(typename TaylorModel<P,FLT>::Iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
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

template<class P, class FLT> auto TaylorModel<P,FLT>::_evaluate(const TaylorModel<P,FLT>& tm, const Vector<IntervalNumericType>& x) -> ArithmeticType<CoefficientType,IntervalNumericType> {
    return horner_evaluate(tm.expansion(),x)+pm(tm.error());
}

template<class P, class FLT> auto TaylorModel<P,FLT>::_evaluate(const TaylorModel<P,FLT>& tm, const Vector<ValidatedNumericType>& x) -> ArithmeticType<CoefficientType,ValidatedNumericType> {
    if constexpr (Same<NumericType,ValidatedNumericType>) {
        return horner_evaluate(tm.expansion(),x)+pm(tm.error());
    } else {
        return horner_evaluate(tm.expansion(),Vector<NumericType>(x))+pm(tm.error());
    }
}

template<class P, class FLT> auto TaylorModel<P,FLT>::_evaluate(const TaylorModel<P,FLT>& tm, const Vector<ApproximateNumericType>& x) -> ArithmeticType<CoefficientType,ApproximateNumericType> {
    return horner_evaluate(tm.expansion(),x);
}


template<class P, class FLT> auto TaylorModel<P,FLT>::_gradient(const TaylorModel<P,FLT>& tm, const Vector<NumericType>& x) -> Covector<NumericType> {
    Vector<Differential<NumericType>> dx=Differential<NumericType>::variables(1u,x);
    Differential<NumericType> df=horner_evaluate(tm.expansion(),dx)+static_cast<NumericType>(pm(tm.error()));
    return gradient(df);
}



template<class P, class FLT> TaylorModel<P,FLT>
TaylorModel<P,FLT>::_compose(TaylorModel<P,FLT> const& x, Vector<TaylorModel<P,FLT>> const& y) {
    return horner_evaluate(x.expansion(),y)+pm(x.error());
}

template<class P, class FLT> Void TaylorModel<P,FLT>::unscale(IntervalDomainType const& codom) {
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
    TaylorModel<P,FLT>& tm=*this;
    auto ivl=convert_interval(codom,this->precision());
    ARIADNE_ASSERT_MSG(decide(ivl.lower_bound()<=ivl.upper_bound()),"Cannot unscale TaylorModel<P,FLT> "<<tm<<" from empty interval "<<ivl);

    if(codom.lower_bound()==codom.upper_bound()) {
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

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::_compose(const Unscaling& u, const TaylorModel<P,FLT>& y) {
    TaylorModel<P,FLT> r=y; r.unscale(u.domain()); return r;
}

template<class P, class FLT> TaylorModel<P,FLT>
TaylorModel<P,FLT>::_compose(const TaylorModel<P,FLT>& x, const VectorUnscaling& u, const Vector<TaylorModel<P,FLT>>& y) {
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

template<class FLT> class Powers<Bounds<FLT>> {
    typedef Bounds<FLT> T;
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



template<class P, class FLT> TaylorModel<P,FLT>
TaylorModel<P,FLT>::_partial_evaluate(const TaylorModel<P,FLT>& x, SizeType k, NumericType c)
{
    const SizeType as=x.argument_size();
    Vector<TaylorModel<P,FLT>> y(as,TaylorModel<P,FLT>(as-1u,x.sweeper()));
    for(SizeType i=0; i!=k; ++i) { y[i]=TaylorModel<P,FLT>::coordinate(as-1u,i,x.sweeper()); }
    y[k]=TaylorModel<P,FLT>::constant(as-1u,c,x.sweeper());
    for(SizeType i=k+1; i!=as; ++i) { y[i]=TaylorModel<P,FLT>::coordinate(as-1u,i-1u,x.sweeper()); }
    return compose(x,y);
}



template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::_embed(SizeType as1, const TaylorModel<P,FLT>& tm2, SizeType as3) {
    return TaylorModel<P,FLT>(embed(as1,tm2.expansion(),as3),tm2.error(),tm2.sweeper());
}


template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::_split(const TaylorModel<P,FLT>& tm, SizeType k, SplitPart h) {
    const DegreeType deg=tm.degree();
    const SizeType as=tm.argument_size();
    SweeperType swp=tm.sweeper();

    TaylorModel<P,FLT> r(tm);

    // Divide all coefficients by 2^a[k]
    // This can be done exactly
    for(typename TaylorModel<P,FLT>::Iterator iter=r.begin(); iter!=r.end(); ++iter) {
        const DegreeType ak=iter->index()[k];
        UniformReference<CoefficientType> c=iter->coefficient();
        c/=pow(two,ak);
    }

    if(h==SplitPart::MIDDLE) { return r; }
    Int tr=( h==SplitPart::UPPER ? +1 : -1 );

    // Replace x[k] with x[k]+tr

    // Split variables by degree in x[k]
    Array<TaylorModel<P,FLT>> ary(deg+1u,TaylorModel<P,FLT>(as,swp));
    for(typename TaylorModel<P,FLT>::ConstIterator iter=r.begin(); iter!=r.end(); ++iter) {
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
            for(typename TaylorModel<P,FLT>::Iterator iter=ary[j].begin(); iter!=ary[j].end(); ++iter) {
                ++iter->index()[k];
            }
         }
    }

    return r;
}



///////////////////////////////////////////////////////////////////////////////

// Banach algebra operations

template<class FLT> TaylorModel<ValidatedTag,UpperInterval<FLT>> operator/(TaylorModel<ValidatedTag,UpperInterval<FLT>> const& x1, FLT const& x2) {
    return x1/UpperInterval<FLT>(x2); }
template<class FLT> TaylorModel<ValidatedTag,UpperInterval<FLT>>& operator/=(TaylorModel<ValidatedTag,UpperInterval<FLT>>& x1, FLT const& x2) {
    return x1/=UpperInterval<FLT>(x2); }
template<class FLT> TaylorModel<ValidatedTag,UpperInterval<FLT>> operator-(TaylorModel<ValidatedTag,UpperInterval<FLT>> x1, FLT const& x2) {
    return x1-UpperInterval<FLT>(x2); }


template<class P, class FLT> auto TaylorModel<P,FLT>::average() const -> ValueType {
    return cast_exact(cast_singleton((*this)[MultiIndex::zero(this->argument_size())]));
}

template<class P, class FLT> auto TaylorModel<P,FLT>::radius() const -> NormType {
    typename TaylorModel<P,FLT>::NormType r(0u,this->precision());
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->index().degree()!=0) {
            r+=mag(iter->coefficient());
        }
    }
    r+=this->error();
    return r;
}

template<class P, class FLT> auto TaylorModel<P,FLT>::norm() const -> NormType {
    NormType r(0u,this->precision());
    for(ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        r+=mag(iter->coefficient());
    }
    r+=this->error();
    return r;
}

template<class P, class FLT> typename TaylorModel<P,FLT>::NormType norm(const TaylorModel<P,FLT>& tm) {
    return tm.norm();
}


template<class P, class FLT> Bool TaylorModel<P,FLT>::_refines(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2)
{
    ARIADNE_ASSERT(tm1.argument_size()==tm2.argument_size());
    TaylorModel<P,FLT> d=tm2;
    d.error()=0u;
    d-=tm1;
    return d.norm().raw() <= tm2.error().raw();
}


template<class P, class FLT> Bool TaylorModel<P,FLT>::_consistent(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    if constexpr(Same<P,ValidatedTag>) {
        return (Ariadne::norm(tm1-tm2).raw() <= (2u*(tm1.error()+tm2.error())).raw());
    } else {
        return true;
    }
}

template<class P, class FLT> Bool TaylorModel<P,FLT>::_inconsistent(const TaylorModel<P,FLT>& tm1, const TaylorModel<P,FLT>& tm2)
{
    ARIADNE_PRECONDITION(tm1.argument_size()==tm2.argument_size());
    return (mag(tm1.value()-tm2.value()).raw() > (tm1.error()+tm2.error()).raw());
}

template<class FLT> inline Bool inconsistent(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2) {
    return definitely(disjoint(x1,x2));
}

template<class P, class FLT> TaylorModel<P,FLT> TaylorModel<P,FLT>::_refinement(const TaylorModel<P,FLT>& x, const TaylorModel<P,FLT>& y) {
    if constexpr (Same<P,ValidatedTag>) {
        TaylorModel<P,FLT> r(x.argument_size(),x.sweeper());

        PrecisionType pr=r.precision();
        ErrorType max_error=nul(r.error());

        const ErrorType& xe=x.error();
        const ErrorType& ye=y.error();

        CoefficientType rv(pr),xv(pr),yv(pr);
        FLT xu(pr),yu(pr),mxl(pr),myl(pr),u(pr),ml(pr);
        MultiIndex a;

        typename TaylorModel<P,FLT>::ConstIterator xiter=x.begin();
        typename TaylorModel<P,FLT>::ConstIterator yiter=y.begin();
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

            if constexpr (AnInterval<CoefficientType>) {
                auto xve=xv+pm(xe); auto yve=yv+pm(ye);
                if (definitely(disjoint(xve,yve))) {
                    ARIADNE_THROW(IntersectionException,"refinement(TaylorModel<ValidatedTag,FLT>,TaylorModel<ValidatedTag,FLT>)",x<<" and "<<y<<" are inconsistent.");
                }
                auto rve=refinement( xve, yve );
                // FIXME: This adds the error to all coefficients,
                // resulting in an overly-large refinement
                r.expansion().append(a,rve);
            } else {
                auto rve=refinement( xv.pm(xe), yv.pm(ye) );
                if(rve.error().raw()<0.0_x) {
                    ARIADNE_THROW(IntersectionException,"refinement(TaylorModel<ValidatedTag,FLT>,TaylorModel<ValidatedTag,FLT>)",x<<" and "<<y<<" are inconsistent.");
                }
                if(rve.value()!=0.0_x) { r.expansion().append(a,rve.value()); }
                max_error=max(max_error,rve.error());
            }
        }
        r.error()=max_error;
        return r;
    } else {
        // Refinement for Approximate models doesn't really make sense
        ARIADNE_NOT_IMPLEMENTED;
    }
}


///////////////////////////////////////////////////////////////////////////////

// Input/output operators

template<class P, class FLT> OutputStream& TaylorModel<P,FLT>::str(OutputStream& os) const {
    TaylorModel<P,FLT> const& tm=*this;

    // Set the variable names to be 'parameter' s0,s1,..
    Array<StringType> variable_names(tm.argument_size());
    for(SizeType j=0; j!=tm.argument_size(); ++j) {
        StringStream sstr;
        sstr << 's' << j;
        variable_names[j]=sstr.str();
    }

    //os << "TaylorModel<P,FLT>";
    os << "TM["<<tm.argument_size()<<"](";
    Expansion<MultiIndex,CoefficientType> e=tm.expansion();
    e.sort(GradedIndexLess());
    e._write(os,variable_names);
    return os << "+/-" << tm.error() << ")";
}

template<class P, class FLT> OutputStream& TaylorModel<P,FLT>::repr(OutputStream& os) const {
    return this->str(os);
}

///////////////////////////////////////////////////////////////////////////////

// Vector-valued named constructors


template<class P, class FLT> Vector<TaylorModel<P,FLT>> TaylorModel<P,FLT>::zeros(SizeType rs, SizeType as, SweeperType swp)
{
    Vector<TaylorModel<P,FLT>> result(rs,TaylorModel<P,FLT>::zero(as,swp));
    return result;
}


template<class P, class FLT> Vector<TaylorModel<P,FLT>> TaylorModel<P,FLT>::constants(SizeType as, const Vector<NumericType>& c, SweeperType swp)
{
    Vector<TaylorModel<P,FLT>> result(c.size(),TaylorModel<P,FLT>::zero(as,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        result[i]=TaylorModel<P,FLT>::constant(as,c[i],swp);
    }
    return result;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> TaylorModel<P,FLT>::coordinates(SizeType as, SweeperType swp)
{
    Vector<TaylorModel<P,FLT>> result(as,TaylorModel<P,FLT>::zero(as,swp));
    for(SizeType i=0; i!=as; ++i) { result[i]=TaylorModel<P,FLT>::coordinate(as,i,swp); }
    return result;
}

template<class P, class FLT> Vector<TaylorModel<P,FLT>> TaylorModel<P,FLT>::scalings(const BoxDomainType& d, SweeperType swp)
{
    Vector<TaylorModel<P,FLT>> result(d.size(),TaylorModel<P,FLT>::zero(d.size(),swp));
    for(SizeType i=0; i!=d.size(); ++i) {
        result[i]=TaylorModel<P,FLT>::scaling(d.size(),i,d[i],swp);
    }
    return result;
}



///////////////////////////////////////////////////////////////////////////////

// Jacobian matrices

// Compute the Jacobian over an arbitrary domain
template<class FLT> Matrix<Bounds<FLT>>
jacobian(const Vector<TaylorModel<ValidatedTag,FLT>>& f, const Vector<Bounds<FLT>>& x) {
    Vector< Differential<Bounds<FLT>> > dx=Differential<Bounds<FLT>>::variables(1u,x);
    Vector< Differential<Bounds<FLT>> > df(f.size(),x.size(),1u);
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<Bounds<FLT>> J=jacobian(df);
    return J;
}

// Compute the Jacobian over an arbitrary domain
template<class FLT> Matrix<Bounds<FLT>>
jacobian(const Vector<TaylorModel<ValidatedTag,FLT>>& f, const Vector<Bounds<FLT>>& x, const Array<SizeType>& p) {
    Vector<Differential<Bounds<FLT>>> dx(x.size(),x.size(),1u);
    for(SizeType j=0; j!=x.size(); ++j) {
        dx[j]=Differential<Bounds<FLT>>::constant(p.size(),1u,x[j]); }
    for(SizeType k=0; k!=p.size(); ++k) {
        SizeType j=p[k];
        dx[j]=Differential<Bounds<FLT>>::variable(p.size(),1u,x[j],k); }
    Vector< Differential<Bounds<FLT>> > df(f.size());
    for(SizeType i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<Bounds<FLT>> J=jacobian(df);
    return J;
}

// Compute the Jacobian at the origin
template<class P, class FLT> Matrix<FLT>
jacobian_value(const Vector<TaylorModel<P,FLT>>& f) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
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
template<class P, class FLT> Matrix<FLT>
jacobian_value(const Vector<TaylorModel<P,FLT>>& f, const Array<SizeType>& p) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
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
template<class P, class FLT> Matrix<typename TaylorModel<P,FLT>::RangeType>
jacobian_range(const Vector<TaylorModel<P,FLT>>& f) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    using RangeType = typename TaylorModel<P,FLT>::RangeType;
    SizeType rs=f.size();
    SizeType as=f.zero_element().argument_size();
    Matrix<RangeType> J(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<P,FLT>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
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
template<class P, class FLT> Matrix<typename TaylorModel<P,FLT>::RangeType>
jacobian_range(const Vector<TaylorModel<P,FLT>>& f, const Array<SizeType>& p) {
    using CoefficientType = typename TaylorModel<P,FLT>::CoefficientType;
    using RangeType = typename TaylorModel<P,FLT>::RangeType;
    SizeType rs=f.size();
    SizeType ps=p.size();
    Matrix<RangeType> J(rs,ps);
    for(SizeType i=0; i!=rs; ++i) {
        for(typename TaylorModel<P,FLT>::ConstIterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
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


template<class FLT> TaylorModel<ValidatedTag,FLT> value_coefficients(TaylorModel<ValidatedTag,Bounds<FLT>> const& tm) {
    TaylorModel<ValidatedTag,FLT> r(tm.argument_size(),tm.sweeper());
    auto& e=r.error();

    for (auto term : tm.expansion()) {
        FLT c=set_err(term.coefficient(),e);
        r.expansion().append(term.index(), c);
    }
    e+=tm.error();
    return r;
}

template<class FLT> TaylorModel<ValidatedTag,Bounds<FLT>> exact_coefficients(TaylorModel<ValidatedTag,Bounds<FLT>> const& tm) {
    TaylorModel<ValidatedTag,Bounds<FLT>> r(tm.argument_size(),tm.sweeper());
    auto& e=r.error();

    for (auto term : tm.expansion()) {
        FLT c=set_err(term.coefficient(),e);
        r.expansion().append(term.index(), c);
    }
    e+=tm.error();
    return r;
}




} //namespace Ariadne


