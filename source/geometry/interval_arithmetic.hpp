/***************************************************************************
 *            geometry/interval_arithmetic.hpp
 *
 *  Copyright  2013-22  Pieter Collins
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

/*! \file geometry/interval_arithmetic.hpp
 *  \brief
 */

#ifndef ARIADNE_INTERVAL_ARITHMETIC_HPP
#define ARIADNE_INTERVAL_ARITHMETIC_HPP

#include "utility/module.hpp"

#include "interval.decl.hpp"

namespace Ariadne {


struct DefineIntervalArithmeticComparisonOperations {
//    template<class U1, class U2> friend decltype(auto) operator==(Interval<U1> const& ivl1, Interval<U2> const& ivl2) {
//        return !(ivl1 != ivl2); }
//    template<class U1, class U2> friend decltype(auto) operator!=(Interval<U1> const& ivl1, Interval<U2> const& ivl2) {
//        return cast_positive_logical(ivl1.lower_bound()> ivl2.upper_bound() || ivl1.upper_bound()< ivl2.lower_bound()); }
    template<class U1, class U2> friend Boolean operator<=(Interval<U1> const& ivl1, Interval<U2> const& ivl2);
    template<class U1, class U2> friend Boolean operator>=(Interval<U1> const& ivl1, Interval<U2> const& ivl2);
    template<class U1, class U2> friend Boolean operator< (Interval<U1> const& ivl1, Interval<U2> const& ivl2);
    template<class U1, class U2> friend Boolean operator> (Interval<U1> const& ivl1, Interval<U2> const& ivl2);
};

template<class U> class DefineIntervalArithmeticOperations;


//template<class U> class DeclareIntervalArithemticOperations {
template<class U> class DefineIntervalArithmeticOperations
    : public DefineIntervalArithmeticComparisonOperations
{
    friend Interval<U> operator+(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> operator-(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> operator*(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> operator/(Interval<U> const& ivl1, Interval<U> const& ivl2);

    friend Interval<U> add(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> sub(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> mul(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> div(Interval<U> const& ivl1, Interval<U> const& ivl2);
};

Boolean cast_positive_logical(Boolean);
ValidatedSierpinskian cast_positive_logical(ValidatedLowerKleenean);
ValidatedNegatedSierpinskian cast_positive_logical(ValidatedUpperKleenean);

template<class U> class IntervalFactory;
template<class F> using UpperIntervalFactory = IntervalFactory<UpperBound<F>>;
template<class F> using ApproximateIntervalFactory = IntervalFactory<Approximation<F>>;

template<class F> class DefineIntervalArithmeticOperations<UpperBound<F>>
    : public DefineIntervalArithmeticComparisonOperations
{
    using U = UpperBound<F>;
  public:
    friend Interval<U> operator+(Interval<U> const& ivl);
    friend Interval<U> operator-(Interval<U> const& ivl);
    friend Interval<U> operator+(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> operator-(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> operator*(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> operator/(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U>& operator+=(Interval<U>& ivl1, Interval<U> const& ivl2);
    friend Interval<U>& operator-=(Interval<U>& ivl1, Interval<U> const& ivl2);
    friend Interval<U>& operator*=(Interval<U>& ivl1, Interval<U> const& ivl2);
    friend Interval<U>& operator/=(Interval<U>& ivl1, Interval<U> const& ivl2);

    friend Interval<U> operator+(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> operator-(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> operator*(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> operator/(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> operator+(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> operator-(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> operator*(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> operator/(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U>& operator+=(Interval<U>& ivl1, Bounds<F> const& x2);
    friend Interval<U>& operator-=(Interval<U>& ivl1, Bounds<F> const& x2);
    friend Interval<U>& operator*=(Interval<U>& ivl1, Bounds<F> const& x2);
    friend Interval<U>& operator/=(Interval<U>& ivl1, Bounds<F> const& x2);

    friend Interval<U> operator+(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> operator-(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> operator*(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> operator/(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> operator+(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> operator-(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> operator*(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> operator/(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U>& operator+=(Interval<U>& ivl1, ValidatedNumber y2);
    friend Interval<U>& operator-=(Interval<U>& ivl1, ValidatedNumber y2);
    friend Interval<U>& operator*=(Interval<U>& ivl1, ValidatedNumber y2);
    friend Interval<U>& operator/=(Interval<U>& ivl1, ValidatedNumber y2);

    friend Boolean eq(Interval<U> const& ivl1, Int n2) {
        return ivl1.lower_bound()==n2 && ivl1.upper_bound()==n2; }
    friend ValidatedKleenean lt(Interval<U> const& ivl1, Int n2) {
        return ivl1.upper_bound()< n2 ? true : ivl1.lower_bound()>=n2 ? false : indeterminate; }
    friend ValidatedKleenean lt(Int n1, Interval<U> const& ivl2) {
        return n1< ivl2.lower_bound() ? true : n1>=ivl2.upper_bound() ? false : indeterminate; }

    friend Boolean operator==(Interval<U> const& ivl1, Int n2);
    friend Boolean operator!=(Interval<U> const& ivl1, Int n2);
    friend ValidatedKleenean operator<=(Interval<U> const& ivl1, Int n2);
    friend ValidatedKleenean operator>=(Interval<U> const& ivl1, Int n2);
    friend ValidatedKleenean operator< (Interval<U> const& ivl1, Int n2);
    friend ValidatedKleenean operator> (Interval<U> const& ivl1, Int n2);
    friend Boolean operator==(Int n1, Interval<U> const& ivl2);
    friend Boolean operator!=(Int n1, Interval<U> const& ivl2);
    friend ValidatedKleenean operator<=(Int n1, Interval<U> const& ivl2);
    friend ValidatedKleenean operator>=(Int n1, Interval<U> const& ivl2);
    friend ValidatedKleenean operator< (Int n1, Interval<U> const& ivl2);
    friend ValidatedKleenean operator> (Int n1, Interval<U> const& ivl2);

/*
    friend Boolean operator==(Interval<U> const& ivl1, Int n2) {
        return eq(ivl1,n2); }
    friend Boolean operator!=(Interval<U> const& ivl1, Int n2) {
        return not eq(ivl1,n2); }
    friend ValidatedKleenean operator<=(Interval<U> const& ivl1, Int n2) {
        return not lt(n2,ivl1); }
    friend ValidatedKleenean operator>=(Interval<U> const& ivl1, Int n2) {
        return not lt(ivl1,n2); }
    friend ValidatedKleenean operator< (Interval<U> const& ivl1, Int n2) {
        return lt(ivl1,n2); }
    friend ValidatedKleenean operator> (Interval<U> const& ivl1, Int n2) {
        return lt(n2,ivl1); }
    friend Boolean operator==(Int n1, Interval<U> const& ivl2) {
        return ivl2==n1; }
    friend Boolean operator!=(Int n1, Interval<U> const& ivl2) {
        return ivl2!=n1; }
    friend ValidatedKleenean operator<=(Int n1, Interval<U> const& ivl2) {
        return n1<=reinterpret_cast<Bounds<F>const&>(ivl2); }
    friend ValidatedKleenean operator>=(Int n1, Interval<U> const& ivl2) {
        return n1>=reinterpret_cast<Bounds<F>const&>(ivl2); }
    friend ValidatedKleenean operator< (Int n1, Interval<U> const& ivl2) {
        return n1< reinterpret_cast<Bounds<F>const&>(ivl2); }
    friend ValidatedKleenean operator> (Int n1, Interval<U> const& ivl2) {
        return n1> reinterpret_cast<Bounds<F>const&>(ivl2); }
*/
    friend Interval<U> nul(Interval<U> const& ivl);
    friend Interval<U> pos(Interval<U> const& ivl);
    friend Interval<U> neg(Interval<U> const& ivl);
    friend Interval<U> sqr(Interval<U> const& ivl);
    friend Interval<U> hlf(Interval<U> const& ivl);
    friend Interval<U> rec(Interval<U> const& ivl);
    friend Interval<U> add(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> sub(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> mul(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> div(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> fma(Interval<U> const& ivl1, Interval<U> const& ivl2, Interval<U> const& ivl3);
    friend Interval<U> pow(Interval<U> const& ivl1, Int n);
    friend Interval<U> sqrt(Interval<U> const& ivl);
    friend Interval<U> exp(Interval<U> const& ivl);
    friend Interval<U> log(Interval<U> const& ivl);
    friend Interval<U> sin(Interval<U> const& ivl);
    friend Interval<U> cos(Interval<U> const& ivl);
    friend Interval<U> tan(Interval<U> const& ivl);
    friend Interval<U> asin(Interval<U> const& ivl);
    friend Interval<U> acos(Interval<U> const& ivl);
    friend Interval<U> atan(Interval<U> const& ivl);

    friend Interval<U> max(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> min(Interval<U> const& ivl1, Interval<U> const& ivl2);
    friend Interval<U> abs(Interval<U> const& ivl);

    friend Interval<U> add(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> sub(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> mul(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> div(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> add(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> sub(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> mul(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> div(Bounds<F> const& x1, Interval<U> const& ivl2);

    friend Interval<U> add(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> sub(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> mul(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> div(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> add(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> sub(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> mul(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> div(ValidatedNumber y1, Interval<U> const& ivl2);

    friend Interval<U> max(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> min(Interval<U> const& ivl1, Bounds<F> const& x2);
    friend Interval<U> max(Bounds<F> const& x1, Interval<U> const& ivl2);
    friend Interval<U> min(Bounds<F> const& x1, Interval<U> const& ivl2);

    friend Interval<U> max(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> min(Interval<U> const& ivl1, ValidatedNumber y2);
    friend Interval<U> max(ValidatedNumber y1, Interval<U> const& ivl2);
    friend Interval<U> min(ValidatedNumber y1, Interval<U> const& ivl2);

    using PrecisionType = typename F::PrecisionType;
    PrecisionType precision() const;

    friend ApproximateInterval<F> operator+(Approximation<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator-(Approximation<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator*(Approximation<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator/(Approximation<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator+(ApproximateInterval<F> const& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator-(ApproximateInterval<F> const& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator*(ApproximateInterval<F> const& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator/(ApproximateInterval<F> const& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator+=(ApproximateInterval<F>& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator-=(ApproximateInterval<F>& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator*=(ApproximateInterval<F>& ivl1, Approximation<F> const& ivl2);
    friend ApproximateInterval<F> operator/=(ApproximateInterval<F>& ivl1, Approximation<F> const& ivl2);

    friend Error<F> mag(UpperInterval<F> const&);

};


template<class F> struct DefineIntervalArithmeticOperations<Approximation<F>>
//    : DeclareNumericOperations<ApproximateInterval<F>>
//    , DeclareComparisonOperations<ApproximateInterval<F>,ApproximateKleenean>
//    , DefineFieldOperators<ApproximateInterval<F>>
//    , DefineComparisonOperators<ApproximateInterval<F>,ApproximateKleenean>
//    , ProvideConvertedFieldOperations<ApproximateInterval<F>,Approximation<F>>
//    , ProvideConvertedComparisonOperations<ApproximateInterval<F>,Approximation<F>,ApproximateInterval<F>,ApproximateKleenean>
//   , ProvideConcreteGenericElementaryOperations<ApproximateInterval<F>,ApproximateNumber>
//    , ProvideConcreteGenericComparisonOperations<ApproximateInterval<F>,ApproximateNumber,ApproximateKleenean>
{
    friend ApproximateInterval<F> operator+(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator-(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator*(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator/(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator+=(ApproximateInterval<F>& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator-=(ApproximateInterval<F>& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator*=(ApproximateInterval<F>& ivl1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator/=(ApproximateInterval<F>& ivl1, ApproximateInterval<F> const& ivl2);

    friend ApproximateInterval<F> operator+(Approximation<F> const& x1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator-(Approximation<F> const& x1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator*(Approximation<F> const& x1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator/(Approximation<F> const& x1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator+(ApproximateInterval<F> const& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator-(ApproximateInterval<F> const& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator*(ApproximateInterval<F> const& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator/(ApproximateInterval<F> const& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator+=(ApproximateInterval<F>& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator-=(ApproximateInterval<F>& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator*=(ApproximateInterval<F>& ivl1, Approximation<F> const& x2);
    friend ApproximateInterval<F> operator/=(ApproximateInterval<F>& ivl1, Approximation<F> const& x2);

    friend ApproximateInterval<F> operator+(ApproximateNumber y1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator-(ApproximateNumber y1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator*(ApproximateNumber y1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator/(ApproximateNumber y1, ApproximateInterval<F> const& ivl2);
    friend ApproximateInterval<F> operator+(ApproximateInterval<F> const& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator-(ApproximateInterval<F> const& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator*(ApproximateInterval<F> const& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator/(ApproximateInterval<F> const& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator+=(ApproximateInterval<F>& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator-=(ApproximateInterval<F>& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator*=(ApproximateInterval<F>& ivl1, ApproximateNumber y2);
    friend ApproximateInterval<F> operator/=(ApproximateInterval<F>& ivl1, ApproximateNumber y2);

    typedef typename F::PrecisionType PR; typedef PR PrecisionType;

    friend ApproximateIntervalFactory<F> factory(ApproximateInterval<F> const& ivl) {
        return ApproximateIntervalFactory<F>(ivl.upper_bound().precision()); }
    friend PositiveApproximation<F> mag(ApproximateInterval<F> const& ivl) {
        return cast_positive(max(-ivl.lower_bound(),ivl.upper_bound())); }
    friend ApproximateInterval<F> add(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2) {
        return ApproximateInterval<F>(ivl1.lower_bound()+ivl2.lower_bound(),ivl1.upper_bound()+ivl2.upper_bound()); }
    friend ApproximateInterval<F> sub(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2) {
        return ApproximateInterval<F>(ivl1.lower_bound()-ivl2.upper_bound(),ivl1.upper_bound()-ivl2.lower_bound()); }
    friend ApproximateInterval<F> mul(ApproximateInterval<F> const& ivl1, ApproximateInterval<F> const& ivl2) {
        auto rll=ivl1.lower_bound()*ivl2.lower_bound(); auto rlu=ivl1.lower_bound()*ivl2.upper_bound();
        auto rul=ivl1.upper_bound()*ivl2.lower_bound(); auto ruu=ivl1.upper_bound()*ivl2.upper_bound();
        return ApproximateInterval<F>(min(min(rll,rlu),min(rul,ruu)),max(max(rll,rlu),max(rul,ruu))); }
    friend ApproximateInterval<F> fma(ApproximateInterval<F> const& x1, ApproximateInterval<F> const& x2, ApproximateInterval<F> y);

    PR precision() const;
};

template<class PR> class DefineIntervalArithmeticOperations<Float<PR>>
    : public DefineIntervalArithmeticOperations<UpperBound<Float<PR>>>
{
  public:
    PR precision() const;

    using F=Float<PR>;
    friend Interval<UpperBound<F>> operator+(Bounds<F> const& x1, Interval<F> const& ivl2);
    friend Interval<UpperBound<F>> operator-(Bounds<F> const& x1, Interval<F> const& ivl2);
    friend Interval<UpperBound<F>> operator*(Bounds<F> const& x1, Interval<F> const& ivl2);
    friend Interval<UpperBound<F>> operator/(Bounds<F> const& x1, Interval<F> const& ivl2);
    friend Interval<UpperBound<F>> operator+(Interval<F> const& ivl1, Bounds<F> const& x2);
    friend Interval<UpperBound<F>> operator-(Interval<F> const& ivl1, Bounds<F> const& x2);
    friend Interval<UpperBound<F>> operator/(Interval<F> const& ivl1, Bounds<F> const& x2);
    friend Interval<UpperBound<F>> operator*(Interval<F> const& ivl1, Bounds<F> const& x2);

};

/*
template<class PR> class DefineIntervalArithmeticOperations<Float<PR>> {
  public:
    PR precision() const;
};
*/


template<class F> auto DefineIntervalArithmeticOperations<UpperBound<F>>::precision() const -> typename F::PrecisionType {
    Interval<UpperBound<F>>const& ivl=static_cast<Interval<UpperBound<F>>const&>(*this);
    return min(ivl.lower_bound().precision(),ivl.upper_bound().precision());
}

} // namespace Ariadne

#endif /* ARIADNE_INTERVAL_ARITHMETIC_HPP */


