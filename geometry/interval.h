/***************************************************************************
 *            geometry/interval.h
 *
 *  Copyright 2013-14  Pieter Collins
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
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file geometry/interval.h
 *  \brief
 */



#ifndef ARIADNE_INTERVAL_H
#define ARIADNE_INTERVAL_H

#include "utility/module.h"

#include "numeric/logical.h"
#include "numeric/float.h"

#include "interval.decl.h"

namespace Ariadne {

//! \defgroup GeometryModule Geometry Module
//!  \brief Classes describing sets in Euclidean space.

class UnitInterval;

struct SizeOne { operator SizeType() const { return 1u; } };

enum class SplitPart : char;

enum class SplitPart : char { LOWER, MIDDLE, UPPER };
InputStream& operator>>(InputStream& is, SplitPart& s);
OutputStream& operator<<(OutputStream& os, const SplitPart& s);

template<class U> struct SplitEndTypedef { typedef U Type; };
template<> struct SplitEndTypedef<UpperFloat64> { typedef ExactFloat64 Type; };
template<> struct SplitEndTypedef<BoundedFloat64> { typedef ExactFloat64 Type; };
template<class U> using SplitEndType = typename SplitEndTypedef<U>::Type;

class UnitInterval;
class EmptyInterval;


//! \ingroup GeometryModule
//! \brief Intervals with upper endoint of type \a U.
//! \details
//! Not intended for use in basic interval arithmetic; represents a \em geometric rather than a \em numerical object.
template<class U> class Interval {
    typedef typename U::Paradigm P;
    typedef decltype(-declval<U>()) L;
    typedef decltype(declval<U>()+declval<L>()) C;
    typedef typename SplitEndTypedef<U>::Type M;
    typedef decltype(declval<U>()-declval<M>()) R;
    typedef decltype(declval<U>()-declval<L>()) W;

  public:
    //! \brief The computational paradigm used by the interval.
    typedef P Paradigm;
    //! \brief The type of the lower bound of the interval.
    typedef L LowerBoundType;
    //! \brief The type of the upper bound of the interval.
    typedef U UpperBoundType;
    //! \brief The type returned by the midpoint() method. The midpoint need not be the exact mid-point
    //! of the interval, but must be convertible to both the upper and lower bound types.
    typedef C CentreType;
    //! \brief The type returned by the midpoint() method. The midpoint need not be the exact mid-point
    //! of the interval, but must be convertible to both the upper and lower bound types.
    typedef M MidpointType;
    //! \brief The type returned by the radius() methods.
    typedef R RadiusType;
    //! \brief The type returned by the width() methods.
    typedef W WidthType;
  public:
    // Deprecated
    typedef Interval<U> NumericType;
  public:
    //! \brief Create the interval which is the single point 0.
    Interval<U> create_zero() const;

    //! \brief Construct an empty interval.
    explicit Interval();
    //! \brief Construct from an empty interval.
    Interval(EmptyInterval);
    //! \brief Construct an interval with the given lower and upper bounds.
    Interval(LowerBoundType l, UpperBoundType u);

    //! \brief Construct an empty interval.
    static Interval<U> empty_interval();
    //! \brief Construct a unit interval.
    static Interval<U> unit_interval();
    //! \brief Construct a biinfinite interval.
    static Interval<U> biinfinite_interval();
    //! \brief Construct a zero interval.
    static Interval<U> singleton_interval(MidpointType x);

    //! \brief Construct a singleton interval from a number.
    //template<class X, EnableIf<And<IsNumber<X>,IsConvertible<X,L>,IsConvertible<X,U>>> = dummy>
    template<class XX, EnableIf<And<IsConstructible<L,XX>,IsConstructible<U,XX>>> = dummy>
        explicit Interval(const XX& x) : _l(x), _u(x) { }

    //! \brief Assign a singleton interval from a number.
    template<class XX, EnableIf<And<IsConvertible<XX,L>,IsConvertible<XX,U>>> = dummy>
        Interval<U>& operator=(const XX& x) { _l=x; _u=x; return *this; }
    //! \brief Convert from an interval of a different type.
    template<class UU, EnableIf<IsConvertible<UU,U>> = dummy>
        Interval(Interval<UU> const& x) : _l(x.lower()), _u(x.upper()) { }

    //! \brief Construct an interval with the given value.
    template<class UU, EnableIf<And<IsConstructible<U,UU>,Not<IsConvertible<UU,U>>>> =dummy>
        explicit Interval(Interval<UU> const& x) : _l(x.lower()), _u(x.upper()) { }

    //! \brief Construct an interval with the given value.
    template<class LL, class UU, EnableIf<And<IsConstructible<L,LL>,IsConstructible<U,UU>,Not<And<IsConvertible<LL,L>,IsConvertible<UU,U>>>>> =dummy>
        // FIXME: Should be explicit
        Interval(const LL& l, const UU& u) : _l(l), _u(u) { }

  public:

    //! \brief The dimension of the set; statically returns size one.
    SizeOne dimension() const;

    //! \brief The lower bound of the interval.
    LowerBoundType const& lower() const;
    LowerBoundType& lower();
    //! \brief The upper bound of the interval.
    UpperBoundType const& upper() const;
    UpperBoundType& upper();
    //! \brief The midpoint of the interval. Need not be the exact midpoint. \sa radius()
    MidpointType midpoint() const;
    //! \brief The centre of the interval, computed with whatever rounding is appropriate.
    CentreType centre() const;
    //! \brief The radius of the interval. In validated computation, must be an upper bound for the both \c u-m and \c m-l, where \c m is the midpoint. \sa midpoint()
    RadiusType radius() const;
    //! \brief The width of the interval.
    WidthType width() const;

    void set_lower(LowerBoundType l);
    void set_upper(UpperBoundType u);
    void set(LowerBoundType l, UpperBoundType u);

    //! Test if the interval is empty.
    auto is_empty() const -> decltype(declval<L>() > declval<U>());
    auto empty() const -> decltype(declval<L>() > declval<U>());
    //! Test if the interval is a singleton.
    auto is_singleton() const -> decltype(declval<L>() == declval<U>());
    //! Test if the interval is bounded.
    auto is_bounded() const -> decltype(declval<U>()<infty);
    auto bounded() const -> decltype(declval<U>()<infty);
  public:
    L _l; U _u;
};

//! \related Interval \brief Write to an output stream.
template<class U> OutputStream& operator<<(OutputStream& os, Interval<U> const& ivl);
//! \related Interval \brief Read from an output stream.  \deprecated
// template<class U> InputStream& operator>>(InputStream& os, Interval<U>& ivl);

template<class U> inline auto midpoint(Interval<U> const& i) -> decltype(i.midpoint());
template<class U> inline auto radius(Interval<U> const& i) -> decltype(i.radius());
template<class U> inline auto width(Interval<U> const& i) -> decltype(i.width());

//! \related Interval \brief Test if the interval is empty.
template<class U> inline auto is_empty(Interval<U> const& i) -> decltype(i.lower()>i.upper());
//! \related Interval \brief Test if the interval is empty.
template<class U> inline auto is_singleton(Interval<U> const& i) -> decltype(i.lower()==i.upper());
//! \related Interval \brief Test if the interval is bounded.
template<class U> inline auto is_bounded(Interval<U> const& i) -> decltype(i.upper()<i.midpoint());

//! \related Interval \brief Test if \a x1 is an element of the interval \a i2.
template<class U, class X> inline auto element(X const& x1, Interval<U> const& i2) -> decltype(i2.lower()<=x1 && i2.upper()>=x1);
//! \related Interval \brief Test if the interval \a i1 contains \a x2.
template<class U, class X> inline auto contains(Interval<U> const& i1, X const& x2) -> decltype(i1.lower()<=x2 && i1.upper()>=x2);
//! \related Interval \brief Test if the interval \a i1 is equal to \a i2.
template<class U> inline auto equal(Interval<U> const& i1, Interval<U> const& i2) -> decltype(i1.upper()==i2.upper());
//! \related Interval \brief Test if the interval \a i1 is a subset of \a i2.
template<class U1, class U2> inline auto subset(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()<=i2.upper());
//! \related Interval \brief Test if the interval \a i1 is a superset of \a i2.
template<class U1, class U2> inline auto superset(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()>=i2.upper());
//! \related Interval \brief Test if the interval \a i1 is disjoint from \a i2. Returns \c false even if the two intervals only have an endpoint in common.
template<class U1, class U2> inline auto disjoint(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()<i2.lower());
//! \related Interval \brief Test if the interval \a i1 intersects \a i2. Returns \c true even if the two intervals only have an endpoint in common.
template<class U1, class U2> inline auto intersect(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()>=i2.lower());

//! \related Interval \brief Test if the closed interval \a i1 is disjoint from the closed interval \a i2.
//! Returns \c false if the two intervals only have an endpoint in common.
template<class U1, class U2> inline auto separated(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()<i2.lower());
//! \related Interval \brief Test if the interval \a i1 overlaps \a i2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
template<class U1, class U2> inline auto overlap(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()>i2.lower());
//! \related Interval \brief Test if the (closed) interval \a i1 is a subset of the interior of \a i2.
template<class U1, class U2> inline auto inside(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()<i2.upper());
//! \related Interval \brief Test if the interior of the interval \a i1 is a superset of the (closed) interval \a i2.
template<class U1, class U2> inline auto covers(Interval<U1> const& i1, Interval<U2> const& i2) -> decltype(i1.upper()>i2.upper());

//! \related Interval \brief The intersection of two intervals.
template<class U1, class U2> inline auto intersection(Interval<U1> const& i1, Interval<U2> const& i2) -> Interval<decltype(max(declval<U1>(),declval<U2>()))>;

//! \related Interval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
template<class U1, class U2> inline auto hull(Interval<U1> const& i1, Interval<U2> const& i2) ->  Interval<decltype(max(declval<U1>(),declval<U2>()))>;

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
template<class U, class X> inline auto hull(Interval<U> const& i1, X x2) -> Interval<decltype(max(declval<U>(),declval<X>()))>;

//! \related Interval \brief The hull of a point and an interval, equal to the smallest interval containing both.
template<class U, class X> inline auto hull(X x1, Interval<U> const& i2) -> Interval<decltype(max(declval<U>(),declval<X>()))>;

//! \related Interval \brief Split an interval into its lower, middle or upper half.
template<class U> inline auto split(Interval<U> const& i1, SplitPart lmu) -> Interval<U>;


// Standard equality operators
//! \related ExactFloatInterval \brief Equality operator.Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
bool operator==(Interval<ExactFloat64> const& i1, Interval<ExactFloat64> const& i2);
//! \related ExactFloatInterval \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true.
bool operator!=(Interval<ExactFloat64> const& i1, Interval<ExactFloat64> const& i2);


//! \related UpperFloatInterval \related ExactFloatInterval \brief Allows the over-approximating interval \a i to be treated as exact.
Interval<ExactFloat64> make_exact(Interval<ApproximateFloat64> const& i);
Interval<ExactFloat64> make_exact_interval(Interval<ApproximateFloat64> const& i);
//! \related UpperFloatInterval \brief Computes a common refinement of \a i1 and \a i2 i.e. the intersection.
Interval<UpperFloat64> refinement(Interval<UpperFloat64> const& i1, Interval<UpperFloat64> const& i2);
//! \related UpperFloatInterval \brief Tests if \a i1 provides a better over-approximation to the exact interval than \a i2.
//! i.e. \a i1 is a subset of \a i2.
bool refines(Interval<UpperFloat64> const& i1, Interval<UpperFloat64> const& i2);
//! \related UpperFloatInterval \brief Computes a common refinement of \a i1 and \a i2 i.e. the intersection.
bool same(Interval<UpperFloat64> const& i1, Interval<UpperFloat64> const& i2);

//! \related UpperFloatInterval \related BoundFloat \brief Allows the over-approximating interval \a i to be treated an over-approximation to a single point.
BoundFloat64 make_singleton(Interval<UpperFloat64> const& i);

//! \related UpperFloatInterval \brief An interval containing the given interval in its interior.
Interval<UpperFloat64> widen(Interval<UpperFloat64> const& i);
Interval<UpperFloat64> widen(Interval<UpperFloat64> const& i, UpperFloat64 e);
//! \related LowerFloatInterval \brief An interval contained in the interior of the given interval.
Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& i);
Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& i, UpperFloat64 e);

InputStream& operator>>(InputStream&, Interval<ExactFloat64>&);

class UnitInterval
    : public ExactInterval
{
  public:
    UnitInterval() : ExactInterval(-1,+1) { }
};

class EmptyInterval { };

} // namespace Ariadne

#include "interval.inl.h"

namespace Ariadne {

inline Bool refines(UpperInterval const& ivl1, UpperInterval const& ivl2) {
    return ivl1.lower().raw()>=ivl2.lower().raw() && ivl1.upper().raw()<=ivl2.upper().raw(); }

inline UpperInterval max(UpperInterval i1, UpperInterval i2) {
    return make_interval(max(make_singleton(i1),make_singleton(i2))); }
inline UpperInterval min(UpperInterval i1, UpperInterval i2) {
    return make_interval(min(make_singleton(i1),make_singleton(i2))); }
inline UpperInterval abs(UpperInterval i) {
    return make_interval(abs(make_singleton(i))); }
inline UpperInterval pos(UpperInterval i) {
    return make_interval(pos(make_singleton(i))); }
inline UpperInterval neg(UpperInterval i) {
    return make_interval(neg(make_singleton(i))); }
inline UpperInterval sqr(UpperInterval i) {
    return make_interval(sqr(make_singleton(i))); }
inline UpperInterval rec(UpperInterval i) {
    return make_interval(rec(make_singleton(i))); }

inline UpperInterval add(UpperInterval i1, UpperInterval i2) {
    return make_interval(add(make_singleton(i1),make_singleton(i2))); }
inline UpperInterval sub(UpperInterval i1, UpperInterval i2) {
    return make_interval(sub(make_singleton(i1),make_singleton(i2))); }
inline UpperInterval mul(UpperInterval i1, UpperInterval i2) {
    return make_interval(mul(make_singleton(i1),make_singleton(i2))); }
inline UpperInterval div(UpperInterval i1, UpperInterval i2) {
    return make_interval(div(make_singleton(i1),make_singleton(i2))); }

inline UpperInterval pow(UpperInterval i, Nat m) {
    return make_interval(pow(make_singleton(i),m)); }
inline UpperInterval pow(UpperInterval i, Int n) {
    return make_interval(pow(make_singleton(i),n)); }

inline UpperInterval sqrt(UpperInterval i) {
    return make_interval(sqrt(make_singleton(i))); }
inline UpperInterval exp(UpperInterval i) {
    return make_interval(exp(make_singleton(i))); }
inline UpperInterval log(UpperInterval i) {
    return make_interval(log(make_singleton(i))); }
inline UpperInterval sin(UpperInterval i) {
    return make_interval(sin(make_singleton(i))); }
inline UpperInterval cos(UpperInterval i) {
    return make_interval(cos(make_singleton(i))); }
inline UpperInterval tan(UpperInterval i) {
    return make_interval(tan(make_singleton(i))); }
inline UpperInterval asin(UpperInterval i) {
    return make_interval(asin(make_singleton(i))); }
inline UpperInterval acos(UpperInterval i) {
    return make_interval(acos(make_singleton(i))); }
inline UpperInterval atan(UpperInterval i) {
    return make_interval(atan(make_singleton(i))); }

inline PositiveUpperFloat64 mag(UpperInterval i) {
    return mag(make_singleton(i)); }
inline LowerFloat64 mig(UpperInterval i) {
    return mig(make_singleton(i)); }

inline UpperInterval operator+(const UpperInterval& i) { return pos(i); }
inline UpperInterval operator-(const UpperInterval& i) { return neg(i); }
inline UpperInterval operator+(const UpperInterval& i1, const UpperInterval& i2) { return add(i1,i2); }
inline UpperInterval operator-(const UpperInterval& i1, const UpperInterval& i2) { return sub(i1,i2); }
inline UpperInterval operator*(const UpperInterval& i1, const UpperInterval& i2) { return mul(i1,i2); }
inline UpperInterval operator/(const UpperInterval& i1, const UpperInterval& i2) { return div(i1,i2); };
inline UpperInterval& operator+=(UpperInterval& i1, const UpperInterval& i2) { i1=add(i1,i2); return i1; }
inline UpperInterval& operator-=(UpperInterval& i1, const UpperInterval& i2) { i1=sub(i1,i2); return i1; }
inline UpperInterval& operator*=(UpperInterval& i1, const UpperInterval& i2) { i1=mul(i1,i2); return i1; }
inline UpperInterval& operator/=(UpperInterval& i1, const UpperInterval& i2) { i1=div(i1,i2); return i1; }

inline Bool same(const UpperInterval& i1, const UpperInterval& i2) {
    return i1.lower().raw() == i2.lower().raw() && i1.upper().raw() == i2.upper().raw(); }

inline Bool operator==(const UpperInterval& i1, const UpperInterval& i2) {
    return i1.lower().raw() == i2.lower().raw() && i1.upper().raw() == i2.upper().raw(); }
inline Tribool operator!=(const UpperInterval& i1, const UpperInterval& i2) {
    return i1.lower().raw() != i2.lower().raw() || i1.upper().raw() != i2.upper().raw(); }
inline Tribool operator<=(UpperInterval i1, UpperInterval i2) {
    return make_singleton(i1) <= make_singleton(i2); }
inline Tribool operator>=(UpperInterval i1, UpperInterval i2) {
    return make_singleton(i1) >= make_singleton(i2); }
inline Tribool operator< (UpperInterval i1, UpperInterval i2) {
    return make_singleton(i1) <  make_singleton(i2); }
inline Tribool operator> (UpperInterval i1, UpperInterval i2) {
    return make_singleton(i1) >  make_singleton(i2); }

// Mixed operations
inline UpperInterval operator+(UpperInterval i1, ValidatedFloat64 x2) { return i1+make_interval(x2); }
inline UpperInterval operator-(UpperInterval i1, ValidatedFloat64 x2) { return i1-make_interval(x2); }
inline UpperInterval operator*(UpperInterval i1, ValidatedFloat64 x2) { return i1*make_interval(x2); }
inline UpperInterval operator/(UpperInterval i1, ValidatedFloat64 x2) { return i1/make_interval(x2); }
inline UpperInterval operator+(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)+i2; }
inline UpperInterval operator-(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)-i2; }
inline UpperInterval operator*(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)*i2; }
inline UpperInterval operator/(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)/i2; }
inline UpperInterval& operator+=(UpperInterval& i1, ValidatedFloat64 x2) { return i1+=make_interval(x2); }
inline UpperInterval& operator-=(UpperInterval& i1, ValidatedFloat64 x2) { return i1-=make_interval(x2); }
inline UpperInterval& operator*=(UpperInterval& i1, ValidatedFloat64 x2) { return i1*=make_interval(x2); }
inline UpperInterval& operator/=(UpperInterval& i1, ValidatedFloat64 x2) { return i1/=make_interval(x2); }
inline Tribool operator==(UpperInterval i1, ValidatedFloat64 x2) { return i1==make_interval(x2); }
inline Tribool operator!=(UpperInterval i1, ValidatedFloat64 x2) { return i1!=make_interval(x2); }
inline Tribool operator<=(UpperInterval i1, ValidatedFloat64 x2) { return i1<=make_interval(x2); }
inline Tribool operator>=(UpperInterval i1, ValidatedFloat64 x2) { return i1>=make_interval(x2); }
inline Tribool operator< (UpperInterval i1, ValidatedFloat64 x2) { return i1< make_interval(x2); }
inline Tribool operator> (UpperInterval i1, ValidatedFloat64 x2) { return i1> make_interval(x2); }


} // namespace Ariadne
#endif

