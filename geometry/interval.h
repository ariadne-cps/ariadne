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

class UnitIntervalType;

struct SizeOne { operator SizeType() const { return 1u; } };

enum class SplitPart : char;

enum class SplitPart : char { LOWER, MIDDLE, UPPER };
InputStream& operator>>(InputStream& is, SplitPart& s);
OutputStream& operator<<(OutputStream& os, const SplitPart& s);

template<class U> struct SplitEndTypedef { typedef U Type; };
template<> struct SplitEndTypedef<UpperFloat64> { typedef ExactFloat64 Type; };
template<> struct SplitEndTypedef<BoundedFloat64> { typedef ExactFloat64 Type; };
template<class U> using SplitEndType = typename SplitEndTypedef<U>::Type;

class UnitIntervalType;
class EmptyIntervalType;


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
    Interval(EmptyIntervalType);
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
    template<class XX, EnableIf<And<IsConstructible<L,XX>,IsConstructible<U,XX>>> = dummy>
        explicit Interval(const XX& x) : _l(x), _u(x) { }
    //! \brief Assign a singleton interval from a number.
    template<class XX, EnableIf<And<IsAssignable<L,XX>,IsAssignable<U,XX>>> = dummy>
        Interval<U>& operator=(const XX& x) { _l=x; _u=x; return *this; }

    //! \brief Convert from an interval of a different type.
    template<class UU, EnableIf<IsConvertible<UU,U>> = dummy>
        Interval(Interval<UU> const& x) : _l(x.lower()), _u(x.upper()) { }
    //! \brief Construct an interval with the given value.
    template<class UU, EnableIf<And<IsConstructible<U,UU>,Not<IsConvertible<UU,U>>>> =dummy>
        explicit Interval(Interval<UU> const& x) : _l(x.lower()), _u(x.upper()) { }

    //! \brief Construct an interval with the given value.
    //! FIXME: Should be explicit, but this would clash with Box constructor from initializer list of double/Float64.
    template<class LL, class UU, EnableIf<And<IsConstructible<L,LL>,IsConstructible<U,UU>,Not<And<IsConvertible<LL,L>,IsConvertible<UU,U>>>>> =dummy>
        Interval(const LL& l, const UU& u) : _l(l), _u(u) { }

  public:
    //! \brief The dimension of the set; statically returns size one.
    SizeOne dimension() const;

    //! \brief The lower bound of the interval.
    LowerBoundType const& lower() const;
    //! \brief The upper bound of the interval.
    UpperBoundType const& upper() const;
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
    //! Test if the interval is a singleton.
    auto is_singleton() const -> decltype(declval<L>() == declval<U>());
    //! Test if the interval is bounded.
    auto is_bounded() const -> decltype(declval<U>()<infty);
  public:
    L _l; U _u;
};

//! \related Interval \brief Write to an output stream.
template<class U> OutputStream& operator<<(OutputStream& os, Interval<U> const& ivl);

template<class U> inline auto lower_bound(Interval<U> const& ivl) -> decltype(ivl.lower());
template<class U> inline auto upper_bound(Interval<U> const& ivl) -> decltype(ivl.upper());
template<class U> inline auto centre(Interval<U> const& ivl) -> decltype(ivl.centre());
template<class U> inline auto midpoint(Interval<U> const& ivl) -> decltype(ivl.midpoint());
template<class U> inline auto radius(Interval<U> const& ivl) -> decltype(ivl.radius());
template<class U> inline auto width(Interval<U> const& ivl) -> decltype(ivl.width());

//! \related Interval \brief Test if the interval is empty.
template<class U> inline auto is_empty(Interval<U> const& ivl) -> decltype(ivl.lower()>ivl.upper());
//! \related Interval \brief Test if the interval is a singleton.
template<class U> inline auto is_singleton(Interval<U> const& ivl) -> decltype(ivl.lower()==ivl.upper());
//! \related Interval \brief Test if the interval is bounded.
template<class U> inline auto is_bounded(Interval<U> const& ivl) -> decltype(ivl.upper()<infty);

//! \related Interval \brief Test if \a x1 is an element of the interval \a ivl2.
template<class U, class X> inline auto element(X const& x1, Interval<U> const& ivl2) -> decltype(ivl2.lower()<=x1 && ivl2.upper()>=x1);
//! \related Interval \brief Test if the interval \a ivl1 contains \a x2.
template<class U, class X> inline auto contains(Interval<U> const& ivl1, X const& x2) -> decltype(ivl1.lower()<=x2 && ivl1.upper()>=x2);
//! \related Interval \brief Test if the interval \a ivl1 is equal to \a ivl2.
template<class U> inline auto equal(Interval<U> const& ivl1, Interval<U> const& ivl2) -> decltype(ivl1.upper()==ivl2.upper());
//! \related Interval \brief Test if the interval \a ivl1 is a subset of \a ivl2.
template<class U1, class U2> inline auto subset(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()<=ivl2.upper());
//! \related Interval \brief Test if the interval \a ivl1 is a superset of \a ivl2.
template<class U1, class U2> inline auto superset(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()>=ivl2.upper());
//! \related Interval \brief Test if the interval \a ivl1 is disjoint from \a ivl2. Returns \c false even if the two intervals only have an endpoint in common.
template<class U1, class U2> inline auto disjoint(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()<ivl2.lower());
//! \related Interval \brief Test if the interval \a ivl1 intersects \a ivl2. Returns \c true even if the two intervals only have an endpoint in common.
template<class U1, class U2> inline auto intersect(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()>=ivl2.lower());

//! \related Interval \brief Test if the closed interval \a ivl1 is disjoint from the closed interval \a ivl2.
//! Returns \c false if the two intervals only have an endpoint in common.
template<class U1, class U2> inline auto separated(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()<ivl2.lower());
//! \related Interval \brief Test if the interval \a ivl1 overlaps \a ivl2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
template<class U1, class U2> inline auto overlap(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()>ivl2.lower());
//! \related Interval \brief Test if the (closed) interval \a ivl1 is a subset of the interior of \a ivl2.
template<class U1, class U2> inline auto inside(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()<ivl2.upper());
//! \related Interval \brief Test if the interior of the interval \a ivl1 is a superset of the (closed) interval \a ivl2.
template<class U1, class U2> inline auto covers(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> decltype(ivl1.upper()>ivl2.upper());

//! \related Interval \brief The intersection of two intervals.
template<class U1, class U2> inline auto intersection(Interval<U1> const& ivl1, Interval<U2> const& ivl2) -> Interval<decltype(min(declval<U1>(),declval<U2>()))>;

//! \related Interval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
template<class U1, class U2> inline auto hull(Interval<U1> const& ivl1, Interval<U2> const& ivl2) ->  Interval<decltype(max(declval<U1>(),declval<U2>()))>;

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
template<class U, class X> inline auto hull(Interval<U> const& ivl1, X x2) -> Interval<decltype(max(declval<U>(),declval<X>()))>;

//! \related Interval \brief The hull of a point and an interval, equal to the smallest interval containing both.
template<class U, class X> inline auto hull(X x1, Interval<U> const& ivl2) -> Interval<decltype(max(declval<U>(),declval<X>()))>;

//! \related Interval \brief Split an interval into its lower, middle or upper half.
template<class U> inline auto split(Interval<U> const& ivl1, SplitPart lmu) -> Interval<U>;


//! \related Interval \brief Equality operator. Tests equality of intervals as geometric objects given information on endpoints.
template<class U> inline auto operator==(Interval<U> const& ivl1, Interval<U> const& ivl2) -> decltype(ivl1.upper()==ivl2.upper());
//! \related Interval \brief Inequality operator.
template<class U> inline auto operator!=(Interval<U> const& ivl1, Interval<U> const& ivl2) -> decltype(ivl1.upper()!=ivl2.upper());


//! \related ApproximateFloat64Interval \related ExactFloat64Interval \brief Allows the over-approximating interval \a ivl to be treated as exact.
Interval<ExactFloat64> cast_exact(Interval<ApproximateFloat64> const& ivl);
Interval<ExactFloat64> cast_exact_interval(Interval<ApproximateFloat64> const& ivl);

//! \related UpperFloat64Interval \brief Computes a common refinement of \a ivl1 and \a ivl2 ivl.e. the intersection.
Interval<UpperFloat64> refinement(Interval<UpperFloat64> const& ivl1, Interval<UpperFloat64> const& ivl2);
//! \related UpperFloat64Interval \brief Tests if \a ivl1 provides a better over-approximation to the exact interval than \a ivl2.
//! ivl.e. \a ivl1 is a subset of \a ivl2.
bool refines(Interval<UpperFloat64> const& ivl1, Interval<UpperFloat64> const& ivl2);
//! \related UpperFloat64Interval \brief Tests if two intervals have the same representation.
bool same(Interval<UpperFloat64> const& ivl1, Interval<UpperFloat64> const& ivl2);

//! \related UpperFloat64Interval \related BoundedFloat \brief Allows the over-approximating interval \a ivl to be treated an over-approximation to a single point.
BoundedFloat64 cast_singleton(Interval<UpperFloat64> const& ivl);

//! \related UpperFloat64Interval \brief An interval containing the given interval in its interior.
Interval<UpperFloat64> widen(Interval<UpperFloat64> const& ivl);
Interval<UpperFloat64> widen(Interval<UpperFloat64> const& ivl, UpperFloat64 e);
Interval<ExactFloat64> widen_domain(Interval<UpperFloat64> const& ivl);
//! \related LowerFloat64Interval \brief An interval contained in the interior of the given interval.
Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& ivl);
Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& ivl, UpperFloat64 e);

//! \related Interval \brief Read from an input stream.
InputStream& operator>>(InputStream&, Interval<ExactFloat64>&);

class UnitIntervalType : public ExactIntervalType {
  public:
    UnitIntervalType() : ExactIntervalType(-1,+1) { }
};

class EmptyIntervalType { };

} // namespace Ariadne

#include "interval.inl.h"

namespace Ariadne {

inline BoundedFloat64 cast_singleton(UpperFloat64Interval const& ivl) {
    return BoundedFloat64(ivl.lower(),ivl.upper()); }
inline UpperFloat64Interval make_interval(BoundedFloat64 const& x) {
    return UpperFloat64Interval(x.lower(),x.upper()); }

inline UpperFloat64Interval max(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return make_interval(max(cast_singleton(ivl1),cast_singleton(ivl2))); }
inline UpperFloat64Interval min(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return make_interval(min(cast_singleton(ivl1),cast_singleton(ivl2))); }
inline UpperFloat64Interval abs(UpperFloat64Interval ivl) {
    return make_interval(abs(cast_singleton(ivl))); }
inline UpperFloat64Interval pos(UpperFloat64Interval ivl) {
    return make_interval(pos(cast_singleton(ivl))); }
inline UpperFloat64Interval neg(UpperFloat64Interval ivl) {
    return make_interval(neg(cast_singleton(ivl))); }
inline UpperFloat64Interval sqr(UpperFloat64Interval ivl) {
    return make_interval(sqr(cast_singleton(ivl))); }
inline UpperFloat64Interval rec(UpperFloat64Interval ivl) {
    return make_interval(rec(cast_singleton(ivl))); }

inline UpperFloat64Interval add(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return make_interval(add(cast_singleton(ivl1),cast_singleton(ivl2))); }
inline UpperFloat64Interval sub(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return make_interval(sub(cast_singleton(ivl1),cast_singleton(ivl2))); }
inline UpperFloat64Interval mul(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return make_interval(mul(cast_singleton(ivl1),cast_singleton(ivl2))); }
inline UpperFloat64Interval div(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return make_interval(div(cast_singleton(ivl1),cast_singleton(ivl2))); }

inline UpperFloat64Interval pow(UpperFloat64Interval ivl, Nat m) {
    return make_interval(pow(cast_singleton(ivl),m)); }
inline UpperFloat64Interval pow(UpperFloat64Interval ivl, Int n) {
    return make_interval(pow(cast_singleton(ivl),n)); }

inline UpperFloat64Interval sqrt(UpperFloat64Interval ivl) {
    return make_interval(sqrt(cast_singleton(ivl))); }
inline UpperFloat64Interval exp(UpperFloat64Interval ivl) {
    return make_interval(exp(cast_singleton(ivl))); }
inline UpperFloat64Interval log(UpperFloat64Interval ivl) {
    return make_interval(log(cast_singleton(ivl))); }
inline UpperFloat64Interval sin(UpperFloat64Interval ivl) {
    return make_interval(sin(cast_singleton(ivl))); }
inline UpperFloat64Interval cos(UpperFloat64Interval ivl) {
    return make_interval(cos(cast_singleton(ivl))); }
inline UpperFloat64Interval tan(UpperFloat64Interval ivl) {
    return make_interval(tan(cast_singleton(ivl))); }
inline UpperFloat64Interval asin(UpperFloat64Interval ivl) {
    return make_interval(asin(cast_singleton(ivl))); }
inline UpperFloat64Interval acos(UpperFloat64Interval ivl) {
    return make_interval(acos(cast_singleton(ivl))); }
inline UpperFloat64Interval atan(UpperFloat64Interval ivl) {
    return make_interval(atan(cast_singleton(ivl))); }

inline PositiveUpperFloat64 mag(UpperFloat64Interval ivl) {
    return mag(cast_singleton(ivl)); }
inline LowerFloat64 mig(UpperFloat64Interval ivl) {
    return mig(cast_singleton(ivl)); }

inline UpperFloat64Interval operator+(const UpperFloat64Interval& ivl) { return pos(ivl); }
inline UpperFloat64Interval operator-(const UpperFloat64Interval& ivl) { return neg(ivl); }
inline UpperFloat64Interval operator+(const UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { return add(ivl1,ivl2); }
inline UpperFloat64Interval operator-(const UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { return sub(ivl1,ivl2); }
inline UpperFloat64Interval operator*(const UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { return mul(ivl1,ivl2); }
inline UpperFloat64Interval operator/(const UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { return div(ivl1,ivl2); };
inline UpperFloat64Interval& operator+=(UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { ivl1=add(ivl1,ivl2); return ivl1; }
inline UpperFloat64Interval& operator-=(UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { ivl1=sub(ivl1,ivl2); return ivl1; }
inline UpperFloat64Interval& operator*=(UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { ivl1=mul(ivl1,ivl2); return ivl1; }
inline UpperFloat64Interval& operator/=(UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) { ivl1=div(ivl1,ivl2); return ivl1; }

inline Bool operator==(const UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) {
    return ivl1.lower().raw() == ivl2.lower().raw() && ivl1.upper().raw() == ivl2.upper().raw(); }
inline Kleenean operator!=(const UpperFloat64Interval& ivl1, const UpperFloat64Interval& ivl2) {
    return ivl1.lower().raw() != ivl2.lower().raw() || ivl1.upper().raw() != ivl2.upper().raw(); }
inline Kleenean operator<=(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return cast_singleton(ivl1) <= cast_singleton(ivl2); }
inline Kleenean operator>=(UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return cast_singleton(ivl1) >= cast_singleton(ivl2); }
inline Kleenean operator< (UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return cast_singleton(ivl1) <  cast_singleton(ivl2); }
inline Kleenean operator> (UpperFloat64Interval ivl1, UpperFloat64Interval ivl2) {
    return cast_singleton(ivl1) >  cast_singleton(ivl2); }

// Mixed operations
inline UpperFloat64Interval operator+(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1+make_interval(x2); }
inline UpperFloat64Interval operator-(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1-make_interval(x2); }
inline UpperFloat64Interval operator*(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1*make_interval(x2); }
inline UpperFloat64Interval operator/(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1/make_interval(x2); }
inline UpperFloat64Interval operator+(BoundedFloat64 x1, UpperFloat64Interval ivl2) { return make_interval(x1)+ivl2; }
inline UpperFloat64Interval operator-(BoundedFloat64 x1, UpperFloat64Interval ivl2) { return make_interval(x1)-ivl2; }
inline UpperFloat64Interval operator*(BoundedFloat64 x1, UpperFloat64Interval ivl2) { return make_interval(x1)*ivl2; }
inline UpperFloat64Interval operator/(BoundedFloat64 x1, UpperFloat64Interval ivl2) { return make_interval(x1)/ivl2; }
inline UpperFloat64Interval& operator+=(UpperFloat64Interval& ivl1, BoundedFloat64 x2) { return ivl1+=make_interval(x2); }
inline UpperFloat64Interval& operator-=(UpperFloat64Interval& ivl1, BoundedFloat64 x2) { return ivl1-=make_interval(x2); }
inline UpperFloat64Interval& operator*=(UpperFloat64Interval& ivl1, BoundedFloat64 x2) { return ivl1*=make_interval(x2); }
inline UpperFloat64Interval& operator/=(UpperFloat64Interval& ivl1, BoundedFloat64 x2) { return ivl1/=make_interval(x2); }
inline Kleenean operator==(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1==make_interval(x2); }
inline Kleenean operator!=(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1!=make_interval(x2); }
inline Kleenean operator<=(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1<=make_interval(x2); }
inline Kleenean operator>=(UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1>=make_interval(x2); }
inline Kleenean operator< (UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1< make_interval(x2); }
inline Kleenean operator> (UpperFloat64Interval ivl1, BoundedFloat64 x2) { return ivl1> make_interval(x2); }


} // namespace Ariadne
#endif

