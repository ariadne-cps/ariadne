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
#include "numeric/number.h"
#include "numeric/float.h"
#include "numeric/arithmetic.h"

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

template<class U> struct CentreTrait { typedef decltype(declval<U>()-declval<U>()) Type; };
template<class PR> struct CentreTrait<FloatUpperBound<PR>> { typedef FloatApproximation<PR> Type; };
template<class PR> struct CentreTrait<FloatLowerBound<PR>> { typedef FloatApproximation<PR> Type; };

template<class U> struct MidpointTrait { typedef decltype(max(-declval<U>(),declval<U>())) Type; };
template<class PR> struct MidpointTrait<FloatUpperBound<PR>> { typedef FloatValue<PR> Type; };
template<class PR> struct MidpointTrait<FloatLowerBound<PR>> { typedef FloatValue<PR> Type; };
template<class PR> struct MidpointTrait<FloatBounds<PR>> { typedef FloatValue<PR> Type; };
template<class PR> struct MidpointTrait<FloatBall<PR>> { typedef FloatValue<PR> Type; };

class UnitIntervalType;
class EmptyIntervalType;

template<class U> struct DeclareIntervalArithmeticOperations { };
template<class PR> struct DeclareIntervalArithmeticOperations<FloatUpperBound<PR>>
    : DeclareNumericOperations<FloatUpperInterval<PR>>
    , DeclareComparisonOperations<FloatUpperInterval<PR>,ValidatedKleenean>
    , DefineFieldOperators<FloatUpperInterval<PR>>
    , DefineComparisonOperators<FloatUpperInterval<PR>,ValidatedKleenean>
    , ProvideConvertedFieldOperations<FloatUpperInterval<PR>,FloatBounds<PR>>
    , ProvideConvertedComparisonOperations<FloatUpperInterval<PR>,FloatBounds<PR>,FloatUpperInterval<PR>,ValidatedKleenean>
    , ProvideConcreteGenericFieldOperations<FloatUpperInterval<PR>,ValidatedNumber>
    , ProvideConcreteGenericComparisonOperations<FloatUpperInterval<PR>,ValidatedNumber,ValidatedKleenean>
{
    typedef PR PrecisionType;
    friend FloatError<PR> mag(FloatUpperInterval<PR> const&);
    FloatUpperInterval<PR> create(ValidatedNumber const& y) const;
};

template<class PR> struct DeclareIntervalArithmeticOperations<FloatValue<PR>> : DeclareIntervalArithmeticOperations<FloatUpperBound<PR>> { };

//! \ingroup GeometryModule
//! \brief Intervals with upper endoint of type \a U.
//! \details
//! Not intended for use in basic interval arithmetic; represents a \em geometric rather than a \em numerical object.
template<class U> class Interval
    : public DeclareIntervalArithmeticOperations<U>
{
    typedef typename U::Paradigm P;
    typedef decltype(-declval<U>()) L;
    typedef typename CentreTrait<U>::Type C;
    typedef typename MidpointTrait<U>::Type M;
    typedef decltype(cast_positive(declval<U>()-declval<M>())) R;
    typedef decltype(cast_positive(declval<U>()-declval<L>())) W;
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
    //! \brief Construct from a unit interval.
    Interval(UnitIntervalType);
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
    template<class V, EnableIf<And<IsConstructible<L,V>,IsConstructible<U,V>>> = dummy>
        explicit Interval(const V& v) : _l(v), _u(v) { }
    //! \brief Assign a singleton interval from a number.
    template<class V, EnableIf<And<IsAssignable<L,V>,IsAssignable<U,V>>> = dummy>
        Interval<U>& operator=(const V& v) { _l=v; _u=v; return *this; }
    //! \brief Construct a singleton interval from a number.
    template<class Y, EnableIf<And<IsConstructible<L,Y,Precision64>,IsConstructible<U,Y,Precision64>>,Not<And<IsConstructible<L,Y>,IsConstructible<U,Y>>>> = dummy>
        explicit Interval(const Y& y) : _l(y,Precision64()), _u(y,Precision64()) { }

    //! \brief Convert from an interval of a different type.
    template<class UU, EnableIf<IsConvertible<UU,U>> = dummy>
        Interval(Interval<UU> const& x) : _l(x.lower()), _u(x.upper()) { }
    //! \brief Construct an interval with the given value.
    template<class UU, EnableIf<And<IsConstructible<U,UU>,Not<IsConvertible<UU,U>>>> =dummy>
        explicit Interval(Interval<UU> const& x) : _l(x.lower()), _u(x.upper()) { }
    //! \brief Construct an interval with the given value.
    template<class UU, class PR, EnableIf<IsConstructible<U,UU,PR>> =dummy>
        explicit Interval(Interval<UU> const& x, PR pr) : _l(x.lower(),pr), _u(x.upper(),pr) { }
    //! \brief Construct an interval with the given value.
    template<class UU, EnableIf<And<IsConstructible<U,UU,Precision64>,Not<IsConstructible<U,UU>>>> =dummy>
        explicit Interval(Interval<UU> const& x) : _l(x.lower(),Precision64()), _u(x.upper(),Precision64()) { }

    //! \brief Construct an interval with the given value.
    //! FIXME: Should be explicit, but this would clash with Box constructor from initializer list of double/Float64.
    template<class LL, class UU, EnableIf<And<IsConstructible<L,LL>,IsConstructible<U,UU>,Not<And<IsConvertible<LL,L>,IsConvertible<UU,U>>>>> =dummy>
        Interval(const LL& l, const UU& u) : _l(l), _u(u) { }
    template<class LL, class UU, EnableIf<And<IsConstructible<L,LL,Precision64>,IsConstructible<U,UU,Precision64>,Not<And<IsConstructible<L,LL>,IsConstructible<U,UU>>>>> =dummy>
        Interval(const LL& l, const UU& u) : _l(l,Precision64()), _u(u,Precision64()) { }

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
    auto is_bounded() const -> decltype(declval<U>()<declval<L>());
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

//! \related Interval \brief Make an interval with the given lower and upper bounds.
template<class L, class U> inline Interval<U> make_interval(L l, U u) { return Interval<U>(l,u); }

template<class U> inline auto is_empty(Interval<U> const& ivl) -> decltype(ivl.lower()>ivl.upper());
//! \related Interval \brief Test if the interval is empty.
template<class U> inline auto is_empty(Interval<U> const& ivl) -> decltype(ivl.lower()>ivl.upper());
//! \related Interval \brief Test if the interval is a singleton.
template<class U> inline auto is_singleton(Interval<U> const& ivl) -> decltype(ivl.lower()==ivl.upper());
//! \related Interval \brief Test if the interval is bounded.
template<class U> inline auto is_bounded(Interval<U> const& ivl) -> decltype(ivl.upper()<ivl.lower());

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


//! \related Float64ApproximationInterval \related Float64ExactInterval \brief Allows the over-approximating interval \a ivl to be treated as exact.
Interval<Float64Value> cast_exact(Interval<Float64Approximation> const& ivl);
Interval<FloatMPValue> cast_exact(Interval<FloatMPApproximation> const& ivl);
Interval<Float64Value> cast_exact_interval(Interval<Float64Approximation> const& ivl);
Interval<FloatMPValue> cast_exact_interval(Interval<FloatMPApproximation> const& ivl);

//! \related Float64UpperInterval \brief Computes a common refinement of \a ivl1 and \a ivl2 ivl.e. the intersection.
template<class PR> Interval<FloatUpperBound<PR>> refinement(Interval<FloatUpperBound<PR>> const& ivl1, Interval<FloatUpperBound<SelfType<PR>>> const& ivl2);
//! \related Float64UpperInterval \brief Tests if \a ivl1 provides a better over-approximation to the exact interval than \a ivl2.
//! ivl.e. \a ivl1 is a subset of \a ivl2.
template<class PR> bool refines(Interval<FloatUpperBound<PR>> const& ivl1, Interval<FloatUpperBound<SelfType<PR>>> const& ivl2);
//! \related Float64UpperInterval \brief Tests if two intervals have the same representation.
template<class PR> bool same(Interval<FloatUpperBound<PR>> const& ivl1, Interval<FloatUpperBound<SelfType<PR>>> const& ivl2);

//! \related Float64UpperInterval \related FloatBounds \brief Allows the over-approximating interval \a ivl to be treated an over-approximation to a single point.
Float64Bounds cast_singleton(Interval<Float64UpperBound> const& ivl);
Float64UpperInterval make_interval(Float64Bounds const& x);
FloatMPBounds cast_singleton(Interval<FloatMPUpperBound> const& ivl);
FloatMPUpperInterval make_interval(FloatMPBounds const& x);

//! \related Float64UpperInterval \brief An interval containing the given interval in its interior.
template<class PR> Interval<FloatUpperBound<PR>> widen(Interval<FloatValue<PR>> const& ivl);
template<class PR> Interval<FloatUpperBound<PR>> widen(Interval<FloatUpperBound<PR>> const& ivl);
template<class PR> Interval<FloatUpperBound<PR>> widen(Interval<FloatUpperBound<PR>> const& ivl, FloatUpperBound<PR> e);
template<class PR> Interval<FloatUpperBound<PR>> widen(Interval<FloatUpperBound<PR>> const& ivl, ValidatedUpperNumber e);
template<class PR> Interval<FloatValue<PR>> widen_domain(Interval<FloatUpperBound<PR>> const& ivl);
//! \related FloatLowerInterval<PR> \brief An interval contained in the interior of the given interval.
template<class PR> Interval<FloatLowerBound<PR>> narrow(Interval<FloatValue<PR>> const& ivl);
template<class PR> Interval<FloatLowerBound<PR>> narrow(Interval<FloatLowerBound<PR>> const& ivl);
template<class PR> Interval<FloatLowerBound<PR>> narrow(Interval<FloatLowerBound<PR>> const& ivl, FloatUpperBound<PR> e);
template<class PR> Interval<FloatUpperBound<PR>> narrow(Interval<FloatLowerBound<PR>> const& ivl, ValidatedUpperNumber e);

//! \related Interval \brief Read from an input stream.
InputStream& operator>>(InputStream&, Interval<Float64Value>&);

class UnitIntervalType : public ExactIntervalType {
  public:
    UnitIntervalType() : ExactIntervalType(-1,+1) { }
};

class EmptyIntervalType { };

} // namespace Ariadne

#include "interval.inl.h"

#endif

