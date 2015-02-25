/***************************************************************************
 *            geometry/interval.inl.h
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

/*! \file geometry/interval.inl.h
 *  \brief
 */

namespace Ariadne {

template<class M> inline M make_split_point(M const& m) { return m; }
template<class PR> inline Float<Exact,PR> make_split_point(Float<Approximate,PR> const& am) { return make_exact(am); }
template<class PR> inline Float<Exact,PR> make_split_point(Float<Bounded,PR> const& vm) { return vm.value(); }
template<class PR> inline Float<Exact,PR> make_split_point(Float<Metric,PR> const& bm) { return bm.value(); }


template<class U> Interval<U>::Interval() : _l(+infty), _u(-infty) { }
template<class U> Interval<U>::Interval(EmptyInterval) : _l(+infty), _u(-infty) { }
template<class U> Interval<U>::Interval(LowerBoundType l, UpperBoundType u) : _l(l), _u(u) { }

template<class U> Interval<U> Interval<U>::create_zero() const { return Interval<U>(0,0); }
template<class U> SizeOne Interval<U>::dimension() const  { return SizeOne(); }
template<class U> auto Interval<U>::lower() const -> LowerBoundType const& { return _l; }
template<class U> auto Interval<U>::upper() const -> UpperBoundType const& { return _u; }
template<class U> auto Interval<U>::lower() -> LowerBoundType& { return _l; }
template<class U> auto Interval<U>::upper() -> UpperBoundType& { return _u; }
template<class U> auto Interval<U>::centre() const -> CentreType { return (_l+_u)/2u; }
template<class U> auto Interval<U>::midpoint() const -> MidpointType { auto m=((_l+_u)/2); return make_split_point(m); }
template<class U> auto Interval<U>::radius() const -> RadiusType { return max(this->upper()-this->midpoint(),this->midpoint()-this->lower()); }
template<class U> auto Interval<U>::width() const -> WidthType { return this->upper()-this->lower(); }

template<class U> Interval<U> Interval<U>::empty_interval() { return Interval<U>(+infty,-infty); }
template<class U> Interval<U> Interval<U>::unit_interval() { return Interval<U>(-1,+1); }

template<class U> auto Interval<U>::empty() const -> decltype(declval<L>()>declval<U>()) { return this->_l > this->_u; }
template<class U> auto Interval<U>::is_empty() const -> decltype(declval<L>()>declval<U>()) { return this->_l > this->_u; }
template<class U> auto Interval<U>::is_bounded() const -> decltype(declval<U>()<infty) { return this->_l > -infty && this->_u < +infty; }
template<class U> auto Interval<U>::is_singleton() const -> decltype(declval<L>() == declval<U>()) { return this->_l == this->_u; }

template<class U> auto Interval<U>::set_lower(LowerBoundType l) -> void { _l=l; }
template<class U> auto Interval<U>::set_upper(UpperBoundType u) -> void { _u=u; }
template<class U> auto Interval<U>::set(LowerBoundType l, UpperBoundType u) -> void { _l=l; _u=u; }

//! \related Interval \brief Write to an output stream
template<class U> inline OutputStream& operator<<(OutputStream& os, Interval<U> const& ivl) {
    return os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
}

template<class UB> inline auto midpoint(Interval<UB> const& i) -> decltype(i.midpoint()) { return i.midpoint(); }
template<class UB> inline auto radius(Interval<UB> const& i) -> decltype(i.radius()) { return i.radius(); }
template<class UB> inline auto width(Interval<UB> const& i) -> decltype(i.width()) { return i.width(); }


//! \related Interval \brief Test if the interval is empty.
template<class UB> inline auto is_empty(Interval<UB> const& i) -> decltype(i.lower()>i.upper()) { return i.lower()>i.upper(); }
//! \related Interval \brief Test if the interval is empty.
template<class UB> inline auto is_singleton(Interval<UB> const& i) -> decltype(i.lower()==i.upper()) { return i.lower()==i.upper(); }
//! \related Interval \brief Test if the interval is bounded.
template<class UB> inline auto is_bounded(Interval<UB> const& i) -> decltype(i.upper()<inf) { return -i.lower()<inf && i.upper()<inf; }

//! \related Interval \brief Test if \a x1 is an element of the interval \a i2.
template<class UB, class X> inline auto element(X const& x1, Interval<UB> const& i2) -> decltype(i2.lower()<=x1 && i2.upper()>=x1) {
    return i2.lower()<=x1 && i2.upper()>=x1; }
//! \related Interval \brief Test if the interval \a i1 contains \a x2.
template<class UB, class X> inline auto contains(Interval<UB> const& i1, X const& x2) -> decltype(i1.lower()<=x2 && i1.upper()>=x2) {
    return i1.lower()<=x2 && i1.upper()>=x2; }
//! \related Interval \brief Test if the interval \a i1 is equal to \a i2.
template<class UB> inline auto equal(Interval<UB> const& i1, Interval<UB> const& i2) -> decltype(i1.upper()==i2.upper()) {
    return i1.lower()<=i2.lower() && i1.upper()==i2.upper(); }
//! \related Interval \brief Test if the interval \a i1 is a subset of \a i2.
template<class UB1, class UB2> inline auto subset(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<=i2.upper()) {
    return i1.lower()>=i2.lower() && i1.upper()<=i2.upper(); }
//! \related Interval \brief Test if the interval \a i1 is a superset of \a i2.
template<class UB1, class UB2> inline auto superset(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>=i2.upper()) {
    return i1.lower()<=i2.lower() && i1.upper()>=i2.upper(); }
//! \related Interval \brief Test if the interval \a i1 is disjoint from \a i2. Returns \c false even if the two intervals only have an endpoint in common.
template<class UB1, class UB2> inline auto disjoint(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<i2.lower()) {
    return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }
//! \related Interval \brief Test if the interval \a i1 intersects \a i2. Returns \c true even if the two intervals only have an endpoint in common.
template<class UB1, class UB2> inline auto intersect(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>=i2.lower()) {
    return i1.lower()<=i2.upper() && i1.upper()>=i2.lower(); }

//! \related Interval \brief Test if the closed interval \a i1 is disjoint from the closed interval \a i2.
//! Returns \c false if the two intervals only have an endpoint in common.
template<class UB1, class UB2> inline auto separated(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<i2.lower()) {
    return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }
//! \related Interval \brief Test if the interval \a i1 overlaps \a i2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
template<class UB1, class UB2> inline auto overlap(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>i2.lower()) {
    return i1.lower()<i2.upper() && i1.upper()>i2.lower(); }
//! \related Interval \brief Test if the (closed) interval \a i1 is a subset of the interior of \a i2.
template<class UB1, class UB2> inline auto inside(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<i2.upper()) {
    return i1.lower()>i2.lower() && i1.upper()<i2.upper(); }
//! \related Interval \brief Test if the interior of the interval \a i1 is a superset of the (closed) interval \a i2.
template<class UB1, class UB2> inline auto covers(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>i2.upper()) {
    return i1.lower()<i2.lower() && i1.upper()>i2.upper(); }

//! \related Interval \brief The intersection of two intervals.
template<class UB1, class UB2> inline Interval<decltype(max(declval<UB1>(),declval<UB2>()))> intersection(Interval<UB1> const& i1, Interval<UB2> const& i2) {
    return Interval<decltype(max(declval<UB1>(),declval<UB2>()))>(max(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}

//! \related Interval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
template<class UB1, class UB2> inline auto
hull(Interval<UB1> const& i1, Interval<UB2> const& i2) -> Interval<decltype(max(declval<UB1>(),declval<UB2>()))> {
    typedef decltype(max(declval<UB1>(),declval<UB2>())) UB0; return Interval<UB0>(min(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
template<class UB, class X> inline auto
hull(Interval<UB> const& i1, X x2) -> Interval<decltype(max(declval<UB>(),declval<X>()))> {
    return Interval<decltype(max(declval<UB>(),declval<X>()))>(min(i1.lower(),x2),max(i1.upper(),x2));
}

//! \related Interval \brief The hull of an interval and a point, equal to the smallest interval containing both.
template<class UB, class X> inline auto
hull(X x1, Interval<UB> const& i2) -> Interval<decltype(max(declval<UB>(),declval<X>()))> {
    return Interval<decltype(max(declval<UB>(),declval<X>()))>(min(x1,i2.lower()),max(x1,i2.upper()));
}


//! \related Interval \brief Split an interval into its lower, middle or upper half.
template<class UB> inline auto split(Interval<UB> const& ivl, SplitPart lmu) -> Interval<UB> {
    auto cc=(ivl.lower()+ivl.upper())/2;
    if(lmu==SplitPart::LOWER) {
        return Interval<UB>(ivl.lower(),make_split_point(cc));
    } else if(lmu==SplitPart::UPPER) {
        return Interval<UB>(make_split_point(cc),ivl.upper());
    } else {
        auto cl=(3*ivl.lower()+ivl.upper())/4;
        auto cu=(ivl.lower()+3*ivl.upper())/4;
        return Interval<UB>(make_split_point(cl),make_split_point(cu));
    }
}

inline BoundFloat64 make_singleton(Interval<UpperFloat64> const& i) {
    return BoundFloat64(i.lower(),i.upper());
}

inline Interval<UpperFloat64> make_interval(BoundFloat64 const& x) {
    return Interval<UpperFloat64>(x.lower(),x.upper());
}

inline Interval<UpperFloat64> widen(Interval<UpperFloat64> const& x, UpperFloat64 e) {
    return Interval<UpperFloat64>(x.lower()-e,x.upper()+e); }

inline Interval<UpperFloat64> widen(Interval<UpperFloat64> const& x) {
    return widen(x,UpperFloat64(Float64::min())); }

inline Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& x, UpperFloat64 e) {
    return Interval<LowerFloat64>(x.lower()+e,x.upper()-e); }

inline Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& x) {
    return narrow(x,UpperFloat64(Float64::min())); }

inline bool operator==(Interval<ExactFloat64> const& i1, Interval<ExactFloat64> const& i2) {
    return equal(i1,i2);
}

inline bool operator!=(Interval<ExactFloat64> const& i1, Interval<ExactFloat64> const& i2) {
    return !(i1==i2);
}

//! \related Interval \brief Test if the interval \a i1 is equal to \a i2.
template<class UB> inline auto operator==(Interval<UB> const& i1, Interval<UB> const& i2) -> decltype(i1.upper()==i2.upper()) {
    return i1.lower()<=i2.lower() && i1.upper()==i2.upper(); }


inline ExactInterval make_exact(ApproximateInterval const& ivl) {
    return reinterpret_cast<ExactInterval const&>(ivl);
}

inline ExactInterval make_exact_interval(ApproximateInterval const& ivl) {
    return reinterpret_cast<ExactInterval const&>(ivl);
}
/*


class RealInterval : public Interval<Real> { public: using Interval<Real>::Interval; };
class Float64Interval : public Interval<Flt64> { public: using Interval<Flt64>::Interval; };
class ExactFloat64Interval : public Interval<ExactFloat64> { public: using Interval<ExactFloat64>::Interval; };
class LowerFloat64Interval : public Interval<LowerFloat64> { public: using Interval<LowerFloat64>::Interval; };
class ApproximateFloat64Interval : public Interval<ApproximateFloat64> { public: using Interval<ApproximateFloat64>::Interval; };

class UpperFloat64Interval : public Interval<UpperFloat64> { public: typedef UpperFloat64Interval NumericType; using Interval<UpperFloat64>::Interval; };
template<> struct IsFloat64<UpperFloat64Interval> : True { };

//! \related Interval \brief Read from an output stream
InputStream& operator>>(InputStream& os, UpperFloat64Interval& ivl);


// Standard equality operators
//! \related UpperFloat64Interval \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
inline bool operator==(const ExactFloat64Interval& i1, const ExactFloat64Interval& i2) { return i1.lower()==i2.lower() && i1.upper()==i2.upper(); }
//! \related UpperFloat64Interval \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true.
inline bool operator!=(const ExactFloat64Interval& i1, const ExactFloat64Interval& i2) { return i1.lower()!=i2.lower() || i1.upper()!=i2.upper(); }


//! \related UpperFloat64Interval \brief Truncate interval to single-precision.
UpperFloat64Interval trunc(UpperFloat64Interval x);
UpperFloat64Interval trunc(UpperFloat64Interval x, uint n);

//! \related UpperFloat64Interval \brief Widen interval by 1ulp.
UpperFloat64Interval widen(UpperFloat64Interval x);
//! \related UpperFloat64Interval \brief Narrow interval by 1ulp.
LowerFloat64Interval narrow(LowerFloat64Interval x);

//! \related UpperFloat64Interval \brief An interval containing the given interval in its interior.
UpperFloat64Interval widen(UpperFloat64Interval i);
//! \related LowerFloat64Interval \brief An interval contained in the interior of the given interval.
LowerFloat64Interval narrow(LowerFloat64Interval i);


//! \related Interval<Real> \brief An under-approximation to an interval.
inline LowerFloat64Interval under_approximation(const RealInterval& rivl) {
    return LowerFloat64Interval(rivl.lower().upper(),rivl.upper().lower());
}
//! \related Interval<Real> \brief An over-approximation to an interval.
inline UpperFloat64Interval over_approximation(const RealInterval& rivl) {
    return ExactFloat64Interval(make_exact(rivl.lower().lower()),make_exact(rivl.upper().upper()));
}
//! \related Interval<Real> \brief An approximation to an interval.
inline ApproximateFloat64Interval approximation(const RealInterval& rivl) {
    return ApproximateFloat64Interval(rivl.lower().approx(),rivl.upper().approx());
}

//! \related Interval<Real> \brief Convert an interval to the Exact paradigm.
inline ExactFloat64Interval const& make_exact(const ApproximateFloat64Interval& ivl) {
    return reinterpret_cast<ExactFloat64Interval const&>(ivl);
}



// Validated paradigm checks

//! \related UpperFloat64Interval \brief Computes a common refinement of \a i1 and \a i2 i.e. the intersection.
inline UpperFloat64Interval refinement(UpperFloat64Interval const& i1, UpperFloat64Interval const& i2) {
    return UpperFloat64Interval(max(i1.lower(),i2.lower()),min(i1.upper(),i2.upper())); }

//! \related UpperFloat64Interval \brief Tests if \a i1 provides a better over-approximation to the exact interval than \a i2.
//! i.e. \a i1 is a subset of \a i2.
inline bool refines(UpperFloat64Interval const& i1, UpperFloat64Interval const& i2) {
    return make_raw(i1.lower())>=make_raw(i2.lower()) && make_raw(i1.upper())<=make_raw(i2.upper()); }

//! \related UpperFloat64Interval \brief Computes if two over-approximating intervals have the same representation.
inline bool same(UpperFloat64Interval const& i1, UpperFloat64Interval const& i2) {
    return make_raw(i1.lower())==make_raw(i2.lower()) && make_raw(i1.upper())==make_raw(i2.upper()); }

//! \related UpperFloat64Interval \brief Restricts \a i to to interval \a r. i.e. Performs an intersection in-place.
inline Void restrict(UpperFloat64Interval& i, UpperFloat64Interval const& r) {
    i.set_lower(max(i.lower(),r.lower())); i.set_upper(min(i.upper(),r.upper()));
}

*/

} // namespace Ariadne
