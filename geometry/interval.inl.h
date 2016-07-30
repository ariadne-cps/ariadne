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


namespace Ariadne {

template<class M> inline M make_split_point(M const& m) { return m; }
template<class PR> inline FloatValue<PR> make_split_point(FloatApproximation<PR> const& am) { return cast_exact(am); }
template<class PR> inline FloatValue<PR> make_split_point(FloatBounds<PR> const& vm) { return vm.value(); }
template<class PR> inline FloatValue<PR> make_split_point(FloatBall<PR> const& bm) { return bm.value(); }


template<class U> Interval<U>::Interval() : _l(+infty), _u(-infty) { }
template<class U> Interval<U>::Interval(EmptyIntervalType) : _l(+infty), _u(-infty) { }
template<class U> Interval<U>::Interval(LowerBoundType l, UpperBoundType u) : _l(l), _u(u) { }

template<class U> Interval<U> Interval<U>::create_zero() const { return Interval<U>(0,0); }
template<class U> SizeOne Interval<U>::dimension() const  { return SizeOne(); }
template<class U> auto Interval<U>::lower() const -> LowerBoundType const& { return _l; }
template<class U> auto Interval<U>::upper() const -> UpperBoundType const& { return _u; }
template<class U> auto Interval<U>::centre() const -> CentreType { return (_l+_u)/2u; }
template<class U> auto Interval<U>::midpoint() const -> MidpointType { auto m=((_l+_u)/2); return make_split_point(m); }
template<class U> auto Interval<U>::radius() const -> RadiusType { return cast_positive(max(this->upper()-this->midpoint(),this->midpoint()-this->lower())); }
template<class U> auto Interval<U>::width() const -> WidthType { return cast_positive(this->upper()-this->lower()); }

template<class U> Interval<U> Interval<U>::empty_interval() { return Interval<U>(+infty,-infty); }
template<class U> Interval<U> Interval<U>::unit_interval() { return Interval<U>(-1,+1); }

template<class U> auto Interval<U>::is_empty() const -> decltype(declval<L>()>declval<U>()) { return this->_l > this->_u; }
template<class U> auto Interval<U>::is_bounded() const -> decltype(declval<U>()<infty) { return this->_l > -infty && this->_u < +infty; }
template<class U> auto Interval<U>::is_singleton() const -> decltype(declval<L>() == declval<U>()) { return this->_l == this->_u; }

template<class U> auto Interval<U>::set_lower(LowerBoundType l) -> void { _l=l; }
template<class U> auto Interval<U>::set_upper(UpperBoundType u) -> void { _u=u; }
template<class U> auto Interval<U>::set(LowerBoundType l, UpperBoundType u) -> void { _l=l; _u=u; }

template<class U> inline OutputStream& operator<<(OutputStream& os, Interval<U> const& ivl) {
    return os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
}


template<class U> inline auto lower_bound(Interval<U> const& ivl) -> decltype(ivl.lower()) { return ivl.lower(); }
template<class U> inline auto upper_bound(Interval<U> const& ivl) -> decltype(ivl.upper()) { return ivl.upper(); }
template<class U> inline auto centre(Interval<U> const& ivl) -> decltype(ivl.centre()) { return ivl.centre(); }
template<class U> inline auto midpoint(Interval<U> const& ivl) -> decltype(ivl.midpoint()) { return ivl.midpoint(); }
template<class U> inline auto radius(Interval<U> const& ivl) -> decltype(ivl.radius()) { return ivl.radius(); }
template<class U> inline auto width(Interval<U> const& ivl) -> decltype(ivl.width()) { return ivl.width(); }

template<class U> inline auto is_empty(Interval<U> const& ivl) -> decltype(ivl.lower()>ivl.upper()) { return ivl.lower()>ivl.upper(); }
template<class U> inline auto is_singleton(Interval<U> const& ivl) -> decltype(ivl.lower()==ivl.upper()) { return ivl.lower()==ivl.upper(); }
template<class U> inline auto is_bounded(Interval<U> const& ivl) -> decltype(ivl.upper()<inf) { return -ivl.lower()<inf && ivl.upper()<inf; }


template<class U, class X> inline auto element(X const& x1, Interval<U> const& ivl2) -> decltype(ivl2.lower()<=x1 && ivl2.upper()>=x1) {
    return ivl2.lower()<=x1 && ivl2.upper()>=x1; }

template<class U, class X> inline auto contains(Interval<U> const& ivl1, X const& x2) -> decltype(ivl1.lower()<=x2 && ivl1.upper()>=x2) {
    return ivl1.lower()<=x2 && ivl1.upper()>=x2; }

template<class U> inline auto equal(Interval<U> const& ivl1, Interval<U> const& ivl2) -> decltype(ivl1.upper()==ivl2.upper()) {
    return ivl1.lower()<=ivl2.lower() && ivl1.upper()==ivl2.upper(); }

template<class UB1, class UB2> inline auto subset(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()<=ivl2.upper()) {
    return ivl1.lower()>=ivl2.lower() && ivl1.upper()<=ivl2.upper(); }

template<class UB1, class UB2> inline auto superset(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()>=ivl2.upper()) {
    return ivl1.lower()<=ivl2.lower() && ivl1.upper()>=ivl2.upper(); }

template<class UB1, class UB2> inline auto disjoint(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()<ivl2.lower()) {
    return ivl1.lower()>ivl2.upper() || ivl1.upper()<ivl2.lower(); }

template<class UB1, class UB2> inline auto intersect(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()>=ivl2.lower()) {
    return ivl1.lower()<=ivl2.upper() && ivl1.upper()>=ivl2.lower(); }


template<class UB1, class UB2> inline auto separated(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()<ivl2.lower()) {
    return ivl1.lower()>ivl2.upper() || ivl1.upper()<ivl2.lower(); }

template<class UB1, class UB2> inline auto overlap(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()>ivl2.lower()) {
    return ivl1.lower()<ivl2.upper() && ivl1.upper()>ivl2.lower(); }

template<class UB1, class UB2> inline auto inside(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()<ivl2.upper()) {
    return ivl1.lower()>ivl2.lower() && ivl1.upper()<ivl2.upper(); }

template<class UB1, class UB2> inline auto covers(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> decltype(ivl1.upper()>ivl2.upper()) {
    return ivl1.lower()<ivl2.lower() && ivl1.upper()>ivl2.upper(); }


template<class UB1, class UB2> inline auto
intersection(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> Interval<decltype(min(declval<UB1>(),declval<UB2>()))> {
    return Interval<decltype(min(declval<UB1>(),declval<UB2>()))>(max(ivl1.lower(),ivl2.lower()),min(ivl1.upper(),ivl2.upper()));
}


template<class UB1, class UB2> inline auto
hull(Interval<UB1> const& ivl1, Interval<UB2> const& ivl2) -> Interval<decltype(max(declval<UB1>(),declval<UB2>()))> {
    typedef decltype(max(declval<UB1>(),declval<UB2>())) UB0; return Interval<UB0>(min(ivl1.lower(),ivl2.lower()),max(ivl1.upper(),ivl2.upper()));
}

template<class U, class X> inline auto
hull(Interval<U> const& ivl1, X x2) -> Interval<decltype(max(declval<U>(),declval<X>()))> {
    return Interval<decltype(max(declval<U>(),declval<X>()))>(min(ivl1.lower(),x2),max(ivl1.upper(),x2));
}

template<class U, class X> inline auto
hull(X x1, Interval<U> const& ivl2) -> Interval<decltype(max(declval<U>(),declval<X>()))> {
    return Interval<decltype(max(declval<U>(),declval<X>()))>(min(x1,ivl2.lower()),max(x1,ivl2.upper()));
}



template<class U> inline auto split(Interval<U> const& ivl, SplitPart lmu) -> Interval<U> {
    auto cc=(ivl.lower()+ivl.upper())/2;
    if(lmu==SplitPart::LOWER) {
        return Interval<U>(ivl.lower(),make_split_point(cc));
    } else if(lmu==SplitPart::UPPER) {
        return Interval<U>(make_split_point(cc),ivl.upper());
    } else {
        auto cl=(3*ivl.lower()+ivl.upper())/4;
        auto cu=(ivl.lower()+3*ivl.upper())/4;
        return Interval<U>(make_split_point(cl),make_split_point(cu));
    }
}



template<class U> inline auto operator==(Interval<U> const& ivl1, Interval<U> const& ivl2) -> decltype(ivl1.upper()==ivl2.upper()) {
    return ivl1.lower()==ivl2.lower() && ivl1.upper()==ivl2.upper(); }
template<class U> inline auto operator!=(Interval<U> const& ivl1, Interval<U> const& ivl2) -> decltype(ivl1.upper()!=ivl2.upper()) {
    return ivl1.lower()!=ivl2.lower() || ivl1.upper()!=ivl2.upper(); }


inline Interval<Float64UpperBound> refinement(Interval<Float64UpperBound> const& ivl1, Interval<Float64UpperBound> const& ivl2) {
    return Interval<Float64UpperBound>(max(ivl1.lower().raw(),ivl2.lower().raw()),min(ivl1.upper().raw(),ivl2.upper().raw())); }
inline Bool refines(Interval<Float64UpperBound> const& ivl1, Interval<Float64UpperBound> const& ivl2) {
    return ivl1.lower().raw()>=ivl2.lower().raw() && ivl1.upper().raw()<=ivl2.upper().raw(); }
inline Bool same(Interval<Float64UpperBound> const& ivl1, Interval<Float64UpperBound> const& ivl2) {
    return ivl1.lower().raw()==ivl2.lower().raw() && ivl1.upper().raw()==ivl2.upper().raw(); }

inline Interval<Float64UpperBound> widen(Interval<Float64UpperBound> const& ivl, Float64UpperBound e) {
    return Interval<Float64UpperBound>(ivl.lower()-e,ivl.upper()+e); }
inline Interval<Float64UpperBound> widen(Interval<Float64UpperBound> const& ivl) {
    return widen(ivl,Float64UpperBound(Float64::min(Precision64()))); }

inline Interval<Float64LowerBound> narrow(Interval<Float64LowerBound> const& ivl, Float64UpperBound e) {
    return Interval<Float64LowerBound>(ivl.lower()+e,ivl.upper()-e); }
inline Interval<Float64LowerBound> narrow(Interval<Float64LowerBound> const& ivl) {
    return narrow(ivl,Float64UpperBound(Float64::min(Precision64()))); }

inline Interval<Float64Value> cast_exact(Interval<Float64Approximation> const& ivl) {
    return reinterpret_cast<ExactIntervalType const&>(ivl); }
inline Interval<Float64Value> cast_exact_interval(Interval<Float64Approximation> const& ivl) {
    return reinterpret_cast<Interval<Float64Value> const&>(ivl); }

} // namespace Ariadne
