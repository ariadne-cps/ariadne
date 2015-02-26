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

template<class U> inline OutputStream& operator<<(OutputStream& os, Interval<U> const& ivl) {
    return os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
}


template<class U> inline auto lower_bound(Interval<U> const& i) -> decltype(i.lower()) { return i.lower(); }
template<class U> inline auto upper_bound(Interval<U> const& i) -> decltype(i.upper()) { return i.upper(); }
template<class UB> inline auto centre(Interval<UB> const& i) -> decltype(i.centre()) { return i.centre(); }
template<class UB> inline auto midpoint(Interval<UB> const& i) -> decltype(i.midpoint()) { return i.midpoint(); }
template<class UB> inline auto radius(Interval<UB> const& i) -> decltype(i.radius()) { return i.radius(); }
template<class UB> inline auto width(Interval<UB> const& i) -> decltype(i.width()) { return i.width(); }

template<class UB> inline auto is_empty(Interval<UB> const& i) -> decltype(i.lower()>i.upper()) { return i.lower()>i.upper(); }
template<class UB> inline auto is_singleton(Interval<UB> const& i) -> decltype(i.lower()==i.upper()) { return i.lower()==i.upper(); }
template<class UB> inline auto is_bounded(Interval<UB> const& i) -> decltype(i.upper()<inf) { return -i.lower()<inf && i.upper()<inf; }


template<class UB, class X> inline auto element(X const& x1, Interval<UB> const& i2) -> decltype(i2.lower()<=x1 && i2.upper()>=x1) {
    return i2.lower()<=x1 && i2.upper()>=x1; }

template<class UB, class X> inline auto contains(Interval<UB> const& i1, X const& x2) -> decltype(i1.lower()<=x2 && i1.upper()>=x2) {
    return i1.lower()<=x2 && i1.upper()>=x2; }

template<class UB> inline auto equal(Interval<UB> const& i1, Interval<UB> const& i2) -> decltype(i1.upper()==i2.upper()) {
    return i1.lower()<=i2.lower() && i1.upper()==i2.upper(); }

template<class UB1, class UB2> inline auto subset(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<=i2.upper()) {
    return i1.lower()>=i2.lower() && i1.upper()<=i2.upper(); }

template<class UB1, class UB2> inline auto superset(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>=i2.upper()) {
    return i1.lower()<=i2.lower() && i1.upper()>=i2.upper(); }

template<class UB1, class UB2> inline auto disjoint(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<i2.lower()) {
    return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }

template<class UB1, class UB2> inline auto intersect(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>=i2.lower()) {
    return i1.lower()<=i2.upper() && i1.upper()>=i2.lower(); }


template<class UB1, class UB2> inline auto separated(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<i2.lower()) {
    return i1.lower()>i2.upper() || i1.upper()<i2.lower(); }

template<class UB1, class UB2> inline auto overlap(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>i2.lower()) {
    return i1.lower()<i2.upper() && i1.upper()>i2.lower(); }

template<class UB1, class UB2> inline auto inside(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()<i2.upper()) {
    return i1.lower()>i2.lower() && i1.upper()<i2.upper(); }

template<class UB1, class UB2> inline auto covers(Interval<UB1> const& i1, Interval<UB2> const& i2) -> decltype(i1.upper()>i2.upper()) {
    return i1.lower()<i2.lower() && i1.upper()>i2.upper(); }


template<class UB1, class UB2> inline auto
intersection(Interval<UB1> const& i1, Interval<UB2> const& i2) -> Interval<decltype(min(declval<UB1>(),declval<UB2>()))> {
    return Interval<decltype(min(declval<UB1>(),declval<UB2>()))>(max(i1.lower(),i2.lower()),min(i1.upper(),i2.upper()));
}


template<class UB1, class UB2> inline auto
hull(Interval<UB1> const& i1, Interval<UB2> const& i2) -> Interval<decltype(max(declval<UB1>(),declval<UB2>()))> {
    typedef decltype(max(declval<UB1>(),declval<UB2>())) UB0; return Interval<UB0>(min(i1.lower(),i2.lower()),max(i1.upper(),i2.upper()));
}

template<class UB, class X> inline auto
hull(Interval<UB> const& i1, X x2) -> Interval<decltype(max(declval<UB>(),declval<X>()))> {
    return Interval<decltype(max(declval<UB>(),declval<X>()))>(min(i1.lower(),x2),max(i1.upper(),x2));
}

template<class UB, class X> inline auto
hull(X x1, Interval<UB> const& i2) -> Interval<decltype(max(declval<UB>(),declval<X>()))> {
    return Interval<decltype(max(declval<UB>(),declval<X>()))>(min(x1,i2.lower()),max(x1,i2.upper()));
}



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



template<class UB> inline auto operator==(Interval<UB> const& i1, Interval<UB> const& i2) -> decltype(i1.upper()==i2.upper()) {
    return i1.lower()==i2.lower() && i1.upper()==i2.upper(); }
template<class UB> inline auto operator!=(Interval<UB> const& i1, Interval<UB> const& i2) -> decltype(i1.upper()!=i2.upper()) {
    return i1.lower()!=i2.lower() || i1.upper()!=i2.upper(); }


inline Interval<UpperFloat64> refinement(Interval<UpperFloat64> const& ivl1, Interval<UpperFloat64> const& ivl2) {
    return Interval<UpperFloat64>(max(ivl1.lower().raw(),ivl2.lower().raw()),min(ivl1.upper().raw(),ivl2.upper().raw())); }
inline Bool refines(UpperInterval const& ivl1, UpperInterval const& ivl2) {
    return ivl1.lower().raw()>=ivl2.lower().raw() && ivl1.upper().raw()<=ivl2.upper().raw(); }
inline Bool same(UpperInterval const& ivl1, UpperInterval const& ivl2) {
    return ivl1.lower().raw()==ivl2.lower().raw() && ivl1.upper().raw()==ivl2.upper().raw(); }

inline Interval<UpperFloat64> widen(Interval<UpperFloat64> const& ivl, UpperFloat64 e) {
    return Interval<UpperFloat64>(ivl.lower()-e,ivl.upper()+e); }
inline Interval<UpperFloat64> widen(Interval<UpperFloat64> const& ivl) {
    return widen(ivl,UpperFloat64(Float64::min())); }

inline Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& ivl, UpperFloat64 e) {
    return Interval<LowerFloat64>(ivl.lower()+e,ivl.upper()-e); }
inline Interval<LowerFloat64> narrow(Interval<LowerFloat64> const& ivl) {
    return narrow(ivl,UpperFloat64(Float64::min())); }

inline ExactInterval make_exact(ApproximateInterval const& ivl) {
    return reinterpret_cast<ExactInterval const&>(ivl); }
inline ExactInterval make_exact_interval(ApproximateInterval const& ivl) {
    return reinterpret_cast<ExactInterval const&>(ivl); }

} // namespace Ariadne
