/***************************************************************************
 *            geometry/interval.inl.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

namespace Ariadne {

template<class M> inline M make_split_point(M const& m) { return m; }
template<class F> inline Value<F> make_split_point(Approximation<F> const& am) { return reinterpret_cast<Value<F>const&>(am); }
template<class F> inline Value<F> make_split_point(Bounds<F> const& bm) { return bm.value(); }
template<class F> inline Value<F> make_split_point(Ball<F> const& bm) { return bm.value(); }


template<class U> Interval<U>::Interval() : Interval(EmptyInterval()) { }
template<class U> Interval<U>::Interval(EmptyInterval const&) {
    if constexpr(IsConstructibleGivenDefaultPrecision<U,Dyadic>::value) {
        _l = L(Dyadic::inf(Sign::POSITIVE),L::RawType::get_default_precision());
        _u = U(Dyadic::inf(Sign::NEGATIVE),U::RawType::get_default_precision());
    } else {
        _l = Dyadic::inf(Sign::POSITIVE);
        _u = Dyadic::inf(Sign::NEGATIVE);
    }
}
template<class U> Interval<U>::Interval(UnitInterval const&) : Interval(-1,+1) { }
template<class U> Interval<U>::Interval(EntireInterval const&) {
    if constexpr(IsConstructibleGivenDefaultPrecision<U,Dyadic>::value) {
        _l = L(Dyadic::inf(Sign::NEGATIVE),L::RawType::get_default_precision());
        _u = U(Dyadic::inf(Sign::POSITIVE),U::RawType::get_default_precision());
    } else {
        _l = Dyadic::inf(Sign::NEGATIVE);
        _u = Dyadic::inf(Sign::POSITIVE);
    }
}
template<class U> Interval<U>::Interval(LowerBoundType l, UpperBoundType u) : _l(l), _u(u) { }

template<class U> Interval<U> Interval<U>::create_zero() const { return Interval<U>(0,0); }
template<class U> SizeOne Interval<U>::dimension() const  { return SizeOne(); }
template<class U> auto Interval<U>::lower() const -> LowerBoundType const& { return _l; }
template<class U> auto Interval<U>::upper() const -> UpperBoundType const& { return _u; }
template<class U> auto Interval<U>::centre() const -> CentreType { return hlf(_l+_u); }
template<class U> auto Interval<U>::midpoint() const -> MidpointType { auto m=(hlf(_l+_u)); return make_split_point(m); }
template<class U> auto Interval<U>::radius() const -> RadiusType { return cast_positive(max(this->upper()-this->midpoint(),this->midpoint()-this->lower())); }
template<class U> auto Interval<U>::width() const -> WidthType { return cast_positive(this->upper()-this->lower()); }

template<class U> Interval<U> Interval<U>::empty_interval() { return Interval<U>(EmptyInterval()); }
template<class U> Interval<U> Interval<U>::unit_interval() { return Interval<U>(-1,+1); }
template<class U> Interval<U> Interval<U>::biinfinite_interval() { return Interval<U>(EntireInterval()); }

template<class U> auto Interval<U>::is_empty() const -> decltype(declval<L>()>declval<U>()) { return this->_l > this->_u; }
template<class U> auto Interval<U>::is_bounded() const -> decltype(declval<U>()<declval<L>()) { return Ariadne::is_bounded(*this); }
template<class U> auto Interval<U>::is_singleton() const -> decltype(declval<L>() == declval<U>()) { return this->_l == this->_u; }

template<class U> auto Interval<U>::set_lower(LowerBoundType l) -> void { _l=l; }
template<class U> auto Interval<U>::set_upper(UpperBoundType u) -> void { _u=u; }
template<class U> auto Interval<U>::set(LowerBoundType l, UpperBoundType u) -> void { _l=l; _u=u; }

template<class U> inline OutputStream& operator<<(OutputStream& os, Interval<U> const& ivl) {
    return os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
}

template<> OutputStream& operator<<(OutputStream& os, Interval<FloatDPValue> const& ivl);

template<class U> inline auto lower_bound(Interval<U> const& ivl) -> decltype(ivl.lower()) { return ivl.lower(); }
template<class U> inline auto upper_bound(Interval<U> const& ivl) -> decltype(ivl.upper()) { return ivl.upper(); }
template<class U> inline auto centre(Interval<U> const& ivl) -> decltype(ivl.centre()) { return ivl.centre(); }
template<class U> inline auto midpoint(Interval<U> const& ivl) -> decltype(ivl.midpoint()) { return ivl.midpoint(); }
template<class U> inline auto radius(Interval<U> const& ivl) -> decltype(ivl.radius()) { return ivl.radius(); }
template<class U> inline auto width(Interval<U> const& ivl) -> decltype(ivl.width()) { return ivl.width(); }

template<class U> inline auto is_empty(Interval<U> const& ivl) -> decltype(ivl.lower()>ivl.upper()) { return ivl.lower()>ivl.upper(); }
template<class U> inline auto is_singleton(Interval<U> const& ivl) -> decltype(ivl.lower()==ivl.upper()) { return ivl.lower()==ivl.upper(); }

template<class U> inline auto is_bounded(Interval<U> const& ivl) -> decltype(ivl.upper()<ivl.lower()) { return -ivl.lower().raw()<inf && ivl.upper().raw()<inf; }
template<> inline auto is_bounded(Interval<Real> const& ivl) -> decltype(ivl.upper()<ivl.lower()) {return -ivl.lower()<infinity && ivl.upper()<infinity; }


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


template<class U1, class U2> inline decltype(auto) refinement(Interval<U1> const& ivl1, Interval<U2> const& ivl2) {
    return make_interval(refinement(ivl1.lower(),ivl2.lower()),refinement(ivl1.upper(),ivl2.upper())); }
template<class U1, class U2> inline bool refines(Interval<U1> const& ivl1, Interval<U2> const& ivl2) {
    return refines(ivl1.lower(),ivl2.lower()) and refines(ivl1.upper(),ivl2.upper()); }
template<class U1, class U2> inline bool same(Interval<U1> const& ivl1, Interval<U2> const& ivl2) {
    return same(ivl1.lower(),ivl2.lower()) and same(ivl1.upper(),ivl2.upper()); }

template<class F> inline Interval<UpperBound<F>> refinement(Interval<UpperBound<F>> const& ivl1, Interval<UpperBound<SelfType<F>>> const& ivl2) {
    return Interval<UpperBound<F>>(max(ivl1.lower().raw(),ivl2.lower().raw()),min(ivl1.upper().raw(),ivl2.upper().raw())); }
template<class F> inline Bool refines(Interval<UpperBound<F>> const& ivl1, Interval<UpperBound<SelfType<F>>> const& ivl2) {
    return ivl1.lower().raw()>=ivl2.lower().raw() && ivl1.upper().raw()<=ivl2.upper().raw(); }
template<class F> inline Bool same(Interval<UpperBound<F>> const& ivl1, Interval<UpperBound<SelfType<F>>> const& ivl2) {
    return ivl1.lower().raw()==ivl2.lower().raw() && ivl1.upper().raw()==ivl2.upper().raw(); }

template<class F> inline Interval<UpperBound<F>> widen(Interval<UpperBound<F>> const& ivl, UpperBound<F> e) {
    return Interval<UpperBound<F>>(ivl.lower()-e,ivl.upper()+e); }
template<class F> inline Interval<UpperBound<F>> widen(Interval<UpperBound<F>> const& ivl) {
    return widen(ivl,UpperBound<F>(F::min(ivl.upper().precision()))); }
template<class F> inline Interval<UpperBound<F>> widen(Interval<Value<F>> const& ivl) {
    return widen(Interval<UpperBound<F>>(ivl),UpperBound<F>(F::min(ivl.upper().precision()))); }

template<class F> inline Interval<LowerBound<F>> narrow(Interval<LowerBound<F>> const& ivl, UpperBound<F> e) {
    return Interval<LowerBound<F>>(ivl.lower()+e,ivl.upper()-e); }
template<class F> inline Interval<LowerBound<F>> narrow(Interval<LowerBound<F>> const& ivl) {
    return narrow(ivl,UpperBound<F>(F::min(ivl.upper().precision()))); }
template<class F> inline Interval<LowerBound<F>> narrow(Interval<Value<F>> const& ivl) {
    return narrow(Interval<LowerBound<F>>(ivl),UpperBound<F>(F::min(ivl.upper().precision()))); }

inline Interval<FloatValue<DP>> cast_exact(Interval<FloatApproximation<DP>> const& ivl) {
    return reinterpret_cast<Interval<FloatValue<DP>> const&>(ivl); }
inline Interval<FloatValue<MP>> cast_exact(Interval<FloatApproximation<MP>> const& ivl) {
    return reinterpret_cast<Interval<FloatValue<MP>> const&>(ivl); }
inline Interval<FloatValue<DP>> cast_exact(Interval<FloatUpperBound<DP>> const& ivl) {
    return reinterpret_cast<Interval<FloatValue<DP>> const&>(ivl); }
inline Interval<FloatValue<MP>> cast_exact(Interval<FloatUpperBound<MP>> const& ivl) {
    return reinterpret_cast<Interval<FloatValue<MP>> const&>(ivl); }
inline Interval<FloatValue<DP>> cast_exact_interval(Interval<FloatApproximation<DP>> const& ivl) {
    return reinterpret_cast<Interval<FloatValue<DP>> const&>(ivl); }
inline Interval<FloatValue<MP>> cast_exact_interval(Interval<FloatApproximation<MP>> const& ivl) {
    return reinterpret_cast<Interval<FloatValue<MP>> const&>(ivl); }

inline FloatDPLowerBound mig(FloatDPUpperInterval const& ivl) { return mig(cast_singleton(ivl)); }
inline FloatMPLowerBound mig(FloatMPUpperInterval const& ivl) { return mig(cast_singleton(ivl)); }

} // namespace Ariadne
