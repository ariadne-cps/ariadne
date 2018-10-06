/***************************************************************************
 *            interval.cpp
 *
 *  Copyright 2013--17  Pieter Collins
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

#include "interval.hpp"
#include "../numeric/dyadic.hpp"

namespace Ariadne {

static const uint GEOMETRY_OUTPUT_PLACES = 4;

template<> OutputStream& operator<<(OutputStream& os, Interval<FloatDPValue> const& ivl) {
    auto places = FloatDPValue::output_places;
    FloatDPValue::output_places=GEOMETRY_OUTPUT_PLACES;
    os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
    FloatDPValue::output_places=places;
    return os;
}

Interval<FloatDPValue> widen_domain(Interval<FloatDPUpperBound> const& ivl) {
    auto rnd=FloatDP::get_rounding_mode();
    FloatDP::set_rounding_mode(upward);
    volatile float min=std::numeric_limits<float>::min();
    volatile double neg_l=(-ivl.lower()).get_d();
    volatile double l=-neg_l;
    volatile double u=ivl.upper().get_d();
    volatile float neg_rl=neg_l;
    volatile float ru=u;
    if(l==u) { neg_rl+=min; ru+=min; }
    if(neg_rl<neg_l) { neg_rl+=min; }
    if(ru<u) { ru+=min; }
    volatile float rl=-neg_l;
    Interval<FloatDPValue> res(rl,ru);
    FloatDP::set_rounding_mode(rnd);
    return res;
}

Interval<FloatDPValue> approximate_domain(Interval<FloatDPUpperBound> const& ivl) {
    auto rnd=FloatDP::get_rounding_mode();
    FloatDP::set_rounding_mode(to_nearest);
    volatile float eps=std::numeric_limits<float>::epsilon();
    volatile double l=ivl.lower().get_d();
    volatile double u=ivl.upper().get_d();
    volatile float rl=l;
    volatile float ru=u;
    if(rl==ru) { rl-=(rl*eps); ru+=(ru*eps); }
    Interval<FloatDPValue> res(rl,ru);
    FloatDP::set_rounding_mode(rnd);
    return res;
}

InputStream&
operator>>(InputStream& is, Interval<FloatDPValue>& ivl)
{
    FloatDP l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(not is.fail());
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(FloatDPValue(l),FloatDPValue(u));
    return is;
}


template<> FloatDPUpperInterval FloatUpperIntervalFactory<DoublePrecision>::create(ValidatedNumber const& y) const {
    return FloatDPUpperInterval(FloatDPBounds(y,this->_precision)); }

FloatDPBounds cast_singleton(FloatDPUpperInterval const& ivl) {
    return FloatDPBounds(ivl.lower(),ivl.upper()); }
FloatDPUpperInterval make_interval(FloatDPBounds const& x) {
    return FloatDPUpperInterval(x.lower(),x.upper()); }

FloatDPUpperInterval max(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return make_interval(max(cast_singleton(ivl1),cast_singleton(ivl2))); }
FloatDPUpperInterval min(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return make_interval(min(cast_singleton(ivl1),cast_singleton(ivl2))); }
FloatDPUpperInterval abs(FloatDPUpperInterval const& ivl) {
    return make_interval(abs(cast_singleton(ivl))); }
FloatDPUpperInterval nul(FloatDPUpperInterval const& ivl) {
    return make_interval(nul(cast_singleton(ivl))); }
FloatDPUpperInterval pos(FloatDPUpperInterval const& ivl) {
    return make_interval(pos(cast_singleton(ivl))); }
FloatDPUpperInterval neg(FloatDPUpperInterval const& ivl) {
    return make_interval(neg(cast_singleton(ivl))); }
FloatDPUpperInterval sqr(FloatDPUpperInterval const& ivl) {
    return make_interval(sqr(cast_singleton(ivl))); }
FloatDPUpperInterval rec(FloatDPUpperInterval const& ivl) {
    return make_interval(rec(cast_singleton(ivl))); }

FloatDPUpperInterval add(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return make_interval(add(cast_singleton(ivl1),cast_singleton(ivl2))); }
FloatDPUpperInterval sub(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return make_interval(sub(cast_singleton(ivl1),cast_singleton(ivl2))); }
FloatDPUpperInterval mul(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return make_interval(mul(cast_singleton(ivl1),cast_singleton(ivl2))); }
FloatDPUpperInterval div(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return make_interval(div(cast_singleton(ivl1),cast_singleton(ivl2))); }

FloatDPUpperInterval pow(FloatDPUpperInterval const& ivl, Nat m) {
    return make_interval(pow(cast_singleton(ivl),m)); }
FloatDPUpperInterval pow(FloatDPUpperInterval const& ivl, Int n) {
    return make_interval(pow(cast_singleton(ivl),n)); }

FloatDPUpperInterval sqrt(FloatDPUpperInterval const& ivl) {
    if(ivl.upper().raw()<0) { return FloatDPUpperInterval::empty_interval(); }
    else { FloatDPUpperBound u=sqrt(ivl.upper()); return FloatDPUpperInterval(-u,+u); } }
FloatDPUpperInterval exp(FloatDPUpperInterval const& ivl) {
    return make_interval(exp(cast_singleton(ivl))); }
FloatDPUpperInterval log(FloatDPUpperInterval const& ivl) {
    if(ivl.upper().raw()<=0) { return FloatDPUpperInterval::empty_interval(); }
    else if(ivl.lower().raw()<=0) { return FloatDPUpperInterval(-inf,log(up,ivl.upper().raw())); }
    else { return make_interval(log(cast_singleton(ivl))); } }
FloatDPUpperInterval sin(FloatDPUpperInterval const& ivl) {
    return make_interval(sin(cast_singleton(ivl))); }
FloatDPUpperInterval cos(FloatDPUpperInterval const& ivl) {
    return make_interval(cos(cast_singleton(ivl))); }
FloatDPUpperInterval tan(FloatDPUpperInterval const& ivl) {
    return make_interval(tan(cast_singleton(ivl))); }
FloatDPUpperInterval asin(FloatDPUpperInterval const& ivl) {
    return make_interval(asin(cast_singleton(ivl))); }
FloatDPUpperInterval acos(FloatDPUpperInterval const& ivl) {
    return make_interval(acos(cast_singleton(ivl))); }
FloatDPUpperInterval atan(FloatDPUpperInterval const& ivl) {
    return make_interval(atan(cast_singleton(ivl))); }

FloatDPError mag(FloatDPUpperInterval const& ivl) {
    return mag(cast_singleton(ivl)); }

ValidatedKleenean eq(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return eq(cast_singleton(ivl1),cast_singleton(ivl2)); }
ValidatedKleenean lt(FloatDPUpperInterval const& ivl1, FloatDPUpperInterval const& ivl2) {
    return lt(cast_singleton(ivl1),cast_singleton(ivl2)); }

FloatMPBounds cast_singleton(FloatMPUpperInterval const& ivl) {
    return FloatMPBounds(ivl.lower(),ivl.upper()); }
FloatMPUpperInterval make_interval(FloatMPBounds const& x) {
    return FloatMPUpperInterval(x.lower(),x.upper()); }
FloatMPError mag(FloatMPUpperInterval const& ivl) {
    return mag(cast_singleton(ivl)); }

template<> String class_name<FloatDPUpperInterval>() { return "FloatDPUpperInterval"; }
template<> String class_name<FloatMPUpperInterval>() { return "FloatMPUpperInterval"; }

} // namespace Ariadne
