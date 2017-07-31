/***************************************************************************
 *            interval.cpp
 *
 *  Copyright 2013--17  Pieter Collins
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

#include "interval.hpp"

namespace Ariadne {

Interval<Float64Value> widen_domain(Interval<Float64UpperBound> const& ivl) {
    auto rnd=Float64::get_rounding_mode();
    Float64::set_rounding_mode(upward);
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
    Interval<Float64Value> res(rl,ru);
    Float64::set_rounding_mode(rnd);
    return res;
}

Interval<Float64Value> approximate_domain(Interval<Float64UpperBound> const& ivl) {
    auto rnd=Float64::get_rounding_mode();
    Float64::set_rounding_mode(to_nearest);
    volatile float eps=std::numeric_limits<float>::epsilon();
    volatile double l=ivl.lower().get_d();
    volatile double u=ivl.upper().get_d();
    volatile float rl=l;
    volatile float ru=u;
    if(rl==ru) { rl-=(rl*eps); ru+=(ru*eps); }
    Interval<Float64Value> res(rl,ru);
    Float64::set_rounding_mode(rnd);
    return res;
}

InputStream&
operator>>(InputStream& is, Interval<Float64Value>& ivl)
{
    Float64 l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(not is.fail());
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(Float64Value(l),Float64Value(u));
    return is;
}


template<> Float64UpperInterval FloatUpperIntervalFactory<Precision64>::create(ValidatedNumber const& y) const {
    return Float64UpperInterval(Float64Bounds(y,this->_precision)); }

Float64Bounds cast_singleton(Float64UpperInterval const& ivl) {
    return Float64Bounds(ivl.lower(),ivl.upper()); }
Float64UpperInterval make_interval(Float64Bounds const& x) {
    return Float64UpperInterval(x.lower(),x.upper()); }

Float64UpperInterval max(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return make_interval(max(cast_singleton(ivl1),cast_singleton(ivl2))); }
Float64UpperInterval min(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return make_interval(min(cast_singleton(ivl1),cast_singleton(ivl2))); }
Float64UpperInterval abs(Float64UpperInterval const& ivl) {
    return make_interval(abs(cast_singleton(ivl))); }
Float64UpperInterval nul(Float64UpperInterval const& ivl) {
    return make_interval(nul(cast_singleton(ivl))); }
Float64UpperInterval pos(Float64UpperInterval const& ivl) {
    return make_interval(pos(cast_singleton(ivl))); }
Float64UpperInterval neg(Float64UpperInterval const& ivl) {
    return make_interval(neg(cast_singleton(ivl))); }
Float64UpperInterval sqr(Float64UpperInterval const& ivl) {
    return make_interval(sqr(cast_singleton(ivl))); }
Float64UpperInterval rec(Float64UpperInterval const& ivl) {
    return make_interval(rec(cast_singleton(ivl))); }

Float64UpperInterval add(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return make_interval(add(cast_singleton(ivl1),cast_singleton(ivl2))); }
Float64UpperInterval sub(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return make_interval(sub(cast_singleton(ivl1),cast_singleton(ivl2))); }
Float64UpperInterval mul(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return make_interval(mul(cast_singleton(ivl1),cast_singleton(ivl2))); }
Float64UpperInterval div(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return make_interval(div(cast_singleton(ivl1),cast_singleton(ivl2))); }

Float64UpperInterval pow(Float64UpperInterval const& ivl, Nat m) {
    return make_interval(pow(cast_singleton(ivl),m)); }
Float64UpperInterval pow(Float64UpperInterval const& ivl, Int n) {
    return make_interval(pow(cast_singleton(ivl),n)); }

Float64UpperInterval sqrt(Float64UpperInterval const& ivl) {
    if(ivl.upper().raw()<0) { return Float64UpperInterval::empty_interval(); }
    else { Float64UpperBound u=sqrt(ivl.upper()); return Float64UpperInterval(-u,+u); } }
Float64UpperInterval exp(Float64UpperInterval const& ivl) {
    return make_interval(exp(cast_singleton(ivl))); }
Float64UpperInterval log(Float64UpperInterval const& ivl) {
    if(ivl.upper().raw()<=0) { return Float64UpperInterval::empty_interval(); }
    else if(ivl.lower().raw()<=0) { return Float64UpperInterval(-inf,log(up,ivl.upper().raw())); }
    else { return make_interval(log(cast_singleton(ivl))); } }
Float64UpperInterval sin(Float64UpperInterval const& ivl) {
    return make_interval(sin(cast_singleton(ivl))); }
Float64UpperInterval cos(Float64UpperInterval const& ivl) {
    return make_interval(cos(cast_singleton(ivl))); }
Float64UpperInterval tan(Float64UpperInterval const& ivl) {
    return make_interval(tan(cast_singleton(ivl))); }
Float64UpperInterval asin(Float64UpperInterval const& ivl) {
    return make_interval(asin(cast_singleton(ivl))); }
Float64UpperInterval acos(Float64UpperInterval const& ivl) {
    return make_interval(acos(cast_singleton(ivl))); }
Float64UpperInterval atan(Float64UpperInterval const& ivl) {
    return make_interval(atan(cast_singleton(ivl))); }

Float64Error mag(Float64UpperInterval const& ivl) {
    return mag(cast_singleton(ivl)); }
Float64LowerBound mig(Float64UpperInterval const& ivl) {
    return mig(cast_singleton(ivl)); }

ValidatedKleenean eq(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return eq(cast_singleton(ivl1),cast_singleton(ivl2)); }
ValidatedKleenean lt(Float64UpperInterval const& ivl1, Float64UpperInterval const& ivl2) {
    return lt(cast_singleton(ivl1),cast_singleton(ivl2)); }

FloatMPBounds cast_singleton(FloatMPUpperInterval const& ivl) {
    return FloatMPBounds(ivl.lower(),ivl.upper()); }
FloatMPUpperInterval make_interval(FloatMPBounds const& x) {
    return FloatMPUpperInterval(x.lower(),x.upper()); }
FloatMPError mag(FloatMPUpperInterval const& ivl) {
    return mag(cast_singleton(ivl)); }
FloatMPLowerBound mig(FloatMPUpperInterval const& ivl) {
    return mig(cast_singleton(ivl)); }

template<> String class_name<Float64UpperInterval>() { return "Float64UpperInterval"; }
template<> String class_name<FloatMPUpperInterval>() { return "FloatMPUpperInterval"; }

} // namespace Ariadne
