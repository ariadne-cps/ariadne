/***************************************************************************
 *            float-user.cc
 *
 *  Copyright 2008-15  Pieter Collins
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

#include "config.h"
#include "utility/macros.h"
#include "utility/exceptions.h"

#include "numeric/float-user.h"

#include "numeric/integer.h"
#include "numeric/dyadic.h"
#include "numeric/decimal.h"
#include "numeric/rational.h"
#include "numeric/real.h"

#include "numeric/number_wrapper.h"

namespace Ariadne {


template<class PR> Nat Float<ApproximateTag,PR>::output_precision = 4;
template<class PR> Nat Float<BoundedTag,PR>::output_precision=8;
template<class PR> Nat Float<ExactTag,PR>::output_precision = 16;

const Float<ExactTag,Precision64> infty = Float<ExactTag,Precision64>(Float64::inf());

Float<ErrorTag,Precision64> operator"" _error(long double lx) {
    double x=lx;
    assert(x==lx);
    return Float<ErrorTag,Precision64>(Float64(x));
}

Float<ExactTag,Precision64> operator"" _exact(long double lx) {
    double x=lx;
    assert(x==lx);
    return Float<ExactTag,Precision64>(x);
}

Float<MetricTag,Precision64> operator"" _near(long double lx) {
    volatile double x=lx;
    volatile long double le=std::abs((long double)x-lx);
    volatile double e=le;
    while(e<le) { e*=(1+std::numeric_limits<double>::epsilon()); }
    return Float<MetricTag,Precision64>(x,e);
}

Float<UpperTag,Precision64> operator"" _upper(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x<lx) { x+=min; }
    while (x<lx) { x+=std::abs(x)*eps; }
    return Float<UpperTag,Precision64>(x);
}

Float<LowerTag,Precision64> operator"" _lower(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x>lx) { x-=min; }
    while (x>lx) { x-=std::abs(x)*eps; }
    return Float<LowerTag,Precision64>(x);
}

Float<ApproximateTag,Precision64> operator"" _approx(long double lx) {
    double x=lx;
    return Float<ApproximateTag,Precision64>(x);
}


TwoExp::operator Float<ExactTag,Precision64> () const {
    return Float<ExactTag,Precision64>(this->get_d());
}


/*
template<class PR> Float<ExactTag,PR>::Float(Rational const& q, PR pr)
    : _v(0.0,pr)
{
    Float<BoundedTag,PR> x(q,pr);
    if(x.lower_raw()==x.upper_raw()) { this->_v==x.value_raw(); }
    else { ARIADNE_THROW(std::runtime_error,"FloatValue(Rational q, Precision pr)","q="<<q<<" cannot be expressed exactly to precision \n"); }
}
*/

template<class PR> Float<ExactTag,PR>::operator Rational() const {
    return Rational(this->get_d());
}

template<class PR> Float<ExactTag,PR>::Float(TwoExp const& t, PR pr)
    : _v(pow(RawFloat<PR>(2.0,pr),t.exponent()))
{
}

template<class PR> Float<ExactTag,PR>::Float(Integer const& z, PR pr)
    : _v(Rational(z),RawFloat<PR>::to_nearest,pr)
{
    Rational q(_v);
    ARIADNE_PRECONDITION(z==q);
}

template<class PR> Float<ExactTag,PR>::operator Number<ExactTag>() const {
    return Number<ExactTag>(new NumberWrapper<Float<ExactTag,PR>>(*this));
}

template<class PR> Float<MetricTag,PR>::Float(Rational const& q, PR pr)
    : _v(RawFloat<PR>(q,RawFloat<PR>::to_nearest,pr)), _e(abs(Rational(_v)-q),RawFloat<PR>::upward,pr) {
}

template<class PR> Float<MetricTag,PR>::Float(Real const& r, PR pr)
    : Float(r.get(pr)) {
}

template<class PR> Float<MetricTag,PR>::Float(Number<ValidatedTag> const& x, PR pr)
    : Float(x.get(MetricTag(),pr)) {
}

template<class PR> Float<MetricTag,PR>::operator Number<ValidatedTag>() const {
    return Number<ValidatedTag>(new NumberWrapper<Float<MetricTag,PR>>(*this));
}


template<class PR> Float<BoundedTag,PR>::Float(Integer const& z, PR pr)
    : Float(RawFloat<PR>(z,RawFloat<PR>::downward,pr),RawFloat<PR>(z,RawFloat<PR>::upward,pr)) {
}

template<class PR> Float<BoundedTag,PR>::Float(Rational const& q, PR pr)
    : Float(RawFloat<PR>(q,RawFloat<PR>::downward,pr),RawFloat<PR>(q,RawFloat<PR>::upward,pr)) {
}

template<class PR> Float<BoundedTag,PR>::Float(Rational const& ql, Rational const& qu, PR pr)
    : Float(RawFloat<PR>(ql,RawFloat<PR>::downward,pr),RawFloat<PR>(qu,RawFloat<PR>::upward,pr)) {
}

template<class PR> Float<BoundedTag,PR>::Float(Real const& x, PR pr)
    : Float(x.get(pr)) {
}

template<class PR> Float<BoundedTag,PR>::Float(Number<ValidatedTag> const& x, PR pr)
    : Float(x.get(BoundedTag(),pr)) {
}

template<class PR> Float<BoundedTag,PR>::operator Number<ValidatedTag>() const {
    return Number<ValidatedTag>(new NumberWrapper<Float<BoundedTag,PR>>(*this));
}

template<class PR> Float<UpperTag,PR>::Float(Rational const& q, PR pr)
    : Float(Float<BoundedTag,PR>(q,pr)) {
}

template<class PR> Float<UpperTag,PR>::Float(Real const& r, PR pr)
    : Float(r.get(pr)) {
}

template<class PR> Float<UpperTag,PR>::Float(Number<ValidatedUpperTag> const& x, PR pr)
    : Float(x.get(UpperTag(),pr)) {
}

template<class PR> Float<UpperTag,PR>::operator Number<ValidatedUpperTag>() const {
    return Number<ValidatedUpperTag>(new NumberWrapper<Float<UpperTag,PR>>(*this));
}

template<class PR> Float<LowerTag,PR>::Float(Rational const& q, PR pr)
    : Float(Float<BoundedTag,PR>(q,pr)) {
}

template<class PR> Float<LowerTag,PR>::Float(Real const& r, PR pr)
    : Float(r.get(pr)) {
}

template<class PR> Float<LowerTag,PR>::Float(Number<ValidatedLowerTag> const& x, PR pr)
    : Float(x.get(LowerTag(),pr)) {
}

template<class PR> Float<LowerTag,PR>::operator Number<ValidatedLowerTag>() const {
    return Number<ValidatedLowerTag>(new NumberWrapper<Float<LowerTag,PR>>(*this));
}

template<class PR> Float<ApproximateTag,PR>::Float(Rational const& q, PR pr)
    : Float(Float<BoundedTag,PR>(q,pr)) {
}

template<class PR> Float<ApproximateTag,PR>::Float(Real const& r, PR pr)
    : Float(r.get(pr)) {
}

template<class PR> Float<ApproximateTag,PR>::Float(Number<ApproximateTag> const& x, PR pr)
    : Float(x.get(ApproximateTag(),pr)) {
}

template<class PR> Float<ApproximateTag,PR>::operator Number<ApproximateTag>() const {
    return Number<ApproximateTag>(new NumberWrapper<Float<ApproximateTag,PR>>(*this));
}



template<class PR> Float<ExactTag,PR>::Float(Integer const& z) : Float(z,RawFloat<PR>::get_default_precision()) { }

template<class PR> Float<MetricTag,PR>::Float(Integer const& z) : Float(Rational(z)) { }
template<class PR> Float<MetricTag,PR>::Float(Rational const& q) : Float(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> Float<MetricTag,PR>::Float(Number<ValidatedTag> const& x) : Float(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> Float<BoundedTag,PR>::Float(Integer const& z) : Float(Rational(z)) { }
template<class PR> Float<BoundedTag,PR>::Float(Dyadic const& b) : Float(Rational(b)) { }
template<class PR> Float<BoundedTag,PR>::Float(Decimal const& d) : Float(Rational(d)) { }
template<class PR> Float<BoundedTag,PR>::Float(Rational const& q) : Float(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> Float<BoundedTag,PR>::Float(Real const& x) : Float(x,RawFloat<PR>::get_default_precision()) { }
template<class PR> Float<BoundedTag,PR>::Float(Number<ValidatedTag> const& x) : Float(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> Float<UpperTag,PR>::Float(Integer const& z) : Float(Rational(z)) { }
template<class PR> Float<UpperTag,PR>::Float(Rational const& q) : Float(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> Float<UpperTag,PR>::Float(Number<ValidatedUpperTag> const& x) : Float(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> Float<LowerTag,PR>::Float(Integer const& z) : Float(Rational(z)) { }
template<class PR> Float<LowerTag,PR>::Float(Rational const& q) : Float(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> Float<LowerTag,PR>::Float(Number<ValidatedLowerTag> const& x) : Float(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> Float<ApproximateTag,PR>::Float(Integer const& z) : Float(Rational(z)) { }
template<class PR> Float<ApproximateTag,PR>::Float(Dyadic const& b) : Float(Rational(b)) { }
template<class PR> Float<ApproximateTag,PR>::Float(Decimal const& d) : Float(Rational(d)) { }
template<class PR> Float<ApproximateTag,PR>::Float(Rational const& q) : Float(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> Float<ApproximateTag,PR>::Float(Number<ApproximateTag> const& x) : Float(x,RawFloat<PR>::get_default_precision()) { }



template<class PR> Float<ApproximateTag,PR>::Float(Float<LowerTag,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<ApproximateTag,PR>::Float(Float<UpperTag,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<ApproximateTag,PR>::Float(Float<BoundedTag,PR> const& x) : _a(x.value_raw()) {
}

template<class PR> Float<ApproximateTag,PR>::Float(Float<MetricTag,PR> const& x) : _a(x.value_raw()) {
}

template<class PR> Float<ApproximateTag,PR>::Float(Float<ExactTag,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<LowerTag,PR>::Float(Float<BoundedTag,PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> Float<LowerTag,PR>::Float(Float<MetricTag,PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> Float<LowerTag,PR>::Float(Float<ExactTag,PR> const& x) : _l(x.raw()) {
}

template<class PR> Float<UpperTag,PR>::Float(Float<BoundedTag,PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> Float<UpperTag,PR>::Float(Float<MetricTag,PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> Float<UpperTag,PR>::Float(Float<ExactTag,PR> const& x) : _u(x.raw()) {
}

template<class PR> Float<BoundedTag,PR>::Float(Float<MetricTag,PR> const& x) : _l(x.lower_raw()), _u(x.upper_raw()) {
}

template<class PR> Float<BoundedTag,PR>::Float(Float<ExactTag,PR> const& x) : _l(x.raw()), _u(x.raw()) {
}

template<class PR> Float<MetricTag,PR>::Float(Float<BoundedTag,PR> const& x) : _v(x.value_raw()), _e(x.error_raw()) {
}

template<class PR> Float<MetricTag,PR>::Float(Float<ExactTag,PR> const& x) : _v(x.raw()), _e(nul(x.raw())) {
}


template<class PR> Float<MetricTag,PR> Float<ExactTag,PR>::pm(Float<ErrorTag,PR> _e) const {
    Float<ExactTag,PR> const& _v=*this; return Float<MetricTag,PR>(_v,_e);
}




template<class PR> Float<ApproximateTag,PR> floor(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(floor(x._a)); }
template<class PR> Float<ApproximateTag,PR> ceil(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(ceil(x._a)); }
template<class PR> Float<ApproximateTag,PR> round(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(round(x._a)); }

template<class PR> Float<ApproximateTag,PR> abs(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(abs_exact(x._a)); }
template<class PR> Float<ApproximateTag,PR> max(Float<ApproximateTag,PR> const& x, Float<ApproximateTag,PR> y) {
    return Float<ApproximateTag,PR>(max_exact(x._a,y._a)); }
template<class PR> Float<ApproximateTag,PR> min(Float<ApproximateTag,PR> const& x, Float<ApproximateTag,PR> y) {
    return Float<ApproximateTag,PR>(min_exact(x._a,y._a)); }
template<class PR> Float<PositiveApproximateTag,PR> mag(Float<ApproximateTag,PR> const& x) {
    return Float<PositiveApproximateTag,PR>(abs(x._a)); }
template<class PR> Float<PositiveApproximateTag,PR> mig(Float<ApproximateTag,PR> const& x) {
    return Float<PositiveApproximateTag,PR>(abs(x._a)); }

template<class PR> Float<ApproximateTag,PR> nul(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(nul_exact(x._a)); }
template<class PR> Float<ApproximateTag,PR> pos(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(pos_exact(x._a)); }
template<class PR> Float<ApproximateTag,PR> neg(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(neg_exact(x._a)); }
template<class PR> Float<ApproximateTag,PR> half(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(half_exact(x._a)); }
template<class PR> Float<ApproximateTag,PR> sqr(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(mul_near(x._a,x._a)); }
template<class PR> Float<ApproximateTag,PR> rec(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(div_near(1.0,x._a)); }

template<class PR> Float<ApproximateTag,PR> add(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return Float<ApproximateTag,PR>(add_near(x1._a,x2._a)); }
template<class PR> Float<ApproximateTag,PR> sub(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return Float<ApproximateTag,PR>(sub_near(x1._a,x2._a)); }
template<class PR> Float<ApproximateTag,PR> mul(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return Float<ApproximateTag,PR>(mul_near(x1._a,x2._a)); }
template<class PR> Float<ApproximateTag,PR> div(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return Float<ApproximateTag,PR>(div_near(x1._a,x2._a)); }

template<class PR> Float<ApproximateTag,PR> pow(Float<ApproximateTag,PR> const& x, Nat m) {
    return Float<ApproximateTag,PR>(pow_approx(x._a,m)); }
template<class PR> Float<ApproximateTag,PR> pow(Float<ApproximateTag,PR> const& x, Int n) {
    return Float<ApproximateTag,PR>(pow_approx(x._a,n)); }

template<class PR> Float<ApproximateTag,PR> sqrt(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(sqrt_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> exp(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(exp_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> log(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(log_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> sin(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(sin_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> cos(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(cos_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> tan(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(tan_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> asin(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(asin_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> acos(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(acos_approx(x._a)); }
template<class PR> Float<ApproximateTag,PR> atan(Float<ApproximateTag,PR> const& x) {
    return Float<ApproximateTag,PR>(atan_approx(x._a)); }

template<class PR> Logical<ApproximateTag> eq(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return x1._a==x2._a; }
template<class PR> Logical<ApproximateTag> leq(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return x1._a<=x2._a; }

template<class PR> Bool same(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) {
    return x1._a==x2._a; }

template<class PR> Float<ApproximateTag,PR> operator+(Float<ApproximateTag,PR> const& x) { return pos(x); }
template<class PR> Float<ApproximateTag,PR> operator-(Float<ApproximateTag,PR> const& x) { return neg(x); }
template<class PR> Float<ApproximateTag,PR> operator+(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<ApproximateTag,PR> operator-(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<ApproximateTag,PR> operator*(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<ApproximateTag,PR> operator/(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<ApproximateTag,PR>& operator+=(Float<ApproximateTag,PR>& x1, Float<ApproximateTag,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<ApproximateTag,PR>& operator-=(Float<ApproximateTag,PR>& x1, Float<ApproximateTag,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<ApproximateTag,PR>& operator*=(Float<ApproximateTag,PR>& x1, Float<ApproximateTag,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<ApproximateTag,PR>& operator/=(Float<ApproximateTag,PR>& x1, Float<ApproximateTag,PR> const& x2) { return x1=x1/x2; }

template<class PR> Fuzzy operator==(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return x1._a==x2._a; }
template<class PR> Fuzzy operator!=(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return x1._a!=x2._a; }
template<class PR> Fuzzy operator<=(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return x1._a<=x2._a; }
template<class PR> Fuzzy operator>=(Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return x1._a>=x2._a; }
template<class PR> Fuzzy operator< (Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return x1._a< x2._a; }
template<class PR> Fuzzy operator> (Float<ApproximateTag,PR> const& x1, Float<ApproximateTag,PR> const& x2) { return x1._a> x2._a; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<ApproximateTag,PR> const& x) {
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_to_nearest();
    os << std::showpoint << std::setprecision(Float<ApproximateTag,PR>::output_precision) << x.raw();
    RawFloat<PR>::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<ApproximateTag,PR>& x) {
    is >> x._a;
    return is;
}



template<class PR> Float<LowerTag,PR> max(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return Float<LowerTag,PR>(max_exact(x1._l,x2._l)); }
template<class PR> Float<LowerTag,PR> min(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return Float<LowerTag,PR>(min_exact(x1._l,x2._l)); }
template<class PR> Float<ApproximateTag,PR> abs(Float<LowerTag,PR> const& x) {
    return abs(Float<ApproximateTag,PR>(x)); }

template<class PR> Float<LowerTag,PR> nul(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(pos_exact(x._l)); }
template<class PR> Float<LowerTag,PR> pos(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(pos_exact(x._l)); }
template<class PR> Float<LowerTag,PR> neg(Float<UpperTag,PR> const& x) {
    return Float<LowerTag,PR>(neg_exact(x._u)); }
template<class PR> Float<LowerTag,PR> half(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(half_exact(x._l)); }

template<class PR> Float<LowerTag,PR> rec(Float<UpperTag,PR> const& x) {
    return Float<LowerTag,PR>(rec_down(x.raw())); }

template<class PR> Float<LowerTag,PR> add(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return Float<LowerTag,PR>(add_down(x1._l,x2._l)); }

template<class PR> Float<ApproximateTag,PR> sub(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return Float<UpperTag,PR>(sub_near(x1._l,x2._l)); }

template<class PR> Float<LowerTag,PR> sub(Float<LowerTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return Float<LowerTag,PR>(sub_down(x1._l,x2._u)); }

template<class PR> Float<LowerTag,PR> mul(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    return Float<LowerTag,PR>(mul_down(x1.raw(),x2.raw())); }

template<class PR> Float<LowerTag,PR> div(Float<LowerTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return Float<LowerTag,PR>(div_down(x1.raw(),x2.raw())); }

template<class PR> Float<LowerTag,PR> pow(Float<LowerTag,PR> const& x, Nat m) {
    ARIADNE_PRECONDITION(x.raw()>=0);
    return Float<LowerTag,PR>(pow_down(x.raw(),m)); }

template<class PR> Float<ApproximateTag,PR> pow(Float<LowerTag,PR> const& x, Int n) {
    return pow(Float<ApproximateTag,PR>(x),n); }

template<class PR> Float<LowerTag,PR> sqrt(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(sqrt_down(x.raw())); }

template<class PR> Float<LowerTag,PR> exp(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(exp_down(x.raw())); }

template<class PR> Float<LowerTag,PR> log(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(log_down(x.raw())); }

template<class PR> Float<LowerTag,PR> atan(Float<LowerTag,PR> const& x) {
    return Float<LowerTag,PR>(atan_down(x.raw())); }

template<class PR> Logical<LowerTag> eq(Float<LowerTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    if(x1._l>x2._u) { return false; }
    else { return Logical<LowerTag>(LogicalValue::INDETERMINATE); }
}

template<class PR> Logical<LowerTag> leq(Float<LowerTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    if(x1._l>x2._u) { return false; }
    else { return Logical<LowerTag>(LogicalValue::LIKELY); }
}

template<class PR> Bool same(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return x1._l==x2._l; }

template<class PR> Float<LowerTag,PR> operator+(Float<LowerTag,PR> const& x) { return pos(x); }
template<class PR> Float<LowerTag,PR> operator-(Float<UpperTag,PR> const& x) { return neg(x); }
template<class PR> Float<LowerTag,PR> operator+(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<LowerTag,PR> operator-(Float<LowerTag,PR> const& x1, Float<UpperTag,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<LowerTag,PR> operator*(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<LowerTag,PR> operator/(Float<LowerTag,PR> const& x1, Float<UpperTag,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<LowerTag,PR>& operator+=(Float<LowerTag,PR>& x1, Float<LowerTag,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<LowerTag,PR>& operator-=(Float<LowerTag,PR>& x1, Float<UpperTag,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<LowerTag,PR>& operator*=(Float<LowerTag,PR>& x1, Float<LowerTag,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<LowerTag,PR>& operator/=(Float<LowerTag,PR>& x1, Float<UpperTag,PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<LowerTag,PR> const& x) {
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_downward();
    os << std::showpoint << std::setprecision(Float<BoundedTag,PR>::output_precision) << x.raw();
    RawFloat<PR>::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<LowerTag,PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
}



template<class PR> Float<UpperTag,PR> max(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return Float<UpperTag,PR>(max_exact(x1._u,x2._u)); }

template<class PR> Float<UpperTag,PR> min(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return Float<UpperTag,PR>(min_exact(x1._u,x2._u)); }

template<class PR> Float<ApproximateTag,PR> abs(Float<UpperTag,PR> const& x) {
    return abs(Float<ApproximateTag,PR>(x)); }

template<class PR> Float<UpperTag,PR> nul(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(pos_exact(x._u)); }

template<class PR> Float<UpperTag,PR> pos(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(pos_exact(x._u)); }

template<class PR> Float<UpperTag,PR> neg(Float<LowerTag,PR> const& x) {
    return Float<UpperTag,PR>(neg_exact(x._l)); }

template<class PR> Float<UpperTag,PR> half(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(half_exact(x._u)); }

template<class PR> Float<UpperTag,PR> sqr(Float<UpperTag,PR> const& x) {
    ARIADNE_ASSERT(false); return Float<UpperTag,PR>(mul_up(x._u,x._u)); }

template<class PR> Float<UpperTag,PR> rec(Float<LowerTag,PR> const& x) {
    return Float<UpperTag,PR>(rec_up(x.raw())); }

template<class PR> Float<UpperTag,PR> add(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return Float<UpperTag,PR>(add_up(x1._u,x2._u)); }

template<class PR> Float<ApproximateTag,PR> sub(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return Float<UpperTag,PR>(sub_near(x1._u,x2._u)); }

template<class PR> Float<UpperTag,PR> sub(Float<UpperTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return Float<UpperTag,PR>(sub_up(x1._u,x2._l)); }

template<class PR> Float<UpperTag,PR> mul(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
//    ARIADNE_WARN("Multiplying FloatUpperBound "<<x1<<" with FloatUpperBound "<<x2<<" is unsafe");
    ARIADNE_PRECONDITION(x1.raw()>=0);
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return Float<UpperTag,PR>(mul_up(x1._u,x2._u)); }

template<class PR> Float<UpperTag,PR> div(Float<UpperTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
//    ARIADNE_WARN("Dividing FloatUpperBound "<<x1<<" by FloatLowerBound "<<x2<<" is unsafe");
    ARIADNE_PRECONDITION(x1.raw()>=0);
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return Float<UpperTag,PR>(div_up(x1._u,x2._l)); }

template<class PR> Float<UpperTag,PR> pow(Float<UpperTag,PR> const& x, Nat m) {
    ARIADNE_PRECONDITION(x.raw()>=0);
    return Float<UpperTag,PR>(pow_up(x._u,m)); }

template<class PR> Float<ApproximateTag,PR> pow(Float<UpperTag,PR> const& x, Int n) {
    return pow(Float<ApproximateTag,PR>(x),n); }

template<class PR> Float<UpperTag,PR> sqrt(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(sqrt_up(x.raw())); }

template<class PR> Float<UpperTag,PR> exp(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(exp_up(x.raw())); }

template<class PR> Float<UpperTag,PR> log(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(log_up(x.raw())); }

template<class PR> Float<UpperTag,PR> atan(Float<UpperTag,PR> const& x) {
    return Float<UpperTag,PR>(atan_up(x.raw())); }

template<class PR> Logical<LowerTag> eq(Float<UpperTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    if(x1._u<x2._l) { return false; }
    else { return Logical<LowerTag>(LogicalValue::INDETERMINATE); }
}

template<class PR> Logical<UpperTag> leq(Float<UpperTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    if(x1._u<=x2._l) { return true; }
    else { return Logical<UpperTag>(LogicalValue::UNLIKELY); }
}

template<class PR> Bool same(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return x1._u==x2._u; }

template<class PR> Integer integer_cast(Float<UpperTag,PR> const& x) { return Integer(static_cast<int>(x._u.get_d())); }

template<class PR> Float<UpperTag,PR> operator+(Float<UpperTag,PR> const& x) { return pos(x); }
template<class PR> Float<UpperTag,PR> operator-(Float<LowerTag,PR> const& x) { return neg(x); }
template<class PR> Float<UpperTag,PR> operator+(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<UpperTag,PR> operator-(Float<UpperTag,PR> const& x1, Float<LowerTag,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<UpperTag,PR> operator*(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<UpperTag,PR> operator/(Float<UpperTag,PR> const& x1, Float<LowerTag,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<UpperTag,PR>& operator+=(Float<UpperTag,PR>& x1, Float<UpperTag,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<UpperTag,PR>& operator-=(Float<UpperTag,PR>& x1, Float<LowerTag,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<UpperTag,PR>& operator*=(Float<UpperTag,PR>& x1, Float<UpperTag,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<UpperTag,PR>& operator/=(Float<UpperTag,PR>& x1, Float<LowerTag,PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<UpperTag,PR> const& x) {
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_upward();
    os << std::showpoint << std::setprecision(Float<BoundedTag,PR>::output_precision) << x.raw();
    RawFloat<PR>::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<UpperTag,PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
}








template<class PR> Float<BoundedTag,PR> round(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(round(x.lower_raw()),round(x.upper_raw()));
}

template<class PR> Float<BoundedTag,PR> max(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
}

template<class PR> Float<BoundedTag,PR> min(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
}


template<class PR> Float<BoundedTag,PR> abs(Float<BoundedTag,PR> const& x) {
    if(x.lower_raw()>=0) {
        return Float<BoundedTag,PR>(x.lower_raw(),x.upper_raw());
    } else if(x.upper_raw()<=0) {
        return Float<BoundedTag,PR>(neg(x.upper_raw()),neg(x.lower_raw()));
    } else {
        return Float<BoundedTag,PR>(static_cast<RawFloat<PR>>(0.0,x.precision()),max(neg(x.lower_raw()),x.upper_raw()));
    }
}

template<class PR> Float<PositiveLowerTag,PR> mig(Float<BoundedTag,PR> const& x) {
    return Float<PositiveLowerTag,PR>(max(0,max(x._l,neg(x._u))));
}

template<class PR> Float<PositiveUpperTag,PR> mag(Float<BoundedTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(max(neg(x._l),x._u));
}

template<class PR> Float<BoundedTag,PR> nul(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(nul(x._l),nul(x._u));
}

template<class PR> Float<BoundedTag,PR> pos(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(pos(x._l),pos(x._u));
}

template<class PR> Float<BoundedTag,PR> neg(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(neg(x._u),neg(x._l));
}

template<class PR> Float<BoundedTag,PR> half(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(half(x._l),half(x._u));
}

template<class PR> Float<BoundedTag,PR> sqr(Float<BoundedTag,PR> const& x) {
    const RawFloat<PR>& xl=x.lower_raw(); const RawFloat<PR>& xu=x.upper_raw();
    RawFloat<PR> rl,ru;
    if(xl>0.0) {
        rl=mul_down(xl,xl); ru=mul_up(xu,xu);
    } else if(xu<0.0) {
        rl=mul_down(xu,xu); ru=mul_up(xl,xl);
    } else {
        rl=nul(xl); ru=max(mul_up(xl,xl),mul_up(xu,xu));
    }
    return Float<BoundedTag,PR>(rl,ru);
}

template<class PR> Float<BoundedTag,PR> rec(Float<BoundedTag,PR> const& x) {
    // IMPORTANT: Need to be careful when one of the bounds is 0, since if xl=-0.0 and xu>0, then 1/xl=-inf
    if(x._l>0 || x._u<0) {
        return Float<BoundedTag,PR>(rec_down(x._u),rec_up(x._l));
    } else {
        RawFloat<PR> inf=RawFloat<PR>::inf(x.precision());
        RawFloat<PR> rl=-inf; RawFloat<PR> ru=+inf;
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds rec(FloatBounds x)","x="<<x);
        return Float<BoundedTag,PR>(-inf,+inf);
    }
}

template<class PR> Float<BoundedTag,PR> add(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(add_down(x1._l,x2._l),add_up(x1._u,x2._u));
}

template<class PR> Float<BoundedTag,PR> sub(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(sub_down(x1._l,x2._u),sub_up(x1._u,x2._l));
}

template<class PR> Float<BoundedTag,PR> mul(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    const RawFloat<PR>& x1l=x1._l; const RawFloat<PR>& x1u=x1._u;
    const RawFloat<PR>& x2l=x2._l; const RawFloat<PR>& x2u=x2._u;
    RawFloat<PR> rl,ru;
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    if(x1l>=0) {
        if(x2l>=0) {
            rl=mul_down(x1l,x2l); ru=mul_up(x1u,x2u);
        } else if(x2u<=0) {
            rl=mul_down(x1u,x2l); ru=mul_up(x1l,x2u);
        } else {
            rl=mul_down(x1u,x2l); ru=mul_up(x1u,x2u);
        }
    }
    else if(x1u<=0) {
        if(x2l>=0) {
            rl=mul_down(x1l,x2u); ru=mul_up(x1u,x2l);
        } else if(x2u<=0) {
            rl=mul_down(x1u,x2u); ru=mul_up(x1l,x2l);
        } else {
            rl=mul_down(x1l,x2u); ru=mul_up(x1l,x2l);
        }
    } else {
        if(x2l>=0) {
            rl=mul_down(x1l,x2u); ru=mul_up(x1u,x2u);
        } else if(x2u<=0) {
            rl=mul_down(x1u,x2l); ru=mul_up(x1l,x2l);
        } else {
            rl=min(mul_down(x1u,x2l),mul_down(x1l,x2u));
            ru=max(mul_up(x1l,x2l),mul_up(x1u,x2u));
        }
    }
    return Float<BoundedTag,PR>(rl,ru);
}

template<class PR> Float<BoundedTag,PR> div(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    const RawFloat<PR>& x1l=x1.lower_raw(); const RawFloat<PR>& x1u=x1.upper_raw();
    const RawFloat<PR>& x2l=x2.lower_raw(); const RawFloat<PR>& x2u=x2.upper_raw();
    RawFloat<PR> rl,ru;

    // IMPORTANT: Need to be careful when one of the bounds is 0, since if x2l=-0.0 and x1u>0, then x2l>=0 but x1u/x2l=-inf
    if(x2l>0) {
        if(x1l>=0) {
            rl=div_down(x1l,x2u); ru=div_up(x1u,x2l);
        } else if(x1u<=0) {
            rl=div_down(x1l,x2l); ru=div_up(x1u,x2u);
        } else {
            rl=div_down(x1l,x2l); ru=div_up(x1u,x2l);
        }
    }
    else if(x2u<0) {
        if(x1l>=0) {
            rl=div_down(x1u,x2u); ru=div_up(x1l,x2l);
        } else if(x1u<=0) {
            rl=div_down(x1u,x2l); ru=div_up(x1l,x2u);
        } else {
            rl=div_down(x1u,x2u); ru=div_up(x1l,x2u);
        }
    }
    else {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
        rl=-RawFloat<PR>::inf();
        ru=+RawFloat<PR>::inf();
    }
    return Float<BoundedTag,PR>(rl,ru);
}






template<class PR> Float<BoundedTag,PR> pow(Float<BoundedTag,PR> const& x, Int n) {
    if(n<0) { return pow(rec(x),Nat(-n)); }
    else return pow(x,Nat(n));
}

template<class PR> Float<BoundedTag,PR> pow(Float<BoundedTag,PR> const& x, Nat m) {
    Float<BoundedTag,PR> y = x;
    if(m%2==0) { y=abs(x); }
    RawFloat<PR> rl=pow_down(y.lower_raw(),m);
    RawFloat<PR> ru=pow_up(y.upper_raw(),m);
    return Float<BoundedTag,PR>(rl,ru);
}


template<class PR> Float<BoundedTag,PR> sqrt(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(sqrt_down(x.lower_raw()),sqrt_up(x.upper_raw()));
}

template<class PR> Float<BoundedTag,PR> exp(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(exp_down(x.lower_raw()),exp_up(x.upper_raw()));
}

template<class PR> Float<BoundedTag,PR> log(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(log_down(x.lower_raw()),log_up(x.upper_raw()));
}


template<class PR> Float<BoundedTag,PR> pi_val(PR pr) { return Float<BoundedTag,PR>(pi_down(pr),pi_up(pr)); }

template<class PR> Float<BoundedTag,PR> sin(Float<BoundedTag,PR> const& x)
{
    return cos(x-half(pi_val<PR>(x.precision())));
}

template<class PR> Float<BoundedTag,PR> cos(Float<BoundedTag,PR> const& x)
{
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    typename RawFloat<PR>::RoundingModeType rnd = RawFloat<PR>::get_rounding_mode();
    PR prec=x.precision();

    static const RawFloat<PR> one(1,prec);
    static const Float<ExactTag,PR> two(2,prec);

    if(x.error().raw()>2*pi_down(prec)) { return Float<BoundedTag,PR>(-one,+one); }

    auto n=floor(x.lower_raw()/(2*pi_approx(prec))+one/2);
    Float<BoundedTag,PR> y=x-two*Float<ExactTag,PR>(n)*pi_val<PR>(x.precision());

    ARIADNE_ASSERT(y.lower_raw()<=pi_up(prec));
    ARIADNE_ASSERT(y.upper_raw()>=-pi_up(prec));

    RawFloat<PR> rl,ru;
    if(y.lower_raw()<=-pi_down(prec)) {
        if(y.upper_raw()<=0.0) { rl=-one; ru=cos_up(y.upper_raw()); }
        else { rl=-one; ru=+one; }
    } else if(y.lower_raw()<=0.0) {
        if(y.upper_raw()<=0.0) { rl=cos_down(y.lower_raw()); ru=cos_up(y.upper_raw()); }
        else if(y.upper_raw()<=pi_down(prec)) { rl=cos_down(max(-y.lower_raw(),y.upper_raw())); ru=+one; }
        else { rl=-one; ru=+one; }
    } else if(y.lower_raw()<=pi_up(prec)) {
        if(y.upper_raw()<=pi_down(prec)) { rl=cos_down(y.upper_raw()); ru=cos_up(y.lower_raw()); }
        else if(y.upper_raw()<=2*pi_down(prec)) { rl=-one; ru=cos_up(min(y.lower_raw(),sub_down(2*pi_down(prec),y.upper_raw()))); }
        else { rl=-one; ru=+one; }
    } else {
        assert(false);
    }

    RawFloat<PR>::set_rounding_mode(rnd);
    return Float<BoundedTag,PR>(rl,ru);
}

template<class PR> Float<BoundedTag,PR> tan(Float<BoundedTag,PR> const& x) {
    return mul(sin(x),rec(cos(x)));
}

template<class PR> Float<BoundedTag,PR> asin(Float<BoundedTag,PR> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class PR> Float<BoundedTag,PR> acos(Float<BoundedTag,PR> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class PR> Float<BoundedTag,PR> atan(Float<BoundedTag,PR> const& x) {
    return Float<BoundedTag,PR>(atan_down(x._l),atan_up(x._u));
}



//! \related Float<BoundedTag,PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> Logical<ValidatedTag> eq(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    if(x1.upper_raw()<x2.lower_raw() || x1.lower_raw()>x2.upper_raw()) { return false; }
    else if(x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw()) { return true; }
    else { return indeterminate; }
}

//! \related Float<BoundedTag,PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> Logical<ValidatedTag> leq(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    if(x1.upper_raw()<=x2.lower_raw()) { return true; }
    else if(x1.lower_raw()> x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

template<class PR> Logical<ValidatedTag> operator==(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return eq(x1,x2);
}

template<class PR> Logical<ValidatedTag> operator!=(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return !eq(x1,x2); }

template<class PR> Logical<ValidatedTag> operator<=(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return leq(x1,x2);
}

template<class PR> Logical<ValidatedTag> operator>=(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return leq(x2,x1);
}

template<class PR> Logical<ValidatedTag> operator< (Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return !leq(x2,x1);
}

template<class PR> Logical<ValidatedTag> operator> (Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return !leq(x1,x2);
}

template<class PR> Float<BoundedTag,PR> widen(Float<BoundedTag,PR> const& x)
{
    typename RawFloat<PR>::RoundingModeType rm=RawFloat<PR>::get_rounding_mode();
    const RawFloat<PR>& xl=x.lower_raw();
    const RawFloat<PR>& xu=x.upper_raw();
    const RawFloat<PR> m=std::numeric_limits<float>::min();
    RawFloat<PR>::set_rounding_upward();
    RawFloat<PR> wu=add(xu,m);
    RawFloat<PR> mwl=add(neg(xl),m);
    RawFloat<PR> wl=neg(mwl);
    RawFloat<PR>::set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return Float<BoundedTag,PR>(wl,wu);
}

template<class PR> Float<BoundedTag,PR> narrow(Float<BoundedTag,PR> const& x)
{
    typename RawFloat<PR>::RoundingModeType rm=RawFloat<PR>::get_rounding_mode();
    const RawFloat<PR>& xl=x.lower_raw();
    const RawFloat<PR>& xu=x.upper_raw();
    const RawFloat<PR> m=std::numeric_limits<float>::min();
    RawFloat<PR>::set_rounding_upward();
    RawFloat<PR> mnu=add(neg(xu),m);
    RawFloat<PR> nu=neg(mnu);
    RawFloat<PR> nl=add(xl,m);
    RawFloat<PR>::set_rounding_mode(rm);
    assert(xl<nl); assert(nu<xu);
    return Float<BoundedTag,PR>(nl,nu);
}

template<class PR> Float<BoundedTag,PR> trunc(Float<BoundedTag,PR> const& x)
{
    typename RawFloat<PR>::RoundingModeType rm=RawFloat<PR>::get_rounding_mode();
    const double& xl=x.lower_raw().get_d();
    const double& xu=x.upper_raw().get_d();
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { RawFloat<PR>::set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { RawFloat<PR>::set_rounding_downward(); tl-=fm; }
    RawFloat<PR>::set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return Float<BoundedTag,PR>(double(tl),double(tu));
}

template<class PR> Float<BoundedTag,PR> trunc(Float<BoundedTag,PR> const& x, Nat n)
{
    Float<BoundedTag,PR> _e=Float<BoundedTag,PR>(std::pow(2.0,52-(Int)n));
    Float<BoundedTag,PR> y=x+_e;
    return y-_e;
}

template<class PR> Integer integer_cast(Float<BoundedTag,PR> const& x) {
    return Integer(static_cast<int>(x.value_raw().get_d()));
}

template<class PR> auto is_zero(Float<BoundedTag,PR> const& x) -> Logical<ValidatedTag> {
    if(x.lower_raw()>0.0 || x.upper_raw()<0.0) { return false; }
    else if(x.lower_raw()==0.0 && x.upper_raw()==0.0) { return true; }
    else { return indeterminate; }
}

template<class PR> auto is_positive(Float<BoundedTag,PR> const& x) -> Logical<ValidatedTag> {
    if(x.lower_raw()>=0.0) { return true; }
    else if(x.upper_raw()<0.0) { return false; }
    else { return indeterminate; }
}

template<class PR> auto is_positive(Float<BoundedTag,PR> const&) -> Logical<ValidatedTag>;

template<class PR> Bool same(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return x1._l==x2._l && x1._u==x2._u; }


template<class PR> OutputStream& operator<<(OutputStream& os, const Float<BoundedTag,PR>& x)
{
    //if(x.lower_raw()==x.upper_raw()) { return os << "{" << std::setprecision(Float<BoundedTag,PR>::output_precision) << x.lower_raw().get_d() << ; }
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    os << '{';
    RawFloat<PR>::set_rounding_downward();
    os << std::showpoint << std::setprecision(Float<BoundedTag,PR>::output_precision) << x.lower().get_d();
    os << ':';
    RawFloat<PR>::set_rounding_upward();
    os << std::showpoint << std::setprecision(Float<BoundedTag,PR>::output_precision) << x.upper().get_d();
    RawFloat<PR>::set_rounding_mode(rnd);
    os << '}';
    return os;

}

template<class PR> InputStream& operator>>(InputStream& is, Float<BoundedTag,PR>& x)
{
    char cl,cm,cr;
    RawFloat<PR> _l,_u;
    auto rnd=RawFloat<PR>::get_rounding_mode();
    is >> cl;
    RawFloat<PR>::set_rounding_downward();
    is >> _l;
    is >> cm;
    RawFloat<PR>::set_rounding_upward();
    is >> _u;
    is >> cr;
    RawFloat<PR>::set_rounding_mode(rnd);
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    x._l=_l; x._u=_u;
    return is;
}



template<class PR> Float<BoundedTag,PR> operator+(Float<BoundedTag,PR> const& x) { return pos(x); }
template<class PR> Float<BoundedTag,PR> operator-(Float<BoundedTag,PR> const& x) { return neg(x); }
template<class PR> Float<BoundedTag,PR> operator+(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<BoundedTag,PR> operator-(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<BoundedTag,PR> operator*(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<BoundedTag,PR> operator/(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<BoundedTag,PR>& operator+=(Float<BoundedTag,PR>& x1, Float<BoundedTag,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<BoundedTag,PR>& operator-=(Float<BoundedTag,PR>& x1, Float<BoundedTag,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<BoundedTag,PR>& operator*=(Float<BoundedTag,PR>& x1, Float<BoundedTag,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<BoundedTag,PR>& operator/=(Float<BoundedTag,PR>& x1, Float<BoundedTag,PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<BoundedTag,PR> const& x);
template<class PR> InputStream& operator>>(InputStream& is, Float<BoundedTag,PR>& x);




template<class PR> Float<MetricTag,PR> nul(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(nul(x._v),nul(x._e));
}

template<class PR> Float<MetricTag,PR> pos(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(pos(x._v),x._e);
}

template<class PR> Float<MetricTag,PR> neg(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(neg(x._v),x._e);
}

template<class PR> Float<MetricTag,PR> half(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(half(x._v),half(x._e));
}

template<class PR> Float<MetricTag,PR> sqr(Float<MetricTag,PR> const& x) {
    Float<MetricTag,PR> r=x*x;
    if(r._e>r._v) {
        r._e=half(add_up(r._e,r._v));
        r._v=r._e;
    }
    return r;
}

template<class PR> Float<MetricTag,PR> rec(Float<MetricTag,PR> const& x) {
    // Use this code to find value same as reciprocal value
    auto rv=rec_approx(x._v);
    auto ru=rec_up(sub_down(x._v,x._e));
    auto rl=rec_down(add_up(x._v,x._e));
    auto re=max(sub_up(ru,rv),sub_up(rv,rl));
    return Float<MetricTag,PR>(rv,re);
#ifdef ARIADNE_UNDEFINED
    // Use this code to get same result as interval computation
    auto ru=rec_up(sub_down(x._v,x._e));
    auto rl=rec_down(add_up(x._v,x._e));
    auto re=half(sub_up(ru,rl));
    auto rv=half(add_near(rl,ru));
    return Float<MetricTag,PR>(rv,re);
#endif
}

template<class PR> Float<MetricTag,PR> add(Float<MetricTag,PR> const& x, Float<MetricTag,PR> y) {
    auto rv=add_near(x._v,y._v);
    auto ru=add_up(x._v,y._v);
    auto rl=add_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return Float<MetricTag,PR>(rv,re);
}

template<class PR> Float<MetricTag,PR> sub(Float<MetricTag,PR> const& x, Float<MetricTag,PR> y) {
    auto rv=sub_near(x._v,y._v);
    auto ru=sub_up(x._v,y._v);
    auto rl=sub_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return Float<MetricTag,PR>(rv,re);
}

template<class PR> Float<MetricTag,PR> mul(Float<MetricTag,PR> const& x, Float<MetricTag,PR> y) {
    auto rv=mul_near(x._v,y._v);
    auto ru=mul_up(x._v,y._v);
    auto rl=mul_down(x._v,y._v);
    auto re1=add_up(half(sub_up(ru,rl)),mul_up(x._e,y._e));
    auto re2=add_up(mul_up(abs(x._v),y._e),mul_up(x._e,abs(y._v)));
    auto re=add_up(re1,re2);
    return Float<MetricTag,PR>(rv,re);
}

template<class PR> Float<MetricTag,PR> div(Float<MetricTag,PR> const& x, Float<MetricTag,PR> y) {
    return x*rec(y);
}

template<class PR> Float<MetricTag,PR> pow(Float<MetricTag,PR> const& x, Nat m) {
    return Float<MetricTag,PR>(pow(Float<BoundedTag,PR>(x),m));
}

template<class PR> Float<MetricTag,PR> pow(Float<MetricTag,PR> const& x, Int n) {
    return Float<MetricTag,PR>(pow(Float<BoundedTag,PR>(x),n));
}

template<class PR> Float<MetricTag,PR> sqrt(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(sqrt(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> exp(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(exp(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> log(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(log(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> sin(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(sin(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> cos(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(cos(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> tan(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(tan(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> asin(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(asin(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> acos(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(acos(Float<BoundedTag,PR>(x)));
}

template<class PR> Float<MetricTag,PR> atan(Float<MetricTag,PR> const& x) {
    return Float<MetricTag,PR>(atan(Float<BoundedTag,PR>(x)));
}


template<class PR> Float<MetricTag,PR> abs(Float<MetricTag,PR> const& x) {
    if(x._e<abs(x._v)) { return x; }
    else { auto rv=half(abs(x._v)+x._e); return Float<MetricTag,PR>(rv,rv); }
}

template<class PR> Float<MetricTag,PR> max(Float<MetricTag,PR> const& x1, Float<MetricTag,PR> const& x2) {
    return half((x1+x2)+abs(x1-x2));
}

template<class PR> Float<MetricTag,PR> min(Float<MetricTag,PR> const& x1, Float<MetricTag,PR> const& x2) {
    return half((x1+x2)-abs(x1-x2));
}

template<class PR> OutputStream& operator<<(OutputStream& os, Float<MetricTag,PR> const& x) {
    return os << x.value() << "\u00b1" << x.error();
}


// Mixed BoundedTag - ExactTag operations
template<class PR> Float<BoundedTag,PR> add(Float<BoundedTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return Float<BoundedTag,PR>(add_down(x1._l,x2._v),add_up(x1._u,x2._v));
}

template<class PR> Float<BoundedTag,PR> add(Float<ExactTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(add_down(x1._v,x2._l),add_down(x1._v,x2._u));
}

template<class PR> Float<BoundedTag,PR> sub(Float<BoundedTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return Float<BoundedTag,PR>(sub_down(x1._l,x2._v),sub_up(x1._u,x2._v));
}

template<class PR> Float<BoundedTag,PR> sub(Float<ExactTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(sub_down(x1._v,x2._u),sub_up(x1._v,x2._l));
}

template<class PR> Float<BoundedTag,PR> mul(Float<BoundedTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    const RawFloat<PR>& x1l=x1.lower_raw(); const RawFloat<PR>& x1u=x1.upper_raw();
    const RawFloat<PR>& x2v=x2.raw();
    RawFloat<PR> rl,ru;
    if(x2v>=0.0) {
        rl=mul_down(x1l,x2v); ru=mul_up(x1u,x2v);
    } else {
        rl=mul_down(x1u,x2v); ru=mul_up(x1l,x2v);
    }
    return Float<BoundedTag,PR>(rl,ru);
}


template<class PR> Float<BoundedTag,PR> mul(Float<ExactTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    const RawFloat<PR>& x1v=x1.raw();
    const RawFloat<PR>& x2l=x2.lower_raw(); const RawFloat<PR>& x2u=x2.upper_raw();
    RawFloat<PR> rl,ru;
    if(x1v>=0.0) {
        rl=mul_down(x1v,x2l); ru=mul_up(x1v,x2u);
    } else {
        rl=mul_down(x1v,x2u); ru=mul_up(x1v,x2l);
    }
    return Float<BoundedTag,PR>(rl,ru);
}

template<class PR> Float<BoundedTag,PR> div(Float<BoundedTag,PR> const& x1, Float<ExactTag,PR> const& x2)
{
    const RawFloat<PR>& x1l=x1.lower_raw();
    const RawFloat<PR>& x1u=x1.upper_raw();
    const RawFloat<PR>& x2v=x2.raw();
    RawFloat<PR> rl,ru;
    if(x2v>0) {
        rl=div_down(x1l,x2v); ru=div_up(x1u,x2v);
    } else if(x2v<0) {
        rl=div_down(x1u,x2v); ru=div_up(x1l,x2v);
    } else {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds const& x1, FloatValue x2)","x1="<<x1<<", x2="<<x2);
        PR pr=min(x1.precision(),x2.precision());
        rl=-RawFloat<PR>::inf(pr);
        ru=+RawFloat<PR>::inf(pr);
    }
    return Float<BoundedTag,PR>(rl,ru);
}


template<class PR> Float<BoundedTag,PR> div(Float<ExactTag,PR> const& x1, Float<BoundedTag,PR> const& x2)
{
    const RawFloat<PR>& x1v=x1.raw();
    const RawFloat<PR>& i2l=x2.lower_raw();
    const RawFloat<PR>& i2u=x2.upper_raw();
    RawFloat<PR> rl,ru;
    if(i2l<=0 && i2u>=0) {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatValue const& x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
        PR pr=min(x1.precision(),x2.precision());
        rl=-RawFloat<PR>::inf(pr);
        ru=+RawFloat<PR>::inf(pr);
    } else if(x1v>=0) {
        rl=div_down(x1v,i2u); ru=div_up(x1v,i2l);
    } else {
        rl=div_down(x1v,i2l); ru=div_up(x1v,i2u);
    }
    return Float<BoundedTag,PR>(rl,ru);
}






template<class PR> Float<ExactTag,PR> max(Float<ExactTag,PR> const& x1,  Float<ExactTag,PR> const& x2) {
    return Float<ExactTag,PR>(max(x1._v,x2._v)); }

template<class PR> Float<ExactTag,PR> min(Float<ExactTag,PR> const& x1,  Float<ExactTag,PR> const& x2) {
    return Float<ExactTag,PR>(min(x1._v,x2._v)); }

template<class PR> Float<ExactTag,PR> abs(Float<ExactTag,PR> const& x) {
    return Float<ExactTag,PR>(abs(x._v)); }

template<class PR> Float<PositiveExactTag,PR> mig(Float<ExactTag,PR> const& x) {
    return Float<PositiveExactTag,PR>(abs(x._v)); }

template<class PR> Float<PositiveExactTag,PR> mag(Float<ExactTag,PR> const& x) {
    return Float<PositiveExactTag,PR>(abs(x._v)); }


template<class PR> Float<ExactTag,PR> nul(Float<ExactTag,PR> const& x) {
    return Float<ExactTag,PR>(nul(x._v)); }

template<class PR> Float<ExactTag,PR> pos(Float<ExactTag,PR> const& x) {
    return Float<ExactTag,PR>(pos(x._v)); }

template<class PR> Float<ExactTag,PR> neg(Float<ExactTag,PR> const& x) {
    return Float<ExactTag,PR>(neg(x._v)); }

template<class PR> Float<ExactTag,PR> half(Float<ExactTag,PR> const& x) {
    return Float<ExactTag,PR>(half(x._v)); }

template<class PR> Float<BoundedTag,PR> sqr(Float<ExactTag,PR> const& x) {
    return Float<BoundedTag,PR>(mul_down(x._v,x._v),mul_up(x._v,x._v)); }

template<class PR> Float<BoundedTag,PR> rec(Float<ExactTag,PR> const& x) {
    return Float<BoundedTag,PR>(rec_down(x._v),rec_up(x._v)); }

template<class PR> Float<BoundedTag,PR> add(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return Float<BoundedTag,PR>(add_down(x1._v,x2._v),add_up(x1._v,x2._v)); }

template<class PR> Float<BoundedTag,PR> sub(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return Float<BoundedTag,PR>(sub_down(x1._v,x2._v),sub_up(x1._v,x2._v)); }

template<class PR> Float<BoundedTag,PR> mul(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return Float<BoundedTag,PR>(mul_down(x1._v,x2._v),mul_up(x1._v,x2._v)); }

template<class PR> Float<BoundedTag,PR> div(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return Float<BoundedTag,PR>(div_down(x1._v,x2._v),div_up(x1._v,x2._v)); }

template<class PR> Float<BoundedTag,PR> pow(Float<ExactTag,PR> const& x, Nat m) {
    return pow(Float<BoundedTag,PR>(x),m); }

template<class PR> Float<BoundedTag,PR> pow(Float<ExactTag,PR> const& x, Int n) {
    return pow(Float<BoundedTag,PR>(x),n); }

template<class PR> Float<BoundedTag,PR> med(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return add(half(x1),half(x2)); }

template<class PR> Float<BoundedTag,PR> rad(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return sub(half(x2),half(x1)); }

template<class PR> Float<BoundedTag,PR> sqrt(Float<ExactTag,PR> const& x) {
    return sqrt(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> exp(Float<ExactTag,PR> const& x) {
    return exp(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> log(Float<ExactTag,PR> const& x) {
    return log(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> sin(Float<ExactTag,PR> const& x) {
    return sin(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> cos(Float<ExactTag,PR> const& x) {
    return cos(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> tan(Float<ExactTag,PR> const& x) {
    return tan(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> asin(Float<ExactTag,PR> const& x) {
    return asin(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> acos(Float<ExactTag,PR> const& x) {
    return acos(Float<BoundedTag,PR>(x)); }

template<class PR> Float<BoundedTag,PR> atan(Float<ExactTag,PR> const& x) {
    return atan(Float<BoundedTag,PR>(x)); }

template<class PR> Logical<ExactTag> eq(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return x1._v == x2._v; }

template<class PR> Logical<ExactTag> leq(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return x1._v <= x2._v; }

template<class PR> Bool same(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return x1._v==x2._v; }

template<class PR> Float<ExactTag,PR> operator+(Float<ExactTag,PR> const& x) { return pos(x); }
template<class PR> Float<ExactTag,PR> operator-(Float<ExactTag,PR> const& x) { return neg(x); }
template<class PR> Float<BoundedTag,PR> operator+(Float<ExactTag,PR> const& x1,  Float<ExactTag,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<BoundedTag,PR> operator-(Float<ExactTag,PR> const& x1,  Float<ExactTag,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<BoundedTag,PR> operator*(Float<ExactTag,PR> const& x1,  Float<ExactTag,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<BoundedTag,PR> operator/(Float<ExactTag,PR> const& x1,  Float<ExactTag,PR> const& x2) { return div(x1,x2); }

template<class PR> Float<ExactTag,PR> operator*(Float<ExactTag,PR> const& x, TwoExp y) {
    Float<ExactTag,PR> yv=y; return Float<ExactTag,PR>(x.raw()*yv.raw()); }
template<class PR> Float<ExactTag,PR> operator/(Float<ExactTag,PR> const& x, TwoExp y) {
    Float<ExactTag,PR> yv=y; return Float<ExactTag,PR>(x.raw()/yv.raw()); }
template<class PR> Float<ExactTag,PR>& operator*=(Float<ExactTag,PR>& x, TwoExp y) {
    Float<ExactTag,PR> yv=y; return x=Float<ExactTag,PR>(x.raw()*yv.raw()); }
template<class PR> Float<ExactTag,PR>& operator/=(Float<ExactTag,PR>& x, TwoExp y) {
    Float<ExactTag,PR> yv=y; return x=Float<ExactTag,PR>(x.raw()/yv.raw()); }

template<class PR> Boolean operator==(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) { return x1.raw()==x2.raw(); }
template<class PR> Boolean operator!=(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) { return x1.raw()!=x2.raw(); }
template<class PR> Boolean operator<=(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) { return x1.raw()<=x2.raw(); }
template<class PR> Boolean operator>=(Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) { return x1.raw()>=x2.raw(); }
template<class PR> Boolean operator< (Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) { return x1.raw()< x2.raw(); }
template<class PR> Boolean operator> (Float<ExactTag,PR> const& x1, Float<ExactTag,PR> const& x2) { return x1.raw()> x2.raw(); }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<ExactTag,PR> const& x) {
    return os << std::setprecision(Float<ExactTag,PR>::output_precision) << x.raw();
}

template<class PR> InputStream& operator>>(InputStream& is, Float<ExactTag,PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
    auto v = nul(x._v);
    is >> v;
    ARIADNE_ASSERT(is);
    x._v=v;
    return is;
}

template<class PR> Integer integer_cast(Float<ExactTag,PR> const& x) {
    Integer z=static_cast<int>(x._v.get_d());
    ARIADNE_ASSERT(z==x);
    return std::move(z);
}


Float<ExactTag,Precision64> cast_exact(Real const& x) {
    return cast_exact(Float<ApproximateTag,Precision64>(x));
}

template<class PR> Bool operator==(Float<ExactTag,PR> const& x, const Rational& q) { return Rational(x)==q; }
template<class PR> Bool operator!=(Float<ExactTag,PR> const& x, const Rational& q) { return Rational(x)!=q; }
template<class PR> Bool operator<=(Float<ExactTag,PR> const& x, const Rational& q) { return Rational(x)<=q; }
template<class PR> Bool operator>=(Float<ExactTag,PR> const& x, const Rational& q) { return Rational(x)>=q; }
template<class PR> Bool operator< (Float<ExactTag,PR> const& x, const Rational& q) { return Rational(x)< q; }
template<class PR> Bool operator> (Float<ExactTag,PR> const& x, const Rational& q) { return Rational(x)> q; }

template<class PR> Bool operator==(const Rational& q, Float<ExactTag,PR> const& x) { return q==Rational(x); }
template<class PR> Bool operator!=(const Rational& q, Float<ExactTag,PR> const& x) { return q!=Rational(x); }
template<class PR> Bool operator<=(const Rational& q, Float<ExactTag,PR> const& x) { return q<=Rational(x); }
template<class PR> Bool operator>=(const Rational& q, Float<ExactTag,PR> const& x) { return q>=Rational(x); }
template<class PR> Bool operator< (const Rational& q, Float<ExactTag,PR> const& x) { return q< Rational(x); }
template<class PR> Bool operator> (const Rational& q, Float<ExactTag,PR> const& x) { return q> Rational(x); }



template<class PR> Bool refines(Float<MetricTag,PR> const& x1, Float<MetricTag,PR> const& x2) {
    return (x1._v>=x2._v ? sub_up(x1._v,x2._v) : sub_up(x2._v,x1._v)) <= sub_down(x2._e, x1._e);
}

template<class PR> Bool models(Float<BoundedTag,PR> const& x1, Float<ExactTag,PR> const& x2) {
    return x1._l<=x2._v && x1._u >= x2._v;
}

template<class PR> Bool refines(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return x1._l>=x2._l && x1._u <= x2._u;
}

template<class PR> Float<BoundedTag,PR> refinement(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return Float<BoundedTag,PR>(max(x1._l,x2._l),min(x1._u,x2._u));
}

template<class PR> Bool consistent(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return x1._l<=x2._u && x1._u >= x2._l;
}

template<class PR> Bool inconsistent(Float<BoundedTag,PR> const& x1, Float<BoundedTag,PR> const& x2) {
    return x1._l>x2._u || x1._u < x2._l;
}


template<class PR> Float<MetricTag,PR> refinement(Float<MetricTag,PR> const& x1, Float<MetricTag,PR> const& x2) {
    return Float<MetricTag,PR>(refinement(Float<BoundedTag,PR>(x1),Float<BoundedTag,PR>(x2)));
}

template<class PR> Bool refines(Float<LowerTag,PR> const& x1, Float<LowerTag,PR> const& x2) {
    return x1._l>=x2._l;
}

template<class PR> Bool refines(Float<UpperTag,PR> const& x1, Float<UpperTag,PR> const& x2) {
    return x1._u <= x2._u;
}


#ifdef ARIADNE_ENABLE_SERIALIZATION
template<class PR, class A> Void serialize(A& _a, Float<BoundedTag,PR>& x, const Nat version) {
    _a & x.lower_raw() & x.upper_raw(); }
#endif




template<class PR> Bool same(Float<PositiveUpperTag,PR> const& x1, Float<PositiveUpperTag,PR> const& x2) {
    return x1._u == x2._u;
}

template<class PR> Float<PositiveUpperTag,PR> max(Float<PositiveUpperTag,PR> const& x1, Float<PositiveUpperTag,PR> const& x2) {
    return Float<PositiveUpperTag,PR>(max(x1._u,x2._u));
}

template<class PR> Float<PositiveUpperTag,PR> min(Float<PositiveUpperTag,PR> const& x1, Float<PositiveUpperTag,PR> const& x2) {
    return Float<PositiveUpperTag,PR>(min(x1._u,x2._u));
}

template<class PR> Float<PositiveUpperTag,PR> abs(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(x._u);
}

template<class PR> Float<PositiveApproximateTag,PR> mig(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveApproximateTag,PR>(x._u);
}

template<class PR> Float<PositiveUpperTag,PR> mag(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(x._u);
}

template<class PR> Float<PositiveUpperTag,PR> pos(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(pos(x._u));
}

template<class PR> Float<LowerTag,PR> neg(Float<PositiveUpperTag,PR> const& x) {
    return Float<LowerTag,PR>(neg(x._u));
}

template<class PR> Float<PositiveUpperTag,PR> rec(Float<PositiveLowerTag,PR> const& x) {
    return Float<PositiveLowerTag,PR>(rec_up(x._l));
}

template<class PR> Float<PositiveUpperTag,PR> half(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(half(x._u));
}

template<class PR> Float<PositiveUpperTag,PR> sqr(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveUpperTag,PR>(mul_up(x._u,x._u));
}

template<class PR> Float<PositiveUpperTag,PR> add(Float<PositiveUpperTag,PR> const& x1, Float<PositiveUpperTag,PR> const& x2) {
    return Float<PositiveUpperTag,PR>(add_up(x1._u,x2._u));
}

template<class PR> Float<PositiveUpperTag,PR> mul(Float<PositiveUpperTag,PR> const& x1, Float<PositiveUpperTag,PR> const& x2) {
    return Float<PositiveUpperTag,PR>(mul_up(x1._u,x2._u));
}

template<class PR> Float<PositiveUpperTag,PR> div(Float<PositiveUpperTag,PR> const& x1, Float<PositiveLowerTag,PR> const& x2) {
    return Float<PositiveUpperTag,PR>(div_up(x1._u,x2._l));
}

template<class PR> Float<PositiveUpperTag,PR> pow(Float<PositiveUpperTag,PR> const& x, Nat m) {
    return Float<PositiveUpperTag,PR>(pow_up(x._u,m));
}

template<class PR> Float<PositiveApproximateTag,PR> pow(Float<PositiveUpperTag,PR> const& x, Int n) {
    return Float<PositiveApproximateTag,PR>(pow_approx(x._u,n));
}

template<class PR> Float<UpperTag,PR> log(Float<PositiveUpperTag,PR> const& x) {
    return Float<UpperTag,PR>(log_up(x._u));
}

template<class PR> OutputStream& operator<<(OutputStream& os, Float<PositiveUpperTag,PR> const& x) {
    return os << static_cast<Float<UpperTag,PR>const&>(x);
}


template<class PR> Float<PositiveLowerTag,PR> mig(Float<PositiveLowerTag,PR> const& x) {
    return Float<PositiveLowerTag,PR>(abs(x._l)); }

template<class PR> Float<PositiveLowerTag,PR> rec(Float<PositiveUpperTag,PR> const& x) {
    return Float<PositiveLowerTag,PR>(rec_down(x._u));
}


template<class PR> Float<PositiveApproximateTag,PR> mul(Float<PositiveApproximateTag,PR> const& x1, Float<PositiveApproximateTag,PR> const& x2) {
    return Float<PositiveApproximateTag,PR>(mul_near(x1._a,x2._a));
}





template<> Nat integer_cast<Nat,Float64Approximation>(Float64Approximation const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,Float64Approximation>(Float64Approximation const& x) {
    return std::round(x.get_d()); }

template<> Nat integer_cast<Nat,Float64LowerBound>(Float64LowerBound const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,Float64LowerBound>(Float64LowerBound const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,Float64Bounds>(Float64Bounds const& x) {
    return std::round((x.lower().get_d()+x.upper().get_d())/2); }



template<class PR, class P> Float<P,PR> max(Float<P,PR> const& x1, Float<P,PR> const& x2) { return max<PR>(x1,x2); }
template<class PR, class P> Float<P,PR> min(Float<P,PR> const& x1, Float<P,PR> const& x2) { return min<PR>(x1,x2); }
template<class PR, class P> Float<Weaker<P,Negated<P>>,PR> abs(Float<P,PR> const& x) { return abs<PR>(x); }
template<class PR, class P> Float<Unsigned<Weaker<P,UpperTag>>,PR> mag(Float<P,PR> const& x) { return mag<PR>(x); }
template<class PR, class P> Float<Unsigned<Weaker<P,LowerTag>>,PR> mig(Float<P,PR> const& x) { return mig<PR>(x); }

template<class PR, class P> Float<P,PR> round(Float<P,PR> const& x) { return round<PR>(x); }

template<class PR, class P> Float<P,PR> nul(Float<P,PR> const& x) { return nul<PR>(x); }
template<class PR, class P> Float<P,PR> pos(Float<P,PR> const& x) { return pos<PR>(x); }
template<class PR, class P> Float<Negated<P>,PR> neg(Float<P,PR> const& x) { return neg<PR>(x); }
template<class PR, class P> Float<P,PR> half(Float<P,PR> const& x) { return half<PR>(x); }
template<class PR, class P> Float<Widen<P>,PR> sqr(Float<P,PR> const& x) { return sqr<PR>(x); }
template<class PR, class P> Float<Widen<Inverted<P>>,PR> rec(Float<P,PR> const& x) { return rec<PR>(x); }

template<class PR, class P> Float<Widen<P>,PR> add(Float<P,PR> const& x1, Float<P,PR> const& x2) { return add<PR>(x1,x2); }
template<class PR, class P> Float<Widen<P>,PR> sub(Float<P,PR> const& x1, Float<Negated<P>,PR> const& x2) { return sub<PR>(x1,x2); }
template<class PR, class P> Float<Widen<P>,PR> mul(Float<P,PR> const& x1, Float<P,PR> const& x2) { return mul<PR>(x1,x2); }
template<class PR, class P> Float<Widen<P>,PR> div(Float<P,PR> const& x1, Float<Inverted<P>,PR> const& x2) { return div<PR>(x1,x2); }

template<class PR, class P> Float<Widen<P>,PR> pow(Float<P,PR> const& x, Nat m) { return pow<PR>(x,m); }
template<class PR, class P> Float<Widen<Unorder<P>>,PR> pow(Float<P,PR> const& x, Int n) { return pow<PR>(x,n); }

template<class PR, class P> Float<Widen<P>,PR> sqrt(Float<P,PR> const& x) { return sqrt<PR>(x); }
template<class PR, class P> Float<Widen<P>,PR> exp(Float<P,PR> const& x) { return exp<PR>(x); }
template<class PR, class P> Float<Widen<Signed<P>>,PR> log(Float<P,PR> const& x) { return log<PR>(x); }
template<class PR, class P> Float<Widen<Unorder<P>>,PR> sin(Float<P,PR> const& x) { return sin<PR>(x); }
template<class PR, class P> Float<Widen<Unorder<P>>,PR> cos(Float<P,PR> const& x) { return cos<PR>(x); }
template<class PR, class P> Float<Widen<Unorder<P>>,PR> tan(Float<P,PR> const& x) { return tan<PR>(x); }
template<class PR, class P> Float<Widen<Unorder<P>>,PR> asin(Float<P,PR> const& x) { return asin<PR>(x); }
template<class PR, class P> Float<Widen<Unorder<P>>,PR> acos(Float<P,PR> const& x) { return acos<PR>(x); }
template<class PR, class P> Float<Widen<P>,PR> atan(Float<P,PR> const& x) { return atan<PR>(x); }

template<class PR, class P> Float<P,PR> operator*(Float<P,PR> const& x1, TwoExp y2) { return operator* <PR>(x1,y2); }
template<class PR, class P> Float<P,PR> operator/(Float<P,PR> const& x1, TwoExp y2) { return operator/ <PR>(x1,y2); }

template<class PR, class P> Integer integer_cast(Float<P,PR> const& x) { return integer_cast<PR>(x); }

template<class PR, class P> Bool same(Float<P,PR> const& x1, Float<P,PR> const& x2) { return same<PR>(x1,x2); }

template<class PR, class P> FloatEqualsType<PR,P,Negated<P>> eq(Float<P,PR> const& x1, Float<Negated<P>,PR> const& x2) { return eq<PR>(x1,x2); }
template<class PR, class P> FloatLessType<PR,P,Negated<P>> leq(Float<P,PR> const& x1, Float<Negated<P>,PR> const& x2) { return leq<PR>(x1,x2); }

template<class PR, class P> OutputStream& operator<<(OutputStream& os, Float<P,PR> const& x) { return operator<<(os,x); }

template<class PR, class P> InputStream& operator>>(InputStream& is, Float<P,PR>& x) { return operator>>(is,x); }

template<class PR> Float<ApproximateTag,PR> make_float(Number<ApproximateTag> x) { return Float<ApproximateTag,PR>(x); }
template<class PR> Float<LowerTag,PR> make_float(Number<ValidatedLowerTag> x) { return Float<LowerTag,PR>(x); }
template<class PR> Float<UpperTag,PR> make_float(Number<ValidatedUpperTag> x) { return Float<UpperTag,PR>(x); }
template<class PR> Float<BoundedTag,PR> make_float(Number<ValidatedTag> x) { return Float<BoundedTag,PR>(x); }
template<class PR> Float<BoundedTag,PR> make_float(Number<EffectiveTag> x) { return Float<BoundedTag,PR>(x); }
template<class PR> Float<BoundedTag,PR> make_float(Number<ExactTag> x) { return Float<BoundedTag,PR>(x); }
template<class PR> Float<BoundedTag,PR> make_float(Real r) { return Float<BoundedTag,PR>(r); }
template<class PR> Float<BoundedTag,PR> make_float(Rational q) { return Float<BoundedTag,PR>(q); }
template<class PR> Float<ExactTag,PR> make_float(Integer z) { return Float<ExactTag,PR>(z); }

template class Float<ApproximateTag,Precision64>;
template class Float<LowerTag,Precision64>;
template class Float<UpperTag,Precision64>;
template class Float<BoundedTag,Precision64>;
template class Float<MetricTag,Precision64>;
template class Float<ExactTag,Precision64>;

template class Float<ApproximateTag,PrecisionMP>;
template class Float<LowerTag,PrecisionMP>;
template class Float<UpperTag,PrecisionMP>;
template class Float<BoundedTag,PrecisionMP>;
template class Float<MetricTag,PrecisionMP>;
template class Float<ExactTag,PrecisionMP>;


template<class T> class Dynamic { public: virtual ~Dynamic() = default; };
template<class T> String mangled_class_name() { Dynamic<T>* p; return typeid(p).name(); }

template<class P, class PR> std::size_t instantiate_float() {
    using std::size_t;
    typedef OutputStream OS;
    typedef OutputStream IS;
    typedef Float<P,PR> X;
    typedef Float<P,PR> const& XCRef;

    typedef X const& Xcr;
    typedef Float<Negated<P>,PR> const& NEGXcr;
    typedef Float<Inverted<P>,PR> const& INVXcr;

    typedef decltype(neg(declval<X>())) NEGX;
    typedef decltype(abs(declval<X>())) ABSX;
    typedef decltype(mig(declval<X>())) MIGX;
    typedef decltype(mag(declval<X>())) MAGX;
    typedef decltype(add(declval<X>(),declval<X>())) WX;
    typedef decltype(sub(declval<X>(),declval<X>())) WUX;

    typedef decltype(pow(declval<X>(),declval<Nat>())) POWNX;
    typedef decltype(pow(declval<X>(),declval<Int>())) POWZX;

    typedef decltype(exp(declval<X>())) PX;
    typedef decltype(log(declval<X>())) SX;

    typedef decltype(eq(declval<X>(),declval<NEGX>())) XE;
    typedef decltype(leq(declval<X>(),declval<NEGX>())) XL;

    return (size_t)(X(*)(Xcr,Xcr))&max<PR,P>
         + (size_t)(X(*)(Xcr,Xcr))&min<PR,P>
         + (size_t)(ABSX(*)(Xcr))&abs<PR,P>
         + (size_t)(MIGX(*)(Xcr))&mig<PR,P>
         + (size_t)(MAGX(*)(Xcr))&mag<PR,P>

         + (size_t)(WX(*)(Xcr,Xcr))&add<PR,P>
         + (size_t)(WX(*)(Xcr,NEGXcr))&sub<PR,P>
         + (size_t)(WX(*)(Xcr,Xcr))&mul<PR,P>
         + (size_t)(WX(*)(Xcr,INVXcr))&div<PR,P>

         + (size_t)(X(*)(Xcr))&nul<PR,P>
         + (size_t)(X(*)(Xcr))&pos<PR,P>
         + (size_t)(NEGX(*)(Xcr))&neg<PR,P>
         + (size_t)(WX(*)(Xcr))&sqr<PR,P>
         + (size_t)(WX(*)(INVXcr))&rec<PR,P>
         + (size_t)(X(*)(Xcr))&half<PR,P>
         + (size_t)(POWNX(*)(Xcr,Nat))&pow<PR,P>
         + (size_t)(POWZX(*)(Xcr,Int))&pow<PR,P>

         + (size_t)(WX(*)(Xcr))&sqrt<PR,P>
         + (size_t)(PX(*)(Xcr))&exp<PR,P>
         + (size_t)(SX(*)(Xcr))&log<PR,P>
         + (size_t)(WUX(*)(Xcr))&sin<PR,P>
         + (size_t)(WUX(*)(Xcr))&cos<PR,P>
         + (size_t)(WUX(*)(Xcr))&tan<PR,P>
         + (size_t)(WUX(*)(Xcr))&asin<PR,P>
         + (size_t)(WUX(*)(Xcr))&acos<PR,P>
         + (size_t)(WX(*)(Xcr))&atan<PR,P>

         + (size_t)(XE(*)(Xcr,NEGXcr))&eq<PR,P>
         + (size_t)(XL(*)(Xcr,NEGXcr))&leq<PR,P>

         + (size_t)(OS&(*)(OS&,Xcr)) &operator<< <PR,P
         > + (size_t)(IS&(*)(IS&,X&)) &operator>> <PR,P>

         + (size_t)(Bool(*)(Xcr,Xcr)) &same<PR,P>
         + (size_t)(Integer(*)(Xcr)) &integer_cast<PR,P>
        ;
}

template<class PR> std::size_t instantiate_floats() {
    return instantiate_float<ApproximateTag,PR>()
        + instantiate_float<LowerTag,PR>()
        + instantiate_float<UpperTag,PR>()
        + instantiate_float<BoundedTag,PR>()
        + instantiate_float<MetricTag,PR>()
        + instantiate_float<ExactTag,PR>()
        + instantiate_float<PositiveUpperTag,PR>()
        ;
}

template std::size_t instantiate_floats<Precision64>();
template std::size_t instantiate_floats<PrecisionMP>();

template Float64Value operator* <Precision64,ExactTag>(Float64Value const&, TwoExp);
template Float64Value operator/ <Precision64,ExactTag>(Float64Value const&, TwoExp);

template PositiveFloat64UpperBound abs(PositiveFloat64UpperBound const&);
template PositiveFloatMPUpperBound abs(PositiveFloatMPUpperBound const&);

template Float64Bounds round<Precision64,BoundedTag>(Float64Bounds const&);
template Float64Approximation round<Precision64,ApproximateTag>(Float64Approximation const&);

template Bool models(Float64Bounds const&, Float64Value const&);
template Bool consistent(Float64Bounds const&, Float64Bounds const&);
template Bool inconsistent(Float64Bounds const&, Float64Bounds const&);
template Bool refines(Float64Bounds const&, Float64Bounds const&);
template Float64Bounds refinement(Float64Bounds const&, Float64Bounds const&);
template Bool same(Float64Bounds const&, Float64Bounds const&);

template Bool refines(Float64Ball const&, Float64Ball const&);
template Float64Ball refinement(Float64Ball const&, Float64Ball const&);

template Bool refines(Float64LowerBound const&, Float64LowerBound const&);
template Bool refines(Float64UpperBound const&, Float64UpperBound const&);

template Bool refines(FloatMPBall const&, FloatMPBall const&);
template Bool refines(FloatMPBounds const&, FloatMPBounds const&);
template Bool refines(FloatMPLowerBound const&, FloatMPLowerBound const&);
template Bool refines(FloatMPUpperBound const&, FloatMPUpperBound const&);

Float64Value midpoint(Float64Bounds const& x) { return x.value(); }

template Float<Widen<PositiveApproximateTag>,Precision64> mul<Precision64,PositiveApproximateTag>(Float<PositiveApproximateTag,Precision64> const& x1, Float<PositiveApproximateTag,Precision64> const& x2);

template<> String class_name<Float64Approximation>() { return "Float64Approximation"; }
template<> String class_name<Float64LowerBound>() { return "Float64LowerBound"; }
template<> String class_name<Float64UpperBound>() { return "Float64UpperBound"; }
template<> String class_name<Float64Bounds>() { return "Float64Bounds"; }
template<> String class_name<Float64Ball>() { return "Float64Ball"; }
template<> String class_name<Float64Value>() { return "Float64Value"; }
template<> String class_name<FloatMPApproximation>() { return "FloatMPApproximation"; }
template<> String class_name<FloatMPLowerBound>() { return "FloatMPLowerBound"; }
template<> String class_name<FloatMPUpperBound>() { return "FloatMPUpperBound"; }
template<> String class_name<FloatMPBounds>() { return "FloatMPBounds"; }
template<> String class_name<FloatMPBall>() { return "FloatMPBall"; }
template<> String class_name<FloatMPValue>() { return "FloatMPValue"; }

} // namespace Ariadne
