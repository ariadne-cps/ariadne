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


template<class PR> Nat FloatError<PR>::output_places = 3;
template<class PR> Nat FloatApproximation<PR>::output_places = 4;
template<class PR> Nat FloatBounds<PR>::output_places=8;
template<class PR> Nat FloatValue<PR>::output_places = 16;

const FloatValue<Precision64> infty = FloatValue<Precision64>(Float64::inf());

FloatError<Precision64> operator"" _error(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatError<Precision64>(Float64(x));
}

FloatValue<Precision64> operator"" _exact(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatValue<Precision64>(x);
}

FloatBall<Precision64> operator"" _near(long double lx) {
    volatile double x=lx;
    volatile long double le=std::abs((long double)x-lx);
    volatile double e=le;
    while(e<le) { e*=(1+std::numeric_limits<double>::epsilon()); }
    return FloatBall<Precision64>(x,e);
}

FloatUpperBound<Precision64> operator"" _upper(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x<lx) { x+=min; }
    while (x<lx) { x+=std::abs(x)*eps; }
    return FloatUpperBound<Precision64>(x);
}

FloatLowerBound<Precision64> operator"" _lower(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x>lx) { x-=min; }
    while (x>lx) { x-=std::abs(x)*eps; }
    return FloatLowerBound<Precision64>(x);
}

FloatApproximation<Precision64> operator"" _approx(long double lx) {
    double x=lx;
    return FloatApproximation<Precision64>(x);
}


TwoExp::operator FloatValue<Precision64> () const {
    return FloatValue<Precision64>(this->get_d());
}


/*
template<class PR> FloatValue<PR>::FloatValue(Rational const& q, PR pr)
    : _v(0.0,pr)
{
    FloatBounds<PR> x(q,pr);
    if(x.lower_raw()==x.upper_raw()) { this->_v==x.value_raw(); }
    else { ARIADNE_THROW(std::runtime_error,"FloatValue(Rational q, Precision pr)","q="<<q<<" cannot be expressed exactly to precision \n"); }
}
*/

template<class PR> FloatValue<PR>::operator Rational() const {
    return Rational(this->get_d());
}

template<class PR> FloatValue<PR>::FloatValue(TwoExp const& t, PR pr)
    : _v(pow(RawFloat<PR>(2.0,pr),t.exponent()))
{
}

template<class PR> FloatValue<PR>::FloatValue(Integer const& z, PR pr)
    : _v(Rational(z),RawFloat<PR>::to_nearest,pr)
{
    Rational q(_v);
    ARIADNE_PRECONDITION(z==q);
}

template<class PR> FloatValue<PR>::operator ExactNumber() const {
    return ExactNumber(new NumberWrapper<FloatValue<PR>>(*this));
}

template<class PR> FloatBall<PR>::FloatBall(Rational const& q, PR pr)
    : _v(RawFloat<PR>(q,RawFloat<PR>::to_nearest,pr)), _e(abs(Rational(_v)-q),RawFloat<PR>::upward,pr) {
}

template<class PR> FloatBall<PR>::FloatBall(Real const& r, PR pr)
    : FloatBall(r.get(pr)) {
}

template<class PR> FloatBall<PR>::FloatBall(ValidatedNumber const& x, PR pr)
    : FloatBall(x.get(MetricTag(),pr)) {
}

template<class PR> FloatBall<PR>::operator ValidatedNumber() const {
    return ValidatedNumber(new NumberWrapper<FloatBall<PR>>(*this));
}


template<class PR> FloatBounds<PR>::FloatBounds(Integer const& z, PR pr)
    : FloatBounds(RawFloat<PR>(z,RawFloat<PR>::downward,pr),RawFloat<PR>(z,RawFloat<PR>::upward,pr)) {
}

template<class PR> FloatBounds<PR>::FloatBounds(Rational const& q, PR pr)
    : FloatBounds(RawFloat<PR>(q,RawFloat<PR>::downward,pr),RawFloat<PR>(q,RawFloat<PR>::upward,pr)) {
}

template<class PR> FloatBounds<PR>::FloatBounds(Rational const& ql, Rational const& qu, PR pr)
    : FloatBounds(RawFloat<PR>(ql,RawFloat<PR>::downward,pr),RawFloat<PR>(qu,RawFloat<PR>::upward,pr)) {
}

template<class PR> FloatBounds<PR>::FloatBounds(Real const& x, PR pr)
    : FloatBounds(x.get(pr)) {
}

template<class PR> FloatBounds<PR>::FloatBounds(Number<ValidatedTag> const& x, PR pr)
    : FloatBounds(x.get(BoundedTag(),pr)) {
}

template<class PR> FloatBounds<PR>::operator Number<ValidatedTag>() const {
    return Number<ValidatedTag>(new NumberWrapper<FloatBounds<PR>>(*this));
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(Rational const& q, PR pr)
    : FloatUpperBound(FloatBounds<PR>(q,pr)) {
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(Real const& r, PR pr)
    : FloatUpperBound(r.get(pr)) {
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(ValidatedUpperNumber const& x, PR pr)
    : FloatUpperBound(x.get(UpperTag(),pr)) {
}

template<class PR> FloatUpperBound<PR>::operator ValidatedUpperNumber() const {
    //return Number<ValidatedUpperTag>(new NumberWrapper<FloatUpperBound<PR>>(*this));
}

template<class PR> FloatLowerBound<PR>::FloatLowerBound(Rational const& q, PR pr)
    : FloatLowerBound(FloatBounds<PR>(q,pr)) {
}

template<class PR> FloatLowerBound<PR>::FloatLowerBound(Real const& r, PR pr)
    : FloatLowerBound(r.get(pr)) {
}

template<class PR> FloatLowerBound<PR>::FloatLowerBound(ValidatedLowerNumber const& x, PR pr)
    : FloatLowerBound(x.get(LowerTag(),pr)) {
}

template<class PR> FloatLowerBound<PR>::operator ValidatedLowerNumber() const {
    // return ValidatedLowerNumber(new NumberWrapper<FloatLowerBound<PR>>(*this));
}

template<class PR> FloatApproximation<PR>::FloatApproximation(Rational const& q, PR pr)
    : FloatApproximation(FloatBounds<PR>(q,pr)) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(Real const& r, PR pr)
    : FloatApproximation(r.get(pr)) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(Number<ApproximateTag> const& x, PR pr)
    : FloatApproximation(x.get(ApproximateTag(),pr)) {
}

template<class PR> FloatApproximation<PR>::operator Number<ApproximateTag>() const {
    return Number<ApproximateTag>(new NumberWrapper<FloatApproximation<PR>>(*this));
}



template<class PR> FloatValue<PR>::FloatValue(Integer const& z) : FloatValue(z,RawFloat<PR>::get_default_precision()) { }

template<class PR> FloatBall<PR>::FloatBall(Integer const& z) : FloatBall(Rational(z)) { }
template<class PR> FloatBall<PR>::FloatBall(Rational const& q) : FloatBall(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> FloatBall<PR>::FloatBall(Number<ValidatedTag> const& x) : FloatBall(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> FloatBounds<PR>::FloatBounds(Integer const& z) : FloatBounds(Rational(z)) { }
template<class PR> FloatBounds<PR>::FloatBounds(Dyadic const& b) : FloatBounds(Rational(b)) { }
template<class PR> FloatBounds<PR>::FloatBounds(Decimal const& d) : FloatBounds(Rational(d)) { }
template<class PR> FloatBounds<PR>::FloatBounds(Rational const& q) : FloatBounds(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> FloatBounds<PR>::FloatBounds(Real const& x) : FloatBounds(x,RawFloat<PR>::get_default_precision()) { }
template<class PR> FloatBounds<PR>::FloatBounds(Number<ValidatedTag> const& x) : FloatBounds(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> FloatUpperBound<PR>::FloatUpperBound(Integer const& z) : FloatUpperBound(Rational(z)) { }
template<class PR> FloatUpperBound<PR>::FloatUpperBound(Rational const& q) : FloatUpperBound(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> FloatUpperBound<PR>::FloatUpperBound(Number<ValidatedUpperTag> const& x) : FloatUpperBound(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> FloatLowerBound<PR>::FloatLowerBound(Integer const& z) : FloatLowerBound(Rational(z)) { }
template<class PR> FloatLowerBound<PR>::FloatLowerBound(Rational const& q) : FloatLowerBound(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> FloatLowerBound<PR>::FloatLowerBound(Number<ValidatedLowerTag> const& x) : FloatLowerBound(x,RawFloat<PR>::get_default_precision()) { }

template<class PR> FloatApproximation<PR>::FloatApproximation(Integer const& z) : FloatApproximation(Rational(z)) { }
template<class PR> FloatApproximation<PR>::FloatApproximation(Dyadic const& b) : FloatApproximation(Rational(b)) { }
template<class PR> FloatApproximation<PR>::FloatApproximation(Decimal const& d) : FloatApproximation(Rational(d)) { }
template<class PR> FloatApproximation<PR>::FloatApproximation(Rational const& q) : FloatApproximation(q,RawFloat<PR>::get_default_precision()) { }
template<class PR> FloatApproximation<PR>::FloatApproximation(Number<ApproximateTag> const& x) : FloatApproximation(x,RawFloat<PR>::get_default_precision()) { }



template<class PR> FloatApproximation<PR>::FloatApproximation(FloatLowerBound<PR> const& x) : _a(x.raw()) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(FloatUpperBound<PR> const& x) : _a(x.raw()) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(FloatBounds<PR> const& x) : _a(x.value_raw()) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(FloatBall<PR> const& x) : _a(x.value_raw()) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(FloatValue<PR> const& x) : _a(x.raw()) {
}

template<class PR> FloatApproximation<PR>::FloatApproximation(FloatError<PR> const& x) : _a(x.raw()) {
}

template<class PR> FloatLowerBound<PR>::FloatLowerBound(FloatBounds<PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> FloatLowerBound<PR>::FloatLowerBound(FloatBall<PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> FloatLowerBound<PR>::FloatLowerBound(FloatValue<PR> const& x) : _l(x.raw()) {
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(FloatBounds<PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(FloatBall<PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(FloatValue<PR> const& x) : _u(x.raw()) {
}

template<class PR> FloatUpperBound<PR>::FloatUpperBound(FloatError<PR> const& x) : _u(x.raw()) {
}

template<class PR> FloatBounds<PR>::FloatBounds(FloatBall<PR> const& x) : _l(x.lower_raw()), _u(x.upper_raw()) {
}

template<class PR> FloatBounds<PR>::FloatBounds(FloatValue<PR> const& x) : _l(x.raw()), _u(x.raw()) {
}

template<class PR> FloatBall<PR>::FloatBall(FloatBounds<PR> const& x) : _v(x.value_raw()), _e(x.error_raw()) {
}

template<class PR> FloatBall<PR>::FloatBall(FloatValue<PR> const& x) : _v(x.raw()), _e(nul(x.raw())) {
}


template<class PR> FloatBall<PR> FloatValue<PR>::pm(FloatError<PR> _e) const {
    FloatValue<PR> const& _v=*this; return FloatBall<PR>(_v,_e);
}




template<class PR> FloatApproximation<PR> _floor(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(floor(x._a)); }
template<class PR> FloatApproximation<PR> _ceil(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(ceil(x._a)); }
template<class PR> FloatApproximation<PR> _round(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(round(x._a)); }

template<class PR> FloatApproximation<PR> _abs(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(abs_exact(x._a)); }
template<class PR> FloatApproximation<PR> _max(FloatApproximation<PR> const& x, FloatApproximation<PR> y) {
    return FloatApproximation<PR>(max_exact(x._a,y._a)); }
template<class PR> FloatApproximation<PR> _min(FloatApproximation<PR> const& x, FloatApproximation<PR> y) {
    return FloatApproximation<PR>(min_exact(x._a,y._a)); }
template<class PR> PositiveFloatApproximation<PR> _mag(FloatApproximation<PR> const& x) {
    return PositiveFloatApproximation<PR>(abs(x._a)); }
template<class PR> PositiveFloatApproximation<PR> _mig(FloatApproximation<PR> const& x) {
    return PositiveFloatApproximation<PR>(abs(x._a)); }

template<class PR> FloatApproximation<PR> _nul(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(nul_exact(x._a)); }
template<class PR> FloatApproximation<PR> _pos(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(pos_exact(x._a)); }
template<class PR> FloatApproximation<PR> _neg(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(neg_exact(x._a)); }
template<class PR> FloatApproximation<PR> _half(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(half_exact(x._a)); }
template<class PR> FloatApproximation<PR> _sqr(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(mul_near(x._a,x._a)); }
template<class PR> FloatApproximation<PR> _rec(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(div_near(1.0,x._a)); }

template<class PR> FloatApproximation<PR> _add(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return FloatApproximation<PR>(add_near(x1._a,x2._a)); }
template<class PR> FloatApproximation<PR> _sub(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return FloatApproximation<PR>(sub_near(x1._a,x2._a)); }
template<class PR> FloatApproximation<PR> _mul(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return FloatApproximation<PR>(mul_near(x1._a,x2._a)); }
template<class PR> FloatApproximation<PR> _div(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return FloatApproximation<PR>(div_near(x1._a,x2._a)); }

template<class PR> FloatApproximation<PR> _pow(FloatApproximation<PR> const& x, Nat m) {
    return FloatApproximation<PR>(pow_approx(x._a,m)); }
template<class PR> FloatApproximation<PR> _pow(FloatApproximation<PR> const& x, Int n) {
    return FloatApproximation<PR>(pow_approx(x._a,n)); }

template<class PR> FloatApproximation<PR> _sqrt(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(sqrt_approx(x._a)); }
template<class PR> FloatApproximation<PR> _exp(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(exp_approx(x._a)); }
template<class PR> FloatApproximation<PR> _log(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(log_approx(x._a)); }
template<class PR> FloatApproximation<PR> _sin(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(sin_approx(x._a)); }
template<class PR> FloatApproximation<PR> _cos(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(cos_approx(x._a)); }
template<class PR> FloatApproximation<PR> _tan(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(tan_approx(x._a)); }
template<class PR> FloatApproximation<PR> _asin(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(asin_approx(x._a)); }
template<class PR> FloatApproximation<PR> _acos(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(acos_approx(x._a)); }
template<class PR> FloatApproximation<PR> _atan(FloatApproximation<PR> const& x) {
    return FloatApproximation<PR>(atan_approx(x._a)); }

template<class PR> ApproximateKleenean _eq(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return x1._a==x2._a; }
template<class PR> ApproximateKleenean _leq(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return x1._a<=x2._a; }

template<class PR> Bool _same(FloatApproximation<PR> const& x1, FloatApproximation<PR> const& x2) {
    return x1._a==x2._a; }

template<class PR> OutputStream& _write(OutputStream& os, FloatApproximation<PR> const& x) {
    return write(os,x.raw(),FloatApproximation<PR>::output_places,RawFloat<PR>::to_nearest);
}

template<class PR> InputStream& _read(InputStream& is, FloatApproximation<PR>& x) {
    is >> x._a;
    return is;
}



template<class PR> FloatLowerBound<PR> _max(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return FloatLowerBound<PR>(max_exact(x1._l,x2._l)); }
template<class PR> FloatLowerBound<PR> _min(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return FloatLowerBound<PR>(min_exact(x1._l,x2._l)); }
template<class PR> FloatApproximation<PR> _abs(FloatLowerBound<PR> const& x) {
    return abs(FloatApproximation<PR>(x)); }

template<class PR> FloatLowerBound<PR> _nul(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(pos_exact(x._l)); }
template<class PR> FloatLowerBound<PR> _pos(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(pos_exact(x._l)); }
template<class PR> FloatLowerBound<PR> _neg(FloatUpperBound<PR> const& x) {
    return FloatLowerBound<PR>(neg_exact(x._u)); }
template<class PR> FloatLowerBound<PR> _half(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(half_exact(x._l)); }

template<class PR> FloatLowerBound<PR> _rec(FloatUpperBound<PR> const& x) {
    return FloatLowerBound<PR>(rec_down(x.raw())); }

template<class PR> FloatLowerBound<PR> _add(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return FloatLowerBound<PR>(add_down(x1._l,x2._l)); }

template<class PR> FloatApproximation<PR> _sub(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return FloatUpperBound<PR>(sub_near(x1._l,x2._l)); }

template<class PR> FloatLowerBound<PR> _sub(FloatLowerBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatLowerBound<PR>(sub_down(x1._l,x2._u)); }

template<class PR> FloatLowerBound<PR> _mul(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    return FloatLowerBound<PR>(mul_down(x1.raw(),x2.raw())); }

template<class PR> FloatLowerBound<PR> _div(FloatLowerBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatLowerBound<PR>(div_down(x1.raw(),x2.raw())); }

template<class PR> FloatLowerBound<PR> _pow(FloatLowerBound<PR> const& x, Nat m) {
    ARIADNE_PRECONDITION(x.raw()>=0);
    return FloatLowerBound<PR>(pow_down(x.raw(),m)); }

template<class PR> FloatApproximation<PR> _pow(FloatLowerBound<PR> const& x, Int n) {
    return pow(FloatApproximation<PR>(x),n); }

template<class PR> FloatLowerBound<PR> _sqrt(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(sqrt_down(x.raw())); }

template<class PR> FloatLowerBound<PR> _exp(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(exp_down(x.raw())); }

template<class PR> FloatLowerBound<PR> _log(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(log_down(x.raw())); }

template<class PR> FloatLowerBound<PR> _atan(FloatLowerBound<PR> const& x) {
    return FloatLowerBound<PR>(atan_down(x.raw())); }

template<class PR> ValidatedNegatedSierpinskian _eq(FloatLowerBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    if(x1._l>x2._u) { return false; }
    else { return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE); }
}

template<class PR> ValidatedNegatedSierpinskian _leq(FloatLowerBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    if(x1._l>x2._u) { return false; }
    else { return ValidatedNegatedSierpinskian(LogicalValue::LIKELY); }
}

template<class PR> Bool _same(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return x1._l==x2._l; }

template<class PR> FloatLowerBound<PR> operator+(FloatLowerBound<PR> const& x) { return pos(x); }
template<class PR> FloatLowerBound<PR> operator-(FloatUpperBound<PR> const& x) { return neg(x); }
template<class PR> FloatLowerBound<PR> operator+(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) { return add(x1,x2); }
template<class PR> FloatLowerBound<PR> operator-(FloatLowerBound<PR> const& x1, FloatUpperBound<PR> const& x2) { return sub(x1,x2); }
template<class PR> FloatLowerBound<PR> operator*(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) { return mul(x1,x2); }
template<class PR> FloatLowerBound<PR> operator/(FloatLowerBound<PR> const& x1, FloatUpperBound<PR> const& x2) { return div(x1,x2); }
template<class PR> FloatLowerBound<PR>& operator+=(FloatLowerBound<PR>& x1, FloatLowerBound<PR> const& x2) { return x1=x1+x2; }
template<class PR> FloatLowerBound<PR>& operator-=(FloatLowerBound<PR>& x1, FloatUpperBound<PR> const& x2) { return x1=x1-x2; }
template<class PR> FloatLowerBound<PR>& operator*=(FloatLowerBound<PR>& x1, FloatLowerBound<PR> const& x2) { return x1=x1*x2; }
template<class PR> FloatLowerBound<PR>& operator/=(FloatLowerBound<PR>& x1, FloatUpperBound<PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& _write(OutputStream& os, FloatLowerBound<PR> const& x) {
    return write(os,x.raw(),FloatBounds<PR>::output_places,RawFloat<PR>::downward);
}

template<class PR> InputStream& _read(InputStream& is, FloatLowerBound<PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
}



template<class PR> FloatUpperBound<PR> _max(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatUpperBound<PR>(max_exact(x1._u,x2._u)); }

template<class PR> FloatUpperBound<PR> _min(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatUpperBound<PR>(min_exact(x1._u,x2._u)); }

template<class PR> FloatApproximation<PR> _abs(FloatUpperBound<PR> const& x) {
    return abs(FloatApproximation<PR>(x)); }

template<class PR> FloatUpperBound<PR> _nul(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(pos_exact(x._u)); }

template<class PR> FloatUpperBound<PR> _pos(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(pos_exact(x._u)); }

template<class PR> FloatUpperBound<PR> _neg(FloatLowerBound<PR> const& x) {
    return FloatUpperBound<PR>(neg_exact(x._l)); }

template<class PR> FloatUpperBound<PR> _half(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(half_exact(x._u)); }

template<class PR> FloatUpperBound<PR> _sqr(FloatUpperBound<PR> const& x) {
    ARIADNE_ASSERT(false); return FloatUpperBound<PR>(mul_up(x._u,x._u)); }

template<class PR> FloatUpperBound<PR> _rec(FloatLowerBound<PR> const& x) {
    return FloatUpperBound<PR>(rec_up(x.raw())); }

template<class PR> FloatUpperBound<PR> _add(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatUpperBound<PR>(add_up(x1._u,x2._u)); }

template<class PR> FloatApproximation<PR> _sub(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatUpperBound<PR>(sub_near(x1._u,x2._u)); }

template<class PR> FloatUpperBound<PR> _sub(FloatUpperBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return FloatUpperBound<PR>(sub_up(x1._u,x2._l)); }

template<class PR> FloatUpperBound<PR> _mul(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
//    ARIADNE_WARN("Multiplying FloatUpperBound "<<x1<<" with FloatUpperBound "<<x2<<" is unsafe");
    ARIADNE_PRECONDITION(x1.raw()>=0);
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return FloatUpperBound<PR>(mul_up(x1._u,x2._u)); }

template<class PR> FloatUpperBound<PR> _div(FloatUpperBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
//    ARIADNE_WARN("Dividing FloatUpperBound "<<x1<<" by FloatLowerBound "<<x2<<" is unsafe");
    ARIADNE_PRECONDITION(x1.raw()>=0);
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return FloatUpperBound<PR>(div_up(x1._u,x2._l)); }

template<class PR> FloatUpperBound<PR> _pow(FloatUpperBound<PR> const& x, Nat m) {
    ARIADNE_PRECONDITION(x.raw()>=0);
    return FloatUpperBound<PR>(pow_up(x._u,m)); }

template<class PR> FloatApproximation<PR> _pow(FloatUpperBound<PR> const& x, Int n) {
    return pow(FloatApproximation<PR>(x),n); }

template<class PR> FloatUpperBound<PR> _sqrt(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(sqrt_up(x.raw())); }

template<class PR> FloatUpperBound<PR> _exp(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(exp_up(x.raw())); }

template<class PR> FloatUpperBound<PR> _log(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(log_up(x.raw())); }

template<class PR> FloatUpperBound<PR> _atan(FloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(atan_up(x.raw())); }

template<class PR> ValidatedNegatedSierpinskian _eq(FloatUpperBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    if(x1._u<x2._l) { return false; }
    else { return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE); }
}

template<class PR> ValidatedSierpinskian _leq(FloatUpperBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    if(x1._u<=x2._l) { return true; }
    else { return ValidatedSierpinskian(LogicalValue::UNLIKELY); }
}

template<class PR> Bool _same(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return x1._u==x2._u; }

template<class PR> Integer integer_cast(FloatUpperBound<PR> const& x) { return Integer(static_cast<int>(x._u.get_d())); }

template<class PR> FloatUpperBound<PR> operator+(FloatUpperBound<PR> const& x) { return pos(x); }
template<class PR> FloatUpperBound<PR> operator-(FloatLowerBound<PR> const& x) { return neg(x); }
template<class PR> FloatUpperBound<PR> operator+(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) { return add(x1,x2); }
template<class PR> FloatUpperBound<PR> operator-(FloatUpperBound<PR> const& x1, FloatLowerBound<PR> const& x2) { return sub(x1,x2); }
template<class PR> FloatUpperBound<PR> operator*(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) { return mul(x1,x2); }
template<class PR> FloatUpperBound<PR> operator/(FloatUpperBound<PR> const& x1, FloatLowerBound<PR> const& x2) { return div(x1,x2); }
template<class PR> FloatUpperBound<PR>& operator+=(FloatUpperBound<PR>& x1, FloatUpperBound<PR> const& x2) { return x1=x1+x2; }
template<class PR> FloatUpperBound<PR>& operator-=(FloatUpperBound<PR>& x1, FloatLowerBound<PR> const& x2) { return x1=x1-x2; }
template<class PR> FloatUpperBound<PR>& operator*=(FloatUpperBound<PR>& x1, FloatUpperBound<PR> const& x2) { return x1=x1*x2; }
template<class PR> FloatUpperBound<PR>& operator/=(FloatUpperBound<PR>& x1, FloatLowerBound<PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& _write(OutputStream& os, FloatUpperBound<PR> const& x) {
    return write(os,x.raw(),FloatBounds<PR>::output_places,RawFloat<PR>::upward);
}

template<class PR> InputStream& _read(InputStream& is, FloatUpperBound<PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
}








template<class PR> FloatBounds<PR> _round(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(round(x.lower_raw()),round(x.upper_raw()));
}

template<class PR> FloatBounds<PR> _max(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
}

template<class PR> FloatBounds<PR> _min(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
}


template<class PR> FloatBounds<PR> _abs(FloatBounds<PR> const& x) {
    if(x.lower_raw()>=0) {
        return FloatBounds<PR>(x.lower_raw(),x.upper_raw());
    } else if(x.upper_raw()<=0) {
        return FloatBounds<PR>(neg(x.upper_raw()),neg(x.lower_raw()));
    } else {
        return FloatBounds<PR>(static_cast<RawFloat<PR>>(0.0,x.precision()),max(neg(x.lower_raw()),x.upper_raw()));
    }
}

template<class PR> PositiveFloatLowerBound<PR> _mig(FloatBounds<PR> const& x) {
    return PositiveFloatLowerBound<PR>(max(0,max(x._l,neg(x._u))));
}

template<class PR> PositiveFloatUpperBound<PR> _mag(FloatBounds<PR> const& x) {
    return PositiveFloatUpperBound<PR>(max(neg(x._l),x._u));
}

template<class PR> FloatBounds<PR> _nul(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(nul(x._l),nul(x._u));
}

template<class PR> FloatBounds<PR> _pos(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(pos(x._l),pos(x._u));
}

template<class PR> FloatBounds<PR> _neg(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(neg(x._u),neg(x._l));
}

template<class PR> FloatBounds<PR> _half(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(half(x._l),half(x._u));
}

template<class PR> FloatBounds<PR> _sqr(FloatBounds<PR> const& x) {
    const RawFloat<PR>& xl=x.lower_raw(); const RawFloat<PR>& xu=x.upper_raw();
    RawFloat<PR> rl,ru;
    if(xl>0.0) {
        rl=mul_down(xl,xl); ru=mul_up(xu,xu);
    } else if(xu<0.0) {
        rl=mul_down(xu,xu); ru=mul_up(xl,xl);
    } else {
        rl=nul(xl); ru=max(mul_up(xl,xl),mul_up(xu,xu));
    }
    return FloatBounds<PR>(rl,ru);
}

template<class PR> FloatBounds<PR> _rec(FloatBounds<PR> const& x) {
    // IMPORTANT: Need to be careful when one of the bounds is 0, since if xl=-0.0 and xu>0, then 1/xl=-inf
    if(x._l>0 || x._u<0) {
        return FloatBounds<PR>(rec_down(x._u),rec_up(x._l));
    } else {
        RawFloat<PR> inf=RawFloat<PR>::inf(x.precision());
        RawFloat<PR> rl=-inf; RawFloat<PR> ru=+inf;
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds rec(FloatBounds x)","x="<<x);
        return FloatBounds<PR>(-inf,+inf);
    }
}

template<class PR> FloatBounds<PR> _add(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(add_down(x1._l,x2._l),add_up(x1._u,x2._u));
}

template<class PR> FloatBounds<PR> _sub(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(sub_down(x1._l,x2._u),sub_up(x1._u,x2._l));
}

template<class PR> FloatBounds<PR> _mul(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
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
    return FloatBounds<PR>(rl,ru);
}

template<class PR> FloatBounds<PR> _div(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
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
    return FloatBounds<PR>(rl,ru);
}






template<class PR> FloatBounds<PR> _pow(FloatBounds<PR> const& x, Int n) {
    if(n<0) { return pow(rec(x),Nat(-n)); }
    else return pow(x,Nat(n));
}

template<class PR> FloatBounds<PR> _pow(FloatBounds<PR> const& x, Nat m) {
    FloatBounds<PR> y = x;
    if(m%2==0) { y=abs(x); }
    RawFloat<PR> rl=pow_down(y.lower_raw(),m);
    RawFloat<PR> ru=pow_up(y.upper_raw(),m);
    return FloatBounds<PR>(rl,ru);
}


template<class PR> FloatBounds<PR> _sqrt(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(sqrt_down(x.lower_raw()),sqrt_up(x.upper_raw()));
}

template<class PR> FloatBounds<PR> _exp(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(exp_down(x.lower_raw()),exp_up(x.upper_raw()));
}

template<class PR> FloatBounds<PR> _log(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(log_down(x.lower_raw()),log_up(x.upper_raw()));
}


template<class PR> FloatBounds<PR> _pi_val(PR pr) { return FloatBounds<PR>(pi_down(pr),pi_up(pr)); }

template<class PR> FloatBounds<PR> _sin(FloatBounds<PR> const& x)
{
    return cos(x-half(_pi_val<PR>(x.precision())));
}

template<class PR> FloatBounds<PR> _cos(FloatBounds<PR> const& x)
{
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    typename RawFloat<PR>::RoundingModeType rnd = RawFloat<PR>::get_rounding_mode();
    PR prec=x.precision();

    static const RawFloat<PR> one(1,prec);
    static const FloatValue<PR> two(2,prec);

    if(x.error().raw()>2*pi_down(prec)) { return FloatBounds<PR>(-one,+one); }

    auto n=floor(x.lower_raw()/(2*pi_approx(prec))+one/2);
    FloatBounds<PR> y=x-two*FloatValue<PR>(n)*_pi_val<PR>(x.precision());

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
    return FloatBounds<PR>(rl,ru);
}

template<class PR> FloatBounds<PR> _tan(FloatBounds<PR> const& x) {
    return mul(sin(x),rec(cos(x)));
}

template<class PR> FloatBounds<PR> _asin(FloatBounds<PR> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class PR> FloatBounds<PR> _acos(FloatBounds<PR> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class PR> FloatBounds<PR> _atan(FloatBounds<PR> const& x) {
    return FloatBounds<PR>(atan_down(x._l),atan_up(x._u));
}



//! \related FloatBounds<PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> Logical<ValidatedTag> _eq(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    if(x1.upper_raw()<x2.lower_raw() || x1.lower_raw()>x2.upper_raw()) { return false; }
    else if(x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw()) { return true; }
    else { return indeterminate; }
}

//! \related FloatBounds<PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> Logical<ValidatedTag> _leq(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    if(x1.upper_raw()<=x2.lower_raw()) { return true; }
    else if(x1.lower_raw()> x2.upper_raw()) { return false; }
    else { return indeterminate; }
}


template<class PR> FloatBounds<PR> _widen(FloatBounds<PR> const& x)
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
    return FloatBounds<PR>(wl,wu);
}

template<class PR> FloatBounds<PR> _narrow(FloatBounds<PR> const& x)
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
    return FloatBounds<PR>(nl,nu);
}

template<class PR> FloatBounds<PR> _trunc(FloatBounds<PR> const& x)
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
    return FloatBounds<PR>(double(tl),double(tu));
}

template<class PR> FloatBounds<PR> _trunc(FloatBounds<PR> const& x, Nat n)
{
    FloatBounds<PR> _e=FloatBounds<PR>(std::pow(2.0,52-(Int)n));
    FloatBounds<PR> y=x+_e;
    return y-_e;
}

template<class PR> Integer integer_cast(FloatBounds<PR> const& x) {
    return Integer(static_cast<int>(x.value_raw().get_d()));
}

template<class PR> auto is_zero(FloatBounds<PR> const& x) -> Logical<ValidatedTag> {
    if(x.lower_raw()>0.0 || x.upper_raw()<0.0) { return false; }
    else if(x.lower_raw()==0.0 && x.upper_raw()==0.0) { return true; }
    else { return indeterminate; }
}

template<class PR> auto is_positive(FloatBounds<PR> const& x) -> Logical<ValidatedTag> {
    if(x.lower_raw()>=0.0) { return true; }
    else if(x.upper_raw()<0.0) { return false; }
    else { return indeterminate; }
}

template<class PR> auto is_positive(FloatBounds<PR> const&) -> Logical<ValidatedTag>;

template<class PR> Bool _same(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return x1._l==x2._l && x1._u==x2._u; }


inline int log10floor(double const& x) { return std::max(std::floor(std::log10(x)),-65280.); }
inline int log10floor(FloatMP const& x) { return log10floor(x.get_d()); }
inline int abslog10floor(double const& x) { return log10floor(std::abs(x)); }


template<class PR> OutputStream& _write(OutputStream& os, const FloatBounds<PR>& x)
{
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    os << '{';
    write(os,x.lower().raw(),FloatBounds<PR>::output_places,RawFloat<PR>::downward);
    os << ':';
    write(os,x.upper().raw(),FloatBounds<PR>::output_places,RawFloat<PR>::upward);
    os << '}';
    return os;

}
template<> OutputStream& _write(OutputStream& os, const FloatBounds<PrecisionMP>& x)
{
    static const double log2ten = 3.3219280948873621817;
    using std::max; using std::min;
    FloatMP const& l=x.lower_raw();
    FloatMP const& u=x.upper_raw();
    double ldbl=l.get_d();
    double udbl=u.get_d();
    if(ldbl==0.0 && udbl==0.0) { return os << "0.0[:]"; }
    int errplc=FloatError<PrecisionMP>::output_places;
    int bndplc=FloatBounds<PrecisionMP>::output_places;
    int precplc=x.precision()/log2ten;
    int log10wdth=log10floor(u-l);
    int log10mag=log10floor(max(-ldbl,udbl));
    int dgtswdth=errplc-(log10wdth+1); // Digits appropriate given width of interval
    int dgtsbnd=bndplc-(log10mag+1); // Digits appropriate given asked-for precision of bounded objects
    int dgtsprec=precplc-(log10mag+1); // Digits appropriate given precision of objects
    int dgts=max(min(dgtswdth,dgtsprec),1);
    DecimalPlaces plcs{dgts};
    String lstr=print(l,plcs,MPFR_RNDD);
    String ustr=print(u,plcs,MPFR_RNDU);
    auto lcstr=lstr.c_str();
    auto ucstr=ustr.c_str();
    size_t cpl=0;
    if(ldbl*udbl>=0 && abslog10floor(ldbl)==abslog10floor(udbl)) {
        while(lcstr[cpl]!='\0' && lcstr[cpl]==ustr[cpl]) { ++cpl; }
    }
    char ocstr[1024];
    ocstr[0]='\0';
    strncat(ocstr,lcstr,cpl);
    strcat(ocstr,"[");
    strcat(ocstr,lcstr+cpl);
    strcat(ocstr,":");
    strcat(ocstr,ucstr+cpl);
    strcat(ocstr,"]");
    return os << ocstr;
}

template<> OutputStream& _write(OutputStream& os, const FloatBounds<Precision64>& x)
{
    PrecisionMP prec(64);
    return os << FloatBounds<PrecisionMP>(FloatMP(x.lower_raw(),prec),FloatMP(x.upper_raw(),prec));
}


template<class PR> InputStream& _read(InputStream& is, FloatBounds<PR>& x)
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



template<class PR> FloatBounds<PR> operator+(FloatBounds<PR> const& x) { return pos(x); }
template<class PR> FloatBounds<PR> operator-(FloatBounds<PR> const& x) { return neg(x); }
template<class PR> FloatBounds<PR> operator+(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) { return add(x1,x2); }
template<class PR> FloatBounds<PR> operator-(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) { return sub(x1,x2); }
template<class PR> FloatBounds<PR> operator*(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) { return mul(x1,x2); }
template<class PR> FloatBounds<PR> operator/(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) { return div(x1,x2); }
template<class PR> FloatBounds<PR>& operator+=(FloatBounds<PR>& x1, FloatBounds<PR> const& x2) { return x1=x1+x2; }
template<class PR> FloatBounds<PR>& operator-=(FloatBounds<PR>& x1, FloatBounds<PR> const& x2) { return x1=x1-x2; }
template<class PR> FloatBounds<PR>& operator*=(FloatBounds<PR>& x1, FloatBounds<PR> const& x2) { return x1=x1*x2; }
template<class PR> FloatBounds<PR>& operator/=(FloatBounds<PR>& x1, FloatBounds<PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& _write(OutputStream& os, FloatBounds<PR> const& x);
template<class PR> InputStream& _read(InputStream& is, FloatBounds<PR>& x);




template<class PR> FloatBall<PR> _nul(FloatBall<PR> const& x) {
    return FloatBall<PR>(nul(x._v),nul(x._e));
}

template<class PR> FloatBall<PR> _pos(FloatBall<PR> const& x) {
    return FloatBall<PR>(pos(x._v),x._e);
}

template<class PR> FloatBall<PR> _neg(FloatBall<PR> const& x) {
    return FloatBall<PR>(neg(x._v),x._e);
}

template<class PR> FloatBall<PR> _half(FloatBall<PR> const& x) {
    return FloatBall<PR>(half(x._v),half(x._e));
}

template<class PR> FloatBall<PR> _sqr(FloatBall<PR> const& x) {
    FloatBall<PR> r=x*x;
    if(r._e>r._v) {
        r._e=half(add_up(r._e,r._v));
        r._v=r._e;
    }
    return r;
}

template<class PR> FloatBall<PR> _rec(FloatBall<PR> const& x) {
    // Use this code to find value same as reciprocal value
    auto rv=rec_approx(x._v);
    auto ru=rec_up(sub_down(x._v,x._e));
    auto rl=rec_down(add_up(x._v,x._e));
    auto re=max(sub_up(ru,rv),sub_up(rv,rl));
    return FloatBall<PR>(rv,re);
#ifdef ARIADNE_UNDEFINED
    // Use this code to get same result as interval computation
    auto ru=rec_up(sub_down(x._v,x._e));
    auto rl=rec_down(add_up(x._v,x._e));
    auto re=half(sub_up(ru,rl));
    auto rv=half(add_near(rl,ru));
    return FloatBall<PR>(rv,re);
#endif
}

template<class PR> FloatBall<PR> _add(FloatBall<PR> const& x, FloatBall<PR> y) {
    auto rv=add_near(x._v,y._v);
    auto ru=add_up(x._v,y._v);
    auto rl=add_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return FloatBall<PR>(rv,re);
}

template<class PR> FloatBall<PR> _sub(FloatBall<PR> const& x, FloatBall<PR> y) {
    auto rv=sub_near(x._v,y._v);
    auto ru=sub_up(x._v,y._v);
    auto rl=sub_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return FloatBall<PR>(rv,re);
}

template<class PR> FloatBall<PR> _mul(FloatBall<PR> const& x, FloatBall<PR> y) {
    auto rv=mul_near(x._v,y._v);
    auto ru=mul_up(x._v,y._v);
    auto rl=mul_down(x._v,y._v);
    auto re1=add_up(half(sub_up(ru,rl)),mul_up(x._e,y._e));
    auto re2=add_up(mul_up(abs(x._v),y._e),mul_up(x._e,abs(y._v)));
    auto re=add_up(re1,re2);
    return FloatBall<PR>(rv,re);
}

template<class PR> FloatBall<PR> _div(FloatBall<PR> const& x, FloatBall<PR> y) {
    return x*rec(y);
}

template<class PR> FloatBall<PR> _pow(FloatBall<PR> const& x, Nat m) {
    return FloatBall<PR>(pow(FloatBounds<PR>(x),m));
}

template<class PR> FloatBall<PR> _pow(FloatBall<PR> const& x, Int n) {
    return FloatBall<PR>(pow(FloatBounds<PR>(x),n));
}

template<class PR> FloatBall<PR> _sqrt(FloatBall<PR> const& x) {
    return FloatBall<PR>(sqrt(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _exp(FloatBall<PR> const& x) {
    return FloatBall<PR>(exp(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _log(FloatBall<PR> const& x) {
    return FloatBall<PR>(log(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _sin(FloatBall<PR> const& x) {
    return FloatBall<PR>(sin(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _cos(FloatBall<PR> const& x) {
    return FloatBall<PR>(cos(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _tan(FloatBall<PR> const& x) {
    return FloatBall<PR>(tan(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _asin(FloatBall<PR> const& x) {
    return FloatBall<PR>(asin(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _acos(FloatBall<PR> const& x) {
    return FloatBall<PR>(acos(FloatBounds<PR>(x)));
}

template<class PR> FloatBall<PR> _atan(FloatBall<PR> const& x) {
    return FloatBall<PR>(atan(FloatBounds<PR>(x)));
}


template<class PR> FloatBall<PR> _abs(FloatBall<PR> const& x) {
    if(x._e<abs(x._v)) { return x; }
    else { auto rv=half(abs(x._v)+x._e); return FloatBall<PR>(rv,rv); }
}

template<class PR> FloatBall<PR> _max(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return half((x1+x2)+abs(x1-x2));
}

template<class PR> FloatBall<PR> _min(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return half((x1+x2)-abs(x1-x2));
}

template<class PR> FloatError<PR> _mag(FloatBall<PR> const& x) {
    return PositiveFloatUpperBound<PR>(max(x._e+x._v,x._e-x._v));
}

template<class PR> PositiveFloatLowerBound<PR> _mig(FloatBall<PR> const& x) {
    return PositiveFloatLowerBound<PR>(max(0,max(x._v-x._e,-x._v-x._e)));
}

//! \related FloatBounds<PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> ValidatedKleenean _eq(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return _eq(FloatBounds<PR>(x1),FloatBounds<PR>(x2));
}

//! \related FloatBounds<PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> ValidatedKleenean _leq(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return _leq(FloatBounds<PR>(x1),FloatBounds<PR>(x2));
}

template<class PR> Bool _same(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return x1._v==x2._v && x1._e==x2._e; }

template<class PR> OutputStream& _write(OutputStream& os, FloatBall<PR> const& x) {
    return os << x.value() << "\u00b1" << x.error();
}

template<> OutputStream& _write(OutputStream& os, FloatBall<PrecisionMP> const& x) {
    // Write based on number of correct digits
    static const double log2ten = 3.3219280948873621817;
    static const char pmstr[] = "\u00b1";
    static const char hlfstr[] = "\u00bd";
    FloatMP const& v=x.value_raw();
    FloatMP const& e=x.error_raw();
    double edbl=e.get_d();
    // Compute the number of decimal places to be displayed
    int errplc = FloatError<PrecisionMP>::output_places;
    int log10err = log10floor(edbl);
    int dgtserr = errplc-(log10err+1);
    int dgtsval = std::floor((x.value().precision()-x.value().raw().exponent())/log2ten);
    int dgts = std::max(std::min(dgtsval,dgtserr),1);
    if(edbl==0.0) { dgts = dgtsval; }
    DecimalPlaces plcs{dgts};

    // Get string version of mpfr values
    String vstr=print(v,plcs,MPFR_RNDN);
    String estr=print(e,plcs,MPFR_RNDU);

    // Find position of first significan digit of error
    auto vcstr=vstr.c_str(); auto ecstr=estr.c_str();
    int cpl=0;
    if(edbl==0.0) {
        cpl=std::strlen(vcstr);
    } else if(edbl<1.0) {
        const char* vptr = std::strchr(vcstr,'.');
        const char* eptr = std::strchr(ecstr,'.');
        ++vptr; ++eptr;
        while((*eptr)=='0') { ++eptr; ++vptr; }
        cpl = vptr-vcstr;
    }

    // Chop and catenate strings
    static const size_t buf_sz = 1024;
    char ocstr[buf_sz];
    ocstr[0]='\0';
    std::strncat(ocstr,vcstr,cpl);
    std::strcat(ocstr,"[");
    std::strcat(ocstr,vcstr+cpl);
    std::strcat(ocstr,pmstr);
    std::strcat(ocstr,ecstr+cpl);
    std::strcat(ocstr,hlfstr);
    std::strcat(ocstr,"]");
    return os << ocstr;

    return os << x.value() << "\u00b1" << x.error();
}

template<> OutputStream& _write(OutputStream& os, FloatBall<Precision64> const& x) {
    PrecisionMP prec(64);
    return os << FloatBall<PrecisionMP>(FloatMP(x.value_raw(),prec),FloatMP(x.error_raw(),prec));
}

template<class PR> InputStream& _read(InputStream& is, FloatBall<PR>& x)
{
    static const char pmstr[] = "\u00b1";
    char cpm[3];
    RawFloat<PR> _v,_e;
    auto rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_to_nearest();
    is >> _v;
    is >> cpm[0] >> cpm[1];
    RawFloat<PR>::set_rounding_upward();
    is >> _e;
    RawFloat<PR>::set_rounding_mode(rnd);
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(std::strcmp(cpm,pmstr));
    x._v=_v; x._e=_e;
    return is;
}


// Mixed BoundedTag - ExactTag operations
template<class PR> FloatBounds<PR> _add(FloatBounds<PR> const& x1, FloatValue<PR> const& x2) {
    return FloatBounds<PR>(add_down(x1._l,x2._v),add_up(x1._u,x2._v));
}

template<class PR> FloatBounds<PR> _add(FloatValue<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(add_down(x1._v,x2._l),add_down(x1._v,x2._u));
}

template<class PR> FloatBounds<PR> _sub(FloatBounds<PR> const& x1, FloatValue<PR> const& x2) {
    return FloatBounds<PR>(sub_down(x1._l,x2._v),sub_up(x1._u,x2._v));
}

template<class PR> FloatBounds<PR> _sub(FloatValue<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(sub_down(x1._v,x2._u),sub_up(x1._v,x2._l));
}

template<class PR> FloatBounds<PR> _mul(FloatBounds<PR> const& x1, FloatValue<PR> const& x2) {
    const RawFloat<PR>& x1l=x1.lower_raw(); const RawFloat<PR>& x1u=x1.upper_raw();
    const RawFloat<PR>& x2v=x2.raw();
    RawFloat<PR> rl,ru;
    if(x2v>=0.0) {
        rl=mul_down(x1l,x2v); ru=mul_up(x1u,x2v);
    } else {
        rl=mul_down(x1u,x2v); ru=mul_up(x1l,x2v);
    }
    return FloatBounds<PR>(rl,ru);
}


template<class PR> FloatBounds<PR> _mul(FloatValue<PR> const& x1, FloatBounds<PR> const& x2) {
    const RawFloat<PR>& x1v=x1.raw();
    const RawFloat<PR>& x2l=x2.lower_raw(); const RawFloat<PR>& x2u=x2.upper_raw();
    RawFloat<PR> rl,ru;
    if(x1v>=0.0) {
        rl=mul_down(x1v,x2l); ru=mul_up(x1v,x2u);
    } else {
        rl=mul_down(x1v,x2u); ru=mul_up(x1v,x2l);
    }
    return FloatBounds<PR>(rl,ru);
}

template<class PR> FloatBounds<PR> _div(FloatBounds<PR> const& x1, FloatValue<PR> const& x2)
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
    return FloatBounds<PR>(rl,ru);
}


template<class PR> FloatBounds<PR> _div(FloatValue<PR> const& x1, FloatBounds<PR> const& x2)
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
    return FloatBounds<PR>(rl,ru);
}






template<class PR> FloatValue<PR> _max(FloatValue<PR> const& x1,  FloatValue<PR> const& x2) {
    return FloatValue<PR>(max(x1._v,x2._v)); }

template<class PR> FloatValue<PR> _min(FloatValue<PR> const& x1,  FloatValue<PR> const& x2) {
    return FloatValue<PR>(min(x1._v,x2._v)); }

template<class PR> FloatValue<PR> _abs(FloatValue<PR> const& x) {
    return FloatValue<PR>(abs(x._v)); }

template<class PR> FloatLowerBound<PR> _mig(FloatValue<PR> const& x) {
    return FloatLowerBound<PR>(abs(x._v)); }

template<class PR> FloatError<PR> _mag(FloatValue<PR> const& x) {
    return FloatError<PR>(abs(x._v)); }


template<class PR> FloatValue<PR> _nul(FloatValue<PR> const& x) {
    return FloatValue<PR>(nul(x._v)); }

template<class PR> FloatValue<PR> _pos(FloatValue<PR> const& x) {
    return FloatValue<PR>(pos(x._v)); }

template<class PR> FloatValue<PR> _neg(FloatValue<PR> const& x) {
    return FloatValue<PR>(neg(x._v)); }

template<class PR> FloatValue<PR> _half(FloatValue<PR> const& x) {
    return FloatValue<PR>(half(x._v)); }

template<class PR> FloatBounds<PR> _sqr(FloatValue<PR> const& x) {
    return FloatBounds<PR>(mul_down(x._v,x._v),mul_up(x._v,x._v)); }

template<class PR> FloatBounds<PR> _rec(FloatValue<PR> const& x) {
    return FloatBounds<PR>(rec_down(x._v),rec_up(x._v)); }

template<class PR> FloatBounds<PR> _add(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return FloatBounds<PR>(add_down(x1._v,x2._v),add_up(x1._v,x2._v)); }

template<class PR> FloatBounds<PR> _sub(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return FloatBounds<PR>(sub_down(x1._v,x2._v),sub_up(x1._v,x2._v)); }

template<class PR> FloatBounds<PR> _mul(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return FloatBounds<PR>(mul_down(x1._v,x2._v),mul_up(x1._v,x2._v)); }

template<class PR> FloatBounds<PR> _div(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return FloatBounds<PR>(div_down(x1._v,x2._v),div_up(x1._v,x2._v)); }

template<class PR> FloatBounds<PR> _pow(FloatValue<PR> const& x, Nat m) {
    return pow(FloatBounds<PR>(x),m); }

template<class PR> FloatBounds<PR> _pow(FloatValue<PR> const& x, Int n) {
    return pow(FloatBounds<PR>(x),n); }

template<class PR> FloatBounds<PR> _med(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return add(half(x1),half(x2)); }

template<class PR> FloatBounds<PR> _rad(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return sub(half(x2),half(x1)); }

template<class PR> FloatBounds<PR> _sqrt(FloatValue<PR> const& x) {
    return sqrt(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _exp(FloatValue<PR> const& x) {
    return exp(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _log(FloatValue<PR> const& x) {
    return log(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _sin(FloatValue<PR> const& x) {
    return sin(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _cos(FloatValue<PR> const& x) {
    return cos(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _tan(FloatValue<PR> const& x) {
    return tan(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _asin(FloatValue<PR> const& x) {
    return asin(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _acos(FloatValue<PR> const& x) {
    return acos(FloatBounds<PR>(x)); }

template<class PR> FloatBounds<PR> _atan(FloatValue<PR> const& x) {
    return atan(FloatBounds<PR>(x)); }

template<class PR> Boolean _eq(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return x1._v == x2._v; }

template<class PR> Boolean _leq(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return x1._v <= x2._v; }

template<class PR> Bool _same(FloatValue<PR> const& x1, FloatValue<PR> const& x2) {
    return x1._v==x2._v; }

template<class PR> FloatValue<PR> operator+(FloatValue<PR> const& x) { return pos(x); }
template<class PR> FloatValue<PR> operator-(FloatValue<PR> const& x) { return neg(x); }
template<class PR> FloatBounds<PR> operator+(FloatValue<PR> const& x1,  FloatValue<PR> const& x2) { return add(x1,x2); }
template<class PR> FloatBounds<PR> operator-(FloatValue<PR> const& x1,  FloatValue<PR> const& x2) { return sub(x1,x2); }
template<class PR> FloatBounds<PR> operator*(FloatValue<PR> const& x1,  FloatValue<PR> const& x2) { return mul(x1,x2); }
template<class PR> FloatBounds<PR> operator/(FloatValue<PR> const& x1,  FloatValue<PR> const& x2) { return div(x1,x2); }

template<class PR> FloatValue<PR> operator*(FloatValue<PR> const& x, TwoExp y) {
    FloatValue<PR> yv=y; return FloatValue<PR>(x.raw()*yv.raw()); }
template<class PR> FloatValue<PR> operator/(FloatValue<PR> const& x, TwoExp y) {
    FloatValue<PR> yv=y; return FloatValue<PR>(x.raw()/yv.raw()); }
template<class PR> FloatValue<PR>& operator*=(FloatValue<PR>& x, TwoExp y) {
    FloatValue<PR> yv=y; return x=FloatValue<PR>(x.raw()*yv.raw()); }
template<class PR> FloatValue<PR>& operator/=(FloatValue<PR>& x, TwoExp y) {
    FloatValue<PR> yv=y; return x=FloatValue<PR>(x.raw()/yv.raw()); }

template<class PR> Boolean operator==(FloatValue<PR> const& x1, FloatValue<PR> const& x2) { return x1.raw()==x2.raw(); }
template<class PR> Boolean operator!=(FloatValue<PR> const& x1, FloatValue<PR> const& x2) { return x1.raw()!=x2.raw(); }
template<class PR> Boolean operator<=(FloatValue<PR> const& x1, FloatValue<PR> const& x2) { return x1.raw()<=x2.raw(); }
template<class PR> Boolean operator>=(FloatValue<PR> const& x1, FloatValue<PR> const& x2) { return x1.raw()>=x2.raw(); }
template<class PR> Boolean operator< (FloatValue<PR> const& x1, FloatValue<PR> const& x2) { return x1.raw()< x2.raw(); }
template<class PR> Boolean operator> (FloatValue<PR> const& x1, FloatValue<PR> const& x2) { return x1.raw()> x2.raw(); }

template<class PR> OutputStream& _write(OutputStream& os, FloatValue<PR> const& x) {
    return write(os,x.raw(),FloatValue<PR>::output_places,RawFloat<PR>::to_nearest);
}

template<class PR> InputStream& _read(InputStream& is, FloatValue<PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
    auto v = nul(x._v);
    is >> v;
    ARIADNE_ASSERT(is);
    x._v=v;
    return is;
}

template<class PR> Integer integer_cast(FloatValue<PR> const& x) {
    Integer z=static_cast<int>(x._v.get_d());
    ARIADNE_ASSERT(z==x);
    return std::move(z);
}


FloatValue<Precision64> cast_exact(Real const& x) {
    return cast_exact(FloatApproximation<Precision64>(x));
}

template<class PR> Bool operator==(FloatValue<PR> const& x, const Rational& q) { return Rational(x)==q; }
template<class PR> Bool operator!=(FloatValue<PR> const& x, const Rational& q) { return Rational(x)!=q; }
template<class PR> Bool operator<=(FloatValue<PR> const& x, const Rational& q) { return Rational(x)<=q; }
template<class PR> Bool operator>=(FloatValue<PR> const& x, const Rational& q) { return Rational(x)>=q; }
template<class PR> Bool operator< (FloatValue<PR> const& x, const Rational& q) { return Rational(x)< q; }
template<class PR> Bool operator> (FloatValue<PR> const& x, const Rational& q) { return Rational(x)> q; }

template<class PR> Bool operator==(const Rational& q, FloatValue<PR> const& x) { return q==Rational(x); }
template<class PR> Bool operator!=(const Rational& q, FloatValue<PR> const& x) { return q!=Rational(x); }
template<class PR> Bool operator<=(const Rational& q, FloatValue<PR> const& x) { return q<=Rational(x); }
template<class PR> Bool operator>=(const Rational& q, FloatValue<PR> const& x) { return q>=Rational(x); }
template<class PR> Bool operator< (const Rational& q, FloatValue<PR> const& x) { return q< Rational(x); }
template<class PR> Bool operator> (const Rational& q, FloatValue<PR> const& x) { return q> Rational(x); }



template<class PR> FloatError<PR> _sqr(FloatError<PR> const& x) {
    return FloatError<PR>(mul_up(x._e,x._e));
}

template<class PR> FloatLowerBound<PR> _rec(FloatError<PR> const& x) {
    return FloatLowerBound<PR>(rec_down(x._e));
}

template<class PR> FloatError<PR> _add(FloatError<PR> const& x1, FloatError<PR> const& x2) {
    return FloatError<PR>(add_up(x1._e,x2._e));
}

template<class PR> FloatError<PR> _mul(FloatError<PR> const& x1, FloatError<PR> const& x2) {
    return FloatError<PR>(mul_up(x1._e,x2._e));
}

template<class PR> FloatError<PR> _pow(FloatError<PR> const& x1, Nat m2) {
    return FloatError<PR>(pow_up(x1._e,m2));
}

template<class PR> FloatUpperBound<PR> _log(FloatError<PR> const& x) {
    return FloatUpperBound<PR>(log_up(x._e));
}

template<class PR> FloatError<PR> _max(FloatError<PR> const& x1, FloatError<PR> const& x2) {
    return FloatError<PR>(max(x1._e,x2._e));
}

template<class PR> FloatError<PR> _min(FloatError<PR> const& x1, FloatError<PR> const& x2) {
    return FloatError<PR>(min(x1._e,x2._e));
}

template<class PR> Bool _same(FloatError<PR> const& x1, FloatError<PR> const& x2) {
    return x1._e == x2._e;
}

template<class PR> OutputStream& _write(OutputStream& os, FloatError<PR> const& x) {
    return write(os,x.raw(),FloatError<PR>::output_places,RawFloat<PR>::upward);
}

template<class PR> InputStream& _read(InputStream& is, FloatError<PR>& x) {
    FloatUpperBound<PR> xu; is >> xu; x=FloatError<PR>(xu);
}





template<class PR> Bool _refines(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return x1._l>=x2._l;
}

template<class PR> FloatLowerBound<PR> _refinement(FloatLowerBound<PR> const& x1, FloatLowerBound<PR> const& x2) {
    return FloatLowerBound<PR>(min(x1._l,x2._l));
}


template<class PR> Bool _refines(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return x1._u <= x2._u;
}

template<class PR> FloatUpperBound<PR> _refinement(FloatUpperBound<PR> const& x1, FloatUpperBound<PR> const& x2) {
    return FloatUpperBound<PR>(max(x1._u,x2._u));
}


template<class PR> Bool _models(FloatBounds<PR> const& x1, FloatValue<PR> const& x2) {
    return x1._l<=x2._v && x1._u >= x2._v;
}

template<class PR> Bool _consistent(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return x1._l<=x2._u && x1._u >= x2._l;
}

template<class PR> Bool _inconsistent(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return x1._l>x2._u || x1._u < x2._l;
}

template<class PR> Bool _refines(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return x1._l>=x2._l && x1._u <= x2._u;
}

template<class PR> FloatBounds<PR> _refinement(FloatBounds<PR> const& x1, FloatBounds<PR> const& x2) {
    return FloatBounds<PR>(max(x1._l,x2._l),min(x1._u,x2._u));
}


template<class PR> Bool _models(FloatBall<PR> const& x1, FloatValue<PR> const& x2) {
    return (x1._v>=x2._v ? sub_up(x1._v,x2._v) : sub_up(x2._v,x1._v)) <= x1._e;
}

template<class PR> Bool _consistent(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return _consistent(FloatBounds<PR>(x1),FloatBounds<PR>(x2));
}

template<class PR> Bool _inconsistent(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return _inconsistent(FloatBounds<PR>(x1),FloatBounds<PR>(x2));
}

template<class PR> Bool _refines(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return (x1._v>=x2._v ? sub_up(x1._v,x2._v) : sub_up(x2._v,x1._v)) <= sub_down(x2._e, x1._e);
}

template<class PR> FloatBall<PR> _refinement(FloatBall<PR> const& x1, FloatBall<PR> const& x2) {
    return FloatBall<PR>(_refinement(FloatBounds<PR>(x1),FloatBounds<PR>(x2)));
}


template<class PR> FloatError<PR> _refinement(FloatError<PR> const& x1, FloatError<PR> const& x2) {
    return FloatError<PR>(max(x1._e,x2._e));
}


#ifdef ARIADNE_ENABLE_SERIALIZATION
template<class PR, class A> Void serialize(A& _a, FloatBounds<PR>& x, const Nat version) {
    _a & x.lower_raw() & x.upper_raw(); }
#endif




template<class PR> Bool _same(PositiveFloatUpperBound<PR> const& x1, PositiveFloatUpperBound<PR> const& x2) {
    return x1._u == x2._u;
}

template<class PR> PositiveFloatUpperBound<PR> _max(PositiveFloatUpperBound<PR> const& x1, PositiveFloatUpperBound<PR> const& x2) {
    return PositiveFloatUpperBound<PR>(max(x1._u,x2._u));
}

template<class PR> PositiveFloatUpperBound<PR> _min(PositiveFloatUpperBound<PR> const& x1, PositiveFloatUpperBound<PR> const& x2) {
    return PositiveFloatUpperBound<PR>(min(x1._u,x2._u));
}

template<class PR> PositiveFloatUpperBound<PR> _abs(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(x._u);
}

template<class PR> PositiveFloatApproximation<PR> _mig(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatApproximation<PR>(x._u);
}

template<class PR> PositiveFloatUpperBound<PR> _mag(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(x._u);
}

template<class PR> PositiveFloatUpperBound<PR> _pos(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(pos(x._u));
}

template<class PR> FloatLowerBound<PR> __neg(PositiveFloatUpperBound<PR> const& x) {
    return FloatLowerBound<PR>(neg(x._u));
}

template<class PR> PositiveFloatUpperBound<PR> _rec(PositiveFloatLowerBound<PR> const& x) {
    return PositiveFloatLowerBound<PR>(rec_up(x._l));
}

template<class PR> PositiveFloatUpperBound<PR> _half(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(half(x._u));
}

template<class PR> PositiveFloatUpperBound<PR> _sqr(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatUpperBound<PR>(mul_up(x._u,x._u));
}

template<class PR> PositiveFloatUpperBound<PR> _add(PositiveFloatUpperBound<PR> const& x1, PositiveFloatUpperBound<PR> const& x2) {
    return PositiveFloatUpperBound<PR>(add_up(x1._u,x2._u));
}

template<class PR> PositiveFloatUpperBound<PR> _mul(PositiveFloatUpperBound<PR> const& x1, PositiveFloatUpperBound<PR> const& x2) {
    return PositiveFloatUpperBound<PR>(mul_up(x1._u,x2._u));
}

template<class PR> PositiveFloatUpperBound<PR> _div(PositiveFloatUpperBound<PR> const& x1, PositiveFloatLowerBound<PR> const& x2) {
    return PositiveFloatUpperBound<PR>(div_up(x1._u,x2._l));
}

template<class PR> PositiveFloatUpperBound<PR> _pow(PositiveFloatUpperBound<PR> const& x, Nat m) {
    return PositiveFloatUpperBound<PR>(pow_up(x._u,m));
}

template<class PR> PositiveFloatApproximation<PR> _pow(PositiveFloatUpperBound<PR> const& x, Int n) {
    return PositiveFloatApproximation<PR>(pow_approx(x._u,n));
}

template<class PR> FloatUpperBound<PR> _log(PositiveFloatUpperBound<PR> const& x) {
    return FloatUpperBound<PR>(log_up(x._u));
}

template<class PR> OutputStream& _write(OutputStream& os, PositiveFloatUpperBound<PR> const& x) {
    return os << static_cast<FloatUpperBound<PR>const&>(x);
}


template<class PR> PositiveFloatLowerBound<PR> _mig(PositiveFloatLowerBound<PR> const& x) {
    return PositiveFloatLowerBound<PR>(abs(x._l)); }

template<class PR> PositiveFloatLowerBound<PR> _rec(PositiveFloatUpperBound<PR> const& x) {
    return PositiveFloatLowerBound<PR>(rec_down(x._u));
}


template<class PR> PositiveFloatApproximation<PR> _mul(PositiveFloatApproximation<PR> const& x1, PositiveFloatApproximation<PR> const& x2) {
    return PositiveFloatApproximation<PR>(mul_near(x1._a,x2._a));
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

template<class PR> FloatApproximation<PR> _make_float(Number<ApproximateTag> x) { return FloatApproximation<PR>(x); }
template<class PR> FloatLowerBound<PR> _make_float(Number<ValidatedLowerTag> x) { return FloatLowerBound<PR>(x); }
template<class PR> FloatUpperBound<PR> _make_float(Number<ValidatedUpperTag> x) { return FloatUpperBound<PR>(x); }
template<class PR> FloatBounds<PR> _make_float(Number<ValidatedTag> x) { return FloatBounds<PR>(x); }
template<class PR> FloatBounds<PR> _make_float(Number<EffectiveTag> x) { return FloatBounds<PR>(x); }
template<class PR> FloatBounds<PR> _make_float(Number<ExactTag> x) { return FloatBounds<PR>(x); }
template<class PR> FloatBounds<PR> _make_float(Real r) { return FloatBounds<PR>(r); }
template<class PR> FloatBounds<PR> _make_float(Rational q) { return FloatBounds<PR>(q); }
template<class PR> FloatValue<PR> _make_float(Integer z) { return FloatValue<PR>(z); }

template class FloatApproximation<Precision64>;
template class FloatLowerBound<Precision64>;
template class FloatUpperBound<Precision64>;
template class FloatBounds<Precision64>;
template class FloatBall<Precision64>;
template class FloatValue<Precision64>;

template class FloatApproximation<PrecisionMP>;
template class FloatLowerBound<PrecisionMP>;
template class FloatUpperBound<PrecisionMP>;
template class FloatBounds<PrecisionMP>;
template class FloatBall<PrecisionMP>;
template class FloatValue<PrecisionMP>;


template<class T> class Dynamic { public: virtual ~Dynamic() = default; };
template<class T> String mangled_class_name() { Dynamic<T>* p; return typeid(p).name(); }



Float64Approximation operator+(Float64Approximation const& x) { return pos(x); }
Float64Approximation operator-(Float64Approximation const& x) { return neg(x); }
Float64Approximation operator+(Float64Approximation const& x1, Float64Approximation const& x2) { return add(x1,x2); }
Float64Approximation operator-(Float64Approximation const& x1, Float64Approximation const& x2) { return sub(x1,x2); }
Float64Approximation operator*(Float64Approximation const& x1, Float64Approximation const& x2) { return mul(x1,x2); }
Float64Approximation operator/(Float64Approximation const& x1, Float64Approximation const& x2) { return div(x1,x2); }
Float64Approximation& operator+=(Float64Approximation& x1, Float64Approximation const& x2) { return x1=x1+x2; }
Float64Approximation& operator-=(Float64Approximation& x1, Float64Approximation const& x2) { return x1=x1-x2; }
Float64Approximation& operator*=(Float64Approximation& x1, Float64Approximation const& x2) { return x1=x1*x2; }
Float64Approximation& operator/=(Float64Approximation& x1, Float64Approximation const& x2) { return x1=x1/x2; }

Float64Approximation operator+(Float64Approximation const& x1, Real const& y2) { return add(x1,y2); }
Float64Approximation operator+(Real const& y1, Float64Approximation const& x2) { return add(y1,x2); }
Float64Approximation operator-(Float64Approximation const& x1, Real const& y2) { return sub(x1,y2); }
Float64Approximation operator-(Real const& y1, Float64Approximation const& x2) { return sub(y1,x2); }
Float64Approximation operator*(Float64Approximation const& x1, Real const& y2) { return mul(x1,y2); }
Float64Approximation operator*(Real const& y1, Float64Approximation const& x2) { return mul(y1,x2); }
Float64Approximation operator/(Float64Approximation const& x1, Real const& y2) { return div(x1,y2); }
Float64Approximation operator/(Real const& y1, Float64Approximation const& x2) { return div(y1,x2); }

ApproximateKleenean operator==(Float64Approximation const& x1, Float64Approximation const& x2) { return _eq(x1,x2); }
ApproximateKleenean operator!=(Float64Approximation const& x1, Float64Approximation const& x2) { return not _eq(x1,x2); }
ApproximateKleenean operator<=(Float64Approximation const& x1, Float64Approximation const& x2) { return _leq(x1,x2); }
ApproximateKleenean operator>=(Float64Approximation const& x1, Float64Approximation const& x2) { return _leq(x2,x1); }
ApproximateKleenean operator< (Float64Approximation const& x1, Float64Approximation const& x2) { return not _leq(x2,x1); }
ApproximateKleenean operator> (Float64Approximation const& x1, Float64Approximation const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, Float64Approximation const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64Approximation& x) { return _read<Precision64>(is,x); }

Float64Approximation nul(Float64Approximation const& x) { return _nul<Precision64>(x); }
Float64Approximation pos(Float64Approximation const& x) { return _pos<Precision64>(x); }
Float64Approximation neg(Float64Approximation const& x) { return _neg<Precision64>(x); }
Float64Approximation half(Float64Approximation const& x) { return _half<Precision64>(x); }
Float64Approximation sqr(Float64Approximation const& x) { return _sqr<Precision64>(x); }
Float64Approximation rec(Float64Approximation const& x) { return _rec<Precision64>(x); }
Float64Approximation add(Float64Approximation const& x1, Float64Approximation const& x2) { return _add<Precision64>(x1,x2); }
Float64Approximation sub(Float64Approximation const& x1, Float64Approximation const& x2) { return _sub<Precision64>(x1,x2); }
Float64Approximation mul(Float64Approximation const& x1, Float64Approximation const& x2) { return _mul<Precision64>(x1,x2); }
Float64Approximation div(Float64Approximation const& x1, Float64Approximation const& x2) { return _div<Precision64>(x1,x2); }
Float64Approximation pow(Float64Approximation const& x1, Nat m2) { return _pow<Precision64>(x1,m2); }
Float64Approximation pow(Float64Approximation const& x1, Int n2) { return _pow<Precision64>(x1,n2); }
Float64Approximation sqrt(Float64Approximation const& x) { return _sqrt<Precision64>(x); }
Float64Approximation exp(Float64Approximation const& x) { return _exp<Precision64>(x); }
Float64Approximation log(Float64Approximation const& x) { return _log<Precision64>(x); }
Float64Approximation sin(Float64Approximation const& x) { return _sin<Precision64>(x); }
Float64Approximation cos(Float64Approximation const& x) { return _cos<Precision64>(x); }
Float64Approximation tan(Float64Approximation const& x) { return _tan<Precision64>(x); }
Float64Approximation asin(Float64Approximation const& x) { return _asin<Precision64>(x); }
Float64Approximation acos(Float64Approximation const& x) { return _acos<Precision64>(x); }
Float64Approximation atan(Float64Approximation const& x) { return _atan<Precision64>(x); }

Float64Approximation max(Float64Approximation const& x1, Float64Approximation const& x2) { return _max<Precision64>(x1,x2); }
Float64Approximation min(Float64Approximation const& x1, Float64Approximation const& x2) { return _min<Precision64>(x1,x2); }
Float64Approximation abs(Float64Approximation const& x) { return _abs<Precision64>(x); }
PositiveFloat64Approximation mag(Float64Approximation const& x) { return _mag<Precision64>(x); }
PositiveFloat64Approximation mig(Float64Approximation const& x) { return _mig<Precision64>(x); }
Float64Approximation round(Float64Approximation const& x) { return _round<Precision64>(x); }
Bool same(Float64Approximation const& x1, Float64Approximation const& x2) { return _same<Precision64>(x1,x2); }

FloatMPApproximation operator+(FloatMPApproximation const& x) { return pos(x); }
FloatMPApproximation operator-(FloatMPApproximation const& x) { return neg(x); }
FloatMPApproximation operator+(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return add(x1,x2); }
FloatMPApproximation operator-(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return sub(x1,x2); }
FloatMPApproximation operator*(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return mul(x1,x2); }
FloatMPApproximation operator/(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return div(x1,x2); }
FloatMPApproximation& operator+=(FloatMPApproximation& x1, FloatMPApproximation const& x2) { return x1=x1+x2; }
FloatMPApproximation& operator-=(FloatMPApproximation& x1, FloatMPApproximation const& x2) { return x1=x1-x2; }
FloatMPApproximation& operator*=(FloatMPApproximation& x1, FloatMPApproximation const& x2) { return x1=x1*x2; }
FloatMPApproximation& operator/=(FloatMPApproximation& x1, FloatMPApproximation const& x2) { return x1=x1/x2; }

FloatMPApproximation operator+(FloatMPApproximation const& x1, Real const& y2) { return add(x1,y2); }
FloatMPApproximation operator+(Real const& y1, FloatMPApproximation const& x2) { return add(y1,x2); }
FloatMPApproximation operator-(FloatMPApproximation const& x1, Real const& y2) { return sub(x1,y2); }
FloatMPApproximation operator-(Real const& y1, FloatMPApproximation const& x2) { return sub(y1,x2); }
FloatMPApproximation operator*(FloatMPApproximation const& x1, Real const& y2) { return mul(x1,y2); }
FloatMPApproximation operator*(Real const& y1, FloatMPApproximation const& x2) { return mul(y1,x2); }
FloatMPApproximation operator/(FloatMPApproximation const& x1, Real const& y2) { return div(x1,y2); }
FloatMPApproximation operator/(Real const& y1, FloatMPApproximation const& x2) { return div(y1,x2); }

ApproximateKleenean operator==(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _eq(x1,x2); }
ApproximateKleenean operator!=(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return not _eq(x1,x2); }
ApproximateKleenean operator<=(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _leq(x1,x2); }
ApproximateKleenean operator>=(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _leq(x2,x1); }
ApproximateKleenean operator< (FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return not _leq(x2,x1); }
ApproximateKleenean operator> (FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, FloatMPApproximation const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPApproximation& x) { return _read<PrecisionMP>(is,x); }

FloatMPApproximation nul(FloatMPApproximation const& x) { return _nul<PrecisionMP>(x); }
FloatMPApproximation pos(FloatMPApproximation const& x) { return _pos<PrecisionMP>(x); }
FloatMPApproximation neg(FloatMPApproximation const& x) { return _neg<PrecisionMP>(x); }
FloatMPApproximation half(FloatMPApproximation const& x) { return _half<PrecisionMP>(x); }
FloatMPApproximation sqr(FloatMPApproximation const& x) { return _sqr<PrecisionMP>(x); }
FloatMPApproximation rec(FloatMPApproximation const& x) { return _rec<PrecisionMP>(x); }
FloatMPApproximation add(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPApproximation sub(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _sub<PrecisionMP>(x1,x2); }
FloatMPApproximation mul(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _mul<PrecisionMP>(x1,x2); }
FloatMPApproximation div(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _div<PrecisionMP>(x1,x2); }
FloatMPApproximation pow(FloatMPApproximation const& x1, Nat m2) { return _pow<PrecisionMP>(x1,m2); }
FloatMPApproximation pow(FloatMPApproximation const& x1, Int n2) { return _pow<PrecisionMP>(x1,n2); }
FloatMPApproximation sqrt(FloatMPApproximation const& x) { return _sqrt<PrecisionMP>(x); }
FloatMPApproximation exp(FloatMPApproximation const& x) { return _exp<PrecisionMP>(x); }
FloatMPApproximation log(FloatMPApproximation const& x) { return _log<PrecisionMP>(x); }
FloatMPApproximation sin(FloatMPApproximation const& x) { return _sin<PrecisionMP>(x); }
FloatMPApproximation cos(FloatMPApproximation const& x) { return _cos<PrecisionMP>(x); }
FloatMPApproximation tan(FloatMPApproximation const& x) { return _tan<PrecisionMP>(x); }
FloatMPApproximation asin(FloatMPApproximation const& x) { return _asin<PrecisionMP>(x); }
FloatMPApproximation acos(FloatMPApproximation const& x) { return _acos<PrecisionMP>(x); }
FloatMPApproximation atan(FloatMPApproximation const& x) { return _atan<PrecisionMP>(x); }

FloatMPApproximation max(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPApproximation min(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _min<PrecisionMP>(x1,x2); }
FloatMPApproximation abs(FloatMPApproximation const& x) { return _abs<PrecisionMP>(x); }
PositiveFloatMPApproximation mag(FloatMPApproximation const& x) { return _mag<PrecisionMP>(x); }
PositiveFloatMPApproximation mig(FloatMPApproximation const& x) { return _mig<PrecisionMP>(x); }
FloatMPApproximation round(FloatMPApproximation const& x) { return _round<PrecisionMP>(x); }
Bool same(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return _same<PrecisionMP>(x1,x2); }



Float64LowerBound operator+(Float64LowerBound const& x) { return pos(x); }
Float64UpperBound operator-(Float64LowerBound const& x) { return neg(x); }
Float64LowerBound operator+(Float64LowerBound const& x1, Float64LowerBound const& x2) { return add(x1,x2); }
Float64LowerBound operator-(Float64LowerBound const& x1, Float64UpperBound const& x2) { return sub(x1,x2); }
Float64LowerBound& operator+=(Float64LowerBound& x1, Float64LowerBound const& x2) { return x1=x1+x2; }
Float64LowerBound& operator-=(Float64LowerBound& x1, Float64UpperBound const& x2) { return x1=x1-x2; }

Float64LowerBound operator+(Float64LowerBound const& x1, Real const& y2) { return add(x1,y2); }
Float64LowerBound operator+(Real const& y1, Float64LowerBound const& x2) { return add(y1,x2); }
Float64LowerBound operator-(Float64LowerBound const& x1, Real const& y2) { return sub(x1,y2); }
Float64LowerBound operator-(Real const& y1, Float64UpperBound const& x2) { return sub(y1,x2); }

ValidatedNegatedSierpinskian operator<=(Float64LowerBound const& x1, Float64UpperBound const& x2) { return _leq(x1,x2); }
ValidatedSierpinskian operator>=(Float64LowerBound const& x1, Float64UpperBound const& x2) { return _leq(x2,x1); }
ValidatedNegatedSierpinskian operator< (Float64LowerBound const& x1, Float64UpperBound const& x2) { return not _leq(x2,x1); }
ValidatedSierpinskian operator> (Float64LowerBound const& x1, Float64UpperBound const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, Float64LowerBound const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64LowerBound& x) { return _read<Precision64>(is,x); }

Float64LowerBound nul(Float64LowerBound const& x) { return _nul<Precision64>(x); }
Float64LowerBound pos(Float64LowerBound const& x) { return _pos<Precision64>(x); }
Float64UpperBound neg(Float64LowerBound const& x) { return _neg<Precision64>(x); }
Float64LowerBound half(Float64LowerBound const& x) { return _half<Precision64>(x); }
Float64LowerBound add(Float64LowerBound const& x1, Float64LowerBound const& x2) { return _add<Precision64>(x1,x2); }
Float64LowerBound sub(Float64LowerBound const& x1, Float64UpperBound const& x2) { return _sub<Precision64>(x1,x2); }
Float64LowerBound sqrt(Float64LowerBound const& x) { return _sqrt<Precision64>(x); }
Float64LowerBound exp(Float64LowerBound const& x) { return _exp<Precision64>(x); }
Float64LowerBound log(Float64LowerBound const& x) { return _log<Precision64>(x); }
Float64LowerBound atan(Float64LowerBound const& x) { return _atan<Precision64>(x); }

Float64LowerBound max(Float64LowerBound const& x1, Float64LowerBound const& x2) { return _max<Precision64>(x1,x2); }
Float64LowerBound min(Float64LowerBound const& x1, Float64LowerBound const& x2) { return _min<Precision64>(x1,x2); }

Bool same(Float64LowerBound const& x1, Float64LowerBound const& x2) { return _same<Precision64>(x1,x2); }
Bool refines(Float64LowerBound const& x1, Float64LowerBound const& x2) { return _refines<Precision64>(x1,x2); }
Float64LowerBound refinement(Float64LowerBound const& x1, Float64LowerBound const& x2) { return _refinement<Precision64>(x1,x2); }


FloatMPLowerBound operator+(FloatMPLowerBound const& x) { return pos(x); }
FloatMPUpperBound operator-(FloatMPLowerBound const& x) { return neg(x); }
FloatMPLowerBound operator+(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return add(x1,x2); }
FloatMPLowerBound operator-(FloatMPLowerBound const& x1, FloatMPUpperBound const& x2) { return sub(x1,x2); }
FloatMPLowerBound& operator+=(FloatMPLowerBound& x1, FloatMPLowerBound const& x2) { return x1=x1+x2; }
FloatMPLowerBound& operator-=(FloatMPLowerBound& x1, FloatMPUpperBound const& x2) { return x1=x1-x2; }

FloatMPLowerBound operator+(FloatMPLowerBound const& x1, Real const& y2) { return add(x1,y2); }
FloatMPLowerBound operator+(Real const& y1, FloatMPLowerBound const& x2) { return add(y1,x2); }
FloatMPLowerBound operator-(FloatMPLowerBound const& x1, Real const& y2) { return sub(x1,y2); }
FloatMPLowerBound operator-(Real const& y1, FloatMPUpperBound const& x2) { return sub(y1,x2); }

ValidatedNegatedSierpinskian operator<=(FloatMPLowerBound const& x1, FloatMPUpperBound const& x2) { return _leq(x1,x2); }
ValidatedSierpinskian operator>=(FloatMPLowerBound const& x1, FloatMPUpperBound const& x2) { return _leq(x2,x1); }
ValidatedNegatedSierpinskian operator< (FloatMPLowerBound const& x1, FloatMPUpperBound const& x2) { return not _leq(x2,x1); }
ValidatedSierpinskian operator> (FloatMPLowerBound const& x1, FloatMPUpperBound const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, FloatMPLowerBound const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPLowerBound& x) { return _read<PrecisionMP>(is,x); }

FloatMPLowerBound nul(FloatMPLowerBound const& x) { return _nul<PrecisionMP>(x); }
FloatMPLowerBound pos(FloatMPLowerBound const& x) { return _pos<PrecisionMP>(x); }
FloatMPUpperBound neg(FloatMPLowerBound const& x) { return _neg<PrecisionMP>(x); }
FloatMPLowerBound half(FloatMPLowerBound const& x) { return _half<PrecisionMP>(x); }
FloatMPLowerBound add(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPLowerBound sub(FloatMPLowerBound const& x1, FloatMPUpperBound const& x2) { return _sub<PrecisionMP>(x1,x2); }
FloatMPLowerBound sqrt(FloatMPLowerBound const& x) { return _sqrt<PrecisionMP>(x); }
FloatMPLowerBound exp(FloatMPLowerBound const& x) { return _exp<PrecisionMP>(x); }
FloatMPLowerBound log(FloatMPLowerBound const& x) { return _log<PrecisionMP>(x); }
FloatMPLowerBound atan(FloatMPLowerBound const& x) { return _atan<PrecisionMP>(x); }

FloatMPLowerBound max(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPLowerBound min(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return _min<PrecisionMP>(x1,x2); }

Bool same(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return _same<PrecisionMP>(x1,x2); }
Bool refines(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return _refines<PrecisionMP>(x1,x2); }
FloatMPLowerBound refinement(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return _refinement<PrecisionMP>(x1,x2); }



Float64UpperBound operator+(Float64UpperBound const& x) { return pos(x); }
Float64LowerBound operator-(Float64UpperBound const& x) { return neg(x); }
Float64UpperBound operator+(Float64UpperBound const& x1, Float64UpperBound const& x2) { return add(x1,x2); }
Float64UpperBound operator-(Float64UpperBound const& x1, Float64LowerBound const& x2) { return sub(x1,x2); }
Float64UpperBound& operator+=(Float64UpperBound& x1, Float64UpperBound const& x2) { return x1=x1+x2; }
Float64UpperBound& operator-=(Float64UpperBound& x1, Float64LowerBound const& x2) { return x1=x1-x2; }

Float64UpperBound operator+(Float64UpperBound const& x1, Real const& y2) { return add(x1,y2); }
Float64UpperBound operator+(Real const& y1, Float64UpperBound const& x2) { return add(y1,x2); }
Float64UpperBound operator-(Float64UpperBound const& x1, Real const& y2) { return sub(x1,y2); }
Float64UpperBound operator-(Real const& y1, Float64LowerBound const& x2) { return sub(y1,x2); }

ValidatedSierpinskian operator<=(Float64UpperBound const& x1, Float64LowerBound const& x2) { return _leq(x1,x2); }
ValidatedNegatedSierpinskian operator>=(Float64UpperBound const& x1, Float64LowerBound const& x2) { return _leq(x2,x1); }
ValidatedSierpinskian operator< (Float64UpperBound const& x1, Float64LowerBound const& x2) { return not _leq(x2,x1); }
ValidatedNegatedSierpinskian operator> (Float64UpperBound const& x1, Float64LowerBound const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, Float64UpperBound const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64UpperBound& x) { return _read<Precision64>(is,x); }

Float64UpperBound nul(Float64UpperBound const& x) { return _nul<Precision64>(x); }
Float64UpperBound pos(Float64UpperBound const& x) { return _pos<Precision64>(x); }
Float64LowerBound neg(Float64UpperBound const& x) { return _neg<Precision64>(x); }
Float64UpperBound half(Float64UpperBound const& x) { return _half<Precision64>(x); }
Float64UpperBound add(Float64UpperBound const& x1, Float64UpperBound const& x2) { return _add<Precision64>(x1,x2); }
Float64UpperBound sub(Float64UpperBound const& x1, Float64LowerBound const& x2) { return _sub<Precision64>(x1,x2); }
Float64UpperBound sqrt(Float64UpperBound const& x) { return _sqrt<Precision64>(x); }
Float64UpperBound exp(Float64UpperBound const& x) { return _exp<Precision64>(x); }
Float64UpperBound log(Float64UpperBound const& x) { return _log<Precision64>(x); }
Float64UpperBound atan(Float64UpperBound const& x) { return _atan<Precision64>(x); }

Float64UpperBound max(Float64UpperBound const& x1, Float64UpperBound const& x2) { return _max<Precision64>(x1,x2); }
Float64UpperBound min(Float64UpperBound const& x1, Float64UpperBound const& x2) { return _min<Precision64>(x1,x2); }

Bool same(Float64UpperBound const& x1, Float64UpperBound const& x2) { return _same<Precision64>(x1,x2); }
Bool refines(Float64UpperBound const& x1, Float64UpperBound const& x2) { return _refines<Precision64>(x1,x2); }
Float64UpperBound refinement(Float64UpperBound const& x1, Float64UpperBound const& x2) { return _refinement<Precision64>(x1,x2); }


FloatMPUpperBound operator+(FloatMPUpperBound const& x) { return pos(x); }
FloatMPLowerBound operator-(FloatMPUpperBound const& x) { return neg(x); }
FloatMPUpperBound operator+(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return add(x1,x2); }
FloatMPUpperBound operator-(FloatMPUpperBound const& x1, FloatMPLowerBound const& x2) { return sub(x1,x2); }
FloatMPUpperBound& operator+=(FloatMPUpperBound& x1, FloatMPUpperBound const& x2) { return x1=x1+x2; }
FloatMPUpperBound& operator-=(FloatMPUpperBound& x1, FloatMPLowerBound const& x2) { return x1=x1-x2; }

FloatMPUpperBound operator+(FloatMPUpperBound const& x1, Real const& y2) { return add(x1,y2); }
FloatMPUpperBound operator+(Real const& y1, FloatMPUpperBound const& x2) { return add(y1,x2); }
FloatMPUpperBound operator-(FloatMPUpperBound const& x1, Real const& y2) { return sub(x1,y2); }
FloatMPUpperBound operator-(Real const& y1, FloatMPLowerBound const& x2) { return sub(y1,x2); }

ValidatedSierpinskian operator<=(FloatMPUpperBound const& x1, FloatMPLowerBound const& x2) { return _leq(x1,x2); }
ValidatedNegatedSierpinskian operator>=(FloatMPUpperBound const& x1, FloatMPLowerBound const& x2) { return _leq(x2,x1); }
ValidatedSierpinskian operator< (FloatMPUpperBound const& x1, FloatMPLowerBound const& x2) { return not _leq(x2,x1); }
ValidatedNegatedSierpinskian operator> (FloatMPUpperBound const& x1, FloatMPLowerBound const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, FloatMPUpperBound const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPUpperBound& x) { return _read<PrecisionMP>(is,x); }

FloatMPUpperBound nul(FloatMPUpperBound const& x) { return _nul<PrecisionMP>(x); }
FloatMPUpperBound pos(FloatMPUpperBound const& x) { return _pos<PrecisionMP>(x); }
FloatMPLowerBound neg(FloatMPUpperBound const& x) { return _neg<PrecisionMP>(x); }
FloatMPUpperBound half(FloatMPUpperBound const& x) { return _half<PrecisionMP>(x); }
FloatMPUpperBound add(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPUpperBound sub(FloatMPUpperBound const& x1, FloatMPLowerBound const& x2) { return _sub<PrecisionMP>(x1,x2); }
FloatMPUpperBound sqrt(FloatMPUpperBound const& x) { return _sqrt<PrecisionMP>(x); }
FloatMPUpperBound exp(FloatMPUpperBound const& x) { return _exp<PrecisionMP>(x); }
FloatMPUpperBound log(FloatMPUpperBound const& x) { return _log<PrecisionMP>(x); }
FloatMPUpperBound atan(FloatMPUpperBound const& x) { return _atan<PrecisionMP>(x); }

FloatMPUpperBound max(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPUpperBound min(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return _min<PrecisionMP>(x1,x2); }

Bool same(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return _same<PrecisionMP>(x1,x2); }
Bool refines(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return _refines<PrecisionMP>(x1,x2); }
FloatMPUpperBound refinement(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return _refinement<PrecisionMP>(x1,x2); }



Float64Bounds operator+(Float64Bounds const& x) { return pos(x); }
Float64Bounds operator-(Float64Bounds const& x) { return neg(x); }
Float64Bounds operator+(Float64Bounds const& x1, Float64Bounds const& x2) { return add(x1,x2); }
Float64Bounds operator-(Float64Bounds const& x1, Float64Bounds const& x2) { return sub(x1,x2); }
Float64Bounds operator*(Float64Bounds const& x1, Float64Bounds const& x2) { return mul(x1,x2); }
Float64Bounds operator/(Float64Bounds const& x1, Float64Bounds const& x2) { return div(x1,x2); }
Float64Bounds& operator+=(Float64Bounds& x1, Float64Bounds const& x2) { return x1=x1+x2; }
Float64Bounds& operator-=(Float64Bounds& x1, Float64Bounds const& x2) { return x1=x1-x2; }
Float64Bounds& operator*=(Float64Bounds& x1, Float64Bounds const& x2) { return x1=x1*x2; }
Float64Bounds& operator/=(Float64Bounds& x1, Float64Bounds const& x2) { return x1=x1/x2; }

Float64Bounds operator+(Float64Bounds const& x1, Real const& y2) { return add(x1,y2); }
Float64Bounds operator+(Real const& y1, Float64Bounds const& x2) { return add(y1,x2); }
Float64Bounds operator-(Float64Bounds const& x1, Real const& y2) { return sub(x1,y2); }
Float64Bounds operator-(Real const& y1, Float64Bounds const& x2) { return sub(y1,x2); }
Float64Bounds operator*(Float64Bounds const& x1, Real const& y2) { return mul(x1,y2); }
Float64Bounds operator*(Real const& y1, Float64Bounds const& x2) { return mul(y1,x2); }
Float64Bounds operator/(Float64Bounds const& x1, Real const& y2) { return div(x1,y2); }
Float64Bounds operator/(Real const& y1, Float64Bounds const& x2) { return div(y1,x2); }

ValidatedKleenean operator==(Float64Bounds const& x1, Float64Bounds const& x2) { return _eq(x1,x2); }
ValidatedKleenean operator!=(Float64Bounds const& x1, Float64Bounds const& x2) { return not _eq(x1,x2); }
ValidatedKleenean operator<=(Float64Bounds const& x1, Float64Bounds const& x2) { return _leq(x1,x2); }
ValidatedKleenean operator>=(Float64Bounds const& x1, Float64Bounds const& x2) { return _leq(x2,x1); }
ValidatedKleenean operator< (Float64Bounds const& x1, Float64Bounds const& x2) { return not _leq(x2,x1); }
ValidatedKleenean operator> (Float64Bounds const& x1, Float64Bounds const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, Float64Bounds const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64Bounds& x) { return _read<Precision64>(is,x); }

Float64Bounds nul(Float64Bounds const& x) { return _nul<Precision64>(x); }
Float64Bounds pos(Float64Bounds const& x) { return _pos<Precision64>(x); }
Float64Bounds neg(Float64Bounds const& x) { return _neg<Precision64>(x); }
Float64Bounds half(Float64Bounds const& x) { return _half<Precision64>(x); }
Float64Bounds sqr(Float64Bounds const& x) { return _sqr<Precision64>(x); }
Float64Bounds rec(Float64Bounds const& x) { return _rec<Precision64>(x); }
Float64Bounds add(Float64Bounds const& x1, Float64Bounds const& x2) { return _add<Precision64>(x1,x2); }
Float64Bounds sub(Float64Bounds const& x1, Float64Bounds const& x2) { return _sub<Precision64>(x1,x2); }
Float64Bounds mul(Float64Bounds const& x1, Float64Bounds const& x2) { return _mul<Precision64>(x1,x2); }
Float64Bounds div(Float64Bounds const& x1, Float64Bounds const& x2) { return _div<Precision64>(x1,x2); }
Float64Bounds pow(Float64Bounds const& x1, Nat m2) { return _pow<Precision64>(x1,m2); }
Float64Bounds pow(Float64Bounds const& x1, Int n2) { return _pow<Precision64>(x1,n2); }
Float64Bounds sqrt(Float64Bounds const& x) { return _sqrt<Precision64>(x); }
Float64Bounds exp(Float64Bounds const& x) { return _exp<Precision64>(x); }
Float64Bounds log(Float64Bounds const& x) { return _log<Precision64>(x); }
Float64Bounds sin(Float64Bounds const& x) { return _sin<Precision64>(x); }
Float64Bounds cos(Float64Bounds const& x) { return _cos<Precision64>(x); }
Float64Bounds tan(Float64Bounds const& x) { return _tan<Precision64>(x); }
Float64Bounds asin(Float64Bounds const& x) { return _asin<Precision64>(x); }
Float64Bounds acos(Float64Bounds const& x) { return _acos<Precision64>(x); }
Float64Bounds atan(Float64Bounds const& x) { return _atan<Precision64>(x); }

Float64Bounds max(Float64Bounds const& x1, Float64Bounds const& x2) { return _max<Precision64>(x1,x2); }
Float64Bounds min(Float64Bounds const& x1, Float64Bounds const& x2) { return _min<Precision64>(x1,x2); }
Float64Bounds abs(Float64Bounds const& x) { return _abs<Precision64>(x); }
Float64Error mag(Float64Bounds const& x) { return _mag<Precision64>(x); }
Float64LowerBound mig(Float64Bounds const& x) { return _mig<Precision64>(x); }
Float64Bounds round(Float64Bounds const& x) { return _round<Precision64>(x); }

Bool same(Float64Bounds const& x1, Float64Bounds const& x2) { return _same<Precision64>(x1,x2); }
Bool models(Float64Bounds const& x1, Float64Value const& x2) { return _models<Precision64>(x1,x2); }
Bool refines(Float64Bounds const& x1, Float64Bounds const& x2) { return _refines<Precision64>(x1,x2); }
Bool consistent(Float64Bounds const& x1, Float64Bounds const& x2) { return _consistent<Precision64>(x1,x2); }
Bool inconsistent(Float64Bounds const& x1, Float64Bounds const& x2) { return _inconsistent<Precision64>(x1,x2); }
Float64Bounds refinement(Float64Bounds const& x1, Float64Bounds const& x2) { return _refinement<Precision64>(x1,x2); }


FloatMPBounds operator+(FloatMPBounds const& x) { return pos(x); }
FloatMPBounds operator-(FloatMPBounds const& x) { return neg(x); }
FloatMPBounds operator+(FloatMPBounds const& x1, FloatMPBounds const& x2) { return add(x1,x2); }
FloatMPBounds operator-(FloatMPBounds const& x1, FloatMPBounds const& x2) { return sub(x1,x2); }
FloatMPBounds operator*(FloatMPBounds const& x1, FloatMPBounds const& x2) { return mul(x1,x2); }
FloatMPBounds operator/(FloatMPBounds const& x1, FloatMPBounds const& x2) { return div(x1,x2); }
FloatMPBounds& operator+=(FloatMPBounds& x1, FloatMPBounds const& x2) { return x1=x1+x2; }
FloatMPBounds& operator-=(FloatMPBounds& x1, FloatMPBounds const& x2) { return x1=x1-x2; }
FloatMPBounds& operator*=(FloatMPBounds& x1, FloatMPBounds const& x2) { return x1=x1*x2; }
FloatMPBounds& operator/=(FloatMPBounds& x1, FloatMPBounds const& x2) { return x1=x1/x2; }

FloatMPBounds operator+(FloatMPBounds const& x1, Real const& y2) { return add(x1,y2); }
FloatMPBounds operator+(Real const& y1, FloatMPBounds const& x2) { return add(y1,x2); }
FloatMPBounds operator-(FloatMPBounds const& x1, Real const& y2) { return sub(x1,y2); }
FloatMPBounds operator-(Real const& y1, FloatMPBounds const& x2) { return sub(y1,x2); }
FloatMPBounds operator*(FloatMPBounds const& x1, Real const& y2) { return mul(x1,y2); }
FloatMPBounds operator*(Real const& y1, FloatMPBounds const& x2) { return mul(y1,x2); }
FloatMPBounds operator/(FloatMPBounds const& x1, Real const& y2) { return div(x1,y2); }
FloatMPBounds operator/(Real const& y1, FloatMPBounds const& x2) { return div(y1,x2); }

OutputStream& operator<<(OutputStream& os, FloatMPBounds const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPBounds& x) { return _read<PrecisionMP>(is,x); }

ValidatedKleenean operator==(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _eq(x1,x2); }
ValidatedKleenean operator!=(FloatMPBounds const& x1, FloatMPBounds const& x2) { return not _eq(x1,x2); }
ValidatedKleenean operator<=(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _leq(x1,x2); }
ValidatedKleenean operator>=(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _leq(x2,x1); }
ValidatedKleenean operator< (FloatMPBounds const& x1, FloatMPBounds const& x2) { return not _leq(x2,x1); }
ValidatedKleenean operator> (FloatMPBounds const& x1, FloatMPBounds const& x2) { return not _leq(x1,x2); }

FloatMPBounds nul(FloatMPBounds const& x) { return _nul<PrecisionMP>(x); }
FloatMPBounds pos(FloatMPBounds const& x) { return _pos<PrecisionMP>(x); }
FloatMPBounds neg(FloatMPBounds const& x) { return _neg<PrecisionMP>(x); }
FloatMPBounds half(FloatMPBounds const& x) { return _half<PrecisionMP>(x); }
FloatMPBounds sqr(FloatMPBounds const& x) { return _sqr<PrecisionMP>(x); }
FloatMPBounds rec(FloatMPBounds const& x) { return _rec<PrecisionMP>(x); }
FloatMPBounds add(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPBounds sub(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _sub<PrecisionMP>(x1,x2); }
FloatMPBounds mul(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _mul<PrecisionMP>(x1,x2); }
FloatMPBounds div(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _div<PrecisionMP>(x1,x2); }
FloatMPBounds pow(FloatMPBounds const& x1, Nat m2) { return _pow<PrecisionMP>(x1,m2); }
FloatMPBounds pow(FloatMPBounds const& x1, Int n2) { return _pow<PrecisionMP>(x1,n2); }
FloatMPBounds sqrt(FloatMPBounds const& x) { return _sqrt<PrecisionMP>(x); }
FloatMPBounds exp(FloatMPBounds const& x) { return _exp<PrecisionMP>(x); }
FloatMPBounds log(FloatMPBounds const& x) { return _log<PrecisionMP>(x); }
FloatMPBounds sin(FloatMPBounds const& x) { return _sin<PrecisionMP>(x); }
FloatMPBounds cos(FloatMPBounds const& x) { return _cos<PrecisionMP>(x); }
FloatMPBounds tan(FloatMPBounds const& x) { return _tan<PrecisionMP>(x); }
FloatMPBounds asin(FloatMPBounds const& x) { return _asin<PrecisionMP>(x); }
FloatMPBounds acos(FloatMPBounds const& x) { return _acos<PrecisionMP>(x); }
FloatMPBounds atan(FloatMPBounds const& x) { return _atan<PrecisionMP>(x); }

FloatMPBounds max(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPBounds min(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _min<PrecisionMP>(x1,x2); }
FloatMPBounds abs(FloatMPBounds const& x) { return _abs<PrecisionMP>(x); }
FloatMPError mag(FloatMPBounds const& x) { return _mag<PrecisionMP>(x); }
FloatMPLowerBound mig(FloatMPBounds const& x) { return _mig<PrecisionMP>(x); }
FloatMPBounds round(FloatMPBounds const& x) { return _round<PrecisionMP>(x); }

Bool same(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _same<PrecisionMP>(x1,x2); }
Bool models(FloatMPBounds const& x1, FloatMPValue const& x2) { return _models<PrecisionMP>(x1,x2); }
Bool refines(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _refines<PrecisionMP>(x1,x2); }
Bool consistent(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _consistent<PrecisionMP>(x1,x2); }
Bool inconsistent(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _inconsistent<PrecisionMP>(x1,x2); }
FloatMPBounds refinement(FloatMPBounds const& x1, FloatMPBounds const& x2) { return _refinement<PrecisionMP>(x1,x2); }


Float64Ball operator+(Float64Ball const& x) { return pos(x); }
Float64Ball operator-(Float64Ball const& x) { return neg(x); }
Float64Ball operator+(Float64Ball const& x1, Float64Ball const& x2) { return add(x1,x2); }
Float64Ball operator-(Float64Ball const& x1, Float64Ball const& x2) { return sub(x1,x2); }
Float64Ball operator*(Float64Ball const& x1, Float64Ball const& x2) { return mul(x1,x2); }
Float64Ball operator/(Float64Ball const& x1, Float64Ball const& x2) { return div(x1,x2); }
Float64Ball& operator+=(Float64Ball& x1, Float64Ball const& x2) { return x1=x1+x2; }
Float64Ball& operator-=(Float64Ball& x1, Float64Ball const& x2) { return x1=x1-x2; }
Float64Ball& operator*=(Float64Ball& x1, Float64Ball const& x2) { return x1=x1*x2; }
Float64Ball& operator/=(Float64Ball& x1, Float64Ball const& x2) { return x1=x1/x2; }

Float64Ball operator+(Float64Ball const& x1, Real const& y2) { return add(x1,y2); }
Float64Ball operator+(Real const& y1, Float64Ball const& x2) { return add(y1,x2); }
Float64Ball operator-(Float64Ball const& x1, Real const& y2) { return sub(x1,y2); }
Float64Ball operator-(Real const& y1, Float64Ball const& x2) { return sub(y1,x2); }
Float64Ball operator*(Float64Ball const& x1, Real const& y2) { return mul(x1,y2); }
Float64Ball operator*(Real const& y1, Float64Ball const& x2) { return mul(y1,x2); }
Float64Ball operator/(Float64Ball const& x1, Real const& y2) { return div(x1,y2); }
Float64Ball operator/(Real const& y1, Float64Ball const& x2) { return div(y1,x2); }

ValidatedKleenean operator==(Float64Ball const& x1, Float64Ball const& x2) { return _eq(x1,x2); }
ValidatedKleenean operator!=(Float64Ball const& x1, Float64Ball const& x2) { return not _eq(x1,x2); }
ValidatedKleenean operator<=(Float64Ball const& x1, Float64Ball const& x2) { return _leq(x1,x2); }
ValidatedKleenean operator>=(Float64Ball const& x1, Float64Ball const& x2) { return _leq(x2,x1); }
ValidatedKleenean operator< (Float64Ball const& x1, Float64Ball const& x2) { return not _leq(x2,x1); }
ValidatedKleenean operator> (Float64Ball const& x1, Float64Ball const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, Float64Ball const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64Ball& x) { return _read<Precision64>(is,x); }

Float64Ball nul(Float64Ball const& x) { return _nul<Precision64>(x); }
Float64Ball pos(Float64Ball const& x) { return _pos<Precision64>(x); }
Float64Ball neg(Float64Ball const& x) { return _neg<Precision64>(x); }
Float64Ball half(Float64Ball const& x) { return _half<Precision64>(x); }
Float64Ball sqr(Float64Ball const& x) { return _sqr<Precision64>(x); }
Float64Ball rec(Float64Ball const& x) { return _rec<Precision64>(x); }
Float64Ball add(Float64Ball const& x1, Float64Ball const& x2) { return _add<Precision64>(x1,x2); }
Float64Ball sub(Float64Ball const& x1, Float64Ball const& x2) { return _sub<Precision64>(x1,x2); }
Float64Ball mul(Float64Ball const& x1, Float64Ball const& x2) { return _mul<Precision64>(x1,x2); }
Float64Ball div(Float64Ball const& x1, Float64Ball const& x2) { return _div<Precision64>(x1,x2); }
Float64Ball pow(Float64Ball const& x1, Nat m2) { return _pow<Precision64>(x1,m2); }
Float64Ball pow(Float64Ball const& x1, Int n2) { return _pow<Precision64>(x1,n2); }
Float64Ball sqrt(Float64Ball const& x) { return _sqrt<Precision64>(x); }
Float64Ball exp(Float64Ball const& x) { return _exp<Precision64>(x); }
Float64Ball log(Float64Ball const& x) { return _log<Precision64>(x); }
Float64Ball sin(Float64Ball const& x) { return _sin<Precision64>(x); }
Float64Ball cos(Float64Ball const& x) { return _cos<Precision64>(x); }
Float64Ball tan(Float64Ball const& x) { return _tan<Precision64>(x); }
Float64Ball asin(Float64Ball const& x) { return _asin<Precision64>(x); }
Float64Ball acos(Float64Ball const& x) { return _acos<Precision64>(x); }
Float64Ball atan(Float64Ball const& x) { return _atan<Precision64>(x); }

Float64Ball max(Float64Ball const& x1, Float64Ball const& x2) { return _max<Precision64>(x1,x2); }
Float64Ball min(Float64Ball const& x1, Float64Ball const& x2) { return _min<Precision64>(x1,x2); }
Float64Ball abs(Float64Ball const& x) { return _abs<Precision64>(x); }
Float64Error mag(Float64Ball const& x) { return _mag<Precision64>(x); }
Float64LowerBound mig(Float64Ball const& x) { return _mig<Precision64>(x); }

Bool same(Float64Ball const& x1, Float64Ball const& x2) { return _same<Precision64>(x1,x2); }
Bool models(Float64Ball const& x1, Float64Value const& x2) { return _models<Precision64>(x1,x2); }
Bool refines(Float64Ball const& x1, Float64Ball const& x2) { return _refines<Precision64>(x1,x2); }
Bool consistent(Float64Ball const& x1, Float64Ball const& x2) { return _consistent<Precision64>(x1,x2); }
Bool inconsistent(Float64Ball const& x1, Float64Ball const& x2) { return _inconsistent<Precision64>(x1,x2); }
Float64Ball refinement(Float64Ball const& x1, Float64Ball const& x2) { return _refinement<Precision64>(x1,x2); }


FloatMPBall operator+(FloatMPBall const& x) { return pos(x); }
FloatMPBall operator-(FloatMPBall const& x) { return neg(x); }
FloatMPBall operator+(FloatMPBall const& x1, FloatMPBall const& x2) { return add(x1,x2); }
FloatMPBall operator-(FloatMPBall const& x1, FloatMPBall const& x2) { return sub(x1,x2); }
FloatMPBall operator*(FloatMPBall const& x1, FloatMPBall const& x2) { return mul(x1,x2); }
FloatMPBall operator/(FloatMPBall const& x1, FloatMPBall const& x2) { return div(x1,x2); }
FloatMPBall& operator+=(FloatMPBall& x1, FloatMPBall const& x2) { return x1=x1+x2; }
FloatMPBall& operator-=(FloatMPBall& x1, FloatMPBall const& x2) { return x1=x1-x2; }
FloatMPBall& operator*=(FloatMPBall& x1, FloatMPBall const& x2) { return x1=x1*x2; }
FloatMPBall& operator/=(FloatMPBall& x1, FloatMPBall const& x2) { return x1=x1/x2; }

FloatMPBall operator+(FloatMPBall const& x1, Real const& y2) { return add(x1,y2); }
FloatMPBall operator+(Real const& y1, FloatMPBall const& x2) { return add(y1,x2); }
FloatMPBall operator-(FloatMPBall const& x1, Real const& y2) { return sub(x1,y2); }
FloatMPBall operator-(Real const& y1, FloatMPBall const& x2) { return sub(y1,x2); }
FloatMPBall operator*(FloatMPBall const& x1, Real const& y2) { return mul(x1,y2); }
FloatMPBall operator*(Real const& y1, FloatMPBall const& x2) { return mul(y1,x2); }
FloatMPBall operator/(FloatMPBall const& x1, Real const& y2) { return div(x1,y2); }
FloatMPBall operator/(Real const& y1, FloatMPBall const& x2) { return div(y1,x2); }

OutputStream& operator<<(OutputStream& os, FloatMPBall const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPBall& x) { return _read<PrecisionMP>(is,x); }

ValidatedKleenean operator==(FloatMPBall const& x1, FloatMPBall const& x2) { return _eq(x1,x2); }
ValidatedKleenean operator!=(FloatMPBall const& x1, FloatMPBall const& x2) { return not _eq(x1,x2); }
ValidatedKleenean operator<=(FloatMPBall const& x1, FloatMPBall const& x2) { return _leq(x1,x2); }
ValidatedKleenean operator>=(FloatMPBall const& x1, FloatMPBall const& x2) { return _leq(x2,x1); }
ValidatedKleenean operator< (FloatMPBall const& x1, FloatMPBall const& x2) { return not _leq(x2,x1); }
ValidatedKleenean operator> (FloatMPBall const& x1, FloatMPBall const& x2) { return not _leq(x1,x2); }

FloatMPBall nul(FloatMPBall const& x) { return _nul<PrecisionMP>(x); }
FloatMPBall pos(FloatMPBall const& x) { return _pos<PrecisionMP>(x); }
FloatMPBall neg(FloatMPBall const& x) { return _neg<PrecisionMP>(x); }
FloatMPBall half(FloatMPBall const& x) { return _half<PrecisionMP>(x); }
FloatMPBall sqr(FloatMPBall const& x) { return _sqr<PrecisionMP>(x); }
FloatMPBall rec(FloatMPBall const& x) { return _rec<PrecisionMP>(x); }
FloatMPBall add(FloatMPBall const& x1, FloatMPBall const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPBall sub(FloatMPBall const& x1, FloatMPBall const& x2) { return _sub<PrecisionMP>(x1,x2); }
FloatMPBall mul(FloatMPBall const& x1, FloatMPBall const& x2) { return _mul<PrecisionMP>(x1,x2); }
FloatMPBall div(FloatMPBall const& x1, FloatMPBall const& x2) { return _div<PrecisionMP>(x1,x2); }
FloatMPBall pow(FloatMPBall const& x1, Nat m2) { return _pow<PrecisionMP>(x1,m2); }
FloatMPBall pow(FloatMPBall const& x1, Int n2) { return _pow<PrecisionMP>(x1,n2); }
FloatMPBall sqrt(FloatMPBall const& x) { return _sqrt<PrecisionMP>(x); }
FloatMPBall exp(FloatMPBall const& x) { return _exp<PrecisionMP>(x); }
FloatMPBall log(FloatMPBall const& x) { return _log<PrecisionMP>(x); }
FloatMPBall sin(FloatMPBall const& x) { return _sin<PrecisionMP>(x); }
FloatMPBall cos(FloatMPBall const& x) { return _cos<PrecisionMP>(x); }
FloatMPBall tan(FloatMPBall const& x) { return _tan<PrecisionMP>(x); }
FloatMPBall asin(FloatMPBall const& x) { return _asin<PrecisionMP>(x); }
FloatMPBall acos(FloatMPBall const& x) { return _acos<PrecisionMP>(x); }
FloatMPBall atan(FloatMPBall const& x) { return _atan<PrecisionMP>(x); }

FloatMPBall max(FloatMPBall const& x1, FloatMPBall const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPBall min(FloatMPBall const& x1, FloatMPBall const& x2) { return _min<PrecisionMP>(x1,x2); }
FloatMPBall abs(FloatMPBall const& x) { return _abs<PrecisionMP>(x); }
FloatMPError mag(FloatMPBall const& x) { return _mag<PrecisionMP>(x); }
FloatMPLowerBound mig(FloatMPBall const& x) { return _mig<PrecisionMP>(x); }

Bool same(FloatMPBall const& x1, FloatMPBall const& x2) { return _same<PrecisionMP>(x1,x2); }
Bool models(FloatMPBall const& x1, FloatMPValue const& x2) { return _models<PrecisionMP>(x1,x2); }
Bool refines(FloatMPBall const& x1, FloatMPBall const& x2) { return _refines<PrecisionMP>(x1,x2); }
Bool consistent(FloatMPBall const& x1, FloatMPBall const& x2) { return _consistent<PrecisionMP>(x1,x2); }
Bool inconsistent(FloatMPBall const& x1, FloatMPBall const& x2) { return _inconsistent<PrecisionMP>(x1,x2); }
FloatMPBall refinement(FloatMPBall const& x1, FloatMPBall const& x2) { return _refinement<PrecisionMP>(x1,x2); }



Float64Value operator+(Float64Value const& x) { return pos(x); }
Float64Value operator-(Float64Value const& x) { return neg(x); }
Float64Bounds operator+(Float64Value const& x1, Float64Value const& x2) { return add(x1,x2); }
Float64Bounds operator-(Float64Value const& x1, Float64Value const& x2) { return sub(x1,x2); }
Float64Bounds operator*(Float64Value const& x1, Float64Value const& x2) { return mul(x1,x2); }
Float64Bounds operator/(Float64Value const& x1, Float64Value const& x2) { return div(x1,x2); }

Float64Bounds operator+(Float64Value const& x1, Real const& y2) { return add(x1,y2); }
Float64Bounds operator+(Real const& y1, Float64Value const& x2) { return add(y1,x2); }
Float64Bounds operator-(Float64Value const& x1, Real const& y2) { return sub(x1,y2); }
Float64Bounds operator-(Real const& y1, Float64Value const& x2) { return sub(y1,x2); }
Float64Bounds operator*(Float64Value const& x1, Real const& y2) { return mul(x1,y2); }
Float64Bounds operator*(Real const& y1, Float64Value const& x2) { return mul(y1,x2); }
Float64Bounds operator/(Float64Value const& x1, Real const& y2) { return div(x1,y2); }
Float64Bounds operator/(Real const& y1, Float64Value const& x2) { return div(y1,x2); }

Boolean operator==(Float64Value const& x1, Float64Value const& x2) { return _eq(x1,x2); }
Boolean operator!=(Float64Value const& x1, Float64Value const& x2) { return not _eq(x1,x2); }
Boolean operator<=(Float64Value const& x1, Float64Value const& x2) { return _leq(x1,x2); }
Boolean operator>=(Float64Value const& x1, Float64Value const& x2) { return _leq(x2,x1); }
Boolean operator< (Float64Value const& x1, Float64Value const& x2) { return not _leq(x2,x1); }
Boolean operator> (Float64Value const& x1, Float64Value const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, Float64Value const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64Value& x) { return _read<Precision64>(is,x); }

Float64Value nul(Float64Value const& x) { return _nul<Precision64>(x); }
Float64Value pos(Float64Value const& x) { return _pos<Precision64>(x); }
Float64Value neg(Float64Value const& x) { return _neg<Precision64>(x); }
Float64Value half(Float64Value const& x) { return _half<Precision64>(x); }
Float64Bounds sqr(Float64Value const& x) { return _sqr<Precision64>(x); }
Float64Bounds rec(Float64Value const& x) { return _rec<Precision64>(x); }
Float64Bounds add(Float64Value const& x1, Float64Value const& x2) { return _add<Precision64>(x1,x2); }
Float64Bounds sub(Float64Value const& x1, Float64Value const& x2) { return _sub<Precision64>(x1,x2); }
Float64Bounds mul(Float64Value const& x1, Float64Value const& x2) { return _mul<Precision64>(x1,x2); }
Float64Bounds div(Float64Value const& x1, Float64Value const& x2) { return _div<Precision64>(x1,x2); }
Float64Bounds pow(Float64Value const& x1, Nat m2) { return _pow<Precision64>(x1,m2); }
Float64Bounds pow(Float64Value const& x1, Int n2) { return _pow<Precision64>(x1,n2); }
Float64Bounds sqrt(Float64Value const& x) { return _sqrt<Precision64>(x); }
Float64Bounds exp(Float64Value const& x) { return _exp<Precision64>(x); }
Float64Bounds log(Float64Value const& x) { return _log<Precision64>(x); }
Float64Bounds sin(Float64Value const& x) { return _sin<Precision64>(x); }
Float64Bounds cos(Float64Value const& x) { return _cos<Precision64>(x); }
Float64Bounds tan(Float64Value const& x) { return _tan<Precision64>(x); }
Float64Bounds asin(Float64Value const& x) { return _asin<Precision64>(x); }
Float64Bounds acos(Float64Value const& x) { return _acos<Precision64>(x); }
Float64Bounds atan(Float64Value const& x) { return _atan<Precision64>(x); }

Float64Value max(Float64Value const& x1, Float64Value const& x2) { return _max<Precision64>(x1,x2); }
Float64Value min(Float64Value const& x1, Float64Value const& x2) { return _min<Precision64>(x1,x2); }
Float64Value abs(Float64Value const& x) { return _abs<Precision64>(x); }
Float64Error mag(Float64Value const& x) { return _mag<Precision64>(x); }


FloatMPValue operator+(FloatMPValue const& x) { return pos(x); }
FloatMPValue operator-(FloatMPValue const& x) { return neg(x); }
FloatMPBounds operator+(FloatMPValue const& x1, FloatMPValue const& x2) { return add(x1,x2); }
FloatMPBounds operator-(FloatMPValue const& x1, FloatMPValue const& x2) { return sub(x1,x2); }
FloatMPBounds operator*(FloatMPValue const& x1, FloatMPValue const& x2) { return mul(x1,x2); }
FloatMPBounds operator/(FloatMPValue const& x1, FloatMPValue const& x2) { return div(x1,x2); }

FloatMPBounds operator+(FloatMPValue const& x1, Real const& y2) { return add(x1,y2); }
FloatMPBounds operator+(Real const& y1, FloatMPValue const& x2) { return add(y1,x2); }
FloatMPBounds operator-(FloatMPValue const& x1, Real const& y2) { return sub(x1,y2); }
FloatMPBounds operator-(Real const& y1, FloatMPValue const& x2) { return sub(y1,x2); }
FloatMPBounds operator*(FloatMPValue const& x1, Real const& y2) { return mul(x1,y2); }
FloatMPBounds operator*(Real const& y1, FloatMPValue const& x2) { return mul(y1,x2); }
FloatMPBounds operator/(FloatMPValue const& x1, Real const& y2) { return div(x1,y2); }
FloatMPBounds operator/(Real const& y1, FloatMPValue const& x2) { return div(y1,x2); }

Boolean operator==(FloatMPValue const& x1, FloatMPValue const& x2) { return _eq(x1,x2); }
Boolean operator!=(FloatMPValue const& x1, FloatMPValue const& x2) { return not _eq(x1,x2); }
Boolean operator<=(FloatMPValue const& x1, FloatMPValue const& x2) { return _leq(x1,x2); }
Boolean operator>=(FloatMPValue const& x1, FloatMPValue const& x2) { return _leq(x2,x1); }
Boolean operator< (FloatMPValue const& x1, FloatMPValue const& x2) { return not _leq(x2,x1); }
Boolean operator> (FloatMPValue const& x1, FloatMPValue const& x2) { return not _leq(x1,x2); }

OutputStream& operator<<(OutputStream& os, FloatMPValue const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPValue& x) { return _read<PrecisionMP>(is,x); }

FloatMPValue nul(FloatMPValue const& x) { return _nul<PrecisionMP>(x); }
FloatMPValue pos(FloatMPValue const& x) { return _pos<PrecisionMP>(x); }
FloatMPValue neg(FloatMPValue const& x) { return _neg<PrecisionMP>(x); }
FloatMPValue half(FloatMPValue const& x) { return _half<PrecisionMP>(x); }
FloatMPBounds sqr(FloatMPValue const& x) { return _sqr<PrecisionMP>(x); }
FloatMPBounds rec(FloatMPValue const& x) { return _rec<PrecisionMP>(x); }
FloatMPBounds add(FloatMPValue const& x1, FloatMPValue const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPBounds sub(FloatMPValue const& x1, FloatMPValue const& x2) { return _sub<PrecisionMP>(x1,x2); }
FloatMPBounds mul(FloatMPValue const& x1, FloatMPValue const& x2) { return _mul<PrecisionMP>(x1,x2); }
FloatMPBounds div(FloatMPValue const& x1, FloatMPValue const& x2) { return _div<PrecisionMP>(x1,x2); }
FloatMPBounds pow(FloatMPValue const& x1, Nat m2) { return _pow<PrecisionMP>(x1,m2); }
FloatMPBounds pow(FloatMPValue const& x1, Int n2) { return _pow<PrecisionMP>(x1,n2); }
FloatMPBounds sqrt(FloatMPValue const& x) { return _sqrt<PrecisionMP>(x); }
FloatMPBounds exp(FloatMPValue const& x) { return _exp<PrecisionMP>(x); }
FloatMPBounds log(FloatMPValue const& x) { return _log<PrecisionMP>(x); }
FloatMPBounds sin(FloatMPValue const& x) { return _sin<PrecisionMP>(x); }
FloatMPBounds cos(FloatMPValue const& x) { return _cos<PrecisionMP>(x); }
FloatMPBounds tan(FloatMPValue const& x) { return _tan<PrecisionMP>(x); }
FloatMPBounds asin(FloatMPValue const& x) { return _asin<PrecisionMP>(x); }
FloatMPBounds acos(FloatMPValue const& x) { return _acos<PrecisionMP>(x); }
FloatMPBounds atan(FloatMPValue const& x) { return _atan<PrecisionMP>(x); }

FloatMPValue max(FloatMPValue const& x1, FloatMPValue const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPValue min(FloatMPValue const& x1, FloatMPValue const& x2) { return _min<PrecisionMP>(x1,x2); }
FloatMPValue abs(FloatMPValue const& x) { return _abs<PrecisionMP>(x); }
FloatMPError mag(FloatMPValue const& x) { return _mag<PrecisionMP>(x); }

Float64Value operator+(TwoExp y) { return Float64Value(y); }
Float64Value operator-(TwoExp y) { return neg(Float64Value(y)); }
Float64Value operator*(Float64Value const& x1, TwoExp y2) { Float64Value x2(y2); return Float64Value(mul_near(x1._v,x2._v)); }
Float64Value operator/(Float64Value const& x1, TwoExp y2) { Float64Value x2(y2); return Float64Value(div_near(x1._v,x2._v)); }
Float64Value& operator*=(Float64Value& x1, TwoExp y2) { return x1=x1*y2; }
Float64Value& operator/=(Float64Value& x1, TwoExp y2) { return x1=x1/y2; }


Float64UpperBound operator+(Float64Error const& x) { return x; }
Float64LowerBound operator-(Float64Error const& x) { return -(+x); }
Float64Error operator+(Float64Error const& x1, Float64Error const& x2) { return _add<Precision64>(x1,x2); }
Float64Error operator*(Float64Error const& x1, Float64Error const& x2) { return _mul<Precision64>(x1,x2); }
Float64Error operator*(Nat y1, Float64Error const& x2) { return Float64Error(y1)*x2; }
Float64Error& operator+=(Float64Error& x1, Float64Error const& x2) { return x1=x1+x2; }
Float64Error& operator*=(Float64Error& x1, Float64Error const& x2) { return x1=x1*x2; }
Float64Error sqr(Float64Error const& x) { return _sqr<Precision64>(x); }
Float64LowerBound rec(Float64Error const& x) { return _rec<Precision64>(x); }
Float64Error pow(Float64Error const& x1, Nat m2) { return _pow<Precision64>(x1,m2); }
Float64UpperBound log(Float64Error const& x) { return _log<Precision64>(x); }
Float64Error max(Float64Error const& x1, Float64Error const& x2) { return _max<Precision64>(x1,x2); }
Float64Error abs(Float64Error const& x) { return x; }
Float64Error mag(Float64Error const& x) { return x; }

OutputStream& operator<<(OutputStream& os, Float64Error const& x) { return _write<Precision64>(os,x); }
InputStream& operator>>(InputStream& is, Float64Error& x) { return _read<Precision64>(is,x); }

Bool same(Float64Error const& x1, Float64Error const& x2) { return _same<Precision64>(x1,x2); }
Bool refines(Float64Error const& x1, Float64Error const& x2) { return _refines<Precision64>(x1,x2); }
Float64Error refinement(Float64Error const& x1, Float64Error const& x2) { return _refinement<Precision64>(x1,x2); }


FloatMPUpperBound operator+(FloatMPError const& x) { return x; }
FloatMPLowerBound operator-(FloatMPError const& x) { return -(+x); }
FloatMPError operator+(FloatMPError const& x1, FloatMPError const& x2) { return _add<PrecisionMP>(x1,x2); }
FloatMPError operator*(FloatMPError const& x1, FloatMPError const& x2) { return _mul<PrecisionMP>(x1,x2); }
FloatMPError operator*(Nat y1, FloatMPError const& x2) { return FloatMPError(y1)*x2; }
FloatMPError& operator+=(FloatMPError& x1, FloatMPError const& x2) { return x1=x1+x2; }
FloatMPError& operator*=(FloatMPError& x1, FloatMPError const& x2) { return x1=x1*x2; }
FloatMPError sqr(FloatMPError const& x) { return _sqr<PrecisionMP>(x); }
FloatMPLowerBound rec(FloatMPError const& x) { return _rec<PrecisionMP>(x); }
FloatMPError pow(FloatMPError const& x1, Nat m2) { return _pow<PrecisionMP>(x1,m2); }
FloatMPUpperBound log(FloatMPError const& x) { return _log<PrecisionMP>(x); }
FloatMPError max(FloatMPError const& x1, FloatMPError const& x2) { return _max<PrecisionMP>(x1,x2); }
FloatMPError abs(FloatMPError const& x) { return x; }
FloatMPError mag(FloatMPError const& x) { return x; }

OutputStream& operator<<(OutputStream& os, FloatMPError const& x) { return _write<PrecisionMP>(os,x); }
InputStream& operator>>(InputStream& is, FloatMPError& x) { return _read<PrecisionMP>(is,x); }

Bool same(FloatMPError const& x1, FloatMPError const& x2) { return _same<PrecisionMP>(x1,x2); }
Bool refines(FloatMPError const& x1, FloatMPError const& x2) { return _refines<PrecisionMP>(x1,x2); }
FloatMPError refinement(FloatMPError const& x1, FloatMPError const& x2) { return _refinement<PrecisionMP>(x1,x2); }





PositiveFloat64Approximation max(PositiveFloat64Approximation const& x1, PositiveFloat64Approximation const& x2) {
    return PositiveFloat64Approximation(max(x1._a,x2._a)); }
PositiveFloat64LowerBound max(PositiveFloat64LowerBound const& x1, PositiveFloat64LowerBound const& x2) {
    return PositiveFloat64LowerBound(max(x1._l,x2._l)); }
PositiveFloat64UpperBound max(PositiveFloat64UpperBound const& x1, PositiveFloat64UpperBound const& x2) {
    return PositiveFloat64UpperBound(max(x1._u,x2._u)); }
PositiveFloat64Bounds max(PositiveFloat64Bounds const& x1, PositiveFloat64Bounds const& x2) {
    return PositiveFloat64Bounds(max(x1._l,x2._l),max(x1._u,x2._u)); }
PositiveFloat64Approximation operator*(PositiveFloat64Approximation const& x1, PositiveFloat64Approximation const& x2) {
    return PositiveFloat64Approximation(mul_near(x1._a,x2._a)); }
PositiveFloat64UpperBound& operator*=(PositiveFloat64UpperBound& x1, PositiveFloat64UpperBound const& x2) {
    x1._u=mul_up(x1._u,x2._u); return x1; }
Float64Error operator/(Float64Error const& x1, PositiveFloat64LowerBound const& x2) {
    return Float64Error(div_up(x1._e,x2._l)); }

Float64UpperBound operator*(Float64UpperBound const& x1, Real const& y2) {
    Float64UpperBound x2(y2,x1.precision()); return Float64UpperBound(mul_up(x1._u,x2._u)); }
Float64UpperBound rec(Float64LowerBound const& x) {
    return Float64UpperBound(rec_up(x._l)); }


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
/*
    return instantiate_float<ApproximateTag,PR>()
        + instantiate_float<LowerTag,PR>()
        + instantiate_float<UpperTag,PR>()
        + instantiate_float<BoundedTag,PR>()
        + instantiate_float<MetricTag,PR>()
        + instantiate_float<ExactTag,PR>()
        + instantiate_float<PositiveUpperTag,PR>()
        ;
*/
}

template std::size_t instantiate_floats<Precision64>();
template std::size_t instantiate_floats<PrecisionMP>();


Float64Value midpoint(Float64Bounds const& x) { return x.value(); }


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
