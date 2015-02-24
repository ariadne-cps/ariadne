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


template<class PR> Nat Float<Approximate,PR>::output_precision = 4;
template<class PR> Nat Float<Bounded,PR>::output_precision=8;
template<class PR> Nat Float<Exact,PR>::output_precision = 16;

const Float<Exact,Precision64> infty = Float<Exact,Precision64>(Float64::inf());

Float<Error,Precision64> operator"" _error(long double lx) {
    double x=lx;
    assert(x==lx);
    return Float<Error,Precision64>(Float64(x));
}

Float<Exact,Precision64> operator"" _exact(long double lx) {
    double x=lx;
    assert(x==lx);
    return Float<Exact,Precision64>(x);
}

Float<Metric,Precision64> operator"" _near(long double lx) {
    volatile double x=lx;
    volatile long double le=std::abs((long double)x-lx);
    volatile double e=le;
    while(e<le) { e*=(1+std::numeric_limits<double>::epsilon()); }
    return Float<Metric,Precision64>(x,e);
}

Float<Upper,Precision64> operator"" _upper(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x<lx) { x+=min; }
    while (x<lx) { x+=std::abs(x)*eps; }
    return Float<Upper,Precision64>(x);
}

Float<Lower,Precision64> operator"" _lower(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x>lx) { x-=min; }
    while (x>lx) { x-=std::abs(x)*eps; }
    return Float<Lower,Precision64>(x);
}

Float<Approximate,Precision64> operator"" _approx(long double lx) {
    double x=lx;
    return Float<Approximate,Precision64>(x);
}


TwoExp::operator Float<Exact,Precision64> () const {
    return Float<Exact,Precision64>(this->get_d());
}


/*
template<class PR> Float<Exact,PR>::Float(Rational const& q, PR pr)
    : _v(0.0,pr)
{
    Float<Bounded,PR> x(q,pr);
    if(x.lower_raw()==x.upper_raw()) { this->_v==x.value_raw(); }
    else { ARIADNE_THROW(std::runtime_error,"ExactFloat(Rational q, Precision pr)","q="<<q<<" cannot be expressed exactly to precision \n"); }
}
*/

template<class PR> Float<Exact,PR>::operator Number<Exact>() const {
    return Number<Exact>(new NumberWrapper<Float<Exact,PR>>(*this));
}

template<class PR> Float<Metric,PR>::Float(Rational const& q, PR pr)
    : _v(RawFloat<PR>(q,RawFloat<PR>::to_nearest,pr)), _e(abs(Rational(_v)-q),RawFloat<PR>::upward,pr) {
}

template<class PR> Float<Metric,PR>::Float(Number<Validated> const& x)
    : Float(x.get(Metric(),RawFloat<PR>::get_default_precision())) {
}

template<class PR> Float<Metric,PR>::Float(Number<Validated> const& x, PR pr)
    : Float(x.get(Metric(),pr)) {
}

template<class PR> Float<Metric,PR>::operator Number<Validated>() const {
    return Number<Validated>(new NumberWrapper<Float<Metric,PR>>(*this));
}

template<class PR> Float<Bounded,PR>::Float(Real const& x)
    : Float(x(FLT::get_default_precision())) {
}

template<class PR> Float<Bounded,PR>::Float(Rational const& q, PR pr)
    : Float(RawFloat<PR>(q,RawFloat<PR>::downward,pr),RawFloat<PR>(q,RawFloat<PR>::upward,pr)) {
}

template<class PR> Float<Bounded,PR>::Float(Real const& x, PR pr)
    : Float(x(pr)) {
}

template<class PR> Float<Bounded,PR>::Float(Number<Validated> const& x)
    : Float(x.get(Bounded(),FLT::get_default_precision())) {
}

template<class PR> Float<Bounded,PR>::Float(Number<Validated> const& x, PR pr)
    : Float(x.get(Bounded(),pr)) {
}

template<class PR> Float<Bounded,PR>::operator Number<Validated>() const {
    return Number<Validated>(new NumberWrapper<Float<Bounded,PR>>(*this));
}

template<class PR> Float<Upper,PR>::Float(Rational const& q)
    : Float(Float<Bounded,PR>(q)) {
}

template<class PR> Float<Upper,PR>::Float(Rational const& q, PR pr)
    : Float(Float<Bounded,PR>(q,pr)) {
}

template<class PR> Float<Upper,PR>::Float(Number<Upper> const& x, PR pr)
    : Float(x.get(Upper(),pr)) {
}

template<class PR> Float<Upper,PR>::operator Number<Upper>() const {
    return Number<Upper>(new NumberWrapper<Float<Upper,PR>>(*this));
}

template<class PR> Float<Lower,PR>::Float(Rational const& q)
    : Float(Float<Bounded,PR>(q)) {
}

template<class PR> Float<Lower,PR>::Float(Rational const& q, PR pr)
    : Float(Float<Bounded,PR>(q,pr)) {
}

template<class PR> Float<Lower,PR>::Float(Number<Lower> const& x, PR pr)
    : Float(x.get(Lower(),pr)) {
}

template<class PR> Float<Lower,PR>::operator Number<Lower>() const {
    return Number<Lower>(new NumberWrapper<Float<Lower,PR>>(*this));
}

template<class PR> Float<Approximate,PR>::Float(Rational const& q, PR pr)
    : Float(Float<Bounded,PR>(q,pr)) {
}

template<class PR> Float<Approximate,PR>::Float(Number<Approximate> const& x, PR pr)
    : Float(x.get(Approximate(),pr)) {
}

template<class PR> Float<Approximate,PR>::operator Number<Approximate>() const {
    return Number<Approximate>(new NumberWrapper<Float<Approximate,PR>>(*this));
}


template<class PR> Float<Exact,PR>::Float(Integer const& z) : Float(z,RawFloatType::get_default_precision()) { }

template<class PR> Float<Exact,PR>::Float(Integer const& z, PR pr)
    : _v(Rational(z),RawFloat<PR>::to_nearest,pr)
{
    Rational q(_v);
    ARIADNE_PRECONDITION(z==q);
}

template<class PR> Float<Bounded,PR>::Float(const Dyadic& b) : Float<Bounded,PR>(Rational(b)) { }

template<class PR> Float<Bounded,PR>::Float(const Decimal& d) : Float<Bounded,PR>(Rational(d)) { }

template<class PR> Float<Bounded,PR>::Float(const Integer& z) : Float<Bounded,PR>(Rational(z)) { }

template<class PR> Float<Bounded,PR>::Float(const Rational& q) : _l(q.get_d()), _u(_l)  {
    while(Rational(_l)>q) { _l=next_down(_l); }
    while(Rational(_u)<q) { _u=next_up(_u); }
}


template<class PR> Float<Metric,PR> Float<Exact,PR>::pm(Float<Error,PR> _e) const {
    Float<Exact,PR> const& _v=*this; return Float<Metric,PR>(_v,_e);
}


template<class PR> Float<Upper,PR>::Float(Number<Upper> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class PR> Float<Lower,PR>::Float(Number<Lower> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class PR> Float<Approximate,PR>::Float(Dyadic const& b) : Float<Approximate,PR>(b.operator Rational()) { }
template<class PR> Float<Approximate,PR>::Float(Decimal const& d) : Float<Approximate,PR>(d.operator Rational()) { }
template<class PR> Float<Approximate,PR>::Float(Rational const& q) : Float<Approximate,PR>(RawFloatType(q,RawFloatType::to_nearest)) { }

template<class PR> Float<Approximate,PR>::Float(Number<Approximate> const& x) { ARIADNE_NOT_IMPLEMENTED; }

template<class PR> Float<Exact,PR>::operator Rational() const {
    return Rational(this->get_d());
}


template<class PR> Float<Approximate,PR>::Float(Float<Lower,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Upper,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Bounded,PR> const& x) : _a(x.value_raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Metric,PR> const& x) : _a(x.value_raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Exact,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<Lower,PR>::Float(Float<Bounded,PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> Float<Lower,PR>::Float(Float<Metric,PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> Float<Lower,PR>::Float(Float<Exact,PR> const& x) : _l(x.raw()) {
}

template<class PR> Float<Upper,PR>::Float(Float<Bounded,PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> Float<Upper,PR>::Float(Float<Metric,PR> const& x) : _u(x.upper_raw()) {
}

template<class PR> Float<Upper,PR>::Float(Float<Exact,PR> const& x) : _u(x.raw()) {
}

template<class PR> Float<Bounded,PR>::Float(Float<Metric,PR> const& x) : _l(x.lower_raw()), _u(x.upper_raw()) {
}

template<class PR> Float<Bounded,PR>::Float(Float<Exact,PR> const& x) : _l(x.raw()), _u(x.raw()) {
}

template<class PR> Float<Metric,PR>::Float(Float<Bounded,PR> const& x) : _v(x.value_raw()), _e(x.error_raw()) {
}

template<class PR> Float<Metric,PR>::Float(Float<Exact,PR> const& x) : _v(x.raw()), _e(0.0) {
}





template<class PR> Float<Approximate,PR> floor(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(floor(x._a)); }
template<class PR> Float<Approximate,PR> ceil(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(ceil(x._a)); }
template<class PR> Float<Approximate,PR> round(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(round(x._a)); }

template<class PR> Float<Approximate,PR> abs(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(abs_exact(x._a)); }
template<class PR> Float<Approximate,PR> max(Float<Approximate,PR> const& x, Float<Approximate,PR> y) {
    return Float<Approximate,PR>(max_exact(x._a,y._a)); }
template<class PR> Float<Approximate,PR> min(Float<Approximate,PR> const& x, Float<Approximate,PR> y) {
    return Float<Approximate,PR>(min_exact(x._a,y._a)); }
template<class PR> Float<PositiveApproximate,PR> mag(Float<Approximate,PR> const& x) {
    return Float<PositiveApproximate,PR>(abs(x._a)); }
template<class PR> Float<PositiveApproximate,PR> mig(Float<Approximate,PR> const& x) {
    return Float<PositiveApproximate,PR>(abs(x._a)); }

template<class PR> Float<Approximate,PR> nul(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(nul_exact(x._a)); }
template<class PR> Float<Approximate,PR> pos(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(pos_exact(x._a)); }
template<class PR> Float<Approximate,PR> neg(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(neg_exact(x._a)); }
template<class PR> Float<Approximate,PR> half(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(half_exact(x._a)); }
template<class PR> Float<Approximate,PR> sqr(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(mul_near(x._a,x._a)); }
template<class PR> Float<Approximate,PR> rec(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(div_near(1.0,x._a)); }

template<class PR> Float<Approximate,PR> add(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return Float<Approximate,PR>(add_near(x1._a,x2._a)); }
template<class PR> Float<Approximate,PR> sub(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return Float<Approximate,PR>(sub_near(x1._a,x2._a)); }
template<class PR> Float<Approximate,PR> mul(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return Float<Approximate,PR>(mul_near(x1._a,x2._a)); }
template<class PR> Float<Approximate,PR> div(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return Float<Approximate,PR>(div_near(x1._a,x2._a)); }

template<class PR> Float<Approximate,PR> pow(Float<Approximate,PR> const& x, Nat m) {
    return Float<Approximate,PR>(pow_approx(x._a,m)); }
template<class PR> Float<Approximate,PR> pow(Float<Approximate,PR> const& x, Int n) {
    return Float<Approximate,PR>(pow_approx(x._a,n)); }

template<class PR> Float<Approximate,PR> sqrt(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(sqrt_approx(x._a)); }
template<class PR> Float<Approximate,PR> exp(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(exp_approx(x._a)); }
template<class PR> Float<Approximate,PR> log(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(log_approx(x._a)); }
template<class PR> Float<Approximate,PR> sin(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(sin_approx(x._a)); }
template<class PR> Float<Approximate,PR> cos(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(cos_approx(x._a)); }
template<class PR> Float<Approximate,PR> tan(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(tan_approx(x._a)); }
template<class PR> Float<Approximate,PR> asin(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(asin_approx(x._a)); }
template<class PR> Float<Approximate,PR> acos(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(acos_approx(x._a)); }
template<class PR> Float<Approximate,PR> atan(Float<Approximate,PR> const& x) {
    return Float<Approximate,PR>(atan_approx(x._a)); }

template<class PR> Logical<Approximate> eq(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return x1._a==x2._a; }
template<class PR> Logical<Approximate> leq(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return x1._a<=x2._a; }

template<class PR> Bool same(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return x1._a==x2._a; }

template<class PR> Float<Approximate,PR> operator+(Float<Approximate,PR> const& x) { return pos(x); }
template<class PR> Float<Approximate,PR> operator-(Float<Approximate,PR> const& x) { return neg(x); }
template<class PR> Float<Approximate,PR> operator+(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<Approximate,PR> operator-(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<Approximate,PR> operator*(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<Approximate,PR> operator/(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<Approximate,PR>& operator+=(Float<Approximate,PR>& x1, Float<Approximate,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<Approximate,PR>& operator-=(Float<Approximate,PR>& x1, Float<Approximate,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<Approximate,PR>& operator*=(Float<Approximate,PR>& x1, Float<Approximate,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<Approximate,PR>& operator/=(Float<Approximate,PR>& x1, Float<Approximate,PR> const& x2) { return x1=x1/x2; }

template<class PR> Fuzzy operator==(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return x1._a==x2._a; }
template<class PR> Fuzzy operator!=(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return x1._a!=x2._a; }
template<class PR> Fuzzy operator<=(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return x1._a<=x2._a; }
template<class PR> Fuzzy operator>=(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return x1._a>=x2._a; }
template<class PR> Fuzzy operator< (Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return x1._a< x2._a; }
template<class PR> Fuzzy operator> (Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) { return x1._a> x2._a; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Approximate,PR> const& x) {
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_to_nearest();
    os << std::showpoint << std::setprecision(Float<Approximate,PR>::output_precision) << x.raw();
    RawFloat<PR>::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<Approximate,PR>& x) {
    is >> x._a;
    return is;
}



template<class PR> Float<Lower,PR> max(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Lower,PR>(max_exact(x1._l,x2._l)); }
template<class PR> Float<Lower,PR> min(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Lower,PR>(min_exact(x1._l,x2._l)); }
template<class PR> Float<Approximate,PR> abs(Float<Lower,PR> const& x) {
    return abs(Float<Approximate,PR>(x)); }

template<class PR> Float<Lower,PR> nul(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(pos_exact(x._l)); }
template<class PR> Float<Lower,PR> pos(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(pos_exact(x._l)); }
template<class PR> Float<Lower,PR> neg(Float<Upper,PR> const& x) {
    return Float<Lower,PR>(neg_exact(x._u)); }
template<class PR> Float<Lower,PR> half(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(half_exact(x._l)); }

template<class PR> Float<Lower,PR> rec(Float<Upper,PR> const& x) {
    return Float<Lower,PR>(rec_down(x.raw())); }

template<class PR> Float<Lower,PR> add(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Lower,PR>(add_down(x1._l,x2._l)); }

template<class PR> Float<Approximate,PR> sub(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Upper,PR>(sub_near(x1._l,x2._l)); }

template<class PR> Float<Lower,PR> sub(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Lower,PR>(sub_down(x1._l,x2._u)); }

template<class PR> Float<Lower,PR> mul(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    return Float<Lower,PR>(mul_down(x1.raw(),x2.raw())); }

template<class PR> Float<Lower,PR> div(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Lower,PR>(div_down(x1.raw(),x2.raw())); }

template<class PR> Float<Lower,PR> pow(Float<Lower,PR> const& x, Nat m) {
    ARIADNE_PRECONDITION(x.raw()>=0);
    return Float<Lower,PR>(pow_down(x.raw(),m)); }

template<class PR> Float<Approximate,PR> pow(Float<Lower,PR> const& x, Int n) {
    return pow(Float<Approximate,PR>(x),n); }

template<class PR> Float<Lower,PR> sqrt(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(sqrt_down(x.raw())); }

template<class PR> Float<Lower,PR> exp(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(exp_down(x.raw())); }

template<class PR> Float<Lower,PR> log(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(log_down(x.raw())); }

template<class PR> Float<Lower,PR> atan(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(atan_down(x.raw())); }

template<class PR> Logical<Lower> eq(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) {
    if(x1._l>x2._u) { return false; }
    else { return Logical<Lower>(LogicalValue::INDETERMINATE); }
}

template<class PR> Logical<Lower> leq(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) {
    if(x1._l>x2._u) { return false; }
    else { return Logical<Lower>(LogicalValue::LIKELY); }
}

template<class PR> Bool same(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return x1._l==x2._l; }

template<class PR> Float<Lower,PR> operator+(Float<Lower,PR> const& x) { return pos(x); }
template<class PR> Float<Lower,PR> operator-(Float<Upper,PR> const& x) { return neg(x); }
template<class PR> Float<Lower,PR> operator+(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<Lower,PR> operator-(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<Lower,PR> operator*(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<Lower,PR> operator/(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<Lower,PR>& operator+=(Float<Lower,PR>& x1, Float<Lower,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<Lower,PR>& operator-=(Float<Lower,PR>& x1, Float<Upper,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<Lower,PR>& operator*=(Float<Lower,PR>& x1, Float<Lower,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<Lower,PR>& operator/=(Float<Lower,PR>& x1, Float<Upper,PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Lower,PR> const& x) {
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_downward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.raw();
    RawFloat<PR>::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<Lower,PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
}



template<class PR> Float<Upper,PR> max(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(max_exact(x1._u,x2._u)); }

template<class PR> Float<Upper,PR> min(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(min_exact(x1._u,x2._u)); }

template<class PR> Float<Approximate,PR> abs(Float<Upper,PR> const& x) {
    return abs(Float<Approximate,PR>(x)); }

template<class PR> Float<Upper,PR> nul(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(pos_exact(x._u)); }

template<class PR> Float<Upper,PR> pos(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(pos_exact(x._u)); }

template<class PR> Float<Upper,PR> neg(Float<Lower,PR> const& x) {
    return Float<Upper,PR>(neg_exact(x._l)); }

template<class PR> Float<Upper,PR> half(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(half_exact(x._u)); }

template<class PR> Float<Upper,PR> sqr(Float<Upper,PR> const& x) {
    ARIADNE_ASSERT(false); return Float<Upper,PR>(mul_up(x._u,x._u)); }

template<class PR> Float<Upper,PR> rec(Float<Lower,PR> const& x) {
    return Float<Upper,PR>(rec_up(x.raw())); }

template<class PR> Float<Upper,PR> add(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(add_up(x1._u,x2._u)); }

template<class PR> Float<Approximate,PR> sub(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(sub_near(x1._u,x2._u)); }

template<class PR> Float<Upper,PR> sub(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Upper,PR>(sub_up(x1._u,x2._l)); }

template<class PR> Float<Upper,PR> mul(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
//    ARIADNE_WARN("Multiplying UpperFloat "<<x1<<" with UpperFloat "<<x2<<" is unsafe");
    ARIADNE_PRECONDITION(x1.raw()>=0);
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return Float<Upper,PR>(mul_up(x1._u,x2._u)); }

template<class PR> Float<Upper,PR> div(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) {
//    ARIADNE_WARN("Dividing UpperFloat "<<x1<<" by LowerFloat "<<x2<<" is unsafe");
    ARIADNE_PRECONDITION(x1.raw()>=0);
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return Float<Upper,PR>(div_up(x1._u,x2._l)); }

template<class PR> Float<Upper,PR> pow(Float<Upper,PR> const& x, Nat m) {
    ARIADNE_PRECONDITION(x.raw()>=0);
    return Float<Upper,PR>(pow_up(x._u,m)); }

template<class PR> Float<Approximate,PR> pow(Float<Upper,PR> const& x, Int n) {
    return pow(Float<Approximate,PR>(x),n); }

template<class PR> Float<Upper,PR> sqrt(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(sqrt_up(x.raw())); }

template<class PR> Float<Upper,PR> exp(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(exp_up(x.raw())); }

template<class PR> Float<Upper,PR> log(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(log_up(x.raw())); }

template<class PR> Float<Upper,PR> atan(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(atan_up(x.raw())); }

template<class PR> Logical<Lower> eq(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) {
    if(x1._u<x2._l) { return false; }
    else { return Logical<Lower>(LogicalValue::INDETERMINATE); }
}

template<class PR> Logical<Upper> leq(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) {
    if(x1._u<=x2._l) { return true; }
    else { return Logical<Upper>(LogicalValue::UNLIKELY); }
}

template<class PR> Bool same(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return x1._u==x2._u; }

template<class PR> Integer integer_cast(Float<Upper,PR> const& x) { return Integer(static_cast<int>(x._u.get_d())); }

template<class PR> Float<Upper,PR> operator+(Float<Upper,PR> const& x) { return pos(x); }
template<class PR> Float<Upper,PR> operator-(Float<Lower,PR> const& x) { return neg(x); }
template<class PR> Float<Upper,PR> operator+(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<Upper,PR> operator-(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<Upper,PR> operator*(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<Upper,PR> operator/(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<Upper,PR>& operator+=(Float<Upper,PR>& x1, Float<Upper,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<Upper,PR>& operator-=(Float<Upper,PR>& x1, Float<Lower,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<Upper,PR>& operator*=(Float<Upper,PR>& x1, Float<Upper,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<Upper,PR>& operator/=(Float<Upper,PR>& x1, Float<Lower,PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Upper,PR> const& x) {
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    RawFloat<PR>::set_rounding_upward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.raw();
    RawFloat<PR>::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<Upper,PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
}








template<class PR> Float<Bounded,PR> round(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(round(x.lower_raw()),round(x.upper_raw()));
}

template<class PR> Float<Bounded,PR> max(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
}

template<class PR> Float<Bounded,PR> min(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
}


template<class PR> Float<Bounded,PR> abs(Float<Bounded,PR> const& x) {
    if(x.lower_raw()>=0) {
        return Float<Bounded,PR>(x.lower_raw(),x.upper_raw());
    } else if(x.upper_raw()<=0) {
        return Float<Bounded,PR>(neg(x.upper_raw()),neg(x.lower_raw()));
    } else {
        return Float<Bounded,PR>(static_cast<RawFloat<PR>>(0.0),max(neg(x.lower_raw()),x.upper_raw()));
    }
}

template<class PR> Float<PositiveLower,PR> mig(Float<Bounded,PR> const& x) {
    return Float<PositiveLower,PR>(max(0,max(x._l,neg(x._u))));
}

template<class PR> Float<PositiveUpper,PR> mag(Float<Bounded,PR> const& x) {
    return Float<PositiveUpper,PR>(max(neg(x._l),x._u));
}

template<class PR> Float<Bounded,PR> nul(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(nul(x._l),nul(x._u));
}

template<class PR> Float<Bounded,PR> pos(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(pos(x._l),pos(x._u));
}

template<class PR> Float<Bounded,PR> neg(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(neg(x._u),neg(x._l));
}

template<class PR> Float<Bounded,PR> half(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(half(x._l),half(x._u));
}

template<class PR> Float<Bounded,PR> sqr(Float<Bounded,PR> const& x) {
    const RawFloat<PR>& xl=x.lower_raw(); const RawFloat<PR>& xu=x.upper_raw();
    RawFloat<PR> rl,ru;
    if(xl>0.0) {
        rl=mul_down(xl,xl); ru=mul_up(xu,xu);
    } else if(xu<0.0) {
        rl=mul_down(xu,xu); ru=mul_up(xl,xl);
    } else {
        rl=nul(xl); ru=max(mul_up(xl,xl),mul_up(xu,xu));
    }
    return Float<Bounded,PR>(rl,ru);
}

template<class PR> Float<Bounded,PR> rec(Float<Bounded,PR> const& x) {
    // IMPORTANT: Need to be careful when one of the bounds is 0, since if xl=-0.0 and xu>0, then 1/xl=-inf
    if(x._l>0 || x._u<0) {
        return Float<Bounded,PR>(rec_down(x._u),rec_up(x._l));
    } else {
        RawFloat<PR> inf=RawFloat<PR>::inf(x.precision());
        RawFloat<PR> rl=-inf; RawFloat<PR> ru=+inf;
        //ARIADNE_THROW(DivideByZeroException,"BoundedFloat rec(BoundedFloat x)","x="<<x);
        return Float<Bounded,PR>(-inf,+inf);
    }
}

template<class PR> Float<Bounded,PR> add(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(add_down(x1._l,x2._l),add_up(x1._u,x2._u));
}

template<class PR> Float<Bounded,PR> sub(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(sub_down(x1._l,x2._u),sub_up(x1._u,x2._l));
}

template<class PR> Float<Bounded,PR> mul(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
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
    return Float<Bounded,PR>(rl,ru);
}

template<class PR> Float<Bounded,PR> div(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
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
        //ARIADNE_THROW(DivideByZeroException,"BoundedFloat div(BoundedFloat x1, BoundedFloat x2)","x1="<<x1<<", x2="<<x2);
        rl=-RawFloat<PR>::inf();
        ru=+RawFloat<PR>::inf();
    }
    return Float<Bounded,PR>(rl,ru);
}






template<class PR> Float<Bounded,PR> pow(Float<Bounded,PR> const& x, Int n) {
    if(n<0) { return pow(rec(x),Nat(-n)); }
    else return pow(x,Nat(n));
}

template<class PR> Float<Bounded,PR> pow(Float<Bounded,PR> const& x, Nat m) {
    Float<Bounded,PR> y = x;
    if(m%2==0) { y=abs(x); }
    RawFloat<PR> rl=pow_down(y.lower_raw(),m);
    RawFloat<PR> ru=pow_up(y.upper_raw(),m);
    return Float<Bounded,PR>(rl,ru);
}


template<class PR> Float<Bounded,PR> sqrt(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(sqrt_down(x.lower_raw()),sqrt_up(x.upper_raw()));
}

template<class PR> Float<Bounded,PR> exp(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(exp_down(x.lower_raw()),exp_up(x.upper_raw()));
}

template<class PR> Float<Bounded,PR> log(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(log_down(x.lower_raw()),log_up(x.upper_raw()));
}


template<class PR> Float<Bounded,PR> pi_val(PR pr) { return Float<Bounded,PR>(pi_down(pr),pi_up(pr)); }

template<class PR> Float<Bounded,PR> sin(Float<Bounded,PR> const& x)
{
    return cos(x-half(pi_val<PR>(x.precision())));
}

template<class PR> Float<Bounded,PR> cos(Float<Bounded,PR> const& x)
{
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    typename RawFloat<PR>::RoundingModeType rnd = RawFloat<PR>::get_rounding_mode();
    PR prec=x.precision();

    static const Float<Exact,PR> two(2);

    if(x.error().raw()>2*pi_down(prec)) { return Float<Bounded,PR>(-1.0,+1.0); }

    auto n=floor(x.lower_raw()/(2*pi_approx(prec))+0.5);
    Float<Bounded,PR> y=x-two*Float<Exact,PR>(n)*pi_val<PR>(x.precision());

    ARIADNE_ASSERT(y.lower_raw()<=pi_up(prec));
    ARIADNE_ASSERT(y.upper_raw()>=-pi_up(prec));

    RawFloat<PR> rl,ru;
    if(y.lower_raw()<=-pi_down(prec)) {
        if(y.upper_raw()<=0.0) { rl=-1.0; ru=cos_up(y.upper_raw()); }
        else { rl=-1.0; ru=+1.0; }
    } else if(y.lower_raw()<=0.0) {
        if(y.upper_raw()<=0.0) { rl=cos_down(y.lower_raw()); ru=cos_up(y.upper_raw()); }
        else if(y.upper_raw()<=pi_down(prec)) { rl=cos_down(max(-y.lower_raw(),y.upper_raw())); ru=+1.0; }
        else { rl=-1.0; ru=+1.0; }
    } else if(y.lower_raw()<=pi_up(prec)) {
        if(y.upper_raw()<=pi_down(prec)) { rl=cos_down(y.upper_raw()); ru=cos_up(y.lower_raw()); }
        else if(y.upper_raw()<=2*pi_down(prec)) { rl=-1.0; ru=cos_up(min(y.lower_raw(),sub_down(2*pi_down(prec),y.upper_raw()))); }
        else { rl=-1.0; ru=+1.0; }
    } else {
        assert(false);
    }

    RawFloat<PR>::set_rounding_mode(rnd);
    return Float<Bounded,PR>(rl,ru);
}

template<class PR> Float<Bounded,PR> tan(Float<Bounded,PR> const& x) {
    return mul(sin(x),rec(cos(x)));
}

template<class PR> Float<Bounded,PR> asin(Float<Bounded,PR> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class PR> Float<Bounded,PR> acos(Float<Bounded,PR> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class PR> Float<Bounded,PR> atan(Float<Bounded,PR> const& x) {
    return Float<Bounded,PR>(atan_down(x._l),atan_up(x._u));
}



//! \related Float<Bounded,PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> Logical<Validated> eq(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    if(x1.upper_raw()<x2.lower_raw() || x1.lower_raw()>x2.upper_raw()) { return false; }
    else if(x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw()) { return true; }
    else { return indeterminate; }
}

//! \related Float<Bounded,PR> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
template<class PR> Logical<Validated> leq(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    if(x1.upper_raw()<=x2.lower_raw()) { return true; }
    else if(x1.lower_raw()> x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

template<class PR> Logical<Validated> operator==(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return eq(x1,x2);
}

template<class PR> Logical<Validated> operator!=(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return !eq(x1,x2); }

template<class PR> Logical<Validated> operator<=(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return leq(x1,x2);
}

template<class PR> Logical<Validated> operator>=(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return leq(x2,x1);
}

template<class PR> Logical<Validated> operator< (Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return !leq(x2,x1);
}

template<class PR> Logical<Validated> operator> (Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return !leq(x1,x2);
}

template<class PR> Float<Bounded,PR> widen(Float<Bounded,PR> const& x)
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
    return Float<Bounded,PR>(wl,wu);
}

template<class PR> Float<Bounded,PR> narrow(Float<Bounded,PR> const& x)
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
    return Float<Bounded,PR>(nl,nu);
}

template<class PR> Float<Bounded,PR> trunc(Float<Bounded,PR> const& x)
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
    return Float<Bounded,PR>(double(tl),double(tu));
}

template<class PR> Float<Bounded,PR> trunc(Float<Bounded,PR> const& x, Nat n)
{
    Float<Bounded,PR> _e=Float<Bounded,PR>(std::pow(2.0,52-(Int)n));
    Float<Bounded,PR> y=x+_e;
    return y-_e;
}

template<class PR> Integer integer_cast(Float<Bounded,PR> const& x) {
    return Integer(static_cast<int>(x.value_raw().get_d()));
}

template<class PR> auto is_zero(Float<Bounded,PR> const& x) -> Logical<Validated> {
    if(x.lower_raw()>0.0 || x.upper_raw()<0.0) { return false; }
    else if(x.lower_raw()==0.0 && x.upper_raw()==0.0) { return true; }
    else { return indeterminate; }
}

template<class PR> auto is_positive(Float<Bounded,PR> const& x) -> Logical<Validated> {
    if(x.lower_raw()>=0.0) { return true; }
    else if(x.upper_raw()<0.0) { return false; }
    else { return indeterminate; }
}

template<class PR> auto is_positive(Float<Bounded,PR> const&) -> Logical<Validated>;

template<class PR> Bool same(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1._l==x2._l && x1._u==x2._u; }


template<class PR> OutputStream& operator<<(OutputStream& os, const Float<Bounded,PR>& x)
{
    //if(x.lower_raw()==x.upper_raw()) { return os << "{" << std::setprecision(Float<Bounded,PR>::output_precision) << x.lower_raw().get_d() << ; }
    typename RawFloat<PR>::RoundingModeType rnd=RawFloat<PR>::get_rounding_mode();
    os << '{';
    RawFloat<PR>::set_rounding_downward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.lower().get_d();
    os << ':';
    RawFloat<PR>::set_rounding_upward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.upper().get_d();
    RawFloat<PR>::set_rounding_mode(rnd);
    os << '}';
    return os;

}

template<class PR> InputStream& operator>>(InputStream& is, Float<Bounded,PR>& x)
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



template<class PR> Float<Bounded,PR> operator+(Float<Bounded,PR> const& x) { return pos(x); }
template<class PR> Float<Bounded,PR> operator-(Float<Bounded,PR> const& x) { return neg(x); }
template<class PR> Float<Bounded,PR> operator+(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<Bounded,PR> operator-(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<Bounded,PR> operator*(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<Bounded,PR> operator/(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) { return div(x1,x2); }
template<class PR> Float<Bounded,PR>& operator+=(Float<Bounded,PR>& x1, Float<Bounded,PR> const& x2) { return x1=x1+x2; }
template<class PR> Float<Bounded,PR>& operator-=(Float<Bounded,PR>& x1, Float<Bounded,PR> const& x2) { return x1=x1-x2; }
template<class PR> Float<Bounded,PR>& operator*=(Float<Bounded,PR>& x1, Float<Bounded,PR> const& x2) { return x1=x1*x2; }
template<class PR> Float<Bounded,PR>& operator/=(Float<Bounded,PR>& x1, Float<Bounded,PR> const& x2) { return x1=x1/x2; }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Bounded,PR> const& x);
template<class PR> InputStream& operator>>(InputStream& is, Float<Bounded,PR>& x);




template<class PR> Float<Metric,PR> nul(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(nul(x._v),nul(x._e));
}

template<class PR> Float<Metric,PR> pos(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(pos(x._v),x._e);
}

template<class PR> Float<Metric,PR> neg(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(neg(x._v),x._e);
}

template<class PR> Float<Metric,PR> half(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(half(x._v),half(x._e));
}

template<class PR> Float<Metric,PR> sqr(Float<Metric,PR> const& x) {
    Float<Metric,PR> r=x*x;
    if(r._e>r._v) {
        r._e=half(add_up(r._e,r._v));
        r._v=r._e;
    }
    return r;
}

template<class PR> Float<Metric,PR> rec(Float<Metric,PR> const& x) {
    auto ru=rec_up(add_down(x._v,x._e));
    auto rl=rec_down(add_up(x._v,x._e));
    auto re=half(sub_up(ru,rl));
    auto rv=half(add_near(rl,ru));
    return Float<Metric,PR>(rv,re);
}

template<class PR> Float<Metric,PR> add(Float<Metric,PR> const& x, Float<Metric,PR> y) {
    auto rv=add_near(x._v,y._v);
    auto ru=add_up(x._v,y._v);
    auto rl=add_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return Float<Metric,PR>(rv,re);
}

template<class PR> Float<Metric,PR> sub(Float<Metric,PR> const& x, Float<Metric,PR> y) {
    auto rv=sub_near(x._v,y._v);
    auto ru=sub_up(x._v,y._v);
    auto rl=sub_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return Float<Metric,PR>(rv,re);
}

template<class PR> Float<Metric,PR> mul(Float<Metric,PR> const& x, Float<Metric,PR> y) {
    auto rv=mul_near(x._v,y._v);
    auto ru=mul_up(x._v,y._v);
    auto rl=mul_down(x._v,y._v);
    auto re1=add_up(half(sub_up(ru,rl)),mul_up(x._e,y._e));
    auto re2=add_up(mul_up(abs(x._v),y._e),mul_up(x._e,abs(y._v)));
    auto re=add_up(re1,re2);
    return Float<Metric,PR>(rv,re);
}

template<class PR> Float<Metric,PR> div(Float<Metric,PR> const& x, Float<Metric,PR> y) {
    return x*rec(y);
}

template<class PR> Float<Metric,PR> pow(Float<Metric,PR> const& x, Nat m) {
    return Float<Metric,PR>(pow(Float<Bounded,PR>(x),m));
}

template<class PR> Float<Metric,PR> pow(Float<Metric,PR> const& x, Int n) {
    return Float<Metric,PR>(pow(Float<Bounded,PR>(x),n));
}

template<class PR> Float<Metric,PR> sqrt(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(sqrt(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> exp(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(exp(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> log(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(log(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> sin(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(sin(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> cos(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(cos(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> tan(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(tan(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> asin(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(asin(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> acos(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(acos(Float<Bounded,PR>(x)));
}

template<class PR> Float<Metric,PR> atan(Float<Metric,PR> const& x) {
    return Float<Metric,PR>(atan(Float<Bounded,PR>(x)));
}


template<class PR> Float<Metric,PR> abs(Float<Metric,PR> const& x) {
    if(x._e<abs(x._v)) { return x; }
    else { auto rv=half(abs(x._v)+x._e); return Float<Metric,PR>(rv,rv); }
}

template<class PR> Float<Metric,PR> max(Float<Metric,PR> const& x1, Float<Metric,PR> const& x2) {
    return half((x1+x2)+abs(x1-x2));
}

template<class PR> Float<Metric,PR> min(Float<Metric,PR> const& x1, Float<Metric,PR> const& x2) {
    return half((x1+x2)-abs(x1-x2));
}

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Metric,PR> const& x) {
    return os << x.value() << "\u00b1" << x.error();
}


// Mixed Bounded - Exact operations
template<class PR> Float<Bounded,PR> add(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2) {
    return Float<Bounded,PR>(add_down(x1._l,x2._v),add_up(x1._u,x2._v));
}

template<class PR> Float<Bounded,PR> add(Float<Exact,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(add_down(x1._v,x2._l),add_down(x1._v,x2._u));
}

template<class PR> Float<Bounded,PR> sub(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2) {
    return Float<Bounded,PR>(sub_down(x1._l,x2._v),sub_up(x1._u,x2._v));
}

template<class PR> Float<Bounded,PR> sub(Float<Exact,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(sub_down(x1._v,x2._u),sub_up(x1._v,x2._l));
}

template<class PR> Float<Bounded,PR> mul(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2) {
    const RawFloat<PR>& x1l=x1.lower_raw(); const RawFloat<PR>& x1u=x1.upper_raw();
    const RawFloat<PR>& x2v=x2.raw();
    RawFloat<PR> rl,ru;
    if(x2v>=0.0) {
        rl=mul_down(x1l,x2v); ru=mul_up(x1u,x2v);
    } else {
        rl=mul_down(x1u,x2v); ru=mul_up(x1l,x2v);
    }
    return Float<Bounded,PR>(rl,ru);
}


template<class PR> Float<Bounded,PR> mul(Float<Exact,PR> const& x1, Float<Bounded,PR> const& x2) {
    const RawFloat<PR>& x1v=x1.raw();
    const RawFloat<PR>& x2l=x2.lower_raw(); const RawFloat<PR>& x2u=x2.upper_raw();
    RawFloat<PR> rl,ru;
    if(x1v>=0.0) {
        rl=mul_down(x1v,x2l); ru=mul_up(x1v,x2u);
    } else {
        rl=mul_down(x1v,x2u); ru=mul_up(x1v,x2l);
    }
    return Float<Bounded,PR>(rl,ru);
}

template<class PR> Float<Bounded,PR> div(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2)
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
        //ARIADNE_THROW(DivideByZeroException,"BoundedFloat div(BoundedFloat const& x1, ExactFloat x2)","x1="<<x1<<", x2="<<x2);
        rl=-RawFloat<PR>::inf();
        ru=+RawFloat<PR>::inf();
    }
    return Float<Bounded,PR>(rl,ru);
}


template<class PR> Float<Bounded,PR> div(Float<Exact,PR> const& x1, Float<Bounded,PR> const& x2)
{
    const RawFloat<PR>& x1v=x1.raw();
    const RawFloat<PR>& i2l=x2.lower_raw();
    const RawFloat<PR>& i2u=x2.upper_raw();
    RawFloat<PR> rl,ru;
    if(i2l<=0 && i2u>=0) {
        //ARIADNE_THROW(DivideByZeroException,"BoundedFloat div(ExactFloat const& x1, BoundedFloat x2)","x1="<<x1<<", x2="<<x2);
        rl=-RawFloat<PR>::inf();
        ru=+RawFloat<PR>::inf();
    } else if(x1v>=0) {
        rl=div_down(x1v,i2u); ru=div_up(x1v,i2l);
    } else {
        rl=div_down(x1v,i2l); ru=div_up(x1v,i2u);
    }
    return Float<Bounded,PR>(rl,ru);
}






template<class PR> Float<Exact,PR> max(Float<Exact,PR> const& x1,  Float<Exact,PR> const& x2) {
    return Float<Exact,PR>(max(x1._v,x2._v)); }

template<class PR> Float<Exact,PR> min(Float<Exact,PR> const& x1,  Float<Exact,PR> const& x2) {
    return Float<Exact,PR>(min(x1._v,x2._v)); }

template<class PR> Float<Exact,PR> abs(Float<Exact,PR> const& x) {
    return Float<Exact,PR>(abs(x._v)); }

template<class PR> Float<PositiveExact,PR> mig(Float<Exact,PR> const& x) {
    return Float<PositiveExact,PR>(abs(x._v)); }

template<class PR> Float<PositiveExact,PR> mag(Float<Exact,PR> const& x) {
    return Float<PositiveExact,PR>(abs(x._v)); }


template<class PR> Float<Exact,PR> nul(Float<Exact,PR> const& x) {
    return Float<Exact,PR>(nul(x._v)); }

template<class PR> Float<Exact,PR> pos(Float<Exact,PR> const& x) {
    return Float<Exact,PR>(pos(x._v)); }

template<class PR> Float<Exact,PR> neg(Float<Exact,PR> const& x) {
    return Float<Exact,PR>(neg(x._v)); }

template<class PR> Float<Exact,PR> half(Float<Exact,PR> const& x) {
    return Float<Exact,PR>(half(x._v)); }

template<class PR> Float<Bounded,PR> sqr(Float<Exact,PR> const& x) {
    return Float<Bounded,PR>(mul_down(x._v,x._v),mul_up(x._v,x._v)); }

template<class PR> Float<Bounded,PR> rec(Float<Exact,PR> const& x) {
    return Float<Bounded,PR>(rec_down(x._v),rec_up(x._v)); }

template<class PR> Float<Bounded,PR> add(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return Float<Bounded,PR>(add_down(x1._v,x2._v),add_up(x1._v,x2._v)); }

template<class PR> Float<Bounded,PR> sub(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return Float<Bounded,PR>(sub_down(x1._v,x2._v),sub_up(x1._v,x2._v)); }

template<class PR> Float<Bounded,PR> mul(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return Float<Bounded,PR>(mul_down(x1._v,x2._v),mul_up(x1._v,x2._v)); }

template<class PR> Float<Bounded,PR> div(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return Float<Bounded,PR>(div_down(x1._v,x2._v),div_up(x1._v,x2._v)); }

template<class PR> Float<Bounded,PR> pow(Float<Exact,PR> const& x, Nat m) {
    return pow(Float<Bounded,PR>(x),m); }

template<class PR> Float<Bounded,PR> pow(Float<Exact,PR> const& x, Int n) {
    return pow(Float<Bounded,PR>(x),n); }

template<class PR> Float<Bounded,PR> med(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return add(half(x1),half(x2)); }

template<class PR> Float<Bounded,PR> rad(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return sub(half(x2),half(x1)); }

template<class PR> Float<Bounded,PR> sqrt(Float<Exact,PR> const& x) {
    return sqrt(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> exp(Float<Exact,PR> const& x) {
    return exp(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> log(Float<Exact,PR> const& x) {
    return log(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> sin(Float<Exact,PR> const& x) {
    return sin(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> cos(Float<Exact,PR> const& x) {
    return cos(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> tan(Float<Exact,PR> const& x) {
    return tan(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> asin(Float<Exact,PR> const& x) {
    return asin(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> acos(Float<Exact,PR> const& x) {
    return acos(Float<Bounded,PR>(x)); }

template<class PR> Float<Bounded,PR> atan(Float<Exact,PR> const& x) {
    return atan(Float<Bounded,PR>(x)); }

template<class PR> Logical<Exact> eq(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return x1._v == x2._v; }

template<class PR> Logical<Exact> leq(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return x1._v <= x2._v; }

template<class PR> Bool same(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return x1._v==x2._v; }

template<class PR> Float<Exact,PR> operator+(Float<Exact,PR> const& x) { return pos(x); }
template<class PR> Float<Exact,PR> operator-(Float<Exact,PR> const& x) { return neg(x); }
template<class PR> Float<Bounded,PR> operator+(Float<Exact,PR> const& x1,  Float<Exact,PR> const& x2) { return add(x1,x2); }
template<class PR> Float<Bounded,PR> operator-(Float<Exact,PR> const& x1,  Float<Exact,PR> const& x2) { return sub(x1,x2); }
template<class PR> Float<Bounded,PR> operator*(Float<Exact,PR> const& x1,  Float<Exact,PR> const& x2) { return mul(x1,x2); }
template<class PR> Float<Bounded,PR> operator/(Float<Exact,PR> const& x1,  Float<Exact,PR> const& x2) { return div(x1,x2); }

template<class PR> Float<Exact,PR> operator*(Float<Exact,PR> const& x, TwoExp y) {
    Float<Exact,PR> yv=y; return Float<Exact,PR>(x.raw()*yv.raw()); }
template<class PR> Float<Exact,PR> operator/(Float<Exact,PR> const& x, TwoExp y) {
    Float<Exact,PR> yv=y; return Float<Exact,PR>(x.raw()/yv.raw()); }
template<class PR> Float<Exact,PR>& operator*=(Float<Exact,PR>& x, TwoExp y) {
    Float<Exact,PR> yv=y; return x=Float<Exact,PR>(x.raw()*yv.raw()); }
template<class PR> Float<Exact,PR>& operator/=(Float<Exact,PR>& x, TwoExp y) {
    Float<Exact,PR> yv=y; return x=Float<Exact,PR>(x.raw()/yv.raw()); }

template<class PR> Boolean operator==(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) { return x1.raw()==x2.raw(); }
template<class PR> Boolean operator!=(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) { return x1.raw()!=x2.raw(); }
template<class PR> Boolean operator<=(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) { return x1.raw()<=x2.raw(); }
template<class PR> Boolean operator>=(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) { return x1.raw()>=x2.raw(); }
template<class PR> Boolean operator< (Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) { return x1.raw()< x2.raw(); }
template<class PR> Boolean operator> (Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) { return x1.raw()> x2.raw(); }

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Exact,PR> const& x) {
    return os << std::setprecision(Float<Exact,PR>::output_precision) << x.raw();
}

template<class PR> InputStream& operator>>(InputStream& is, Float<Exact,PR>& x) {
    ARIADNE_NOT_IMPLEMENTED;
    auto v = nul(x._v);
    is >> v;
    ARIADNE_ASSERT(is);
    x._v=v;
    return is;
}

template<class PR> Integer integer_cast(Float<Exact,PR> const& x) {
    Integer z=static_cast<int>(x._v.get_d());
    ARIADNE_ASSERT(z==x);
    return std::move(z);
}


Float<Exact,Precision64> make_exact(Real const& x) {
    return make_exact(Float<Approximate,Precision64>(x));
}

template<class PR> Bool operator==(Float<Exact,PR> const& x, const Rational& q) { return Rational(x)==q; }
template<class PR> Bool operator!=(Float<Exact,PR> const& x, const Rational& q) { return Rational(x)!=q; }
template<class PR> Bool operator<=(Float<Exact,PR> const& x, const Rational& q) { return Rational(x)<=q; }
template<class PR> Bool operator>=(Float<Exact,PR> const& x, const Rational& q) { return Rational(x)>=q; }
template<class PR> Bool operator< (Float<Exact,PR> const& x, const Rational& q) { return Rational(x)< q; }
template<class PR> Bool operator> (Float<Exact,PR> const& x, const Rational& q) { return Rational(x)> q; }

template<class PR> Bool operator==(const Rational& q, Float<Exact,PR> const& x) { return q==Rational(x); }
template<class PR> Bool operator!=(const Rational& q, Float<Exact,PR> const& x) { return q!=Rational(x); }
template<class PR> Bool operator<=(const Rational& q, Float<Exact,PR> const& x) { return q<=Rational(x); }
template<class PR> Bool operator>=(const Rational& q, Float<Exact,PR> const& x) { return q>=Rational(x); }
template<class PR> Bool operator< (const Rational& q, Float<Exact,PR> const& x) { return q< Rational(x); }
template<class PR> Bool operator> (const Rational& q, Float<Exact,PR> const& x) { return q> Rational(x); }



template<class PR> Bool refines(Float<Metric,PR> const& x1, Float<Metric,PR> const& x2) {
    return (x1._v>=x2._v ? sub_up(x1._v,x2._v) : sub_up(x2._v,x1._v)) <= sub_down(x2._e, x1._e);
}

template<class PR> Bool models(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2) {
    return x1._l<=x2._v && x1._u >= x2._v;
}

template<class PR> Bool refines(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1._l>=x2._l && x1._u <= x2._u;
}

template<class PR> Float<Bounded,PR> refinement(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(max(x1._l,x2._l),min(x1._u,x2._u));
}

template<class PR> Bool consistent(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1._l<=x2._u && x1._u >= x2._l;
}

template<class PR> Bool inconsistent(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1._l>x2._u || x1._u < x2._l;
}


template<class PR> Float<Metric,PR> refinement(Float<Metric,PR> const& x1, Float<Metric,PR> const& x2) {
    return Float<Metric,PR>(refinement(Float<Bounded,PR>(x1),Float<Bounded,PR>(x2)));
}

template<class PR> Bool refines(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return x1._l>=x2._l;
}

template<class PR> Bool refines(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return x1._u <= x2._u;
}


#ifdef ARIADNE_ENABLE_SERIALIZATION
template<class PR, class A> Void serialize(A& _a, Float<Bounded,PR>& x, const Nat version) {
    _a & x.lower_raw() & x.upper_raw(); }
#endif




template<class PR> Bool same(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return x1._u == x2._u;
}

template<class PR> Float<PositiveUpper,PR> max(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(max(x1._u,x2._u));
}

template<class PR> Float<PositiveUpper,PR> min(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(min(x1._u,x2._u));
}

template<class PR> Float<PositiveUpper,PR> abs(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveUpper,PR>(x._u);
}

template<class PR> Float<PositiveApproximate,PR> mig(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveApproximate,PR>(x._u);
}

template<class PR> Float<PositiveUpper,PR> mag(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveUpper,PR>(x._u);
}

template<class PR> Float<PositiveUpper,PR> pos(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveUpper,PR>(pos(x._u));
}

template<class PR> Float<Lower,PR> neg(Float<PositiveUpper,PR> const& x) {
    return Float<Lower,PR>(neg(x._u));
}

template<class PR> Float<PositiveUpper,PR> rec(Float<PositiveLower,PR> const& x) {
    return Float<PositiveLower,PR>(rec_up(x._l));
}

template<class PR> Float<PositiveUpper,PR> half(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveUpper,PR>(half(x._u));
}

template<class PR> Float<PositiveUpper,PR> sqr(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveUpper,PR>(mul_up(x._u,x._u));
}

template<class PR> Float<PositiveUpper,PR> add(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(add_up(x1._u,x2._u));
}

template<class PR> Float<PositiveUpper,PR> mul(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(mul_up(x1._u,x2._u));
}

template<class PR> Float<PositiveUpper,PR> div(Float<PositiveUpper,PR> const& x1, Float<PositiveLower,PR> const& x2) {
    return Float<PositiveUpper,PR>(div_up(x1._u,x2._l));
}

template<class PR> Float<PositiveUpper,PR> pow(Float<PositiveUpper,PR> const& x, Nat m) {
    return Float<PositiveUpper,PR>(pow_up(x._u,m));
}

template<class PR> Float<PositiveApproximate,PR> pow(Float<PositiveUpper,PR> const& x, Int n) {
    return Float<PositiveApproximate,PR>(pow_approx(x._u,n));
}

template<class PR> Float<Upper,PR> log(Float<PositiveUpper,PR> const& x) {
    return Float<Upper,PR>(log_up(x._u));
}

template<class PR> OutputStream& operator<<(OutputStream& os, Float<PositiveUpper,PR> const& x) {
    return os << static_cast<Float<Upper,PR>const&>(x);
}


template<class PR> Float<PositiveLower,PR> mig(Float<PositiveLower,PR> const& x) {
    return Float<PositiveLower,PR>(abs(x._l)); }

template<class PR> Float<PositiveLower,PR> rec(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveLower,PR>(rec_down(x._u));
}


template<class PR> Float<PositiveApproximate,PR> mul(Float<PositiveApproximate,PR> const& x1, Float<PositiveApproximate,PR> const& x2) {
    return Float<PositiveApproximate,PR>(mul_near(x1._a,x2._a));
}





template<> Nat integer_cast<Nat,ApproximateFloat64>(ApproximateFloat64 const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,ApproximateFloat64>(ApproximateFloat64 const& x) {
    return std::round(x.get_d()); }

template<> Nat integer_cast<Nat,LowerFloat64>(LowerFloat64 const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,LowerFloat64>(LowerFloat64 const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,BoundedFloat64>(BoundedFloat64 const& x) {
    return std::round((x.lower().get_d()+x.upper().get_d())/2); }



template<class PR, class P> Float<P,PR> max(Float<P,PR> const& x1, Float<P,PR> const& x2) { return max<PR>(x1,x2); }
template<class PR, class P> Float<P,PR> min(Float<P,PR> const& x1, Float<P,PR> const& x2) { return min<PR>(x1,x2); }
template<class PR, class P> Float<Weaker<P,Negated<P>>,PR> abs(Float<P,PR> const& x) { return abs<PR>(x); }
template<class PR, class P> Float<Unsigned<Weaker<P,Upper>>,PR> mag(Float<P,PR> const& x) { return mag<PR>(x); }
template<class PR, class P> Float<Unsigned<Weaker<P,Lower>>,PR> mig(Float<P,PR> const& x) { return mig<PR>(x); }

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

template<class PR> Float<Approximate,PR> make_float(Number<Approximate> x) { return Float<Approximate,PR>(x); }
template<class PR> Float<Lower,PR> make_float(Number<Lower> x) { return Float<Lower,PR>(x); }
template<class PR> Float<Upper,PR> make_float(Number<Upper> x) { return Float<Upper,PR>(x); }
template<class PR> Float<Bounded,PR> make_float(Number<Validated> x) { return Float<Bounded,PR>(x); }
template<class PR> Float<Bounded,PR> make_float(Number<Effective> x) { return Float<Bounded,PR>(x); }
template<class PR> Float<Bounded,PR> make_float(Number<Exact> x) { return Float<Bounded,PR>(x); }
template<class PR> Float<Bounded,PR> make_float(Real r) { return Float<Bounded,PR>(r); }
template<class PR> Float<Bounded,PR> make_float(Rational q) { return Float<Bounded,PR>(q); }
template<class PR> Float<Exact,PR> make_float(Integer z) { return Float<Exact,PR>(z); }

template class Float<Approximate,Precision64>;
template class Float<Lower,Precision64>;
template class Float<Upper,Precision64>;
template class Float<Bounded,Precision64>;
template class Float<Metric,Precision64>;
template class Float<Exact,Precision64>;

template class Float<Approximate,PrecisionMP>;
template class Float<Lower,PrecisionMP>;
template class Float<Upper,PrecisionMP>;
template class Float<Bounded,PrecisionMP>;
template class Float<Metric,PrecisionMP>;
template class Float<Exact,PrecisionMP>;


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
    return instantiate_float<Approximate,PR>()
        + instantiate_float<Lower,PR>()
        + instantiate_float<Upper,PR>()
        + instantiate_float<Bounded,PR>()
        + instantiate_float<Metric,PR>()
        + instantiate_float<Exact,PR>()
        + instantiate_float<PositiveUpper,PR>()
        ;
}

template std::size_t instantiate_floats<Precision64>();
template std::size_t instantiate_floats<PrecisionMP>();

template ExactFloat64 operator* <Precision64,Exact>(ExactFloat64 const&, TwoExp);
template ExactFloat64 operator/ <Precision64,Exact>(ExactFloat64 const&, TwoExp);

template PositiveUpperFloat64 abs(PositiveUpperFloat64 const&);
template PositiveUpperFloatMP abs(PositiveUpperFloatMP const&);

template BoundedFloat64 round<Precision64,Bounded>(BoundedFloat64 const&);
template ApproximateFloat64 round<Precision64,Approximate>(ApproximateFloat64 const&);

template Bool models(BoundedFloat64 const&, ExactFloat64 const&);
template Bool consistent(BoundedFloat64 const&, BoundedFloat64 const&);
template Bool inconsistent(BoundedFloat64 const&, BoundedFloat64 const&);
template Bool refines(BoundedFloat64 const&, BoundedFloat64 const&);
template BoundedFloat64 refinement(BoundedFloat64 const&, BoundedFloat64 const&);
template Bool same(BoundedFloat64 const&, BoundedFloat64 const&);

template Bool refines(MetricFloat64 const&, MetricFloat64 const&);
template MetricFloat64 refinement(MetricFloat64 const&, MetricFloat64 const&);

template Bool refines(LowerFloat64 const&, LowerFloat64 const&);
template Bool refines(UpperFloat64 const&, UpperFloat64 const&);

template Bool refines(MetricFloatMP const&, MetricFloatMP const&);
template Bool refines(BoundedFloatMP const&, BoundedFloatMP const&);
template Bool refines(LowerFloatMP const&, LowerFloatMP const&);
template Bool refines(UpperFloatMP const&, UpperFloatMP const&);

ExactFloat64 midpoint(BoundedFloat64 const& x) { return x.value(); }

template Float<Widen<PositiveApproximate>,Precision64> mul<Precision64,PositiveApproximate>(Float<PositiveApproximate,Precision64> const& x1, Float<PositiveApproximate,Precision64> const& x2);

template<> String class_name<ApproximateFloat64>() { return "ApproximateFloat64"; }
template<> String class_name<LowerFloat64>() { return "LowerFloat64"; }
template<> String class_name<UpperFloat64>() { return "UpperFloat64"; }
template<> String class_name<BoundedFloat64>() { return "BoundedFloat64"; }
template<> String class_name<MetricFloat64>() { return "MetricFloat64"; }
template<> String class_name<ExactFloat64>() { return "ExactFloat64"; }
template<> String class_name<ApproximateFloatMP>() { return "ApproximateFloatMP"; }
template<> String class_name<LowerFloatMP>() { return "LowerFloatMP"; }
template<> String class_name<UpperFloatMP>() { return "UpperFloatMP"; }
template<> String class_name<BoundedFloatMP>() { return "BoundedFloatMP"; }
template<> String class_name<MetricFloatMP>() { return "MetricFloatMP"; }
template<> String class_name<ExactFloatMP>() { return "ExactFloatMP"; }

} // namespace Ariadne
