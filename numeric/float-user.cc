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

namespace Ariadne {

template<> Nat ApproximateFloat64::output_precision = 4;
template<> Nat ValidatedFloat64::output_precision = 8;
template<> Nat ExactFloat64::output_precision = 16;


template<class PR> Float<Exact,PR>::Float(Integer const& z)
    : _v(z.get_si())
{
    int n=z.get_si();
    ARIADNE_PRECONDITION(z==n);
}

template<class PR> Float<Exact,PR>::Float(Integer const& z, PR pr)
    : _v(z.get_si(),pr)
{
    int n=z.get_si();
    ARIADNE_PRECONDITION(z==n);
}

template<class PR> Float<Bounded,PR>::Float(Number<Validated> const& x)
    : Float(x.get(Bounded(),FLT::get_default_precision())) {
}

template<class PR> Float<Bounded,PR>::Float(Number<Validated> const& x, PR pr)
    : Float(x.get(Bounded(),pr)) {
}

template<class PR> Float<Upper,PR>::Float(Number<Upper> const& x, PR pr)
    : Float(x.get(Upper(),pr)) {
}

template<class PR> Float<Lower,PR>::Float(Number<Lower> const& x, PR pr)
    : Float(x.get(Lower(),pr)) {
}

template<class PR> Float<Approximate,PR>::Float(Number<Approximate> const& x, PR pr)
    : Float(x.get(Approximate(),pr)) {
}


template<class PR> Float<Bounded,PR>::Float(const Dyadic& b) : Float<Bounded,PR>(b.operator Rational()) { }

template<class PR> Float<Bounded,PR>::Float(const Decimal& d) : Float<Bounded,PR>(d.operator Rational()) { }


template<class PR> Float<Bounded,PR>::Float(const Integer& z) : Float<Bounded,PR>(Rational(z)) {
}

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

//Float<Exact,PR> inf = Float<Exact,PR>(std::numeric_limits< double >::infinity());
template<class PR> Float<Approximate,PR>::Float(Dyadic const& b) : Float<Approximate,PR>(b.operator Rational()) { }
template<class PR> Float<Approximate,PR>::Float(Decimal const& d) : Float<Approximate,PR>(d.operator Rational()) { }

template<class PR> Float<Approximate,PR>::Float(Number<Approximate> const& x) { ARIADNE_NOT_IMPLEMENTED; }

template<class PR> Float<Exact,PR>::operator Rational() const {
    return Rational(this->get_d());
}

template<class PR> Float<Approximate,PR>::Float(Rational const& q) : Float<Approximate,PR>(q.get_d()) {
}






ApproximateFloat64 floor(ApproximateFloat64 const& x) { return ApproximateFloat64(floor(x._a)); }
ApproximateFloat64 ceil(ApproximateFloat64 const& x) { return ApproximateFloat64(ceil(x._a)); }

ApproximateFloat64 abs(ApproximateFloat64 const& x) { return ApproximateFloat64(abs_exact(x._a)); }
ApproximateFloat64 max(ApproximateFloat64 const& x, ApproximateFloat64 y) { return ApproximateFloat64(max_exact(x._a,y._a)); }
ApproximateFloat64 min(ApproximateFloat64 const& x, ApproximateFloat64 y) { return ApproximateFloat64(min_exact(x._a,y._a)); }

ApproximateFloat64 nul(ApproximateFloat64 const& x) { return ApproximateFloat64(nul_exact(x._a)); }
ApproximateFloat64 pos(ApproximateFloat64 const& x) { return ApproximateFloat64(pos_exact(x._a)); }
ApproximateFloat64 neg(ApproximateFloat64 const& x) { return ApproximateFloat64(neg_exact(x._a)); }
ApproximateFloat64 half(ApproximateFloat64 const& x) { return ApproximateFloat64(half_exact(x._a)); }
ApproximateFloat64 sqr(ApproximateFloat64 const& x) { return ApproximateFloat64(mul_near(x._a,x._a)); }
ApproximateFloat64 rec(ApproximateFloat64 const& x) { return ApproximateFloat64(div_near(1.0,x._a)); }

ApproximateFloat64 add(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(add_near(x1._a,x2._a)); }
ApproximateFloat64 sub(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(sub_near(x1._a,x2._a)); }
ApproximateFloat64 mul(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(mul_near(x1._a,x2._a)); }
ApproximateFloat64 div(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return ApproximateFloat64(div_near(x1._a,x2._a)); }

ApproximateFloat64 pow(ApproximateFloat64 const& x, Nat m) { return ApproximateFloat64(pow_approx(x._a,m)); }
ApproximateFloat64 pow(ApproximateFloat64 const& x, Int n) { return ApproximateFloat64(pow_approx(x._a,n)); }

ApproximateFloat64 sqrt(ApproximateFloat64 const& x) { return ApproximateFloat64(sqrt_approx(x._a)); }
ApproximateFloat64 exp(ApproximateFloat64 const& x) { return ApproximateFloat64(exp_approx(x._a)); }
ApproximateFloat64 log(ApproximateFloat64 const& x) { return ApproximateFloat64(log_approx(x._a)); }
ApproximateFloat64 sin(ApproximateFloat64 const& x) { return ApproximateFloat64(sin_approx(x._a)); }
ApproximateFloat64 cos(ApproximateFloat64 const& x) { return ApproximateFloat64(cos_approx(x._a)); }
ApproximateFloat64 tan(ApproximateFloat64 const& x) { return ApproximateFloat64(tan_approx(x._a)); }
ApproximateFloat64 asin(ApproximateFloat64 const& x) { return ApproximateFloat64(asin_approx(x._a)); }
ApproximateFloat64 acos(ApproximateFloat64 const& x) { return ApproximateFloat64(acos_approx(x._a)); }
ApproximateFloat64 atan(ApproximateFloat64 const& x) { return ApproximateFloat64(atan_approx(x._a)); }

ApproximateFloat64 operator+(ApproximateFloat64 const& x) { return pos(x); }
ApproximateFloat64 operator-(ApproximateFloat64 const& x) { return neg(x); }
ApproximateFloat64 operator+(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return add(x1,x2); }
ApproximateFloat64 operator-(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return sub(x1,x2); }
ApproximateFloat64 operator*(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return mul(x1,x2); }
ApproximateFloat64 operator/(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return div(x1,x2); }
ApproximateFloat64& operator+=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1._a+=x2._a; return x1; }
ApproximateFloat64& operator-=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1._a-=x2._a; return x1; }
ApproximateFloat64& operator*=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1._a*=x2._a; return x1; }
ApproximateFloat64& operator/=(ApproximateFloat64& x1, ApproximateFloat64 const& x2) { x1._a/=x2._a; return x1; }

Bool operator==(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1._a==x2._a; }
Bool operator!=(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1._a!=x2._a; }
Bool operator<=(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1._a<=x2._a; }
Bool operator>=(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1._a>=x2._a; }
Bool operator< (ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1._a< x2._a; }
Bool operator> (ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) { return x1._a> x2._a; }

OutputStream& operator<<(OutputStream& os, ApproximateFloat64 const& x);
InputStream& operator>>(InputStream& is, ApproximateFloat64& x);



LowerFloat64 max(LowerFloat64 const& x1, LowerFloat64 const& x2) { return LowerFloat64(max_exact(x1._l,x2._l)); }
LowerFloat64 min(LowerFloat64 const& x1, LowerFloat64 const& x2) { return LowerFloat64(min_exact(x1._l,x2._l)); }

LowerFloat64 nul(LowerFloat64 const& x) { return LowerFloat64(pos_exact(x._l)); }
LowerFloat64 pos(LowerFloat64 const& x) { return LowerFloat64(pos_exact(x._l)); }
LowerFloat64 neg(UpperFloat64 const& x) { return LowerFloat64(neg_exact(x._u)); }
LowerFloat64 half(LowerFloat64 const& x) { return LowerFloat64(half_exact(x._l)); }
LowerFloat64 sqr(LowerFloat64 const& x);
LowerFloat64 rec(UpperFloat64 const& x);

LowerFloat64 add(LowerFloat64 const& x1, LowerFloat64 const& x2) { return LowerFloat64(add_down(x1._l,x2._l)); }
LowerFloat64 sub(LowerFloat64 const& x1, UpperFloat64 const& x2) { return LowerFloat64(sub_down(x1._l,x2._u)); }
LowerFloat64 mul(LowerFloat64 const& x1, LowerFloat64 const& x2);
LowerFloat64 div(LowerFloat64 const& x1, UpperFloat64 const& x2);
LowerFloat64 pow(LowerFloat64 const& x, Nat m);

LowerFloat64 sqrt(LowerFloat64 const& x);
LowerFloat64 exp(LowerFloat64 const& x);
LowerFloat64 log(LowerFloat64 const& x);
LowerFloat64 atan(LowerFloat64 const& x);

LowerFloat64 operator+(LowerFloat64 const& x) { return pos(x); }
LowerFloat64 operator-(UpperFloat64 const& x) { return neg(x); }
LowerFloat64 operator+(LowerFloat64 const& x1, LowerFloat64 const& x2) { return add(x1,x2); }
LowerFloat64 operator-(LowerFloat64 const& x1, UpperFloat64 const& x2) { return sub(x1,x2); }
LowerFloat64 operator*(LowerFloat64 const& x1, LowerFloat64 const& x2) { return mul(x1,x2); }
LowerFloat64 operator/(LowerFloat64 const& x1, UpperFloat64 const& x2) { return div(x1,x2); }
LowerFloat64& operator+=(LowerFloat64& x1, LowerFloat64 const& x2) { return x1=x1+x2; }
LowerFloat64& operator-=(LowerFloat64& x1, UpperFloat64 const& x2) { return x1=x1-x2; }
LowerFloat64& operator*=(LowerFloat64& x1, LowerFloat64 const& x2) { return x1=x1*x2; }
LowerFloat64& operator/=(LowerFloat64& x1, UpperFloat64 const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, LowerFloat64 const& x);
InputStream& operator>>(InputStream& is, LowerFloat64& x);

UpperFloat64 max(UpperFloat64 const& x1, UpperFloat64 const& x2) { return UpperFloat64(max_exact(x1._u,x2._u)); }
UpperFloat64 min(UpperFloat64 const& x1, UpperFloat64 const& x2) { return UpperFloat64(min_exact(x1._u,x2._u)); }

UpperFloat64 nul(UpperFloat64 const& x) { return UpperFloat64(pos_exact(x._u)); }
UpperFloat64 pos(UpperFloat64 const& x) { return UpperFloat64(pos_exact(x._u)); }
UpperFloat64 neg(LowerFloat64 const& x) { return UpperFloat64(neg_exact(x._l)); }
UpperFloat64 half(UpperFloat64 const& x) { return UpperFloat64(half_exact(x._u)); }
UpperFloat64 sqr(UpperFloat64 const& x);
UpperFloat64 rec(LowerFloat64 const& x);

UpperFloat64 add(UpperFloat64 const& x1, UpperFloat64 const& x2) { return UpperFloat64(add_up(x1._u,x2._u)); }
UpperFloat64 sub(UpperFloat64 const& x1, LowerFloat64 const& x2) { return UpperFloat64(sub_up(x1._u,x2._l)); }
UpperFloat64 mul(UpperFloat64 const& x1, UpperFloat64 const& x2);
UpperFloat64 div(UpperFloat64 const& x1, LowerFloat64 const& x2);
UpperFloat64 pow(UpperFloat64 const& x, Nat m);

UpperFloat64 sqrt(UpperFloat64 const& x);
UpperFloat64 exp(UpperFloat64 const& x);
UpperFloat64 log(UpperFloat64 const& x);
UpperFloat64 atan(UpperFloat64 const& x);

UpperFloat64 operator+(UpperFloat64 const& x) { return pos(x); }
UpperFloat64 operator-(LowerFloat64 const& x) { return neg(x); }
UpperFloat64 operator+(UpperFloat64 const& x1, UpperFloat64 const& x2) { return add(x1,x2); }
UpperFloat64 operator-(UpperFloat64 const& x1, LowerFloat64 const& x2) { return sub(x1,x2); }
UpperFloat64 operator*(UpperFloat64 const& x1, UpperFloat64 const& x2) { return mul(x1,x2); }
UpperFloat64 operator/(UpperFloat64 const& x1, LowerFloat64 const& x2) { return div(x1,x2); }
UpperFloat64& operator+=(UpperFloat64& x1, UpperFloat64 const& x2) { return x1=x1+x2; }
UpperFloat64& operator-=(UpperFloat64& x1, LowerFloat64 const& x2) { return x1=x1-x2; }
UpperFloat64& operator*=(UpperFloat64& x1, UpperFloat64 const& x2) { return x1=x1*x2; }
UpperFloat64& operator/=(UpperFloat64& x1, LowerFloat64 const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, UpperFloat64 const& x);
InputStream& operator>>(InputStream& is, UpperFloat64& x);





ValidatedFloat64 max(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 min(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 abs(ValidatedFloat64 const& x);

ValidatedFloat64 nul(ValidatedFloat64 const& x);
ValidatedFloat64 pos(ValidatedFloat64 const& x);
ValidatedFloat64 neg(ValidatedFloat64 const& x);
ValidatedFloat64 half(ValidatedFloat64 const& x);
ValidatedFloat64 sqr(ValidatedFloat64 const& x);
ValidatedFloat64 rec(ValidatedFloat64 const& x);

ValidatedFloat64 add(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 sub(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 mul(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 div(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2);
ValidatedFloat64 pow(ValidatedFloat64 const& x, Nat m);
ValidatedFloat64 pow(ValidatedFloat64 const& x, Int m);

ValidatedFloat64 sqrt(ValidatedFloat64 const& x);
ValidatedFloat64 exp(ValidatedFloat64 const& x);
ValidatedFloat64 log(ValidatedFloat64 const& x);
ValidatedFloat64 sin(ValidatedFloat64 const& x);
ValidatedFloat64 cos(ValidatedFloat64 const& x);
ValidatedFloat64 tan(ValidatedFloat64 const& x);
ValidatedFloat64 asin(ValidatedFloat64 const& x);
ValidatedFloat64 acos(ValidatedFloat64 const& x);
ValidatedFloat64 atan(ValidatedFloat64 const& x);

ValidatedFloat64 operator+(ValidatedFloat64 const& x) { return pos(x); }
ValidatedFloat64 operator-(ValidatedFloat64 const& x) { return neg(x); }
ValidatedFloat64 operator+(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return add(x1,x2); }
ValidatedFloat64 operator-(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return sub(x1,x2); }
ValidatedFloat64 operator*(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return mul(x1,x2); }
ValidatedFloat64 operator/(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return div(x1,x2); }
ValidatedFloat64& operator+=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1+x2; }
ValidatedFloat64& operator-=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1-x2; }
ValidatedFloat64& operator*=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1*x2; }
ValidatedFloat64& operator/=(ValidatedFloat64& x1, ValidatedFloat64 const& x2) { return x1=x1/x2; }

OutputStream& operator<<(OutputStream& os, ValidatedFloat64 const& x);
InputStream& operator>>(InputStream& is, ValidatedFloat64& x);



extern const ExactFloat64 infty;
ExactFloat64 operator"" _exact(long double lx) { double x=lx; assert(x==lx); return ExactFloat64(x); }
TwoExp::operator ExactFloat64 () const { return ExactFloat64(this->get_d()); }

ExactFloat64 max(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return ExactFloat64(max(x1._v,x2._v)); }
ExactFloat64 min(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return ExactFloat64(min(x1._v,x2._v)); }
ExactFloat64 abs(ExactFloat64 const& x) { return ExactFloat64(abs(x._v)); }

ExactFloat64 nul(ExactFloat64 const& x) { return ExactFloat64(nul(x._v)); }
ExactFloat64 pos(ExactFloat64 const& x) { return ExactFloat64(pos(x._v)); }
ExactFloat64 neg(ExactFloat64 const& x) { return ExactFloat64(neg(x._v)); }
ExactFloat64 half(ExactFloat64 const& x) { return ExactFloat64(half(x._v)); }

ValidatedFloat64 sqr(ExactFloat64 const& x);
ValidatedFloat64 rec(ExactFloat64 const& x);
ValidatedFloat64 add(ExactFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 sub(ExactFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 mul(ExactFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 div(ExactFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 pow(ExactFloat64 const& x, Int n);

ValidatedFloat64 sqrt(ExactFloat64 const& x);
ValidatedFloat64 exp(ExactFloat64 const& x);
ValidatedFloat64 log(ExactFloat64 const& x);
ValidatedFloat64 sin(ExactFloat64 const& x);
ValidatedFloat64 cos(ExactFloat64 const& x);
ValidatedFloat64 tan(ExactFloat64 const& x);
ValidatedFloat64 atan(ExactFloat64 const& x);

ValidatedFloat64 rad(ExactFloat64 const& x1, ExactFloat64 const& x2);
ValidatedFloat64 med(ExactFloat64 const& x1, ExactFloat64 const& x2);

ExactFloat64 operator+(ExactFloat64 const& x) { return pos(x); }
ExactFloat64 operator-(ExactFloat64 const& x) { return neg(x); }
ValidatedFloat64 operator+(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return add(x1,x2); }
ValidatedFloat64 operator-(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return sub(x1,x2); }
ValidatedFloat64 operator*(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return mul(x1,x2); }
ValidatedFloat64 operator/(ExactFloat64 const& x1,  ExactFloat64 const& x2) { return div(x1,x2); }

/*
ValidatedFloat64 operator+(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 operator-(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 operator*(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 operator/(ValidatedFloat64 const& x1,  ExactFloat64 const& x2);
ValidatedFloat64 operator+(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
ValidatedFloat64 operator-(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
ValidatedFloat64 operator*(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
ValidatedFloat64 operator/(ExactFloat64 const& x1,  ValidatedFloat64 const& x2);
*/

ExactFloat64 operator*(ExactFloat64 const& x, TwoExp y) { ExactFloat64 yv=y; return ExactFloat64(x.raw()*yv.raw()); }
ExactFloat64 operator/(ExactFloat64 const& x, TwoExp y) { ExactFloat64 yv=y; return ExactFloat64(x.raw()/yv.raw()); }
ExactFloat64& operator*=(ExactFloat64& x, TwoExp y) { ExactFloat64 yv=y; return x=ExactFloat64(x.raw()*yv.raw()); }
ExactFloat64& operator/=(ExactFloat64& x, TwoExp y) { ExactFloat64 yv=y; return x=ExactFloat64(x.raw()/yv.raw()); }

Boolean operator==(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()==x2.raw(); }
Boolean operator!=(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()!=x2.raw(); }
Boolean operator<=(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()<=x2.raw(); }
Boolean operator>=(ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()>=x2.raw(); }
Boolean operator< (ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()< x2.raw(); }
Boolean operator> (ExactFloat64 const& x1, ExactFloat64 const& x2) { return x1.raw()> x2.raw(); }

OutputStream& operator<<(OutputStream& os, ExactFloat64 const& x);


ExactFloat64 make_exact(const Real& x);

Bool operator==(ExactFloat64 const& x, const Rational& q) { return Rational(x)==q; }
Bool operator!=(ExactFloat64 const& x, const Rational& q) { return Rational(x)!=q; }
Bool operator<=(ExactFloat64 const& x, const Rational& q) { return Rational(x)<=q; }
Bool operator>=(ExactFloat64 const& x, const Rational& q) { return Rational(x)>=q; }
Bool operator< (ExactFloat64 const& x, const Rational& q) { return Rational(x)< q; }
Bool operator> (ExactFloat64 const& x, const Rational& q) { return Rational(x)> q; }

Bool operator==(const Rational& q, ExactFloat64 const& x) { return q==Rational(x); }
Bool operator!=(const Rational& q, ExactFloat64 const& x) { return q!=Rational(x); }
Bool operator<=(const Rational& q, ExactFloat64 const& x) { return q<=Rational(x); }
Bool operator>=(const Rational& q, ExactFloat64 const& x) { return q>=Rational(x); }
Bool operator< (const Rational& q, ExactFloat64 const& x) { return q< Rational(x); }
Bool operator> (const Rational& q, ExactFloat64 const& x) { return q> Rational(x); }


PositiveExactFloat64 mag(ExactFloat64 const& x) {
    return PositiveExactFloat64(abs(x.raw())); }
// FIXME: Unsafe since x may be negative
PositiveUpperFloat64 mag(UpperFloat64 const& x) {
    return PositiveUpperFloat64(abs(x.raw())); }
PositiveUpperFloat64 mag(ValidatedFloat64 const& x) {
    return PositiveUpperFloat64(max(neg(x.lower_raw()),x.upper_raw())); }
PositiveLowerFloat64 mig(ValidatedFloat64 const& x) {
    Float64 r=max(x.lower_raw(),neg(x.upper_raw()));
    return PositiveLowerFloat64(max(r,nul(r))); }
PositiveApproximateFloat64 mag(ApproximateFloat64 const& x) {
    return PositiveApproximateFloat64(abs(x.raw())); }


ValidatedFloat64 make_bounds(PositiveUpperFloat64 const& _e) {
    return ValidatedFloat64(-_e.raw(),+_e.raw());
}

ExactFloat64 value(ValidatedFloat64 const& x) {
    return ExactFloat64(half_exact(add_near(x.lower_raw(),x.upper_raw())));
}

PositiveUpperFloat64 error(ValidatedFloat64 const& x) {
    return PositiveUpperFloat64(half_exact(sub_up(x.upper_raw(),x.lower_raw())));
}



template<class PR> Float<Approximate,PR>::Float(Float<Lower,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Upper,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Bounded,PR> const& x) : _a(x.value_raw()) {
}

template<class PR> Float<Approximate,PR>::Float(Float<Exact,PR> const& x) : _a(x.raw()) {
}

template<class PR> Float<Lower,PR>::Float(Float<Bounded,PR> const& x) : _l(x.lower_raw()) {
}

template<class PR> Float<Lower,PR>::Float(Float<Exact,PR> const& x) : _l(x.raw()) {
}

template<class PR> Float<Upper,PR>::Float(Float<Bounded,PR> const& x) : _u(x.upper_raw()) {
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


Bool same(ApproximateFloat64 const& x1, ApproximateFloat64 const& x2) {
    return x1.raw()==x2.raw(); }

Bool same(LowerFloat64 const& x1, LowerFloat64 const& x2) {
    return x1.raw()==x2.raw(); }

Bool same(UpperFloat64 const& x1, UpperFloat64 const& x2) {
    return x1.raw()==x2.raw(); }

Bool same(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }

Bool same(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    return x1.raw()==x2.raw(); }





ExactFloat64 midpoint(ValidatedFloat64 const& x) { return x.value(); }






ValidatedFloat64 max(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    return ValidatedFloat64(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
}

ValidatedFloat64 min(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    return ValidatedFloat64(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
}


ValidatedFloat64 abs(ValidatedFloat64 const& x)
{
    if(x.lower_raw()>=0) {
        return ValidatedFloat64(x.lower_raw(),x.upper_raw());
    } else if(x.upper_raw()<=0) {
        return ValidatedFloat64(neg(x.upper_raw()),neg(x.lower_raw()));
    } else {
        return ValidatedFloat64(static_cast<Float64>(0.0),max(neg(x.lower_raw()),x.upper_raw()));
    }
}

ValidatedFloat64 pos(ValidatedFloat64 const& x)
{
    return ValidatedFloat64(pos(x.lower_raw()),pos(x.upper_raw()));
}

ValidatedFloat64 neg(ValidatedFloat64 const& x)
{
    return ValidatedFloat64(neg(x.upper_raw()),neg(x.lower_raw()));
}

ValidatedFloat64 half(ValidatedFloat64 const& x) {
    return ValidatedFloat64(half(x.lower_raw()),half(x.upper_raw()));
}

ValidatedFloat64 sqr(ExactFloat64 const& x)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& xv=x.raw();
    Float64::set_rounding_downward();
    Float64 rl=mul(xv,xv);
    Float64::set_rounding_upward();
    Float64 ru=mul(xv,xv);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 rec(ExactFloat64 const& x)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& xv=x.raw();
    Float64::set_rounding_downward();
    Float64 rl=1.0/xv;
    Float64::set_rounding_upward();
    Float64 ru=1.0/xv;
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}



ValidatedFloat64 add(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2l=x2.lower_raw();
    Float64 const& x2u=x2.upper_raw();
    Float64::set_rounding_downward();
    Float64 rl=add(x1l,x2l);
    Float64::set_rounding_upward();
    Float64 ru=add(x1u,x2u);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 add(ValidatedFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=add(x1l,x2v);
    Float64::set_rounding_upward();
    Float64 ru=add(x1u,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 add(ExactFloat64 const& x1, ValidatedFloat64 const& x2)
{
    return add(x2,x1);
}

ValidatedFloat64 add(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=add(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=add(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 sub(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2l=x2.lower_raw();
    Float64 const& x2u=x2.upper_raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1l,x2u);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1u,x2l);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 sub(ValidatedFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1l=x1.lower_raw();
    Float64 const& x1u=x1.upper_raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1l,x2v);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1u,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 sub(ExactFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2l=x2.lower_raw();
    Float64 const& x2u=x2.upper_raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1v,x2u);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1v,x2l);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 sub(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=sub(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=sub(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 mul(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=mul(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=mul(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 div(ExactFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64 const& x1v=x1.raw();
    Float64 const& x2v=x2.raw();
    Float64::set_rounding_downward();
    Float64 rl=div(x1v,x2v);
    Float64::set_rounding_upward();
    Float64 ru=div(x1v,x2v);
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 pow(ExactFloat64 const& x1, Int n2) {
    return pow(ValidatedFloat64(x1),n2);
}

ValidatedFloat64 med(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    return add(half(x1),half(x2));
}

ValidatedFloat64 rad(ExactFloat64 const& x1, ExactFloat64 const& x2) {
    return sub(half(x2),half(x1));
}

ValidatedFloat64 sqrt(ExactFloat64 const& x) {
    return sqrt(ValidatedFloat64(x));
}

ValidatedFloat64 exp(ExactFloat64 const& x) {
    return exp(ValidatedFloat64(x));
}

ValidatedFloat64 log(ExactFloat64 const& x) {
    return log(ValidatedFloat64(x));
}

ValidatedFloat64 sin(ExactFloat64 const& x) {
    return sin(ValidatedFloat64(x));
}

ValidatedFloat64 cos(ExactFloat64 const& x) {
    return cos(ValidatedFloat64(x));
}



ValidatedFloat64 med(ValidatedFloat64 const& x);

ValidatedFloat64 rad(ValidatedFloat64 const& x);


/*
ValidatedFloat64 operator+(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return add(x1,x2); }
ValidatedFloat64 operator-(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return sub(x1,x2); }
ValidatedFloat64 operator*(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return mul(x1,x2); }
ValidatedFloat64 operator/(ValidatedFloat64 const& x1, ExactFloat64 const& x2) { return div(x1,x2); }
ValidatedFloat64 operator+(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return add(x2,x1); }
ValidatedFloat64 operator-(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return sub(x1,x2); }
ValidatedFloat64 operator*(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return mul(x2,x1); }
ValidatedFloat64 operator/(ExactFloat64 const& x1, ValidatedFloat64 const& x2) { return div(x1,x2); }
*/

// Standard equality operators
//! \related ValidatedFloat64 \brief Tests if \_a x1 provides tighter bounds than \_a x2.
Bool refines(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return x1.lower_raw()>=x2.lower_raw() && x1.upper_raw()<=x2.upper_raw(); }

//! \related ValidatedFloat64 \brief The common refinement of \_a x1 and \_a x2.
ValidatedFloat64 refinement(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return ValidatedFloat64(max(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw())); }

//! \related ValidatedFloat64 \brief Tests if \_a x1 and \_a x2 are consistent with representing the same number.
Bool consistent(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return x1.lower_raw()<=x2.upper_raw() && x1.upper_raw()>=x2.lower_raw(); }

//! \related ValidatedFloat64 \brief  Tests if \_a x1 and \_a x2 are inconsistent with representing the same number.
Bool inconsistent(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    return not consistent(x1,x2); }

//! \related ValidatedFloat64 \brief  Tests if \_a x1 is a model for the exact value \_a x2. number.
Bool models(ValidatedFloat64 const& x1, ExactFloat64 const& x2) {
    return x1.lower_raw()<=x2.raw() && x1.upper_raw()>=x2.raw(); }

// Standard equality operators
//! \related ValidatedFloat64 \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
Bool operator==(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }
//! \related ValidatedFloat64 \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
//! even though the intervals possibly represent the same exact real value.
Bool operator!=(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) { return x1.lower_raw()!=x2.lower_raw() || x1.upper_raw()!=x2.upper_raw(); }



//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
//! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
Tribool operator> (ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.lower_raw()> x2.upper_raw()) { return true; }
    else if(x1.upper_raw()<=x2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
Tribool operator< (ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.upper_raw()< x2.lower_raw()) { return true; }
    else if(x1.lower_raw()>=x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
Tribool operator>=(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.lower_raw()>=x2.upper_raw()) { return true; }
    else if(x1.upper_raw()< x2.lower_raw()) { return false; }
    else { return indeterminate; }
}

//! \related ValidatedFloat64 \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
Tribool operator<=(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2) {
    if(x1.upper_raw()<=x2.lower_raw()) { return true; }
    else if(x1.lower_raw()> x2.upper_raw()) { return false; }
    else { return indeterminate; }
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> Void serialize(A& _a, ValidatedFloat64& ivl, const Nat version) {
    _a & ivl.lower_raw() & ivl.upper_raw(); }
#endif

OutputStream& operator<<(OutputStream&, ValidatedFloat64 const&);
InputStream& operator>>(InputStream&, ValidatedFloat64&);


PositiveUpperFloat64 operator+(PositiveUpperFloat64 const& x1, PositiveUpperFloat64 const& x2);
PositiveUpperFloat64 operator-(PositiveUpperFloat64 const& x1, LowerFloat64 const& x2);
PositiveUpperFloat64 operator*(PositiveUpperFloat64 const& x1, PositiveUpperFloat64 const& x2);
PositiveUpperFloat64 operator/(PositiveUpperFloat64 const& x1, LowerFloat64 const& x2);
PositiveUpperFloat64 pow(PositiveUpperFloat64 const& x, Nat m);
PositiveUpperFloat64 half(PositiveUpperFloat64 const& x);


ErrorFloat64 operator"" _error(long double lx) { double x=lx; assert(x==lx); return ErrorFloat64(Float64(x)); }

const ExactFloat64 infty = ExactFloat64(Float64::inf());


ExactFloat64 make_exact(Real const& r) {
    ApproximateFloat64 _a(r); return ExactFloat64(_a.raw());
}

OutputStream& operator<<(OutputStream& os, ExactFloat64 const& x) {
    os << std::showpoint << std::setprecision(ExactFloat64::output_precision) << x.raw();
    return os;
}


ValidatedFloat64 widen(ValidatedFloat64 const& x)
{
    Float64::RoundingModeType rm=Float64::get_rounding_mode();
    const Float64& xl=x.lower_raw();
    const Float64& xu=x.upper_raw();
    const Float64 m=std::numeric_limits<float>::min();
    Float64::set_rounding_upward();
    Float64 wu=add(xu,m);
    Float64 mwl=add(neg(xl),m);
    Float64 wl=neg(mwl);
    Float64::set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return ValidatedFloat64(wl,wu);
}

ValidatedFloat64 narrow(ValidatedFloat64 const& x)
{
    Float64::RoundingModeType rm=Float64::get_rounding_mode();
    const Float64& xl=x.lower_raw();
    const Float64& xu=x.upper_raw();
    const Float64 m=std::numeric_limits<float>::min();
    Float64::set_rounding_upward();
    Float64 mnu=add(neg(xu),m);
    Float64 nu=neg(mnu);
    Float64 nl=add(xl,m);
    Float64::set_rounding_mode(rm);
    assert(xl<nl); assert(nu<xu);
    return ValidatedFloat64(nl,nu);
}

ValidatedFloat64 trunc(ValidatedFloat64 const& x)
{
    Float64::RoundingModeType rm=Float64::get_rounding_mode();
    const double& xl=x.lower_raw().get_d();
    const double& xu=x.upper_raw().get_d();
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { Float64::set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { Float64::set_rounding_downward(); tl-=fm; }
    Float64::set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return ValidatedFloat64(double(tl),double(tu));
}

ValidatedFloat64 trunc(ValidatedFloat64 const& x, Nat n)
{
    ValidatedFloat64 _e=ValidatedFloat64(std::pow(2.0,52-(Int)n));
    ValidatedFloat64 y=x+_e;
    return y-_e;
}

ValidatedFloat64 rec(ValidatedFloat64 const& x)
{
    const Float64& xl=x.lower_raw();
    const Float64& xu=x.upper_raw();
    Float64 rl,ru;
    if(xl>0 || xu<0) {
        rl=rec_down(xu);
        ru=rec_up(xl);
    } else {
        rl=-inf;
        ru=+inf;
        ARIADNE_THROW(DivideByZeroException,"ValidatedFloat64 rec(ValidatedFloat64 ivl)","ivl="<<x);
    }
    return ValidatedFloat64(rl,ru);
}


ValidatedFloat64 mul(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    const Float64& x1l=x1.lower_raw();
    const Float64& x1u=x1.upper_raw();
    const Float64& i2l=x2.lower_raw();
    const Float64& i2u=x2.upper_raw();
    Float64 rl,ru;
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    if(x1l>=0) {
        if(i2l>=0) {
            rl=mul_down(x1l,i2l); ru=mul_up(x1u,i2u);
        } else if(i2u<=0) {
            rl=mul_down(x1u,i2l); ru=mul_up(x1l,i2u);
        } else {
            rl=mul_down(x1u,i2l); ru=mul_up(x1u,i2u);
        }
    }
    else if(x1u<=0) {
        if(i2l>=0) {
            rl=mul_down(x1l,i2u); ru=mul_up(x1u,i2l);
        } else if(i2u<=0) {
            rl=mul_down(x1u,i2u); ru=mul_up(x1l,i2l);
        } else {
            rl=mul_down(x1l,i2u); ru=mul_up(x1l,i2l);
        }
    } else {
        if(i2l>=0) {
            rl=mul_down(x1l,i2u); ru=mul_up(x1u,i2u);
        } else if(i2u<=0) {
            rl=mul_down(x1u,i2l); ru=mul_up(x1l,i2l);
        } else {
            Float64::set_rounding_downward();
            rl=min(mul(x1u,i2l),mul(x1l,i2u));
            Float64::set_rounding_upward();
            ru=max(mul(x1l,i2l),mul(x1u,i2u));
        }
    }
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}


ValidatedFloat64 mul(ValidatedFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    const Float64& x1l=x1.lower_raw();
    const Float64& x1u=x1.upper_raw();
    const Float64& x2v=x2.raw();
    Float64 rl,ru;
    if(x2>=0) {
        rl=mul_down(x1l,x2v); ru=mul_up(x1u,x2v);
    } else {
        rl=mul_down(x1u,x2v); ru=mul_up(x1l,x2v);
    }
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}


ValidatedFloat64 mul(ExactFloat64 const& x1, ValidatedFloat64 const& x2) {
    return mul(x2,x1);
}


ValidatedFloat64 div(ValidatedFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    const Float64& x1l=x1.lower_raw();
    const Float64& x1u=x1.upper_raw();
    const Float64& i2l=x2.lower_raw();
    const Float64& i2u=x2.upper_raw();
    Float64 rl,ru;
    if(i2l>0) {
        if(x1l>=0) {
            rl=div_down(x1l,i2u); ru=div_up(x1u,i2l);
        } else if(x1u<=0) {
            rl=div_down(x1l,i2l); ru=div_up(x1u,i2u);
        } else {
            rl=div_down(x1l,i2l); ru=div_up(x1u,i2l);
        }
    }
    else if(i2u<0) {
        if(x1l>=0) {
            rl=div_down(x1u,i2u); ru=div_up(x1l,i2l);
        } else if(x1u<=0) {
            rl=div_down(x1u,i2l); ru=div_up(x1l,i2u);
        } else {
            rl=div_down(x1u,i2u); ru=div_up(x1l,i2u);
        }
    }
    else {
        // ARIADNE_THROW(DivideByZeroException,"ValidatedFloat64 div(ValidatedFloat64 ivl1, ValidatedFloat64 ivl2)","ivl1="<<x1<<", ivl2="<<x2);
        rl=-Float64::inf();
        ru=+Float64::inf();
    }
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}



ValidatedFloat64 div(ValidatedFloat64 const& x1, ExactFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    const Float64& x1l=x1.lower_raw();
    const Float64& x1u=x1.upper_raw();
    const Float64& x2v=x2.raw();
    Float64 rl,ru;
    if(x2v>0) {
        rl=div_down(x1l,x2v); ru=div_up(x1u,x2v);
    } else if(x2v<0) {
        rl=div_down(x1u,x2v); ru=div_up(x1l,x2v);
    } else {
        rl=-Float64::inf();
        ru=+Float64::inf();
    }
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}


ValidatedFloat64 div(ExactFloat64 const& x1, ValidatedFloat64 const& x2)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    const Float64& x1v=x1.raw();
    const Float64& i2l=x2.lower_raw();
    const Float64& i2u=x2.upper_raw();
    Float64 rl,ru;
    if(i2l<=0 && i2u>=0) {
        ARIADNE_THROW(DivideByZeroException,"ValidatedFloat64 div(Float64 const& x1, ValidatedFloat64 ivl2)","x1="<<x1<<", ivl2="<<x2);
        rl=-Float64::inf();
        ru=+Float64::inf();
    } else if(x1v>=0) {
        rl=div_down(x1v,i2u); ru=div_up(x1v,i2l);
    } else {
        rl=div_down(x1v,i2l); ru=div_up(x1v,i2u);
    }
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 sqr(ValidatedFloat64 const& x)
{
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    const Float64& xl=x.lower_raw();
    const Float64& xu=x.upper_raw();
    Float64 rl,ru;
    if(xl>0.0) {
        rl=mul_down(xl,xl);
        ru=mul_up(xu,xu);
    } else if(xu<0.0) {
        rl=mul_down(xu,xu);
        ru=mul_up(xl,xl);
    } else {
        rl=nul(xl);
        Float64 ru1=mul_up(xl,xl);
        Float64 ru2=mul_up(xu,xu);
        ru=max(ru1,ru2);
    }
    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}




ValidatedFloat64 pow(ValidatedFloat64 const& x, Int n) {
    if(n<0) { return pow(rec(x),Nat(-n)); }
    else return pow(x,Nat(n));
}

ValidatedFloat64 pow(ValidatedFloat64 const& x, Nat m) {
    ValidatedFloat64 y = x;
    if(m%2==0) { y=abs(x); }
    Float64 rl=pow_down(y.lower_raw(),m);
    Float64 ru=pow_up(y.upper_raw(),m);
    return ValidatedFloat64(rl,ru);
}



ValidatedFloat64 sqrt(ValidatedFloat64 const& x) {
    return ValidatedFloat64(sqrt_down(x.lower_raw()),sqrt_up(x.upper_raw()));
}

ValidatedFloat64 exp(ValidatedFloat64 const& x) {
    return ValidatedFloat64(exp_down(x.lower_raw()),exp_up(x.upper_raw()));
}

ValidatedFloat64 log(ValidatedFloat64 const& x) {
    return ValidatedFloat64(log_down(x.lower_raw()),log_up(x.upper_raw()));
}


ValidatedFloat64 pi_val() { return ValidatedFloat64(pi_down(),pi_up()); }


ValidatedFloat64 sin(ValidatedFloat64 const& x)
{
    return cos(x-half(pi_val()));
}

ValidatedFloat64 cos(ValidatedFloat64 const& x)
{
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    Float64::RoundingModeType rnd = Float64::get_rounding_mode();

    static const ExactFloat64 two(2);

    if(x.error().raw()>2*pi_down()) { return ValidatedFloat64(-1.0,+1.0); }

    Float64 n=floor(x.lower_raw()/(2*pi_approx())+0.5);
    ValidatedFloat64 y=x-two*ExactFloat64(n)*pi_val();

    ARIADNE_ASSERT(y.lower_raw()<=pi_up());
    ARIADNE_ASSERT(y.upper_raw()>=-pi_up());

    Float64 rl,ru;
    if(y.lower_raw()<=-pi_down()) {
        if(y.upper_raw()<=0.0) { rl=-1.0; ru=cos_up(y.upper_raw()); }
        else { rl=-1.0; ru=+1.0; }
    } else if(y.lower_raw()<=0.0) {
        if(y.upper_raw()<=0.0) { rl=cos_down(y.lower_raw()); ru=cos_up(y.upper_raw()); }
        else if(y.upper_raw()<=pi_down()) { rl=cos_down(max(-y.lower_raw(),y.upper_raw())); ru=+1.0; }
        else { rl=-1.0; ru=+1.0; }
    } else if(y.lower_raw()<=pi_up()) {
        if(y.upper_raw()<=pi_down()) { rl=cos_down(y.upper_raw()); ru=cos_up(y.lower_raw()); }
        else if(y.upper_raw()<=2*pi_down()) { rl=-1.0; ru=cos_up(min(y.lower_raw(),sub_down(2*pi_down(),y.upper_raw()))); }
        else { rl=-1.0; ru=+1.0; }
    } else {
        assert(false);
    }

    Float64::set_rounding_mode(rnd);
    return ValidatedFloat64(rl,ru);
}

ValidatedFloat64 tan(ValidatedFloat64 const& x) {
    return mul(sin(x),rec(cos(x)));
}

ValidatedFloat64 asin(ValidatedFloat64 const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat64 acos(ValidatedFloat64 const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat64 atan(ValidatedFloat64 const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}


OutputStream&
operator<<(OutputStream& os, const ValidatedFloat64& ivl)
{
    //if(ivl.lower_raw()==ivl.upper_raw()) { return os << "{" << std::setprecision(ValidatedFloat64::output_precision) << ivl.lower_raw().get_d() << ; }
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    os << '{';
    Float64::set_rounding_downward();
    os << std::showpoint << std::setprecision(ValidatedFloat64::output_precision) << ivl.lower().get_d();
    os << ':';
    Float64::set_rounding_upward();
    os << std::showpoint << std::setprecision(ValidatedFloat64::output_precision) << ivl.upper().get_d();
    Float64::set_rounding_mode(rnd);
    os << '}';
    return os;

}

InputStream&
operator>>(InputStream& is, ValidatedFloat64& x)
{
    Float64 _l,_u;
    char cl,cm,cr;
    is >> cl >> _l >> cm >> _u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    x._l=_l; x._u=_u;
    return is;
}


OutputStream& operator<<(OutputStream& os, UpperFloat64 const& x) {
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64::set_rounding_upward();
    os << std::showpoint << std::setprecision(ValidatedFloat64::output_precision) << x.raw();
    Float64::set_rounding_mode(rnd);
    return os;
}

UpperFloat64 sqr(UpperFloat64 const& x) {
    ARIADNE_PRECONDITION(x.raw()>=0.0);
    return UpperFloat64(mul_up(x.raw(),x.raw()));
}

UpperFloat64 mul(UpperFloat64 const& x1, UpperFloat64 const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0.0 && x2.raw() >= 0.0);
    return UpperFloat64(mul_up(x1.raw(),x2.raw()));
}

UpperFloat64 div(UpperFloat64 const& x1, LowerFloat64 const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0.0 && x2.raw() > 0.0);
    return UpperFloat64(div_up(x1.raw(),x2.raw()));
}

UpperFloat64 pow(UpperFloat64 const& x, Nat n) {
    ARIADNE_PRECONDITION(x.raw()>=0.0);
    return UpperFloat64(pow_up(x.raw(),n));
}


UpperFloat64 rec(LowerFloat64 const& x) {
    return UpperFloat64(rec_up(x.raw()));
}

UpperFloat64 sqrt(UpperFloat64 const& x) {
    return UpperFloat64(sqrt_up(x.raw()));
}

UpperFloat64 exp(UpperFloat64 const& x) {
    return UpperFloat64(exp_up(x.raw()));
}

UpperFloat64 log(UpperFloat64 const& x) {
    return UpperFloat64(log_up(x.raw()));
}

template<> Int integer_cast(UpperFloat64 const& x) { return static_cast<Int>(x._u.get_d()); }
template<> Nat integer_cast(UpperFloat64 const& x) { return static_cast<Nat>(x._u.get_d()); }




PositiveUpperFloat64 operator+(PositiveUpperFloat64 const& x1, PositiveUpperFloat64 const& x2) {
    return PositiveUpperFloat64(add_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat64 operator-(PositiveUpperFloat64 const& x1, LowerFloat64 const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=x2.raw());
    return PositiveUpperFloat64(sub_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat64 operator*(PositiveUpperFloat64 const& x1, PositiveUpperFloat64 const& x2) {
    return PositiveUpperFloat64(mul_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat64 operator/(PositiveUpperFloat64 const& x1, LowerFloat64 const& x2) {
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return PositiveUpperFloat64(div_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat64 pow(PositiveUpperFloat64 const& x, Nat m) {
    return PositiveUpperFloat64(pow_up(x.raw(),m));
}

PositiveUpperFloat64 half(PositiveUpperFloat64 const& x) {
    return PositiveUpperFloat64(half(x.raw()));
}



LowerFloat64 mul(LowerFloat64 const& x1, LowerFloat64 const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    return LowerFloat64(mul_down(x1.raw(),x2.raw()));
}

LowerFloat64 div(LowerFloat64 const& x1, UpperFloat64 const& x2) {
    //ARIADNE_PRECONDITION_MSG(x1.raw()>=0 && x2.raw()>=0,"x1="<<x1<<", x2="<<x2);
    return LowerFloat64(div_down(x1.raw(),x2.raw()));
}

LowerFloat64 rec(UpperFloat64 const& x) {
    return LowerFloat64(rec_down(x.raw()));
}

LowerFloat64 sqrt(LowerFloat64 const& x) {
    return LowerFloat64(sqrt_down(x.raw()));
}

LowerFloat64 exp(LowerFloat64 const& x) {
    return LowerFloat64(exp_down(x.raw()));
}

LowerFloat64 log(LowerFloat64 const& x) {
    return LowerFloat64(log_down(x.raw()));
}

OutputStream& operator<<(OutputStream& os, LowerFloat64 const& x) {
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64::set_rounding_downward();
    os << std::showpoint << std::setprecision(ValidatedFloat64::output_precision) << x.raw();
    Float64::set_rounding_mode(rnd);
    return os;
}

template<> Int integer_cast(LowerFloat64 const& x) { return static_cast<Int>(x._l.get_d()); }
template<> Nat integer_cast(LowerFloat64 const& x) { return static_cast<Nat>(x._l.get_d()); }


OutputStream& operator<<(OutputStream& os, ApproximateFloat64 const& x) {
    return os << std::showpoint << std::setprecision(ApproximateFloat64::output_precision) << x.raw();
}

template<> Int integer_cast(ApproximateFloat64 const& x) { return static_cast<Int>(x._a.get_d()); }
template<> Nat integer_cast(ApproximateFloat64 const& x) { return static_cast<Nat>(x._a.get_d()); }




ApproximateFloat64 make_float(Number<Approximate> x) { return ApproximateFloat64(x); }
LowerFloat64 make_float(Number<Lower> x) { return LowerFloat64(x); }
UpperFloat64 make_float(Number<Upper> x) { return UpperFloat64(x); }
ValidatedFloat64 make_float(Number<Validated> x) { return ValidatedFloat64(x); }
ValidatedFloat64 make_float(Number<Effective> x) { return ValidatedFloat64(x); }
ValidatedFloat64 make_float(Number<Exact> x) { return ValidatedFloat64(x); }
ValidatedFloat64 make_float(Real r) { return ValidatedFloat64(r); }
ValidatedFloat64 make_float(Rational q) { return ValidatedFloat64(q); }
ExactFloat64 make_float(Integer z) { return ExactFloat64(z); }

template class Float<Approximate,Precision64>;
template class Float<Lower,Precision64>;
template class Float<Upper,Precision64>;
template class Float<Bounded,Precision64>;
template class Float<Metric,Precision64>;
template class Float<Exact,Precision64>;

template Float<Lower,PrecisionMP>::Float(Float<Bounded,PrecisionMP>const&);
template Float<Upper,PrecisionMP>::Float(Float<Bounded,PrecisionMP>const&);

} // namespace Ariadne
