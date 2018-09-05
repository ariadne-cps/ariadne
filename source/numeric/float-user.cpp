/***************************************************************************
 *            float-user.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../config.hpp"
#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"

#include "../numeric/float-user.hpp"

#include "../numeric/integer.hpp"
#include "../numeric/dyadic.hpp"
#include "../numeric/decimal.hpp"
#include "../numeric/rational.hpp"
#include "../numeric/real.hpp"

#include "../numeric/number_wrapper.hpp"

namespace Ariadne {

template<class FE, class FLT, DisableIf<IsSame<FE,FLT>> =dummy> inline FE _make_error(FLT const& x) {
    typename FE::PrecisionType pre; return FE(Dyadic(x),upward,pre); }
template<class FE, class FLT, EnableIf<IsSame<FE,FLT>> =dummy> inline FE _make_error(FLT const& x) {
    return x; }

template<class FE, class FLT, class PRE> inline FE _make_error(FLT const& x, PRE pre) {
    return FE(Dyadic(x),upward,pre); }

FloatDP set(RoundUpward rnd, FloatMP const& x, DoublePrecision pr) { return FloatDP(Dyadic(x),rnd,pr); }
FloatDP set(RoundDownward rnd, FloatMP const& x, DoublePrecision pr) { return FloatDP(Dyadic(x),rnd,pr); }


template<class F> Nat Error<F>::output_places = 3;
template<class F> Nat Approximation<F>::output_places = 4;
template<class F> Nat Bounds<F>::output_places=8;
template<class F> Nat Value<F>::output_places = 16;

const FloatDPValue infty = FloatDPValue(FloatDP::inf(dp));

OutputStream& operator<<(OutputStream& os, Rounding const& rnd) {
    return os << ( rnd._rbp == ROUND_TO_NEAREST ? "near" : (rnd._rbp == ROUND_DOWNWARD ? "down" : "up") ); }

FloatError<DoublePrecision> operator"" _error(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatError<DoublePrecision>(FloatDP(x));
}

FloatValue<DoublePrecision> operator"" _exact(long double lx) {
    double x=lx;
    assert(x==lx);
    return FloatValue<DoublePrecision>(x);
}

FloatBall<DoublePrecision> operator"" _near(long double lx) {
    volatile double x=lx;
    volatile long double le=std::abs((long double)x-lx);
    volatile double e=le;
    while(e<le) { e*=(1+std::numeric_limits<double>::epsilon()); }
    return FloatBall<DoublePrecision>(x,e);
}

FloatUpperBound<DoublePrecision> operator"" _upper(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x<lx) { x+=min; }
    while (x<lx) { x+=std::abs(x)*eps; }
    return FloatUpperBound<DoublePrecision>(x);
}

FloatLowerBound<DoublePrecision> operator"" _lower(long double lx) {
    static const double eps = std::numeric_limits<double>::epsilon();
    static const double min = std::numeric_limits<double>::min();
    double x=lx;
    if(x>lx) { x-=min; }

    while (x>lx) { x-=std::abs(x)*eps; }
    return FloatLowerBound<DoublePrecision>(x);
}

FloatApproximation<DoublePrecision> operator"" _approx(long double lx) {
    double x=lx;
    return FloatApproximation<DoublePrecision>(x);
}



template<class F> Value<F>::Value(ExactDouble d, PR pr)
    : _v(d.get_d(),pr)
{
}

template<class F> Value<F>::Value(TwoExp const& t, PR pr)
    : _v(t,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==Dyadic(t),"Number 2^"<<t.exponent()<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Integer const& z, PR pr)
    : _v(z,to_nearest,pr)
{
    Rational q(_v);
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==z,"Integer "<<z<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Dyadic const& w, PR pr)
    : _v(w,pr)
{
    ARIADNE_ASSERT_MSG(Dyadic(this->_v)==w,"Dyadic number "<<w<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

template<class F> Value<F>::Value(Value<F> const& x, PR pr)
    : _v(x._v,to_nearest,pr)
{
    ARIADNE_ASSERT_MSG(*this==x,"Exact FloatValue "<<x<<" cannot be converted exactly to a floating-point number with precision "<<pr<<"; nearest is "<<(*this));
}

/*
template<class F> Value<F>::Value(Rational const& q, PR pr)
    : _v(0.0,pr)
{
    Bounds<F> x(q,pr);
    if(x.lower_raw()==x.upper_raw()) { this->_v==x.value_raw(); }
    else { ARIADNE_THROW(std::runtime_error,"FloatValue(Rational q, Precision pr)","q="<<q<<" cannot be expressed exactly to precision \n"); }
}
*/

template<class F> Value<F>::operator Dyadic() const {
    return this->_v.operator Dyadic();
}

template<class F> Value<F>::operator Rational() const {
    return Rational(this->operator Dyadic());
}

template<class F> Value<F>& Value<F>::operator=(Dyadic const& w) {
    _v=F(w,this->precision());
    ARIADNE_ASSERT_MSG(Dyadic(_v)==w,"Dyadic number "<<w<<" cannot be assigned exactly to a floating-point number with precision "<<this->precision()<<"; nearest is "<<(*this));
    return *this;
}

template<class F> auto Value<F>::create(ValidatedNumber const& y) const -> Ball<F> {
    return Ball<F>(y,this->precision());
}
/*
template<class F> Value<F>::operator ExactNumber() const {
    return ExactNumber(new NumberWrapper<Value<F>>(*this));
}
*/
template<class F, class FE> Ball<F,FE>::Ball(ExactDouble d, PR pr)
    : _v(d.get_d(),to_nearest,pr), _e(0,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(TwoExp t, PR pr)
    : _v(t,pr), _e(0u,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Integer const& z, PR pr) : Ball<F,FE>(Rational(z),pr) {
}

template<class F, class FE> Ball<F,FE>::Ball(Dyadic const& w, PR pr)
    : _v(F(w,to_nearest,pr)), _e(abs(Dyadic(_v)-w),upward,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Decimal const& d, PR pr)
    : Ball(Rational(d),pr) {
}

template<class F, class FE> Ball<F,FE>::Ball(Rational const& q, PR pr)
    : _v(F(q,to_nearest,pr)), _e(abs(Rational(_v)-q),upward,_error_precision<PRE>(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Real const& r, PR pr)
    : Ball(r.get(pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Ball<F,FE> const& x, PR pr)
    : _v(x._v,to_nearest,pr), _e(x._e,upward,_error_precision<PRE>(pr))
{
    F d = (this->_v>=x._v) ? sub(up,this->_v,x._v) : sub(up,x._v,this->_v);
    _e=add(up,_e,_make_error<FE>(d));
}

template<class F, class FE> Ball<F,FE>::Ball(ValidatedNumber const& y, PR pr)
    : Ball(y.get(MetricTag(),pr)) {
}

template<class F, class FE> Ball<F,FE>::Ball(Bounds<F> const& x, PRE pre)
    : _v(x.value_raw()) , _e(_make_error<FE>(x.error_raw(),pre)) {
}



template<class F, class FE> Ball<F,FE>::Ball(Real const& r, PR pr, PRE pre)
    : Ball(r.get(pr),pre) {
}

template<class F, class FE> Ball<F,FE>::Ball(Rational const& q, PR pr, PRE pre)
    : _v(F(q,to_nearest,pr)), _e(abs(Rational(_v)-q),upward,pre) {
}

template<class F, class FE> Ball<F,FE>::Ball(ValidatedNumber const& y, PR pr, PRE pre)
    : Ball(y.get(MetricTag(),pr,pre)) {
}

template<class F, class FE> Ball<F,FE> Ball<F,FE>::create(ValidatedNumber const& y) const {
    return Ball<F,FE>(y,this->precision());
}

template<class F, class FE> Ball<F,FE>& Ball<F,FE>::operator=(ValidatedNumber const& y)
{
    return *this = Ball(y,this->precision(),this->error_precision());
}

/*
template<class F, class FE> Ball<F,FE>::operator ValidatedNumber() const {
    return ValidatedNumber(new NumberWrapper<Ball<F,FE>>(*this));
}
*/

template<class F> Bounds<F>::Bounds(ExactDouble d, PR pr)
    : _l(d.get_d(),downward,pr),_u(d.get_d(),upward,pr) {
}

template<class F> Bounds<F>::Bounds(TwoExp t, PR pr)
    : _l(t,pr),_u(t,pr) {
}

template<class F> Bounds<F>::Bounds(Integer const& z, PR pr)
    : _l(z,downward,pr),_u(z,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Dyadic const& w, PR pr)
    : _l(w,downward,pr),_u(w,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Decimal const& d, PR pr)
    : Bounds(Rational(d),pr) {
}

template<class F> Bounds<F>::Bounds(Rational const& q, PR pr)
    : _l(q,downward,pr),_u(q,upward,pr) {
}

template<class F> Bounds<F>::Bounds(ExactDouble const& dl, ExactDouble const& du, PR pr)
    : _l(dl.get_d(),downward,pr),_u(du.get_d(),upward,pr) {
}

template<class F> Bounds<F>::Bounds(Dyadic const& wl, Dyadic const& wu, PR pr)
    : _l(wl,downward,pr),_u(wu,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Bounds<F> const& x, PR pr)
    : _l(x._l,downward,pr), _u(x._u,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Rational const& ql, Rational const& qu, PR pr)
    : _l(ql,downward,pr),_u(qu,upward,pr) {
}

template<class F> Bounds<F>::Bounds(Real const& x, PR pr)
    : Bounds(x.get(pr)) {
}

template<class F> Bounds<F>::Bounds(LowerBound<F> const& lower, ValidatedUpperNumber const& upper)
    : Bounds<F>(lower,lower.create(upper)) { }

template<class F> Bounds<F>::Bounds(ValidatedLowerNumber const& lower, UpperBound<F> const& upper)
    : Bounds<F>(upper.create(lower),upper) { }

template<class F> Bounds<F>::Bounds(ValidatedLowerNumber const& lower, ValidatedUpperNumber const& upper, PR pr)
    : Bounds<F>(lower.get(LowerTag(),pr),upper.get(UpperTag(),pr)) { }

template<class F> Bounds<F>::Bounds(ValidatedNumber const& y, PR pr)
    : Bounds(y.get(BoundedTag(),pr)) {
}

template<class F> Bounds<F> Bounds<F>::create(ValidatedNumber const& y) const {
    return Bounds<F>(y,this->precision());
}

template<class F> Bounds<F>& Bounds<F>::operator=(ValidatedNumber const& y) {
    return *this = Bounds<F>(y,this->precision());
}
/*
template<class F> Bounds<F>::operator ValidatedNumber() const {
    return ValidatedNumber(new NumberWrapper<Bounds<F>>(*this));
}
*/
template<class F> UpperBound<F>::UpperBound(ExactDouble d, PR pr)
    : _u(d.get_d(),upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(TwoExp t, PR pr)
    : _u(t,pr) {
}

template<class F> UpperBound<F>::UpperBound(Integer const& z, PR pr)
    : _u(z,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(Dyadic const& w, PR pr)
    : _u(w,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(Decimal const& d, PR pr)
    : UpperBound(Rational(d),pr) {
}

template<class F> UpperBound<F>::UpperBound(Rational const& q, PR pr)
    : _u(q,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(Real const& r, PR pr)
    : UpperBound(r.get(pr)) {
}

template<class F> UpperBound<F>::UpperBound(UpperBound<F> const& x, PR pr)
    : _u(x._u,upward,pr) {
}

template<class F> UpperBound<F>::UpperBound(ValidatedUpperNumber const& y, PR pr)
    : UpperBound(y.get(UpperTag(),pr)) {
}

template<class F> UpperBound<F>& UpperBound<F>::operator=(ValidatedUpperNumber const& y) {
    return *this = UpperBound<F>(y,this->precision());
}

template<class F> UpperBound<F>::operator ValidatedUpperNumber() const {
    ARIADNE_NOT_IMPLEMENTED;
    // return ValidatedUpperNumber(new NumberWrapper<UpperBound<F>>(*this));
}

template<class F> LowerBound<F> UpperBound<F>::create(ValidatedLowerNumber const& y) const {
    return LowerBound<F>(y,this->precision());
}

template<class F> UpperBound<F> UpperBound<F>::create(ValidatedUpperNumber const& y) const {
    return UpperBound<F>(y,this->precision());
}

template<class F> LowerBound<F>::LowerBound(ExactDouble d, PR pr)
    : _l(d.get_d(),pr) {
}

template<class F> LowerBound<F>::LowerBound(Integer const& z, PR pr)
    : _l(z,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(Dyadic const& w, PR pr)
    : _l(w,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(Decimal const& d, PR pr)
    : LowerBound(Rational(d),pr) {
}

template<class F> LowerBound<F>::LowerBound(Rational const& q, PR pr)
    : _l(q,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(Real const& r, PR pr)
    : LowerBound(r.get(pr)) {
}

template<class F> LowerBound<F>::LowerBound(LowerBound<F> const& x, PR pr)
    : _l(x._l,downward,pr) {
}

template<class F> LowerBound<F>::LowerBound(ValidatedLowerNumber const& y, PR pr)
    : LowerBound(y.get(LowerTag(),pr)) {
}

template<class F> LowerBound<F>& LowerBound<F>::operator=(ValidatedLowerNumber const& y) {
    return *this=LowerBound<F>(y,this->precision());
}

template<class F> LowerBound<F>::operator ValidatedLowerNumber() const {
    ARIADNE_NOT_IMPLEMENTED;
    //return ValidatedLowerNumber(new NumberWrapper<LowerBound<F>>(*this));
}

template<class F> LowerBound<F> LowerBound<F>::create(ValidatedLowerNumber const& y) const {
    return LowerBound<F>(y,this->precision());
}

template<class F> UpperBound<F> LowerBound<F>::create(ValidatedUpperNumber const& y) const {
    return UpperBound<F>(y,this->precision());
}

template<class F> Approximation<F>::Approximation(double d, PR pr)
    : _a(d,to_nearest,pr)
{
}

template<class F> Approximation<F>::Approximation(ExactDouble d, PR pr)
    : _a(d.get_d(),to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(TwoExp t, PR pr)
    : _a(t,pr) {
}

template<class F> Approximation<F>::Approximation(Integer const& z, PR pr)
    : _a(z,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Dyadic const& w, PR pr)
    : _a(w,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Decimal const& d, PR pr)
    : Approximation<F>(Rational(d),pr) {
}

template<class F> Approximation<F>::Approximation(Rational const& q, PR pr)
    : _a(q,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Approximation<F> const& x, PR pr)
    : _a(x._a,to_nearest,pr) {
}

template<class F> Approximation<F>::Approximation(Real const& r, PR pr)
    : Approximation<F>(r.get(pr)) {
}

template<class F> Approximation<F>::Approximation(ApproximateNumber const& y, PR pr)
    : Approximation<F>(y.get(ApproximateTag(),pr)) {
}

template<class F> Approximation<F>& Approximation<F>::operator=(ApproximateNumber const& y) {
    return *this=Approximation<F>(y,this->precision());
}

template<class F> Approximation<F> Approximation<F>::create(ApproximateNumber const& y) const {
    return Approximation<F>(y,this->precision());
}
/*
template<class F> Approximation<F>::operator ApproximateNumber() const {
    return ApproximateNumber(new NumberWrapper<Approximation<F>>(*this));
}
*/




template<class F> Approximation<F>::Approximation(LowerBound<F> const& x) : _a(x.raw()) {
}

template<class F> Approximation<F>::Approximation(UpperBound<F> const& x) : _a(x.raw()) {
}

template<class F> Approximation<F>::Approximation(Bounds<F> const& x) : _a(x.value_raw()) {
}

template<class F> Approximation<F>::Approximation(Value<F> const& x) : _a(x.raw()) {
}

template<class F> Approximation<F>::Approximation(Error<F> const& x) : _a(x.raw()) {
}

template<class F> LowerBound<F>::LowerBound(Bounds<F> const& x) : _l(x.lower_raw()) {
}

template<class F> LowerBound<F>::LowerBound(Ball<F> const& x) : _l(x.lower_raw()) {
}

template<class F> LowerBound<F>::LowerBound(Value<F> const& x) : _l(x.raw()) {
}

template<class F> UpperBound<F>::UpperBound(Bounds<F> const& x) : _u(x.upper_raw()) {
}

template<class F> UpperBound<F>::UpperBound(Ball<F> const& x) : _u(x.upper_raw()) {
}

template<class F> UpperBound<F>::UpperBound(Value<F> const& x) : _u(x.raw()) {
}

template<class F> UpperBound<F>::UpperBound(Error<F> const& x) : _u(x.raw()) {
}

template<class F> Bounds<F>::Bounds(Value<F> const& x) : _l(x.raw()), _u(x.raw()) {
}

template<class F, class FE> Ball<F,FE>::Ball(Bounds<F> const& x)
    : _v(x.value_raw()) , _e(_make_error<FE>(x.error_raw())) {
}

template<class F, class FE> Ball<F,FE>::Ball(Value<F> const& x) : _v(x.raw()), _e(_error_precision<PRE>(x.precision())) {
}


template<class F, class FE> Ball<F,FE> Ball<F,FE>::pm(Error<FE> const& e) const {
    return Ball<F,FE>(this->_v,add(up,this->_e,e._e));
}

template<class F> Bounds<F> Bounds<F>::pm(Error<F> e) const {
    return Bounds<F>(sub(down,this->_l,e._e),add(up,this->_u,e._e));
}



template<class F> struct Operations<Approximation<F>> {
    static Approximation<F> _floor(Approximation<F> const& x) {
        return Approximation<F>(floor(x._a)); }
    static Approximation<F> _ceil(Approximation<F> const& x) {
        return Approximation<F>(ceil(x._a)); }
    static Approximation<F> _round(Approximation<F> const& x) {
        return Approximation<F>(round(x._a)); }

    static Approximation<F> _abs(Approximation<F> const& x) {
        return Approximation<F>(abs(x._a)); }
    static Approximation<F> _max(Approximation<F> const& x, Approximation<F> const& y) {
        return Approximation<F>(max(x._a,y._a)); }
    static Approximation<F> _min(Approximation<F> const& x, Approximation<F> const& y) {
        return Approximation<F>(min(x._a,y._a)); }
    static PositiveApproximation<F> _mag(Approximation<F> const& x) {
        return PositiveApproximation<F>(abs(x._a)); }
    static PositiveApproximation<F> _mig(Approximation<F> const& x) {
        return PositiveApproximation<F>(abs(x._a)); }

    static Approximation<F> _nul(Approximation<F> const& x) {
        return Approximation<F>(nul(x._a)); }
    static Approximation<F> _pos(Approximation<F> const& x) {
        return Approximation<F>(pos(x._a)); }
    static Approximation<F> _neg(Approximation<F> const& x) {
        return Approximation<F>(neg(x._a)); }
    static Approximation<F> _hlf(Approximation<F> const& x) {
        return Approximation<F>(hlf(x._a)); }
    static Approximation<F> _sqr(Approximation<F> const& x) {
        return Approximation<F>(mul(near,x._a,x._a)); }
    static Approximation<F> _rec(Approximation<F> const& x) {
        return Approximation<F>(div(near,1.0,x._a)); }

    static Approximation<F> _add(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(add(near,x1._a,x2._a)); }
    static Approximation<F> _sub(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(sub(near,x1._a,x2._a)); }
    static Approximation<F> _mul(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(mul(near,x1._a,x2._a)); }
    static Approximation<F> _div(Approximation<F> const& x1, Approximation<F> const& x2) {
        return Approximation<F>(div(near,x1._a,x2._a)); }

    static Approximation<F> _pow(Approximation<F> const& x, Nat m) {
        return Approximation<F>(pow(approx,x._a,static_cast<Int>(m))); }
    static Approximation<F> _pow(Approximation<F> const& x, Int n) {
        return Approximation<F>(pow(approx,x._a,n)); }

    static Approximation<F> _sqrt(Approximation<F> const& x) {
        return Approximation<F>(sqrt(approx,x._a)); }
    static Approximation<F> _exp(Approximation<F> const& x) {
        return Approximation<F>(exp(approx,x._a)); }
    static Approximation<F> _log(Approximation<F> const& x) {
        return Approximation<F>(log(approx,x._a)); }
    static Approximation<F> _sin(Approximation<F> const& x) {
        return Approximation<F>(sin(approx,x._a)); }
    static Approximation<F> _cos(Approximation<F> const& x) {
        return Approximation<F>(cos(approx,x._a)); }
    static Approximation<F> _tan(Approximation<F> const& x) {
        return Approximation<F>(tan(approx,x._a)); }
    static Approximation<F> _asin(Approximation<F> const& x) {
        return Approximation<F>(asin(approx,x._a)); }
    static Approximation<F> _acos(Approximation<F> const& x) {
        return Approximation<F>(acos(approx,x._a)); }
    static Approximation<F> _atan(Approximation<F> const& x) {
        return Approximation<F>(atan(approx,x._a)); }

    static ApproximateKleenean _eq(Approximation<F> const& x1, Approximation<F> const& x2) {
        return x1._a==x2._a; }
    static ApproximateKleenean _lt(Approximation<F> const& x1, Approximation<F> const& x2) {
        return x1._a< x2._a; }

    static Bool _same(Approximation<F> const& x1, Approximation<F> const& x2) {
        return x1._a==x2._a; }

    static OutputStream& _write(OutputStream& os, Approximation<F> const& x) {
        return write(os,x.raw(),Approximation<F>::output_places,to_nearest);
    }

    static InputStream& _read(InputStream& is, Approximation<F>& x) {
        is >> x._a;
        return is;
    }
};


template<class F> struct Operations<LowerBound<F>> {

    static LowerBound<F> _max(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(max(x1._l,x2._l)); }
    static LowerBound<F> _min(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(min(x1._l,x2._l)); }
    static Approximation<F> _abs(LowerBound<F> const& x) {
        return abs(Approximation<F>(x)); }

    static LowerBound<F> _nul(LowerBound<F> const& x) {
        return LowerBound<F>(pos(x._l)); }
    static LowerBound<F> _pos(LowerBound<F> const& x) {
        return LowerBound<F>(pos(x._l)); }
    static UpperBound<F> _neg(LowerBound<F> const& x) {
        return UpperBound<F>(neg(x._l)); }
    static LowerBound<F> _hlf(LowerBound<F> const& x) {
        return LowerBound<F>(hlf(x._l)); }

    static LowerBound<F> _rec(UpperBound<F> const& x) {
        return LowerBound<F>(rec(down,x.raw())); }

    static LowerBound<F> _add(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(add(down,x1._l,x2._l)); }

    static Approximation<F> _sub(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return UpperBound<F>(sub(near,x1._l,x2._l)); }

    static LowerBound<F> _sub(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        return LowerBound<F>(sub(down,x1._l,x2._u)); }

    static LowerBound<F> _mul(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
        return LowerBound<F>(mul(down,x1.raw(),x2.raw())); }

    static LowerBound<F> _div(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        return LowerBound<F>(div(down,x1.raw(),x2.raw())); }

    static LowerBound<F> _pow(LowerBound<F> const& x, Nat m) {
        ARIADNE_PRECONDITION(x.raw()>=0);
        return LowerBound<F>(pow(down,x.raw(),static_cast<int>(m))); }

    static Approximation<F> _pow(LowerBound<F> const& x, Int n) {
        return pow(Approximation<F>(x),n); }

    static LowerBound<F> _sqrt(LowerBound<F> const& x) {
        return LowerBound<F>(sqrt(down,x.raw())); }

    static LowerBound<F> _exp(LowerBound<F> const& x) {
        return LowerBound<F>(exp(down,x.raw())); }

    static LowerBound<F> _log(LowerBound<F> const& x) {
        return LowerBound<F>(log(down,x.raw())); }

    static LowerBound<F> _atan(LowerBound<F> const& x) {
        return LowerBound<F>(atan(down,x.raw())); }

    static ValidatedNegatedSierpinskian _eq(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        if(x1._l>x2._u) { return false; }
        else { return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE); }
    }

    static ValidatedNegatedSierpinskian _lt(LowerBound<F> const& x1, UpperBound<F> const& x2) {
        if(x1._l>=x2._u) { return false; }
        else { return ValidatedNegatedSierpinskian(LogicalValue::LIKELY); }
    }

    static Bool _same(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return x1._l==x2._l;
    }

    static Bool _refines(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return x1._l>=x2._l;
    }

    static LowerBound<F> _refinement(LowerBound<F> const& x1, LowerBound<F> const& x2) {
        return LowerBound<F>(max(x1._l,x2._l));
    }


    static OutputStream& _write(OutputStream& os, LowerBound<F> const& x) {
        return write(os,x.raw(),Bounds<F>::output_places,downward);
    }

    static InputStream& _read(InputStream& is, LowerBound<F>& x) {
        ARIADNE_NOT_IMPLEMENTED;
    }
};

template<class F> struct Operations<UpperBound<F>> {
    static UpperBound<F> _max(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(max(x1._u,x2._u)); }

    static UpperBound<F> _min(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(min(x1._u,x2._u)); }

    static Approximation<F> _abs(UpperBound<F> const& x) {
        return abs(Approximation<F>(x)); }

    static UpperBound<F> _nul(UpperBound<F> const& x) {
        return UpperBound<F>(pos(x._u)); }

    static UpperBound<F> _pos(UpperBound<F> const& x) {
        return UpperBound<F>(pos(x._u)); }

    static LowerBound<F> _neg(UpperBound<F> const& x) {
        return LowerBound<F>(neg(x._u)); }

    static UpperBound<F> _hlf(UpperBound<F> const& x) {
        return UpperBound<F>(hlf(x._u)); }

    static UpperBound<F> _sqr(UpperBound<F> const& x) {
        ARIADNE_ASSERT(false); return UpperBound<F>(mul(up,x._u,x._u)); }

    static UpperBound<F> _rec(LowerBound<F> const& x) {
        return UpperBound<F>(rec(up,x.raw())); }

    static UpperBound<F> _add(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(add(up,x1._u,x2._u)); }

    static Approximation<F> _sub(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(sub(near,x1._u,x2._u)); }

    static UpperBound<F> _sub(UpperBound<F> const& x1, LowerBound<F> const& x2) {
        return UpperBound<F>(sub(up,x1._u,x2._l)); }

    static UpperBound<F> _mul(UpperBound<F> const& x1, UpperBound<F> const& x2) {
    //    ARIADNE_WARN("Multiplying FloatUpperBound "<<x1<<" with FloatUpperBound "<<x2<<" is unsafe");
        ARIADNE_PRECONDITION(x1.raw()>=0);
        ARIADNE_PRECONDITION(x2.raw()>=0);
        return UpperBound<F>(mul(up,x1._u,x2._u)); }

    static UpperBound<F> _div(UpperBound<F> const& x1, LowerBound<F> const& x2) {
    //    ARIADNE_WARN("Dividing FloatUpperBound "<<x1<<" by FloatLowerBound "<<x2<<" is unsafe");
        ARIADNE_PRECONDITION(x1.raw()>=0);
        ARIADNE_PRECONDITION(x2.raw()>=0);
        return UpperBound<F>(div(up,x1._u,x2._l)); }

    static UpperBound<F> _pow(UpperBound<F> const& x, Nat m) {
        ARIADNE_PRECONDITION(x.raw()>=0);
        return UpperBound<F>(pow(up,x._u,static_cast<Int>(m))); }

    static Approximation<F> _pow(UpperBound<F> const& x, Int n) {
        return pow(Approximation<F>(x),n); }

    static UpperBound<F> _sqrt(UpperBound<F> const& x) {
        return UpperBound<F>(sqrt(up,x.raw())); }

    static UpperBound<F> _exp(UpperBound<F> const& x) {
        return UpperBound<F>(exp(up,x.raw())); }

    static UpperBound<F> _log(UpperBound<F> const& x) {
        return UpperBound<F>(log(up,x.raw())); }

    static UpperBound<F> _atan(UpperBound<F> const& x) {
        return UpperBound<F>(atan(up,x.raw())); }

    static ValidatedNegatedSierpinskian _eq(UpperBound<F> const& x1, LowerBound<F> const& x2) {
        if(x1._u<x2._l) { return false; }
        else { return ValidatedNegatedSierpinskian(LogicalValue::INDETERMINATE); }
    }

    static ValidatedSierpinskian _lt(UpperBound<F> const& x1, LowerBound<F> const& x2) {
        if(x1._u< x2._l) { return true; }
        else { return ValidatedSierpinskian(LogicalValue::UNLIKELY); }
    }

    static Bool _same(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return x1._u==x2._u;
    }

    static Bool _refines(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return x1._u <= x2._u;
    }

    static UpperBound<F> _refinement(UpperBound<F> const& x1, UpperBound<F> const& x2) {
        return UpperBound<F>(min(x1._u,x2._u));
    }


    static Integer integer_cast(UpperBound<F> const& x) { return Integer(static_cast<int>(x._u.get_d())); }

    static OutputStream& _write(OutputStream& os, UpperBound<F> const& x) {
        return write(os,x.raw(),Bounds<F>::output_places,upward);
    }

    static InputStream& _read(InputStream& is, UpperBound<F>& x) {
        ARIADNE_NOT_IMPLEMENTED;
    }
};






template<class F> struct Operations<Bounds<F>> {
    typedef typename F::PrecisionType PR;

    static Bounds<F> _round(Bounds<F> const& x) {
        return Bounds<F>(round(x.lower_raw()),round(x.upper_raw()));
    }

    static Bounds<F> _max(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(max(x1.lower_raw(),x2.lower_raw()),max(x1.upper_raw(),x2.upper_raw()));
    }

    static Bounds<F> _min(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(min(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw()));
    }


    static Bounds<F> _abs(Bounds<F> const& x) {
        if(x.lower_raw()>=0) {
            return Bounds<F>(x.lower_raw(),x.upper_raw());
        } else if(x.upper_raw()<=0) {
            return Bounds<F>(neg(x.upper_raw()),neg(x.lower_raw()));
        } else {
            return Bounds<F>(F(0.0,x.precision()),max(neg(x.lower_raw()),x.upper_raw()));
        }
    }

    static PositiveLowerBound<F> _mig(Bounds<F> const& x) {
        return PositiveLowerBound<F>(max(0,max(x._l,neg(x._u))));
    }

    static PositiveUpperBound<F> _mag(Bounds<F> const& x) {
        return PositiveUpperBound<F>(max(neg(x._l),x._u));
    }

    static Bounds<F> _nul(Bounds<F> const& x) {
        return Bounds<F>(nul(x._l),nul(x._u));
    }

    static Bounds<F> _pos(Bounds<F> const& x) {
        return Bounds<F>(pos(x._l),pos(x._u));
    }

    static Bounds<F> _neg(Bounds<F> const& x) {
        return Bounds<F>(neg(x._u),neg(x._l));
    }

    static Bounds<F> _hlf(Bounds<F> const& x) {
        return Bounds<F>(hlf(x._l),hlf(x._u));
    }

    static Bounds<F> _sqr(Bounds<F> const& x) {
        const F& xl=x.lower_raw(); const F& xu=x.upper_raw();
        F rl,ru;
        if(xl>0.0) {
            rl=mul(down,xl,xl); ru=mul(up,xu,xu);
        } else if(xu<0.0) {
            rl=mul(down,xu,xu); ru=mul(up,xl,xl);
        } else {
            rl=nul(xl); ru=max(mul(up,xl,xl),mul(up,xu,xu));
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _rec(Bounds<F> const& x) {
        // IMPORTANT: Need to be careful when one of the bounds is 0, since if xl=-0.0 and xu>0, then 1/xl=-inf
        if(x._l>0 || x._u<0) {
            return Bounds<F>(rec(down,x._u),rec(up,x._l));
        } else {
            F inf_=F::inf(x.precision());
            //ARIADNE_THROW(DivideByZeroException,"FloatBounds rec(FloatBounds x)","x="<<x);
            return Bounds<F>(-inf_,+inf_);
        }
    }

    static Bounds<F> _add(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(add(down,x1._l,x2._l),add(up,x1._u,x2._u));
    }

    static Bounds<F> _sub(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(sub(down,x1._l,x2._u),sub(up,x1._u,x2._l));
    }

    static Bounds<F> _mul(Bounds<F> const& x1, Bounds<F> const& x2) {
        const F& x1l=x1._l; const F& x1u=x1._u;
        const F& x2l=x2._l; const F& x2u=x2._u;
        F rl,ru;
        if(x1l>=0) {
            if(x2l>=0) {
                rl=mul(down,x1l,x2l); ru=mul(up,x1u,x2u);
            } else if(x2u<=0) {
                rl=mul(down,x1u,x2l); ru=mul(up,x1l,x2u);
            } else {
                rl=mul(down,x1u,x2l); ru=mul(up,x1u,x2u);
            }
        }
        else if(x1u<=0) {
            if(x2l>=0) {
                rl=mul(down,x1l,x2u); ru=mul(up,x1u,x2l);
            } else if(x2u<=0) {
                rl=mul(down,x1u,x2u); ru=mul(up,x1l,x2l);
            } else {
                rl=mul(down,x1l,x2u); ru=mul(up,x1l,x2l);
            }
        } else {
            if(x2l>=0) {
                rl=mul(down,x1l,x2u); ru=mul(up,x1u,x2u);
            } else if(x2u<=0) {
                rl=mul(down,x1u,x2l); ru=mul(up,x1l,x2l);
            } else {
                rl=min(mul(down,x1u,x2l),mul(down,x1l,x2u));
                ru=max(mul(up,x1l,x2l),mul(up,x1u,x2u));
            }
        }
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _div(Bounds<F> const& x1, Bounds<F> const& x2) {
        const F& x1l=x1.lower_raw(); const F& x1u=x1.upper_raw();
        const F& x2l=x2.lower_raw(); const F& x2u=x2.upper_raw();
        F rl,ru;

        // IMPORTANT: Need to be careful when one of the bounds is 0, since if x2l=-0.0 and x1u>0, then x2l>=0 but x1u/x2l=-inf
        if(x2l>0) {
            if(x1l>=0) {
                rl=div(down,x1l,x2u); ru=div(up,x1u,x2l);
            } else if(x1u<=0) {
                rl=div(down,x1l,x2l); ru=div(up,x1u,x2u);
            } else {
                rl=div(down,x1l,x2l); ru=div(up,x1u,x2l);
            }
        }
        else if(x2u<0) {
            if(x1l>=0) {
                rl=div(down,x1u,x2u); ru=div(up,x1l,x2l);
            } else if(x1u<=0) {
                rl=div(down,x1u,x2l); ru=div(up,x1l,x2u);
            } else {
                rl=div(down,x1u,x2u); ru=div(up,x1l,x2u);
            }
        }
        else {
            //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
            PR pr=max(x1.precision(),x2.precision());
            rl=-F::inf(pr);
            ru=+F::inf(pr);
        }
        return Bounds<F>(rl,ru);
    }






    static Bounds<F> _pow(Bounds<F> const& x, Int n) {
        if(n<0) { return pow(rec(x),Nat(-n)); }
        else return pow(x,Nat(n));
    }

    static Bounds<F> _pow(Bounds<F> const& x, Nat m) {
        Bounds<F> y = x;
        if(m%2==0) { y=abs(x); }
        F rl=pow(down,y.lower_raw(),static_cast<Int>(m));
        F ru=pow(up,y.upper_raw(),static_cast<Int>(m));
        return Bounds<F>(rl,ru);
    }


    static Bounds<F> _sqrt(Bounds<F> const& x) {
        return Bounds<F>(sqrt(down,x.lower_raw()),sqrt(up,x.upper_raw()));
    }

    static Bounds<F> _exp(Bounds<F> const& x) {
        return Bounds<F>(exp(down,x.lower_raw()),exp(up,x.upper_raw()));
    }

    static Bounds<F> _log(Bounds<F> const& x) {
        return Bounds<F>(log(down,x.lower_raw()),log(up,x.upper_raw()));
    }


    static Bounds<F> _pi_val(PR pr) { return Bounds<F>(F::pi(down,pr),F::pi(up,pr)); }

    static Bounds<F> _sin(Bounds<F> const& x)
    {
        return cos(x-hlf(_pi_val(x.precision())));
    }

    static Bounds<F> _cos(Bounds<F> const& x)
    {
        ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
        typename F::RoundingModeType rnd = F::get_rounding_mode();
        PR prec=x.precision();

        const F one(1,prec);
        const Value<F> _two(2,prec);
        const Bounds<F> _pi=_pi_val(prec);
        const Bounds<F> two_pi=2*_pi_val(prec);
        if(x.error().raw()>two_pi.lower().raw()) { return Bounds<F>(-one,+one); }

        Value<F> n(round(div(near,x.value_raw(),(two_pi.value_raw()))));
        Bounds<F> y=x-_two*(n*_pi);

        ARIADNE_ASSERT(y.lower_raw()<=_pi.upper_raw());
        ARIADNE_ASSERT(y.upper_raw()>=-_pi.upper_raw());

        F rl,ru;
        if(y.lower_raw()<=-_pi.lower_raw()) {
            if(y.upper_raw()<=0.0) { rl=-one; ru=cos(up,y.upper_raw()); }
            else { rl=-one; ru=+one; }
        } else if(y.lower_raw()<=0.0) {
            if(y.upper_raw()<=0.0) { rl=cos(down,y.lower_raw()); ru=cos(up,y.upper_raw()); }
            else if(y.upper_raw()<=_pi.lower_raw()) { rl=cos(down,max(-y.lower_raw(),y.upper_raw())); ru=+one; }
            else { rl=-one; ru=+one; }
        } else if(y.lower_raw()<=_pi.upper_raw()) {
            if(y.upper_raw()<=_pi.lower_raw()) { rl=cos(down,y.upper_raw()); ru=cos(up,y.lower_raw()); }
            else if(y.upper_raw()<=two_pi.lower_raw()) { rl=-one; ru=cos(up,min(y.lower_raw(),sub(down,two_pi.lower_raw(),y.upper_raw()))); }
            else { rl=-one; ru=+one; }
        } else {
            assert(false);
        }

        F::set_rounding_mode(rnd);
        return Bounds<F>(rl,ru);
    }

    static Bounds<F> _tan(Bounds<F> const& x) {
        return mul(sin(x),rec(cos(x)));
    }

    static Bounds<F> _asin(Bounds<F> const& x) {
        ARIADNE_NOT_IMPLEMENTED;
    }

    static Bounds<F> _acos(Bounds<F> const& x) {
        ARIADNE_NOT_IMPLEMENTED;
    }

    static Bounds<F> _atan(Bounds<F> const& x) {
        return Bounds<F>(atan(down,x._l),atan(up,x._u));
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static LogicalType<ValidatedTag> _eq(Bounds<F> const& x1, Bounds<F> const& x2) {
        if(x1.upper_raw()<x2.lower_raw() || x1.lower_raw()>x2.upper_raw()) { return false; }
        else if(x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw()) { return true; }
        else { return indeterminate; }
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static LogicalType<ValidatedTag> _lt(Bounds<F> const& x1, Bounds<F> const& x2) {
        if(x1.upper_raw()< x2.lower_raw()) { return true; }
        else if(x1.lower_raw()>=x2.upper_raw()) { return false; }
        else { return indeterminate; }
    }


    static Bounds<F> _widen(Bounds<F> const& x)
    {
        const F& xl=x.lower_raw();
        const F& xu=x.upper_raw();
        const F m=std::numeric_limits<float>::min();
        F wl=sub(down,xl,m);
        F wu=add(up,xu,m);
        assert(wl<xl); assert(wu>xu);
        return Bounds<F>(wl,wu);
    }

    static Bounds<F> _narrow(Bounds<F> const& x)
    {
        const F& xl=x.lower_raw();
        const F& xu=x.upper_raw();
        const F m=std::numeric_limits<float>::min();
        F nu=add(down,xu,m);
        F nl=add(up,xl,m);
        assert(xl<nl); assert(nu<xu);
        return Bounds<F>(nl,nu);
    }

    static Bounds<F> _trunc(Bounds<F> const& x)
    {
        typename F::RoundingModeType rm=F::get_rounding_mode();
        const double& xl=x.lower_raw().get_d();
        const double& xu=x.upper_raw().get_d();
        // Use machine epsilon instead of minimum to move away from zero
        const float fm=std::numeric_limits<float>::epsilon();
        volatile float tu=xu;
        if(tu<xu) { F::set_rounding_upward(); tu+=fm; }
        volatile float tl=xl;
        if(tl>xl) { F::set_rounding_downward(); tl-=fm; }
        F::set_rounding_mode(rm);
        assert(tl<=xl); assert(tu>=xu);
        return Bounds<F>(double(tl),double(tu));
    }

    static Bounds<F> _trunc(Bounds<F> const& x, Nat n)
    {
        Bounds<F> _e=Bounds<F>(std::pow(2.0,52-(Int)n));
        Bounds<F> y=x+_e;
        return y-_e;
    }

    static Integer integer_cast(Bounds<F> const& x) {
        return Integer(static_cast<int>(x.value_raw().get_d()));
    }

    static auto is_zero(Bounds<F> const& x) -> LogicalType<ValidatedTag> {
        if(x.lower_raw()>0.0 || x.upper_raw()<0.0) { return false; }
        else if(x.lower_raw()==0.0 && x.upper_raw()==0.0) { return true; }
        else { return indeterminate; }
    }

    static auto is_positive(Bounds<F> const& x) -> LogicalType<ValidatedTag> {
        if(x.lower_raw()>=0.0) { return true; }
        else if(x.upper_raw()<0.0) { return false; }
        else { return indeterminate; }
    }

    static Bool _same(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l==x2._l && x1._u==x2._u; }

    static Bool _models(Bounds<F> const& x1, Value<F> const& x2) {
        return x1._l<=x2._v && x1._u >= x2._v; }

    static Bool _consistent(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l<=x2._u && x1._u >= x2._l; }

    static Bool _inconsistent(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l>x2._u || x1._u < x2._l; }

    static Bool _refines(Bounds<F> const& x1, Bounds<F> const& x2) {
        return x1._l>=x2._l && x1._u <= x2._u; }

    static Bounds<F> _refinement(Bounds<F> const& x1, Bounds<F> const& x2) {
        return Bounds<F>(max(x1._l,x2._l),min(x1._u,x2._u)); }



    static OutputStream& _write(OutputStream& os, const Bounds<F>& x) {
        typename F::RoundingModeType rnd=F::get_rounding_mode();
        os << '{';
        write(os,x.lower().raw(),Bounds<F>::output_places,downward);
        os << ':';
        write(os,x.upper().raw(),Bounds<F>::output_places,upward);
        os << '}';
        return os;

    }

    static InputStream& _read(InputStream& is, Bounds<F>& x) {
        char cl,cm,cr;
        F _l,_u;
        auto rnd=F::get_rounding_mode();
        is >> cl;
        F::set_rounding_downward();
        is >> _l;
        is >> cm;
        F::set_rounding_upward();
        is >> _u;
        is >> cr;
        F::set_rounding_mode(rnd);
        ARIADNE_ASSERT(not is.fail());
        ARIADNE_ASSERT(cl=='[' || cl=='(');
        ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
        ARIADNE_ASSERT(cr==']' || cr==')');
        x._l=_l; x._u=_u;
        return is;
    }
};

template<class F> auto is_positive(Bounds<F> const&) -> LogicalType<ValidatedTag>;

inline int log10floor(double const& x) { return std::max(std::floor(std::log10(x)),-65280.); }
inline int log10floor(FloatMP const& x) { return log10floor(x.get_d()); }
inline int abslog10floor(double const& x) { return log10floor(std::abs(x)); }


template<> OutputStream& Operations<FloatBounds<MultiplePrecision>>::_write(OutputStream& os, const FloatBounds<MultiplePrecision>& x)
{
    static const double log2ten = 3.3219280948873621817;
    using std::max; using std::min;
    FloatMP const& l=x.lower_raw();
    FloatMP const& u=x.upper_raw();
    double ldbl=l.get_d();
    double udbl=u.get_d();
    if(ldbl==0.0 && udbl==0.0) { return os << "0.0[:]"; }
    int errplc=static_cast<int>(FloatError<MultiplePrecision>::output_places);
    //int bndplc=FloatBounds<MultiplePrecision>::output_places;
    int precplc=x.precision()/log2ten;
    int log10wdth=log10floor(sub(to_nearest,u,l));
    int log10mag=log10floor(max(-ldbl,udbl));
    int dgtswdth=errplc-(log10wdth+1); // Digits appropriate given width of interval
    //int dgtsbnd=bndplc-(log10mag+1); // Digits appropriate given asked-for precision of bounded objects
    int dgtsprec=precplc-(log10mag+1); // Digits appropriate given precision of objects
    Nat dgts=static_cast<Nat>(max(min(dgtswdth,dgtsprec),1));
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

template<> OutputStream& Operations<FloatBounds<DoublePrecision>>::_write(OutputStream& os, const FloatBounds<DoublePrecision>& x)
{
    MultiplePrecision prec(64);
    return os << FloatBounds<MultiplePrecision>(FloatMP(x.lower_raw(),prec),FloatMP(x.upper_raw(),prec));
}

template<class F, class FE> struct Operations<Ball<F,FE>> {

    static Ball<F,FE> _nul(Ball<F,FE> const& x) {
        return Ball<F,FE>(nul(x._v),nul(x._e));
    }

    static Ball<F,FE> _pos(Ball<F,FE> const& x) {
        return Ball<F,FE>(pos(x._v),x._e);
    }

    static Ball<F,FE> _neg(Ball<F,FE> const& x) {
        return Ball<F,FE>(neg(x._v),x._e);
    }

    static Ball<F,FE> _hlf(Ball<F,FE> const& x) {
        return Ball<F,FE>(hlf(x._v),hlf(x._e));
    }

    static Ball<F,FE> _sqr(Ball<F,FE> const& x) {
        Ball<F,FE> r=x*x;
        if(r._e>r._v) {
            r._e=hlf(add(up,r._e,_make_error<FE>(r._v)));
            r._v=F(Dyadic(r._e),upward,x.precision());
        }
        return r;
    }

    static Ball<F,FE> _rec(Ball<F,FE> const& x) {
        // Use this code to find value same as reciprocal value
        auto rv=rec(approx,x._v);
        auto ru=rec(up,sub(down,x._v,x._e));
        auto rl=rec(down,add(up,x._v,x._e));
        auto re=max(sub(up,ru,rv),sub(up,rv,rl));
        return Ball<F,FE>(rv,_make_error<FE>(re));
    #ifdef ARIADNE_UNDEFINED
        // Use this code to get same result as interval computation
        auto ru=rec(up,sub(down,x._v,x._e));
        auto rl=rec(down,add(up,x._v,x._e));
        auto re=hlf(sub(up,ru,rl));
        auto rv=hlf(add(near,rl,ru));
        return Ball<F,FE>(rv,re);
    #endif
    }

    static Ball<F,FE> _add(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        auto rv=add(near,x._v,y._v);
        auto ru=add(up,x._v,y._v);
        auto rl=add(down,x._v,y._v);
        auto ae=_make_error<FE>(hlf(sub(up,ru,rl)));
        auto re=add(up,ae,add(up,x._e,y._e));
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _sub(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        auto rv=sub(near,x._v,y._v);
        auto ru=sub(up,x._v,y._v);
        auto rl=sub(down,x._v,y._v);
        auto ae=_make_error<FE>(hlf(sub(up,ru,rl)));
        auto re=add(up,ae,add(up,x._e,y._e));
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _mul(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        auto rv=mul(near,x._v,y._v);
        auto ru=mul(up,x._v,y._v);
        auto rl=mul(down,x._v,y._v);
        auto re0=_make_error<FE>(hlf(sub(up,ru,rl)));
        auto re1=add(up,re0,mul(up,x._e,y._e));
        auto re2=add(up,mul(up,_make_error<FE>(abs(x._v)),y._e),mul(up,x._e,_make_error<FE>(abs(y._v))));
        auto re=add(up,re1,re2);
        return Ball<F,FE>(rv,re);
    }

    static Ball<F,FE> _div(Ball<F,FE> const& x, Ball<F,FE> const& y) {
        return x*rec(y);
    }

    static Ball<F,FE> _pow(Ball<F,FE> const& x, Nat m) {
        return Ball<F,FE>(pow(Bounds<F>(x),m));
    }

    static Ball<F,FE> _pow(Ball<F,FE> const& x, Int n) {
        return Ball<F,FE>(pow(Bounds<F>(x),n));
    }

    static Ball<F,FE> _sqrt(Ball<F,FE> const& x) {
        return Ball<F,FE>(sqrt(Bounds<F>(x)));
    }

    static Ball<F,FE> _exp(Ball<F,FE> const& x) {
        return Ball<F,FE>(exp(Bounds<F>(x)));
    }

    static Ball<F,FE> _log(Ball<F,FE> const& x) {
        return Ball<F,FE>(log(Bounds<F>(x)));
    }

    static Ball<F,FE> _sin(Ball<F,FE> const& x) {
        return Ball<F,FE>(sin(Bounds<F>(x)));
    }

    static Ball<F,FE> _cos(Ball<F,FE> const& x) {
        return Ball<F,FE>(cos(Bounds<F>(x)));
    }

    static Ball<F,FE> _tan(Ball<F,FE> const& x) {
        return Ball<F,FE>(tan(Bounds<F>(x)));
    }

    static Ball<F,FE> _asin(Ball<F,FE> const& x) {
        return Ball<F,FE>(asin(Bounds<F>(x)));
    }

    static Ball<F,FE> _acos(Ball<F,FE> const& x) {
        return Ball<F,FE>(acos(Bounds<F>(x)));
    }

    static Ball<F,FE> _atan(Ball<F,FE> const& x) {
        return Ball<F,FE>(atan(Bounds<F>(x)));
    }


    static Ball<F,FE> _abs(Ball<F,FE> const& x) {
        if(x._e<abs(x._v)) { return x; }
        else { auto rv=hlf(add(up,abs(x._v),x._e)); return Ball<F,FE>(rv,_make_error<FE>(rv)); }
    }

    static Ball<F,FE> _max(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return hlf((x1+x2)+abs(x1-x2));
    }

    static Ball<F,FE> _min(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return hlf((x1+x2)-abs(x1-x2));
    }

    static Error<F> _mag(Ball<F,FE> const& x) {
        return PositiveUpperBound<F>(add(up,abs(x._v),x._e));
    }

    static PositiveLowerBound<F> _mig(Ball<F,FE> const& x) {
        return PositiveLowerBound<F>(max(0,sub(down,abs(x._v),x._e)));
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static ValidatedKleenean _eq(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Bounds<F>(x1) == Bounds<F>(x2);
    }

    //! \related Bounds<F> \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    static ValidatedKleenean _lt(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Bounds<F>(x1) <  Bounds<F>(x2);
    }

    static Bool _same(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return x1._v==x2._v && x1._e==x2._e;
    }

    static Bool _models(Ball<F,FE> const& x1, Value<F> const& x2) {
        return (x1._v>=x2._v ? sub(up,x1._v,x2._v) : sub(up,x2._v,x1._v)) <= x1._e;
    }

    static Bool _consistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return consistent(Bounds<F>(x1),Bounds<F>(x2));
    }

    static Bool _inconsistent(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return inconsistent(Bounds<F>(x1),Bounds<F>(x2));
    }

    static Bool _refines(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return (x1._v>=x2._v ? sub(up,x1._v,x2._v) : sub(up,x2._v,x1._v)) <= sub(down,x2._e, x1._e);
    }

    static Ball<F,FE> _refinement(Ball<F,FE> const& x1, Ball<F,FE> const& x2) {
        return Ball<F,FE>(refinement(Bounds<F>(x1),Bounds<F>(x2)));
    }

    static OutputStream& _write(OutputStream& os, Ball<F,FE> const& x) {
        return os << x.value() << "\u00b1" << x.error();
    }

    static InputStream& _read(InputStream& is, Ball<F,FE>& x) {
        static const char pmstr[] = "\u00b1";
        char cpm[3];
        F _v; FE _e;
        auto rnd=F::get_rounding_mode();
        F::set_rounding_to_nearest();
        is >> _v;
        is >> cpm[0] >> cpm[1];
        F::set_rounding_mode(rnd);
        auto rnde=FE::get_rounding_mode();
        FE::set_rounding_upward();
        is >> _e;
        FE::set_rounding_mode(rnde);
        ARIADNE_ASSERT(not is.fail());
        ARIADNE_ASSERT(std::strcmp(cpm,pmstr));
        x._v=_v; x._e=_e;
        return is;
    }

};

template<> OutputStream& Operations<FloatBall<MultiplePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision> const& x) {
    // Write based on number of correct digits
    static const double log2ten = 3.3219280948873621817;
    static const char pmstr[] = "\u00b1";
    static const char hlfstr[] = "\u00bd";
    FloatMP const& v=x.value_raw();
    FloatMP const& e=x.error_raw();
    double edbl=e.get_d();
    // Compute the number of decimal places to be displayed
    Nat errplc = static_cast<Nat>(FloatError<MultiplePrecision>::output_places);
    Nat log10err = static_cast<Nat>(log10floor(edbl));
    Nat dgtserr = errplc-(log10err+1);
    Nat dgtsval = (x.value().raw()==0) ? dgtserr : std::floor((x.value().precision()+1-x.value().raw().exponent())/log2ten);
    Nat dgts = std::max(std::min(dgtsval,dgtserr),errplc);
    if(edbl==0.0) { dgts = dgtsval; }
    DecimalPlaces plcs{dgts};
    // Get string version of mpfr values
    String vstr=print(v,plcs,MPFR_RNDN);
    String estr=print(e,plcs,MPFR_RNDU);

    // Find position of first significan digit of error
    auto vcstr=vstr.c_str(); auto ecstr=estr.c_str();
    size_t cpl=0;
    if(edbl==0.0) {
        cpl=std::strlen(vcstr);
    } else if(edbl<1.0) {
        const char* vptr = std::strchr(vcstr,'.');
        const char* eptr = std::strchr(ecstr,'.');
        ++vptr; ++eptr;
        while((*eptr)=='0') { ++eptr; ++vptr; }
        cpl = static_cast<size_t>(vptr-vcstr);
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

template<> OutputStream& Operations<FloatBall<MultiplePrecision,DoublePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision,DoublePrecision> const& x) {
    MultiplePrecision prec(64);
    return os << FloatBall<MultiplePrecision,MultiplePrecision>(x.value_raw(),FloatMP(x.error_raw(),prec));
}

template<> OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, FloatBall<DoublePrecision> const& x) {
    MultiplePrecision prec(64);
    return os << FloatBall<MultiplePrecision>(FloatMP(x.value_raw(),prec),FloatMP(x.error_raw(),prec));
}





// Mixed BoundedTag - ExactTag operations
template<class F> Bounds<F> _add(Bounds<F> const& x1, Value<F> const& x2) {
    return Bounds<F>(add(down,x1._l,x2._v),add(up,x1._u,x2._v));
}

template<class F> Bounds<F> _add(Value<F> const& x1, Bounds<F> const& x2) {
    return Bounds<F>(add(down,x1._v,x2._l),add(down,x1._v,x2._u));
}

template<class F> Bounds<F> _sub(Bounds<F> const& x1, Value<F> const& x2) {
    return Bounds<F>(sub(down,x1._l,x2._v),sub(up,x1._u,x2._v));
}

template<class F> Bounds<F> _sub(Value<F> const& x1, Bounds<F> const& x2) {
    return Bounds<F>(sub(down,x1._v,x2._u),sub(up,x1._v,x2._l));
}

template<class F> Bounds<F> _mul(Bounds<F> const& x1, Value<F> const& x2) {
    const F& x1l=x1.lower_raw(); const F& x1u=x1.upper_raw();
    const F& x2v=x2.raw();
    F rl,ru;
    if(x2v>=0.0) {
        rl=mul(down,x1l,x2v); ru=mul(up,x1u,x2v);
    } else {
        rl=mul(down,x1u,x2v); ru=mul(up,x1l,x2v);
    }
    return Bounds<F>(rl,ru);
}


template<class F> Bounds<F> _mul(Value<F> const& x1, Bounds<F> const& x2) {
    const F& x1v=x1.raw();
    const F& x2l=x2.lower_raw(); const F& x2u=x2.upper_raw();
    F rl,ru;
    if(x1v>=0.0) {
        rl=mul(down,x1v,x2l); ru=mul(up,x1v,x2u);
    } else {
        rl=mul(down,x1v,x2u); ru=mul(up,x1v,x2l);
    }
    return Bounds<F>(rl,ru);
}

template<class F> Bounds<F> _div(Bounds<F> const& x1, Value<F> const& x2)
{
    const F& x1l=x1.lower_raw();
    const F& x1u=x1.upper_raw();
    const F& x2v=x2.raw();
    F rl,ru;
    if(x2v>0) {
        rl=div(down,x1l,x2v); ru=div(up,x1u,x2v);
    } else if(x2v<0) {
        rl=div(down,x1u,x2v); ru=div(up,x1l,x2v);
    } else {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatBounds const& x1, FloatValue x2)","x1="<<x1<<", x2="<<x2);
        auto pr=min(x1.precision(),x2.precision());
        rl=-F::inf(pr);
        ru=+F::inf(pr);
    }
    return Bounds<F>(rl,ru);
}


template<class F> Bounds<F> _div(Value<F> const& x1, Bounds<F> const& x2)
{
    const F& x1v=x1.raw();
    const F& i2l=x2.lower_raw();
    const F& i2u=x2.upper_raw();
    F rl,ru;
    if(i2l<=0 && i2u>=0) {
        //ARIADNE_THROW(DivideByZeroException,"FloatBounds div(FloatValue const& x1, FloatBounds x2)","x1="<<x1<<", x2="<<x2);
        auto pr=min(x1.precision(),x2.precision());
        rl=-F::inf(pr);
        ru=+F::inf(pr);
    } else if(x1v>=0) {
        rl=div(down,x1v,i2u); ru=div(up,x1v,i2l);
    } else {
        rl=div(down,x1v,i2l); ru=div(up,x1v,i2u);
    }
    return Bounds<F>(rl,ru);
}




template<class F> struct Operations<Value<F>> {
    static Value<F> _max(Value<F> const& x1,  Value<F> const& x2) {
        return Value<F>(max(x1._v,x2._v)); }

    static Value<F> _min(Value<F> const& x1,  Value<F> const& x2) {
        return Value<F>(min(x1._v,x2._v)); }

    static Value<F> _abs(Value<F> const& x) {
        return Value<F>(abs(x._v)); }

    static LowerBound<F> _mig(Value<F> const& x) {
        return LowerBound<F>(abs(x._v)); }

    static Error<F> _mag(Value<F> const& x) {
        return Error<F>(abs(x._v)); }


    static Value<F> _nul(Value<F> const& x) {
        return Value<F>(nul(x._v)); }

    static Value<F> _pos(Value<F> const& x) {
        return Value<F>(pos(x._v)); }

    static Value<F> _neg(Value<F> const& x) {
        return Value<F>(neg(x._v)); }

    static Value<F> _hlf(Value<F> const& x) {
        return Value<F>(hlf(x._v)); }

    static Bounds<F> _sqr(Value<F> const& x) {
        return Bounds<F>(mul(down,x._v,x._v),mul(up,x._v,x._v)); }

    static Bounds<F> _rec(Value<F> const& x) {
        return Bounds<F>(rec(down,x._v),rec(up,x._v)); }

    static Bounds<F> _add(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(add(down,x1._v,x2._v),add(up,x1._v,x2._v)); }

    static Bounds<F> _sub(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(sub(down,x1._v,x2._v),sub(up,x1._v,x2._v)); }

    static Bounds<F> _mul(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(mul(down,x1._v,x2._v),mul(up,x1._v,x2._v)); }

    static Bounds<F> _div(Value<F> const& x1, Value<F> const& x2) {
        return Bounds<F>(div(down,x1._v,x2._v),div(up,x1._v,x2._v)); }

    static Value<F> _mul(Value<F> const& x, TwoExp y) {
        Value<F> yv(y,x.precision()); return Value<F>(mul(near,x.raw(),yv.raw())); }

    static Value<F> _div(Value<F> const& x, TwoExp y) {
        Value<F> yv(y,x.precision()); return Value<F>(div(near,x.raw(),yv.raw())); }

    static Bounds<F> _pow(Value<F> const& x, Nat m) {
        return pow(Bounds<F>(x),m); }

    static Bounds<F> _pow(Value<F> const& x, Int n) {
        return pow(Bounds<F>(x),n); }

    static Bounds<F> _med(Value<F> const& x1, Value<F> const& x2) {
        return add(hlf(x1),hlf(x2)); }

    static Bounds<F> _rad(Value<F> const& x1, Value<F> const& x2) {
        return sub(hlf(x2),hlf(x1)); }

    static Bounds<F> _sqrt(Value<F> const& x) {
        return sqrt(Bounds<F>(x)); }

    static Bounds<F> _exp(Value<F> const& x) {
        return exp(Bounds<F>(x)); }

    static Bounds<F> _log(Value<F> const& x) {
        return log(Bounds<F>(x)); }

    static Bounds<F> _sin(Value<F> const& x) {
        return sin(Bounds<F>(x)); }

    static Bounds<F> _cos(Value<F> const& x) {
        return cos(Bounds<F>(x)); }

    static Bounds<F> _tan(Value<F> const& x) {
        return tan(Bounds<F>(x)); }

    static Bounds<F> _asin(Value<F> const& x) {
        return asin(Bounds<F>(x)); }

    static Bounds<F> _acos(Value<F> const& x) {
        return acos(Bounds<F>(x)); }

    static Bounds<F> _atan(Value<F> const& x) {
        return atan(Bounds<F>(x)); }

    static Boolean _eq(Value<F> const& x1, Value<F> const& x2) {
        return x1._v == x2._v; }

    static Boolean _lt(Value<F> const& x1, Value<F> const& x2) {
        return x1._v <  x2._v; }

    static Bool _same(Value<F> const& x1, Value<F> const& x2) {
        return x1._v==x2._v; }

    static OutputStream& _write(OutputStream& os, Value<F> const& x) {
        return write(os,x.raw(),Value<F>::output_places,to_nearest);
    }

    static InputStream& _read(InputStream& is, Value<F>& x) {
        ARIADNE_NOT_IMPLEMENTED;
        auto v = nul(x._v);
        is >> v;
        ARIADNE_ASSERT(not is.fail());
        x._v=v;
        return is;
    }

    static Integer integer_cast(Value<F> const& x) {
        Dyadic w(x);
        Integer z=round(w);
        ARIADNE_ASSERT(z==w);
        return z;
    }

};

Rational cast_exact(Real const& x) {
    return Rational(cast_exact(FloatApproximation<DoublePrecision>(x,dp)));
}


template<class F> struct Operations<PositiveApproximation<F>> {
    static PositiveApproximation<F> _nul(PositiveApproximation<F> const& x) {
        return PositiveApproximation<F>(nul(x._a)); }
    static PositiveApproximation<F> _sqr(PositiveApproximation<F> const& x) {
        return PositiveApproximation<F>(mul(near,x._a,x._a)); }
    static PositiveApproximation<F> _rec(PositiveApproximation<F> const& x) {
        return PositiveApproximation<F>(rec(near,x._a)); }
    static PositiveApproximation<F> _add(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(add(near,x1._a,x2._a)); }
    static PositiveApproximation<F> _mul(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(mul(near,x1._a,x2._a)); }
    static PositiveApproximation<F> _div(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(div(near,x1._a,x2._a)); }
    static PositiveApproximation<F> _pow(PositiveApproximation<F> const& x1, Int n2) {
        return PositiveApproximation<F>(pow(approx,x1._a,n2)); }
    static Approximation<F> _log(PositiveApproximation<F> const& x) {
        return Approximation<F>(log(approx,x._a)); }
    static PositiveApproximation<F> _max(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(max(x1._a,x2._a)); }
    static PositiveApproximation<F> _min(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return PositiveApproximation<F>(min(x1._a,x2._a)); }
    static PositiveApproximation<F> _abs(PositiveApproximation<F> const& x) {
        return PositiveApproximation<F>(x._a); }
    static Bool _same(PositiveApproximation<F> const& x1, PositiveApproximation<F> const& x2) {
        return x1._a == x2._a; }
    static OutputStream& _write(OutputStream& os, PositiveApproximation<F> const& x) {
        return write(os,x.raw(),Approximation<F>::output_places,upward); }
    static InputStream& _read(InputStream& is, PositiveApproximation<F>& x) {
        Approximation<F> xa; is >> xa; x=PositiveApproximation<F>(xa); return is; }
};

template<class F> struct Operations<PositiveUpperBound<F>> {
    static PositiveUpperBound<F> _nul(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(nul(x._u));
    }

    static PositiveUpperBound<F> _hlf(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(hlf(x._u));
    }

    static PositiveUpperBound<F> _sqr(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(mul(up,x._u,x._u));
    }

    static PositiveLowerBound<F> _rec(PositiveUpperBound<F> const& x) {
        return PositiveLowerBound<F>(rec(down,x._u));
    }

    static PositiveUpperBound<F> _rec(PositiveLowerBound<F> const& x) {
        ARIADNE_ASSERT_MSG(x._l>=0.0,x); return PositiveUpperBound<F>(rec(up,x._l));
    }

    static PositiveUpperBound<F> _add(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(add(up,x1._u,x2._u));
    }

    static PositiveUpperBound<F> _mul(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(mul(up,x1._u,x2._u));
    }

    static PositiveUpperBound<F> _div(PositiveUpperBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveUpperBound<F>(div(up,x1._u,x2._l));
    }

    static PositiveUpperBound<F> _div(PositiveUpperBound<F> const& x1, LowerBound<F> const& x2) {
        ARIADNE_ASSERT_MSG(x2._l>=0.0,x2); return PositiveUpperBound<F>(div(up,x1._u,x2._l));
    }

    static PositiveUpperBound<F> _pow(PositiveUpperBound<F> const& x1, Nat m2) {
        return PositiveUpperBound<F>(pow(up,x1._u,static_cast<Int>(m2)));
    }

    static UpperBound<F> _log(PositiveUpperBound<F> const& x) {
        return UpperBound<F>(log(up,x._u));
    }

    static PositiveUpperBound<F> _max(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(max(x1._u,x2._u));
    }

    static PositiveUpperBound<F> _min(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(min(x1._u,x2._u));
    }

    static PositiveUpperBound<F> _abs(PositiveUpperBound<F> const& x) {
        return PositiveUpperBound<F>(x._u);
    }

    static Bool _same(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return x1._u == x2._u;
    }

    static Bool _refines(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return x1._u <= x2._u;
    }

    static PositiveUpperBound<F> _refinement(PositiveUpperBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveUpperBound<F>(min(x1._u,x2._u));
    }

    static OutputStream& _write(OutputStream& os, PositiveUpperBound<F> const& x) {
        return write(os,x.raw(),Bounds<F>::output_places,upward);
    }

    static InputStream& _read(InputStream& is, PositiveUpperBound<F>& x) {
        UpperBound<F> xu; is >> xu; x=PositiveUpperBound<F>(xu); return is;
    }

};



template<class F> struct Operations<PositiveLowerBound<F>> {
    static PositiveLowerBound<F> _nul(PositiveLowerBound<F> const& x) {
        return PositiveLowerBound<F>(nul(x._l));
    }

    static PositiveLowerBound<F> _sqr(PositiveLowerBound<F> const& x) {
        return PositiveLowerBound<F>(mul(down,x._l,x._l));
    }

    static PositiveUpperBound<F> _rec(PositiveLowerBound<F> const& x) {
        return PositiveUpperBound<F>(rec(up,x._l));
    }

    static PositiveLowerBound<F> _rec(PositiveUpperBound<F> const& x) {
        return PositiveLowerBound<F>(rec(down,x._u));
    }

    static PositiveLowerBound<F> _add(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(add(down,x1._l,x2._l));
    }

    static PositiveLowerBound<F> _mul(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(mul(down,x1._l,x2._l));
    }

    static PositiveLowerBound<F> _div(PositiveLowerBound<F> const& x1, PositiveUpperBound<F> const& x2) {
        return PositiveLowerBound<F>(div(down,x1._l,x2._u));
    }

    static PositiveLowerBound<F> _pow(PositiveLowerBound<F> const& x1, Nat m2) {
        return PositiveLowerBound<F>(pow(down,x1._l,static_cast<Int>(m2)));
    }

    static LowerBound<F> _log(PositiveLowerBound<F> const& x) {
        return LowerBound<F>(log(down,x._l));
    }

    static PositiveLowerBound<F> _max(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(max(x1._l,x2._l));
    }

    static PositiveLowerBound<F> _min(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(min(x1._l,x2._l));
    }

    static PositiveLowerBound<F> _abs(PositiveLowerBound<F> const& x) {
        return PositiveLowerBound<F>(x._l);
    }

    static Bool _same(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return x1._l == x2._l;
    }

    static Bool _refines(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return x1._l >= x2._l;
    }

    static PositiveLowerBound<F> _refinement(PositiveLowerBound<F> const& x1, PositiveLowerBound<F> const& x2) {
        return PositiveLowerBound<F>(max(x1._l,x2._l));
    }

    static OutputStream& _write(OutputStream& os, PositiveLowerBound<F> const& x) {
        return write(os,x.raw(),Bounds<F>::output_places,upward);
    }

    static InputStream& _read(InputStream& is, PositiveLowerBound<F>& x) {
        LowerBound<F> xu; is >> xu; x=PositiveLowerBound<F>(xu); return is;
    }

};

template<class F> struct Operations<PositiveBounds<F>> {
    static PositiveBounds<F> _nul(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(nul(x._l),nul(x._u)); }
    static PositiveBounds<F> _sqr(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(mul(down,x._l,x._l),mul(up,x._u,x._u)); }
    static PositiveBounds<F> _rec(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(rec(down,x._u),rec(up,x._l)); }
    static PositiveBounds<F> _add(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(add(down,x1._l,x2._l),add(up,x1._u,x2._u)); }
    static PositiveBounds<F> _mul(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(mul(down,x1._l,x2._l),mul(up,x1._u,x2._u)); }
    static PositiveBounds<F> _div(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(div(down,x1._l,x2._u),div(up,x1._u,x2._l)); }
    static PositiveBounds<F> _pow(PositiveBounds<F> const& x1, Nat m2) {
        return PositiveBounds<F>(pow(down,x1._l,static_cast<Int>(m2)),pow(up,x1._u,static_cast<Int>(m2))); }
    static PositiveBounds<F> _pow(PositiveBounds<F> const& x1, Int n2) {
        if(n2>=0) { return _pow(x1,Nat(n2)); } else { return _rec(_pow(x1,Nat(-n2))); } }
    static Bounds<F> _log(PositiveBounds<F> const& x) {
        return Bounds<F>(log(down,x._l),log(up,x._u)); }
    static PositiveBounds<F> _max(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(max(x1._l,x2._l),max(x1._u,x2._u)); }
    static PositiveBounds<F> _min(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return PositiveBounds<F>(min(x1._l,x2._l),min(x1._u,x2._u)); }
    static PositiveBounds<F> _abs(PositiveBounds<F> const& x) {
        return PositiveBounds<F>(x._l,x._u); }
    static Bool _same(PositiveBounds<F> const& x1, PositiveBounds<F> const& x2) {
        return x1._l == x2._l && x1._u == x2._u; }
    static OutputStream& _write(OutputStream& os, PositiveBounds<F> const& x) {
        return os << static_cast<Bounds<F>const&>(x); }
    static InputStream& _read(InputStream& is, PositiveBounds<F>& x) {
        Bounds<F> xb; is >> xb; x=PositiveBounds<F>(xb); return is; }
};

template<class F> struct Operations<Error<F>> {
    static OutputStream& _write(OutputStream& os, Error<F> const& x) {
        return write(os,x.raw(),DecimalPrecision{Error<F>::output_places},upward);
    }

    static InputStream& _read(InputStream& is, Error<F>& x) {
        UpperBound<F> xu; is >> xu; x=Error<F>(xu); return is;
    }
};





template<> Nat integer_cast<Nat,FloatDPApproximation>(FloatDPApproximation const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,FloatDPApproximation>(FloatDPApproximation const& x) {
    return std::round(x.get_d()); }

template<> Nat integer_cast<Nat,FloatDPLowerBound>(FloatDPLowerBound const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,FloatDPLowerBound>(FloatDPLowerBound const& x) {
    return std::round(x.get_d()); }

template<> Int integer_cast<Int,FloatDPBounds>(FloatDPBounds const& x) {
    return std::round((x.lower().get_d()+x.upper().get_d())/2); }

template<> Nat integer_cast<Nat,FloatMPApproximation>(FloatMPApproximation const& x) {
    return std::round(x.get_d()); }
template<> Int integer_cast<Int,FloatMPApproximation>(FloatMPApproximation const& x) {
    return std::round(x.get_d()); }




template<class F> Approximation<F> _make_float(Number<ApproximateTag> x) { return Approximation<F>(x); }
template<class F> LowerBound<F> _make_float(Number<ValidatedLowerTag> x) { return LowerBound<F>(x); }
template<class F> UpperBound<F> _make_float(Number<ValidatedUpperTag> x) { return UpperBound<F>(x); }
template<class F> Bounds<F> _make_float(Number<ValidatedTag> x) { return Bounds<F>(x); }
template<class F> Bounds<F> _make_float(Number<EffectiveTag> x) { return Bounds<F>(x); }
template<class F> Bounds<F> _make_float(Number<ExactTag> x) { return Bounds<F>(x); }
template<class F> Bounds<F> _make_float(Real r) { return Bounds<F>(r); }
template<class F> Bounds<F> _make_float(Rational q) { return Bounds<F>(q); }
template<class F> Value<F> _make_float(Integer z) { return Value<F>(z); }


template class Approximation<FloatDP>;
template class LowerBound<FloatDP>;
template class UpperBound<FloatDP>;
template class Bounds<FloatDP>;
template class Ball<FloatDP>;
template class Value<FloatDP>;

template class Approximation<FloatMP>;
template class LowerBound<FloatMP>;
template class UpperBound<FloatMP>;
template class Bounds<FloatMP>;
template class Ball<FloatMP>;
template class Value<FloatMP>;

template class Ball<FloatMP, FloatDP>;

template class Operations<FloatDPApproximation>;
template class Operations<FloatDPLowerBound>;
template class Operations<FloatDPUpperBound>;
template class Operations<FloatDPBounds>;
template class Operations<FloatDPBall>;
template class Operations<FloatDPValue>;
template class Operations<PositiveFloatDPApproximation>;
template class Operations<PositiveFloatDPLowerBound>;
template class Operations<PositiveFloatDPUpperBound>;
template class Operations<PositiveFloatDPBounds>;
template class Operations<FloatDPError>;

template class Operations<FloatMPApproximation>;
template class Operations<FloatMPLowerBound>;
template class Operations<FloatMPUpperBound>;
template class Operations<FloatMPBounds>;
template class Operations<FloatMPBall>;
template class Operations<FloatMPValue>;
template class Operations<PositiveFloatMPApproximation>;
template class Operations<PositiveFloatMPLowerBound>;
template class Operations<PositiveFloatMPUpperBound>;
template class Operations<PositiveFloatMPBounds>;
template class Operations<FloatMPError>;

template class Operations<FloatMDPBall>;

PositiveFloatDPApproximation mag(FloatDPApproximation const& x) { return Operations<FloatDPApproximation>::_mag(x); }
PositiveFloatDPApproximation mig(FloatDPApproximation const& x) { return Operations<FloatDPApproximation>::_mig(x); }
FloatDPApproximation round(FloatDPApproximation const& x) { return Operations<FloatDPApproximation>::_round(x); }
Bool same(FloatDPApproximation const& x1, FloatDPApproximation const& x2) { return Operations<FloatDPApproximation>::_same(x1,x2); }


PositiveFloatMPApproximation mag(FloatMPApproximation const& x) { return Operations<FloatMPApproximation>::_mag(x); }
PositiveFloatMPApproximation mig(FloatMPApproximation const& x) { return Operations<FloatMPApproximation>::_mig(x); }
FloatMPApproximation round(FloatMPApproximation const& x) { return Operations<FloatMPApproximation>::_round(x); }
Bool same(FloatMPApproximation const& x1, FloatMPApproximation const& x2) { return Operations<FloatMPApproximation>::_same(x1,x2); }



Bool same(FloatDPLowerBound const& x1, FloatDPLowerBound const& x2) { return Operations<FloatDPLowerBound>::_same(x1,x2); }
Bool refines(FloatDPLowerBound const& x1, FloatDPLowerBound const& x2) { return Operations<FloatDPLowerBound>::_refines(x1,x2); }
FloatDPLowerBound refinement(FloatDPLowerBound const& x1, FloatDPLowerBound const& x2) { return Operations<FloatDPLowerBound>::_refinement(x1,x2); }


Bool same(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return Operations<FloatMPLowerBound>::_same(x1,x2); }
Bool refines(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return Operations<FloatMPLowerBound>::_refines(x1,x2); }
FloatMPLowerBound refinement(FloatMPLowerBound const& x1, FloatMPLowerBound const& x2) { return Operations<FloatMPLowerBound>::_refinement(x1,x2); }




Bool same(FloatDPUpperBound const& x1, FloatDPUpperBound const& x2) { return Operations<FloatDPUpperBound>::_same(x1,x2); }
Bool refines(FloatDPUpperBound const& x1, FloatDPUpperBound const& x2) { return Operations<FloatDPUpperBound>::_refines(x1,x2); }
FloatDPUpperBound refinement(FloatDPUpperBound const& x1, FloatDPUpperBound const& x2) { return Operations<FloatDPUpperBound>::_refinement(x1,x2); }


Bool same(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return Operations<FloatMPUpperBound>::_same(x1,x2); }
Bool refines(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return Operations<FloatMPUpperBound>::_refines(x1,x2); }
FloatMPUpperBound refinement(FloatMPUpperBound const& x1, FloatMPUpperBound const& x2) { return Operations<FloatMPUpperBound>::_refinement(x1,x2); }



PositiveFloatDPUpperBound mag(FloatDPBounds const& x) { return Operations<FloatDPBounds>::_mag(x); }
PositiveFloatDPLowerBound mig(FloatDPBounds const& x) { return Operations<FloatDPBounds>::_mig(x); }
FloatDPBounds round(FloatDPBounds const& x) { return Operations<FloatDPBounds   >::_round(x); }

Bool same(FloatDPBounds const& x1, FloatDPBounds const& x2) { return Operations<FloatDPBounds>::_same(x1,x2); }
Bool models(FloatDPBounds const& x1, FloatDPValue const& x2) { return Operations<FloatDPBounds>::_models(x1,x2); }
Bool refines(FloatDPBounds const& x1, FloatDPBounds const& x2) { return Operations<FloatDPBounds>::_refines(x1,x2); }
Bool consistent(FloatDPBounds const& x1, FloatDPBounds const& x2) { return Operations<FloatDPBounds>::_consistent(x1,x2); }
Bool inconsistent(FloatDPBounds const& x1, FloatDPBounds const& x2) { return Operations<FloatDPBounds>::_inconsistent(x1,x2); }
FloatDPBounds refinement(FloatDPBounds const& x1, FloatDPBounds const& x2) { return Operations<FloatDPBounds>::_refinement(x1,x2); }


PositiveFloatMPUpperBound mag(FloatMPBounds const& x) { return Operations<FloatMPBounds>::_mag(x); }
PositiveFloatMPLowerBound mig(FloatMPBounds const& x) { return Operations<FloatMPBounds>::_mig(x); }
FloatMPBounds round(FloatMPBounds const& x) { return Operations<FloatMPBounds>::_round(x); }

Bool same(FloatMPBounds const& x1, FloatMPBounds const& x2) { return Operations<FloatMPBounds>::_same(x1,x2); }
Bool models(FloatMPBounds const& x1, FloatMPValue const& x2) { return Operations<FloatMPBounds>::_models(x1,x2); }
Bool refines(FloatMPBounds const& x1, FloatMPBounds const& x2) { return Operations<FloatMPBounds>::_refines(x1,x2); }
Bool consistent(FloatMPBounds const& x1, FloatMPBounds const& x2) { return Operations<FloatMPBounds>::_consistent(x1,x2); }
Bool inconsistent(FloatMPBounds const& x1, FloatMPBounds const& x2) { return Operations<FloatMPBounds>::_inconsistent(x1,x2); }
FloatMPBounds refinement(FloatMPBounds const& x1, FloatMPBounds const& x2) { return Operations<FloatMPBounds>::_refinement(x1,x2); }



PositiveFloatDPUpperBound mag(FloatDPBall const& x) { return Operations<FloatDPBall>::_mag(x); }
PositiveFloatDPLowerBound mig(FloatDPBall const& x) { return Operations<FloatDPBall>::_mig(x); }

Bool same(FloatDPBall const& x1, FloatDPBall const& x2) { return Operations<FloatDPBall>::_same(x1,x2); }
Bool models(FloatDPBall const& x1, FloatDPValue const& x2) { return Operations<FloatDPBall>::_models(x1,x2); }
Bool refines(FloatDPBall const& x1, FloatDPBall const& x2) { return Operations<FloatDPBall>::_refines(x1,x2); }
Bool consistent(FloatDPBall const& x1, FloatDPBall const& x2) { return Operations<FloatDPBall>::_consistent(x1,x2); }
Bool inconsistent(FloatDPBall const& x1, FloatDPBall const& x2) { return Operations<FloatDPBall>::_inconsistent(x1,x2); }
FloatDPBall refinement(FloatDPBall const& x1, FloatDPBall const& x2) { return Operations<FloatDPBall>::_refinement(x1,x2); }


PositiveFloatMPUpperBound mag(FloatMPBall const& x) { return Operations<FloatMPBall>::_mag(x); }
PositiveFloatMPLowerBound mig(FloatMPBall const& x) { return Operations<FloatMPBall>::_mig(x); }

Bool same(FloatMPBall const& x1, FloatMPBall const& x2) { return Operations<FloatMPBall>::_same(x1,x2); }
Bool models(FloatMPBall const& x1, FloatMPValue const& x2) { return Operations<FloatMPBall>::_models(x1,x2); }
Bool refines(FloatMPBall const& x1, FloatMPBall const& x2) { return Operations<FloatMPBall>::_refines(x1,x2); }
Bool consistent(FloatMPBall const& x1, FloatMPBall const& x2) { return Operations<FloatMPBall>::_consistent(x1,x2); }
Bool inconsistent(FloatMPBall const& x1, FloatMPBall const& x2) { return Operations<FloatMPBall>::_inconsistent(x1,x2); }
FloatMPBall refinement(FloatMPBall const& x1, FloatMPBall const& x2) { return Operations<FloatMPBall>::_refinement(x1,x2); }



FloatDPError mag(FloatDPValue const& x) { return Operations<FloatDPValue>::_mag(x); }
FloatMPError mag(FloatMPValue const& x) { return Operations<FloatMPValue>::_mag(x); }
Bool same(FloatDPValue const& x1, FloatDPValue const& x2) { return Operations<FloatDPValue>::_same(x1,x2); }
Bool same(FloatMPValue const& x1, FloatMPValue const& x2) { return Operations<FloatMPValue>::_same(x1,x2); }


PositiveFloatDPValue hlf(PositiveFloatDPValue const& x) { return PositiveFloatDPValue(hlf(x._v)); }
PositiveFloatMPValue hlf(PositiveFloatMPValue const& x) { return PositiveFloatMPValue(hlf(x._v)); }

FloatDPError operator/(FloatDPError const& x1, PositiveFloatDPLowerBound const& x2) {
    return FloatDPError(div(up,x1._e,x2._l)); }

FloatDPUpperBound operator*(FloatDPUpperBound const& x1, Real const& y2) {
    FloatDPUpperBound x2(y2,x1.precision()); return FloatDPUpperBound(mul(up,x1._u,x2._u)); }


FloatDPValue midpoint(FloatDPBounds const& x) { return x.value(); }

template<> String class_name<FloatDPApproximation>() { return "FloatDPApproximation"; }
template<> String class_name<FloatDPLowerBound>() { return "FloatDPLowerBound"; }
template<> String class_name<FloatDPUpperBound>() { return "FloatDPUpperBound"; }
template<> String class_name<FloatDPBounds>() { return "FloatDPBounds"; }
template<> String class_name<FloatDPBall>() { return "FloatDPBall"; }
template<> String class_name<FloatDPValue>() { return "FloatDPValue"; }
template<> String class_name<FloatMPApproximation>() { return "FloatMPApproximation"; }
template<> String class_name<FloatMPLowerBound>() { return "FloatMPLowerBound"; }
template<> String class_name<FloatMPUpperBound>() { return "FloatMPUpperBound"; }
template<> String class_name<FloatMPBounds>() { return "FloatMPBounds"; }
template<> String class_name<FloatMPBall>() { return "FloatMPBall"; }
template<> String class_name<FloatMPValue>() { return "FloatMPValue"; }

template<> String class_name<FloatMDPBall>() { return "FloatMDPBall"; }

} // namespace Ariadne
