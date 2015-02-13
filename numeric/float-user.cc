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

template<class PR> Nat Float<Approximate,PR>::output_precision = 4;
template<class PR> Nat Float<Bounded,PR>::output_precision;
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


template<class PR> Float<Exact,PR>::Float(Integer const& z) : Float(z,RawFloatType::get_default_precision()) { }

template<class PR> Float<Exact,PR>::Float(Integer const& z, PR pr)
    : _v(Rational(z),RawFloatType::to_nearest)
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

template<class PR> Fuzzy eq(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return x1._a==x2._a; }
template<class PR> Fuzzy leq(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return x1._a<=x2._a; }

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

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Approximate,PR> const& x);
template<class PR> InputStream& operator>>(InputStream& is, Float<Approximate,PR>& x);



template<class PR> Float<Lower,PR> max(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Lower,PR>(max_exact(x1._l,x2._l)); }
template<class PR> Float<Lower,PR> min(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Lower,PR>(min_exact(x1._l,x2._l)); }

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

template<class PR> Float<Lower,PR> sub(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Lower,PR>(sub_down(x1._l,x2._u)); }

template<class PR> Float<Lower,PR> mul(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    return Float<Lower,PR>(mul_down(x1.raw(),x2.raw())); }

template<class PR> Float<Lower,PR> div(Float<Lower,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Lower,PR>(div_down(x1.raw(),x2.raw())); }

template<class PR> Float<Lower,PR> pow(Float<Lower,PR> const& x, Nat m);

template<class PR> Float<Lower,PR> sqrt(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(sqrt_down(x.raw())); }

template<class PR> Float<Lower,PR> exp(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(exp_down(x.raw())); }

template<class PR> Float<Lower,PR> log(Float<Lower,PR> const& x) {
    return Float<Lower,PR>(log_down(x.raw())); }

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
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64::set_rounding_downward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.raw();
    Float64::set_rounding_mode(rnd);
    return os;
}

template<class PR> InputStream& operator>>(InputStream& is, Float<Lower,PR>& x);



template<class PR> Float<Upper,PR> max(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(max_exact(x1._u,x2._u)); }

template<class PR> Float<Upper,PR> min(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(min_exact(x1._u,x2._u)); }

template<class PR> Float<Upper,PR> abs(Float<Upper,PR> const& x) = delete;

template<class PR> Float<Upper,PR> nul(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(pos_exact(x._u)); }

template<class PR> Float<Upper,PR> pos(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(pos_exact(x._u)); }

template<class PR> Float<Upper,PR> neg(Float<Lower,PR> const& x) {
    return Float<Upper,PR>(neg_exact(x._l)); }

template<class PR> Float<Upper,PR> half(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(half_exact(x._u)); }

template<class PR> Float<Upper,PR> sqr(Float<Upper,PR> const& x);

template<class PR> Float<Upper,PR> rec(Float<Lower,PR> const& x) {
    return Float<Upper,PR>(rec_up(x.raw())); }

template<class PR> Float<Upper,PR> add(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return Float<Upper,PR>(add_up(x1._u,x2._u)); }

template<class PR> Float<Upper,PR> sub(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) {
    return Float<Upper,PR>(sub_up(x1._u,x2._l)); }

template<class PR> Float<Upper,PR> mul(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2);

template<class PR> Float<Upper,PR> div(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2);

template<class PR> Float<Upper,PR> pow(Float<Upper,PR> const& x, Nat m);

template<class PR> Float<Upper,PR> sqrt(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(sqrt_up(x.raw())); }

template<class PR> Float<Upper,PR> exp(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(exp_up(x.raw())); }

template<class PR> Float<Upper,PR> log(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(log_up(x.raw())); }

template<class PR> Float<Upper,PR> atan(Float<Upper,PR> const& x) {
    return Float<Upper,PR>(atan_up(x.raw())); }

template<class PR> Int integer_cast(Float<Upper,PR> const& x) { return static_cast<Int>(x._u.get_d()); }
template<class PR> Nat integer_cast(Float<Upper,PR> const& x) { return static_cast<Nat>(x._u.get_d()); }

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

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Upper,PR> const& x);
template<class PR> InputStream& operator>>(InputStream& is, Float<Upper,PR>& x);







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
        return Float<Bounded,PR>(static_cast<Float64>(0.0),max(neg(x.lower_raw()),x.upper_raw()));
    }
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
    const Float64& xl=x.lower_raw(); const Float64& xu=x.upper_raw();
    Float64 rl,ru;
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
    if(x._l>0 || x._u<0) {
        return Float<Bounded,PR>(rec_down(x._u),rec_up(x._l));
    } else {
        RawFloat<PR> rl=-inf; RawFloat<PR> ru=+inf;
        ARIADNE_THROW(DivideByZeroException,"BoundedFloat rec(BoundedFloat x)","x="<<x);
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
    const Float64& x1l=x1._l; const Float64& x1u=x1._u;
    const Float64& x2l=x2._l; const Float64& x2u=x2._u;
    Float64 rl,ru;
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
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
    const Float64& x1l=x1.lower_raw(); const Float64& x1u=x1.upper_raw();
    const Float64& i2l=x2.lower_raw(); const Float64& i2u=x2.upper_raw();
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
        ARIADNE_THROW(DivideByZeroException,"BoundedFloat div(BoundedFloat x1, BoundedFloat x2)","x1="<<x1<<", x2="<<x2);
        rl=-Float64::inf();
        ru=+Float64::inf();
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
    Float64 rl=pow_down(y.lower_raw(),m);
    Float64 ru=pow_up(y.upper_raw(),m);
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


template<class PR> Float<Bounded,PR> pi_val() { return Float<Bounded,PR>(pi_down(),pi_up()); }


template<class PR> Float<Bounded,PR> sin(Float<Bounded,PR> const& x)
{
    return cos(x-half(pi_val<PR>()));
}

template<class PR> Float<Bounded,PR> cos(Float<Bounded,PR> const& x)
{
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    Float64::RoundingModeType rnd = Float64::get_rounding_mode();

    static const Float<Exact,PR> two(2);

    if(x.error().raw()>2*pi_down()) { return Float<Bounded,PR>(-1.0,+1.0); }

    Float64 n=floor(x.lower_raw()/(2*pi_approx())+0.5);
    Float<Bounded,PR> y=x-two*Float<Exact,PR>(n)*pi_val<PR>();

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
    ARIADNE_NOT_IMPLEMENTED;
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
    return Float<Bounded,PR>(wl,wu);
}

template<class PR> Float<Bounded,PR> narrow(Float<Bounded,PR> const& x)
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
    return Float<Bounded,PR>(nl,nu);
}

template<class PR> Float<Bounded,PR> trunc(Float<Bounded,PR> const& x)
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
    return Float<Bounded,PR>(double(tl),double(tu));
}

template<class PR> Float<Bounded,PR> trunc(Float<Bounded,PR> const& x, Nat n)
{
    Float<Bounded,PR> _e=Float<Bounded,PR>(std::pow(2.0,52-(Int)n));
    Float<Bounded,PR> y=x+_e;
    return y-_e;
}

template<class PR> Int integer_cast(Float<Bounded,PR> const& x) { return static_cast<Int>(x.value_raw().get_d()); }
template<class PR> Nat integer_cast(Float<Bounded,PR> const& x) { return static_cast<Nat>(x.value_raw().get_d()); }

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
template<class PR> auto same(Float<Bounded,PR> const&, Float<Bounded,PR> const&) -> Bool;


template<class PR> OutputStream& operator<<(OutputStream& os, const Float<Bounded,PR>& x)
{
    //if(x.lower_raw()==x.upper_raw()) { return os << "{" << std::setprecision(Float<Bounded,PR>::output_precision) << x.lower_raw().get_d() << ; }
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    os << '{';
    Float64::set_rounding_downward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.lower().get_d();
    os << ':';
    Float64::set_rounding_upward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.upper().get_d();
    Float64::set_rounding_mode(rnd);
    os << '}';
    return os;

}

template<class PR> InputStream& operator>>(InputStream& is, Float<Bounded,PR>& x)
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
    const Float64& x1l=x1.lower_raw(); const Float64& x1u=x1.upper_raw();
    const Float64& x2v=x2.raw();
    Float64 rl,ru;
    if(x2v>=0.0) {
        rl=mul_down(x1l,x2v); ru=mul_up(x1u,x2v);
    } else {
        rl=mul_down(x1u,x2v); ru=mul_up(x1l,x2v);
    }
    return Float<Bounded,PR>(rl,ru);
}


template<class PR> Float<Bounded,PR> mul(Float<Exact,PR> const& x1, Float<Bounded,PR> const& x2) {
    const Float64& x1v=x1.raw();
    const Float64& x2l=x2.lower_raw(); const Float64& x2u=x2.upper_raw();
    Float64 rl,ru;
    if(x1v>=0.0) {
        rl=mul_down(x1v,x2l); ru=mul_up(x1v,x2u);
    } else {
        rl=mul_down(x1v,x2u); ru=mul_up(x1v,x2l);
    }
    return Float<Bounded,PR>(rl,ru);
}

template<class PR> Float<Bounded,PR> div(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2)
{
    const Float64& x1l=x1.lower_raw();
    const Float64& x1u=x1.upper_raw();
    const Float64& x2v=x2.raw();
    Float64 rl,ru;
    if(x2v>0) {
        rl=div_down(x1l,x2v); ru=div_up(x1u,x2v);
    } else if(x2v<0) {
        rl=div_down(x1u,x2v); ru=div_up(x1l,x2v);
    } else {
        ARIADNE_THROW(DivideByZeroException,"BoundedFloat div(BoundedFloat const& x1, ExactFloat x2)","x1="<<x1<<", x2="<<x2);
        rl=-Float64::inf();
        ru=+Float64::inf();
    }
    return Float<Bounded,PR>(rl,ru);
}


template<class PR> Float<Bounded,PR> div(Float<Exact,PR> const& x1, Float<Bounded,PR> const& x2)
{
    const Float64& x1v=x1.raw();
    const Float64& i2l=x2.lower_raw();
    const Float64& i2u=x2.upper_raw();
    Float64 rl,ru;
    if(i2l<=0 && i2u>=0) {
        ARIADNE_THROW(DivideByZeroException,"BoundedFloat div(ExactFloat const& x1, BoundedFloat x2)","x1="<<x1<<", x2="<<x2);
        rl=-Float64::inf();
        ru=+Float64::inf();
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

template<class PR> Float<Bounded,PR> pow(Float<Exact,PR> const& x1, Nat m2) {
    return pow(Float<Bounded,PR>(x1),m2); }

template<class PR> Float<Bounded,PR> pow(Float<Exact,PR> const& x1, Int n2) {
    return pow(Float<Bounded,PR>(x1),n2); }

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

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Exact,PR> const& x);
template<class PR> InputStream& operator>>(InputStream& is, Float<Exact,PR>& x);

template<class PR> Float<Exact,PR> make_exact(const Real& x);

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



template<class PR> Float<PositiveExact,PR> mag(Float<Exact,PR> const& x) {
    return Float<PositiveExact,PR>(abs(x.raw())); }
// FIXME: Unsafe since x may be negative
template<class PR> Float<PositiveUpper,PR> mag(Float<Upper,PR> const& x) {
    return Float<PositiveUpper,PR>(abs(x.raw())); }
template<class PR> Float<PositiveUpper,PR> mag(Float<Bounded,PR> const& x) {
    return Float<PositiveUpper,PR>(max(neg(x.lower_raw()),x.upper_raw())); }
template<class PR> Float<PositiveApproximate,PR> mag(Float<Approximate,PR> const& x) {
    return Float<PositiveApproximate,PR>(abs(x.raw())); }

template<class PR> Float<PositiveLower,PR> mig(Float<Bounded,PR> const& x) {
    Float64 r=max(x.lower_raw(),neg(x.upper_raw()));
    return Float<PositiveLower,PR>(max(r,nul(r))); }
template<class PR> Float<PositiveLower,PR> mig(Float<Lower,PR> const& x) {
    Float64 r=min(x.raw(),neg(x.raw()));
    return Float<PositiveLower,PR>(max(r,nul(r))); }

template<class PR> Float<Bounded,PR> make_bounds(Float<PositiveUpper,PR> const& _e) {
    return Float<Bounded,PR>(-_e.raw(),+_e.raw());
}

template<class PR> Float<Exact,PR> value(Float<Bounded,PR> const& x) {
    return Float<Exact,PR>(half_exact(add_near(x.lower_raw(),x.upper_raw())));
}

template<class PR> Float<PositiveUpper,PR> error(Float<Bounded,PR> const& x) {
    return Float<PositiveUpper,PR>(half_exact(sub_up(x.upper_raw(),x.lower_raw())));
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


template<class PR> Bool same(Float<Approximate,PR> const& x1, Float<Approximate,PR> const& x2) {
    return x1.raw()==x2.raw(); }

template<class PR> Bool same(Float<Lower,PR> const& x1, Float<Lower,PR> const& x2) {
    return x1.raw()==x2.raw(); }

template<class PR> Bool same(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    return x1.raw()==x2.raw(); }

template<class PR> Bool same(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1.lower_raw()==x2.lower_raw() && x1.upper_raw()==x2.upper_raw(); }

template<class PR> Bool same(Float<Exact,PR> const& x1, Float<Exact,PR> const& x2) {
    return x1.raw()==x2.raw(); }



#ifdef ARIADNE_ENABLE_SERIALIZATION
template<class PR, class A> Void serialize(A& _a, Float<Bounded,PR>& x, const Nat version) {
    _a & x.lower_raw() & x.upper_raw(); }
#endif

template<class PR> OutputStream& operator<<(OutputStream&, Float<Bounded,PR> const&);
template<class PR> InputStream& operator>>(InputStream&, Float<Bounded,PR>&);


template<class PR> Float<PositiveUpper,PR> operator+(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2);
template<class PR> Float<PositiveUpper,PR> operator-(Float<PositiveUpper,PR> const& x1, Float<Lower,PR> const& x2);
template<class PR> Float<PositiveUpper,PR> operator*(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2);
template<class PR> Float<PositiveUpper,PR> operator/(Float<PositiveUpper,PR> const& x1, Float<Lower,PR> const& x2);
template<class PR> Float<PositiveUpper,PR> pow(Float<PositiveUpper,PR> const& x, Nat m);
template<class PR> Float<PositiveUpper,PR> half(Float<PositiveUpper,PR> const& x);
template<class PR> Float<PositiveUpper,PR> max(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(max(x1.raw(),x2.raw())); }



template<class PR> Float<Exact,PR> make_exact(Real const& r) {
    Float<Approximate,PR> _a(r); return Float<Exact,PR>(_a.raw());
}

template<class PR> OutputStream& operator<<(OutputStream& os, Float<Exact,PR> const& x) {
    os << std::showpoint << std::setprecision(Float<Exact,PR>::output_precision) << x.raw();
    return os;
}



template<class PR> Float<PositiveUpper,PR> operator+(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(add_up(x1.raw(),x2.raw()));
}

template<class PR> Float<PositiveUpper,PR> operator-(Float<PositiveUpper,PR> const& x1, Float<Lower,PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=x2.raw());
    return Float<PositiveUpper,PR>(sub_up(x1.raw(),x2.raw()));
}

template<class PR> Float<PositiveUpper,PR> operator*(Float<PositiveUpper,PR> const& x1, Float<PositiveUpper,PR> const& x2) {
    return Float<PositiveUpper,PR>(mul_up(x1.raw(),x2.raw()));
}

template<class PR> Float<PositiveUpper,PR> operator/(Float<PositiveUpper,PR> const& x1, Float<Lower,PR> const& x2) {
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return Float<PositiveUpper,PR>(div_up(x1.raw(),x2.raw()));
}

template<class PR> Float<PositiveUpper,PR> pow(Float<PositiveUpper,PR> const& x, Nat m) {
    return Float<PositiveUpper,PR>(pow_up(x.raw(),m));
}

template<class PR> Float<PositiveUpper,PR> half(Float<PositiveUpper,PR> const& x) {
    return Float<PositiveUpper,PR>(half(x.raw()));
}




template<class PR> OutputStream& operator<<(OutputStream& os, Float<Upper,PR> const& x) {
    Float64::RoundingModeType rnd=Float64::get_rounding_mode();
    Float64::set_rounding_upward();
    os << std::showpoint << std::setprecision(Float<Bounded,PR>::output_precision) << x.raw();
    Float64::set_rounding_mode(rnd);
    return os;
}

template<class PR> Float<Upper,PR> sqr(Float<Upper,PR> const& x) {
    ARIADNE_PRECONDITION(x.raw()>=0.0);
    return Float<Upper,PR>(mul_up(x.raw(),x.raw()));
}

template<class PR> Float<Upper,PR> mul(Float<Upper,PR> const& x1, Float<Upper,PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0.0 && x2.raw() >= 0.0);
    return Float<Upper,PR>(mul_up(x1.raw(),x2.raw()));
}

template<class PR> Float<Upper,PR> div(Float<Upper,PR> const& x1, Float<Lower,PR> const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0.0 && x2.raw() > 0.0);
    return Float<Upper,PR>(div_up(x1.raw(),x2.raw()));
}

template<class PR> Float<Upper,PR> pow(Float<Upper,PR> const& x, Nat n) {
    ARIADNE_PRECONDITION(x.raw()>=0.0);
    return Float<Upper,PR>(pow_up(x.raw(),n));
}


template<class PR> Int integer_cast(Float<Lower,PR> const& x) { return static_cast<Int>(x._l.get_d()); }
template<class PR> Nat integer_cast(Float<Lower,PR> const& x) { return static_cast<Nat>(x._l.get_d()); }


template<class PR> OutputStream& operator<<(OutputStream& os, Float<Approximate,PR> const& x) {
    return os << std::showpoint << std::setprecision(Float<Approximate,PR>::output_precision) << x.raw();
}

template<class PR> Int integer_cast(Float<Approximate,PR> const& x) { return static_cast<Int>(x._a.get_d()); }
template<class PR> Nat integer_cast(Float<Approximate,PR> const& x) { return static_cast<Nat>(x._a.get_d()); }



// Operations relating to validated operations

//! \related Float<Bounded,PR> \brief Tests if \_a x1 provides tighter bounds than \_a x2.
template<class PR> Bool refines(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1.lower_raw()>=x2.lower_raw() && x1.upper_raw()<=x2.upper_raw(); }

//! \related Float<Bounded,PR> \brief The common refinement of \_a x1 and \_a x2.
template<class PR> Float<Bounded,PR> refinement(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return Float<Bounded,PR>(max(x1.lower_raw(),x2.lower_raw()),min(x1.upper_raw(),x2.upper_raw())); }

//! \related Float<Bounded,PR> \brief Tests if \_a x1 and \_a x2 are consistent with representing the same number.
template<class PR> Bool consistent(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return x1.lower_raw()<=x2.upper_raw() && x1.upper_raw()>=x2.lower_raw(); }

//! \related Float<Bounded,PR> \brief  Tests if \_a x1 and \_a x2 are inconsistent with representing the same number.
template<class PR> Bool inconsistent(Float<Bounded,PR> const& x1, Float<Bounded,PR> const& x2) {
    return not consistent(x1,x2); }

//! \related Float<Bounded,PR> \brief  Tests if \_a x1 is a model for the exact value \_a x2. number.
template<class PR> Bool models(Float<Bounded,PR> const& x1, Float<Exact,PR> const& x2) {
    return x1.lower_raw()<=x2.raw() && x1.upper_raw()>=x2.raw(); }



template<class P, class PR> Float<P,PR> max(Float<P,PR> const& x1, Float<P,PR> const& x2) { return max<PR>(x1,x2); }
template<class P, class PR> Float<P,PR> min(Float<P,PR> const& x1, Float<P,PR> const& x2) { return min<PR>(x1,x2); }
template<class P, class PR> Float<P,PR> abs(Float<P,PR> const& x) { return abs<PR>(x); }

template<class P, class PR> Float<Widen<P>,PR> add(Float<P,PR> const& x1, Float<P,PR> const& x2) { return add<PR>(x1,x2); }
template<class P, class PR> Float<Widen<P>,PR> sub(Float<P,PR> const& x1, Float<Negated<P>,PR> const& x2) { return sub<PR>(x1,x2); }
template<class P, class PR> Float<Widen<P>,PR> mul(Float<P,PR> const& x1, Float<P,PR> const& x2) { return mul<PR>(x1,x2); }
template<class P, class PR> Float<Widen<P>,PR> div(Float<P,PR> const& x1, Float<Inverted<P>,PR> const& x2) { return div<PR>(x1,x2); }

template<class P, class PR> Float<Widen<P>,PR> pow(Float<P,PR> const& x, Nat m) { return pow<PR>(x,m); }
template<class P, class PR> Float<Widen<P>,PR> pow(Float<P,PR> const& x, Int n) { return pow<PR>(x,n); }

template<class P, class PR> Float<Widen<P>,PR> sqrt(Float<P,PR> const& x) { return sqrt<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> exp(Float<P,PR> const& x) { return exp<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> log(Float<P,PR> const& x) { return log<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> sin(Float<P,PR> const& x) { return sin<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> cos(Float<P,PR> const& x) { return cos<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> tan(Float<P,PR> const& x) { return tan<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> asin(Float<P,PR> const& x) { return asin<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> acos(Float<P,PR> const& x) { return acos<PR>(x); }
template<class P, class PR> Float<Widen<P>,PR> atan(Float<P,PR> const& x) { return atan<PR>(x); }

template<class P, class PR> OutputStream& operator<<(OutputStream& os, Float<P,PR> const& x) {
    return operator<<(os,x); }



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

template Float<Lower,PrecisionMP>::Float(Float<Bounded,PrecisionMP>const&);
template Float<Upper,PrecisionMP>::Float(Float<Bounded,PrecisionMP>const&);


template<class P, class PR> std::size_t instantiate_float() {
    using std::size_t;
    typedef OutputStream OS;
    typedef Float<P,PR> X;
    typedef Float<P,PR> const& XCRef;

    typedef Float<Negated<P>,PR> NEGX;
    typedef Float<Inverted<P>,PR> INVX;
    typedef decltype(abs(declval<X>())) ABSX;
    typedef decltype(add(declval<X>(),declval<X>())) WX;

    typedef decltype(mag(declval<X>())) MAGX;

    typedef decltype(operator==(declval<X>(),declval<X>())) XE;
    typedef decltype(operator!=(declval<X>(),declval<X>())) XNE;
    typedef decltype(operator<=(declval<X>(),declval<X>())) XL;

    return (size_t)(X(*)(XCRef,XCRef))&max + (size_t)(X(*)(XCRef,XCRef))&min + (size_t)(ABSX(*)(XCRef))&abs
         + (size_t)(WX(*)(X const&,X const&))&add + (size_t)(WX(*)(X const&,NEGX const&))&sub
         + (size_t)(WX(*)(X const&,X const&))&mul + (size_t)(WX(*)(X const&,INVX const&))&div
         + (size_t)(X(*)(XCRef))&nul + (size_t)(X(*)(XCRef))&pos + (size_t)(NEGX(*)(XCRef))&neg
         + (size_t)(WX(*)(XCRef))&sqr + (size_t)(WX(*)(XCRef))&rec + (size_t)(INVX(*)(XCRef))&half
         + (size_t)(WX(*)(XCRef,Int))&pow + (size_t)(WX(*)(XCRef,Nat))&pow
         + (size_t)(WX(*)(XCRef))&sqrt + (size_t)(WX(*)(XCRef))&exp  + (size_t)(WX(*)(XCRef))&log
         + (size_t)(WX(*)(XCRef))&sin  + (size_t)(WX(*)(XCRef))&cos  + (size_t)(WX(*)(XCRef))&tan
         + (size_t)(WX(*)(XCRef))&asin + (size_t)(WX(*)(XCRef))&acos + (size_t)(WX(*)(XCRef))&atan
         + (size_t)(MAGX(*)(XCRef))&mag
         + ((size_t)(XE(*)(XCRef,XCRef))&operator==) + ((size_t)(XNE(*)(XCRef,XCRef))&operator!=)
         + ((size_t)(XL(*)(XCRef,XCRef))&operator<=) + ((size_t)(XL(*)(XCRef,XCRef))&operator>=)
         + ((size_t)(XL(*)(XCRef,XCRef))&operator< ) + ((size_t)(XL(*)(XCRef,XCRef))&operator> )
         + (size_t)(OS&(*)(OS&,XCRef)) operator<<;

}

template<class PR> std::size_t instantiate_floats() {
    return instantiate_float<Approximate,PR>()
        + instantiate_float<Bounded,PR>()
//        + instantiate_float<Metric,PR>()
        + instantiate_float<Exact,PR>();
        ;
}

template std::size_t instantiate_floats<Precision64>();
template OutputStream& operator<<(OutputStream&, ExactFloat64 const&);
template InputStream& operator>>(InputStream&, BoundedFloat64&);
template Bool same(BoundedFloat64 const&, BoundedFloat64 const&);


} // namespace Ariadne
