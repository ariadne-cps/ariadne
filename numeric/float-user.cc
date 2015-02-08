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

const ExactFloat64 infty = ExactFloat64(Float64::inf());

template<class PR> Float<Exact,PR>::Float(Integer const& z)
    : v(z.get_si())
{
    int n=z.get_si();
    ARIADNE_PRECONDITION(z==n);
}

template<class PR> Float<Validated,PR>::Float(Number<Validated> const& x)
    : Float(x.get(Validated(),FLT::get_default_precision())) {
}

template<class PR> Float<Validated,PR>::Float(Number<Validated> const& x, PR pr)
    : Float(x.get(Validated(),pr)) {
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


template<class PR> Float<Validated,PR>::Float(const Dyadic& b) : Float<Validated,PR>(b.operator Rational()) { }

template<class PR> Float<Validated,PR>::Float(const Decimal& d) : Float<Validated,PR>(d.operator Rational()) { }


template<class PR> Float<Validated,PR>::Float(const Integer& z) : Float<Validated,PR>(Rational(z)) {
}

template<class PR> Float<Validated,PR>::Float(const Rational& q) : Float<Validated,PR>(q,q) {
}

template<class PR> Float<Validated,PR>::Float(const Rational& ql, const Rational& qu) : l(ql.get_d()), u(qu.get_d())  {
    while(Rational(l)>ql) { l=next_down(l); }
    while(Rational(u)<qu) { u=next_up(u); }
}

template<class PR> Float<Validated,PR> Float<Exact,PR>::pm(Float<Error,PR> e) const {
    Float<Exact,PR> const& v=*this; return Float<Validated,PR>(v-e,v+e);
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



ExactFloat64 make_exact(Real const& r) {
    ApproximateFloat64 a(r); return ExactFloat64(a.raw());
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
    ValidatedFloat64 e=ValidatedFloat64(std::pow(2.0,52-(Int)n));
    ValidatedFloat64 y=x+e;
    return y-e;
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
    Float64 l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    x.l=l; x.u=u;
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

template<> Int integer_cast(UpperFloat64 const& x) { return static_cast<Int>(x.u.get_d()); }
template<> Nat integer_cast(UpperFloat64 const& x) { return static_cast<Nat>(x.u.get_d()); }




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

template<> Int integer_cast(LowerFloat64 const& x) { return static_cast<Int>(x.l.get_d()); }
template<> Nat integer_cast(LowerFloat64 const& x) { return static_cast<Nat>(x.l.get_d()); }


OutputStream& operator<<(OutputStream& os, ApproximateFloat64 const& x) {
    return os << std::showpoint << std::setprecision(ApproximateFloat64::output_precision) << x.raw();
}

template<> Int integer_cast(ApproximateFloat64 const& x) { return static_cast<Int>(x.a.get_d()); }
template<> Nat integer_cast(ApproximateFloat64 const& x) { return static_cast<Nat>(x.a.get_d()); }




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
template class Float<Validated,Precision64>;
template class Float<Exact,Precision64>;


} // namespace Ariadne
