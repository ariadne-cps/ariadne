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

Nat ApproximateFloat::output_precision = 4;
Nat ValidatedFloat::output_precision = 8;
Nat ExactFloat::output_precision = 16;

const ExactFloat infty = ExactFloat(Float::inf());


ExactFloat::ExactFloat(Integer const& z)
    : v(z.get_si())
{
    int n=z.get_si();
    ARIADNE_PRECONDITION(z==n);
}

ExactFloat make_exact(Real const& r) {
    ApproximateFloat a(r); return ExactFloat(a.raw());
}

OutputStream& operator<<(OutputStream& os, ExactFloat const& x) {
    os << std::showpoint << std::setprecision(ExactFloat::output_precision) << x.raw();
    return os;
}

namespace {
inline Float cos_down(Float const& x) { Float::set_rounding_downward(); Float y=cos_rnd(x); set_default_rounding(); return y; }
inline Float cos_up(Float const& x) { Float::set_rounding_upward(); Float y=cos_rnd(x); set_default_rounding(); return y; }
} // namespace

ValidatedFloat::ValidatedFloat(Number<Validated> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat widen(ValidatedFloat const& x)
{
    Float::RoundingModeType rm=Float::get_rounding_mode();
    const Float& xl=x.lower_raw();
    const Float& xu=x.upper_raw();
    const Float m=std::numeric_limits<float>::min();
    Float::set_rounding_upward();
    Float wu=add(xu,m);
    Float mwl=add(neg(xl),m);
    Float wl=neg(mwl);
    Float::set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return ValidatedFloat(wl,wu);
}

ValidatedFloat narrow(ValidatedFloat const& x)
{
    Float::RoundingModeType rm=Float::get_rounding_mode();
    const Float& xl=x.lower_raw();
    const Float& xu=x.upper_raw();
    const Float m=std::numeric_limits<float>::min();
    Float::set_rounding_upward();
    Float mnu=add(neg(xu),m);
    Float nu=neg(mnu);
    Float nl=add(xl,m);
    Float::set_rounding_mode(rm);
    assert(xl<nl); assert(nu<xu);
    return ValidatedFloat(nl,nu);
}

ValidatedFloat trunc(ValidatedFloat const& x)
{
    Float::RoundingModeType rm=Float::get_rounding_mode();
    const double& xl=x.lower_raw().get_d();
    const double& xu=x.upper_raw().get_d();
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { Float::set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { Float::set_rounding_downward(); tl-=fm; }
    Float::set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return ValidatedFloat(double(tl),double(tu));
}

ValidatedFloat trunc(ValidatedFloat const& x, Nat n)
{
    ValidatedFloat e=ValidatedFloat(std::pow(2.0,52-(Int)n));
    ValidatedFloat y=x+e;
    return y-e;
}

ValidatedFloat rec(ValidatedFloat const& x)
{
    const Float& xl=x.lower_raw();
    const Float& xu=x.upper_raw();
    Float rl,ru;
    if(xl>0 || xu<0) {
        rl=rec_down(xu);
        ru=rec_up(xl);
    } else {
        rl=-inf;
        ru=+inf;
        ARIADNE_THROW(DivideByZeroException,"ValidatedFloat rec(ValidatedFloat ivl)","ivl="<<x);
    }
    return ValidatedFloat(rl,ru);
}


ValidatedFloat mul(ValidatedFloat const& x1, ValidatedFloat const& x2)
{
    const Float& x1l=x1.lower_raw();
    const Float& x1u=x1.upper_raw();
    const Float& i2l=x2.lower_raw();
    const Float& i2u=x2.upper_raw();
    Float rl,ru;
    Float::RoundingModeType rnd=Float::get_rounding_mode();
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
            Float::set_rounding_downward();
            rl=min(mul(x1u,i2l),mul(x1l,i2u));
            Float::set_rounding_upward();
            ru=max(mul(x1l,i2l),mul(x1u,i2u));
        }
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat mul(ValidatedFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& x1l=x1.lower_raw();
    const Float& x1u=x1.upper_raw();
    const Float& x2v=x2.raw();
    Float rl,ru;
    if(x2>=0) {
        rl=mul_down(x1l,x2v); ru=mul_up(x1u,x2v);
    } else {
        rl=mul_down(x1u,x2v); ru=mul_up(x1l,x2v);
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat mul(ExactFloat const& x1, ValidatedFloat const& x2) {
    return mul(x2,x1);
}


ValidatedFloat div(ValidatedFloat const& x1, ValidatedFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& x1l=x1.lower_raw();
    const Float& x1u=x1.upper_raw();
    const Float& i2l=x2.lower_raw();
    const Float& i2u=x2.upper_raw();
    Float rl,ru;
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
        // ARIADNE_THROW(DivideByZeroException,"ValidatedFloat div(ValidatedFloat ivl1, ValidatedFloat ivl2)","ivl1="<<x1<<", ivl2="<<x2);
        rl=-Float::inf();
        ru=+Float::inf();
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



ValidatedFloat div(ValidatedFloat const& x1, ExactFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& x1l=x1.lower_raw();
    const Float& x1u=x1.upper_raw();
    const Float& x2v=x2.raw();
    Float rl,ru;
    if(x2v>0) {
        rl=div_down(x1l,x2v); ru=div_up(x1u,x2v);
    } else if(x2v<0) {
        rl=div_down(x1u,x2v); ru=div_up(x1l,x2v);
    } else {
        rl=-Float::inf();
        ru=+Float::inf();
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat div(ExactFloat const& x1, ValidatedFloat const& x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& x1v=x1.raw();
    const Float& i2l=x2.lower_raw();
    const Float& i2u=x2.upper_raw();
    Float rl,ru;
    if(i2l<=0 && i2u>=0) {
        ARIADNE_THROW(DivideByZeroException,"ValidatedFloat div(Float const& x1, ValidatedFloat ivl2)","x1="<<x1<<", ivl2="<<x2);
        rl=-Float::inf();
        ru=+Float::inf();
    } else if(x1v>=0) {
        rl=div_down(x1v,i2u); ru=div_up(x1v,i2l);
    } else {
        rl=div_down(x1v,i2l); ru=div_up(x1v,i2u);
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat sqr(ValidatedFloat const& x)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& xl=x.lower_raw();
    const Float& xu=x.upper_raw();
    Float rl,ru;
    if(xl>0.0) {
        rl=mul_down(xl,xl);
        ru=mul_up(xu,xu);
    } else if(xu<0.0) {
        rl=mul_down(xu,xu);
        ru=mul_up(xl,xl);
    } else {
        rl=nul(xl);
        Float ru1=mul_up(xl,xl);
        Float ru2=mul_up(xu,xu);
        ru=max(ru1,ru2);
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}




ValidatedFloat pow(ValidatedFloat const& x, Int n) {
    if(n<0) { return pow(rec(x),Nat(-n)); }
    else return pow(x,Nat(n));
}

ValidatedFloat pow(ValidatedFloat const& x, Nat m)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    ValidatedFloat sx = x;
    if(m%2==0) { sx=abs(x); }
    Float rl=pow_down(sx.lower_raw(),m);
    Float ru=pow_up(sx.upper_raw(),m);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



ValidatedFloat sqrt(ValidatedFloat const& x)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    Float::set_rounding_downward();
    Float rl=sqrt_rnd(x.lower_raw());
    Float::set_rounding_upward();
    Float ru=sqrt_rnd(x.upper_raw());
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat exp(ValidatedFloat const& x)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    Float::set_rounding_downward();
    Float rl=exp_rnd(x.lower_raw());
    Float::set_rounding_upward();
    Float ru=exp_rnd(x.upper_raw());
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat log(ValidatedFloat const& x)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    Float::set_rounding_downward();
    Float rl=log_rnd(x.lower_raw());
    Float::set_rounding_upward();
    Float ru=log_rnd(x.upper_raw());
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat pi_val() { return ValidatedFloat(pi_down(),pi_up()); }


ValidatedFloat sin(ValidatedFloat const& x)
{
    return cos(x-half(pi_val()));
}

ValidatedFloat cos(ValidatedFloat const& x)
{
    ARIADNE_ASSERT(x.lower_raw()<=x.upper_raw());
    Float::RoundingModeType rnd = Float::get_rounding_mode();

    static const ExactFloat two(2);

    if(x.error().raw()>2*pi_down()) { return ValidatedFloat(-1.0,+1.0); }

    Float n=floor(x.lower_raw()/(2*pi_approx())+0.5);
    ValidatedFloat y=x-two*ExactFloat(n)*pi_val();

    ARIADNE_ASSERT(y.lower_raw()<=pi_up());
    ARIADNE_ASSERT(y.upper_raw()>=-pi_up());

    Float rl,ru;
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

    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat tan(ValidatedFloat const& x) {
    return mul(sin(x),rec(cos(x)));
}

ValidatedFloat asin(ValidatedFloat const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat acos(ValidatedFloat const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat atan(ValidatedFloat const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}


ValidatedFloat::ValidatedFloat(const Dyadic& b) : ValidatedFloat(b.operator Rational()) { }

ValidatedFloat::ValidatedFloat(const Decimal& d) : ValidatedFloat(d.operator Rational()) { }


ValidatedFloat::ValidatedFloat(const Integer& z) : ValidatedFloat(Rational(z)) {
}

ValidatedFloat::ValidatedFloat(const Rational& q) : ValidatedFloat(q,q) {
}

ValidatedFloat::ValidatedFloat(const Rational& ql, const Rational& qu) : l(ql.get_d()), u(qu.get_d())  {
    static const double min_dbl=std::numeric_limits<double>::min();
    Float::RoundingModeType rounding_mode=Float::get_rounding_mode();
    Float::set_rounding_downward();
    while(Rational(l)>ql) { l=sub_rnd(l,min_dbl); }
    Float::set_rounding_upward();
    while(Rational(u)<qu) { u=add_rnd(u,min_dbl); }
    Float::set_rounding_mode(rounding_mode);
}

ValidatedFloat ExactFloat::pm(ErrorFloat e) const {
    ExactFloat const& v=*this; return ValidatedFloat(v-e,v+e);
}

OutputStream&
operator<<(OutputStream& os, const ValidatedFloat& ivl)
{
    //if(ivl.lower_raw()==ivl.upper_raw()) { return os << "{" << std::setprecision(ValidatedFloat::output_precision) << ivl.lower_raw().get_d() << ; }
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    os << '{';
    Float::set_rounding_downward();
    os << std::showpoint << std::setprecision(ValidatedFloat::output_precision) << ivl.lower().get_d();
    os << ':';
    Float::set_rounding_upward();
    os << std::showpoint << std::setprecision(ValidatedFloat::output_precision) << ivl.upper().get_d();
    Float::set_rounding_mode(rnd);
    os << '}';
    return os;

}

InputStream&
operator>>(InputStream& is, ValidatedFloat& x)
{
    Float l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    x.l=l; x.u=u;
    return is;
}


UpperFloat::UpperFloat(Number<Upper> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}


OutputStream& operator<<(OutputStream& os, UpperFloat const& x) {
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float::set_rounding_upward();
    os << std::showpoint << std::setprecision(ValidatedFloat::output_precision) << x.raw();
    Float::set_rounding_mode(rnd);
    return os;
}

UpperFloat sqr(UpperFloat const& x) {
    ARIADNE_PRECONDITION(x.raw()>=0.0);
    return UpperFloat(mul_up(x.raw(),x.raw()));
}

UpperFloat mul(UpperFloat const& x1, UpperFloat const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0.0 && x2.raw() >= 0.0);
    return UpperFloat(mul_up(x1.raw(),x2.raw()));
}

UpperFloat div(UpperFloat const& x1, LowerFloat const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0.0 && x2.raw() > 0.0);
    return UpperFloat(div_up(x1.raw(),x2.raw()));
}

UpperFloat pow(UpperFloat const& x, Nat n) {
    ARIADNE_PRECONDITION(x.raw()>=0.0);
    return UpperFloat(pow_up(x.raw(),n));
}


UpperFloat rec(LowerFloat const& x) {
    return UpperFloat(rec_up(x.raw()));
}

UpperFloat sqrt(UpperFloat const& x) {
    return UpperFloat(sqrt_up(x.raw()));
}

UpperFloat exp(UpperFloat const& x) {
    return UpperFloat(exp_up(x.raw()));
}

UpperFloat log(UpperFloat const& x) {
    return UpperFloat(log_up(x.raw()));
}

template<> Int integer_cast(UpperFloat const& x) { return static_cast<Int>(x.u.get_d()); }
template<> Nat integer_cast(UpperFloat const& x) { return static_cast<Nat>(x.u.get_d()); }




PositiveUpperFloat operator+(PositiveUpperFloat const& x1, PositiveUpperFloat const& x2) {
    return PositiveUpperFloat(add_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat operator-(PositiveUpperFloat const& x1, LowerFloat const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=x2.raw());
    return PositiveUpperFloat(sub_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat operator*(PositiveUpperFloat const& x1, PositiveUpperFloat const& x2) {
    return PositiveUpperFloat(mul_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat operator/(PositiveUpperFloat const& x1, LowerFloat const& x2) {
    ARIADNE_PRECONDITION(x2.raw()>=0);
    return PositiveUpperFloat(div_up(x1.raw(),x2.raw()));
}

PositiveUpperFloat pow(PositiveUpperFloat const& x, Nat m) {
    return PositiveUpperFloat(pow_up(x.raw(),m));
}

PositiveUpperFloat half(PositiveUpperFloat const& x) {
    return PositiveUpperFloat(half(x.raw()));
}


LowerFloat::LowerFloat(Number<Lower> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

LowerFloat mul(LowerFloat const& x1, LowerFloat const& x2) {
    ARIADNE_PRECONDITION(x1.raw()>=0 && x2.raw()>=0);
    return LowerFloat(mul_down(x1.raw(),x2.raw()));
}

LowerFloat div(LowerFloat const& x1, UpperFloat const& x2) {
    //ARIADNE_PRECONDITION_MSG(x1.raw()>=0 && x2.raw()>=0,"x1="<<x1<<", x2="<<x2);
    return LowerFloat(div_down(x1.raw(),x2.raw()));
}

LowerFloat rec(UpperFloat const& x) {
    return LowerFloat(rec_down(x.raw()));
}

LowerFloat sqrt(LowerFloat const& x) {
    return LowerFloat(sqrt_down(x.raw()));
}

LowerFloat exp(LowerFloat const& x) {
    return LowerFloat(exp_down(x.raw()));
}

LowerFloat log(LowerFloat const& x) {
    return LowerFloat(log_down(x.raw()));
}

OutputStream& operator<<(OutputStream& os, LowerFloat const& x) {
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    Float::set_rounding_downward();
    os << std::showpoint << std::setprecision(ValidatedFloat::output_precision) << x.raw();
    Float::set_rounding_mode(rnd);
    return os;
}

template<> Int integer_cast(LowerFloat const& x) { return static_cast<Int>(x.l.get_d()); }
template<> Nat integer_cast(LowerFloat const& x) { return static_cast<Nat>(x.l.get_d()); }



//ExactFloat inf = ExactFloat(std::numeric_limits< double >::infinity());
ApproximateFloat::ApproximateFloat(Dyadic const& b) : ApproximateFloat(b.operator Rational()) { }
ApproximateFloat::ApproximateFloat(Decimal const& d) : ApproximateFloat(d.operator Rational()) { }

ApproximateFloat::ApproximateFloat(Number<Approximate> const& x) { ARIADNE_NOT_IMPLEMENTED; }

ExactFloat::operator Rational() const {
    return Rational(this->get_d());
}

ApproximateFloat::ApproximateFloat(Rational const& q) : ApproximateFloat(q.get_d()) {
}

OutputStream& operator<<(OutputStream& os, ApproximateFloat const& x) {
    return os << std::showpoint << std::setprecision(ApproximateFloat::output_precision) << x.raw();
}

template<> Int integer_cast(ApproximateFloat const& x) { return static_cast<Int>(x.a.get_d()); }
template<> Nat integer_cast(ApproximateFloat const& x) { return static_cast<Nat>(x.a.get_d()); }




ApproximateFloat create_float(Number<Approximate> x) { return ApproximateFloat(x); }
LowerFloat create_float(Number<Lower> x) { return LowerFloat(x); }
UpperFloat create_float(Number<Upper> x) { return UpperFloat(x); }
ValidatedFloat create_float(Number<Validated> x) { return ValidatedFloat(x); }
ValidatedFloat create_float(Number<Effective> x) { return ValidatedFloat(x); }
ValidatedFloat create_float(Number<Exact> x) { return ValidatedFloat(x); }
ValidatedFloat create_float(Real r) { return ValidatedFloat(r); }
ValidatedFloat create_float(Rational q) { return ValidatedFloat(q); }
ExactFloat create_float(Integer z) { return ExactFloat(z); }



} // namespace Ariadne
