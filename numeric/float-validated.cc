/***************************************************************************
 *            float-validated.cc
 *
 *  Copyright 2008-14  Pieter Collins
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

#include "utility/standard.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include "utility/container.h"



#include "config.h"
#include "utility/typedefs.h"
#include "utility/macros.h"
#include "utility/exceptions.h"
#include "numeric/integer.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "numeric/rational.h"
#include "numeric/real.h"
#include "numeric/number.h"
#include "numeric/float.h"
#include "numeric/float-exact.h"
#include "numeric/float-approximate.h"

#include "numeric/float-validated.h"

namespace Ariadne {

namespace {
inline Float cos_down(Float x) { Float::set_rounding_downward(); Float y=cos_rnd(x); set_default_rounding(); return y; }
inline Float cos_up(Float x) { Float::set_rounding_upward(); Float y=cos_rnd(x); set_default_rounding(); return y; }
} // namespace
const ValidatedFloat pi_val=ValidatedFloat(pi_down(Precision64()),pi_up(Precision64()));

ValidatedFloat::ValidatedFloat(Number<Validated> const& x) {
    ARIADNE_NOT_IMPLEMENTED;
}

Nat ValidatedFloat::output_precision = 6;

ValidatedFloat::ValidatedFloat(const ExactFloat& x) : ValidatedFloat(x.raw(),x.raw()) {
}

ValidatedFloat widen(ValidatedFloat x)
{
    Float::RoundingModeType rm=Float::get_rounding_mode();
    const Float& xl=x.lower_raw();
    const Float& xu=x.upper_raw();
    const Float m=std::numeric_limits<float>::min();
    Float::set_rounding_upward();
    Float wu=xu+m;
    Float mwl=-xl+m;
    Float wl=-mwl;
    Float::set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return ValidatedFloat(wl,wu);
}

ValidatedFloat narrow(ValidatedFloat x)
{
    Float::RoundingModeType rm=Float::get_rounding_mode();
    const Float& xl=x.lower_raw();
    const Float& xu=x.upper_raw();
    const Float m=std::numeric_limits<float>::min();
    Float::set_rounding_upward();
    Float mnu=-xu+m;
    Float nu=-mnu;
    Float nl=xl+m;
    Float::set_rounding_mode(rm);
    assert(xl<nl); assert(nu<xu);
    return ValidatedFloat(nl,nu);
}

ValidatedFloat trunc(ValidatedFloat x)
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

ValidatedFloat trunc(ValidatedFloat x, Nat n)
{
    ValidatedFloat e=ValidatedFloat(std::pow(2.0,52-(Int)n));
    ValidatedFloat y=x+e;
    return y-e;
}

ValidatedFloat rec(ValidatedFloat i)
{
    const Float& il=i.lower_raw();
    const Float& iu=i.upper_raw();
    Float rl,ru;
    if(il>0 || iu<0) {
        Float::RoundingModeType rnd=Float::get_rounding_mode();
        rl=div_down(1.0,iu);
        ru=div_up(1.0,il);
        Float::set_rounding_mode(rnd);
    } else {
        rl=-inf;
        ru=+inf;
        ARIADNE_THROW(DivideByZeroException,"ValidatedFloat rec(ValidatedFloat ivl)","ivl="<<i);
    }
    return ValidatedFloat(rl,ru);
}


ValidatedFloat mul(ValidatedFloat i1, ValidatedFloat i2)
{
    const Float& i1l=i1.lower_raw();
    const Float& i1u=i1.upper_raw();
    const Float& i2l=i2.lower_raw();
    const Float& i2u=i2.upper_raw();
    Float rl,ru;
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    if(i1l>=0) {
        if(i2l>=0) {
            rl=mul_down(i1l,i2l); ru=mul_up(i1u,i2u);
        } else if(i2u<=0) {
            rl=mul_down(i1u,i2l); ru=mul_up(i1l,i2u);
        } else {
            rl=mul_down(i1u,i2l); ru=mul_up(i1u,i2u);
        }
    }
    else if(i1u<=0) {
        if(i2l>=0) {
            rl=mul_down(i1l,i2u); ru=mul_up(i1u,i2l);
        } else if(i2u<=0) {
            rl=mul_down(i1u,i2u); ru=mul_up(i1l,i2l);
        } else {
            rl=mul_down(i1l,i2u); ru=mul_up(i1l,i2l);
        }
    } else {
        if(i2l>=0) {
            rl=mul_down(i1l,i2u); ru=mul_up(i1u,i2u);
        } else if(i2u<=0) {
            rl=mul_down(i1u,i2l); ru=mul_up(i1l,i2l);
        } else {
            Float::set_rounding_downward();
            rl=min(i1u*i2l,i1l*i2u);
            Float::set_rounding_upward();
            ru=max(i1l*i2l,i1u*i2u);
        }
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat mul(ValidatedFloat i1, ExactFloat x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& i1l=i1.lower_raw();
    const Float& i1u=i1.upper_raw();
    const Float& x2v=x2.raw();
    Float rl,ru;
    if(x2>=0) {
        rl=mul_down(i1l,x2v); ru=mul_up(i1u,x2v);
    } else {
        rl=mul_down(i1u,x2v); ru=mul_up(i1l,x2v);
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat mul(ExactFloat x1, ValidatedFloat i2)
{
    return mul(i2,x1);
}


ValidatedFloat div(ValidatedFloat i1, ValidatedFloat i2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& i1l=i1.lower_raw();
    const Float& i1u=i1.upper_raw();
    const Float& i2l=i2.lower_raw();
    const Float& i2u=i2.upper_raw();
    Float rl,ru;
    if(i2l>0) {
        if(i1l>=0) {
            rl=div_down(i1l,i2u); ru=div_up(i1u,i2l);
        } else if(i1u<=0) {
            rl=div_down(i1l,i2l); ru=div_up(i1u,i2u);
        } else {
            rl=div_down(i1l,i2l); ru=div_up(i1u,i2l);
        }
    }
    else if(i2u<0) {
        if(i1l>=0) {
            rl=div_down(i1u,i2u); ru=div_up(i1l,i2l);
        } else if(i1u<=0) {
            rl=div_down(i1u,i2l); ru=div_up(i1l,i2u);
        } else {
            rl=div_down(i1u,i2u); ru=div_up(i1l,i2u);
        }
    }
    else {
        // ARIADNE_THROW(DivideByZeroException,"ValidatedFloat div(ValidatedFloat ivl1, ValidatedFloat ivl2)","ivl1="<<i1<<", ivl2="<<i2);
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



ValidatedFloat div(ValidatedFloat i1, ExactFloat x2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& i1l=i1.lower_raw();
    const Float& i1u=i1.upper_raw();
    const Float& x2v=x2.raw();
    Float rl,ru;
    if(x2v>0) {
        rl=div_down(i1l,x2v); ru=div_up(i1u,x2v);
    } else if(x2v<0) {
        rl=div_down(i1u,x2v); ru=div_up(i1l,x2v);
    } else {
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}


ValidatedFloat div(ExactFloat x1, ValidatedFloat i2)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& x1v=x1.raw();
    const Float& i2l=i2.lower_raw();
    const Float& i2u=i2.upper_raw();
    Float rl,ru;
    if(i2l<=0 && i2u>=0) {
        ARIADNE_THROW(DivideByZeroException,"ValidatedFloat div(Float x1, ValidatedFloat ivl2)","x1="<<x1<<", ivl2="<<i2);
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    } else if(x1v>=0) {
        rl=div_down(x1v,i2u); ru=div_up(x1v,i2l);
    } else {
        rl=div_down(x1v,i2l); ru=div_up(x1v,i2u);
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat sqr(ValidatedFloat i)
{
    Float::RoundingModeType rnd=Float::get_rounding_mode();
    const Float& il=i.lower_raw();
    const Float& iu=i.upper_raw();
    Float rl,ru;
    if(il>0.0) {
        Float::set_rounding_downward();
        rl=il*il;
        Float::set_rounding_upward();
        ru=iu*iu;
    } else if(iu<0.0) {
        Float::set_rounding_downward();
        rl=iu*iu;
        Float::set_rounding_upward();
        ru=il*il;
    } else {
        rl=0.0;
        Float::set_rounding_upward();
        Float ru1=il*il;
        Float ru2=iu*iu;
        ru=max(ru1,ru2);
    }
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}




ValidatedFloat pow(ValidatedFloat i, Int n)
{
    if(n<0) { return pow(rec(i),Nat(-n)); }
    else return pow(i,Nat(n));
}

ValidatedFloat pow(ValidatedFloat i, Nat m)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    const ValidatedFloat& nvi=i;
    if(m%2==0) { i=abs(nvi); }
    Float::set_rounding_downward();
    Float rl=pow_rnd(i.lower_raw(),m);
    Float::set_rounding_upward();
    Float ru=pow_rnd(i.upper_raw(),m);
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



ValidatedFloat sqrt(ValidatedFloat i)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    Float::set_rounding_downward();
    Float rl=sqrt_rnd(i.lower_raw());
    Float::set_rounding_upward();
    Float ru=sqrt_rnd(i.upper_raw());
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat exp(ValidatedFloat i)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    Float::set_rounding_downward();
    Float rl=exp_rnd(i.lower_raw());
    Float::set_rounding_upward();
    Float ru=exp_rnd(i.upper_raw());
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat log(ValidatedFloat i)
{
    Float::RoundingModeType rnd = Float::get_rounding_mode();
    Float::set_rounding_downward();
    Float rl=log_rnd(i.lower_raw());
    Float::set_rounding_upward();
    Float ru=log_rnd(i.upper_raw());
    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}



ValidatedFloat sin(ValidatedFloat i)
{
    return cos(i-half(pi_val));
}

ValidatedFloat cos(ValidatedFloat i)
{
    ARIADNE_ASSERT(i.lower_raw()<=i.upper_raw());
    Float::RoundingModeType rnd = Float::get_rounding_mode();

    static const ExactFloat two(2);

    if(i.radius().raw()>2*pi_down()) { return ValidatedFloat(-1.0,+1.0); }

    Float n=floor(i.lower_raw()/(2*pi_approx())+0.5);
    i=i-two*ExactFloat(n)*pi_val;

    ARIADNE_ASSERT(i.lower_raw()<=pi_up());
    ARIADNE_ASSERT(i.upper_raw()>=-pi_up());

    Float rl,ru;
    if(i.lower_raw()<=-pi_down()) {
        if(i.upper_raw()<=0.0) { rl=-1.0; ru=cos_up(i.upper_raw()); }
        else { rl=-1.0; ru=+1.0; }
    } else if(i.lower_raw()<=0.0) {
        if(i.upper_raw()<=0.0) { rl=cos_down(i.lower_raw()); ru=cos_up(i.upper_raw()); }
        else if(i.upper_raw()<=pi_down()) { rl=cos_down(max(-i.lower_raw(),i.upper_raw())); ru=+1.0; }
        else { rl=-1.0; ru=+1.0; }
    } else if(i.lower_raw()<=pi_up()) {
        if(i.upper_raw()<=pi_down()) { rl=cos_down(i.upper_raw()); ru=cos_up(i.lower_raw()); }
        else if(i.upper_raw()<=2*pi_down()) { rl=-1.0; ru=cos_up(min(i.lower_raw(),sub_down(2*pi_down(),i.upper_raw()))); }
        else { rl=-1.0; ru=+1.0; }
    } else {
        assert(false);
    }

    Float::set_rounding_mode(rnd);
    return ValidatedFloat(rl,ru);
}

ValidatedFloat tan(ValidatedFloat i)
{
    return mul(sin(i),rec(cos(i)));
}

ValidatedFloat asin(ValidatedFloat i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat acos(ValidatedFloat i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedFloat atan(ValidatedFloat i)
{
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

/*
OutputStream&
operator<<(OutputStream& os, const ValidatedFloat& ivl)
{
    return os << '[' << ivl.l << ':' << ivl.u << ']';
}
*/

/*
OutputStream&
operator<<(OutputStream& os, const ValidatedFloat& ivl)
{
    if(ivl.lower_raw()==ivl.upper_raw()) {
        return os << std::setprecision(18) << ivl.lower_raw();
    }

    StringStream iss,uss;
    iss << std::setprecision(18) << ivl.lower_raw();
    uss << std::setprecision(18) << ivl.upper_raw();

    StringType lstr,ustr;
    iss >> lstr; uss >> ustr;

    // Test if one endpoint is an integer and the other is not
    // If this is the case, append ".0" to to integer value
    if( (lstr.find('.')==lstr.size()) xor (ustr.find('.')==lstr.size()) ) {
        if(lstr.find('.')==lstr.size()) {
            lstr+=".0";
        } else {
            ustr+=".0";
        }
    }

    // Write common head
    Nat i;
    for(i=0; (i<std::min(lstr.size(),ustr.size()) && lstr[i]==ustr[i]); ++i) {
        os << lstr[i];
    }

    os << "[";
    if(i==lstr.size()) {
        os << "0";
    }
    for(Nat li=i; li != lstr.size(); ++li) {
        os << lstr[li];
    }
    os << ":";
    if(i==ustr.size()) {
        os << "0";
    }
    for(Nat ui=i; ui != ustr.size(); ++ui) {
        os << ustr[ui];
    }
    os << "]";
    return os;

}
*/
InputStream&
operator>>(InputStream& is, ValidatedFloat& ivl)
{
    double l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(l,u);
    return is;
}


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

