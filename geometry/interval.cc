/***************************************************************************
 *            interval.cc
 *
 *  Copyright 2008-10  Pieter Collins
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
#include "config.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include "utility/container.h"

#include "utility/macros.h"
#include "utility/exceptions.h"
#include "numeric/integer.h"
#include "numeric/float.h"
#include "numeric/rational.h"
#include "numeric/float-exact.h"
#include "numeric/decimal.h"
#include "numeric/dyadic.h"
#include "geometry/interval.h"


namespace Ariadne {

namespace {
static const double _pi_up  =3.1415926535897936;
static const double _pi_down=3.1415926535897931;
static const double _pi_approx=3.1415926535897931;
inline double _add_down(volatile double x, volatile double y) { set_rounding_downward(); return x+y; }
inline double _add_up(volatile double x, volatile double y) { set_rounding_upward(); return x+y; }
inline double _sub_down(volatile double x, volatile double y) { set_rounding_downward(); return x-y; }
inline double _sub_up(volatile double x, volatile double y) { set_rounding_upward(); return x-y; }
inline double _mul_down(volatile double x, volatile double y) { set_rounding_downward(); return x*y; }
inline double _mul_up(volatile double x, volatile double y) { set_rounding_upward(); return x*y; }
inline double _div_down(volatile double x, volatile double y) { set_rounding_downward(); return x/y; }
inline double _div_up(volatile double x, volatile double y) { set_rounding_upward(); return x/y; }
inline Float cos_down(Float x) { set_rounding_downward(); Float y=cos_rnd(x); return y; }
inline Float cos_up(Float x) { set_rounding_upward(); Float y=cos_rnd(x); return y; }
} // namespace



const UpperInterval pi_ivl=ExactInterval(pi_down,pi_up);

uint ExactInterval::output_precision = 6;

ExactInterval widen(ExactInterval x)
{
    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower().value());
    const double& xu=internal_cast<const double&>(x.upper().value());
    const double m=std::numeric_limits<float>::min();
    set_rounding_upward();
    volatile double wu=xu+m;
    volatile double mwl=-xl+m;
    volatile double wl=-mwl;
    set_rounding_mode(rm);
    assert(wl<xl); assert(wu>xu);
    return ExactInterval(wl,wu);
}

ExactInterval narrow(ExactInterval x)
{
    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower().value());
    const double& xu=internal_cast<const double&>(x.upper().value());
    const double m=std::numeric_limits<float>::min();
    set_rounding_upward();
    volatile double mnu=-xu+m;
    volatile double nu=-mnu;
    volatile double nl=xl+m;
    set_rounding_mode(rm);
    assert(xl<nl); assert(nu<xu);
    return ExactInterval(nl,nu);
}

ExactInterval trunc(ExactInterval x)
{

    rounding_mode_t rm=get_rounding_mode();
    const double& xl=internal_cast<const double&>(x.lower().value());
    const double& xu=internal_cast<const double&>(x.upper().value());
    // Use machine epsilon instead of minimum to move away from zero
    const float fm=std::numeric_limits<float>::epsilon();
    volatile float tu=xu;
    if(tu<xu) { set_rounding_upward(); tu+=fm; }
    volatile float tl=xl;
    if(tl>xl) { set_rounding_downward(); tl-=fm; }
    set_rounding_mode(rm);
    assert(tl<=xl); assert(tu>=xu);
    return ExactInterval(double(tl),double(tu));
}

ExactInterval trunc(ExactInterval x, uint n)
{
    ExactInterval e=ExactInterval(std::pow(2.0,52-(int)n));
    UpperInterval y=x+e;
    UpperInterval r=y-e;
    return ExactInterval(r.lower_raw(),r.upper_raw());
}

UpperInterval rec(UpperInterval i)
{
    volatile double& il=internal_cast<volatile double&>(i.lower().value());
    volatile double& iu=internal_cast<volatile double&>(i.upper().value());
    volatile double rl,ru;
    if(il>0 || iu<0) {
        rounding_mode_t rnd=get_rounding_mode();
        rl=_div_down(1.0,iu);
        ru=_div_up(1.0,il);
        set_rounding_mode(rnd);
    } else {
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
        ARIADNE_THROW(DivideByZeroException,"UpperInterval rec(UpperInterval ivl)","ivl="<<i);
    }
    return UpperInterval(rl,ru);
}


UpperInterval mul(UpperInterval i1, UpperInterval i2)
{
    volatile double& i1l=internal_cast<volatile double&>(i1.lower().value());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper().value());
    volatile double& i2l=internal_cast<volatile double&>(i2.lower().value());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper().value());
    volatile double rl,ru;
    rounding_mode_t rnd=get_rounding_mode();
    if(i1l>=0) {
        if(i2l>=0) {
            rl=_mul_down(i1l,i2l); ru=_mul_up(i1u,i2u);
        } else if(i2u<=0) {
            rl=_mul_down(i1u,i2l); ru=_mul_up(i1l,i2u);
        } else {
            rl=_mul_down(i1u,i2l); ru=_mul_up(i1u,i2u);
        }
    }
    else if(i1u<=0) {
        if(i2l>=0) {
            rl=_mul_down(i1l,i2u); ru=_mul_up(i1u,i2l);
        } else if(i2u<=0) {
            rl=_mul_down(i1u,i2u); ru=_mul_up(i1l,i2l);
        } else {
            rl=_mul_down(i1l,i2u); ru=_mul_up(i1l,i2l);
        }
    } else {
        if(i2l>=0) {
            rl=_mul_down(i1l,i2u); ru=_mul_up(i1u,i2u);
        } else if(i2u<=0) {
            rl=_mul_down(i1u,i2l); ru=_mul_up(i1l,i2l);
        } else {
            set_rounding_mode(downward);
            rl=std::min(i1u*i2l,i1l*i2u);
            set_rounding_mode(upward);
            ru=std::max(i1l*i2l,i1u*i2u);
        }
    }
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}


UpperInterval mul(UpperInterval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower().value());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper().value());
    volatile double& x2v=internal_cast<volatile double&>(x2);
    volatile double rl,ru;
    if(x2>=0) {
        rl=_mul_down(i1l,x2v); ru=_mul_up(i1u,x2v);
    } else {
        rl=_mul_down(i1u,x2v); ru=_mul_up(i1l,x2v);
    }
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}


UpperInterval mul(Float x1, UpperInterval i2)
{
    return mul(i2,x1);
}


UpperInterval div(UpperInterval i1, UpperInterval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower().value());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper().value());
    volatile double& i2l=internal_cast<volatile double&>(i2.lower().value());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper().value());
    volatile double rl,ru;
    if(i2l>0) {
        if(i1l>=0) {
            rl=_div_down(i1l,i2u); ru=_div_up(i1u,i2l);
        } else if(i1u<=0) {
            rl=_div_down(i1l,i2l); ru=_div_up(i1u,i2u);
        } else {
            rl=_div_down(i1l,i2l); ru=_div_up(i1u,i2l);
        }
    }
    else if(i2u<0) {
        if(i1l>=0) {
            rl=_div_down(i1u,i2u); ru=_div_up(i1l,i2l);
        } else if(i1u<=0) {
            rl=_div_down(i1u,i2l); ru=_div_up(i1l,i2u);
        } else {
            rl=_div_down(i1u,i2u); ru=_div_up(i1l,i2u);
        }
    }
    else {
        // ARIADNE_THROW(DivideByZeroException,"UpperInterval div(UpperInterval ivl1, UpperInterval ivl2)","ivl1="<<i1<<", ivl2="<<i2);
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    }
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}



UpperInterval div(UpperInterval i1, Float x2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=internal_cast<volatile double&>(i1.lower().value());
    volatile double& i1u=internal_cast<volatile double&>(i1.upper().value());
    volatile double& x2v=internal_cast<volatile double&>(x2);
    volatile double rl,ru;
    if(x2v>0) {
        rl=_div_down(i1l,x2v); ru=_div_up(i1u,x2v);
    } else if(x2v<0) {
        rl=_div_down(i1u,x2v); ru=_div_up(i1l,x2v);
    } else {
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    }
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}


UpperInterval div(Float x1, UpperInterval i2)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=internal_cast<volatile double&>(x1);
    volatile double& i2l=internal_cast<volatile double&>(i2.lower().value());
    volatile double& i2u=internal_cast<volatile double&>(i2.upper().value());
    volatile double rl,ru;
    if(i2l<=0 && i2u>=0) {
        ARIADNE_THROW(DivideByZeroException,"UpperInterval div(Float x1, UpperInterval ivl2)","x1="<<x1<<", ivl2="<<i2);
        rl=-std::numeric_limits<double>::infinity();
        ru=+std::numeric_limits<double>::infinity();
    } else if(x1v>=0) {
        rl=_div_down(x1v,i2u); ru=_div_up(x1v,i2l);
    } else {
        rl=_div_down(x1v,i2l); ru=_div_up(x1v,i2u);
    }
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

UpperInterval sqr(UpperInterval i)
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& il=internal_cast<volatile double&>(i.lower().value());
    volatile double& iu=internal_cast<volatile double&>(i.upper().value());
    volatile double rl,ru;
    if(il>0.0) {
        set_rounding_mode(downward);
        rl=il*il;
        set_rounding_mode(upward);
        ru=iu*iu;
    } else if(iu<0.0) {
        set_rounding_mode(downward);
        rl=iu*iu;
        set_rounding_mode(upward);
        ru=il*il;
    } else {
        rl=0.0;
        set_rounding_mode(upward);
        volatile double ru1=il*il;
        volatile double ru2=iu*iu;
        ru=max(ru1,ru2);
    }
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}




UpperInterval pow(UpperInterval i, int n)
{
    if(n<0) { return pow(rec(i),uint(-n)); }
    else return pow(i,uint(n));
}

UpperInterval pow(UpperInterval i, uint m)
{
    rounding_mode_t rnd = get_rounding_mode();
    const UpperInterval& nvi=i;
    if(m%2==0) { i=abs(nvi); }
    set_rounding_mode(downward);
    Float rl=pow_rnd(i.lower().value(),m);
    set_rounding_mode(upward);
    Float ru=pow_rnd(i.upper().value(),m);
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}



UpperInterval sqrt(UpperInterval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    Float rl=sqrt_rnd(i.lower().value());
    set_rounding_upward();
    Float ru=sqrt_rnd(i.upper().value());
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

UpperInterval exp(UpperInterval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    Float rl=exp_rnd(i.lower().value());
    set_rounding_upward();
    Float ru=exp_rnd(i.upper().value());
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

UpperInterval log(UpperInterval i)
{
    rounding_mode_t rnd = get_rounding_mode();
    set_rounding_downward();
    Float rl=log_rnd(i.lower().value());
    set_rounding_upward();
    Float ru=log_rnd(i.upper().value());
    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}



UpperInterval sin(UpperInterval i)
{
    return cos(i-pi_ivl/2);
}

UpperInterval cos(UpperInterval i)
{
    ARIADNE_ASSERT(i.lower().value()<=i.upper().value());
    rounding_mode_t rnd = get_rounding_mode();

    if(i.radius().value().raw()>2*pi_down) { return UpperInterval(-1.0,+1.0); }

    Float n=floor(i.lower().value()/(2*_pi_approx)+0.5);
    i=i-(2*n)*pi_ivl;

    ARIADNE_ASSERT(i.lower().value()<=pi_up);
    ARIADNE_ASSERT(i.upper().value()>=-pi_up);

    Float rl,ru;
    if(i.lower().value()<=-pi_down) {
        if(i.upper().value()<=0.0) { rl=-1.0; ru=cos_up(i.upper().value()); }
        else { rl=-1.0; ru=+1.0; }
    } else if(i.lower().value()<=0.0) {
        if(i.upper().value()<=0.0) { rl=cos_down(i.lower().value()); ru=cos_up(i.upper().value()); }
        else if(i.upper().value()<=pi_down) { rl=cos_down(max(-i.lower().value(),i.upper().value())); ru=+1.0; }
        else { rl=-1.0; ru=+1.0; }
    } else if(i.lower().value()<=pi_up) {
        if(i.upper().value()<=pi_down) { rl=cos_down(i.upper().value()); ru=cos_up(i.lower().value()); }
        else if(i.upper().value()<=2*pi_down) { rl=-1.0; ru=cos_up(min(i.lower().value(),sub_down(2*pi_down,i.upper().value()))); }
        else { rl=-1.0; ru=+1.0; }
    } else {
        assert(false);
    }

    set_rounding_mode(rnd);
    return UpperInterval(rl,ru);
}

UpperInterval tan(UpperInterval i)
{
    return mul(sin(i),rec(cos(i)));
}

UpperInterval asin(UpperInterval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

UpperInterval acos(UpperInterval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}

UpperInterval atan(UpperInterval i)
{
    ARIADNE_NOT_IMPLEMENTED;
}


ExactInterval::ExactInterval(const Dyadic& b) : ExactInterval(b.operator Rational()) { }

ExactInterval::ExactInterval(const Decimal& d) : ExactInterval(d.operator Rational()) { }


#ifdef HAVE_GMPXX_H

ExactInterval::ExactInterval(const Integer& z) : ExactInterval(Rational(z)) {
}

ExactInterval::ExactInterval(const Rational& q) : ExactInterval(q,q) {
}

ExactInterval::ExactInterval(const Rational& ql, const Rational& qu) : l(ql.get_d()), u(qu.get_d())  {
    static const double min_dbl=std::numeric_limits<double>::min();
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(downward);
    while(l.get_d()>static_cast<mpq_class const&>(ql)) { l=sub_rnd(l,min_dbl); }
    set_rounding_mode(upward);
    while(u.get_d()<static_cast<mpq_class const&>(qu)) { u=add_rnd(u,min_dbl); }
    set_rounding_mode(rounding_mode);
}

ExactInterval::ExactInterval(const Real& lower, const Real& upper)
    : l(ValidatedFloat(lower).lower().value()), u(ValidatedFloat(upper).upper().value()) {
}

ExactInterval& ExactInterval::operator=(const Rational& q) {
    return *this = ExactInterval(q);
}


#endif // HAVE_GMPXX_H


std::ostream&
operator<<(std::ostream& os, const ExactInterval& ivl)
{
    //if(ivl.lower()==ivl.upper().value()) { return os << "{" << std::setprecision(ExactInterval::output_precision) << ivl.lower().value().get_d() << ; }
    rounding_mode_t rnd=get_rounding_mode();
    os << '{';
    set_rounding_downward();
    os << std::showpoint << std::setprecision(ExactInterval::output_precision) << ivl.lower().value().get_d();
    os << ':';
    set_rounding_upward();
    os << std::showpoint << std::setprecision(ExactInterval::output_precision) << ivl.upper().value().get_d();
    set_rounding_mode(rnd);
    os << '}';
    return os;

}

/*
std::ostream&
operator<<(std::ostream& os, const ExactInterval& ivl)
{
    return os << '[' << ivl.l << ':' << ivl.u << ']';
}
*/

/*
std::ostream&
operator<<(std::ostream& os, const ExactInterval& ivl)
{
    if(ivl.lower().value()==ivl.upper().value()) {
        return os << std::setprecision(18) << ivl.lower().value();
    }

    std::stringstream iss,uss;
    iss << std::setprecision(18) << ivl.lower().value();
    uss << std::setprecision(18) << ivl.upper().value();

    std::string lstr,ustr;
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
    uint i;
    for(i=0; (i<std::min(lstr.size(),ustr.size()) && lstr[i]==ustr[i]); ++i) {
        os << lstr[i];
    }

    os << "[";
    if(i==lstr.size()) {
        os << "0";
    }
    for(uint li=i; li != lstr.size(); ++li) {
        os << lstr[li];
    }
    os << ":";
    if(i==ustr.size()) {
        os << "0";
    }
    for(uint ui=i; ui != ustr.size(); ++ui) {
        os << ustr[ui];
    }
    os << "]";
    return os;

}
*/
std::istream&
operator>>(std::istream& is, ExactInterval& ivl)
{
    Float l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(l,u);
    return is;
}



} // namespace Ariadne

