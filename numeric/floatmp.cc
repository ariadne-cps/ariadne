/***************************************************************************
 *            floatmp.cc
 *
 *  Copyright 2013--17  Pieter Collins
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
 *  You should have received _a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file floatmp.cc
 *  \brief
 */



#include "utility/module.h"
#include "logical.h"
#include "floatmp.h"
#include "float64.h"
#include "dyadic.h"
#include "rational.h"

namespace Ariadne {

// Rule for combining mixed precision
inline PrecisionMP cmb(PrecisionMP pr1, PrecisionMP pr2) { return min(pr1,pr2); }
inline PrecisionMP cmb(PrecisionMP pr1, PrecisionMP pr2, PrecisionMP pr3) { return cmb(cmb(pr1,pr2),pr3); }

FloatMP::~FloatMP() {
    mpfr_clear(_mpfr);
}

FloatMP::FloatMP() {
    mpfr_init_set_si(_mpfr,0l,get_rounding_mode());
}

FloatMP::FloatMP(NoInit) {
    mpfr_init(_mpfr);
}

/*
FloatMP::FloatMP(const mpfr_t x) : FloatMP(NoInit()) {
    mpfr_set_prec(this->_mpfr,mpfr_get_prec(x));
    mpfr_set(this->_mpfr,x,MPFR_RNDN);
}
*/

FloatMP::FloatMP(double d) : FloatMP(d,get_default_precision()) {
}

FloatMP::FloatMP(double d, PrecisionMP pr) : FloatMP(d,MPFR_RNDN,pr) {
    ARIADNE_ASSERT(d==this->get_d());
}

FloatMP::FloatMP(Float64 x, PrecisionMP pr) : FloatMP(x.get_d(),pr) {
}

FloatMP::FloatMP(PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_si(_mpfr,0l,get_rounding_mode());
}

FloatMP::FloatMP(PrecisionMP pr, NoInit) {
    mpfr_init2(_mpfr,pr);
}

FloatMP::FloatMP(Int32 n, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_si(_mpfr,n.get_si(),get_rounding_mode());
}

FloatMP::FloatMP(Dyadic const& w, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_f(_mpfr,w.get_mpf(),get_rounding_mode());
    ARIADNE_ASSERT(Dyadic(*this)==w);
}

FloatMP::FloatMP(double d, RoundingModeType rnd, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_d(_mpfr,d,rnd);
}

FloatMP::FloatMP(Float64 x, RoundingModeType rnd, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_d(_mpfr,x.get_d(),rnd);
}

FloatMP::FloatMP(Integer const& z, RoundingModeType rnd, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_z(_mpfr,z.get_mpz(),rnd);
}

FloatMP::FloatMP(Dyadic const& w, RoundingModeType rnd, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_f(_mpfr,w.get_mpf(),rnd);
}

FloatMP::FloatMP(Rational const& q, RoundingModeType rnd, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_q(_mpfr,q.get_mpq(),rnd);
}

FloatMP::FloatMP(FloatMP const& x, RoundingModeType rnd, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set(_mpfr,x._mpfr,rnd);
}

FloatMP::FloatMP(const FloatMP& x) {
    mpfr_init2(_mpfr,mpfr_get_prec(x._mpfr));
    mpfr_set(_mpfr,x._mpfr,get_rounding_mode());
}

FloatMP::FloatMP(FloatMP&& x) {
    mpfr_init(_mpfr);
    mpfr_swap(_mpfr,x._mpfr);
}

FloatMP& FloatMP::operator=(const FloatMP& x) {
    // TODO: Decide whether equality changes precision
    // NOTE: mpfr_set_prec clears a number, even if precision does not change
    //   Hence we should check for self-assignment explicitly, and/or only
    //   change precision when necessary
    if(this->precision()!=x.precision()) {
        mpfr_set_prec(_mpfr,x.precision());
    }
    mpfr_set(_mpfr,x._mpfr,get_rounding_mode());
    return *this;
}

FloatMP& FloatMP::operator=(FloatMP&& x) {
        mpfr_swap(_mpfr,x._mpfr);
    return *this;
}

FloatMP::operator Dyadic() const {
    Dyadic res;
    mpfr_get_f(res._mpf,this->_mpfr, MPFR_RNDN);
    return res;
}

FloatMP::operator Rational() const {
    return Rational(Dyadic(*this));
    mpz_t num; mpz_init(num);
    mpfr_exp_t exp = mpfr_get_z_2exp (num, this->_mpfr);
    mpq_t res; mpq_init(res); mpq_set_z(res,num);
    if(exp>=0) { mpq_mul_2exp(res,res,exp); }
    else { mpq_div_2exp(res,res,-exp); }
    return Rational(res);
}

const FloatMP::RoundingModeType FloatMP::upward = MPFR_RNDU;

const FloatMP::RoundingModeType FloatMP::downward = MPFR_RNDD;

const FloatMP::RoundingModeType FloatMP::to_nearest = MPFR_RNDN;

const FloatMP::RoundingModeType FloatMP::toward_zero = MPFR_RNDZ;

FloatMP::RoundingModeType FloatMP::get_rounding_mode() {
    return mpfr_get_default_rounding_mode();
}

Void FloatMP::set_rounding_mode(FloatMP::RoundingModeType rnd) {
    mpfr_set_default_rounding_mode(rnd);
}

Void FloatMP::set_rounding_upward() {
    set_rounding_mode(FloatMP::upward);
}

Void FloatMP::set_rounding_downward() {
    set_rounding_mode(FloatMP::downward);
}

Void FloatMP::set_rounding_to_nearest() {
    set_rounding_mode(FloatMP::to_nearest);
}

Void FloatMP::set_rounding_toward_zero() {
    set_rounding_mode(FloatMP::toward_zero);
}

/*
Void FloatMP::set_rounding_mode(RoundingModeMP rnd) {
    switch(rnd) {
        case RoundingModeMP::NEAREST: set_rounding_mode(MPFR_RNDN; break);
        case RoundingModeMP::DOWNWARD: set_rounding_mode(MPFR_RNDD; break);
        case RoundingModeMP::UPWARD: set_rounding_mode(MPFR_RNDU; break);
        default: assert(false);
    }
}
*/

Void FloatMP::set_default_precision(PrecisionType pr) {
    mpfr_set_default_prec(pr);
}

FloatMP::PrecisionType FloatMP::get_default_precision() {
    return PrecisionType(mpfr_get_default_prec());
}

Void FloatMP::set_precision(PrecisionMP pr) {
    mpfr_t tmp;
    mpfr_init2(tmp,pr);
    mpfr_set(tmp,_mpfr,get_rounding_mode());
    mpfr_swap(_mpfr,tmp);
}

FloatMP::ExponentType FloatMP::exponent() const {
    return mpfr_get_exp(this->_mpfr);
}

FloatMP::PrecisionType FloatMP::precision() const {
    return PrecisionType(mpfr_get_prec(this->_mpfr));
}

double FloatMP::get_d() const {
    return mpfr_get_d(this->_mpfr,get_rounding_mode());
}

FloatMP FloatMP::nan(PrecisionMP pr) {
    FloatMP x(pr);
    mpfr_set_nan(x._mpfr);
    return x;
}

FloatMP FloatMP::inf(PrecisionMP pr) {
    FloatMP x(pr);
    mpfr_set_inf(x._mpfr,+1);
    return x;
}

FloatMP FloatMP::eps(PrecisionMP pr) {
    FloatMP x(pr);
    mpfr_set_ui_2exp(x._mpfr,1u,1-pr.bits(),FloatMP::to_nearest);
    return x;
}

FloatMP FloatMP::min(PrecisionMP pr) {
    mpfr_exp_t emin=mpfr_get_emin();
    FloatMP x(pr);
    mpfr_set_ui_2exp(x._mpfr,1u,emin-1,FloatMP::to_nearest);
    return x;
}

FloatMP FloatMP::max(PrecisionMP pr) {
    mpfr_exp_t emax=mpfr_get_emax();
    FloatMP x(2,pr);
    x=sub(x,FloatMP::min(pr),MPFR_RNDD);
    assert(x<2);
    FloatMP e(pr);
    mpfr_set_ui_2exp(e._mpfr,1u,(emax-1),FloatMP::to_nearest);
    return x*e;
}

Bool is_nan(FloatMP const& x) {
    return mpfr_nan_p(x._mpfr);
}

Bool is_inf(FloatMP const& x) {
    return mpfr_inf_p(x._mpfr);
}


FloatMP nul(FloatMP const& x) { return nul(x,FloatMP::get_rounding_mode()); }
FloatMP hlf(FloatMP const& x) { return hlf(x,FloatMP::get_rounding_mode()); }
FloatMP pos(FloatMP const& x) { return pos(x,FloatMP::get_rounding_mode()); }
FloatMP neg(FloatMP const& x) { return neg(x,FloatMP::get_rounding_mode()); }
FloatMP sqr(FloatMP const& x) { return sqr(x,FloatMP::get_rounding_mode()); }
FloatMP rec(FloatMP const& x) { return rec(x,FloatMP::get_rounding_mode()); }
FloatMP add(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP sub(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP mul(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP div(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3) { return fma(x1,x2,x3,FloatMP::get_rounding_mode()); }
FloatMP pow(FloatMP const& x, Int n) { return pow(x,n,FloatMP::get_rounding_mode()); }
FloatMP sqrt(FloatMP const& x) { return sqrt(x,FloatMP::get_rounding_mode()); }
FloatMP exp(FloatMP const& x) { return exp(x,FloatMP::get_rounding_mode()); }
FloatMP log(FloatMP const& x) { return log(x,FloatMP::get_rounding_mode()); }
FloatMP sin(FloatMP const& x) { return sin(x,FloatMP::get_rounding_mode()); }
FloatMP cos(FloatMP const& x) { return cos(x,FloatMP::get_rounding_mode()); }
FloatMP tan(FloatMP const& x) { return tan(x,FloatMP::get_rounding_mode()); }
FloatMP asin(FloatMP const& x) { return asin(x,FloatMP::get_rounding_mode()); }
FloatMP acos(FloatMP const& x) { return acos(x,FloatMP::get_rounding_mode()); }
FloatMP atan(FloatMP const& x) { return atan(x,FloatMP::get_rounding_mode()); }
FloatMP FloatMP::pi(PrecisionMP pr) { return pi(pr,FloatMP::get_rounding_mode()); }

FloatMP max(FloatMP const& x1, FloatMP const& x2) { return max(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP min(FloatMP const& x1, FloatMP const& x2) { return min(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP abs(FloatMP const& x) { return abs(x,FloatMP::get_rounding_mode()); }
FloatMP mag(FloatMP const& x) { return mag(x,FloatMP::get_rounding_mode()); }

    // Mixed operations
FloatMP add(FloatMP const& x1, Dbl x2) { return add(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP sub(FloatMP const& x1, Dbl x2) { return sub(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP mul(FloatMP const& x1, Dbl x2) { return mul(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP div(FloatMP const& x1, Dbl x2) { return div(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP add(Dbl x1, FloatMP const& x2) { return add(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP sub(Dbl x1, FloatMP const& x2) { return sub(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP mul(Dbl x1, FloatMP const& x2) { return mul(x1,x2,FloatMP::get_rounding_mode()); }
FloatMP div(Dbl x1, FloatMP const& x2) { return div(x1,x2,FloatMP::get_rounding_mode()); }

FloatMP operator+(FloatMP const& x) {
    return x;
}

FloatMP operator-(FloatMP const& x1) {
    return neg(x1,FloatMP::get_rounding_mode());
}

FloatMP operator+(FloatMP const& x1, FloatMP const& x2) {
    return add(x1,x2,FloatMP::get_rounding_mode());
}

FloatMP operator-(FloatMP const& x1, FloatMP const& x2) {
    return sub(x1,x2,FloatMP::get_rounding_mode());
}

FloatMP operator*(FloatMP const& x1, FloatMP const& x2) {
    return mul(x1,x2,FloatMP::get_rounding_mode());
}

FloatMP operator/(FloatMP const& x1, FloatMP const& x2) {
    return div(x1,x2,FloatMP::get_rounding_mode());
}


FloatMP& operator+=(FloatMP& x1, FloatMP const& x2) {
    mpfr_add(x1._mpfr,x1._mpfr,x2._mpfr,FloatMP::get_rounding_mode()); return x1;
}

FloatMP& operator-=(FloatMP& x1, FloatMP const& x2) {
    mpfr_sub(x1._mpfr,x1._mpfr,x2._mpfr,FloatMP::get_rounding_mode()); return x1;
}

FloatMP& operator*=(FloatMP& x1, FloatMP const& x2) {
    mpfr_mul(x1._mpfr,x1._mpfr,x2._mpfr,FloatMP::get_rounding_mode()); return x1;
}

FloatMP& operator/=(FloatMP& x1, FloatMP const& x2) {
    mpfr_div(x1._mpfr,x1._mpfr,x2._mpfr,FloatMP::get_rounding_mode()); return x1;
}



FloatMP operator+(FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_add_d(r._mpfr,x1._mpfr,x2,FloatMP::get_rounding_mode()); return r;
}
FloatMP operator-(FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_sub_d(r._mpfr,x1._mpfr,x2,FloatMP::get_rounding_mode()); return r;
}
FloatMP operator*(FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_mul_d(r._mpfr,x1._mpfr,x2,FloatMP::get_rounding_mode()); return r;
}
FloatMP operator/(FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_div_d(r._mpfr,x1._mpfr,x2,FloatMP::get_rounding_mode()); return r;
}

FloatMP operator+(Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_add_d(r._mpfr,x2._mpfr,x1,FloatMP::get_rounding_mode()); return r;
}
FloatMP operator-(Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_d_sub(r._mpfr,x1,x2._mpfr,FloatMP::get_rounding_mode()); return r;
}
FloatMP operator*(Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_mul_d(r._mpfr,x2._mpfr,x1,FloatMP::get_rounding_mode()); return r;
}
FloatMP operator/(Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_d_div(r._mpfr,x1,x2._mpfr,FloatMP::get_rounding_mode()); return r;
}

inline int log10floor(double const& x) { return std::max(std::floor(std::log10(x)),-65280.); }
inline int log10floor(FloatMP const& x) { return log10floor(x.get_d()); }
inline int abslog10floor(double const& x) { return log10floor(std::abs(x)); }


String print(const mpfr_t x, int zdgts, int fdgts, mpfr_rnd_t rnd) {
    char fmt[16];
    std::strcpy(fmt,"%.");
    std::sprintf(fmt+2, "%d", fdgts);
    std::strcat(fmt,"R*f");
    static const uint buf_size=1024;
    char cstr[buf_size];
    // uint buf_size = zdgts+fdgts+8;
    mpfr_snprintf(cstr,buf_size,fmt,rnd,x);
    cstr[1023]='\0';
    return String(cstr);
};

String print(FloatMP const& x, DecimalPrecision figs, RoundingModeMP rnd) {
    static const double log2ten = 3.3219280948873621817;
    int pdgts = std::ceil(x.precision()/log2ten);
    int zdgts = std::max(log10floor(x),0)+1;
    int fdgts = pdgts-zdgts;
    return print(x._mpfr,zdgts,fdgts,rnd);
}

String print(FloatMP const& x, DecimalPlaces plcs, RoundingModeMP rnd) {
    static const double log2ten = 3.3219280948873621817;
    int zdgts = std::max(log10floor(x),0)+1;
    int fdgts = plcs;
    return print(x._mpfr,zdgts,fdgts,rnd);
};

OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPlaces plcs, RoundingModeMP rnd) {
    return os << print(x,plcs,rnd);
}


OutputStream& operator<<(OutputStream& os, FloatMP const& x) {
    static const double log2ten = 3.3219280948873621817;
    return os << print(x,DecimalPlaces(std::max((int)x.precision()-(int)x.exponent(),0)/log2ten),FloatMP::get_rounding_mode());
}

InputStream& operator>>(InputStream& is, FloatMP& x) {
    char c;
    std::string str;
    c=is.get();
    while( c==' ' or c=='\t' or c=='\n' ) {
        c=is.get();
    }
    while( (c>='0' and c<='9') or c=='+' or c=='-' or c=='.' or c=='e' ) {
        str.push_back(c);
        c=is.get();
    }
    is.putback(c);
    mpfr_set_str(x._mpfr,str.c_str(),0,MPFR_RNDN);
    return is;
}

//Integer floor(FloatMP const& x) {
//    Integer z; mpfr_get_z(z._mpz,x._mpfr,MPFR_RNDD); return std::move(z);
//}
//Integer ceil(FloatMP const& x) {
//    Integer z; mpfr_get_z(z._mpz,x._mpfr,MPFR_RNDU); return std::move(z);
//}

FloatMP floor(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_floor(r._mpfr,x._mpfr); return r;
};
FloatMP ceil(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_ceil(r._mpfr,x._mpfr); return r;
};
FloatMP round(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_round(r._mpfr,x._mpfr); return r;
};

FloatMP abs(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_abs(r._mpfr,x._mpfr,MPFR_RNDN); return r;
}

FloatMP max(FloatMP const& x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_max(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}

FloatMP min(FloatMP const& x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_min(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}

FloatMP mag(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_abs(r._mpfr,x._mpfr,MPFR_RNDN); return r;
}

FloatMP nul(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_set_si(r._mpfr,0,rnd); return r;
}

FloatMP pos(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_set(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP neg(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_neg(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP sqr(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_sqr(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP hlf(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_div_si(r._mpfr,x._mpfr,2,rnd); return r;
}

FloatMP rec(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_ui_div(r._mpfr,1u,x._mpfr,rnd); return r;
}

FloatMP sqrt(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_sqrt(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP exp(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_exp(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP log(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_log(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP sin(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_sin(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP cos(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_cos(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP tan(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_tan(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP asin(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_asin(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP acos (FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_acos (r._mpfr,x._mpfr,rnd); return r;
}

FloatMP atan(FloatMP const& x, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_atan(r._mpfr,x._mpfr,rnd); return r;
}


FloatMP FloatMP::pi(PrecisionMP pr, FloatMP::RoundingModeType rnd) {
    FloatMP r(pr); mpfr_const_pi(r._mpfr,rnd); return r;
}

FloatMP add(FloatMP const& x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_add(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP sub(FloatMP const& x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_sub(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP mul(FloatMP const& x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_mul(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP div(FloatMP const& x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_div(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3, FloatMP::RoundingModeType rnd) {
    FloatMP r(cmb(x1.precision(),x2.precision(),x3.precision()),NoInit()); mpfr_fma(r._mpfr,x1._mpfr,x2._mpfr,x3._mpfr,rnd); return r;
}
FloatMP pow(FloatMP const& x, Int n, FloatMP::RoundingModeType rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_pow_si(r._mpfr,x._mpfr,n,rnd); return r;
}

Comparison cmp(FloatMP const& x1, FloatMP const& x2) {
    auto c=mpfr_cmp(x1._mpfr,x2._mpfr);
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

Comparison cmp(FloatMP const& x1, Dbl x2) {
    auto c=mpfr_cmp_d(x1._mpfr,x2);
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

Comparison cmp(FloatMP const& x1, Rational const& q2) {
    auto c=mpfr_cmp_q(x1._mpfr,q2.get_mpq());
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

FloatMP add(FloatMP const& x1, Dbl x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x1.precision(),NoInit()); mpfr_add_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP sub(FloatMP const& x1, Dbl x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x1.precision(),NoInit()); mpfr_sub_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP mul(FloatMP const& x1, Dbl x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x1.precision(),NoInit()); mpfr_mul_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP div(FloatMP const& x1, Dbl x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x1.precision(),NoInit()); mpfr_div_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP add(Dbl x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x2.precision(),NoInit()); mpfr_add_d(r._mpfr,x2._mpfr,x1,rnd); return r;
}
FloatMP sub(Dbl x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x2.precision(),NoInit()); mpfr_d_sub(r._mpfr,x1,x2._mpfr,rnd); return r;
}
FloatMP mul(Dbl x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x2.precision(),NoInit()); mpfr_mul_d(r._mpfr,x2._mpfr,x1,rnd); return r;
}
FloatMP div(Dbl x1, FloatMP const& x2, FloatMP::RoundingModeType rnd) {
    FloatMP r(x2.precision(),NoInit()); mpfr_d_div(r._mpfr,x1,x2._mpfr,rnd); return r;
}

Bool operator==(FloatMP const& x1, FloatMP const& x2) {
    return mpfr_equal_p(x1._mpfr,x2._mpfr);
}
Bool operator!=(FloatMP const& x1, FloatMP const& x2) {
    return not mpfr_equal_p(x1._mpfr,x2._mpfr);
}
Bool operator<=(FloatMP const& x1, FloatMP const& x2) {
    return mpfr_lessequal_p(x1._mpfr,x2._mpfr);
}
Bool operator>=(FloatMP const& x1, FloatMP const& x2) {
    return mpfr_greaterequal_p(x1._mpfr,x2._mpfr);
}
Bool operator< (FloatMP const& x1, FloatMP const& x2) {
    return mpfr_less_p(x1._mpfr,x2._mpfr);
}
Bool operator> (FloatMP const& x1, FloatMP const& x2) {
    return mpfr_greater_p(x1._mpfr,x2._mpfr);
}


Bool operator==(FloatMP const& x1, Dbl x2) {
    return cmp(x1,x2)==Comparison::EQUAL;
}
Bool operator!=(FloatMP const& x1, Dbl x2) {
    return cmp(x1,x2)!=Comparison::EQUAL;
}
Bool operator<=(FloatMP const& x1, Dbl x2) {
    return cmp(x1,x2)<=Comparison::EQUAL;
}
Bool operator>=(FloatMP const& x1, Dbl x2) {
    return cmp(x1,x2)>=Comparison::EQUAL;
}
Bool operator< (FloatMP const& x1, Dbl x2) {
    return cmp(x1,x2)< Comparison::EQUAL;
}
Bool operator> (FloatMP const& x1, Dbl x2) {
    return cmp(x1,x2)> Comparison::EQUAL;
}

FloatMP pos_exact(FloatMP const& x) { return pos(x,MPFR_RNDN); }
FloatMP neg_exact(FloatMP const& x) { return neg(x,MPFR_RNDN); }
FloatMP abs_exact(FloatMP const& x) { return abs(x,MPFR_RNDN); }
FloatMP half_exact(FloatMP const& x) { return hlf(x,MPFR_RNDN); }

FloatMP next_down(FloatMP const& x) { FloatMP r(x); mpfr_nextbelow(r._mpfr); return r; }
FloatMP next_up(FloatMP const& x) { FloatMP r(x); mpfr_nextabove(r._mpfr); return r; }

FloatMP neg_up(FloatMP const& x) { return neg(x,MPFR_RNDU); }
FloatMP rec_up(FloatMP const& x) { return rec(x,MPFR_RNDU); }
FloatMP add_up(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDU); }
FloatMP sub_up(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDU); }
FloatMP mul_up(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDU); }
FloatMP div_up(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDU); }
FloatMP pow_up(FloatMP const& x, Int n) { return pow(x,n,MPFR_RNDU); }

FloatMP neg_down(FloatMP const& x) { return neg(x,MPFR_RNDD); }
FloatMP rec_down(FloatMP const& x) { return rec(x,MPFR_RNDD); }
FloatMP add_down(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDD); }
FloatMP sub_down(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDD); }
FloatMP mul_down(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDD); }
FloatMP div_down(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDD); }
FloatMP pow_down(FloatMP const& x, Int n) { return pow(x,n,MPFR_RNDD); }

FloatMP neg_near(FloatMP const& x) { return neg(x,MPFR_RNDN); }
// FloatMP rec_near(FloatMP const& x) { return rec(x,MPFR_RNDN); }
FloatMP add_near(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDN); }
FloatMP sub_near(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDN); }
FloatMP mul_near(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDN); }
FloatMP div_near(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDN); }
FloatMP pow_near(FloatMP const& x, Int n) { return pow(x,n,MPFR_RNDN); }

FloatMP neg_approx(FloatMP const& x) { return neg(x,MPFR_RNDN); }
FloatMP rec_approx(FloatMP const& x) { return rec(x,MPFR_RNDN); }
FloatMP add_approx(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDN); }
FloatMP sub_approx(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDN); }
FloatMP mul_approx(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDN); }
FloatMP div_approx(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDN); }
FloatMP pow_approx(FloatMP const& x, Int n) { return pow(x,n,MPFR_RNDN); }
FloatMP sqrt_approx(FloatMP const& x) { return sqrt(x,MPFR_RNDN); }
FloatMP exp_approx(FloatMP const& x) { return exp(x,MPFR_RNDN); }
FloatMP log_approx(FloatMP const& x) { return log(x,MPFR_RNDN); }
FloatMP sin_approx(FloatMP const& x) { return sin(x,MPFR_RNDN); }
FloatMP cos_approx(FloatMP const& x) { return cos(x,MPFR_RNDN); }
FloatMP tan_approx(FloatMP const& x) { return tan(x,MPFR_RNDN); }
FloatMP atan_approx(FloatMP const& x) { return atan(x,MPFR_RNDN); }

FloatMP med_approx(FloatMP const& x1, FloatMP const& x2) { return hlf(add(x1,x2,MPFR_RNDN),MPFR_RNDN); }
FloatMP rad_up(FloatMP const& x1, FloatMP const& x2) { return hlf(abs(sub(x2,x1,MPFR_RNDU),MPFR_RNDU),MPFR_RNDU); }

FloatMP sqr_rnd(FloatMP const& x) { return sqr(x,FloatMP::get_rounding_mode()); }
FloatMP add_rnd(FloatMP const& x, FloatMP const& y) { return add(x,y,FloatMP::get_rounding_mode()); }
FloatMP sub_rnd(FloatMP const& x, FloatMP const& y) { return sub(x,y,FloatMP::get_rounding_mode()); }
FloatMP mul_rnd(FloatMP const& x, FloatMP const& y) { return mul(x,y,FloatMP::get_rounding_mode()); }
FloatMP div_rnd(FloatMP const& x, FloatMP const& y) { return div(x,y,FloatMP::get_rounding_mode()); }
FloatMP pow_rnd(FloatMP const& x, Int n) { return pow(x,n,FloatMP::get_rounding_mode()); }
FloatMP sqrt_rnd(FloatMP const& x) { return sqrt(x,FloatMP::get_rounding_mode()); }
FloatMP exp_rnd(FloatMP const& x) { return exp(x,FloatMP::get_rounding_mode()); }
FloatMP log_rnd(FloatMP const& x) { return log(x,FloatMP::get_rounding_mode()); }
FloatMP sin_rnd(FloatMP const& x) { return sin(x,FloatMP::get_rounding_mode()); }
FloatMP cos_rnd(FloatMP const& x) { return cos(x,FloatMP::get_rounding_mode()); }
FloatMP tan_rnd(FloatMP const& x) { return tan(x,FloatMP::get_rounding_mode()); }
FloatMP atan_rnd(FloatMP const& x) { return atan(x,FloatMP::get_rounding_mode()); }

FloatMP add_opp(FloatMP const& x, FloatMP const& y);
FloatMP sub_opp(FloatMP const& x, FloatMP const& y);
FloatMP mul_opp(FloatMP const& x, FloatMP const& y);
FloatMP div_opp(FloatMP const& x, FloatMP const& y);

FloatMP sqrt_near(FloatMP const& x) { return sqrt(x,MPFR_RNDN); }
FloatMP exp_near(FloatMP const& x) { return exp(x,MPFR_RNDN); }
FloatMP log_near(FloatMP const& x) { return log(x,MPFR_RNDN); }
FloatMP sin_near(FloatMP const& x) { return sin(x,MPFR_RNDN); }
FloatMP cos_near(FloatMP const& x) { return cos(x,MPFR_RNDN); }
FloatMP tan_near(FloatMP const& x) { return tan(x,MPFR_RNDN); }
FloatMP atan_near(FloatMP const& x) { return atan(x,MPFR_RNDN); }

FloatMP sqrt_up(FloatMP const& x) { return sqrt(x,MPFR_RNDU); }
FloatMP exp_up(FloatMP const& x) { return exp(x,MPFR_RNDU); }
FloatMP log_up(FloatMP const& x) { return log(x,MPFR_RNDU); }
FloatMP sin_up(FloatMP const& x) { return sin(x,MPFR_RNDU); }
FloatMP cos_up(FloatMP const& x) { return cos(x,MPFR_RNDU); }
FloatMP tan_up(FloatMP const& x) { return tan(x,MPFR_RNDU); }
FloatMP atan_up(FloatMP const& x) { return atan(x,MPFR_RNDU); }

FloatMP sqrt_down(FloatMP const& x) { return sqrt(x,MPFR_RNDD); }
FloatMP exp_down(FloatMP const& x) { return exp(x,MPFR_RNDD); }
FloatMP log_down(FloatMP const& x) { return log(x,MPFR_RNDD); }
FloatMP sin_down(FloatMP const& x) { return sin(x,MPFR_RNDD); }
FloatMP cos_down(FloatMP const& x) { return cos(x,MPFR_RNDD); }
FloatMP tan_down(FloatMP const& x) { return tan(x,MPFR_RNDD); }
FloatMP atan_down(FloatMP const& x) { return atan(x,MPFR_RNDD); }

FloatMP pi_near(PrecisionMP pr) { return FloatMP::pi(pr,MPFR_RNDN); }
FloatMP pi_down(PrecisionMP pr) { return FloatMP::pi(pr,MPFR_RNDD); }
FloatMP pi_up(PrecisionMP pr) { return FloatMP::pi(pr,MPFR_RNDU); }

} // namespace Ariadne
