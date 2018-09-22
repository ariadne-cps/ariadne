/***************************************************************************
 *            floatmp.cpp
 *
 *  Copyright 2013--17  Pieter Collins
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

/*! \file floatmp.cpp
 *  \brief
 */



#include "../utility/module.hpp"
#include "logical.hpp"
#include "floatmp.hpp"
#include "floatdp.hpp"
#include "dyadic.hpp"
#include "rational.hpp"

namespace Ariadne {

// Rule for combining mixed precision
inline MultiplePrecision cmb(MultiplePrecision pr1, MultiplePrecision pr2) { return min(pr1,pr2); }
inline MultiplePrecision cmb(MultiplePrecision pr1, MultiplePrecision pr2, MultiplePrecision pr3) { return cmb(cmb(pr1,pr2),pr3); }

FloatMP::~FloatMP() {
    mpfr_clear(_mpfr);
}

FloatMP::FloatMP() {
    mpfr_init_set_si(_mpfr,0l,get_rounding_mode());
}

FloatMP::FloatMP(NoInit) {
    mpfr_init(_mpfr);
}

FloatMP::FloatMP(const mpfr_t x, RawPtr) : FloatMP(NoInit()) {
    mpfr_set_prec(this->_mpfr,mpfr_get_prec(x));
    mpfr_set(this->_mpfr,x,MPFR_RNDN);
}


FloatMP::FloatMP(double d) : FloatMP(d,get_default_precision()) {
}

FloatMP::FloatMP(double d, MultiplePrecision pr) : FloatMP(d,MPFR_RNDN,pr) {
    ARIADNE_ASSERT(d==this->get_d());
}

FloatMP::FloatMP(FloatDP const& x, MultiplePrecision pr) : FloatMP(x.get_d(),pr) {
}

FloatMP::FloatMP(MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_si(_mpfr,0l,get_rounding_mode());
}

FloatMP::FloatMP(MultiplePrecision pr, NoInit) {
    mpfr_init2(_mpfr,pr);
}

FloatMP::FloatMP(Int32 n, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_si(_mpfr,n.get_si(),get_rounding_mode());
}

FloatMP::FloatMP(Dyadic const& w, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_f(_mpfr,w.get_mpf(),get_rounding_mode());
    ARIADNE_ASSERT(Dyadic(*this)==w);
}

FloatMP::FloatMP(double d, RoundingModeType rnd, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_d(_mpfr,d,rnd);
}

FloatMP::FloatMP(FloatDP const& x, RoundingModeType rnd, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_d(_mpfr,x.get_d(),rnd);
}

FloatMP::FloatMP(Integer const& z, RoundingModeType rnd, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_z(_mpfr,z.get_mpz(),rnd);
}

FloatMP::FloatMP(Dyadic const& w, RoundingModeType rnd, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_f(_mpfr,w.get_mpf(),rnd);
}

FloatMP::FloatMP(Rational const& q, RoundingModeType rnd, MultiplePrecision pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_q(_mpfr,q.get_mpq(),rnd);
}

FloatMP::FloatMP(FloatMP const& x, RoundingModeType rnd, MultiplePrecision pr) {
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
    if(exp>=0) { mpq_mul_2exp(res,res,static_cast<mp_bitcnt_t>(exp)); }
    else { mpq_div_2exp(res,res,static_cast<mp_bitcnt_t>(-exp)); }
    return Rational(res);
}

const FloatMP::RoundingModeType FloatMP::ROUND_TO_NEAREST = MPFR_RNDN;

const FloatMP::RoundingModeType FloatMP::ROUND_UPWARD = MPFR_RNDU;

const FloatMP::RoundingModeType FloatMP::ROUND_DOWNWARD = MPFR_RNDD;

const FloatMP::RoundingModeType FloatMP::ROUND_TOWARD_ZERO = MPFR_RNDZ;

FloatMP::RoundingModeType FloatMP::get_rounding_mode() {
    return mpfr_get_default_rounding_mode();
}

Void FloatMP::set_rounding_mode(FloatMP::RoundingModeType rnd) {
    mpfr_set_default_rounding_mode(rnd);
}

Void FloatMP::set_rounding_to_nearest() {
    set_rounding_mode(to_nearest);
}

Void FloatMP::set_rounding_downward() {
    set_rounding_mode(downward);
}

Void FloatMP::set_rounding_upward() {
    set_rounding_mode(upward);
}

Void FloatMP::set_rounding_toward_zero() {
    set_rounding_mode(toward_zero);
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

Void FloatMP::set_precision(MultiplePrecision pr) {
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

FloatMP FloatMP::nan(MultiplePrecision pr) {
    FloatMP x(pr);
    mpfr_set_nan(x._mpfr);
    return x;
}

FloatMP FloatMP::inf(MultiplePrecision pr) {
    FloatMP x(pr);
    mpfr_set_inf(x._mpfr,+1);
    return x;
}

FloatMP FloatMP::inf(Sign sgn, MultiplePrecision pr) {
    FloatMP x(pr);
    switch (sgn) {
    case Sign::POSITIVE: mpfr_set_inf(x._mpfr,+1); break;
    case Sign::NEGATIVE: mpfr_set_inf(x._mpfr,-1); break;
    default: mpfr_set_nan(x._mpfr);
    }
    return x;
}

FloatMP FloatMP::eps(MultiplePrecision pr) {
    FloatMP x(pr);
    mpfr_set_ui_2exp(x._mpfr,1u,1-mpfr_exp_t(pr.bits()),to_nearest);
    return x;
}

FloatMP FloatMP::min(MultiplePrecision pr) {
    mpfr_exp_t emin=mpfr_get_emin();
    FloatMP x(pr);
    mpfr_set_ui_2exp(x._mpfr,1u,emin-1,to_nearest);
    return x;
}

FloatMP FloatMP::max(MultiplePrecision pr) {
    mpfr_exp_t emax=mpfr_get_emax();
    FloatMP x(2,pr);
    x=sub(down,x,FloatMP::min(pr));
    assert(x<2);
    FloatMP e(pr);
    mpfr_set_ui_2exp(e._mpfr,1u,(emax-1),to_nearest);
    return mul(to_nearest,e,x);
}

Bool is_nan(FloatMP const& x) {
    return mpfr_nan_p(x._mpfr);
}

Bool is_inf(FloatMP const& x) {
    return mpfr_inf_p(x._mpfr);
}

Bool is_finite(FloatMP const& x) {
    return mpfr_number_p(x._mpfr);
}

Bool is_zero(FloatMP const& x) {
    return mpfr_zero_p(x._mpfr);
}

FloatMP next(RoundUpward rnd, FloatMP const& x) { return add(rnd,x,FloatMP::min(x.precision())); }
FloatMP next(RoundDownward rnd, FloatMP const& x) { return sub(rnd,x,FloatMP::min(x.precision())); }


FloatMP nul(FloatMP const& x) { return nul(FloatMP::get_rounding_mode(),x); }
FloatMP hlf(FloatMP const& x) { return hlf(FloatMP::get_rounding_mode(),x); }
FloatMP pos(FloatMP const& x) { return pos(FloatMP::get_rounding_mode(),x); }
FloatMP neg(FloatMP const& x) { return neg(FloatMP::get_rounding_mode(),x); }
FloatMP sqr(FloatMP const& x) { return sqr(FloatMP::get_rounding_mode(),x); }
FloatMP rec(FloatMP const& x) { return rec(FloatMP::get_rounding_mode(),x); }
FloatMP add(FloatMP const& x1, FloatMP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(FloatMP const& x1, FloatMP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(FloatMP const& x1, FloatMP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(FloatMP const& x1, FloatMP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP fma(FloatMP const& x1, FloatMP const& x2, FloatMP const& x3) { return fma(FloatMP::get_rounding_mode(),x1,x2,x3); }
FloatMP pow(FloatMP const& x, Int n) { return pow(FloatMP::get_rounding_mode(),x,n); }
FloatMP sqrt(FloatMP const& x) { return sqrt(FloatMP::get_rounding_mode(),x); }
FloatMP exp(FloatMP const& x) { return exp(FloatMP::get_rounding_mode(),x); }
FloatMP log(FloatMP const& x) { return log(FloatMP::get_rounding_mode(),x); }
FloatMP sin(FloatMP const& x) { return sin(FloatMP::get_rounding_mode(),x); }
FloatMP cos(FloatMP const& x) { return cos(FloatMP::get_rounding_mode(),x); }
FloatMP tan(FloatMP const& x) { return tan(FloatMP::get_rounding_mode(),x); }
FloatMP asin(FloatMP const& x) { return asin(FloatMP::get_rounding_mode(),x); }
FloatMP acos(FloatMP const& x) { return acos(FloatMP::get_rounding_mode(),x); }
FloatMP atan(FloatMP const& x) { return atan(FloatMP::get_rounding_mode(),x); }
// FIXME FloatMP FloatMP::pi(MultiplePrecision pr) { return pi(FloatMP::get_rounding_mode(),pr); }

FloatMP max(FloatMP const& x1, FloatMP const& x2) { return max(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP min(FloatMP const& x1, FloatMP const& x2) { return min(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP abs(FloatMP const& x) { return abs(FloatMP::get_rounding_mode(),x); }
FloatMP mag(FloatMP const& x) { return mag(FloatMP::get_rounding_mode(),x); }

// Mixed operations
FloatMP add(FloatMP const& x1, Dbl x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(FloatMP const& x1, Dbl x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(FloatMP const& x1, Dbl x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(FloatMP const& x1, Dbl x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP add(Dbl x1, FloatMP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub(Dbl x1, FloatMP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul(Dbl x1, FloatMP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div(Dbl x1, FloatMP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }

FloatMP operator+(FloatMP const& x) {
    return x;
}

FloatMP operator-(FloatMP const& x) {
    return neg(FloatMP::get_rounding_mode(),x);
}

FloatMP operator+(FloatMP const& x1, FloatMP const& x2) {
    return add(FloatMP::get_rounding_mode(),x1,x2);
}

FloatMP operator-(FloatMP const& x1, FloatMP const& x2) {
    return sub(FloatMP::get_rounding_mode(),x1,x2);
}

FloatMP operator*(FloatMP const& x1, FloatMP const& x2) {
    return mul(FloatMP::get_rounding_mode(),x1,x2);
}

FloatMP operator/(FloatMP const& x1, FloatMP const& x2) {
    return div(FloatMP::get_rounding_mode(),x1,x2);
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
}

String print(FloatMP const& x, DecimalPrecision figs, RoundingModeMP rnd) {
    static const double log2ten = 3.3219280948873621817;
    int pdgts = std::ceil(x.precision()/log2ten);
    int zdgts = std::max(log10floor(x),0)+1;
    int fdgts = pdgts-zdgts;
    return print(x._mpfr,zdgts,fdgts,rnd);
}

String print(FloatMP const& x, DecimalPlaces plcs, RoundingModeMP rnd) {
    int zdgts = std::max(log10floor(x),0)+1;
    int fdgts = static_cast<int>(plcs);
    return print(x._mpfr,zdgts,fdgts,rnd);
}

OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPlaces plcs, RoundingModeMP rnd) {
    return os << print(x,plcs,rnd);
}

OutputStream& write(OutputStream& os, FloatMP const& x, DecimalPrecision plcs, RoundingModeMP rnd) {
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
}
FloatMP ceil(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_ceil(r._mpfr,x._mpfr); return r;
}
FloatMP round(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_round(r._mpfr,x._mpfr); return r;
}

FloatMP abs(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_abs(r._mpfr,x._mpfr,MPFR_RNDN); return r;
}

FloatMP max(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_max(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}

FloatMP min(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_min(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}

FloatMP mag(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_abs(r._mpfr,x._mpfr,MPFR_RNDN); return r;
}

FloatMP nul(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_set_si(r._mpfr,0,rnd); return r;
}

FloatMP pos(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_set(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP neg(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_neg(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP sqr(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_sqr(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP hlf(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_div_2ui(r._mpfr,x._mpfr,1u,rnd); return r;
}

FloatMP shft(FloatMP::RoundingModeType rnd, FloatMP const& x, Int n) {
    FloatMP r(x.precision(),NoInit());
    mpfr_mul_2si(r._mpfr,x._mpfr,n,rnd);
    return r;
}

FloatMP rec(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_ui_div(r._mpfr,1u,x._mpfr,rnd); return r;
}

FloatMP sqrt(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_sqrt(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP exp(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_exp(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP log(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_log(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP sin(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_sin(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP cos(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_cos(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP tan(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_tan(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP asin(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_asin(r._mpfr,x._mpfr,rnd); return r;
}

FloatMP acos (FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_acos (r._mpfr,x._mpfr,rnd); return r;
}

FloatMP atan(FloatMP::RoundingModeType rnd, FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_atan(r._mpfr,x._mpfr,rnd); return r;
}


FloatMP FloatMP::pi(FloatMP::RoundingModeType rnd, MultiplePrecision pr) {
    FloatMP r(pr); mpfr_const_pi(r._mpfr,rnd); return r;
}

FloatMP add(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_add(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP sub(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_sub(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP mul(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_mul(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP div(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2) {
    FloatMP r(cmb(x1.precision(),x2.precision()),NoInit()); mpfr_div(r._mpfr,x1._mpfr,x2._mpfr,rnd); return r;
}
FloatMP fma(FloatMP::RoundingModeType rnd, FloatMP const& x1, FloatMP const& x2, FloatMP const& x3) {
    FloatMP r(cmb(x1.precision(),x2.precision(),x3.precision()),NoInit()); mpfr_fma(r._mpfr,x1._mpfr,x2._mpfr,x3._mpfr,rnd); return r;
}
FloatMP pow(FloatMP::RoundingModeType rnd, FloatMP const& x, Int n) {
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

Comparison cmp(Dbl x1, FloatMP const& x2) {
    auto c=mpfr_cmp_d(x2._mpfr,x1);
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::LESS:Comparison::GREATER);
}

Comparison cmp(FloatMP const& x1, FloatDP const& x2) {
    return cmp(x1,x2.get_d());
}

Comparison cmp(FloatDP const& x1, FloatMP const& x2) {
    return cmp(x1.get_d(),x2);
}

Comparison cmp(FloatMP const& x1, Rational const& q2) {
    auto c=mpfr_cmp_q(x1._mpfr,q2.get_mpq());
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

Comparison cmp(Rational const& q1, FloatMP const& x2) {
    auto c=mpfr_cmp_q(x2._mpfr,q1.get_mpq());
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::LESS:Comparison::GREATER);
}

FloatMP add(FloatMP::RoundingModeType rnd, FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_add_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP sub(FloatMP::RoundingModeType rnd, FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_sub_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP mul(FloatMP::RoundingModeType rnd, FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_mul_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP div(FloatMP::RoundingModeType rnd, FloatMP const& x1, Dbl x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_div_d(r._mpfr,x1._mpfr,x2,rnd); return r;
}
FloatMP add(FloatMP::RoundingModeType rnd, Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_add_d(r._mpfr,x2._mpfr,x1,rnd); return r;
}
FloatMP sub(FloatMP::RoundingModeType rnd, Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_d_sub(r._mpfr,x1,x2._mpfr,rnd); return r;
}
FloatMP mul(FloatMP::RoundingModeType rnd, Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_mul_d(r._mpfr,x2._mpfr,x1,rnd); return r;
}
FloatMP div(FloatMP::RoundingModeType rnd, Dbl x1, FloatMP const& x2) {
    FloatMP r(x2.precision(),NoInit()); mpfr_d_div(r._mpfr,x1,x2._mpfr,rnd); return r;
}

FloatMP add(RoundUpward rnd, FloatMP const& x1, FloatDP const& x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_add_d(r._mpfr,x1._mpfr,x2.get_d(),rnd); return r;
}
FloatMP sub(RoundDownward rnd, FloatMP const& x1, FloatDP const& x2) {
    FloatMP r(x1.precision(),NoInit()); mpfr_sub_d(r._mpfr,x1._mpfr,x2.get_d(),rnd); return r;
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




FloatMP sqr_rnd(FloatMP const& x) { return sqr(FloatMP::get_rounding_mode(),x); }
FloatMP add_rnd(FloatMP const& x1, FloatMP const& x2) { return add(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP sub_rnd(FloatMP const& x1, FloatMP const& x2) { return sub(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP mul_rnd(FloatMP const& x1, FloatMP const& x2) { return mul(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP div_rnd(FloatMP const& x1, FloatMP const& x2) { return div(FloatMP::get_rounding_mode(),x1,x2); }
FloatMP pow_rnd(FloatMP const& x, Int n) { return pow(FloatMP::get_rounding_mode(),x,n); }
FloatMP sqrt_rnd(FloatMP const& x) { return sqrt(FloatMP::get_rounding_mode(),x); }
FloatMP exp_rnd(FloatMP const& x) { return exp(FloatMP::get_rounding_mode(),x); }
FloatMP log_rnd(FloatMP const& x) { return log(FloatMP::get_rounding_mode(),x); }
FloatMP sin_rnd(FloatMP const& x) { return sin(FloatMP::get_rounding_mode(),x); }
FloatMP cos_rnd(FloatMP const& x) { return cos(FloatMP::get_rounding_mode(),x); }
FloatMP tan_rnd(FloatMP const& x) { return tan(FloatMP::get_rounding_mode(),x); }
FloatMP atan_rnd(FloatMP const& x) { return atan(FloatMP::get_rounding_mode(),x); }

FloatMP add_opp(FloatMP const& x, FloatMP const& y);
FloatMP sub_opp(FloatMP const& x, FloatMP const& y);
FloatMP mul_opp(FloatMP const& x, FloatMP const& y);
FloatMP div_opp(FloatMP const& x, FloatMP const& y);

template<> String class_name<FloatMP>() { return "FloatMP"; }

} // namespace Ariadne
