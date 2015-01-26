/***************************************************************************
 *            floatmp.cc
 *
 *  Copyright 2013-14  Pieter Collins
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

namespace Ariadne {

mpfr_rnd_t FloatMP::current_rounding_mode=MPFR_RNDU;

FloatMP::~FloatMP() {
    mpfr_clear(_mpfr);
}

FloatMP::FloatMP() {
    mpfr_init_set_si(_mpfr,0l,current_rounding_mode);
}

FloatMP::FloatMP(NoInit) {
    mpfr_init(_mpfr);
}

FloatMP::FloatMP(PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_si(_mpfr,0l,current_rounding_mode);
}

FloatMP::FloatMP(PrecisionMP pr, NoInit) {
    mpfr_init2(_mpfr,pr);
}

FloatMP::FloatMP(Int32 n) {
    mpfr_init_set_si(_mpfr,n.get_si(),current_rounding_mode);
}

FloatMP::FloatMP(Int32 n, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_si(_mpfr,n.get_si(),current_rounding_mode);
}

FloatMP::FloatMP(Integer const& z, PrecisionMP pr, RoundingModeType rnd) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_z(_mpfr,z.get_mpz(),rnd);
}

FloatMP::FloatMP(Rational const& q, PrecisionMP pr, RoundingModeType rnd) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_q(_mpfr,q.get_mpq(),rnd);
}

FloatMP::FloatMP(double d) {
    mpfr_init(_mpfr);
    mpfr_set_d(_mpfr,d,current_rounding_mode);
//    mpfr_init_set_d(_mpfr,d,current_rounding_mode);
}

FloatMP::FloatMP(double d, PrecisionMP pr) {
    mpfr_init2(_mpfr,pr);
    mpfr_set_d(_mpfr,d,current_rounding_mode);
}

//FloatMP::FloatMP(const mpfr_t x) {
//    mpfr_init_set(_mpfr,x,current_rounding_mode);
//}

FloatMP::FloatMP(const FloatMP& x) {
    mpfr_init2(_mpfr,mpfr_get_prec(x._mpfr));
    mpfr_set(_mpfr,x._mpfr,current_rounding_mode);
}

FloatMP::FloatMP(FloatMP&& x) {
    mpfr_init(_mpfr);
    mpfr_swap(_mpfr,x._mpfr);
}

FloatMP& FloatMP::operator=(const FloatMP& x) {
    // TOOD: Decide whether equality changes precision
    mpfr_set_prec(_mpfr,x.precision());
    mpfr_set(_mpfr,x._mpfr,current_rounding_mode);
    return *this;
}

FloatMP& FloatMP::operator=(FloatMP&& x) {
    mpfr_swap(_mpfr,x._mpfr);
    return *this;
}

FloatMP& FloatMP::operator=(double d) {
    mpfr_set_d(_mpfr,d,current_rounding_mode);
    return *this;
}

FloatMP::operator Rational() const {
    mpz_t num; mpz_init(num);
    mpfr_exp_t exp = mpfr_get_z_2exp (num, this->_mpfr);
    mpq_t res; mpq_init(res); mpq_set_z(res,num);
    mpq_mul_2exp(res,res,exp);
    return Rational(res);
}

Void FloatMP::set_rounding_mode(RoundingModeMP rnd) {
    switch(rnd) {
        case RoundingModeMP::NEAREST: current_rounding_mode=MPFR_RNDN; break;
        case RoundingModeMP::DOWNWARD: current_rounding_mode=MPFR_RNDD; break;
        case RoundingModeMP::UPWARD: current_rounding_mode=MPFR_RNDU; break;
        default: assert(false);
    }
}

Void FloatMP::set_default_precision(PrecisionType pr) {
    mpfr_set_default_prec(pr);
}

FloatMP::PrecisionType FloatMP::get_default_precision() {
    return mpfr_get_default_prec();
}

Void FloatMP::set_precision(PrecisionMP pr) {
    mpfr_t tmp;
    mpfr_init2(tmp,pr);
    mpfr_set(tmp,_mpfr,current_rounding_mode);
    mpfr_swap(_mpfr,tmp);
}

FloatMP::PrecisionType FloatMP::precision() const {
    return mpfr_get_prec(this->_mpfr);
}

double FloatMP::get_d() const {
    return mpfr_get_d(this->_mpfr,current_rounding_mode);
}

FloatMP operator+(FloatMP const& x) {
    return x;
}

FloatMP operator-(FloatMP const& x) {
    return neg(x,FloatMP::current_rounding_mode);
}

FloatMP operator+(FloatMP const& x, FloatMP const& y) {
    return add(x,y,FloatMP::current_rounding_mode);
}

FloatMP operator-(FloatMP const& x, FloatMP const& y) {
    return sub(x,y,FloatMP::current_rounding_mode);
}

FloatMP operator*(FloatMP const& x, FloatMP const& y) {
    return mul(x,y,FloatMP::current_rounding_mode);
}

FloatMP operator/(FloatMP const& x, FloatMP const& y) {
    return div(x,y,FloatMP::current_rounding_mode);
}

String sprintf(FloatMP const& x, RoundingModeMP rnd);
String sprintf(FloatMP const&, PrecisionMP, RoundingModeMP);

//   mpfr_get_str (char *str, mpfr_exp_t *expptr, Int b, SizeType n, mpfr_t op, mpfr_rnd_t rnd)
// If str is not _a null pointer, it should point to _a block of storage large enough for the significand,
// i._e., at least max(n + 2, 7). The extra two bytes are for _a possible minus sign,
// and for the terminating null character, and the value 7 accounts for -@Inf@ plus the terminating null character.
OutputStream& operator<<(OutputStream& os, FloatMP const& x) {
    char str[1024];
//    mpfr_snprintf (str, 1024, "%R*e", MPFR_RNDN, x._mpfr);
    mpfr_snprintf (str, 1024, "%.100R*f", MPFR_RNDN, x._mpfr);
    return os << str;
}

//Integer floor(FloatMP const& x) {
//    Integer z; mpfr_get_z(z._mpz,x._mpfr,MPFR_RNDD); return std::move(z);
//}
//Integer ceil(FloatMP const& x) {
//    Integer z; mpfr_get_z(z._mpz,x._mpfr,MPFR_RNDU); return std::move(z);
//}

FloatMP floor(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_floor(r._mpfr,x._mpfr); return std::move(r);
};
FloatMP ceil(FloatMP const& x) {
    FloatMP r(x.precision()); mpfr_ceil(r._mpfr,x._mpfr); return std::move(r);
};

FloatMP nul(FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_set_si(r._mpfr,0,MPFR_RNDN); return std::move(r);
}

FloatMP pos(FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_set(r._mpfr,x._mpfr,MPFR_RNDN); return std::move(r);
}

FloatMP neg(FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_neg(r._mpfr,x._mpfr,MPFR_RNDN); return std::move(r);
}

FloatMP abs(FloatMP const& x) {
    FloatMP r(x.precision(),NoInit()); mpfr_abs(r._mpfr,x._mpfr,MPFR_RNDN); return std::move(r);
}

FloatMP max(FloatMP const& x, FloatMP const& y) {
    FloatMP r(std::max(x.precision(),y.precision()),NoInit()); mpfr_max(r._mpfr,x._mpfr,y._mpfr,MPFR_RNDN); return std::move(r);
}

FloatMP min(FloatMP const& x, FloatMP const& y) {
    FloatMP r(std::max(x.precision(),y.precision()),NoInit()); mpfr_min(r._mpfr,x._mpfr,y._mpfr,MPFR_RNDN); return std::move(r);
}

FloatMP half(FloatMP&& x) {
    mpfr_div_si(x._mpfr,x._mpfr,2,MPFR_RNDU); return std::move(x);
}

FloatMP pos(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_set(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP neg(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_neg(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP rec(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_ui_div(r._mpfr,1u,x._mpfr,rnd); return std::move(r);
}

FloatMP sqrt(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_sqrt(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP exp(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_exp(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP log(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_log(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP sin(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_sin(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP cos(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_cos(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP tan(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_tan(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP atan(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_atan(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP abs(FloatMP const& x, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_abs(r._mpfr,x._mpfr,rnd); return std::move(r);
}

FloatMP const_pi(PrecisionMP pr, mpfr_rnd_t rnd) {
    FloatMP r(pr); mpfr_const_pi(r._mpfr,rnd); return std::move(r);
}

FloatMP add(FloatMP const& x1, FloatMP const& x2, mpfr_rnd_t rnd) {
    FloatMP r(std::min(x1.precision(),x2.precision()),NoInit()); mpfr_add(r._mpfr,x1._mpfr,x2._mpfr,rnd); return std::move(r);
}

FloatMP sub(FloatMP const& x1, FloatMP const& x2, mpfr_rnd_t rnd) {
    FloatMP r(std::min(x1.precision(),x2.precision()),NoInit()); mpfr_sub(r._mpfr,x1._mpfr,x2._mpfr,rnd); return std::move(r);
}
FloatMP mul(FloatMP const& x1, FloatMP const& x2, mpfr_rnd_t rnd) {
    FloatMP r(std::min(x1.precision(),x2.precision()),NoInit()); mpfr_mul(r._mpfr,x1._mpfr,x2._mpfr,rnd); return std::move(r);
}
FloatMP div(FloatMP const& x1, FloatMP const& x2, mpfr_rnd_t rnd) {
    FloatMP r(std::min(x1.precision(),x2.precision()),NoInit()); mpfr_div(r._mpfr,x1._mpfr,x2._mpfr,rnd); return std::move(r);
}
FloatMP pow(FloatMP const& x, Int n, mpfr_rnd_t rnd) {
    FloatMP r(x.precision(),NoInit()); mpfr_pow_si(r._mpfr,x._mpfr,n,rnd); return std::move(r);
}

Comparison cmp(FloatMP const& x1, FloatMP const& x2) {
    return Comparison(mpfr_cmp(x1._mpfr,x2._mpfr));
}

Comparison cmp(FloatMP const& x1, long int n2) {
    return Comparison(mpfr_cmp_si(x1._mpfr,n2));
}

Comparison cmp(FloatMP const& x1, double x2) {
    return Comparison(mpfr_cmp_d(x1._mpfr,x2));
}

Bool operator==(FloatMP const& x1, FloatMP const& x2) {
    return cmp(x1,x2)==Comparison::EQUAL;
}
Bool operator!=(FloatMP const& x1, FloatMP const& x2) {
    return cmp(x1,x2)!=Comparison::EQUAL;
}
Bool operator<=(FloatMP const& x1, FloatMP const& x2) {
    return cmp(x1,x2)<=Comparison::EQUAL;
}
Bool operator>=(FloatMP const& x1, FloatMP const& x2) {
    return cmp(x1,x2)>=Comparison::EQUAL;
}
Bool operator< (FloatMP const& x1, FloatMP const& x2) {
    return cmp(x1,x2)< Comparison::EQUAL;
}
Bool operator> (FloatMP const& x1, FloatMP const& x2) {
    return cmp(x1,x2)> Comparison::EQUAL;
}


Bool operator==(FloatMP const& x1, double x2) {
    return cmp(x1,x2)==Comparison::EQUAL;
}
Bool operator!=(FloatMP const& x1, double x2) {
    return cmp(x1,x2)!=Comparison::EQUAL;
}
Bool operator<=(FloatMP const& x1, double x2) {
    return cmp(x1,x2)<=Comparison::EQUAL;
}
Bool operator>=(FloatMP const& x1, double x2) {
    return cmp(x1,x2)>=Comparison::EQUAL;
}
Bool operator< (FloatMP const& x1, double x2) {
    return cmp(x1,x2)< Comparison::EQUAL;
}
Bool operator> (FloatMP const& x1, double x2) {
    return cmp(x1,x2)> Comparison::EQUAL;
}


FloatMP next_down(FloatMP x) { ARIADNE_NOT_IMPLEMENTED; }
FloatMP next_up(FloatMP x) { ARIADNE_NOT_IMPLEMENTED; }

FloatMP neg_up(FloatMP const& x) { return neg(x,MPFR_RNDU); }
FloatMP rec_up(FloatMP const& x) { return rec(x,MPFR_RNDU); }
FloatMP add_up(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDU); }
FloatMP sub_up(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDU); }
FloatMP mul_up(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDU); }
FloatMP div_up(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDU); }

FloatMP neg_down(FloatMP const& x) { return neg(x,MPFR_RNDD); }
FloatMP rec_down(FloatMP const& x) { return rec(x,MPFR_RNDD); }
FloatMP add_down(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDD); }
FloatMP sub_down(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDD); }
FloatMP mul_down(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDD); }
FloatMP div_down(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDD); }

FloatMP neg_near(FloatMP const& x) { return neg(x,MPFR_RNDN); }
FloatMP rec_near(FloatMP const& x) { return rec(x,MPFR_RNDN); }
FloatMP add_near(FloatMP const& x1, FloatMP const& x2) { return add(x1,x2,MPFR_RNDN); }
FloatMP sub_near(FloatMP const& x1, FloatMP const& x2) { return sub(x1,x2,MPFR_RNDN); }
FloatMP mul_near(FloatMP const& x1, FloatMP const& x2) { return mul(x1,x2,MPFR_RNDN); }
FloatMP div_near(FloatMP const& x1, FloatMP const& x2) { return div(x1,x2,MPFR_RNDN); }
FloatMP pow_near(FloatMP const& x, Int n) { return pow(x,n,MPFR_RNDN); }


FloatMP add_rnd(FloatMP const& x, FloatMP const& y) { return add(x,y,FloatMP::current_rounding_mode); }
FloatMP sub_rnd(FloatMP const& x, FloatMP const& y) { return sub(x,y,FloatMP::current_rounding_mode); }
FloatMP mul_rnd(FloatMP const& x, FloatMP const& y) { return mul(x,y,FloatMP::current_rounding_mode); }
FloatMP div_rnd(FloatMP const& x, FloatMP const& y) { return div(x,y,FloatMP::current_rounding_mode); }
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

FloatMP pi_near(PrecisionMP pr) { return const_pi(pr,MPFR_RNDN); }
FloatMP pi_down(PrecisionMP pr) { return const_pi(pr,MPFR_RNDD); }
FloatMP pi_up(PrecisionMP pr) { return const_pi(pr,MPFR_RNDU); }

/************ FloatMPTemplate<Error> ***************************************************/

FloatMPTemplate<Error>::FloatMPTemplate(unsigned int e) : _e(e) { }

FloatMPTemplate<Error>::FloatMPTemplate(double e) : _e(e) { }

FloatMPTemplate<Error>::FloatMPTemplate(double e, PrecisionMP pr) : _e(pr) { _e=e; }

FloatMPTemplate<Error>::FloatMPTemplate(FloatMP e) : _e(e) { }

FloatMP const& FloatMPTemplate<Error>::get_flt() const { return this->_e; }

FloatMPTemplate<Error> operator+(FloatMPTemplate<Error> x, FloatMPTemplate<Error> y) {
    return FloatMPTemplate<Error>(add_up(x._e,y._e));
}

FloatMPTemplate<Error> operator*(FloatMPTemplate<Error> x, FloatMPTemplate<Error> y) {
    return FloatMPTemplate<Error>(mul_up(x._e,y._e));
}

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Error> const& x) {
    char str[1024];
    mpfr_snprintf (str, 1024, "%.2RUe", x._e._mpfr);
    return os << "\u00b1" << str;
}


/************ FloatMPTemplate<Exact> ***************************************************/

FloatMPTemplate<Exact>::FloatMPTemplate(double v) : _v(v) { }

FloatMPTemplate<Exact>::FloatMPTemplate(FloatMP v) : _v(v) { }

FloatMPTemplate<Exact>::operator FloatMPTemplate<Metrc>() const {
    return FloatMPTemplate<Metrc>(_v);
}
FloatMPTemplate<Exact>::operator FloatMPTemplate<Bound>() const {
    return FloatMPTemplate<Bound>(_v);
}
FloatMPTemplate<Exact>::operator FloatMPTemplate<Apprx>() const {
    return FloatMPTemplate<Apprx>(_v);
}

Void FloatMPTemplate<Exact>::set_precision(PrecisionMP pr) {
    _v.set_precision(pr);
}

PrecisionMP FloatMPTemplate<Exact>::precision() const {
    return _v.precision();
}

FloatMP const& FloatMPTemplate<Exact>::get_flt() const {
    return this->_v;
}

FloatMPTemplate<Metrc> FloatMPTemplate<Exact>::pm(FloatMPTemplate<Error> e) const {
    return FloatMPTemplate<Metrc>(this->_v,e.get_flt());
}

FloatMPTemplate<Exact> pos(FloatMPTemplate<Exact> x) {
    return FloatMPTemplate<Exact>(-x._v);
}

FloatMPTemplate<Exact> neg(FloatMPTemplate<Exact> x) {
    return FloatMPTemplate<Exact>(-x._v);
}

FloatMPTemplate<Bound> add(FloatMPTemplate<Exact> x, FloatMPTemplate<Exact> y) {
    return FloatMPTemplate<Bound>(add_down(x._v,y._v),add_up(x._v,y._v));
}

FloatMPTemplate<Bound> sub(FloatMPTemplate<Exact> x, FloatMPTemplate<Exact> y) {
    return FloatMPTemplate<Bound>(sub_down(x._v,y._v),sub_up(x._v,y._v));
}

FloatMPTemplate<Bound> mul(FloatMPTemplate<Exact> x, FloatMPTemplate<Exact> y) {
    return FloatMPTemplate<Bound>(mul_down(x._v,y._v),mul_up(x._v,y._v));
}

FloatMPTemplate<Bound> div(FloatMPTemplate<Exact> x, FloatMPTemplate<Exact> y) {
    return FloatMPTemplate<Bound>(div_down(x._v,y._v),div_up(x._v,y._v));
}

FloatMPTemplate<Bound> rec(FloatMPTemplate<Exact> x) {
    FloatMP one(1,x.precision());
    auto y=FloatMPTemplate<Bound>(div_down(one,x._v),div_up(one,x._v));
    std::cerr<<x.precision()<<" "<<one.precision()<<" "<<y.precision()<<"\n";
    return y;
}

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Exact> const& x) {
    char str[1024];
    mpfr_snprintf (str, 1024, "%.R*e", MPFR_RNDN, x._v._mpfr);
    return os << str;
}



/************ FloatMPTemplate<Metrc> ***************************************************/

FloatMPTemplate<Metrc>::FloatMPTemplate(Int32 n) : _v(n), _e(Int32(0)) {
}

FloatMPTemplate<Metrc>::FloatMPTemplate(Integer const& z, PrecisionMP pr) : FloatMPTemplate<Metrc>(Rational(z),pr) {
}

FloatMPTemplate<Metrc>::FloatMPTemplate(Rational const& q, PrecisionMP pr) : _v(q.get_d()), _e(Int32(0)) {
}

FloatMPTemplate<Metrc>::FloatMPTemplate(FloatMP v) : _v(v), _e(v) {
}

FloatMPTemplate<Metrc>::FloatMPTemplate(FloatMP v, FloatMP e) : _v(v), _e(e) {
}

FloatMPTemplate<Metrc>::operator FloatMPTemplate<Bound>() const {
    return FloatMPTemplate<Bound>(sub_down(_v,_e),add_up(_v,_e));
}
FloatMPTemplate<Metrc>::operator FloatMPTemplate<Upper>() const {
    return FloatMPTemplate<Upper>(add_up(_v,_e));
}
FloatMPTemplate<Metrc>::operator FloatMPTemplate<Lower>() const {
    return FloatMPTemplate<Lower>(add_down(_v,_e));
}
FloatMPTemplate<Metrc>::operator FloatMPTemplate<Apprx>() const {
    return FloatMPTemplate<Apprx>(_v);
}

Void FloatMPTemplate<Metrc>::set_precision(PrecisionMP pr) {
    _v.set_precision(pr);
    _e.set_precision(pr);
}

PrecisionMP FloatMPTemplate<Metrc>::precision() const {
    return _v.precision();
}

FloatMPTemplate<Exact> const& FloatMPTemplate<Metrc>::value() const {
    return reinterpret_cast<FloatMPTemplate<Exact>const&>(this->_v);
}

FloatMPTemplate<Error> const& FloatMPTemplate<Metrc>::error() const {
    return reinterpret_cast<FloatMPTemplate<Error>const&>(this->_e);
}

FloatMPTemplate<Metrc> nul(FloatMPTemplate<Metrc> x) {
    return FloatMPTemplate<Metrc>(nul(x._v),nul(x._e));
}

FloatMPTemplate<Metrc> pos(FloatMPTemplate<Metrc> x) {
    return FloatMPTemplate<Metrc>(+x._v,x._e);
}

FloatMPTemplate<Metrc> sqr(FloatMPTemplate<Metrc> x) {
    FloatMPTemplate<Metrc> r=x*x;
    if(r._e>r._v) {
        r._e=half(add_up(r._e,r._v));
        r._v=r._e;
    }
    return r;
}

FloatMPTemplate<Metrc> neg(FloatMPTemplate<Metrc> x) {
    return FloatMPTemplate<Metrc>(-x._v,x._e);
}

FloatMPTemplate<Metrc> rec(FloatMPTemplate<Metrc> x) {
    auto ru=rec_up(add_down(x._v,x._e));
    auto rl=rec_down(add_down(x._v,x._e));
    auto re=half(sub_up(ru,rl));
    auto rv=half(add_near(rl,ru));
    return FloatMPTemplate<Metrc>(rv,re);
}

FloatMPTemplate<Metrc> add(FloatMPTemplate<Metrc> x, FloatMPTemplate<Metrc> y) {
    auto rv=add_near(x._v,y._v);
    auto ru=add_up(x._v,y._v);
    auto rl=add_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return FloatMPTemplate<Metrc>(rv,re);
}

FloatMPTemplate<Metrc> sub(FloatMPTemplate<Metrc> x, FloatMPTemplate<Metrc> y) {
    auto rv=sub_near(x._v,y._v);
    auto ru=sub_up(x._v,y._v);
    auto rl=sub_down(x._v,y._v);
    auto re=add_up(half(sub_up(ru,rl)),add_up(x._e,y._e));
    return FloatMPTemplate<Metrc>(rv,re);
}

FloatMPTemplate<Metrc> mul(FloatMPTemplate<Metrc> x, FloatMPTemplate<Metrc> y) {
    auto rv=mul_near(x._v,y._v);
    auto ru=mul_up(x._v,y._v);
    auto rl=mul_down(x._v,y._v);
    auto re1=add_up(half(sub_up(ru,rl)),mul_up(x._e,y._e));
    auto re2=add_up(mul_up(abs(x._v),y._e),mul_up(x._e,abs(y._v)));
    auto re=add_up(re1,re2);
    return FloatMPTemplate<Metrc>(rv,re);
}

FloatMPTemplate<Metrc> div(FloatMPTemplate<Metrc> x, FloatMPTemplate<Metrc> y) {
    return x*rec(y);
}

FloatMPTemplate<Metrc> pow(FloatMPTemplate<Metrc> x, Int n) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> sqrt(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> exp(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> log(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> sin(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> cos(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> tan(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}

FloatMPTemplate<Metrc> atan(FloatMPTemplate<Metrc> x) {
    ARIADNE_NOT_IMPLEMENTED;
}


FloatMPTemplate<Metrc> abs(FloatMPTemplate<Metrc> x) {
    if(x._e<abs(x._v)) { return x; }
    else { FloatMP rv=half(abs(x._v)+x._e); return MetrcFloatMP(rv,rv); }
}

MetrcFloatMP max(MetrcFloatMP x1, MetrcFloatMP x2) {
    return ((x1+x2)+abs(x1-x2))/2;
}

MetrcFloatMP min(MetrcFloatMP x1, MetrcFloatMP x2) {
    return ((x1+x2)-abs(x1-x2))/2;
}

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Metrc> const& x) {
    return os << x.value() << x.error();
}



/************ FloatMPTemplate<Bound> ***************************************************/

FloatMP operator*(double d, FloatMP x) { mpfr_mul_d(x._mpfr,x._mpfr,d,MPFR_RNDN); return std::move(x); }
FloatMP operator+(FloatMP x, double d) { mpfr_add_d(x._mpfr,x._mpfr,d,MPFR_RNDN); return std::move(x); }


FloatMPTemplate<Bound>::FloatMPTemplate() : _l(0.0), _u(0.0) {
}

FloatMPTemplate<Bound>::FloatMPTemplate(Int32 n) : _l(n), _u(n) {
}

FloatMPTemplate<Bound>::FloatMPTemplate(const Integer& z, PrecisionMP pr) : FloatMPTemplate<Bound>(Rational(z),pr) {
}

FloatMPTemplate<Bound>::FloatMPTemplate(const Rational& q, PrecisionMP pr) : _l(q.get_d(),pr), _u(q.get_d(),pr) {
    // FIXME: Change this to increment up
    while(Rational(_l)>q) { _l=next_down(_l); }
    while(Rational(_u)<q) { _u=next_up(_u); }
}

//FloatMPTemplate<Bound>::FloatMPTemplate(double _v) : _l(_v), _u(_v) {
//}

//FloatMPTemplate<Bound>::FloatMPTemplate(double l, double u) : _l(l), _u(u) {
//}

FloatMPTemplate<Bound>::FloatMPTemplate(FloatMP _v) : _l(_v), _u(_v) {
}

FloatMPTemplate<Bound>::FloatMPTemplate(FloatMP l, FloatMP u) : _l(l), _u(u) {
}

FloatMPTemplate<Bound>::operator FloatMPTemplate<Metrc>() const {
    FloatMP v=half(add_near(_l,_u));
    FloatMP e=max(sub_up(v,_l),sub_up(_u,v));
    return FloatMPTemplate<Metrc>(v,e);
}

FloatMPTemplate<Bound>::operator FloatMPTemplate<Upper>() const {
    return FloatMPTemplate<Upper>(_u);
}

FloatMPTemplate<Bound>::operator FloatMPTemplate<Lower>() const {
    return FloatMPTemplate<Lower>(_l);
}

FloatMPTemplate<Bound>::operator FloatMPTemplate<Apprx>() const {
    return FloatMPTemplate<Apprx>(add_near(_l,_u));
}

FloatMPTemplate<Lower> const& FloatMPTemplate<Bound>::lower() const {
    return reinterpret_cast<FloatMPTemplate<Lower>const&>(this->_l);
}

FloatMPTemplate<Upper> const& FloatMPTemplate<Bound>::upper() const {
    return reinterpret_cast<FloatMPTemplate<Upper>const&>(this->_u);
}

FloatMPTemplate<Exact> FloatMPTemplate<Bound>::value() const {
    return FloatMPTemplate<Exact>(half(add_near(_l,_u)));
}

FloatMPTemplate<Error> FloatMPTemplate<Bound>::error() const {
    return FloatMPTemplate<Error>(half(sub_up(_u,_l)));
}

Void FloatMPTemplate<Bound>::set_precision(PrecisionMP pr) {
    _l.set_precision(pr);
    _u.set_precision(pr);
}

PrecisionMP FloatMPTemplate<Bound>::precision() const {
    return _l.precision();
}

BoundFloatMP pos(BoundFloatMP x) {
    return BoundFloatMP(x._l,x._u);
}

BoundFloatMP neg(BoundFloatMP x) {
    return BoundFloatMP(-x._u,-x._l);
}

BoundFloatMP rec(BoundFloatMP x) {
    FloatMP one(1,x.precision());
    if(x._l>0 or x._u<0) {
        return BoundFloatMP{div_down(one,x._u),div_up(one,x._l)};
    } else {
        FloatMP inf(1e300*1e300);
        return BoundFloatMP{-inf,+inf};
    }
}

BoundFloatMP sqr(BoundFloatMP x) {
    if(x._l>0) {
        return BoundFloatMP{mul_down(x._l,x._l),mul_up(x._u,x._u)};
    } else if(x._u<0) {
        return BoundFloatMP{mul_down(x._u,x._u),mul_up(x._l,x._l)};
    } else {
        FloatMP zero(0,x.precision());
        return BoundFloatMP{zero,std::max(mul_up(x._l,x._l),mul_up(x._u,x._u))};
    }
}


BoundFloatMP pow(BoundFloatMP x, Nat m) {
    ARIADNE_NOT_IMPLEMENTED;
}

BoundFloatMP pow(BoundFloatMP x, Int n) {
//    if(n<0) { return pow(rec(x),Nat(-n)); }
    if(n<0) { return rec(pow(x,Nat(-n))); }
    else return pow(x,Nat(n));
}


BoundFloatMP sqrt(BoundFloatMP x) {
    return BoundFloatMP(Ariadne::sqrt(x._l,MPFR_RNDD),Ariadne::sqrt(x._u,MPFR_RNDU));
}
BoundFloatMP exp(BoundFloatMP x) {
    return BoundFloatMP(Ariadne::exp(x._l,MPFR_RNDD),Ariadne::exp(x._u,MPFR_RNDU));
}

BoundFloatMP log(BoundFloatMP x) {
    return BoundFloatMP(Ariadne::log(x._l,MPFR_RNDD),Ariadne::log(x._u,MPFR_RNDU));
}

BoundFloatMP const_pi(PrecisionMP pr) {
    return BoundFloatMP(pi_down(pr),pi_up(pr));
}

BoundFloatMP sin(BoundFloatMP x) {
    return cos(x-pi(x.precision())/ExactFloatMP(2.0));
}

BoundFloatMP cos(BoundFloatMP x) {

    BoundFloatMP pi_bnd=pi(x.precision());
    auto pi_down=pi_bnd.lower().get_flt();
    auto pi_up=pi_bnd.lower().get_flt();
    auto pi_approx=pi_bnd.value().get_flt();

    ARIADNE_ASSERT(x._l<=x._u);

    FloatMP r_l(x.precision());
    FloatMP r_u(x.precision());

    if(x.error().get_flt()>2.0*pi_down) {
        r_l=-1.0; r_u=+1.0; return BoundedFloatMP(r_l,r_u);
    }

    FloatMP n=floor(x._l/(2*pi_approx)+0.5);
    x=x-ExactFloatMP(2*n)*pi_bnd;
    ARIADNE_ASSERT(x._l<=x._u);

    ARIADNE_ASSERT(x._l<=pi_up);
    ARIADNE_ASSERT(x._u>=-pi_up);

    if(x._l<=-pi_down) {
        if(x._u<=0.0) { r_l=-1.0; r_u=cos_up(x._u); }
        else { r_l=-1.0; r_u=+1.0; }
    } else if(x._l<=0.0) {
        if(x._u<=0.0) { r_l=cos_down(x._l); r_u=cos_up(x._u); }
        else if(x._u<=pi_down) { r_l=cos_down(max(-x._l,x._u)); r_u=+1.0; }
        else { r_l=-1.0; r_u=+1.0; }
    } else if(x._l<=pi_up) {
        if(x._u<=pi_down) { r_l=cos_down(x._u); r_u=cos_up(x._l); }
        else if(x._u<=2*pi_down) { r_l=-1.0; r_u=cos_up(min(x._l,sub_down(2*pi_down,x._u))); }
        else { r_l=-1.0; r_u=+1.0; }
    } else {
        assert(false);
    }
    return BoundedFloatMP(r_l,r_u);
}

BoundFloatMP tan(BoundFloatMP x) {
    return div(sin(x),cos(x));
}

BoundFloatMP atan(BoundFloatMP x) {
    return BoundFloatMP(Ariadne::atan(x._l,MPFR_RNDD),Ariadne::atan(x._u,MPFR_RNDU));
}


BoundFloatMP min(BoundFloatMP x, BoundFloatMP y) {
    using std::min;
    return BoundFloatMP(min(x._l,y._l),min(x._u,y._u));
}

BoundFloatMP max(BoundFloatMP x, BoundFloatMP y) {
    using std::max;
    return BoundFloatMP(max(x._l,y._l),max(x._u,y._u));
}

BoundFloatMP abs(BoundFloatMP x) {
    using std::max;
    FloatMP zero(x.precision());
    return BoundFloatMP(max(max(x._l,-x._u),zero),max(-x._l,x._u));
}


FloatMPTemplate<Bound> add(FloatMPTemplate<Bound> x,FloatMPTemplate<Bound> y) {
    return FloatMPTemplate<Bound>(add_down(x._l,y._l),add_up(x._u,y._u));
}

FloatMPTemplate<Bound> sub(FloatMPTemplate<Bound> x,FloatMPTemplate<Bound> y) {
    return FloatMPTemplate<Bound>(sub_down(x._l,y._u),sub_up(x._u,y._l));
}

FloatMPTemplate<Bound> mul(FloatMPTemplate<Bound> x,FloatMPTemplate<Bound> y) {
    if(y._l>0) {
        if(x._l>0) {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._l),mul_up(x._u,y._u)};
        } else if(x._u<0) {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._u),mul_up(x._u,y._l)};
        } else {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._u),mul_up(x._u,y._u)};
        }
    } else if(y._u<0) {
        if(x._l>0) {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._u),mul_up(x._u,y._l)};
        } else if(x._u<0) {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._l),mul_up(x._u,y._u)};
        } else {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._l),mul_up(x._u,y._l)};
        }
    } else {
        if(x._l>0) {
            return FloatMPTemplate<Bound>{mul_down(x._u,y._l),mul_up(x._u,y._l)};
        } else if(x._u<0) {
            return FloatMPTemplate<Bound>{mul_down(x._l,y._l),mul_up(x._u,y._u)};
        } else {
            return FloatMPTemplate<Bound>{std::min(mul_down(x._l,y._u),mul_down(x._u,y._l)),
                              std::max(mul_up(x._l,y._l),mul_up(x._u,y._u))};
        }
        return FloatMPTemplate<Bound>();
    }
}

FloatMPTemplate<Bound> div(FloatMPTemplate<Bound> x,FloatMPTemplate<Bound> y) {
    if(y._l>0) {
        if(x._l>0) {
            return FloatMPTemplate<Bound>{div_down(x._l,y._u),div_up(x._u,y._l)};
        } else if(x._u<0) {
            return FloatMPTemplate<Bound>{div_down(x._l,y._l),div_up(x._u,y._u)};
        } else {
            return FloatMPTemplate<Bound>{div_down(x._l,y._l),div_up(x._u,y._u)};
        }
    } else if(y._u<0) {
        if(x._l>0) {
            return FloatMPTemplate<Bound>{div_down(x._l,y._l),div_up(x._u,y._u)};
        } else if(x._u<0) {
            return FloatMPTemplate<Bound>{div_down(x._l,y._u),div_up(x._u,y._l)};
        } else {
            return FloatMPTemplate<Bound>{div_down(x._l,y._u),div_up(x._u,y._l)};
        }
    } else {
        assert(y._l>0 || y._u<0);
        return FloatMPTemplate<Bound>();
    }
}

Bool operator==(FloatMPTemplate<Bound> x, Int n) {
    return x._l==n && x._u==n;
}

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Bound> const& x) {
    return os << "{" << x.lower() << ":" << x.upper() << "}";
}



/************ FloatMPTemplate<Upper> ***************************************************/

FloatMPTemplate<Upper>::FloatMPTemplate() : _u(0.0) { }
FloatMPTemplate<Upper>::FloatMPTemplate(FloatMP f) : _u(f) { }
FloatMP const& FloatMPTemplate<Upper>::get_flt() const { return this->_u; }

FloatMPTemplate<Upper> pos(FloatMPTemplate<Upper> x) { return +x._u; }
FloatMPTemplate<Lower> neg(FloatMPTemplate<Upper> x) { return -x._u; }
FloatMPTemplate<Upper> add(FloatMPTemplate<Upper> x,FloatMPTemplate<Upper> y) { return add_near(x._u,y._u); }
FloatMPTemplate<Upper> sub(FloatMPTemplate<Upper> x,FloatMPTemplate<Lower> y) { return sub_down(x._u,y._l); }
FloatMPTemplate<Upper> mul(FloatMPTemplate<Upper> x,FloatMPTemplate<Upper> y) { return mul_near(x._u,y._u); }
FloatMPTemplate<Upper> div(FloatMPTemplate<Upper> x,FloatMPTemplate<Lower> y) { return div_down(x._u,y._l); }

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Upper> const& x) {
    char str[1024];
    mpfr_snprintf (str, 1024, "%.RUe", x._u._mpfr);
    return os << str;
}


/************ FloatMPTemplate<Lower> ***************************************************/

FloatMPTemplate<Lower>::FloatMPTemplate() : _l(0.0) { }
FloatMPTemplate<Lower>::FloatMPTemplate(FloatMP f) : _l(f) { }
FloatMP const& FloatMPTemplate<Lower>::get_flt() const { return this->_l; }

FloatMPTemplate<Lower> pos(FloatMPTemplate<Lower> x) { return +x._l; }
FloatMPTemplate<Upper> neg(FloatMPTemplate<Lower> x) { return -x._l; }
FloatMPTemplate<Lower> add(FloatMPTemplate<Lower> x,FloatMPTemplate<Lower> y) { return add_near(x._l,y._l); }
FloatMPTemplate<Lower> sub(FloatMPTemplate<Lower> x,FloatMPTemplate<Upper> y) { return sub_down(x._l,y._u); }
FloatMPTemplate<Lower> mul(FloatMPTemplate<Lower> x,FloatMPTemplate<Lower> y) { return mul_near(x._l,y._l); }
FloatMPTemplate<Lower> div(FloatMPTemplate<Lower> x,FloatMPTemplate<Upper> y) { return div_down(x._l,y._u); }

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Lower> const& x) {
    char str[1024];
    mpfr_snprintf (str, 1024, "%.RUe", x._l._mpfr);
    return os << str;
}


/************ FloatMPTemplate<Apprx> ***************************************************/

FloatMPTemplate<Apprx>::FloatMPTemplate() : _a(0.0) { }
FloatMPTemplate<Apprx>::FloatMPTemplate(double d) : _a(d) { }
FloatMPTemplate<Apprx>::FloatMPTemplate(FloatMP f) : _a(f) { }
Void FloatMPTemplate<Apprx>::set_precision(PrecisionMP pr) { _a.set_precision(pr); }
PrecisionMP FloatMPTemplate<Apprx>::precision() const { return _a.precision(); }
FloatMP const& FloatMPTemplate<Apprx>::get_flt() const { return this->_a; }

FloatMPTemplate<Apprx>::FloatMPTemplate(Integer const& z, PrecisionMP pr) : FloatMPTemplate(Rational(z),pr) { }
FloatMPTemplate<Apprx>::FloatMPTemplate(Rational const& q, PrecisionMP pr) : FloatMPTemplate(FloatMP(q,pr,MPFR_RNDN)) { }

FloatMPTemplate<Apprx> nul(FloatMPTemplate<Apprx> x);
FloatMPTemplate<Apprx> pos(FloatMPTemplate<Apprx> x) { return +x._a; }
FloatMPTemplate<Apprx> sqr(FloatMPTemplate<Apprx> x) { return x._a*x._a; }
FloatMPTemplate<Apprx> neg(FloatMPTemplate<Apprx> x) { return -x._a; }
FloatMPTemplate<Apprx> rec(FloatMPTemplate<Apprx> x) { return rec(x._a); }
FloatMPTemplate<Apprx> add(FloatMPTemplate<Apprx> x,FloatMPTemplate<Apprx> y) { return add_near(x._a,y._a); }
FloatMPTemplate<Apprx> sub(FloatMPTemplate<Apprx> x,FloatMPTemplate<Apprx> y) { return sub_near(x._a,y._a); }
FloatMPTemplate<Apprx> mul(FloatMPTemplate<Apprx> x,FloatMPTemplate<Apprx> y) { return mul_near(x._a,y._a); }
FloatMPTemplate<Apprx> div(FloatMPTemplate<Apprx> x,FloatMPTemplate<Apprx> y) { return div_near(x._a,y._a); }
FloatMPTemplate<Apprx> abs(FloatMPTemplate<Apprx> x) { return abs(x._a); }
FloatMPTemplate<Apprx> max(FloatMPTemplate<Apprx> x,FloatMPTemplate<Apprx> y) { return max(x._a,y._a); }
FloatMPTemplate<Apprx> min(FloatMPTemplate<Apprx> x,FloatMPTemplate<Apprx> y) { return min(x._a,y._a); }
FloatMPTemplate<Apprx> pow(FloatMPTemplate<Apprx> x, Int n) { return pow_near(x.get_flt(),n); }
FloatMPTemplate<Apprx> sqrt(FloatMPTemplate<Apprx> x) { return sqrt_near(x.get_flt()); }
FloatMPTemplate<Apprx> exp(FloatMPTemplate<Apprx> x) { return exp_near(x.get_flt()); }
FloatMPTemplate<Apprx> log(FloatMPTemplate<Apprx> x) { return log_near(x.get_flt()); }
FloatMPTemplate<Apprx> sin(FloatMPTemplate<Apprx> x) { return sin_near(x.get_flt()); }
FloatMPTemplate<Apprx> cos(FloatMPTemplate<Apprx> x) { return cos_near(x.get_flt()); }
FloatMPTemplate<Apprx> tan(FloatMPTemplate<Apprx> x) { return tan_near(x.get_flt()); }
FloatMPTemplate<Apprx> atan(FloatMPTemplate<Apprx> x) { return atan_near(x.get_flt()); }

OutputStream& operator<<(OutputStream& os, FloatMPTemplate<Apprx> const& x) {
    char str[1024];
    mpfr_snprintf (str, 1024, "%.RNe", x._a._mpfr);
    return os << str;
}



/************ Class name **********************************************/

template<> String class_name<FloatMP>() { return "FloatMPTemplate"; }
template<> String class_name<ExactFloatMP>() { return "ExactFloatMP"; }
template<> String class_name<ErrorFloatMP>() { return "ErrorFloatMP"; }
template<> String class_name<MetrcFloatMP>() { return "MetrcatedFloatMP"; }
template<> String class_name<BoundFloatMP>() { return "BoundedFloatMP"; }
template<> String class_name<LowerFloatMP>() { return "LowerFloatMP"; }
template<> String class_name<UpperFloatMP>() { return "UpperFloatMP"; }
template<> String class_name<ApprxFloatMP>() { return "ApproximateFloatMP"; }



} // namespace Ariadne
