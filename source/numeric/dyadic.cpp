/***************************************************************************
 *            dyadic.cpp
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file dyadic.cpp
 *  \brief
 */



#include "utility/stdlib.hpp"

#include "dyadic.hpp"

#include "utility/macros.hpp"
#include "utility/string.hpp"
#include "numeric/logical.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/builtin.hpp"
#include "numeric/rational.hpp"

#include <limits>

namespace Ariadne {

static const mp_bitcnt_t maximum_precision = 65535;

Dyadic::~Dyadic() {
    mpf_clear(_mpf);
}

Dyadic::Dyadic() {
    mpf_init2(_mpf,maximum_precision);
}

Dyadic::Dyadic(mpf_t mpf) {
    mpf_init2(_mpf,maximum_precision);
    mpf_set(_mpf,_mpf);
}

Dyadic::Dyadic(Integer const& p, Nat q) {
    mpf_init2(_mpf,maximum_precision);
    mpf_set_z(_mpf,p._mpz);
    mpf_div_2exp(_mpf,_mpf,q);
    //if(q>=0) { mpf_div_2exp(_mpf,_mpf,q); } else { mpf_mul_2exp(_mpf,_mpf,-q); }
}

Dyadic::Dyadic(Dbl d) : Dyadic(ExactDouble(d)) {
}

Dyadic::Dyadic(ExactDouble const& x) {
    static bool give_inf_warning=true;
    Dbl d=x.get_d();
    if(std::isinf(d)) {
        if(give_inf_warning) {
            std::cerr<<"WARNING: Converting Double inf to Dyadic is not supported by GMP; returning numeric_limits<double>::max()\n";
            give_inf_warning=false;
        }
        const double max=std::numeric_limits<double>::max();
        d=(d>0?+max:-max);
    }
    ARIADNE_ASSERT(std::isfinite(d));
    mpf_init2(_mpf,maximum_precision);
    mpf_set_d(_mpf,d);
}

Dyadic::Dyadic(Integer const& z) {
    mpf_init2(_mpf,maximum_precision);
    mpf_set_z(_mpf,z._mpz);
}

Dyadic::Dyadic(TwoExp const& w) : Dyadic(1u) {
    const int q=w.exponent();
    if(q>=0) { mpf_mul_2exp(_mpf,_mpf,q); }
    else { mpf_div_2exp(_mpf,_mpf,-q); }
}

Dyadic::Dyadic(const Dyadic& x) {
    mpf_init2(_mpf,maximum_precision);
    mpf_set(_mpf,x._mpf);
}

Dyadic::Dyadic(Dyadic&& x) {
    mpf_init(_mpf);
    mpf_swap(_mpf,x._mpf);
}

Dyadic& Dyadic::operator=(const Dyadic& x) {
    mpf_set(_mpf,x._mpf);
    return *this;
}

Dyadic& Dyadic::operator=(Dyadic&& x) {
    mpf_swap(_mpf,x._mpf);
    return *this;
}

Integer Dyadic::mantissa() const {
    return Rational(*this).numerator();
}

Int Dyadic::exponent() const {
    return log2floor(Rational(*this).denominator());
}

mpf_t const& Dyadic::get_mpf() const {
    return this->_mpf;
}

double Dyadic::get_d() const {
    return mpf_get_d(this->_mpf);
}

Dyadic operator+(TwoExp y) {
    return +Dyadic(y); 
}

Dyadic operator-(TwoExp y) {
    return -Dyadic(y);
}


/*
Dyadic& operator+=(Dyadic& x1, Dyadic const& x2) {
    mpf_add(x1._mpf,x1._mpf,x2._mpf);
    return x1;
}

Dyadic& operator-=(Dyadic& x1, Dyadic const& x2) {
    mpf_sub(x1._mpf,x1._mpf,x2._mpf);
    return x1;
}

Dyadic& operator*=(Dyadic& x1, Dyadic const& x2) {
    mpf_mul(x1._mpf,x1._mpf,x2._mpf);
    return x1;
}
*/

Integer round(Dyadic const& x) {
    Integer r;
    mpz_set_f(r._mpz,x._mpf);
    return r;
}


Dyadic nul(Dyadic const& x) {
    Dyadic r;
    mpf_set_si(r._mpf,0);
    return r;
}

Dyadic pos(Dyadic const& x) {
    Dyadic r;
    mpf_set(r._mpf,x._mpf);
    return r;
}

Dyadic neg(Dyadic const& x) {
    Dyadic r;
    mpf_neg(r._mpf,x._mpf);
    return r;
}

Dyadic sqr(Dyadic const& x) {
    return x*x;
}

Dyadic add(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r;
    mpf_add(r._mpf,x1._mpf,x2._mpf);
    return r;
}

Dyadic sub(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r;
    mpf_sub(r._mpf,x1._mpf,x2._mpf);
    return r;
}

Dyadic mul(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r;
    mpf_mul(r._mpf,x1._mpf,x2._mpf);
    return r;
}

Dyadic hlf(Dyadic const& x) {
    Dyadic r;
    mpf_div_2exp(r._mpf,x._mpf,1);
    return r;
}

Dyadic pow(Dyadic const& x, Nat m) {
    unsigned long int lm=m;
    Dyadic r;
    mpf_pow_ui(r._mpf,x._mpf,lm);
    return r;
}


Dyadic abs(Dyadic const& x) {
    Dyadic r;
    mpf_abs(r._mpf,x._mpf);
    return r;
}

Dyadic min(Dyadic const& x1,Dyadic const& x2) {
    return (x1<x2)?x1:x2;
}

Dyadic max(Dyadic const& x1,Dyadic const& x2) {
    return (x1>x2)?x1:x2;
}


Comparison cmp(Dyadic const& x1, Dyadic const& x2) {
    auto c=mpf_cmp(x1._mpf,x2._mpf);
    return c==0 ? Comparison::EQUAL : (c>0?Comparison::GREATER:Comparison::LESS);
}

Boolean eq(Dyadic const& x1, Dyadic const& x2) {
    return mpf_cmp(x1._mpf,x2._mpf)==0;
}

Boolean lt(Dyadic const& x1, Dyadic const& x2) {
    return mpf_cmp(x1._mpf,x2._mpf)<0;
}


//   mpf_get_str (char *str, mp_exp_t *expptr, int base, size_t n_digits, const mpf_t op)
OutputStream& operator<<(OutputStream& os, Dyadic const& x) {
    Rational q;
    mpq_set_f (q._mpq,x._mpf);
    os << q.numerator();
    Int exp = log2floor(q.denominator());
    if (exp!=0) { if(exp==1) { os << "/2"; } else { os << "/2^" << exp; } }
    return os;
}

Dyadic make_dyadic(unsigned long long int n) {
    static const unsigned int max=std::numeric_limits<Int>::max();
    static const unsigned long long int m=max;
    unsigned long long int q = n / m;
    unsigned long long int r = n % m;
    unsigned int rem = static_cast<unsigned int>(r);
    assert(n==q*m+r);
    if(q==0) {
        return Dyadic(rem);
    } else {
        return make_dyadic(q)*Dyadic(max)+Dyadic(rem);
    }
}

template<> String class_name<Dyadic>() { return "Dyadic"; }

} // namespace Ariadne
