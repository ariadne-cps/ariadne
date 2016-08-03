/***************************************************************************
 *            dyadic.cc
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
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file dyadic.cc
 *  \brief
 */



#include "utility/stdlib.h"

#include "dyadic.h"

#include "utility/macros.h"
#include "utility/string.h"
#include "numeric/logical.h"
#include "numeric/rational.h"

#include <limits>

namespace Ariadne {

static const mp_bitcnt_t maximum_precision = 65535;

Dyadic::~Dyadic() {
    mpf_clear(_mpf);
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

Dyadic::Dyadic() : Dyadic(Integer(0)) {
}

Dyadic::Dyadic(Integer const& z) : Dyadic(z,1u) {
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

Dyadic::Dyadic(Dbl x) {
    mpf_init2(_mpf,maximum_precision);
    mpf_set_d(_mpf,x);
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

Dyadic nul(Dyadic const& x) {
    Dyadic r;
    mpf_set_si(r._mpf,0);
    return std::move(r);
}

Dyadic pos(Dyadic const& x) {
    Dyadic r;
    mpf_set(r._mpf,x._mpf);
    return std::move(r);
}

Dyadic neg(Dyadic const& x) {
    Dyadic r;
    mpf_neg(r._mpf,x._mpf);
    return std::move(r);
}

Dyadic sqr(Dyadic const& x) {
    return x*x;
}

Dyadic add(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r;
    mpf_add(r._mpf,x1._mpf,x2._mpf);
    return std::move(r);
}

Dyadic sub(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r;
    mpf_sub(r._mpf,x1._mpf,x2._mpf);
    return std::move(r);
}

Dyadic mul(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r;
    mpf_mul(r._mpf,x1._mpf,x2._mpf);
    return std::move(r);
}

Dyadic hlf(Dyadic const& x) {
    Dyadic r;
    mpf_div_2exp(r._mpf,x._mpf,1);
    return std::move(r);
}

Dyadic pow(Dyadic const& x, Nat m) {
    unsigned long int lm=m;
    Dyadic r;
    mpf_pow_ui(r._mpf,x._mpf,lm);
    return std::move(r);
}


Dyadic abs(Dyadic const& x) {
    Dyadic r;
    mpf_abs(r._mpf,x._mpf);
    return std::move(r);
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
