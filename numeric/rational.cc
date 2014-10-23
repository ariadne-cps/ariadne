/***************************************************************************
 *            rational.cc
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

/*! \file rational.cc
 *  \brief
 */



#include "utility/stdlib.h"
#include "utility/macros.h"

#include "rational.h"
#include "logical.h"
#include "integer.h"
#include "sign.h"
#include <limits>

#include <limits>
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace Ariadne {

class InvalidRationalLiteralException {
  public:
    InvalidRationalLiteralException(std::string what) { }
};

// Shortened version of raw float classes sufficient for comparison operator
struct Float { double dbl; };
struct Flt64 { volatile double d; };

Rational rec(Integer const& z) {
    return Rational(1,z);
}

Rational pow(Integer const& z, Int n) {
    if(n>=0) { return Rational(pow(z,Nat(n))); } else { return Rational(1,pow(z,Nat(-n))); }
}

Rational::~Rational() {
    mpq_clear(_mpq);
}

Rational::Rational() {
    mpq_init(_mpq);
    mpq_set_si(_mpq,0,1);
}

/*
Rational::Rational(uint m) {
    mpq_init(_mpq);
    mpq_set_ui(_mpq,m,1u);
}

Rational::Rational(int n) {
    mpq_init(_mpq);
    mpq_set_si(_mpq,n,1);
}
*/

Rational::Rational(Integer const& z) {
    mpq_init(_mpq);
    mpq_set_z(_mpq,z._mpz);
}

Rational::Rational(Integer const& znum, Integer const& zden) {
    mpq_init(_mpq);
    mpq_set_num(_mpq,znum._mpz);
    mpq_set_den(_mpq,zden._mpz);
    mpq_canonicalize(_mpq);
}

Rational::Rational(Int64 n) : Rational(Integer(n)) {
}

Rational::Rational(double d) {
    mpq_init(_mpq);
    mpq_set_d(_mpq,d);
    mpq_canonicalize(_mpq);
}

Rational::Rational(Float x) : Rational(x.dbl) {
}

Rational::Rational(Float64 x) : Rational(x.d) {
}

Rational::Rational(ExactFloat const& x) : Rational(reinterpret_cast<Float const&>(x)) {
}

Rational::Rational(ExactFloat64 const& x) : Rational(reinterpret_cast<Float64 const&>(x)) {
}

Rational::Rational(const String& s) {
    int base=10;
    mpq_init(_mpq);
    mpq_set_str(_mpq,s.c_str(),base);
    mpq_canonicalize(_mpq);
}

Rational::Rational(const mpq_t q) {
    mpq_init(_mpq);
    mpq_set(_mpq,q);
    mpq_canonicalize(_mpq);
}

Rational::Rational(const Rational& q) {
    mpq_init(_mpq);
    mpq_set(_mpq,q._mpq);
}

Rational::Rational(Rational&& q) {
    mpq_init(_mpq);
    mpq_swap(_mpq,q._mpq);
}

Rational& Rational::operator=(const Rational& q) {
    mpq_set(_mpq,q._mpq);
    return *this;
}

Rational& Rational::operator=(Rational&& q) {
    mpq_swap(_mpq,q._mpq);
    return *this;
}


mpq_t const& Rational::get_mpq() const {
    return this->_mpq;
}

double Rational::get_d() const {
    return mpq_get_d(this->_mpq);
}

Integer Rational::get_num() const {
    Integer z;
    mpq_get_num(z._mpz,this->_mpq);
    return z;
}

Integer Rational::get_den() const {
    Integer z;
    mpq_get_den(z._mpz,this->_mpq);
    return z;
}

Rational operator+(Rational const& q) {
    return Rational(q);
}

Rational operator-(Rational const& q) {
    Rational r; mpq_neg(r._mpq,q._mpq);
    return std::move(r);
}

Rational operator+(Rational const& q1, Rational const& q2) {
    Rational r; mpq_add(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational operator-(Rational const& q1, Rational const& q2) {
    Rational r; mpq_sub(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational operator*(Rational const& q1, Rational const& q2) {
    Rational r; mpq_mul(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational operator/(Rational const& q1, Rational const& q2) {
    Rational r; mpq_div(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational operator/(Integer const& z1, Integer const& z2) {
    return Rational(z1,z2);
}

Rational& operator+=(Rational& q1, Rational const& q2) {
    mpq_add(q1._mpq,q1._mpq,q2._mpq);
    return q1;
}

Rational& operator-=(Rational& q1, Rational const& q2) {
    mpq_sub(q1._mpq,q1._mpq,q2._mpq);
    return q1;
}

Rational& operator*=(Rational& q1, Rational const& q2) {
    mpq_mul(q1._mpq,q1._mpq,q2._mpq);
    return q1;
}

Rational& operator/=(Rational& q1, Rational const& q2) {
    mpq_div(q1._mpq,q1._mpq,q2._mpq);
    return q1;
}

Rational max(Rational const& q1, Rational const& q2) {
    return q1>q2 ? q1 :  q2;
}

Rational min(Rational const& q1, Rational const& q2) {
    return q1<q2 ? q1 :  q2;
}

Rational abs(Rational const& q) {
    Rational r; mpq_abs(r._mpq,q._mpq);
    return std::move(r);
}

Rational pos(Rational const& q) {
    Rational r; mpq_set(r._mpq,q._mpq);
    return std::move(r);
}

Rational neg(Rational const& q) {
    Rational r; mpq_neg(r._mpq,q._mpq);
    return std::move(r);
}

Rational sqr(Rational const& q) {
    Rational r; mpq_mul(r._mpq,q._mpq,q._mpq);
    return std::move(r);
}

Rational rec(Rational const& q) {
    Rational r; mpq_inv(r._mpq,q._mpq);
    return std::move(r);
}

Rational add(Rational const& q1, Rational const& q2) {
    Rational r; mpq_add(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational sub(Rational const& q1, Rational const& q2) {
    Rational r; mpq_sub(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational mul(Rational const& q1, Rational const& q2) {
    Rational r; mpq_mul(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational div(Rational const& q1, Rational const& q2) {
    Rational r; mpq_div(r._mpq,q1._mpq,q2._mpq);
    return std::move(r);
}

Rational div(Integer const& z1, Integer const& z2) {
    return Rational(z1,z2);
}

Rational pow(Rational const& q, Nat m) {
    Rational r=1; Rational p=q;
    while(m!=0) { if(m%2==1) { r=r*p; } p=p*p; m/=2; }
    return std::move(r);
}

Rational pow(Rational const& q, Int n) {
    if(n<0) { return rec(pow(q,Nat(-n))); }
    else { return pow(q,Nat(n)); }
}

Boolean eq(Rational const& q1, Rational const& q2) {
    return Boolean(mpq_equal(q1._mpq,q2._mpq));
}

Comparison cmp(Rational const& q1, Rational const& q2) {
    return Comparison(mpq_cmp(q1._mpq,q2._mpq));
}


Rational operator"" _q(long double x) {
    static const uint32_t max_cf_coef = std::numeric_limits<uint32_t>::max();
    static const std::size_t N=11;

    // See if the long double value is an exact single-precision value
    volatile float sx=static_cast<float>(x);
    if(static_cast<long double>(sx)==x) {
        return Rational(static_cast<double>(sx));
    }
    // Compute the continued fraction expansion of x, storing coefficients in cf
    // Stop if a coefficient is larger than max_cf_coef, since then the result
    //   is accurate to within 2^{-32}
    long double t=x;
    int s=1;
    if(x<0) { s=-1; t=-x; }
    long double cf[N+1];
    cf[N]=0;
    std::size_t i;
    for(i=0; i!=N; ++i) {
        cf[i]=std::floor(t);
        t=1/(t-cf[i]);
        if(t>max_cf_coef) { break; }
    }
    if(i==N) {
        ARIADNE_THROW(InvalidRationalLiteralException,"Rational operator"" _q(long double)",
                      "x="<<x<<" is not a sufficiently close approximation to a simple rational number.");
    }
    // Compute the result from the continued fraction coefficients
    Rational q = Rational(static_cast<double>(cf[i]));
    while(i!=0) {
        --i;
        q=Rational(cf[i])+rec(q);
    }
    if(s==-1) { q=-q; }
    double xd=static_cast<double>(x);
    Rational xq=Rational(xd);
    //volatile double qd=q.get_d();
    double ae=std::abs((q-xq).get_d());
    double re=ae/std::max(1.0,std::abs(xd));
    if(re>std::numeric_limits<double>::epsilon()) {
        ARIADNE_THROW(InvalidRationalLiteralException,"Rational operator"" _q(long double)",
                      "Rational approximation q="<<q<<" to x="<<x<<"="<<xd<<" has error "<<ae<<" and relative error "<<re<<" while is larger than machine epsilon");
    }
    mpq_canonicalize(q._mpq);
    return q;
}

OutputStream& write(OutputStream& os, mpz_t const z) {
    char str[255];
    str[254]='\0';
    mpz_get_str (str, 10, z);
    assert(str[254]=='\0');
    return os << str;
}

InputStream& operator>>(InputStream& is, Rational& q1) {
    ARIADNE_NOT_IMPLEMENTED;
}

//   mpq_get_str (char *str, mpq_eq1p_t *eq1pptr, int b, size_t n, mpq_t op, mpq_rnd_t rnd)
// If str is not a null pointer, it should point to a block of storage large enough for the significand,
// i.e., at least maq1(n + 2, 7). The eq1tra two bytes are for a possible minus sign,
// and for the terminating null character, and the value 7 accounts for -@Inf@ plus the terminating null character.
OutputStream& operator<<(OutputStream& os, Rational const& q1) {
    mpz_t num; mpz_t den;
    mpz_init(num); mpz_init(den);
    mpq_get_num(num,q1._mpq);
    mpq_get_den(den,q1._mpq);
    write(os,num);
    if(mpz_cmp_si(den,1)!=0) { os << "/"; write(os,den); }
    mpz_clear(num); mpz_clear(den);
    return os;
}

Boolean operator==(Rational const& q1, Rational const& q2) { return cmp(q1,q2)==Comparison::EQUAL; }
Boolean operator!=(Rational const& q1, Rational const& q2) { return cmp(q1,q2)!=Comparison::EQUAL; }
Boolean operator<=(Rational const& q1, Rational const& q2) { return cmp(q1,q2)!=Comparison::GREATER; }
Boolean operator>=(Rational const& q1, Rational const& q2) { return cmp(q1,q2)!=Comparison::LESS; }
Boolean operator< (Rational const& q1, Rational const& q2) { return cmp(q1,q2)==Comparison::LESS; }
Boolean operator> (Rational const& q1, Rational const& q2) { return cmp(q1,q2)==Comparison::GREATER; }

Boolean eq(Rational const& q1, Flt64 const& x2) {
    assert(false);
    //return eq(q1,Rational(x2.get_d()));
}

Comparison cmp(Rational const& q1, Flt64 const& x2) {
    assert(false);
    //return cmp(q1,Rational(x2.get_d()));
}


template<> String class_name<Rational>() { return "Rational"; }

} // namespace Ariadne
