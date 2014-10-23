/***************************************************************************
 *            numeric/rational.h
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

/*! \file numeric/rational.h
 *  \brief
 */



#ifndef ARIADNE_RATIONAL_H
#define ARIADNE_RATIONAL_H

#include "external/gmp.h"
#include "utility/typedefs.h"
#include "utility/metaprogramming.h"
#include "utility/string.h"
#include "numeric/integer.h"

namespace Ariadne {

struct Exact;

class Boolean;

class Int64;

class Integer;
template<class P> class Number;
class Float;
class ExactFloat;

template<class X> struct IsNumber;

enum class Comparison : char;


/************ Rational *******************************************************/

//! \ingroup UserNumberSubModule
//! \brief %Rational numbers.
class Rational
{
  public:
    mpq_t _mpq;
  public:
    typedef Exact Paradigm;
  public:
    ~Rational();
    Rational();
    Rational(const Integer&, const Integer&);
    template<class N, EnableIf<IsIntegral<N>> = dummy> Rational(N n);
    Rational(Int64);
    explicit Rational(Float);
    explicit Rational(Float64);
    Rational(const Integer&);
    explicit Rational(const String&);
    explicit Rational(const ExactFloat&);
    explicit Rational(const ExactFloat64&);
    explicit Rational(const mpq_t);
    Rational(const Rational&);
    Rational(Rational&&);
    Rational& operator=(const Rational&);
    Rational& operator=(Rational&&);
    operator Number<Exact> () const;
    Integer get_num() const;
    Integer get_den() const;
    friend Rational operator+(Rational const& q);
    friend Rational operator-(Rational const& q);
    friend Rational operator+(Rational const& q1, Rational const& q2);
    friend Rational operator-(Rational const& q1, Rational const& q2);
    friend Rational operator*(Rational const& q1, Rational const& q2);
    friend Rational operator/(Rational const& q1, Rational const& q2);
    friend Rational operator/(Integer const& z1, Integer const& z2);
    friend Rational& operator+=(Rational& q1, Rational const& q2);
    friend Rational& operator-=(Rational& q1, Rational const& q2);
    friend Rational& operator*=(Rational& q1, Rational const& q2);
    friend Rational& operator/=(Rational& q1, Rational const& q2);
    friend Rational max(Rational const& q1, Rational const& q2);
    friend Rational min(Rational const& q1, Rational const& q2);
    friend Rational abs(Rational const& q);
    friend Rational pos(Rational const& q);
    friend Rational neg(Rational const& q);
    friend Rational sqr(Rational const& q);
    friend Rational rec(Rational const& q);
    friend Rational add(Rational const& q1, Rational const& q2);
    friend Rational sub(Rational const& q1, Rational const& q2);
    friend Rational mul(Rational const& q1, Rational const& q2);
    friend Rational div(Rational const& q1, Rational const& q2);
    friend Rational pow(Rational const& q, Nat m);
    friend Rational pow(Rational const& q, Int n);

    friend Boolean eq(Rational const& q1, Rational const& q2);
    friend Comparison cmp(Rational const& q1, Rational const& q2);

    friend Boolean operator==(Rational const& q1, Rational const& q2);
    friend Boolean operator!=(Rational const& q1, Rational const& q2);
    friend Boolean operator>=(Rational const& q1, Rational const& q2);
    friend Boolean operator<=(Rational const& q1, Rational const& q2);
    friend Boolean operator> (Rational const& q1, Rational const& q2);
    friend Boolean operator< (Rational const& q1, Rational const& q2);

    friend OutputStream& operator<<(OutputStream& os, Rational const& q);
    friend InputStream& operator>>(InputStream& os, Rational& q);
    friend Rational operator"" _q(long double x);
  public:
    double get_d() const;
    mpq_t const& get_mpq() const;
  private:
    friend class Dyadic;
  public:
    explicit Rational(double);
};
template<> struct IsNumber<Rational> : True { };
Rational operator"" _q(long double x);

template<class N, EnableIf<IsIntegral<N>>> inline Rational::Rational(N n) : Rational(Int64(n)) { }


} // namespace Ariadne

#endif
