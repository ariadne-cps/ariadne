/***************************************************************************
 *            numeric/integer.h
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

/*! \file numeric/integer.h
 *  \brief
 */



#ifndef ARIADNE_INTEGER_H
#define ARIADNE_INTEGER_H

#include "external/gmp.h"

#include <cassert>

#include "utility/typedefs.h"
#include "utility/metaprogramming.h"
#include "numeric/sign.h"
#include "numeric/logical.h"
#include "numeric/number.decl.h"

namespace Ariadne {

/************  Ints ********************************************************/

template<class X> struct IsNumber;

class Boolean;
class Integer;
class Rational;

template<class P> class Number;

class Nat32 {
    uint32_t _n;
  public:
    template<class N, EnableIf<IsIntegral<N>> = dummy> Nat32(N n) : _n(n) { assert(_n==n); }
    uint32_t get_ui() const { return _n; }
};

class Nat64 {
    uint64_t _n;
  public:
    template<class N, EnableIf<IsIntegral<N>> = dummy> Nat64(N n) : _n(n) { assert(_n==n); }
    uint64_t get_ui() const { return _n; }
};

class Int32 {
    int32_t _n;
  public:
    template<class N, EnableIf<IsIntegral<N>> = dummy> Int32(N n) : _n(n) { assert(_n==n); }
    int32_t get_si() const { return _n; }
};

class Int64 {
    int64_t _n;
  public:
    Int64() : _n(0) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> Int64(N n) : _n(n) { assert(_n==n); }
    int64_t get_si() const { return _n; }
};

uint32_t fac(uint8_t n);
uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint32_t bin(uint8_t n, uint8_t k);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);

struct Exact;
class Integer;
template<> struct IsNumber<Nat> : True { };
template<> struct IsNumber<Int> : True { };
template<> struct IsNumber<Integer> : True { };

//! \ingroup UserNumberSubModule
//! \brief Arbitrarily-sized integers.
class Integer
{
  public:
    mpz_t _mpz;
  public:
    typedef Exact Paradigm;
  public:
    ~Integer();
    Integer();
    template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>> = dummy> Integer(M m);
    template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>> = dummy> Integer(N n);
    explicit Integer(const mpz_t);
    Integer(const Integer&);
    Integer(Integer&&);
    Integer& operator=(const Integer&);
    Integer& operator=(Integer&&);
    operator Number<Exact> () const;

    friend Integer operator+(Integer const& z);
    friend Integer operator-(Integer const& z);
    friend Integer operator+(Integer const& z1, Integer const& z2);
    friend Integer operator-(Integer const& z1, Integer const& z2);
    friend Integer operator*(Integer const& z1, Integer const& z2);
    friend Rational operator/(Integer const& z1, Integer const& z2);
    friend Integer& operator++(Integer& z);
    friend Integer& operator--(Integer& z);
    friend Integer& operator+=(Integer& z1, Integer const& z2);
    friend Integer& operator*=(Integer& z1, Integer const& z2);
    friend Integer const& max(Integer const& z1, Integer const& z2);
    friend Integer const& min(Integer const& z1, Integer const& z2);
    friend Integer abs(Integer const& z);
    friend Integer pos(Integer const& z);
    friend Integer neg(Integer const& z);
    friend Integer sqr(Integer const& z);
    friend Integer add(Integer const& z1, Integer const& z2);
    friend Integer sub(Integer const& z1, Integer const& z2);
    friend Integer mul(Integer const& z1, Integer const& z2);
    friend Integer pow(Integer const& z, Nat m);

    friend Rational rec(Integer const& z);
    friend Rational div(Integer const& z1, Integer const& z2);
    friend Rational pow(Integer const& z, Int n);

    friend Boolean operator==(Integer const& z1, Integer const& z2);
    friend Boolean operator!=(Integer const& z1, Integer const& z2);
    friend Boolean operator>=(Integer const& z1, Integer const& z2);
    friend Boolean operator<=(Integer const& z1, Integer const& z2);
    friend Boolean operator> (Integer const& z1, Integer const& z2);
    friend Boolean operator< (Integer const& z1, Integer const& z2);

    friend OutputStream& operator<<(OutputStream& os, Integer const& z);
    friend Integer operator"" _z(unsigned long long int n);
  public:
    long int get_si() const;
    mpz_t const& get_mpz() const;
  private:
  public:
    Integer(Nat32 m);
    Integer(Int32 n);
    Integer(Nat64 m);
    Integer(Int64 n);
};

template<class M, EnableIf<And<IsIntegral<M>,IsUnsigned<M>>>> inline Integer::Integer(M m) : Integer(Nat64(m)) { }
template<class N, EnableIf<And<IsIntegral<N>,IsSigned<N>>>> inline Integer::Integer(N n) : Integer(Int64(n)) { }
Integer operator"" _z(unsigned long long int n);

// Comparisons with arbitary ints go through Int64
template<class N, EnableIf<IsIntegral<N>> = dummy> inline auto operator==(Integer const& x, N n) -> decltype(x==Int64(n)) { return x==Int64(n); }
template<class N, EnableIf<IsIntegral<N>> = dummy> inline auto operator!=(Integer const& x, N n) -> decltype(x!=Int64(n)) { return x!=Int64(n); }
template<class N, EnableIf<IsIntegral<N>> = dummy> inline auto operator< (Integer const& x, N n) -> decltype(x!=Int64(n)) { return x< Int64(n); }
template<class N, EnableIf<IsIntegral<N>> = dummy> inline auto operator> (Integer const& x, N n) -> decltype(x!=Int64(n)) { return x> Int64(n); }
template<class N, EnableIf<IsIntegral<N>> = dummy> inline auto operator<=(Integer const& x, N n) -> decltype(x!=Int64(n)) { return x<=Int64(n); }
template<class N, EnableIf<IsIntegral<N>> = dummy> inline auto operator>=(Integer const& x, N n) -> decltype(x!=Int64(n)) { return x>=Int64(n); }


} // namespace Ariadne

#endif
