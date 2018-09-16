/***************************************************************************
 *            numeric/rational.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file numeric/rational.hpp
 *  \brief
 */



#ifndef ARIADNE_RATIONAL_HPP
#define ARIADNE_RATIONAL_HPP

#include "../external/gmp.hpp"
#include "../utility/typedefs.hpp"
#include "../utility/metaprogramming.hpp"
#include "../utility/string.hpp"
#include "../numeric/integer.hpp"
#include "../numeric/arithmetic.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"

namespace Ariadne {

class Int64;
class FloatDP;

enum class Comparison : char;


/************ Rational *******************************************************/

//! \ingroup NumericModule
//! \brief %Rational numbers.
class Rational
    : DeclareFieldOperations<Rational>
    , DeclareLatticeOperations<Rational,Rational>
    , DeclareComparisonOperations<Rational,Boolean,Boolean>
    , DefineFieldOperators<Rational>
    , DefineComparisonOperators<Rational,Boolean,Boolean>
    , DeclareTranscendentalOperations<Real>
{
  public:
    mpq_t _mpq;
  public:
    typedef ExactTag Paradigm;
    typedef Rational NumericType;
  public:
    ~Rational();
    Rational();
    Rational(const Integer&, const Integer&);
    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> Rational(N n);
    Rational(Int64);
    explicit Rational(FloatDP const&);
    Rational(const ExactDouble&);
    Rational(const Integer&);
    Rational(const Dyadic&);
    explicit Rational(const String&);
    explicit Rational(const FloatDPValue&);
    explicit Rational(const mpq_t);
    Rational(const Rational&);
    Rational(Rational&&);
    Rational& operator=(const Rational&);
    Rational& operator=(Rational&&);
    operator Number<ExactTag> () const;

    static Rational inf(Sign sgn);
    static Rational inf();
    static Rational nan();

    Integer get_num() const;
    Integer get_den() const;
    Integer numerator() const;
    Natural denominator() const;
    friend Rational operator/(Integer const& z1, Integer const& z2);

    friend Real sqrt(Real const&);
    friend Real exp(Real const&);
    friend Real log(Real const&);
    friend Real sin(Real const&);
    friend Real cos(Real const&);
    friend Real tan(Real const&);
    friend Real atan(Real const&);

    friend Sign sgn(Rational const& q);

    friend Bool is_nan(Rational const& q);
    friend Bool is_inf(Rational const& q);
    friend Bool is_finite(Rational const& q);
    friend Bool is_zero(Rational const& q);

    friend Comparison cmp(Rational const& q1, Rational const& q2);
    friend Comparison cmp(Rational const& q1, ExactDouble const& d2);
    friend Comparison cmp(ExactDouble const& d1, Rational const& q2);

    friend OutputStream& operator<<(OutputStream& os, Rational const& q);
    friend InputStream& operator>>(InputStream& os, Rational& q);
    friend Rational operator"" _q(long double x);
  public:
    double get_d() const;
    mpq_t const& get_mpq() const;
  private:
    friend class Dyadic;
};
template<> struct IsNumericType<Rational> : True { };
Rational operator"" _q(unsigned long long int n);
Rational operator"" _q(long double x);

template<class N, EnableIf<IsBuiltinIntegral<N>>> inline Rational::Rational(N n) : Rational(Int64(n)) { }

OutputStream& write(OutputStream& os, mpz_t const z);
InputStream& operator>>(InputStream& is, Rational& q1);


} // namespace Ariadne

#endif
