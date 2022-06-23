/***************************************************************************
 *            numeric/decimal.hpp
 *
 *  Copyright  2014-20  Pieter Collins
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

/*! \file numeric/decimal.hpp
 *  \brief ExactTag decimal numbers, useful for user input.
 */

#ifndef ARIADNE_DECIMAL_HPP
#define ARIADNE_DECIMAL_HPP

#include <string>
#include <iostream>
#include "utility/typedefs.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/integer.hpp"

namespace Ariadne {

//! \ingroup NumericModule
//! \sa Dyadic, Rational, Real
//! \brief A decimal number.
class Decimal
    : public DefineComparisonOperators<Decimal,Boolean,Boolean>
{
  public:
    static const Integer _ten;
    Integer _p; Nat _q;
  public:
    typedef ExactTag Paradigm;
    //! \brief Default constructor creates the number 0 (zero).
    Decimal();
    //! \brief Construct the number p/10^q.
    explicit Decimal(Integer p, Nat q);
    explicit Decimal(Integer p, Int q) = delete;
    //! \brief Construct from a builtin integral type.
    template<BuiltinIntegral N> Decimal(N n);
    //! \brief Construct from a double-precision floating-point number representation.
    explicit Decimal(double d);
    //! \brief Construct from a string representation.
    explicit Decimal(String const&);
    //! \brief Convert from a dyadic.
    Decimal(Dyadic const&);
    //! \brief Convert to a rational number.
    operator Rational () const;
    //! \brief Convert to an generic number.
    operator ExactNumber () const;
    //! \brief A string literal
    String literal() const;
    //! \brief Unary plus of a decimal value.
    friend Decimal operator+(Decimal const& d);
    //! \brief Negation of a decimal value.
    friend Decimal operator-(Decimal const& d);
    //! \brief Addition of two decimal values.
    friend Decimal operator+(Decimal const& d1, Decimal const& d2);
    //! \brief Subtraction of two decimal values.
    friend Decimal operator-(Decimal const& d1, Decimal const& d2);
    //! \brief Multiplication of two decimal values.
    friend Decimal operator*(Decimal const& d1, Decimal const& d2);
    //! \brief Division of two decimal values yields a rational.
    friend Rational operator/(Decimal const& d1, Decimal const& d2);
    //! \brief Inplace addition of two decimal values.
    friend Decimal& operator+=(Decimal& d1, Decimal const& d2);
    //! \brief Squared value of a decimal.
    friend Decimal nul(Decimal const& d);
    //! \brief Zeroed value of a decimal.
    friend Decimal sqr(Decimal const& d);
    //! \brief Half the value of a decimal.
    friend Decimal hlf(Decimal const& d);
    //! \brief Absolute value of a decimal.
    friend Decimal abs(Decimal const& d);
    //! \brief Maximum of two decimal values.
    friend Decimal max(Decimal const& d1, Decimal const& d2);
    //! \brief Minimum of two decimal values.
    friend Decimal min(Decimal const& d1, Decimal const& d2);
    //! \brief Comparison of two decimal values.
    friend Comparison cmp(Decimal const& d1, Decimal const& d2);
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, Decimal const& d);
    //! \brief Construct from an integer literal.
    friend Decimal operator"" _decimal (unsigned long long int n);
    //! \brief Construct from a floating-point literal.
    friend Decimal operator"" _decimal (long double dbl);
    //! \brief Construct from a string literal.
    friend Decimal operator"" _decimal (const char* str, std::size_t);

    //! \brief Shorthand for operator""_decimal.
    friend Decimal operator"" _dec (long double dbl);

    //! \brief Alternative for operator""_decimal for use in Python interface.
    friend Decimal dec_(long double dbl);

    void canonicalize();
};
Decimal operator"" _dec (unsigned long long int n);
Decimal operator"" _dec (long double dbl);
Decimal operator"" _dec (const char* str, std::size_t);
Decimal operator"" _decimal (unsigned long long int n);
Decimal operator"" _decimal (long double dbl);
Decimal operator"" _decimal (const char* str, std::size_t);

template<BuiltinIntegral N> inline Decimal::Decimal(N n) : Decimal(n,0u) { }

template<> class Bounds<Decimal> {
    Decimal _l, _u;
  public:
    Bounds(Decimal w) : _l(w), _u(w) { }
    template<ConvertibleTo<Decimal> X> Bounds(X const& x)
        : Bounds(Decimal(x)) { }
    Bounds(Decimal l, Decimal u) : _l(l), _u(u) { }
    template<class X> requires Constructible<Decimal,X> Bounds(Bounds<X> const& x)
        : Bounds(Decimal(x.lower_raw()),Decimal(x.upper_raw())) { }
    operator ValidatedNumber() const;
    Decimal lower() const { return _l; }
    Decimal upper() const { return _u; }
    Decimal lower_raw() const { return _l; }
    Decimal upper_raw() const { return _u; }
    friend OutputStream& operator<<(OutputStream& os, Bounds<Decimal> y) { return os << "[" << y._l << ":" << y._u << "]"; }
};


} // namespace Ariadne

#endif
