/***************************************************************************
 *            integer.class.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file numeric/integer.class.h
 *  \brief Definition of integer class. 
 */

#ifndef ARIADNE_NUMERIC_INTEGER_CLASS_H
#define ARIADNE_NUMERIC_INTEGER_CLASS_H

#include <gmp.h>

namespace Ariadne {
  
  

    template<class E> class Expression;
    template<class X> class Value;
  
    /*!\ingroup Numeric
     * \brief An integer of arbitrary size.
     *  
     * An element of the ring of integers.
     * Must allow denotation of any integer, including arbitrarily large values.
     * Integer quotient and remainder must be supported.
     *
     * Currently implemented using mpz_t from the GNU Multiple Precision Library.
     */
    class Integer
      : public Value<Integer>
    { 
     public:
      mpz_t _value; 
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs the integer 0. */
      ~Integer();
      /*! \brief Default constructor constructs the integer 0. */
      Integer();
      /*! \brief Copy constructor. */
      Integer(const Integer& z);
      /*! \brief Copy assignment operator. */
      Integer& operator=(const Integer& z);
      
      /*! \brief Construct from a string literal. */
      explicit Integer(const std::string& n);
      /*! \brief Convert from a built-in integer. */
      Integer(const int& n);
      Integer(const long int& n);
      /*! \brief Convert from an unsigned built-in integer. */
      Integer(const unsigned int& n);
      Integer(const unsigned long int& n);
      /*! \brief Conversion assignment operator from a built-in integer. */
      Integer& operator=(const int& n);
      Integer& operator=(const long int& n);
      /*! \brief Conversion assignment operator from an unsigned built-in  integer. */
      Integer& operator=(const unsigned int& n);
      Integer& operator=(const unsigned long int& n);

      /*! \brief Convert from a GMP integer. */
      Integer(mpz_srcptr z);
      /*! \brief Conversion assignment operator from an GMP integer. */
      Integer& operator=(mpz_srcptr z);

      /*! \brief Convert from a numerical expression. */
      template<class E> Integer(const Expression<E>& e);
      /*! \brief Assign from a numerical expression. */
      template<class E> Integer& operator=(const Expression<E>& e);
      //@}

      //@{
      //! \name Accessor methods
      /*! \brief A reference to the internal value. */
      mpz_ptr value() { return this->_value; }
      /*! \brief A constant reference to the internal value. */
      mpz_srcptr value() const { return this->_value; }
      /*! \brief Convert to an built-in integer value. A check is performed to ensure that the conversion is be performed exactly. */
      operator long int() const;
      //@}

#ifdef DOXYGEN
      //@{
      //! \name Arithmetic operations
      /*! \brief The minimum of z1 and z2. */
      friend Integer min(const Integer& z1, const Integer& z2);
      /*! \brief The maximum of z1 and z2. */
      friend Integer max(const Integer& z1, const Integer& z2);
      /*! \brief The absolute value \a z. */
      friend Integer abs(const Integer& z);
      
      /*! \brief Unary plus. */
      friend Integer pos(const Integer& z);
      /*! \brief Negation. */
      friend Integer neg(const Integer& z);
      /*! \brief Addition. */
      friend Integer add(const Integer& z1, const Integer& z2);
      /*! \brief Subtraction. */
      friend Integer sub(const Integer& z1, const Integer& z2);
      /*! \brief Multiplication. */
      friend Integer mul(const Integer& z1, const Integer& z2);
      /*! \brief Multiplication. */
      friend Rational div(const Integer& z1, const Integer& z2);

      /*! \brief Quotient. */
      friend Integer quot(const Integer& z1, const Integer& z2);
      /*! \brief Remainder. */
      friend Integer rem(const Integer& z1, const Integer& z2);
      /*! \brief Power by a positive integer. */
      friend Integer pow(const Integer& z, const unsigned int& n);

      /*! \brief In-place increment. */
      friend Integer& operator++(Integer& z);
      /*! \brief In-place decrement. */
      friend Integer& operator--(Integer& z);
      /*! \brief In-place addition. */
      friend Integer& operator+=(Integer& z1, const Integer& z2);
      /*! \brief In-place subtraction of a number. */
      friend Integer& operator-=(Integer& z1, const Integer& z2);
      /*! \brief In-place multiplication. */
      friend Integer& operator*=(Integer& z1, const Integer& z2);

      /*! \brief Unary plus. */
      friend Integer operator+(const Integer& z);
      /*! \brief Negation. */
      friend Integer operator-(const Integer& z);
      /*! \brief Addition. */
      friend Integer operator+(const Integer& z1, const Integer& z2);
      /*! \brief Subtraction. */
      friend Integer operator-(const Integer& z1, const Integer& z2);
      /*! \brief Multiplication. */
      friend Integer operator*(const Integer& z1, const Integer& z2);
      /*! \brief Division. */
      friend Rational operator/(const Integer& z1, const Integer& z2);
      /*! \brief Remainder. */
      friend Integer operator%(const Integer& z1, const Integer& z2);
      //@}
      
      
      //@{
      //! \name Comparison operators.
      /*! \brief Equality operator. */
      friend bool operator==(const Integer& z1, const Integer& z2); 
      /*! \brief Inezuality operator. */
      friend bool operator!=(const Integer& z1, const Integer& z2); 
      /*! \brief Less than operator. */
      friend bool operator<(const Integer& z1, const Integer& z2);  
      /*! \brief Greater than operator. */
      friend bool operator>(const Integer& z1, const Integer& z2);
      /*! \brief Less than or ezual to operator. */
      friend bool operator<=(const Integer& z1, const Integer& z2);
      /*! \brief Greater than or ezual to operator. */
      friend bool operator>=(const Integer& z1, const Integer& z2);
      //@}

      //@{
      //! \name %Integer functions
      /*! \brief Factorial. */
      friend Integer fac(const Integer& n);

      /*! \brief The number of ways of choosing \a k objects from \a n. */
      friend Integer bin(const Integer& n, const Integer& k);

      /*! \brief Greatest common divisor. */
      friend Integer gcd(const Integer& n1, const Integer& n2);

      /*! \brief Least common multiple. */
      friend Integer lcm(const Integer& n1, const Integer& n2);
      //@}

      //! \name %Integer bit-shift functions
      //@{
      /*! \brief The integer power \f$2^n\f$. */
      friend Integer exp2(const Integer& n);

      /*! \brief The floor of the logarithm of \a n in base 2. */
      friend Integer log2_floor(const Integer& n);

      /*! \brief The ceiling of the logarithm of \a n in base 2. */
      friend Integer log2_ceil(const Integer& n);
      //@}

      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const Integer& z);
      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, Integer& z);
      //@}
#endif
    };

  
}

#endif /* ARIADNE_NUMERIC_INTEGER_CLASS_H */
