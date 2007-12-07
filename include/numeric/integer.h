/***************************************************************************
 *            numeric/integer.h
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
 
/*! \file numeric/integer.h
 *  \brief Multiple-precision integer type and interger functions.
 */

#ifndef ARIADNE_NUMERIC_INTEGER_H
#define ARIADNE_NUMERIC_INTEGER_H

//#include <gmpxx.h>
#include <gmp.h>
#include <iosfwd>
#include <stdexcept>

#include "numeric/traits.h"
#include "numeric/expression.h"
#include "numeric/operators.h"

namespace Ariadne {
  namespace Numeric {
  
    /*!\ingroup Numeric
     * \brief An integer of arbitrary size.
     *  
     * An element of the ring of integers.
     * Must allow denotation of any integer, including arbitrarily large values.
     * Integer quotient and remainder must be supported.
     *
     * Currently implemented using mpz_class from the GNU Multiple Precision Library.
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
      /*! \brief Convert from an unsigned built-in integer. */
      Integer(const unsigned int& n);
      /*! \brief Conversion assignment operator from a built-in integer. */
      Integer& operator=(const int& n);
      /*! \brief Conversion assignment operator from an unsigned built-in  integer. */
      Integer& operator=(const unsigned int& n);

      /*! \brief Convert from a GMP integer. */
      Integer(const mpz_class& z);
      /*! \brief Conversion assignment operator from an GMP integer. */
      Integer& operator=(const mpz_class& z);

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
      operator int() const;
      //@}

#ifdef DOXYGEN
      //@{
      //! \name Arithmetic operations
      /*! \brief The min_imum of z1 and z2. */
      friend Integer min_(const Integer& z1, const Integer& z2);
      /*! \brief The max_imum of z1 and z2. */
      friend Integer max_(const Integer& z1, const Integer& z2);
      /*! \brief The abs_olute value \a z. */
      friend Integer abs_(const Integer& z);
      
      /*! \brief Power by a pos_itive integer. */
      friend Integer pow_(const Integer& z, const unsigned int& n);
      /*! \brief Quotient. */
      friend Integer quot(const Integer& z1, const Integer& z2);
      /*! \brief Remainder. */
      friend Integer rem(const Integer& z1, const Integer& z2);

      /*! \brief In-place increment. */
      friend Integer& operator++(Integer& z);
      /*! \brief In-place decrement. */
      friend Integer& operator--(Integer& z);
      /*! \brief In-place add_ition. */
      friend Integer& operator+=(Integer& z1, const Integer& z2);
      /*! \brief In-place sub_traction of a number. */
      friend Integer& operator-=(Integer& z1, const Integer& z2);
      /*! \brief In-place mul_tiplication. */
      friend Integer& operator*=(Integer& z1, const Integer& z2);

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
      friend Integer factorial(const Integer& n);

      /*! \brief The number of ways of choosing \a k objects from \a n. */
      friend Integer choose(const Integer& n, const Integer& k);

      /*! \brief Greatest common div_isor. */
      friend Integer gcd(const Integer& n1, const Integer& n2);

      /*! \brief Least common mul_tiple. */
      friend Integer lcm(const Integer& n1, const Integer& n2);
      //@}

      //! \name %Integer bit-shift functions
      //@{
      /*! \brief The integer pow_er \f$2^n\f$. */
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

    // Declare stream i/o operators
    std::ostream& operator<<(std::ostream& os, const Integer& n);
    std::istream& operator>>(std::istream& is, Integer& n);

    // Operators on builtin types
    uint factorial(uint n);
    uint choose(uint n, uint k);
    int choose(int n, int k);


  }
}

#include "integer.inline.h"

#endif /* ARIADNE_NUMERIC_INTEGER_H */
