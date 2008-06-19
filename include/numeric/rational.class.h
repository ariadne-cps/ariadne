/***************************************************************************
 *            numeric/rational.class.h
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
 
/*! \file numeric/rational.class.h
 *  \brief Definition of rational number class.
 */

#ifndef ARIADNE_NUMERIC_RATIONAL_CLASS_H
#define ARIADNE_NUMERIC_RATIONAL_CLASS_H

#include <gmp.h>
#include <gmpxx.h>


namespace Ariadne {
  
  
    template<class E> class Expression;
    template<class X> class Value;

    class Integer;
    class Rational;
    template<class T> class Float;
  
 
   /*!\ingroup Numeric
    * \brief A rational number.
    * 
    * An element of the field of rationals.
    * Must allow denotation of any rational.
    * May be created without loss of precision from any integral or floating point type, and from a dyadic.
    *
    * Currently implemented using mpq_class from the GNU Multiple Precision library.
    */
    class Rational
      : public Value<Rational>
    {
     public:
      mpq_t _value;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Destructor. */
      ~Rational();
      /*! \brief Default constructor constructs the rational 0/1. */
      Rational();
#ifdef DOXYGEN
      /*! \brief Construct from a numerator and a denomin_ator. */
      template<class N1, class N2> Rational(const N1& n, const N2& d);
#else
      Rational(const int& n, const int& d);
      Rational(const Integer& n, const Integer& d);
#endif

      /*! \brief Convert from an int. */
      Rational(const int& n);
      Rational(const long int& n);
      /*! \brief Convert from an unsigned int. */
      Rational(const unsigned int& n);
      Rational(const unsigned long int& n);
      /*! \brief Convert from a double. */
      Rational(const double& x);

      /*! \brief Convert from an integer. */
      Rational(const Integer& z);
      /*! \brief Convert from a Float. */
      template<class T> Rational(const Float<T>& x);
      /*! \brief Copy constructor. */
      Rational(const Rational& q);

      /*! \brief Convert from a raw mpz_class. */
      Rational(const mpz_class& z);
      /*! \brief Convert from a raw mpf_class. */
      Rational(const mpf_class& x);
      /*! \brief Convert from a raw mpq_class. */
      Rational(const mpq_class& q);

      /*! \brief Construct from a string literal. */
      explicit Rational(const std::string& n);

      /*! \brief Assign from an int. */
      Rational& operator=(const int& n);
      Rational& operator=(const long int& n);
      /*! \brief Conversion assignment operator from anunsigned built-in integer. */
      Rational& operator=(const unsigned int& n);
      Rational& operator=(const unsigned long int& n);
      /*! \brief Conversion assignment operator from a double. */
      Rational& operator=(const double& x);
      /*! \brief Assign from an integer. */
      Rational& operator=(const Integer& q);
      /*! \brief Assign from an float. */
      template<class T> Rational& operator=(const Float<T>& x);
      /*! \brief Copy assignment operator. */
      Rational& operator=(const Rational& q);

      /*! \brief Convert from a numerical expression. */
      template<class E> Rational(const Expression<E>& e);
      /*! \brief Assign from a numerical expression. */
      template<class E> Rational& operator=(const Expression<E>& e);

      // Convert from a numerical expression with a rounding mode. (For convenience only). */
      template<class E, class Rnd> Rational(const Expression<E>& e, Rnd);
      template<class X, class Rnd> Rational(const X& e, Rnd);

      //@}
      
      //@{
      //! \name Accessor methods
      /*! \brief A reference to the internal value. */
      mpq_ptr value() { return this->_value; }
      /*! \brief A constant reference to the internal value. */
      mpq_srcptr value() const { return this->_value; }
      //@}

      //@{
      //! \name Data access
      /*! \brief The numerator. */
      Integer numerator() const;
      /*! \brief The denomin_ator. */
      Integer denominator() const;
      
      double get_d() const;
      //@}
     public:
      void canonicalize();

#ifdef DOXYGEN
      //@{
      //! \name Arithmetic operations
      /*! \brief The min_imum of q1 and q2. */
      friend Rational min_(const Rational& q1, const Rational& q2);
      /*! \brief The max_imum of q1 and q2. */
      friend Rational max_(const Rational& q1, const Rational& q2);
      /*! \brief The abs_olute value \a q. */
      friend Rational abs_(const Rational& q);
      
      /*! \brief In-place add_ition. */
      friend Rational& operator+=(Rational& q1, const Rational& q2);
      /*! \brief In-place sub_traction of a number. */
      friend Rational& operator-=(Rational& q1, const Rational& q2);
      /*! \brief In-place mul_tiplication. */
      friend Rational& operator*=(Rational& q1, const Rational& q2);
      /*! \brief In-place div_ision. */
      friend Rational& operator/=(Rational& q1, const Rational& q2);

      /*! \brief Negation. */
      friend Rational operator-(const Rational& q);
      /*! \brief Addition. */
      friend Rational operator+(const Rational& q1, const Rational& q2);
      /*! \brief Subtraction. */
      friend Rational operator-(const Rational& q1, const Rational& q2);
      /*! \brief Multiplication. */
      friend Rational operator*(const Rational& q1, const Rational& q2);
      /*! \brief Division. */
      friend Rational operator/(const Rational& q1, const Rational& q2);
      /*! \brief %Integer pow_er. */
      friend Rational pow_(const Rational& q, const Integer& n);
      //@}
      
      
      //@{
      //! \name Comparison operators.
      /*! \brief Equality operator. */
      friend bool operator==(const Rational& q1, const Rational& q2); 
      /*! \brief Inequality operator. */
      friend bool operator!=(const Rational& q1, const Rational& q2); 
      /*! \brief Less than operator. */
      friend bool operator<(const Rational& q1, const Rational& q2);  
      /*! \brief Greater than operator. */
      friend bool operator>(const Rational& q1, const Rational& q2);
      /*! \brief Less than or equal to operator. */
      friend bool operator<=(const Rational& q1, const Rational& q2);
      /*! \brief Greater than or equal to operator. */
      friend bool operator>=(const Rational& q1, const Rational& q2);
      //@}

      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const Rational& q);
      /*! \brief Stream extraction operator. */
      friend std::istream& operator>>(std::istream& is, Rational& q);
      //@}
#endif
    };

    // Declare stream i/o operators
    std::ostream& operator<<(std::ostream& os, const Rational& q);
    std::istream& operator>>(std::istream& is, Rational& q);
      
  
} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_RATIONAL_CLASS_H */
