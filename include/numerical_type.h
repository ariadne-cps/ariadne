/***************************************************************************
 *            numerical_type.h
 *
 *  Thu Oct 02 16:27:05 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numerical_type.h
 *  \brief Type definitions and conversion operators for fundamental Ariadne types.
 */

#ifndef _NUMERICAL_TYPE_H
#define _NUMERICAL_TYPE_H

#include <gmpxx.h>
#include <string>

namespace Ariadne {

  /*! \brief An integer
  *
  * An element of the ring of integers.
  * Must allow denotation of any integer, including arbitrarily large values.
  * Integer quotient and remainder must be supported.
  *
  * Currently implemented using mpz_class from the GNU Multiple Precision Library.
  */
  typedef mpz_class Integer;

  /*! \brief A dyadic rational (i.e. of form @a m/2^n).
  *
  * An element of the ring of dyadic rationals.
  * Must allow denotation of any dyadic rational.
  * May be created without loss of precision from any integral or floating point type,
  * or from any rational of the form m/2^n.
  *
  * Currently implemented using mpf_class from the GNU Multiple Precision library.
  */
  typedef mpf_class Dyadic;

  /*! \brief A rational number.
  *
  * An element of the field of rationals.
  * Must allow denotation of any rational.
  * May be created without loss of precision from any integral or floating point type, and from a dyadic.
  *
  * Currently implemented using mpq_class from the GNU Multiple Precision library.
  */
  typedef mpq_class Rational;


  inline Rational transform_into_rational(
        const Rational &num){ return num; }

  inline Integer numerator(const Rational 
      &num){ return num.get_num(); }

  inline Integer denumerator(const Rational
      &num){ return num.get_den();}


  inline Rational epsilon(const Rational &a) {

    Rational e(1,1000);

    return e*e*e;
  }


  template <typename NumType>
  NumType gcd(const NumType &a, const NumType &b){

    NumType c=a%b;

    if (c==0) return b;

    return (gcd(b, c));
  }

  template <typename NumType>
  inline NumType lcm(const NumType &a, const NumType &b) {
    return ((a*b)/gcd(a,b));
  }

  template <typename T> std::string name();

  template<> std::string name<double>() { return "double"; }
  template<> std::string name<Rational>() { return "Rational"; }
  template<> std::string name<Dyadic>() { return "Dyadic"; }

}

#endif
