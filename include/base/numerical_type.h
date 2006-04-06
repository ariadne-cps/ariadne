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

#ifndef _ARIADNE_NUMERICAL_TYPE_H
#define _ARIADNE_NUMERICAL_TYPE_H

#include <gmpxx.h>
#include <string>
#include <assert.h>

#include "../base/dyadic.h"

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

  /*! \brief A 64-bit fixed-precision floating point number.
   *
   * Standard operations are not exact, but must support interval arithmetic.
   *
   * Currently implemented by the built-in type double.
   */
  typedef double Float64;
  
  /*! \brief A dyadic rational (i.e. of form \f$m/2^n\f$).
  *
  * An element of the ring of dyadic rationals.
  * Must allow denotation of any dyadic rational.
  * May be created without loss of precision from any integral or floating point type,
  * or from any rational of the form m/2^n.
  * May be created without loss of precision from any integral or floating point type,
  * or from any rational of the form m/2^n.
  *
  * Currently implemented using a modification of the Synaps dyadic class.
  *
  * FIXME: mpf_class does not implement addition, subtraction and multiplication exactly.
  */
  //typedef Synaps::dyadic Dyadic;
  
  /*! \brief A multiple-precision floating-point type.
   *
   * Currently implemented using mpf_class from the GNU Multiple Precision library.
   */
  typedef mpf_class MPFloat;
  
  /*! \brief A rational number.
  *
  * An element of the field of rationals.
  * Must allow denotation of any rational.
  * May be created without loss of precision from any integral or floating point type, and from a dyadic.
  *
  * Currently implemented using mpq_class from the GNU Multiple Precision library.
  */
  typedef mpq_class Rational;


  inline Integer precision(const MPFloat& num) {
    return mpf_get_prec(num.get_mpf_t());
  }

  inline Integer exponent(const MPFloat& num) {
    long int res; 
    mpf_get_d_2exp(&res,num.get_mpf_t()); 
    return res-1; 
  }

  inline MPFloat mantissa(const MPFloat& num) {
    long int exp; 
    mpf_class res; 
    mpf_get_d_2exp(&exp,num.get_mpf_t()); 
    exp=exp-1;
    mpf_div_2exp(res.get_mpf_t(),num.get_mpf_t(),exp); 
    return res; 
  }

  inline Dyadic mantissa(const Dyadic& num) {
    return num.mantissa();
  }

  inline int exponent(const Dyadic& num) {
    return num.exponent();
  }

  inline int precision(const Dyadic& num) {
    return num.precision();
  }

  inline Integer numerator(const Dyadic& num) {
    return num.numerator();
  }

  inline Integer denominator(const Dyadic& num) {
    return num.denominator();
  }


  inline Integer numerator(const Rational& num){ 
    return num.get_num(); }

  inline Integer denominator(const Rational& num){ 
    return num.get_den();}


  /* numerical traits */
  /*! \brief Tags a class representing a ring. */
  class ring_tag { };
  /*! \brief Tags a class representing a field. */
  class field_tag { };
    
  /*! \brief Typedef's describing a numerical type. */
  template<typename T> class numerical_traits;
    
  template<> class numerical_traits<double> {
   public:
    typedef field_tag algebraic_category;
    typedef double field_extension_type;
  };

  template<> class numerical_traits<MPFloat> {
   public:
    typedef ring_tag algebraic_category;
    typedef Rational field_extension_type;
  };
    
  template<> class numerical_traits<Dyadic> {
   public:
    typedef ring_tag algebraic_category;
    typedef Rational field_extension_type;
  };
    
  template<> class numerical_traits<Rational> {
   public:
    typedef field_tag algebraic_category;
    typedef Rational field_extension_type;
  };

  /* names (for diagnostic output) */
  template <typename T> std::string name();
  template<> inline std::string name<Float64>() { return "Float64"; }
  template<> inline std::string name<Rational>() { return "Rational"; }
  template<> inline std::string name<Dyadic>() { return "Dyadic"; }


  /*! \brief Convert from type \p Arg to type \p Res. */
  template<typename Res, typename Arg> Res convert_to(const Arg& x) { return Res(x); }
  
  template<> inline double convert_to(const Rational& x) { return x.get_d(); }
  template<> inline double convert_to(const MPFloat& x) { return x.get_d(); }
  template<> inline double convert_to(const Dyadic& x) { return x.get_d(); }

  template<> inline int convert_to(const Integer& n) { return n.get_si(); }
  template<> inline long convert_to(const Integer& n) { return n.get_si(); }
    
}

#endif /* _ARIADNE_NUMERICAL_TYPE */
