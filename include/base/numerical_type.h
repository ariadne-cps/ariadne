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
  
  /*! \brief A dyadic rational (i.e. of form @a m/2^n).
  *
  * An element of the ring of dyadic rationals.
  * Must allow denotation of any dyadic rational.
  * May be created without loss of precision from any integral or floating point type,
  * or from any rational of the form m/2^n.
  * May be created without loss of precision from any integral or floating point type,
  * or from any rational of the form m/2^n.
  *
  * Currently implemented using mpf_class from the GNU Multiple Precision library.
  *
  * FIXME: mpf_class does not implement addition, subtraction and multiplication exactly.
  */
  typedef mpf_class Dyadic;
  //typedef Synaps::Dyadic Dyadic;
  
/* Hand-coded wrapped (currently not used)
  class Dyadic {
   public:
    Dyadic(const mpf_class& x) : _rep(x) { }
    Dyadic(const double& x) : _rep(x) { }
    Dyadic(const int& x) : _rep(x) { }
    Dyadic(const Dyadic& x) : _rep(x._rep) { }
    
    Dyadic(const double& x, uint n) : _rep(x,n) { } 
    Dyadic(const Dyadic& x, uint n) : _rep(x._rep,n) { }

    uint precision() const { return _rep.get_prec(); }
    int exponent() const { long int res; mpf_get_d_2exp(&res,_rep.get_mpf_t()); return res-1; }
    Dyadic mantissa() const { mpf_class res; mpf_div_2exp(res.get_mpf_t(),_rep.get_mpf_t(),this->exponent()); return res; }
    double mantissa_approx() const { long int tmp; return mpf_get_d_2exp(&tmp,_rep.get_mpf_t()); }
    
    Dyadic operator-() const { Dyadic res; mpf_neg(res._rep.get_mpf_t(),this->_rep.get_mpf_t()); return res; }
    Dyadic operator+(const Dyadic& y) const { Dyadic res; mpf_add(res._rep.get_mpf_t(),this->_rep.get_mpf_t(),y._rep.get_mpf_t()); return res; }
    Dyadic operator-(const Dyadic& y) const { Dyadic res; mpf_sub(res._rep.get_mpf_t(),this->_rep.get_mpf_t(),y._rep.get_mpf_t()); return res; }
    Dyadic operator*(const Dyadic& y) const { Dyadic res; mpf_mul(res._rep.get_mpf_t(),this->_rep.get_mpf_t(),y._rep.get_mpf_t()); return res; }
    
    bool operator==(const Dyadic& y) const { return *this==y; }
    bool operator!=(const Dyadic& y) const { return *this!=y; }
    bool operator<(const Dyadic& y) const { return *this<y; }
    bool operator<=(const Dyadic& y) const { return *this<=y; }
    
    Dyadic operator-() const { return mpf_class(-this->_rep); }
    Dyadic operator+(const Dyadic& y) const { return mpf_class(this->_rep + y._rep); }
    Dyadic operator-(const Dyadic& y) const { return mpf_class(this->_rep - y._rep); }
    Dyadic operator*(const Dyadic& y) const { return mpf_class(this->_rep * y._rep); }
    
    friend Dyadic div_approx(const Dyadic& x, const Dyadic& y, const Dyadic& e);
    friend Dyadic div_approx(const Dyadic& x, const Dyadic& y, const uint& n);
    friend std::ostream& operator<<(std::ostream&, const Dyadic&);
    friend std::istream& operator>>(std::istream&, Dyadic&);
    Dyadic() { }
   private:
    mpf_class _rep;
  };
 
  inline std::ostream& operator<<(std::ostream& os, const Dyadic& x) {
    mpf_class tmp(x._rep);
    os << tmp;
    return os;
  }
  
  inline std::istream& operator>>(std::istream& is, Dyadic& x) {
    mpf_class tmp;
    is >> tmp;
    x=Dyadic(tmp);
    return is;
  }
*/
  
  /*! \brief A rational number.
  *
  * An element of the field of rationals.
  * Must allow denotation of any rational.
  * May be created without loss of precision from any integral or floating point type, and from a dyadic.
  *
  * Currently implemented using mpq_class from the GNU Multiple Precision library.
  */
  typedef mpq_class Rational;

  inline Integer precision(const Dyadic& num) {
    return mpf_get_prec(num.get_mpf_t());
  }

  inline Integer exponent(const Dyadic& num) {
    long int res; 
    mpf_get_d_2exp(&res,num.get_mpf_t()); 
    return res-1; 
  }

  inline Dyadic mantissa(const Dyadic& num) {
    long int exp; 
    mpf_class res; 
    mpf_get_d_2exp(&exp,num.get_mpf_t()); 
    exp=exp-1;
    mpf_div_2exp(res.get_mpf_t(),num.get_mpf_t(),exp); 
    return res; 
  }

  inline Integer numerator(const Rational& num){ 
    return num.get_num(); }

  inline Integer denominator(const Rational& num){ 
    return num.get_den();}

  inline Integer floor(const Rational& num){ 
    return numerator(num)/denominator(num);
  }

  /* numerical traits */
  class ring_tag { };
  class field_tag { };
    
  template<typename T> class numerical_traits;
    
  template<> class numerical_traits<double> {
   public:
    typedef field_tag algebraic_category;
    typedef double field_extension_type;
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

  template<typename Res, typename Arg> inline Res convert_to(const Arg& x) { return Res(x); }
  
  template<> inline Rational convert_to(const Rational& q) { return q; }
  template<> inline Rational convert_to(const Dyadic& d) { return Rational(d); }
  template<> inline Rational convert_to(const double& x) { return Rational(x); }
  
  template<> inline double convert_to(const Rational& x) { return x.get_d(); }
  template<> inline double convert_to(const Dyadic& x) { return x.get_d(); }

  template<> inline int convert_to(const Integer& n) { return n.get_si(); }
  template<> inline long convert_to(const Integer& n) { return n.get_si(); }
  
  template<typename Res, typename Arg>
  Res approximate(const Arg& x, const int& p);

  
  inline Rational approximate(const Rational& q, const Rational& e) {
    //Rational r=floor(q);
    return q;
  }
  
  /* names (for diagnostic output) */
  template <typename T> std::string name();
  template<> inline std::string name<double>() { return "double"; }
  template<> inline std::string name<Rational>() { return "Rational"; }
  template<> inline std::string name<Dyadic>() { return "Dyadic"; }
    
}

#endif /* _ARIADNE_NUMERICAL_TYPE */
