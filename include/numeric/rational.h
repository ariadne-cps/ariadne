/***************************************************************************
 *            rational.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file rational.h
 *  \brief Type definitions and conversion operators for rational numbers.
 */

#ifndef ARIADNE_RATIONAL_H
#define ARIADNE_RATIONAL_H

#include <gmpxx.h>
#include <iostream>

#include "../numeric/numerical_traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

#include "../numeric/integer.h"

namespace Ariadne {
  namespace Numeric {

   /*!\ingroup Numeric
    * \brief A rational number.
    * 
    * An element of the field of rationals.
    * Must allow denotation of any rational.
    * May be created without loss of precision from any integral or floating point type, and from a dyadic.
    *
    * Currently implemented using mpq_class from the GNU Multiple Precision library.
    */
    class Rational : public mpq_class 
    {
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs the rational 0/1. */
      Rational() : mpq_class() { }
      /*! \brief Construct from a numerator and a denominator. */
      Rational(const Integer& n, const Integer& d)
        : mpq_class(n,d) { this->mpq_class::canonicalize(); }
      /*! \brief Copy constructor. */
      Rational(const Rational& q)
        : mpq_class(q) { this->mpq_class::canonicalize(); }
      /*! \brief Copy assignment operator. */
      Rational& operator=(const Rational& q) {
        this->mpq_class::operator=(q); return *this; }

      /*! \brief Construct from a numerator and a denominator of arbitrary types. */
      template<class R1,class R2> Rational(const R1& n, const R2& d)
        : mpq_class(n,d) { this->mpq_class::canonicalize(); }
      /*! \brief Convert from another numerical type. */
      template<class R> Rational(const R& x)
        : mpq_class(x) { this->mpq_class::canonicalize(); }
      /*! \brief Conversion assignment operator from another numerical type. */
      template<class R> Rational& operator=(const R& x) {
        (*this)=Rational(x); return *this; }
      //@}
      
      //@{
      //! \name Data access
      /*! \brief The numerator. */
      Integer numerator() const { return this->get_num(); }
      /*! \brief The denominator. */
      Integer denominator() const { return this->get_den();}
      //@}
     public:
      mpq_class get_base() const { return *this; }

#ifdef DOXYGEN
      //@{
      //! \name Arithmetic operations
      /*! \brief The minimum of q1 and q2. */
      friend Rational min(const Rational& q1, const Rational& q2);
      /*! \brief The maximum of q1 and q2. */
      friend Rational max(const Rational& q1, const Rational& q2);
      /*! \brief The absolute value \a q. */
      friend Rational abs(const Rational& q);
      
      /*! \brief In-place addition. */
      friend Rational& operator+=(Rational& q1, const Rational& q2);
      /*! \brief In-place subtraction of a number. */
      friend Rational& operator-=(Rational& q1, const Rational& q2);
      /*! \brief In-place multiplication. */
      friend Rational& operator*=(Rational& q1, const Rational& q2);
      /*! \brief In-place division. */
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
      /*! \brief %Integer power. */
      friend Rational pow(const Rational& q, const Integer& n);
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


    std::istream& 
    operator>>(std::istream& is, Rational& q);
      
    inline Integer numerator(const Rational& num){ 
      return num.get_num(); }
  
    inline Integer denominator(const Rational& num){ 
      return num.get_den();}
  
  
    template<> inline Rational min(const Rational& x1, const Rational& x2) {
      return (x1<=x2) ? x1 : x2; }
    template<> inline Rational max(const Rational& x1, const Rational& x2) {
      return (x1>=x2) ? x1 : x2; }
    template<> inline Rational abs(const Rational& x) {
      return (x>=0) ? x : static_cast<Rational>(-x); }

    template<> inline Rational neg(const Rational& x) {
      return -x; }
    template<> inline Rational add(const Rational& x1, const Rational& x2) {
      return x1+x2; }
    template<> inline Rational sub(const Rational& x1, const Rational& x2) {
      return x1-x2; }
    template<> inline Rational mul(const Rational& x1, const Rational& x2) {
      return x1*x2; }
    template<> inline Rational div(const Rational& x1, const Rational& x2) {
      return x1/x2; }

    template<> inline Rational pow(const Rational& q, const int& n) {
      if(n<0) { return pow(q,-n); }
      Rational r=1; Rational p=q; uint e=1; uint un=n;
      while(e<=un) { if(e&un) { r*=p; } p*=p; e*=2; }
      return r; 
    }      
    
    template<> inline Rational pow(const Rational& q, const uint& n) {
      return pow(q,int(n));
    }      
    
    template<> inline Rational floor(const Rational& x) { 
      return Rational((x.get_num()+x.get_den()-1)/x.get_den()); }
    template<> inline Rational ceil(const Rational& x) { 
      return Rational(x.get_num()/x.get_den()); }

    template<> inline Integer int_down(const Rational& x) { 
      return Integer((x.get_num()+x.get_den()-1)/x.get_den()); }
    template<> inline Integer int_up(const Rational& x) { 
      return Integer(x.get_num()/x.get_den()); }
    template<> inline int int_down(const Rational& x) { 
      return int_down<Integer>(x).get_si(); }
    template<> inline int int_up(const Rational& x) { 
      return int_up<Integer>(x).get_si(); }
    
      
    template<> inline std::string name<Numeric::Rational>() { return "Rational"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Rational> >() { return "Interval<Rational>"; }
    
    template<> inline double conv_approx(const Rational& x) { return x.get_d(); }
 
    template<> inline Rational conv_exact(const int& n) { return Rational(n); }
    template<> inline Rational conv_approx(const int& n) { return conv_exact<Rational>(n); }
    template<> inline Rational conv_down(const int& n) { return conv_exact<Rational>(n); }
    template<> inline Rational conv_up(const int& n) { return conv_exact<Rational>(n); }
 
    template<> inline Rational conv_exact(const double& x) { return Rational(x); }
    template<> inline Rational conv_approx(const double& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_down(const double& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_up(const double& x) { return conv_exact<Rational>(x); }
 
    template<> inline Rational conv_exact(const Rational& x) { return x; }
    template<> inline Rational conv_approx(const Rational& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_down(const Rational& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_up(const Rational& x) { return conv_exact<Rational>(x); }
 
    template<> inline Rational min_exact(const Rational& x1, const Rational& x2) {
      return (x1<=x2) ? x1 : x2; }
    template<> inline Rational min_approx(const Rational& x1, const Rational& x2) { 
      return min_exact(x1,x2); }
    template<> inline Rational min_down(const Rational& x1, const Rational& x2) { 
      return min_exact(x1,x2); }
    template<> inline Rational min_up(const Rational& x1, const Rational& x2) { 
      return min_exact(x1,x2); }
  
    template<> inline Rational max_exact(const Rational& x1, const Rational& x2) {
      return (x1>=x2) ? x1 : x2; }
    template<> inline Rational max_approx(const Rational& x1, const Rational& x2) { 
      return max_exact(x1,x2); }
    template<> inline Rational max_down(const Rational& x1, const Rational& x2) { 
      return max_exact(x1,x2); }
    template<> inline Rational max_up(const Rational& x1, const Rational& x2) { 
      return max_exact(x1,x2); }
  
    
    template<> inline Rational neg_exact(const Rational& x) { return neg(x); }
    template<> inline Rational neg_approx(const Rational& x) { return neg_exact(x); }
    template<> inline Rational neg_down(const Rational& x) { return neg_exact(x); }
    template<> inline Rational neg_up(const Rational& x) { return neg_exact(x); }
    
    template<> inline Rational abs_exact(const Rational& x) { return abs(x); }
    template<> inline Rational abs_approx(const Rational& x) { return abs_exact(x); }
    template<> inline Rational abs_down(const Rational& x) { return abs_exact(x); }
    template<> inline Rational abs_up(const Rational& x) { return abs_exact(x); }
    
    template<> inline Rational add_exact(const Rational& x1, const Rational& x2) { return x1+x2; }
    template<> inline Rational add_down(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
    template<> inline Rational add_up(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
    template<> inline Rational add_approx(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }

    template<> inline Rational sub_exact(const Rational& x1, const Rational& x2) { return x1-x2; }
    template<> inline Rational sub_down(const Rational& x1, const Rational& x2) { return sub_exact(x1,x2); }
    template<> inline Rational sub_up(const Rational& x1, const Rational& x2) { return sub_exact(x1,x2); }
    template<> inline Rational sub_approx(const Rational& x1, const Rational& x2) { return sub_exact(x1,x2); }
    
    template<> inline Rational mul_exact(const Rational& x1, const Rational& x2) { return x1*x2; }
    template<> inline Rational mul_down(const Rational& x1, const Rational& x2) { return mul_exact(x1,x2); }
    template<> inline Rational mul_up(const Rational& x1, const Rational& x2) { return mul_exact(x1,x2); }
    template<> inline Rational mul_approx(const Rational& x1, const Rational& x2) { return mul_exact(x1,x2); }
    
    template<> inline Rational div_exact(const Rational& x1, const Rational& x2) { return x1/x2; }
    template<> inline Rational div_down(const Rational& x1, const Rational& x2) { return div_exact(x1,x2); }
    template<> inline Rational div_up(const Rational& x1, const Rational& x2) { return div_exact(x1,x2); }
    template<> inline Rational div_approx(const Rational& x1, const Rational& x2) { return div_exact(x1,x2); }
  
  
  }

}

#endif /* ARIADNE_RATIONAL_H */
