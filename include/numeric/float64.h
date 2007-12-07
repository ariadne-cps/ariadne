/***************************************************************************
 *            numeric/float64.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file numeric/float64.h
 *  \brief Type definitions and conversion operators for 64-bit fixed precision floating point numbers.
 */

#ifndef ARIADNE_NUMERIC_FLOAT64_H
#define ARIADNE_NUMERIC_FLOAT64_H

#include <cmath>
#include <limits>
#include <mpfr.h>
#include <boost/numeric/interval/rounded_arith.hpp>
#include <boost/numeric/interval/rounded_transc.hpp>
#include <boost/numeric/interval/hw_rounding.hpp>

#include "numeric/traits.h"
#include "numeric/rounding.h"

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/interval.class.h"


namespace Ariadne {
  namespace Numeric {

    template<class T> class Float;
    typedef Float<double> Float64;
    typedef Interval<Float64> Interval64;

#ifdef DOXYGEN
    /*!\ingroup Numeric
     * \brief A 64-bit fixed-precision floating point number.
     *  
     * Standard operations are not exact, but must support interval arithmetic.
     *
     * Currently implemented by the built-in type double.
     */
    template<> class Float<double> { };
#else
    //typedef double Float64;
#endif
    template<>
    class Float<double> 
      : public Value< Float<double> >
    {
     public:
      typedef boost::numeric::interval_lib::rounded_arith_std<double> rounded_arith;
      typedef boost::numeric::interval_lib::rounded_transc_std<double> rounding;
     public:
      Float();
      Float(const int& n);
      Float(const unsigned int& n);
      Float(const double& x);
      Float(const Float64& x);
      
      explicit Float(const std::string& x);
      
      Float64& operator=(const int& n);
      Float64& operator=(const double& x);
      Float64& operator=(const Float64& x);
      
      template<class E> Float(const Expression<E>& e);
      template<class E> Float64& operator=(const Expression<E>& e);

      template<class X, class Rnd> Float(const X& x, Rnd rnd);
      template<class E, class Rnd> Float(const Expression<E>& x, Rnd rnd);
     public:
      double _value;
    };

    inline std::ostream& operator<<(std::ostream& os, const Rational& x) { 
      return os << x._value; }
    
    inline std::ostream& operator<<(std::ostream& os, const Float64& x) { 
      return os << x._value; }
    
    inline std::istream& operator>>(std::istream& is, Float64& x) {
      return is >> x._value; }

    template<class Rnd> void set_rounding_mode();
    template<> inline void set_rounding_mode<RoundUp>() { }
    template<> inline void set_rounding_mode<RoundDown>() { }
    template<> inline void set_rounding_mode<RoundApprox>() { }

    // Declare input/output operators
    std::ostream& operator<<(std::ostream& os, const Float64& x);
    std::istream& operator>>(std::istream& is, Float64& x);

    inline Float64::Float() : _value() { }
    inline Float64::Float(const int& n) : _value(n) { }
    inline Float64::Float(const uint& n) : _value(n) { }
    inline Float64::Float(const double& x) : _value(x) { }
    inline Float64::Float(const Float<double>& x) : _value(x._value) { }
      
    inline Float64& Float64::operator=(const int& n) { this->_value=n; return *this; }
    inline Float64& Float64::operator=(const double& x) { this->_value=x; return *this; }
    inline Float64& Float64::operator=(const Float64& x) { this->_value=x._value; return *this; }
      
    template<class E> inline Float64::Float(const Expression<E>& e) : _value() { e.assign_to(*this); }
    template<class E> inline Float64& Float64::operator=(const Expression<E>& e) { e.assign_to(*this); return *this; }

    template<class X, class Rnd> inline Float64::Float(const X& x, Rnd rnd) { set_(*this,x,rnd); }
    template<class E, class Rnd> inline Float64::Float(const Expression<E>& e, Rnd rnd) { e.assign_to(*this); }

    template<class Rnd, class T> inline Float<T> next(const Float<T>& x) {
      Float<T> r; next_(r,x,Rnd()); return r; }

    inline void nan_(Float64& r) { r=std::numeric_limits<double>::quiet_NaN(); }
    inline void inf_(Float64& r) { r=std::numeric_limits<double>::infinity(); }

    template<class Rnd> inline void set_(double& r, const Float64& x, Rnd) { 
      set_rounding_mode<Rnd>(); r=x._value; }

    inline void set_(Float64& r, const int& n, RoundDown) { r._value=n; }
    inline void set_(Float64& r, const int& n, RoundUp) { r._value=n; }
    inline void set_(Float64& r, const int& n, RoundApprox) { r._value=n; }

    inline void set_(Float64& r, const unsigned int& n, RoundDown) { r._value=n; }
    inline void set_(Float64& r, const unsigned int& n, RoundUp) { r._value=n; }
    inline void set_(Float64& r, const unsigned int& n, RoundApprox) { r._value=n; }

    inline void set_(Float64& r, const double& x, RoundDown) { r._value=x; }
    inline void set_(Float64& r, const double& n, RoundUp) { r._value=n; }
    inline void set_(Float64& r, const double& n, RoundApprox) { r._value=n; }

    inline void set_(Float64& r, const Integer& z, RoundDown) { r._value=mpz_get_d(z._value); }
    inline void set_(Float64& r, const Integer& z, RoundUp) { r._value=mpz_get_d(z._value); }
    inline void set_(Float64& r, const Integer& z, RoundApprox) { r._value=mpz_get_d(z._value); }

    template<> inline std::string name<Numeric::Float64>() { return "Float64"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Float64> >() { return "Interval<Float64>"; }
   
    inline bool operator==(const Float64& x1, const Float64& x2) { return x1._value==x2._value; }
    inline bool operator!=(const Float64& x1, const Float64& x2) { return x1._value!=x2._value; }
    inline bool operator<=(const Float64& x1, const Float64& x2) { return x1._value<=x2._value; }
    inline bool operator>=(const Float64& x1, const Float64& x2) { return x1._value>=x2._value; }
    inline bool operator< (const Float64& x1, const Float64& x2) { return x1._value< x2._value; }
    inline bool operator> (const Float64& x1, const Float64& x2) { return x1._value> x2._value; }
    
    inline bool operator==(const Float64& x1, const int& x2) { return x1._value==x2; }
    inline bool operator!=(const Float64& x1, const int& x2) { return x1._value!=x2; }
    inline bool operator<=(const Float64& x1, const int& x2) { return x1._value<=x2; }
    inline bool operator>=(const Float64& x1, const int& x2) { return x1._value>=x2; }
    inline bool operator< (const Float64& x1, const int& x2) { return x1._value< x2; }
    inline bool operator> (const Float64& x1, const int& x2) { return x1._value> x2; }    

    inline bool operator==(const int& x1, const Float64& x2) { return x1==x2._value; }
    inline bool operator!=(const int& x1, const Float64& x2) { return x1!=x2._value; }
    inline bool operator<=(const int& x1, const Float64& x2) { return x1<=x2._value; }
    inline bool operator>=(const int& x1, const Float64& x2) { return x1>=x2._value; }
    inline bool operator< (const int& x1, const Float64& x2) { return x1< x2._value; }
    inline bool operator> (const int& x1, const Float64& x2) { return x1> x2._value; }
    
    inline bool operator==(const Float64& x1, const uint& x2) { return x1._value==x2; }
    inline bool operator!=(const Float64& x1, const uint& x2) { return x1._value!=x2; }
    inline bool operator<=(const Float64& x1, const uint& x2) { return x1._value<=x2; }
    inline bool operator>=(const Float64& x1, const uint& x2) { return x1._value>=x2; }
    inline bool operator< (const Float64& x1, const uint& x2) { return x1._value< x2; }
    inline bool operator> (const Float64& x1, const uint& x2) { return x1._value> x2; }    

    inline bool operator==(const uint& x1, const Float64& x2) { return x1==x2._value; }
    inline bool operator!=(const uint& x1, const Float64& x2) { return x1!=x2._value; }
    inline bool operator<=(const uint& x1, const Float64& x2) { return x1<=x2._value; }
    inline bool operator>=(const uint& x1, const Float64& x2) { return x1>=x2._value; }
    inline bool operator< (const uint& x1, const Float64& x2) { return x1< x2._value; }
    inline bool operator> (const uint& x1, const Float64& x2) { return x1> x2._value; }
    
    inline bool operator==(const Float64& x1, const double& x2) { return x1._value==x2; }
    inline bool operator!=(const Float64& x1, const double& x2) { return x1._value!=x2; }
    inline bool operator<=(const Float64& x1, const double& x2) { return x1._value<=x2; }
    inline bool operator>=(const Float64& x1, const double& x2) { return x1._value>=x2; }
    inline bool operator< (const Float64& x1, const double& x2) { return x1._value< x2; }
    inline bool operator> (const Float64& x1, const double& x2) { return x1._value> x2; }    

    inline bool operator==(const double& x1, const Float64& x2) { return x1==x2._value; }
    inline bool operator!=(const double& x1, const Float64& x2) { return x1!=x2._value; }
    inline bool operator<=(const double& x1, const Float64& x2) { return x1<=x2._value; }
    inline bool operator>=(const double& x1, const Float64& x2) { return x1>=x2._value; }
    inline bool operator< (const double& x1, const Float64& x2) { return x1< x2._value; }
    inline bool operator> (const double& x1, const Float64& x2) { return x1> x2._value; }    

    inline bool operator==(const Float64& x1, const Integer& x2) { return mpz_cmp_d(x2._value,x1._value)==0; }
    inline bool operator!=(const Float64& x1, const Integer& x2) { return mpz_cmp_d(x2._value,x1._value)!=0; }
    inline bool operator<=(const Float64& x1, const Integer& x2) { return mpz_cmp_d(x2._value,x1._value)>=0; }
    inline bool operator>=(const Float64& x1, const Integer& x2) { return mpz_cmp_d(x2._value,x1._value)<=0; }
    inline bool operator< (const Float64& x1, const Integer& x2) { return mpz_cmp_d(x2._value,x1._value)> 0; }
    inline bool operator> (const Float64& x1, const Integer& x2) { return mpz_cmp_d(x2._value,x1._value)< 0; }
      
    inline bool operator==(const Integer& x1, const Float64& x2) { return mpz_cmp_d(x1._value,x2._value)==0; }
    inline bool operator!=(const Integer& x1, const Float64& x2) { return mpz_cmp_d(x1._value,x2._value)!=0; }
    inline bool operator<=(const Integer& x1, const Float64& x2) { return mpz_cmp_d(x1._value,x2._value)<=0; }
    inline bool operator>=(const Integer& x1, const Float64& x2) { return mpz_cmp_d(x1._value,x2._value)>=0; }
    inline bool operator< (const Integer& x1, const Float64& x2) { return mpz_cmp_d(x1._value,x2._value)< 0; }
    inline bool operator> (const Integer& x1, const Float64& x2) { return mpz_cmp_d(x1._value,x2._value)> 0; }
      
    inline bool operator==(const Float64& x1, const Rational& x2) { return Rational(x1)==x2; }
    inline bool operator!=(const Float64& x1, const Rational& x2) { return Rational(x1)!=x2; }
    inline bool operator<=(const Float64& x1, const Rational& x2) { return Rational(x1)<=x2; }
    inline bool operator>=(const Float64& x1, const Rational& x2) { return Rational(x1)>=x2; }
    inline bool operator< (const Float64& x1, const Rational& x2) { return Rational(x1)< x2; }
    inline bool operator> (const Float64& x1, const Rational& x2) { return Rational(x1)> x2; }    
      
    inline bool operator==(const Rational& x1, const Float64& x2) { return x1==Rational(x2); }
    inline bool operator!=(const Rational& x1, const Float64& x2) { return x1!=Rational(x2); }
    inline bool operator<=(const Rational& x1, const Float64& x2) { return x1<=Rational(x2); }
    inline bool operator>=(const Rational& x1, const Float64& x2) { return x1>=Rational(x2); }
    inline bool operator< (const Rational& x1, const Float64& x2) { return x1< Rational(x2); }
    inline bool operator> (const Rational& x1, const Float64& x2) { return x1> Rational(x2); }    
      
      
     
    
      
 
    // Operations which may be performed exactly
    inline void min_(Float64& r, const Float64& x1, const Float64& x2) { r = (x1<=x2) ? x1 : x2; }
    inline void max_(Float64& r, const Float64& x1, const Float64& x2) { r = (x1>=x2) ? x1 : x2; }
    inline void pos_(Float64& r, const Float64& x) { r._value = x._value; }
    inline void neg_(Float64& r, const Float64& x) { r._value = -x._value; }
    inline void abs_(Float64& r, const Float64& x) { if(x>=0) { r._value=x._value; } else { r._value=-x._value; } }

    inline Float64 abs(const Float64& x) { if(x>=0) { return x; } else { return -x; } }



  }
}

#include "arithmetic64.h"
#include "transcendental64.h"

#endif /* ARIADNE_NUMERIC_FLOAT64_H */
