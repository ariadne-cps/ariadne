/***************************************************************************
 *            numeric/float64.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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

inline std::ostream& operator<<(std::ostream& os, const Float64& x) { 
  return os << x._value; }

inline std::istream& operator>>(std::istream& is, Float64& x) {
  return is >> x._value; }


template<class Rnd> void set_rounding_mode();
template<> inline void set_rounding_mode<RoundUp>() { }
template<> inline void set_rounding_mode<RoundDown>() { }
template<> inline void set_rounding_mode<RoundApprox>() { }

inline Float64::~Float() { }
inline Float64::Float() : _value() { }
inline Float64::Float(const int& n) : _value(n) { }
inline Float64::Float(const unsigned int& n) : _value(n) { }
inline Float64::Float(const double& x) : _value(x) { }
inline Float64::Float(const Integer& n) { _value=mpz_get_d(n._value); assert(*this==n); }
inline Float64::Float(const Float<double>& x) : _value(x._value) { }

inline Float64& Float64::operator=(const int& n) { 
  this->_value=n; return *this; }
inline Float64& Float64::operator=(const unsigned int& n) { 
  this->_value=n; return *this; }
inline Float64& Float64::operator=(const double& x) {
  this->_value=x; return *this; }
inline Float64& Float64::operator=(const Integer& n) { 
  _value=mpz_get_d(n._value); assert(*this==n); return *this; }
inline Float64& Float64::operator=(const Float64& x) { 
  this->_value=x._value; return *this; }

template<class E> 
inline Float64::Float(const Expression<E>& e) : _value() { 
  e.assign_to(*this); }
template<class E> 
inline Float64& Float64::operator=(const Expression<E>& e) {
  e.assign_to(*this); return *this; }

template<class X, class Rnd> 
inline Float64::Float(const X& x, Rnd rnd) { 
  set_(*this,x,rnd); }
template<class E, class Rnd> 
inline Float64::Float(const Expression<E>& e, Rnd rnd) {
  e.assign_to(*this); }

inline double Float64::get_d() const {
  return this->_value; }

template<> 
inline std::string name<Float64>() { 
  return "Float64"; }
template<> 
inline std::string name< Interval<Float64> >() { 
  return "IntervalMP"; }





inline void nan_(Float64& r) { r=std::numeric_limits<double>::quiet_NaN(); }
inline void inf_(Float64& r) { r=std::numeric_limits<double>::infinity(); }
inline void eps_(Float64& r) { r=std::numeric_limits<double>::min(); }

template<class Rnd> 
inline void get_(double& r, const Float64& x, Rnd) { 
  set_rounding_mode<Rnd>(); r=x._value; }

inline void set_(Rational& q, const Float64& x) { 
  mpq_set_d(q._value,x._value); mpq_canonicalize(q._value); }

template<> inline Rational& Rational::operator=(const Float64& x) { 
  set_(*this,x); return *this; }

inline void set_(Float64& r, const int& n) { r._value=n; }
inline void set_(Float64& r, const unsigned int& n) { r._value=n; }
inline void set_(Float64& r, const double& x) { r._value=x; }
inline void set_(Float64& r, const Integer& z) { r._value=mpz_get_d(z._value); assert(r==z); }

// Conversion from any integer type should be able to be performed exactly
inline void set_(Float64& r, const int& n, RoundDown) { r._value=n; }
inline void set_(Float64& r, const int& n, RoundUp) { r._value=n; }
inline void set_(Float64& r, const int& n, RoundApprox) { r._value=n; }

inline void set_(Float64& r, const long int& n, RoundDown) { r._value=n; }
inline void set_(Float64& r, const long int& n, RoundUp) { r._value=n; }
inline void set_(Float64& r, const long int& n, RoundApprox) { r._value=n; }

inline void set_(Float64& r, const unsigned int& n, RoundDown) { r._value=n; }
inline void set_(Float64& r, const unsigned int& n, RoundUp) { r._value=n; }
inline void set_(Float64& r, const unsigned int& n, RoundApprox) { r._value=n; }

inline void set_(Float64& r, const unsigned long int& n, RoundDown) { r._value=n; }
inline void set_(Float64& r, const unsigned long int& n, RoundUp) { r._value=n; }
inline void set_(Float64& r, const unsigned long int& n, RoundApprox) { r._value=n; }

inline void set_(Float64& r, const double& x, RoundDown) { r._value=x; }
inline void set_(Float64& r, const double& n, RoundUp) { r._value=n; }
inline void set_(Float64& r, const double& n, RoundApprox) { r._value=n; }

inline void set_(Float64& r, const Integer& z, RoundDown) { 
  r._value=mpz_get_d(z._value); }
inline void set_(Float64& r, const Integer& z, RoundUp) { 
  r._value=mpz_get_d(z._value); }
inline void set_(Float64& r, const Integer& z, RoundApprox) { 
  r._value=mpz_get_d(z._value); }




// Operations which may be performed exactly
inline void min_(Float64& r, const Float64& x1, const Float64& x2) { 
  r = (x1<=x2) ? x1 : x2; }
inline void max_(Float64& r, const Float64& x1, const Float64& x2) { 
  r = (x1>=x2) ? x1 : x2; }
inline void pos_(Float64& r, const Float64& x) { 
  r._value = x._value; }
inline void neg_(Float64& r, const Float64& x) { 
  r._value = -x._value; }
inline void abs_(Float64& r, const Float64& x) { 
  if(x>=0) { r._value=x._value; } else { r._value=-x._value; } }

	  



inline void floor_(Float64& r, const Float64& x) {
  r=std::floor(x._value); }
inline void floor_(Integer& r, const Float64& x) {
  r=(long int)(std::floor(x._value)); }
inline void floor_(long int& r, const Float64& x) {
  r=(long int)(std::floor(x._value)); }

inline void ceil_(Float64& r, const Float64& x) {
  r=std::ceil(x._value); }
inline void ceil_(Integer& r, const Float64& x) {
  r=(long int)(std::ceil(x._value)); }
inline void ceil_(long int& r, const Float64& x) {
  r=(long int)(std::ceil(x._value)); }
	  
inline long int floor(const Float64& x) {
  long int r; floor_(r,x); return r; }
inline long int ceil(const Float64& x) {
  long int r; ceil_(r,x); return r; }




template<class Rnd> 
inline void add_(Float64& r, const Float64& x1,const Float64& x2, Rnd) {
  set_rounding_mode<Rnd>(); r._value=x1._value+x2._value; }

inline void add_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
  r._value=Float64::rounding().add_down(x1._value,x2._value); }
inline void add_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
  r._value=Float64::rounding().add_up(x1._value,x2._value); }
inline void add_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
  r._value=x1._value+x2._value; }

    
inline void sub_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
  r._value=Float64::rounding().sub_down(x1._value,x2._value); }
inline void sub_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
  r._value=Float64::rounding().sub_up(x1._value,x2._value); }
inline void sub_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
  r._value=x1._value-x2._value; }
    
inline void mul_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
  r._value=Float64::rounding().mul_down(x1._value,x2._value); }
inline void mul_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
  r._value=Float64::rounding().mul_up(x1._value,x2._value); }
inline void mul_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
  r._value=x1._value*x2._value; }
     
inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
  r._value=Float64::rounding().div_down(x1._value,x2._value); }
inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
  r._value=Float64::rounding().div_up(x1._value,x2._value); }
inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
  r._value=x1._value/x2._value; }

inline void med_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
  r._value=(x1._value+x2._value)/2; }
inline void rad_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
  r._value=Float64::rounding().div_up(
             Float64::rounding().sub_up(x2._value,x1._value),2); }
    

	  
	  
// Mixed-mode arithmetic
inline void mul_(Float64& r, const Float64& x, const int& n, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const Float64& x,const int& n, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const Float64& x, const int& n, RoundApprox) {
  r._value=x._value*n; }
     
inline void mul_(Float64& r, const int& n, const Float64& x, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const int& n, const Float64& x, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const int& n, const Float64& x, RoundApprox) {
  r._value=x._value*n; }

inline void mul_(Float64& r, const Float64& x, const long int& n, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const Float64& x,const long int& n, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const Float64& x, const long int& n, RoundApprox) {
  r._value=x._value*n; }
     
inline void mul_(Float64& r, const long int& n, const Float64& x, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const long int& n, const Float64& x, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const long int& n, const Float64& x, RoundApprox) {
  r._value=x._value*n; }

inline void mul_(Float64& r, const Float64& x, const unsigned int& n, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const Float64& x,const unsigned int& n, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const Float64& x, const unsigned int& n, RoundApprox) {
  r._value=x._value*n; }
     
inline void mul_(Float64& r, const unsigned int& n, const Float64& x, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const unsigned int& n, const Float64& x, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const unsigned int& n, const Float64& x, RoundApprox) {
  r._value=x._value*n; }

inline void mul_(Float64& r, const Float64& x, const unsigned long int& n, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const Float64& x,const unsigned long int& n, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const Float64& x, const unsigned long int& n, RoundApprox) {
  r._value=x._value*n; }
     
inline void mul_(Float64& r, const unsigned long int& n, const Float64& x, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,n); }
inline void mul_(Float64& r, const unsigned long int& n, const Float64& x, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,n); }
inline void mul_(Float64& r, const unsigned long int& n, const Float64& x, RoundApprox) {
  r._value=x._value*n; }

inline void mul_(Float64& r, const Float64& x, const double& d, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,d); }
inline void mul_(Float64& r, const Float64& x,const double& d, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,d); }
inline void mul_(Float64& r, const Float64& x, const double& d, RoundApprox) {
  r._value=x._value*d; }

inline void mul_(Float64& r, const double& d, const Float64& x, RoundDown) {
  r._value=Float64::rounding().mul_down(x._value,d); }
inline void mul_(Float64& r, const double& d, const Float64& x, RoundUp) {
  r._value=Float64::rounding().mul_up(x._value,d); }
inline void mul_(Float64& r, const double& d, const Float64& x, RoundApprox) {
  r._value=x._value*d; }


inline void div_(Float64& r, const Float64& x,const int& n, RoundDown) {
  r._value=Float64::rounding().div_down(x._value,n); }
inline void div_(Float64& r, const Float64& x,const int& n, RoundUp) {
  r._value=Float64::rounding().div_up(x._value,n); }
inline void div_(Float64& r, const Float64& x,const int& n, RoundApprox) {
  r._value=x._value/n; }

inline void div_(Float64& r, const Float64& x,const long int& n, RoundDown) {
  r._value=Float64::rounding().div_down(x._value,n); }
inline void div_(Float64& r, const Float64& x,const long int& n, RoundUp) {
  r._value=Float64::rounding().div_up(x._value,n); }
inline void div_(Float64& r, const Float64& x,const long int& n, RoundApprox) {
  r._value=x._value/n; }

inline void div_(Float64& r, const Float64& x,const double& y, RoundDown) {
  r._value=Float64::rounding().div_down(x._value,y); }
inline void div_(Float64& r, const Float64& x,const double& y, RoundUp) {
  r._value=Float64::rounding().div_up(x._value,y); }
inline void div_(Float64& r, const Float64& x,const double& y, RoundApprox) {
  r._value=x._value/y; }


template<class Rnd>
inline void pow_(Float64& r, const Float64& x,const uint& n, Rnd rnd) {
  Float64 p=x; uint m=n; r=1; 
  while(m) { if(m%2) { mul_(r,r,p,rnd); } mul_(p,p,p,rnd); m/=2; }
}

inline void next_(Float64& r, const Float64& x, RoundDown) { 
  Float64 min=std::numeric_limits<double>::min();
  sub_(r,x,min,round_down); }
inline void next_(Float64& r, const Float64& x, RoundUp) { 
  Float64 min=std::numeric_limits<double>::min();
  add_(r,x,min,round_up); }

template<class Rnd, class T> 
inline Float<T> next(const Float<T>& x) {
  Float<T> r; next_(r,x,Rnd()); return r; }

inline Float64 next(Float64& x, RoundDown) { 
  Float64 r; next_(r,x,round_down); return r; }
inline Float64 next(Float64& x, RoundUp) { 
  Float64 r; next_(r,x,round_down); return r; }


// Convert from a rational. These need the next_ function
inline void set_(Float64& r, const Rational& q, RoundDown) { 
  r._value=mpq_get_d(q._value); next_(r,r,round_down); }
inline void set_(Float64& r, const Rational& q, RoundUp) { 
  r._value=mpq_get_d(q._value); next_(r,r,round_up); }
inline void set_(Float64& r, const Rational& q, RoundApprox) { 
  r._value=mpq_get_d(q._value); }




inline void sqrt_(Float64& r, const Float64& x, RoundDown) { 
  r._value=Float64::rounding().sqrt_down(x._value); }
inline void sqrt_(Float64& r, const Float64& x, RoundUp) { 
  r._value=Float64::rounding().sqrt_up(x._value); }
inline void sqrt_(Float64& r, const Float64& x, RoundApprox) { 
  r._value=std::sqrt(x._value); }


inline void exp_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::exp(x._value); }
inline void exp_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().exp_down(x._value); }
inline void exp_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().exp_up(x._value); }

inline void log_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::log(x._value); }
inline void log_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().log_down(x._value); }
inline void log_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().log_up(x._value); }


inline void pi_(Float64& r, RoundApprox) {
  r._value=3.1415926535897931; }
inline void pi_(Float64& r, RoundDown) {
  r._value=3.1415926535897928; }
inline void pi_(Float64& r, RoundUp) {
  r._value=3.1415926535897934; }

inline void sin_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::sin(x._value); }
inline void sin_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().sin_down(x._value); }
inline void sin_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().sin_up(x._value); }

inline void cos_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::cos(x._value); }
inline void cos_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().cos_down(x._value); }
inline void cos_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().cos_up(x._value); }

inline void tan_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::tan(x._value); }
inline void tan_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().tan_down(x._value); }
inline void tan_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().tan_up(x._value); }

inline void asin_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::asin(x._value); }
inline void asin_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().asin_down(x._value); }
inline void asin_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().asin_up(x._value); }

inline void acos_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::acos(x._value); }
inline void acos_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().acos_down(x._value); }
inline void acos_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().acos_up(x._value); }

inline void atan_(Float64& r, const Float64& x, RoundApprox) {
  r._value=std::atan(x._value); }
inline void atan_(Float64& r, const Float64& x, RoundDown) {
  r._value=Float64::rounding().atan_down(x._value); }
inline void atan_(Float64& r, const Float64& x, RoundUp) {
  r._value=Float64::rounding().atan_up(x._value); }





// Comparison operators
inline bool operator==(const Float64& x1, const Float64& x2) {
  return x1._value==x2._value; }
inline bool operator!=(const Float64& x1, const Float64& x2) {
  return x1._value!=x2._value; }
inline bool operator<=(const Float64& x1, const Float64& x2) {
  return x1._value<=x2._value; }
inline bool operator>=(const Float64& x1, const Float64& x2) {
  return x1._value>=x2._value; }
inline bool operator< (const Float64& x1, const Float64& x2) {
  return x1._value< x2._value; }
inline bool operator> (const Float64& x1, const Float64& x2) {
  return x1._value> x2._value; }

inline bool operator==(const Float64& x1, const int& x2) {
  return x1._value==x2; }
inline bool operator!=(const Float64& x1, const int& x2) {
  return x1._value!=x2; }
inline bool operator<=(const Float64& x1, const int& x2) {
  return x1._value<=x2; }
inline bool operator>=(const Float64& x1, const int& x2) {
  return x1._value>=x2; }
inline bool operator< (const Float64& x1, const int& x2) {
  return x1._value< x2; }
inline bool operator> (const Float64& x1, const int& x2) {
  return x1._value> x2; }    

inline bool operator==(const int& x1, const Float64& x2) {
  return x1==x2._value; }
inline bool operator!=(const int& x1, const Float64& x2) {
  return x1!=x2._value; }
inline bool operator<=(const int& x1, const Float64& x2) {
  return x1<=x2._value; }
inline bool operator>=(const int& x1, const Float64& x2) {
  return x1>=x2._value; }
inline bool operator< (const int& x1, const Float64& x2) {
  return x1< x2._value; }
inline bool operator> (const int& x1, const Float64& x2) {
  return x1> x2._value; }

inline bool operator==(const Float64& x1, const long int& x2) {
  return x1._value==x2; }
inline bool operator!=(const Float64& x1, const long int& x2) {
  return x1._value!=x2; }
inline bool operator<=(const Float64& x1, const long int& x2) {
  return x1._value<=x2; }
inline bool operator>=(const Float64& x1, const long int& x2) {
  return x1._value>=x2; }
inline bool operator< (const Float64& x1, const long int& x2) {
  return x1._value< x2; }
inline bool operator> (const Float64& x1, const long int& x2) {
  return x1._value> x2; }    

inline bool operator==(const long int& x1, const Float64& x2) {
  return x1==x2._value; }
inline bool operator!=(const long int& x1, const Float64& x2) {
  return x1!=x2._value; }
inline bool operator<=(const long int& x1, const Float64& x2) {
  return x1<=x2._value; }
inline bool operator>=(const long int& x1, const Float64& x2) {
  return x1>=x2._value; }
inline bool operator< (const long int& x1, const Float64& x2) {
  return x1< x2._value; }
inline bool operator> (const long int& x1, const Float64& x2) {
  return x1> x2._value; }

inline bool operator==(const Float64& x1, const unsigned int& x2) {
  return x1._value==x2; }
inline bool operator!=(const Float64& x1, const unsigned int& x2) {
  return x1._value!=x2; }
inline bool operator<=(const Float64& x1, const unsigned int& x2) {
  return x1._value<=x2; }
inline bool operator>=(const Float64& x1, const unsigned int& x2) {
  return x1._value>=x2; }
inline bool operator< (const Float64& x1, const unsigned int& x2) {
  return x1._value< x2; }
inline bool operator> (const Float64& x1, const unsigned int& x2) {
  return x1._value> x2; }    

inline bool operator==(const unsigned int& x1, const Float64& x2) {
  return x1==x2._value; }
inline bool operator!=(const unsigned int& x1, const Float64& x2) {
  return x1!=x2._value; }
inline bool operator<=(const unsigned int& x1, const Float64& x2) {
  return x1<=x2._value; }
inline bool operator>=(const unsigned int& x1, const Float64& x2) {
  return x1>=x2._value; }
inline bool operator< (const unsigned int& x1, const Float64& x2) {
  return x1< x2._value; }
inline bool operator> (const unsigned int& x1, const Float64& x2) {
  return x1> x2._value; }

inline bool operator==(const Float64& x1, const unsigned long int& x2) {
  return x1._value==x2; }
inline bool operator!=(const Float64& x1, const unsigned long int& x2) {
  return x1._value!=x2; }
inline bool operator<=(const Float64& x1, const unsigned long int& x2) {
  return x1._value<=x2; }
inline bool operator>=(const Float64& x1, const unsigned long int& x2) {
  return x1._value>=x2; }
inline bool operator< (const Float64& x1, const unsigned long int& x2) {
  return x1._value< x2; }
inline bool operator> (const Float64& x1, const unsigned long int& x2) {
  return x1._value> x2; }    

inline bool operator==(const unsigned long int& x1, const Float64& x2) {
  return x1==x2._value; }
inline bool operator!=(const unsigned long int& x1, const Float64& x2) {
  return x1!=x2._value; }
inline bool operator<=(const unsigned long int& x1, const Float64& x2) {
  return x1<=x2._value; }
inline bool operator>=(const unsigned long int& x1, const Float64& x2) {
  return x1>=x2._value; }
inline bool operator< (const unsigned long int& x1, const Float64& x2) {
  return x1< x2._value; }
inline bool operator> (const unsigned long int& x1, const Float64& x2) {
  return x1> x2._value; }

inline bool operator==(const Float64& x1, const double& x2) {
  return x1._value==x2; }
inline bool operator!=(const Float64& x1, const double& x2) {
  return x1._value!=x2; }
inline bool operator<=(const Float64& x1, const double& x2) {
  return x1._value<=x2; }
inline bool operator>=(const Float64& x1, const double& x2) {
  return x1._value>=x2; }
inline bool operator< (const Float64& x1, const double& x2) {
  return x1._value< x2; }
inline bool operator> (const Float64& x1, const double& x2) {
  return x1._value> x2; }    

inline bool operator==(const double& x1, const Float64& x2) {
  return x1==x2._value; }
inline bool operator!=(const double& x1, const Float64& x2) {
  return x1!=x2._value; }
inline bool operator<=(const double& x1, const Float64& x2) {
  return x1<=x2._value; }
inline bool operator>=(const double& x1, const Float64& x2) {
  return x1>=x2._value; }
inline bool operator< (const double& x1, const Float64& x2) {
  return x1< x2._value; }
inline bool operator> (const double& x1, const Float64& x2) {
  return x1> x2._value; }    

inline bool operator==(const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value)==0; }
inline bool operator!=(const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value)!=0; }
inline bool operator<=(const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value)>=0; }
inline bool operator>=(const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value)<=0; }
inline bool operator< (const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value)> 0; }
inline bool operator> (const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value)< 0; }

inline bool operator==(const Integer& x1, const Float64& x2) {
  return mpz_cmp_d(x1._value,x2._value)==0; }
inline bool operator!=(const Integer& x1, const Float64& x2) {
  return mpz_cmp_d(x1._value,x2._value)!=0; }
inline bool operator<=(const Integer& x1, const Float64& x2) {
  return mpz_cmp_d(x1._value,x2._value)<=0; }
inline bool operator>=(const Integer& x1, const Float64& x2) {
  return mpz_cmp_d(x1._value,x2._value)>=0; }
inline bool operator< (const Integer& x1, const Float64& x2) {
  return mpz_cmp_d(x1._value,x2._value)< 0; }
inline bool operator> (const Integer& x1, const Float64& x2) {
  return mpz_cmp_d(x1._value,x2._value)> 0; }

inline bool operator==(const Float64& x1, const Rational& x2) {
  return Rational(x1)==x2; }
inline bool operator!=(const Float64& x1, const Rational& x2) {
  return Rational(x1)!=x2; }
inline bool operator<=(const Float64& x1, const Rational& x2) {
  return Rational(x1)<=x2; }
inline bool operator>=(const Float64& x1, const Rational& x2) {
  return Rational(x1)>=x2; }
inline bool operator< (const Float64& x1, const Rational& x2) {
  return Rational(x1)< x2; }
inline bool operator> (const Float64& x1, const Rational& x2) {
  return Rational(x1)> x2; }    

inline bool operator==(const Rational& x1, const Float64& x2) {
  return x1==Rational(x2); }
inline bool operator!=(const Rational& x1, const Float64& x2) {
  return x1!=Rational(x2); }
inline bool operator<=(const Rational& x1, const Float64& x2) {
  return x1<=Rational(x2); }
inline bool operator>=(const Rational& x1, const Float64& x2) {
  return x1>=Rational(x2); }
inline bool operator< (const Rational& x1, const Float64& x2) {
  return x1< Rational(x2); }
inline bool operator> (const Rational& x1, const Float64& x2) {
  return x1> Rational(x2); }    


// Comparison functions
inline int cmp(const Float64& x1, const Float64& x2) {
  return (x1==x2) ? 0 : (x1>x2) ? +1 : -1; }
inline int cmp(const Float64& x1, const int& x2) {
  return (x1==x2) ? 0 : (x1>x2) ? +1 : -1; }
inline int cmp(const Float64& x1, const long int& x2) {
  return (x1==x2) ? 0 : (x1>x2) ? +1 : -1; }
inline int cmp(const Float64& x1, const unsigned int& x2) {
  return (x1==x2) ? 0 : (x1>x2) ? +1 : -1; }
inline int cmp(const Float64& x1, const unsigned long int& x2) {
  return (x1==x2) ? 0 : (x1>x2) ? +1 : -1; }
inline int cmp(const Float64& x1, const double& x2) {
  return (x1==x2) ? 0 : (x1>x2) ? +1 : -1; }
inline int cmp(const Float64& x1, const Integer& x2) {
  return -mpz_cmp_d(x2._value,x1._value); }
inline int cmp(const Float64& x1, const Rational& x2) {
  return cmp(Rational(x1),x2); }




} // namespace Numeric
} // namespace Ariadne


