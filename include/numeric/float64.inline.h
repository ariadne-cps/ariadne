/****************************************************************************
 *            numeric/float64.inline.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file numeric/float64.inline.h
 *  \brief Definition of 64-bit floating point.
 */

#include <iosfwd>
#include <cmath>
#include "numeric/expression.h"
#include "numeric/rounding.h"
#include "numeric/double.h"
#include "numeric/float64.class.h"


#include "numeric/traits.h"
#include "numeric/macros.h"

#include "numeric/integer.h"
#include "numeric/rational.h"

namespace Ariadne {
namespace Numeric {



inline Float64::~Float() { }
inline Float64::Float() : _value() { }
inline Float64::Float(const int& n) : _value(n) { }
inline Float64::Float(const unsigned int& n) : _value(n) { }
inline Float64::Float(const double& x) : _value(x) { }
inline Float64::Float(const Float64& x) : _value(x._value) { }

inline Float64& Float64::operator=(const int& n) { 
  _value=n; return *this; }
inline Float64& Float64::operator=(const unsigned int& n) { 
  _value=n; return *this; }
inline Float64& Float64::operator=(const double& x) { 
  _value=x; return *this; }
inline Float64& Float64::operator=(const Float64& x) { 
  _value=x._value; return *this; }

template<class E> 
inline Float64::Float(const Expression<E>& e) { 
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


// FIXME: Rounding mode!
inline Float64::Float(const Integer& z)
  : _value(mpz_get_d(z._value)) { ARIADNE_ASSERT(*this==z); }
inline Float64& Float64::operator=(const Integer& z) { 
  _value=mpz_get_d(z._value); ARIADNE_ASSERT(*this==z); return *this; }

template<class Rnd> inline Float64::Float(const Rational& q, Rnd rnd) { 
  set_(*this,q,rnd); }
template<class Rnd> inline void Float64::set(const Rational& q, Rnd rnd) {
  set_(*this,q,rnd); }

inline double Float64::get_d() const {
  return _value; }

inline unsigned int precision(const Float64& x) {
  return 64; }
inline void set_precision(Float64& x, unsigned int p) { 
  ARIADNE_ASSERT(false); }

  
inline void nan_(Float64& r) { r._value=std::numeric_limits<double>::quiet_NaN(); }
inline void inf_(Float64& r) { r._value=std::numeric_limits<double>::infinity(); }
inline void eps_(Float64& r) { r._value=std::numeric_limits<double>::epsilon(); }
inline void min_(Float64& r) { r._value=std::numeric_limits<double>::min(); }

inline void next_(Float64& r, const Float64& x, RoundUp rnd) { 
  succ_(r._value,x._value); }
inline void next_(Float64& r, const Float64& x, RoundDown rnd) { 
  prec_(r._value,x._value); }

inline void get_(double& r, const Float64& x, RoundApprox rnd) { 
  r=x._value; }
inline void set_(Rational& r, const Float64& x) { 
  mpq_set_d(r._value,x._value); }

inline void set_(Float64& r, const int& n) { r._value=n; }
inline void set_(Float64& r, const unsigned int& n) { r._value=n; }
inline void set_(Float64& r, const double& x) { r._value=x; }
inline void set_(Float64& r, const Integer& z) { r._value=mpz_get_d(z._value); ARIADNE_ASSERT(r==z); }
inline void set_(Float64& r, const Float64& x) { r._value=x._value; }


template<class Rnd> inline void set_(Float64& r, const int& n, Rnd) {
  r._value=n; }
template<class Rnd> inline void set_(Float64& r, const unsigned int& n, Rnd) {
  r._value=n; }
template<class Rnd> inline void set_(Float64& r, const double& x, Rnd) {
  r._value=x; }
template<class Rnd> inline void set_(Float64& r, const Integer& z, Rnd) { 
  r._value=mpz_get_d(z._value); assert(r==z); }
template<class Rnd> inline void set_(Float64& r, const Float64& z, Rnd) { 
  r._value=z._value; }
template<class Rnd> inline void set_(Float64& r, const Rational& q, Rnd rnd) { 
  r._value=mpq_get_d(q._value); next_(r,r,rnd); }


inline void floor_(int& r, const Float64& x) { r=(int)std::floor(x._value); }
inline void floor_(long int& r, const Float64& x) { r=(long int)std::floor(x._value); }
inline void floor_(Integer& r, const Float64& x);
inline void floor_(Float64& r, const Float64& x) { r._value=std::floor(x._value); }

inline void ceil_(int& r, const Float64& x) { r=(int)std::ceil(x._value); }
inline void ceil_(long int& r, const Float64& x) { r=(long int)std::ceil(x._value); }
inline void ceil_(Integer& r, const Float64& x);
inline void ceil_(Float64& r, const Float64& x) { r._value=std::ceil(x._value); }

	  

// Operations which may be performed exactly
inline void min_(Float64& r, const Float64& x, const Float64& y) {
  r._value = (x._value<=y._value ? x._value : y._value); }
inline void max_(Float64& r, const Float64& x, const Float64& y) {
  r._value = (x._value>=y._value ? x._value : y._value); }
inline void pos_(Float64& r, const Float64& x) {
  r._value=x._value; }
inline void neg_(Float64& r, const Float64& x) {
  r._value=-x._value; }
inline void abs_(Float64& r, const Float64& x) { 
  r._value = (x._value>=0 ? x._value : -x._value); }



// Rounded arithmetic operations
template<class Rnd> 
inline void pos_(Float64& r, const Float64& x, Rnd) {
  r._value=x._value; }
template<class Rnd> 
inline void neg_(Float64& r, const Float64& x, Rnd) {
  r._value=-x._value; }


template<class Rnd> 
inline void add_(Float64& r, const Float64& x, const Float64& y, Rnd rnd) {
  add_(r._value,x._value,y._value,rnd); }

template<class Rnd> 
inline void sub_(Float64& r, const Float64& x, const Float64& y, Rnd rnd) {
  sub_(r._value,x._value,y._value,rnd); }

template<class Rnd> 
inline void mul_(Float64& r, const Float64& x, const Float64& y, Rnd rnd) {
  mul_(r._value,x._value,y._value,rnd); }

template<class Rnd> 
inline void div_(Float64& r, const Float64& x, const Float64& y, Rnd rnd) {
  div_(r._value,x._value,y._value,rnd); }

inline void med_(Float64& r, const Float64& x, const Float64& y, RoundApprox rnd) {
  med_(r._value,x._value,y._value,rnd); }
inline void rad_(Float64& r, const Float64& x, const Float64& y, RoundUp rnd) {
  rad_(r._value,x._value,y._value,rnd); }


// Mixed-mode arithmetic
template<class Rnd> inline 
void mul_(Float64& r, const Float64& x, const int& y, Rnd rnd) {
  mul_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const int& x, const Float64& y, Rnd rnd) {
  mul_(r._value,x,y._value,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const Float64& x, const long int& y, Rnd rnd) {
  mul_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const long int& x, const Float64& y, Rnd rnd) {
  mul_(r._value,x,y._value,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const Float64& x, const unsigned int& y, Rnd rnd) {
  mul_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const unsigned int& x, const Float64& y, Rnd rnd) {
  mul_(r._value,x,y._value,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const Float64& x, const unsigned long int& y, Rnd rnd) {
  mul_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const unsigned long int& x, const Float64& y, Rnd rnd) {
  mul_(r._value,x,y._value,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const Float64& x, const double& y, Rnd rnd) {
  mul_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void mul_(Float64& r, const double& x, const Float64& y, Rnd rnd) {
  mul_(r._value,x,y._value,rnd); }


template<class Rnd> inline 
void div_(Float64& r, const Float64& x, const int& y, Rnd rnd) {
  div_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void div_(Float64& r, const Float64& x, const long int& y, Rnd rnd) {
  div_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void div_(Float64& r, const Float64& x, const unsigned int& y, Rnd rnd) {
  div_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void div_(Float64& r, const Float64& x, const unsigned long int& y, Rnd rnd) {
  div_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void div_(Float64& r, const Float64& x, const double& y, Rnd rnd) {
  div_(r._value,x._value,y,rnd); }

template<class Rnd> inline 
void div_(Float64& r, const int& x, const Float64& y, Rnd rnd) {
  div_(r._value,x,y._value,rnd); }



template<class Rnd> 
inline void pow_(Float64& r, const Float64& x, const unsigned int& n, Rnd rnd) {
  pow_(r._value,x._value,n,rnd); }

template<class Rnd> 
void pow_(Float64& r, const Float64& x, const int& n, Rnd rnd) {
  pow_(r._value,x._value,n,rnd); }



template<class Rnd> 
inline void sqrt_(Float64& r, const Float64& x, Rnd rnd) { 
  sqrt_(r._value,x._value,rnd); }

template<class Rnd> 
inline void hypot_(Float64& r, const Float64& x, const Float64& y, Rnd rnd) { 
  hypot_(r._value,x._value,rnd); }

template<class Rnd> 
inline void exp_(Float64& r, const Float64& x, Rnd rnd) { 
  exp_(r._value,x._value,rnd); }

template<class Rnd> 
inline void log_(Float64& r, const Float64& x, Rnd rnd) { 
  log_(r._value,x._value,rnd); }

  

template<class Rnd> void pi_(Float64& r, Rnd rnd) {
  pi_(r._value,rnd); }

/*
template<> 
inline void pi_(Float64& r, RoundApprox) {
  r._value=3.1415926535897932; }
template<> 
inline void pi_(Float64& r, RoundDown) {
  r._value=3.1415926535897931; }
template<> 
inline void pi_(Float64& r, RoundUp) {
  r._value=3.1415926535897936; }
*/

template<class Rnd> 
inline void sin_(Float64& r, const Float64& x, Rnd rnd) {
  sin_(r._value,x._value,rnd); }
template<class Rnd> 
inline void cos_(Float64& r, const Float64& x, Rnd rnd) {
  cos_(r._value,x._value,rnd); }
template<class Rnd> 
inline void tan_(Float64& r, const Float64& x, Rnd rnd) {
  tan_(r._value,x._value,rnd); }

template<class Rnd> 
inline void asin_(Float64& r, const Float64& x, Rnd rnd) {
  asin_(r._value,x._value,rnd); }
template<class Rnd> 
inline void acos_(Float64& r, const Float64& x, Rnd rnd) {
  acos_(r._value,x._value,rnd); }
template<class Rnd> 
inline void atan_(Float64& r, const Float64& x, Rnd rnd) {
  atan_(r._value,x._value,rnd); }

template<class Rnd> 
inline void sinh_(Float64& r, const Float64& x, Rnd rnd) {
  sinh_(r._value,x._value,rnd); }
template<class Rnd> 
inline void cosh_(Float64& r, const Float64& x, Rnd rnd) {
  cosh_(r._value,x._value,rnd); }
template<class Rnd> 
inline void tanh_(Float64& r, const Float64& x, Rnd rnd) {
  tanh_(r._value,x._value,rnd); }

template<class Rnd> 
inline void asinh_(Float64& r, const Float64& x, Rnd rnd);
template<class Rnd> 
inline void acosh_(Float64& r, const Float64& x, Rnd rnd);
template<class Rnd> 
inline void atanh_(Float64& r, const Float64& x, Rnd rnd);


// Comparison functions
template<class X1, class X2> 
inline int builtin_cmp(const X1& x1, const X2& x2) {
  return ( x1==x2 ? 0 : (x1<x2) ? -1 : +1 ); }

inline int cmp(const Float64& x1, const Float64& x2) {
  return builtin_cmp(x1._value,x2._value); }

inline int cmp(const Float64& x1, const int& x2) {
  return builtin_cmp(x1._value,x2); }
inline int cmp(const Float64& x1, const long int& x2) {
  return builtin_cmp(x1._value,x2); }
inline int cmp(const Float64& x1, const unsigned int& x2) {
  return builtin_cmp(x1._value,x2); }
inline int cmp(const Float64& x1, const unsigned long int& x2) {
  return builtin_cmp(x1._value,x2); }
inline int cmp(const Float64& x1, const double& x2) {
  return builtin_cmp(x1._value,x2); }
inline int cmp(const Float64& x1, const Integer& x2) {
  return mpz_cmp_d(x2._value,x1._value); }
inline int cmp(const Float64& x1, const Rational& x2) {
  return cmp(Rational(x1),x2); }

inline int cmp(const int& x1, const Float64& x2) {
  return builtin_cmp(x1,x2._value); }
inline int cmp(const long int& x1, const Float64& x2) {
  return builtin_cmp(x1,x2._value); }
inline int cmp(const unsigned int& x1, const Float64& x2) {
  return builtin_cmp(x1,x2._value); }
inline int cmp(const unsigned long int& x1, const Float64& x2) {
  return builtin_cmp(x1,x2._value); }
inline int cmp(const double& x1, const Float64& x2) {
  return builtin_cmp(x1,x2._value); }
inline int cmp(const Integer& x1, const Float64& x2) {
  return -mpz_cmp_d(x1._value,x2._value); }
inline int cmp(const Rational& x1, const Float64& x2) {
  return cmp(x1,Rational(x2)); }



// Comparison operators
ARIADNE_VALUE_VALUE_COMPARISON(bool,Float64,Float64);

ARIADNE_VALUE_BUILTIN_COMPARISON(bool,Float64,int);
ARIADNE_VALUE_BUILTIN_COMPARISON(bool,Float64,long int);
ARIADNE_VALUE_BUILTIN_COMPARISON(bool,Float64,unsigned int);
ARIADNE_VALUE_BUILTIN_COMPARISON(bool,Float64,unsigned long int);
ARIADNE_VALUE_BUILTIN_COMPARISON(bool,Float64,double);
ARIADNE_FUNCTION_COMPARISON(bool,Float64,Integer);
ARIADNE_FUNCTION_COMPARISON(bool,Float64,Rational);

ARIADNE_BUILTIN_VALUE_COMPARISON(bool,int,Float64);
ARIADNE_BUILTIN_VALUE_COMPARISON(bool,long int,Float64);
ARIADNE_BUILTIN_VALUE_COMPARISON(bool,unsigned int,Float64);
ARIADNE_BUILTIN_VALUE_COMPARISON(bool,unsigned long int,Float64);
ARIADNE_BUILTIN_VALUE_COMPARISON(bool,double,Float64);
ARIADNE_FUNCTION_COMPARISON(bool,Integer,Float64);
ARIADNE_FUNCTION_COMPARISON(bool,Rational,Float64);

inline std::ostream& operator<<(std::ostream& os, const Float64& x) { return os << x._value; }
inline std::istream& operator>>(std::istream& is, Float64& x) { return is >> x._value; }
inline Float64::Float(const std::string& str) { std::stringstream ss(str); ss>>*this; }
template<class R> std::string name();
template<> inline std::string name<Float64>() { return "Float64"; }
template<> inline std::string name< Interval<Float64> >() { return "Interval64"; }
    


inline Float64 pos(const Float64& x) { return x; }
inline Float64 neg(const Float64& x) { return -x._value; }

inline Float64 operator+(const Float64& x) { return pos(x); }
inline Float64 operator-(const Float64& x) { return neg(x); }

} // namespace Numeric
} // namespace Ariadne
