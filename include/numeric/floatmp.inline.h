/***************************************************************************
 *            numeric/floatmp.inline.h
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
 
#include <mpfr.h>

#include "numeric/traits.h"
#include "numeric/macros.h"
#include "numeric/rounding.h"

#include "numeric/integer.h"
#include "numeric/rational.h"

namespace Ariadne {




inline FloatMP::~Float() { 
  mpfr_clear(this->_value); }
inline FloatMP::Float() { 
  mpfr_init_set_si(this->_value,0,GMP_RNDN); }
inline FloatMP::Float(const int& n) {
  mpfr_init_set_si(this->_value,n,GMP_RNDN); }
inline FloatMP::Float(const unsigned int& n) { 
  mpfr_init_set_ui(this->_value,n,GMP_RNDN); }
inline FloatMP::Float(const double& x) { 
  mpfr_init_set_d(this->_value,x,GMP_RNDN); }
inline FloatMP::Float(const Integer& n) { 
  mpfr_init_set_z(this->_value,n._value,GMP_RNDN); 
  assert(*this==n); }
inline FloatMP::Float(const FloatMP& x) { 
  mpfr_init_set(this->_value,x._value,GMP_RNDN); }

inline FloatMP& FloatMP::operator=(const int& n) { 
  mpfr_set_si(this->_value,n,GMP_RNDN); return *this; }
inline FloatMP& FloatMP::operator=(const unsigned int& n) { 
  mpfr_set_ui(this->_value,n,GMP_RNDN); return *this; }
inline FloatMP& FloatMP::operator=(const double& x) { 
  mpfr_set_d(this->_value,x,GMP_RNDN); return *this; }
inline FloatMP& FloatMP::operator=(const Integer& n) { 
  mpfr_set_z(this->_value,n._value,GMP_RNDN); 
  assert(*this==n); return *this; }
inline FloatMP& FloatMP::operator=(const FloatMP& x) { 
  if(this!=&x) { mpfr_set(this->_value,x._value,GMP_RNDN); } return *this; }

template<class E> 
inline FloatMP::Float(const Expression<E>& e) { 
  mpfr_init(_value); e.assign_to(*this); }
template<class E> 
inline FloatMP& FloatMP::operator=(const Expression<E>& e) { 
  e.assign_to(*this); return *this; }

template<class X, class Rnd> 
inline FloatMP::Float(const X& x, Rnd rnd) { 
  mpfr_init(_value); set_(*this,x,rnd); }
template<class E, class Rnd> 
inline FloatMP::Float(const Expression<E>& e, Rnd rnd) {
  mpfr_init(_value); e.assign_to(*this); }


inline unsigned int FloatMP::default_precision() { 
  return mpfr_get_default_prec(); }
inline void FloatMP::set_default_precision(unsigned int p) { 
  mpfr_set_default_prec(p); }
inline unsigned int FloatMP::precision() const { 
  return mpfr_get_prec(this->_value); }
inline void FloatMP::set_precision(unsigned int p) { 
  mpfr_set_prec(this->_value,p); }

template<class Rnd> inline FloatMP::Float(const Rational& q, Rnd rnd) {
  mpfr_init_set_q(this->_value,q._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void FloatMP::set(const Rational& q, Rnd rnd) {
  mpfr_set_q(this->_value,q._value,mpfr_rounding_mode(rnd)); }

inline double FloatMP::get_d() const {
  return mpfr_get_d(this->_value, GMP_RNDN); }

inline unsigned int precision(const FloatMP& x) {
  return mpfr_get_prec(x._value); }
inline void set_precision(FloatMP& x, unsigned int p) { 
  mpfr_set_prec(x._value,p); }

  
template<> 
inline std::string name<FloatMP>() { 
  return "FloatMP"; }
template<> 
inline std::string name< Interval<FloatMP> >() { 
  return "IntervalMP"; }
    




inline void nan_(FloatMP& r) { mpfr_set_nan(r._value); }
inline void inf_(FloatMP& r) { mpfr_set_inf(r._value,+1); }
inline void eps_(FloatMP& r) { mpfr_set_si(r._value,0,GMP_RNDN); mpfr_nextabove(r._value); }

inline void next_(FloatMP& r, const FloatMP& x, RoundUp) { 
  mpfr_set(r._value,x._value,GMP_RNDU); mpfr_nextabove(r._value); }
inline void next_(FloatMP& r, const FloatMP& x, RoundDown) { 
  mpfr_set(r._value,x._value,GMP_RNDD); mpfr_nextbelow(r._value); }



inline void get_(double& r, const FloatMP& x, RoundApprox) { 
  r=mpfr_get_d(x._value,GMP_RNDN); }
inline void set_(Rational& r, const FloatMP& x) { 
  mpf_t f; mpf_init2(f,mpfr_get_prec(x._value)+128); mpfr_get_f(f,x._value,GMP_RNDN); 
  mpq_set_f(r._value,f); mpq_canonicalize(r._value);
  mpf_clear(f);
  if(mpfr_cmp_q(x._value,r._value)!=0) {
    std::cerr << "ERROR: conversion from FloatMP x to Rational q failed; x="<<x<<", q="<<r<<std::endl;
    ARIADNE_ASSERT(mpfr_cmp_q(x._value,r._value)==0); 
  }
}

inline void set_(FloatMP& r, const int& n) { mpfr_set_si(r._value,n,GMP_RNDN); }
inline void set_(FloatMP& r, const unsigned int& n) { mpfr_set_ui(r._value,n,GMP_RNDN); }
inline void set_(FloatMP& r, const double& x) { mpfr_set_d(r._value,x,GMP_RNDN); }
inline void set_(FloatMP& r, const Integer& n) { mpfr_set_z(r._value,n._value,GMP_RNDN); assert(r==n); }


template<class Rnd> inline void set_(FloatMP& r, const int& x, Rnd rnd) {
  mpfr_set_si(r._value,x,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void set_(FloatMP& r, const long int& x, Rnd rnd) {
  mpfr_set_si(r._value,x,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void set_(FloatMP& r, const unsigned int& n, Rnd rnd) {
  mpfr_set_ui(r._value,n,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void set_(FloatMP& r, const unsigned long int& n, Rnd rnd) {
  mpfr_set_ui(r._value,n,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void set_(FloatMP& r, const double& x, Rnd rnd) {
  mpfr_set_d(r._value,x,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void set_(FloatMP& r, const Integer& x, Rnd rnd) {
  mpfr_set_z(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> inline void set_(FloatMP& r, const Rational& x, Rnd rnd) {
  mpfr_set_q(r._value,x._value,mpfr_rounding_mode(rnd)); }




inline void floor_(int& r, const FloatMP& x) { 
  r=mpfr_get_si(x._value,GMP_RNDD); }
inline void floor_(long int& r, const FloatMP& x) { 
  r=mpfr_get_si(x._value,GMP_RNDD); }
inline void floor_(Integer& r, const FloatMP& x) { 
  mpfr_get_z(r._value,x._value,GMP_RNDD); }
inline void floor_(FloatMP& r, const FloatMP& x) { 
  mpfr_floor(r._value,x._value); }

inline void ceil_(int& r, const FloatMP& x) { 
  r=mpfr_get_si(x._value,GMP_RNDU); }
inline void ceil_(long int& r, const FloatMP& x) { 
  r=mpfr_get_si(x._value,GMP_RNDU); }
inline void ceil_(Integer& r, const FloatMP& x) { 
  mpfr_get_z(r._value,x._value,GMP_RNDU); }
inline void ceil_(FloatMP& r, const FloatMP& x) { 
  mpfr_ceil(r._value,x._value); }


	  

// Operations which may be performed exactly
inline void min_(FloatMP& r, const FloatMP& x, const FloatMP& y) {
  assert(r.precision()>=std::max(x.precision(),y.precision())); 
  r = (x<=y ? x : y); }
inline void max_(FloatMP& r, const FloatMP& x, const FloatMP& y) {
  assert(r.precision()>=std::max(x.precision(),y.precision())); 
  r = (x>=y ? x : y); }
inline void pos_(FloatMP& r, const FloatMP& x) {
  assert(r.precision()>=x.precision()); 
  mpfr_set(r._value,x._value,GMP_RNDN); }
inline void neg_(FloatMP& r, const FloatMP& x) {
  assert(r.precision()>=x.precision()); 
  mpfr_neg(r._value,x._value,GMP_RNDN); }
inline void abs_(FloatMP& r, const FloatMP& x) { 
  assert(r.precision()>=x.precision()); 
  mpfr_abs(r._value,x._value,GMP_RNDN); }




// Rounded operations which may be performed exactly
template<class Rnd>
inline void min_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) { 
  mpfr_min(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }

template<class Rnd>
inline void max_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) { 
  mpfr_max(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }

template<class Rnd>
inline void abs_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) { 
  mpfr_abs(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd>
inline void pos_(FloatMP& r, const FloatMP& x, Rnd rnd) { 
  mpfr_set(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd>
inline void neg_(FloatMP& r, const FloatMP& x, Rnd rnd) { 
  mpfr_neg(r._value,x._value,mpfr_rounding_mode(rnd)); }



// Rounded arithmetic operations
template<class Rnd> 
inline void add_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) {
  mpfr_add(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }


template<class Rnd> 
inline void sub_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) {
  mpfr_sub(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void mul_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) {
  mpfr_mul(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void div_(FloatMP& r, const FloatMP& x, const FloatMP& y, Rnd rnd) {
  mpfr_div(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }


inline void med_(FloatMP& r, const FloatMP& x, const FloatMP& y, RoundApprox rnd) {
  mpfr_add(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); 
  mpfr_div_ui(r._value,r._value,2u,mpfr_rounding_mode(rnd)); }
inline void rad_(FloatMP& r, const FloatMP& x, const FloatMP& y, RoundUp rnd) {
  mpfr_sub(r._value,y._value,x._value,mpfr_rounding_mode(rnd)); 
  mpfr_div_ui(r._value,r._value,2u,mpfr_rounding_mode(rnd)); }




// Mixed-mode arithmetic
template<class Rnd> inline 
void mul_(FloatMP& r, const FloatMP& x, const int& y, Rnd rnd) {
  mpfr_mul_si(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const int& x, const FloatMP& y, Rnd rnd) {
  mpfr_mul_si(r._value,y._value,x,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const FloatMP& x, const long int& y, Rnd rnd) {
  mpfr_mul_si(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const long int& x, const FloatMP& y, Rnd rnd) {
  mpfr_mul_si(r._value,y._value,x,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const FloatMP& x, const unsigned int& y, Rnd rnd) {
  mpfr_mul_ui(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const unsigned int& x, const FloatMP& y, Rnd rnd) {
  mpfr_mul_ui(r._value,y._value,x,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const FloatMP& x, const unsigned long int& y, Rnd rnd) {
  mpfr_mul_ui(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const unsigned long int& x, const FloatMP& y, Rnd rnd) {
  mpfr_mul_ui(r._value,y._value,x,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const FloatMP& x, const double& y, Rnd rnd) {
  mpfr_t t; mpfr_init_set_d(t,y,mpfr_rounding_mode(rnd)); 
  mpfr_mul(r._value,x._value,t,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void mul_(FloatMP& r, const double& x, const FloatMP& y, Rnd rnd) {
  mpfr_t t; mpfr_init_set_d(t,x,mpfr_rounding_mode(rnd)); 
  mpfr_mul(r._value,y._value,t,mpfr_rounding_mode(rnd)); }


template<class Rnd> inline 
void div_(FloatMP& r, const FloatMP& x, const int& y, Rnd rnd) {
  mpfr_div_si(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void div_(FloatMP& r, const FloatMP& x, const long int& y, Rnd rnd) {
  mpfr_div_si(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void div_(FloatMP& r, const FloatMP& x, const unsigned int& y, Rnd rnd) {
  mpfr_div_ui(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void div_(FloatMP& r, const FloatMP& x, const unsigned long int& y, Rnd rnd) {
  mpfr_div_ui(r._value,x._value,y,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void div_(FloatMP& r, const FloatMP& x, const double& y, Rnd rnd) {
  mpfr_t t; mpfr_init_set_d(t,y,mpfr_rounding_mode(rnd)); 
  mpfr_div(r._value,x._value,t,mpfr_rounding_mode(rnd)); }

template<class Rnd> inline 
void div_(FloatMP& r, const int& x, const FloatMP& y, Rnd rnd) {
  mpfr_si_div(r._value,x,y._value,mpfr_rounding_mode(rnd)); }


template<class Rnd> 
inline void pow_(FloatMP& r, const FloatMP& x, const unsigned int& n, Rnd rnd) {
  mpfr_pow_ui(r._value,x._value,n,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void pow_(FloatMP& r, const FloatMP& x, const int& n, Rnd rnd) {
  mpfr_pow_si(r._value,x._value,n,mpfr_rounding_mode(rnd)); }



template<class Rnd> 
inline void sqrt_(FloatMP& r, const FloatMP& x, Rnd rnd) { 
  mpfr_sqrt(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void hypot_(FloatMP& r, const FloatMP& x,const FloatMP& y, Rnd rnd) { 
  mpfr_hypot(r._value,x._value,y._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void exp_(FloatMP& r, const FloatMP& x, Rnd rnd) { 
  mpfr_exp(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void log_(FloatMP& r, const FloatMP& x, Rnd rnd) { 
  mpfr_log(r._value,x._value,mpfr_rounding_mode(rnd)); }




template<class Rnd> 
inline void pi_(FloatMP& r, Rnd rnd) {
  mpfr_const_pi(r._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void sin_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_sin(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void cos_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_cos(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void tan_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_tan(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void asin_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_asin(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void acos_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_acos(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void atan_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_atan(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void sinh_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_sinh(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void cosh_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_cosh(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void tanh_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_tanh(r._value,x._value,mpfr_rounding_mode(rnd)); }

template<class Rnd> 
inline void asinh_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_asinh(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void acosh_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_acosh(r._value,x._value,mpfr_rounding_mode(rnd)); }
template<class Rnd> 
inline void atanh_(FloatMP& r, const FloatMP& x, Rnd rnd) {
  mpfr_atanh(r._value,x._value,mpfr_rounding_mode(rnd)); }




// Comparison functions
inline int cmp(const FloatMP& x1, const FloatMP& x2) {
  return mpfr_cmp(x1._value,x2._value); }
inline int cmp(const FloatMP& x1, const int& x2) {
  return mpfr_cmp_si(x1._value,x2); }
inline int cmp(const FloatMP& x1, const long int& x2) {
  return mpfr_cmp_si(x1._value,x2); }
inline int cmp(const FloatMP& x1, const unsigned int& x2) {
  return mpfr_cmp_ui(x1._value,x2); }
inline int cmp(const FloatMP& x1, const unsigned long int& x2) {
  return mpfr_cmp_ui(x1._value,x2); }
inline int cmp(const FloatMP& x1, const double& x2) {
  return mpfr_cmp_d(x1._value,x2); }
inline int cmp(const FloatMP& x1, const Integer& x2) {
  return mpfr_cmp_z(x1._value,x2._value); }
inline int cmp(const FloatMP& x1, const Rational& x2) {
  return mpfr_cmp_q(x1._value,x2._value); }



// Comparison operators
inline bool operator==(const FloatMP& x1, const FloatMP& x2) {
  return mpfr_equal_p(x1._value,x2._value); }
inline bool operator!=(const FloatMP& x1, const FloatMP& x2) {
  return mpfr_lessgreater_p(x1._value,x2._value); }
inline bool operator<=(const FloatMP& x1, const FloatMP& x2) {
  return mpfr_lessequal_p(x1._value,x2._value); }
inline bool operator>=(const FloatMP& x1, const FloatMP& x2) {
  return mpfr_greaterequal_p(x1._value,x2._value); }
inline bool operator< (const FloatMP& x1, const FloatMP& x2) {
  return mpfr_less_p(x1._value,x2._value); }
inline bool operator> (const FloatMP& x1, const FloatMP& x2) {
  return mpfr_greater_p(x1._value,x2._value); }

ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,int);
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,long int);
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,unsigned int);
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,unsigned long int);
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,double);
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,Integer);
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,FloatMP,Rational);



} // namespace Ariadne

