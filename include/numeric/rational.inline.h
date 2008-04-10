/***************************************************************************
 *            rational.inline.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
#include <gmp.h>
#include <iostream>
#include <cassert>

#include "numeric/rounding.h"
#include "numeric/operators.h"

namespace Ariadne {
namespace Numeric {
  
// These functions are not provided by GMP
inline void mpq_min(mpq_ptr rop, mpq_srcptr op1, mpq_srcptr op2) { 
  if(mpq_cmp(op1,op2)<=0) { mpq_set(rop,op1); } else { mpq_set(rop,op2); } 
}
inline void mpq_max(mpq_ptr rop, mpq_srcptr op1, mpq_srcptr op2) { 
  if(mpq_cmp(op1,op2)>=0) { mpq_set(rop,op1); } else { mpq_set(rop,op2); } 
}
inline void mpq_pow_ui(mpq_ptr rop, mpq_srcptr op1, unsigned long int op2) { 
  mpz_pow_ui(mpq_numref(rop),mpq_numref(op1),op2); 
  mpz_pow_ui(mpq_denref(rop),mpq_denref(op1),op2); 
}
inline void mpq_pow_si(mpq_ptr rop, mpq_srcptr op1, long int op2) { 
  if(op2==0) { mpq_set_ui(rop,1u,1u); }
  else if(op2>0) { mpz_pow_ui(mpq_numref(rop),mpq_numref(op1),op2); mpz_pow_ui(mpq_denref(rop),mpq_denref(op1),op2); }
  else { mpz_pow_ui(mpq_numref(rop),mpq_denref(op1),-op2); mpz_pow_ui(mpq_denref(rop),mpq_numref(op1),-op2); } 
}

inline void mpq_pow_z(mpq_ptr rop, mpq_srcptr op1, mpz_srcptr op2) {
  long int tmp=mpz_get_si(op2);
  assert(mpz_cmp_si(op2,tmp)==0);
  mpq_pow_si(rop,op1,tmp);
}


inline Rational::~Rational() { 
  mpq_clear(this->_value); }
inline Rational::Rational() { 
  mpq_init(this->_value); 
  mpq_set_si(this->_value,0,1); }
inline Rational::Rational(const int& n, const int& d) { 
  mpq_init(this->_value); 
  mpq_set_si(this->_value,n,d); 
  mpq_canonicalize(this->_value); }
inline Rational::Rational(const Integer& n, const Integer& d) { 
  mpq_init(this->_value); 
  mpz_set(mpq_numref(this->_value),n._value); 
  mpz_set(mpq_denref(this->_value),d._value); }

inline Rational::Rational(const int& n) { 
  mpq_init(this->_value); 
  mpq_set_si(this->_value,n,1); }
inline Rational::Rational(const long int& n) { 
  mpq_init(this->_value); 
  mpq_set_si(this->_value,n,1); }
inline Rational::Rational(const unsigned int& n) { 
  mpq_init(this->_value); 
  mpq_set_ui(this->_value,n,1u); }
inline Rational::Rational(const unsigned long int& n) { 
  mpq_init(this->_value); 
  mpq_set_ui(this->_value,n,1u); }
inline Rational::Rational(const double& x) { 
  mpq_init(this->_value); 
  mpq_set_d(this->_value,x); }
inline Rational::Rational(const mpz_class& z) { 
  mpq_init(this->_value); 
  mpq_set_z(this->_value,z.get_mpz_t()); }
inline Rational::Rational(const mpf_class& x) { 
  mpq_init(this->_value); 
  mpq_set_f(this->_value,x.get_mpf_t()); }
inline Rational::Rational(const mpq_class& q) { 
  mpq_init(this->_value); 
  mpq_set(this->_value,q.get_mpq_t()); }
inline Rational::Rational(const Integer& z) { 
  mpq_init(this->_value); 
  mpq_set_z(this->_value,z._value); }
inline Rational::Rational(const Rational& q) { 
  mpq_init(this->_value); 
  mpq_set(this->_value,q._value); }
template<class T> 
inline Rational::Rational(const Float<T>& x) { 
  mpq_init(this->_value); 
  set_(*this,x); }

inline Rational& Rational::operator=(const int& n) { 
  mpq_set_si(this->_value,n,1); return *this; }
inline Rational& Rational::operator=(const long int& n) { 
  mpq_set_si(this->_value,n,1); return *this; }
inline Rational& Rational::operator=(const unsigned int& n) { 
  mpq_set_ui(this->_value,n,1u); return *this; }
inline Rational& Rational::operator=(const unsigned long int& n) { 
  mpq_set_ui(this->_value,n,1u); return *this; }
inline Rational& Rational::operator=(const double& x) { 
  mpq_set_d(this->_value,x); return *this; }
inline Rational& Rational::operator=(const Integer& z) { 
  mpq_set_z(this->_value,z._value); return *this; }
inline Rational& Rational::operator=(const Rational& q) { 
  mpq_set(this->_value,q._value); return *this; }
template<class T> 
inline Rational& Rational::operator=(const Float<T>& x) { 
  set_(*this,x); return *this; }

template<class E> 
inline Rational::Rational(const Expression<E>& e) { 
  mpq_init(this->_value); 
  e.assign_to(*this); 
  mpq_canonicalize(this->_value); }
template<class E> 
inline Rational& Rational::operator=(const Expression<E>& e) { 
  e.assign_to(*this); 
  return *this; }

      
inline Integer Rational::numerator() const { 
  Integer z; mpz_set(z._value,mpq_numref(this->_value)); return z; }
inline Integer Rational::denominator() const { 
  Integer z; mpz_set(z._value,mpq_denref(this->_value)); return z; }
inline void Rational::canonicalize() { 
  mpq_canonicalize(this->_value); }


inline void set_(double& r, const Rational& x, RoundApprox) { 
  r=mpq_get_d(x._value); }

inline void set_(Rational& r, const int& x) { 
  mpq_set_si(r._value,x,1); }
inline void set_(Rational& r, const unsigned int& x) { 
  mpq_set_ui(r._value,x,1u); }
inline void set_(Rational& r, const double& x) { 
  mpq_set_d(r._value,x); }
inline void set_(Rational& r, const Integer& x) { 
  mpq_set_z(r._value,x._value); }
inline void set_(Rational& r, const Rational& x) { 
  mpq_set(r._value,x._value); }

template<class E>
inline void set_(Rational& r, const Expression<E>& x) { 
  x.assign_to(r); }


inline void min_(Rational& r, const Rational& x1, const Rational& x2) { 
  r=std::min(x1,x2); }
inline void max_(Rational& r, const Rational& x1, const Rational& x2) { 
  r=std::max(x1,x2); }

inline void pos_(Rational& r, const Rational& x) { 
  mpq_set(r._value,x._value); }
inline void neg_(Rational& r, const Rational& x) { 
  mpq_neg(r._value,x._value); }
inline void abs_(Rational& r, const Rational& x) { 
  mpq_abs(r._value,x._value); }
inline Rational abs(const Rational& x) { 
  Rational r; abs_(r,x); return r; }


inline void inf_(Rational& r) { 
  r=Rational(1,0); }

inline void add_(Rational& r, const Rational& x1, const Rational& x2) { 
  mpq_add(r._value,x1._value,x2._value); }
inline void sub_(Rational& r, const Rational& x1, const Rational& x2) { 
  mpq_sub(r._value,x1._value,x2._value); }
inline void mul_(Rational& r, const Rational& x1, const Rational& x2) { 
  mpq_mul(r._value,x1._value,x2._value); }
inline void div_(Rational& r, const Rational& x1, const Rational& x2) { 
  mpq_div(r._value,x1._value,x2._value); }
inline void med_(Rational& r, const Rational& x1, const Rational& x2) { 
  mpq_add(r._value,x1._value,x2._value); 
  mpq_div_2exp(r._value,r._value,1); }
inline void rad_(Rational& r, const Rational& x1, const Rational& x2) { 
  mpq_sub(r._value,x2._value,x1._value); 
  mpq_div_2exp(r._value,r._value,1); }

// Mixed-mode arithmetic involving division by an Integer
inline void div_(Rational& r, const Integer& x1, const Integer& x2) {
  r=Rational(x1,x2); }
inline void div_(Rational& r, const Integer& x1, const int& x2) {
  r=Rational(x1,Integer(x2)); }
inline void div_(Rational& r, const Integer& x1, const unsigned int& x2) {
  r=Rational(x1,Integer(x2)); }

// Mixed-mode arithmetic with int, long int, unsigned int, long unsigned int, double and Integer
inline void add_(Rational& r, const Rational& x1, const int& x2) {
  add_(r,x1,Rational(x2)); }
inline void sub_(Rational& r, const Rational& x1, const int& x2) {
  sub_(r,x1,Rational(x2)); }
inline void mul_(Rational& r, const Rational& x1, const int& x2) {
  mul_(r,x1,Rational(x2)); }
inline void div_(Rational& r, const Rational& x1, const int& x2) {
  div_(r,x1,Rational(x2)); }

inline void add_(Rational& r, const int& x1, const Rational& x2) {
  add_(r,Rational(x1),x2); }
inline void sub_(Rational& r, const int& x1, const Rational& x2) {
  sub_(r,Rational(x1),x2); }
inline void mul_(Rational& r, const int& x1, const Rational& x2) {
  mul_(r,Rational(x1),x2); }
inline void div_(Rational& r, const int& x1, const Rational& x2) {
  div_(r,Rational(x1),x2); }

inline void add_(Rational& r, const Rational& x1, const long int& x2) {
  add_(r,x1,Rational(x2)); }
inline void sub_(Rational& r, const Rational& x1, const long int& x2) {
  sub_(r,x1,Rational(x2)); }
inline void mul_(Rational& r, const Rational& x1, const long int& x2) {
  mul_(r,x1,Rational(x2)); }
inline void div_(Rational& r, const Rational& x1, const long int& x2) {
  div_(r,x1,Rational(x2)); }

inline void add_(Rational& r, const long int& x1, const Rational& x2) {
  add_(r,Rational(x1),x2); }
inline void sub_(Rational& r, const long int& x1, const Rational& x2) {
  sub_(r,Rational(x1),x2); }
inline void mul_(Rational& r, const long int& x1, const Rational& x2) {
  mul_(r,Rational(x1),x2); }
inline void div_(Rational& r, const long int& x1, const Rational& x2) {
  div_(r,Rational(x1),x2); }

inline void add_(Rational& r, const Rational& x1, const unsigned int& x2) {
  add_(r,x1,Rational(x2)); }
inline void sub_(Rational& r, const Rational& x1, const unsigned int& x2) {
  sub_(r,x1,Rational(x2)); }
inline void mul_(Rational& r, const Rational& x1, const unsigned int& x2) {
  mul_(r,x1,Rational(x2)); }
inline void div_(Rational& r, const Rational& x1, const unsigned int& x2) {
  div_(r,x1,Rational(x2)); }

inline void add_(Rational& r, const unsigned int& x1, const Rational& x2) {
  add_(r,Rational(x1),x2); }
inline void sub_(Rational& r, const unsigned int& x1, const Rational& x2) {
  sub_(r,Rational(x1),x2); }
inline void mul_(Rational& r, const unsigned int& x1, const Rational& x2) {
  mul_(r,Rational(x1),x2); }
inline void div_(Rational& r, const unsigned int& x1, const Rational& x2) {
  div_(r,Rational(x1),x2); }

inline void add_(Rational& r, const Rational& x1, const unsigned long int& x2) {
  add_(r,x1,Rational(x2)); }
inline void sub_(Rational& r, const Rational& x1, const unsigned long int& x2) {
  sub_(r,x1,Rational(x2)); }
inline void mul_(Rational& r, const Rational& x1, const unsigned long int& x2) {
  mul_(r,x1,Rational(x2)); }
inline void div_(Rational& r, const Rational& x1, const unsigned long int& x2) {
  div_(r,x1,Rational(x2)); }

inline void add_(Rational& r, const unsigned long int& x1, const Rational& x2) {
  add_(r,Rational(x1),x2); }
inline void sub_(Rational& r, const unsigned long int& x1, const Rational& x2) {
  sub_(r,Rational(x1),x2); }
inline void mul_(Rational& r, const unsigned long int& x1, const Rational& x2) {
  mul_(r,Rational(x1),x2); }
inline void div_(Rational& r, const unsigned long int& x1, const Rational& x2) {
  div_(r,Rational(x1),x2); }

inline void add_(Rational& r, const Rational& x1, const double& x2) {
  add_(r,x1,Rational(x2)); }
inline void sub_(Rational& r, const Rational& x1, const double& x2) {
  sub_(r,x1,Rational(x2)); }
inline void mul_(Rational& r, const Rational& x1, const double& x2) {
  mul_(r,x1,Rational(x2)); }
inline void div_(Rational& r, const Rational& x1, const double& x2) {
  div_(r,x1,Rational(x2)); }

inline void add_(Rational& r, const double& x1, const Rational& x2) {
  add_(r,Rational(x1),x2); }
inline void sub_(Rational& r, const double& x1, const Rational& x2) {
  sub_(r,Rational(x1),x2); }
inline void mul_(Rational& r, const double& x1, const Rational& x2) {
  mul_(r,Rational(x1),x2); }
inline void div_(Rational& r, const double& x1, const Rational& x2) {
  div_(r,Rational(x1),x2); }

inline void add_(Rational& r, const Rational& x1, const Integer& x2) {
  add_(r,x1,Rational(x2)); }
inline void sub_(Rational& r, const Rational& x1, const Integer& x2) {
  sub_(r,x1,Rational(x2)); }
inline void mul_(Rational& r, const Rational& x1, const Integer& x2) {
  mul_(r,x1,Rational(x2)); }
inline void div_(Rational& r, const Rational& x1, const Integer& x2) {
  div_(r,x1,Rational(x2)); }

inline void add_(Rational& r, const Integer& x1, const Rational& x2) {
  add_(r,Rational(x1),x2); }
inline void sub_(Rational& r, const Integer& x1, const Rational& x2) {
  sub_(r,Rational(x1),x2); }
inline void mul_(Rational& r, const Integer& x1, const Rational& x2) {
  mul_(r,Rational(x1),x2); }
inline void div_(Rational& r, const Integer& x1, const Rational& x2) {
  div_(r,Rational(x1),x2); }

inline void pow_(Rational& r, const Rational& x1, const unsigned int& x2) { 
  mpq_pow_ui(r._value,x1._value,x2); }
inline void pow_(Rational& r, const Rational& x1, const int& x2) { 
  mpq_pow_si(r._value,x1._value,x2); }
inline void pow_(Rational& r, const Rational& x1, const Integer& x2) { 
  mpq_pow_z(r._value,x1._value,x2._value); }

inline Rational med(const Rational& x1, const Rational& x2) {
  Rational r; med_(r,x1,x2); return r; }
inline Rational rad(const Rational& x1, const Rational& x2) {
  Rational r; rad_(r,x1,x2); return r; }


inline void floor_(Rational& r, const Rational& x) { 
  mpz_t tmp; mpz_init(tmp); 
  mpz_fdiv_q(tmp,mpq_numref(x._value),mpq_denref(x._value)); 
  mpq_set_z(r._value,tmp); mpz_clear(tmp); }
inline void ceil_(Rational& r, const Rational& x) { 
  mpz_t tmp; mpz_init(tmp); 
  mpz_cdiv_q(tmp,mpq_numref(x._value),mpq_denref(x._value)); 
  mpq_set_z(r._value,tmp); mpz_clear(tmp); }

inline void floor_(Integer& r, const Rational& x) { 
  mpz_cdiv_q(r._value,mpq_numref(x._value),mpq_denref(x._value)); }
inline void ceil_(Integer& r, const Rational& x) { 
  mpz_fdiv_q(r._value,mpq_numref(x._value),mpq_denref(x._value)); }

inline void floor_(long int& r, const Rational& x) { 
  mpz_t tmp; 
  mpz_init(tmp); 
  mpz_fdiv_q(tmp,mpq_numref(x._value),mpq_denref(x._value)); 
  r=mpz_get_si(tmp); 
  mpz_clear(tmp); }
inline void ceil_(long int& r, const Rational& x) { 
  mpz_t tmp; 
  mpz_init(tmp); 
  mpz_cdiv_q(tmp,mpq_numref(x._value),mpq_denref(x._value)); 
  r=mpz_get_si(tmp); 
  mpz_clear(tmp); }

inline long int floor(const Rational& x) {
  long int r; floor_(r,x); return r; }
inline long int ceil(const Rational& x) {
  long int r; ceil_(r,x); return r; }



// Rounded arithmetic needed for interval class.
template<class Rnd, class X>
inline void set_(Rational& r, const X& x, Rnd) {
  set_(r,x); }

template<class Rnd> 
inline void add_(Rational& r, const Rational& x, const Rational& y, Rnd) { 
  add_(r,x,y); }
template<class Rnd> 
inline void sub_(Rational& r, const Rational& x, const Rational& y, Rnd) { 
  sub_(r,x,y); }
template<class Rnd> 
inline void mul_(Rational& r, const Rational& x, const Rational& y, Rnd) { 
  mul_(r,x,y); }
template<class Rnd> 
inline void div_(Rational& r, const Rational& x, const Rational& y, Rnd) { 
  div_(r,x,y); }
template<class Rnd> 
inline void med_(Rational& r, const Rational& x, const Rational& y, Rnd) { 
  med_(r,x,y); }
template<class Rnd> 
inline void rad_(Rational& r, const Rational& x, const Rational& y, Rnd) { 
  rad_(r,x,y); }




// Comparison functions
inline int cmp(const Rational& x, const Rational& y) {
  return mpq_cmp(x._value,y._value); }

inline int cmp(const Rational& x, const int& y) {
  return mpq_cmp_si(x._value,y,1); }
inline int cmp(const Rational& x, const unsigned int& y) {
  return mpq_cmp_ui(x._value,y,1u); }
inline int cmp(const Rational& x, const double& y) {
  mpq_t t; mpq_init(t); mpq_set_d(t,y); 
  int r=mpq_cmp(x._value,t); mpq_clear(t); return r; }
inline int cmp(const Rational& x, const Integer& y) {
  mpq_t t; mpq_init(t); mpq_set_z(t,y._value); 
  int r=mpq_cmp(x._value,t); mpq_clear(t); return r; }

template<class E> 
inline int cmp(const Rational& x, const Expression<E>& y) {
  return cmp(x,Rational(y)); }
template<class E> 
inline int cmp(const Expression<E>& x, const Rational& y) {
  return cmp(Rational(x),y); }


// Comparison operators
inline bool operator==(const Rational& x1, const Rational& x2) {
  return mpq_equal(x1._value,x2._value); }
inline bool operator!=(const Rational& x1, const Rational& x2) {
  return !mpq_equal(x1._value,x2._value); }
inline bool operator<=(const Rational& x1, const Rational& x2) {
  return mpq_cmp(x1._value,x2._value)<=0; }
inline bool operator>=(const Rational& x1, const Rational& x2) {
  return mpq_cmp(x1._value,x2._value)>=0; }
inline bool operator<(const Rational& x1, const Rational& x2) {
  return mpq_cmp(x1._value,x2._value)<0; }
inline bool operator>(const Rational& x1, const Rational& x2) {
  return mpq_cmp(x1._value,x2._value)>0; }


// Mixed comparison operators
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,Rational,int)
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,Rational,uint)
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,Rational,double)
ARIADNE_MIXED_FUNCTION_COMPARISON(bool,Rational,Integer)

template<class E>
inline bool operator==(const Rational& x1, const Expression<E>& x2) {
  return x1==Rational(x2); }
template<class E>
inline bool operator!=(const Rational& x1, const Expression<E>& x2) {
  return x1!=Rational(x2); }
template<class E>
inline bool operator<=(const Rational& x1, const Expression<E>& x2) {
  return x1<=Rational(x2); }
template<class E>
inline bool operator>=(const Rational& x1, const Expression<E>& x2) {
  return x1>=Rational(x2); }
template<class E>
inline bool operator<(const Rational& x1, const Expression<E>& x2) {
  return x1<Rational(x2); }
template<class E>
inline bool operator>(const Rational& x1, const Expression<E>& x2) {
  return x1>Rational(x2); }

template<class E>
inline bool operator==(const Expression<E>& x1, const Rational& x2) {
  return Rational(x1)==x2; }
template<class E>
inline bool operator!=(const Expression<E>& x1, const Rational& x2) {
  return Rational(x1)!=x2; }
template<class E>
inline bool operator<=(const Expression<E>& x1, const Rational& x2) {
  return Rational(x1)<=x2; }
template<class E>
inline bool operator>=(const Expression<E>& x1, const Rational& x2) {
  return Rational(x1)>=x2; }
template<class E>
inline bool operator<(const Expression<E>& x1, const Rational& x2) {
  return Rational(x1)<x2; }
template<class E>
inline bool operator>(const Expression<E>& x1, const Rational& x2) {
  return Rational(x1)>x2; }




// Explicit rounding operators for Rational numbers (provided as a convenience)

template<class E, class Rnd> 
inline Rational::Rational(const Expression<E>& e, Rnd) { 
  mpq_init(this->_value); set_(*this,e); }

template<class X, class Rnd> 
inline Rational::Rational(const X& x, Rnd) { 
  mpq_init(this->_value); set_(*this,x); }


inline Rational& operator+=(Rational& q, const int& x) {
  add_(q,q,x); return q; }

inline Rational& operator+=(Rational& q, const long int& x) {
  add_(q,q,x); return q; }

inline Rational& operator+=(Rational& q, const double& x) {
  add_(q,q,x); return q; }


inline Rational add_approx(const Rational& x, const Rational& y) {
  Rational r; add_(r,x,y); return r; }

inline Rational sub_approx(const Rational& x, const Rational& y) {
  Rational r; sub_(r,x,y); return r; }

inline Rational mul_approx(const Rational& x, const Rational& y) {
  Rational r; mul_(r,x,y); return r; }

inline Rational div_approx(const Rational& x, const Rational& y) {
  Rational r; div_(r,x,y); return r; }

inline Rational div_up(const Rational& x, const long int& y) {
  Rational r; div_(r,x,y); return r; }


inline void quot_(Integer& r, const Rational& x1, const Rational& x2) { 
  Rational q; Integer n,d; 
  mpq_div(q._value,x1._value,x2._value); 
  mpq_get_num(n._value,q._value);
  mpq_get_den(d._value,q._value);
  mpz_tdiv_q(r._value,d._value,n._value);
}
   
inline void rem_(Rational& r, const Rational& x1, const Rational& x2) { 
  Integer q; quot_(q,x1,x2); r=x1-x2*q;
}
   
template<> inline std::string name<Rational>() { 
  return "Rational"; }
    

  }

}

