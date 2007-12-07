/***************************************************************************
 *            integer.inline.h
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
 
#include "debug.h"

#define ARIADNE_DIRECT_COMPARISON(Arg1,Arg2,Cmp)                            \
  inline bool operator==(const Arg1& x1, const Arg2& x2) { return Cmp==0; } \
  inline bool operator!=(const Arg1& x1, const Arg2& x2) { return Cmp!=0; } \
  inline bool operator<=(const Arg1& x1, const Arg2& x2) { return Cmp<=0; } \
  inline bool operator>=(const Arg1& x1, const Arg2& x2) { return Cmp>=0; } \
  inline bool operator< (const Arg1& x1, const Arg2& x2) { return Cmp< 0; } \
  inline bool operator> (const Arg1& x1, const Arg2& x2) { return Cmp> 0; } 

#define ARIADNE_REVERSE_COMPARISON(Arg1,Arg2,Cmp)                            \
  inline bool operator==(const Arg1& x1, const Arg2& x2) { return Cmp==0; } \
  inline bool operator!=(const Arg1& x1, const Arg2& x2) { return Cmp!=0; } \
  inline bool operator<=(const Arg1& x1, const Arg2& x2) { return Cmp>=0; } \
  inline bool operator>=(const Arg1& x1, const Arg2& x2) { return Cmp<=0; } \
  inline bool operator< (const Arg1& x1, const Arg2& x2) { return Cmp> 0; } \
  inline bool operator> (const Arg1& x1, const Arg2& x2) { return Cmp< 0; } 

namespace Ariadne {
namespace Numeric {
	
// Functions on built-in types
//inline void add_(uint& r, const uint& x, const uint& y) { r=x+y; } 
//inline void mul_(uint& r, const uint& x, const uint& y) { r=x*y; } 
	
// These functions are not provided by GMP
inline void mpz_min_(mpz_ptr rop, mpz_srcptr op1, mpz_srcptr op2) { if(mpz_cmp(op1,op2)<=0) { mpz_set(rop,op1); } else { mpz_set(rop,op2); } }
inline void mpz_max_(mpz_ptr rop, mpz_srcptr op1, mpz_srcptr op2) { if(mpz_cmp(op1,op2)>=0) { mpz_set(rop,op1); } else { mpz_set(rop,op2); } }
inline void mpz_addsi(mpz_ptr rop, mpz_srcptr op1, long int op2) { if(op2>=0) { mpz_add_ui(rop,op1,op2); } else { mpz_sub_ui(rop,op1,(uint)(-op2)); } }
inline void mpz_subsi(mpz_ptr rop, mpz_srcptr op1, long int op2) { if(op2>=0) { mpz_sub_ui(rop,op1,op2); } else { mpz_add_ui(rop,op1,(uint)(-op2)); } }
inline void mpz_si_sub(mpz_ptr rop, long int op1, mpz_srcptr op2) { if(op1>=0) { mpz_ui_sub(rop,op1,op2); } else { mpz_add_ui(rop,op2,-op1); mpz_neg(rop,rop); } }

// Class methods
inline Integer::~Integer() { mpz_clear(this->_value); }
inline Integer::Integer() :  _value() { mpz_init_set_si(this->_value,0); }
inline Integer::Integer(const int& n) : _value() { mpz_init_set_si(this->_value,n); }
inline Integer::Integer(const unsigned int& n) : _value() { mpz_init_set_ui(this->_value,n); }
inline Integer::Integer(const mpz_class& z) : _value() { mpz_init_set(this->_value,z.get_mpz_t()); }
inline Integer::Integer(const Integer& z) :  _value() { mpz_init_set(this->_value,z._value); }
inline Integer& Integer::operator=(const int& n) { mpz_set_si(this->_value,n); return *this; }
inline Integer& Integer::operator=(const unsigned int& n) { mpz_set_ui(this->_value,n); return *this; }
inline Integer& Integer::operator=(const mpz_class& z) { mpz_set(this->_value,z.get_mpz_t()); return *this; }
inline Integer& Integer::operator=(const Integer& z) { mpz_set(this->_value,z._value); return *this; }

template<class E> inline Integer::Integer(const Expression<E>& e) :  _value() { mpz_init(this->_value); e.assign_to(*this); }
template<class E> inline Integer& Integer::operator=(const Expression<E>& e) { e.assign_to(*this); return *this; }

// Conversion to built-in types
inline void conv(int& r, const Integer& x) { r=mpz_get_si(x._value); }
inline void conv(long& r, const Integer& x) { r=mpz_get_si(x._value); }
inline Integer::operator int() const { int r=mpz_get_si(this->_value); assert(mpz_cmp_si(_value,r)==0); return r; }

// Minimum and max_imum
inline void min_(Integer& r, const Integer& x1, const Integer& x2) { mpz_min_(r._value,x1._value,x2._value); }
inline void max_(Integer& r, const Integer& x1, const Integer& x2) { mpz_max_(r._value,x1._value,x2._value); }
  
// Increment and decrement 
inline void incr(Integer& r, const Integer& x) { mpz_addsi(r._value,x._value,+1); }
inline void decr(Integer& r, const Integer& x) { mpz_addsi(r._value,x._value,-1); }

// Unary arithmetic operators
inline void pos_(Integer& r, const Integer& x) { mpz_set(r._value,x._value); }
inline void neg_(Integer& r, const Integer& x) { mpz_neg(r._value,x._value); }
inline void abs_(Integer& r, const Integer& x) { mpz_abs(r._value,x._value); }

// Binary arithmetic 
inline void add_(Integer& r, const Integer& x1, const Integer& x2) { mpz_add(r._value,x1._value,x2._value); }
inline void sub_(Integer& r, const Integer& x1, const Integer& x2) { mpz_sub(r._value,x1._value,x2._value); }
inline void mul_(Integer& r, const Integer& x1, const Integer& x2) { mpz_mul(r._value,x1._value,x2._value); }

// Mixed binary arithmetic
inline void add_(Integer& r, const Integer& x1, const int& x2) { mpz_addsi(r._value,x1._value,x2); }
inline void sub_(Integer& r, const Integer& x1, const int& x2) { mpz_subsi(r._value,x1._value,x2); }
inline void mul_(Integer& r, const Integer& x1, const int& x2) { mpz_mul_si(r._value,x1._value,x2); }
inline void add_(Integer& r, const int& x1, const Integer& x2) { mpz_addsi(r._value,x2._value,x1); }
inline void sub_(Integer& r, const int& x1, const Integer& x2) { mpz_si_sub(r._value,x1,x2._value); }
inline void mul_(Integer& r, const int& x1, const Integer& x2) { mpz_mul_si(r._value,x2._value,x1); }

// Arithmetic pow_er functions
inline void pow_(Integer& r, const Integer& x1, const uint& x2) { mpz_pow_ui(r._value,x1._value,x2); }

// Integer operations
inline void quot_(Integer& r, const Integer& x1, const Integer& x2) { mpz_tdiv_q(r._value,x1._value,x2._value); }
inline void rem_(Integer& r, const Integer& x1, const Integer& x2) { mpz_tdiv_r(r._value,x1._value,x2._value); }
inline void lcm_(Integer& r, const Integer& x1, const Integer& x2) { mpz_lcm(r._value,x1._value,x2._value); }
inline void gcd_(Integer& r, const Integer& x1, const Integer& x2) { mpz_gcd(r._value,x1._value,x2._value); }
inline void fac_(Integer& r, unsigned long int x1) { mpz_fac_ui(r._value,x1); }
inline void bin_(Integer& r, const Integer& x1, unsigned long int x2) { mpz_bin_ui(r._value,x1._value,x2); }
inline void bin_(Integer& r, unsigned long int x1, unsigned long int x2) { mpz_bin_uiui(r._value,x1,x2); }

inline int cmp(const Integer& x, const Integer& y) {
  return mpz_cmp(x._value,y._value); }
template<class E> inline int cmp(const Integer& z, const Expression<E>& e) {
  return cmp(z,Integer(e)); }
template<class E> inline int cmp(const Expression<E>& e, const Integer& z) {
  return cmp(Integer(e),z); }

// Comparison operators
ARIADNE_DIRECT_COMPARISON(Integer,Integer,mpz_cmp(x1._value,x2._value));
ARIADNE_DIRECT_COMPARISON(Integer,int,mpz_cmp_si(x1._value,x2));
ARIADNE_REVERSE_COMPARISON(int,Integer,mpz_cmp_si(x2._value,x1));
ARIADNE_DIRECT_COMPARISON(Integer,uint,mpz_cmp_si(x1._value,x2));
ARIADNE_REVERSE_COMPARISON(uint,Integer,mpz_cmp_si(x2._value,x1));
ARIADNE_DIRECT_COMPARISON(Integer,long,mpz_cmp_si(x1._value,x2));
ARIADNE_REVERSE_COMPARISON(long,Integer,mpz_cmp_si(x2._value,x1));
ARIADNE_DIRECT_COMPARISON(Integer,ulong,mpz_cmp_ui(x1._value,x2));
ARIADNE_REVERSE_COMPARISON(ulong,Integer,mpz_cmp_ui(x2._value,x1));
ARIADNE_DIRECT_COMPARISON(Integer,double,mpz_cmp_d(x1._value,x2));
ARIADNE_REVERSE_COMPARISON(double,Integer,mpz_cmp_d(x2._value,x1));

template<class E> 
inline bool operator==(const Integer& x, const Expression<E>& y) {
  return x==Integer(y); }

template<class E> 
inline bool operator==(const Expression<E>& x, const Integer& y) {
  return Integer(x)==y; }

inline Integer operator%(const Integer& n1, const Integer& n2) { Integer r; rem_(r,n1,n2); return r; }
inline Integer& operator++(Integer& n) { add_(n,n,1); return n; }
inline Integer& operator--(Integer& n) { sub_(n,n,1); return n; }

static uint factorials[13]={ 
	  1, 1, 2, 6, 
	  24, 120, 720, 5040, 
	  40320, 362880, 3628800, 39916800, 
	  479001600 
};

template<class R, class A> inline
void factorial_(R& r, const A& n) {
  if(n<=1) { r=1; return; }
  if(n<=12) { int m(n); r=factorials[m]; return; }
  A i=n; r=479001600;
  while(i!=12) { r*=i; --i; }	
}


template<class R, class A1, class A2> inline 
void choose_(R& r, const A1& n, const A2& k) 
{
  //std::cerr << "choose(" << n << "," << k << ")=" << std::flush;
  if(k==0 || k==n) { r=1; return; }
  if(k<0 || k>n) { r=0; return; }
  A2 m=(n-k < k) ? k : static_cast<A2>(n-k);
  R result=1;
  for(A1 i=n; i!=n-m; --i) { result*=i; }
  for(A1 i=m; i!=1; --i) { result/=i; }
  //std::cerr << result << std::endl;
  r=result;
}

   
   

template<class R, class A1, class A2> inline 
void gcd_(R& r, const A1& a, const A2& b) {
  R aa=a; R bb=b; R cc=aa%bb;
  while(cc!=0) { aa=bb; bb=cc; cc=aa%bb; }
  return bb;
}

template<class R, class A1, class A2> inline 
void lcm_(R& r, const A1& a, const A2& b) {
  R res; quot(res,a*b,gcd(a,b)); return res;
}


template<class N> inline
N exp2(const N& n) {
  return 1<<n;
}


template<class N> inline
N log2_floor(const N& n) {
  if(n<1) { throw std::invalid_argument(__PRETTY_FUNCTION__); }
  N r=0;
  N y=n;
  while(y>=n) {
    y/=2;
    r+=1;
  }
  return r;
}

template<class N> inline
N log2_ceil(const N& n) {
  if(n<1) { throw std::invalid_argument(__PRETTY_FUNCTION__); }
  N r=0;
  N y=n;
  while(y>1) {
    y/=2;
    r+=1;
  }
  return r;
}




inline uint factorial(uint n) { 
	ARIADNE_ASSERT(n<13) 
	return factorials[n];
}

inline uint choose(uint n, uint k) { 
	ARIADNE_ASSERT(n<13);
	ARIADNE_ASSERT(k<=n);
	return factorial(n)/(factorial(k)*factorial(uint(n-k)));
}

inline int choose(int n, int k) { return choose(uint(n),uint(k)); }

template<class N> inline N factorial(uint n) {
  N r; factorial_(r,n); return r; }
template<class N> inline N choose(uint n, uint k) {
  N r; choose_(r,n,k); return r; }


}}

#undef ARIADNE_DIRECT_COMPARISON
#undef ARIADNE_REVERSE_COMPARISON
