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
 
#include "macros/assert.h"
#include "macros/throw.h"
#include "numeric/exceptions.h"

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

	
using std::min;
using std::max;

// Class methods
inline Integer::~Integer() { 
  mpz_clear(this->_value); }
inline Integer::Integer() :  _value() {
  mpz_init_set_si(this->_value,0); }
inline Integer::Integer(const int& n) : _value() {
  mpz_init_set_si(this->_value,n); }
inline Integer::Integer(const long int& n) : _value() {
  mpz_init_set_si(this->_value,n); }
inline Integer::Integer(const unsigned int& n) : _value() {
  mpz_init_set_ui(this->_value,n); }
inline Integer::Integer(const unsigned long int& n) : _value() {
  mpz_init_set_ui(this->_value,n); }
inline Integer::Integer(mpz_srcptr z) : _value() {
  mpz_init_set(this->_value,z); }
inline Integer::Integer(const Integer& z) :  _value() {
  mpz_init_set(this->_value,z._value); }
inline Integer& Integer::operator=(const int& n) {
  mpz_set_si(this->_value,n); return *this; }
inline Integer& Integer::operator=(const long int& n) {
  mpz_set_si(this->_value,n); return *this; }
inline Integer& Integer::operator=(const unsigned int& n) {
  mpz_set_ui(this->_value,n); return *this; }
inline Integer& Integer::operator=(const unsigned long int& n) {
  mpz_set_ui(this->_value,n); return *this; }
inline Integer& Integer::operator=(mpz_srcptr z) {
  mpz_set(this->_value,z); return *this; }
inline Integer& Integer::operator=(const Integer& z) {
  mpz_set(this->_value,z._value); return *this; }

template<class E> inline Integer::Integer(const Expression<E>& e)
  :  _value() { mpz_init(this->_value); e.assign_to(*this); }
template<class E> inline Integer& Integer::operator=(const Expression<E>& e) {
  e.assign_to(*this); return *this; }

inline Integer::operator long int() const { 
  ARIADNE_ASSERT(mpz_fits_sint_p(this->_value));
  return mpz_get_si(this->_value); }





// Comparison operators
inline int cmp(const Integer& x, const Integer& y) {
  return mpz_cmp(x._value,y._value); }
template<class E> inline int cmp(const Integer& z, const Expression<E>& e) {
  return cmp(z,Integer(e)); }
template<class E> inline int cmp(const Expression<E>& e, const Integer& z) {
  return cmp(Integer(e),z); }

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
  return cmp(x,y)==0; }
template<class E>
inline bool operator!=(const Integer& x, const Expression<E>& y) {
  return cmp(x,y)!=0; }
template<class E>
inline bool operator>=(const Integer& x, const Expression<E>& y) {
  return cmp(x,y)>=0; }
template<class E>
inline bool operator<=(const Integer& x, const Expression<E>& y) {
  return cmp(x,y)<=0; }
template<class E>
inline bool operator> (const Integer& x, const Expression<E>& y) {
  return cmp(x,y)> 0; }
template<class E>
inline bool operator< (const Integer& x, const Expression<E>& y) {
  return cmp(x,y)< 0; }

template<class E>
inline bool operator==(const Expression<E>& x, const Integer& y) {
  return cmp(x,y)==0; }
template<class E>
inline bool operator!=(const Expression<E>& x, const Integer& y) {
  return cmp(x,y)!=0; }
template<class E>
inline bool operator>=(const Expression<E>& x, const Integer& y) {
  return cmp(x,y)>=0; }
template<class E>
inline bool operator<=(const Expression<E>& x, const Integer& y) {
  return cmp(x,y)<=0; }
template<class E>
inline bool operator> (const Expression<E>& x, const Integer& y) {
  return cmp(x,y)> 0; }
template<class E>
inline bool operator< (const Expression<E>& x, const Integer& y) {
  return cmp(x,y)< 0; }



// Conversion to built-in types
inline void set_(int& r, const Integer& x) { 
  ARIADNE_ASSERT(mpz_fits_sint_p(x._value));
  r=mpz_get_si(x._value); }

// Conversion from built-in types
inline void set_(Integer& r, const int& x) { 
  mpz_set_si(r._value,x); }
inline void set_(Integer& r, const uint& x) { 
  mpz_set_ui(r._value,x); }


// Increment and decrement 
inline void incr_(Integer& r, const Integer& x) { 
  mpz_add_ui(r._value,x._value,1u); }
inline void decr_(Integer& r, const Integer& x) { 
  mpz_sub_ui(r._value,x._value,1u); }


// Minimum and maximum
inline void min_(Integer& r, const Integer& x1, const Integer& x2) { 
  r=std::min(x1,x2); }
inline void max_(Integer& r, const Integer& x1, const Integer& x2) { 
  r=std::max(x1,x2); }
  
// Unary arithmetic operators
inline void pos_(Integer& r, const Integer& x) { 
  mpz_set(r._value,x._value); }
inline void neg_(Integer& r, const Integer& x) { 
  mpz_neg(r._value,x._value); }
inline void abs_(Integer& r, const Integer& x) { 
  mpz_abs(r._value,x._value); }

// Binary arithmetic 
inline void add_(Integer& r, const Integer& x1, const Integer& x2) {
  mpz_add(r._value,x1._value,x2._value); }
inline void sub_(Integer& r, const Integer& x1, const Integer& x2) { 
  mpz_sub(r._value,x1._value,x2._value); }
inline void mul_(Integer& r, const Integer& x1, const Integer& x2) { 
  mpz_mul(r._value,x1._value,x2._value); }

// Mixed binary arithmetic
inline void add_(Integer& r, const Integer& x1, const uint& x2) { 
  mpz_add_ui(r._value,x1._value,x2); }
inline void sub_(Integer& r, const Integer& x1, const uint& x2) { 
  mpz_sub_ui(r._value,x1._value,x2); } 
inline void mul_(Integer& r, const Integer& x1, const uint& x2) { 
  mpz_mul_ui(r._value,x1._value,x2); }

inline void add_(Integer& r, const uint& x1, const Integer& x2) { 
  mpz_add_ui(r._value,x2._value,x1); }
inline void sub_(Integer& r, const uint& x1, const Integer& x2) { 
  mpz_ui_sub(r._value,x1,x2._value); }
inline void mul_(Integer& r, const uint& x1, const Integer& x2) { 
  mpz_mul_ui(r._value,x2._value,x1); }

inline void add_(Integer& r, const Integer& x1, const int& x2) { 
  if(x2>=0) { mpz_add_ui(r._value,x1._value,x2); } 
  else { mpz_sub_ui(r._value,x1._value,-x2); } }
inline void sub_(Integer& r, const Integer& x1, const int& x2) { 
  if(x2>=0) { mpz_sub_ui(r._value,x1._value,x2); } 
  else { mpz_add_ui(r._value,x1._value,-x2); } }
inline void mul_(Integer& r, const Integer& x1, const int& x2) { 
  mpz_mul_si(r._value,x1._value,x2); }

inline void add_(Integer& r, const int& x1, const Integer& x2) { 
  if(x1>=0) { mpz_add_ui(r._value,x2._value,x1); }
  else { mpz_sub_ui(r._value,x2._value,-x1); } }
inline void sub_(Integer& r, const int& x1, const Integer& x2) { 
  if(x1>=0) { mpz_ui_sub(r._value,x1,x2._value); }
  else { mpz_ui_sub(r._value,-x1,x2._value); } }
inline void mul_(Integer& r, const int& x1, const Integer& x2) { 
  mpz_mul_si(r._value,x2._value,x1); }

// Arithmetic power functions
inline void pow_(Integer& r, const Integer& x1, const uint& x2) { 
  mpz_pow_ui(r._value,x1._value,x2); }

// Integer operations
inline void quot_(Integer& r, const Integer& x1, const Integer& x2) { 
  mpz_tdiv_q(r._value,x1._value,x2._value); }
inline void quot_(Integer& r, const Integer& x1, const unsigned long int& x2) { 
  mpz_tdiv_q_ui(r._value,x1._value,x2); }
inline void quot_(Integer& r, const long int& x1, const unsigned long int& x2) { 
  r=x1/x2; }
inline void quot_(int& r, const int& x1, const int& x2) { 
  r=x1/x2; }
inline void quot_(long int& r, const long int& x1, const long int& x2) { 
  r=x1/x2; }
inline void quot_(unsigned int& r, const unsigned int& x1, const unsigned int& x2) { 
  r=x1/x2; }
inline void quot_(unsigned long int& r, const unsigned long int& x1, const unsigned long int& x2) { 
  r=x1/x2; }
inline void rem_(Integer& r, const Integer& x1, const Integer& x2) { 
  mpz_tdiv_r(r._value,x1._value,x2._value); }

inline void lcm_(Integer& r, const Integer& x1, const Integer& x2) { 
  mpz_lcm(r._value,x1._value,x2._value); }
inline void gcd_(Integer& r, const Integer& x1, const Integer& x2) { 
  mpz_gcd(r._value,x1._value,x2._value); }

inline void fac_(Integer& r, Integer& x) { 
  ARIADNE_ASSERT(mpz_fits_ulong_p(x._value)); 
  mpz_fac_ui(r._value,mpz_get_ui(x._value)); } 
inline void fac_(Integer& r, unsigned long int x1) { 
  mpz_fac_ui(r._value,x1); }

inline void bin_(Integer& r, const Integer& x1, const Integer& x2) { 
  ARIADNE_ASSERT(mpz_fits_ulong_p(x2._value)); 
  mpz_bin_ui(r._value,x1._value,mpz_get_ui(x2._value)); }
inline void bin_(Integer& r, const Integer& x1, unsigned long int x2) { 
  mpz_bin_ui(r._value,x1._value,x2); }
inline void bin_(Integer& r, unsigned long int x1, unsigned long int x2) { 
  mpz_bin_uiui(r._value,x1,x2); }





static uint factorials[13]={ 
  1, 1, 2, 6, 24, 120, 720, 5040, 
  40320, 362880, 3628800, 39916800, 479001600 
};


template<class R, class A> inline
void exp2_(R& r, const A& n) {
  r = 1<<n; }

template<class R, class N> 
void fac_(R& r, const N& n);

template<class R, class N1, class N2> 
void bin_(R& r, const N1& n, const N2& k); 

template<class R, class N1, class N2>  
void gcd_(R& r, const N1& a, const N2& b);

template<class R, class N1, class N2>  
void lcm_(R& r, const N1& a, const N2& b);

template<class R, class N> 
void log2_floor_(R& r, const N& n);

template<class R, class N> 
  void log2_ceil_(R& r, const N& n);


inline int fac(int n) { 
#ifndef NDEBUG 
  if(n>=13) {     
    ARIADNE_THROW(OverflowException,"int fac(int n)"," with n="<<n);
    //ARIADNE_THROW(OverflowException,__FUNCTION__," with n="<<n);
  }
#endif
  return factorials[n]; 
}

inline uint fac(uint n) { 
#ifndef NDEBUG 
  if(n>=13) { 
    ARIADNE_THROW(OverflowException,"uint fac(uint n)"," with n="<<n);
    //ARIADNE_THROW(OverflowException,__FUNCTION__," with n="<<n);
  }
#endif
  return factorials[n]; 
}

inline unsigned long int fac(unsigned long int n) { 
#ifndef NDEBUG 
  if(n>=13) { 
    ARIADNE_THROW(OverflowException,"ulong fac(ulong n)"," with n="<<n);
    //ARIADNE_THROW(OverflowException,__FUNCTION__," with n="<<n);
  }
#endif
  return factorials[n]; 
}


inline int bin(int n, int k) { 
  static const int kmax=16;
  static const int nmax[kmax]={ 
    2147483647, 2147483647, 46341, 1626, 338, 140, 82, 58, 
    46, 39, 35, 33, 31, 30, 30, 29 };
  if(k>n || k<0) { return 0; }
  int c=std::min(k,n-k);
  if(c>=kmax || n>nmax[c]) {
    ARIADNE_THROW(OverflowException,"int bin(int n, int k)","with n="<<n<<", k="<<k);
  }
  uint r; 
  bin_(r,n,c); 
  return r; 
}

inline uint bin(unsigned int n, unsigned int k) { 
  static const uint kmax=16u;
  static const uint nmax[kmax]={ 
    4294967295u, 4294967295u, 65536u, 2049u, 402u, 161u, 92u, 63u, 
    49u, 42u, 37u, 34u, 33u, 31u, 31u, 30u };
  if(k>n) { return 0; }
  uint c=std::min(k,n-k);
  
  if(c>=kmax || n>nmax[c]) {
    ARIADNE_THROW(OverflowException,"uint bin(uint n,uint k)","with n="<<n<<", k="<<k);
  }
  
  uint r; 
  bin_(r,n,c); 
  return r; 
}

inline long unsigned int bin(long unsigned int n, long unsigned int k) { 
  static const unsigned long int kmax=16u;
  static const unsigned long int nmax[kmax]={ 
    4294967295u, 4294967295u, 65536u, 2049u, 402u, 161u, 92u, 63u, 
    49u, 42u, 37u, 34u, 33u, 31u, 31u, 30u };
  /*
  static const unsigned long int kmax=32;
  static const unsigned long int nmax[kmax]={
    18446744073709551615ul, 4294967296ul, 3329022ul, 102570ul, 13467ul, 3612ul, 1449ul, 746ul, 
    453ul, 308ul, 227ul, 178ul, 147ul, 125ul, 110ul, 99ul, 
    90ul, 84ul, 79ul, 75ul, 72ul, 69ul, 68ul, 66ul, 
    65ul, 64ul, 63ul, 63ul, 62ul, 62ul, 62ul, 62ul };
  */
  if(k>n) { return 0; }
  long unsigned int c=std::min(k,n-k);
  if(c>=kmax || n>nmax[c]) {
    ARIADNE_THROW(OverflowException,"ulong bin(ulong n, ulong k)","with n="<<n<<", k="<<k);
  }
  long unsigned int r; 
  bin_(r,n,c); 
  return r; 
}

inline uint exp2(const uint& n) {
  ARIADNE_ASSERT(n<32);
  uint r; exp2_(r,n); return r; }
inline uint log2_floor(const uint& n) {
  uint r; log2_floor_(r,n); return r; }
inline uint log2_ceil(const uint& n) {
  uint r; log2_ceil_(r,n); return r; }


inline Integer& operator++(Integer& r) {
  incr_(r,r); return r; }
inline Integer& operator--(Integer& r) {
  decr_(r,r); return r; }

template<> inline std::string name<Integer>() { 
  return "Integer"; }

} // namespace Ariadne

#undef ARIADNE_DIRECT_COMPARISON
#undef ARIADNE_REVERSE_COMPARISON
