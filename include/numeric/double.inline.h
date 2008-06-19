/***************************************************************************
 *            double.inline.h
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


#include <iostream>
#include <limits>
#include "config.h"
#include "numeric/rounding.h"
#include "numeric/double.h"

#if defined ARIADNE_STANDARD_ROUNDING 

#define ARIADNE_ROUNDED_ARITHMETIC(R,X,Y) \
  template<class RM> inline void add_(R& r, const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); (volatile R&)r = (volatile X&)x + (volatile Y&)y; } \
  template<class RM> inline void sub_(R& r, const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); (volatile R&)r = (volatile X&)x - (volatile Y&)y; } \
  template<class RM> inline void mul_(R& r, const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); (volatile R&)r = (volatile X&)x * (volatile Y&)y; } \
  template<class RM> inline void div_(R& r, const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); (volatile R&)r = (volatile X&)x / (volatile Y&)y; } \
  \
  template<class RM> inline R add(const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); return (volatile X&)x + (volatile Y&)y; } \
  template<class RM> inline R sub(const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd);  return (volatile X&)x - (volatile Y&)y; } \
  template<class RM> inline R mul(const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); return (volatile X&)x * (volatile Y&)y; } \
  template<class RM> inline R div(const X& x, const Y& y, Round<RM> rnd) { \
    set_rounding_mode(rnd); return (volatile X&)x / (volatile Y&)y; } \

#elif defined ARIADNE_OPPOSITE_ROUNDING 

#ifdef ARIADNE_STD_TRANSCENDENTAL 
#error "Cannot use std transcendental functions with opposite arithmetic"
#endif

#define ARIADNE_ROUNDED_ARITHMETIC(R,X,Y) \
  inline void add_(R& r, const X& x, const Y& y, RoundApprox) { r=x+y; } \
  inline void add_(R& r, const X& x, const Y& y, RoundUp)     { r=(volatile X&)x+(volatile Y&)y; } \
  inline void add_(R& r, const X& x, const Y& y, RoundDown)   { volatile X t=-x; r=-(t-y); } \
  inline void sub_(R& r, const X& x, const Y& y, RoundApprox) { r=x-y; } \
  inline void sub_(R& r, const X& x, const Y& y, RoundUp)     { r=(volatile X&)x-(volatile Y&)y; } \
  inline void sub_(R& r, const X& x, const Y& y, RoundDown)   { volatile X t=-x; r=-(t+y); } \
  inline void mul_(R& r, const X& x, const Y& y, RoundApprox) { r=x*y; } \
  inline void mul_(R& r, const X& x, const Y& y, RoundUp)     { r=(volatile X&)x*(volatile Y&)y; } \
  inline void mul_(R& r, const X& x, const Y& y, RoundDown)   { volatile X t=-x; r=-(t*y); } \
  inline void div_(R& r, const X& x, const Y& y, RoundApprox) { r=x/y; } \
  inline void div_(R& r, const X& x, const Y& y, RoundUp)     { r=(volatile X&)x/(volatile Y&)y; } \
  inline void div_(R& r, const X& x, const Y& y, RoundDown)   { volatile X t=-x; r=-(t/y); } \
  \
  inline R add(const X& x, const Y& y, RoundApprox) { return x+y; } \
  inline R add(const X& x, const Y& y, RoundUp)     { return (volatile X&)x+(volatile Y&)y; } \
  inline R add(const X& x, const Y& y, RoundDown)   { volatile X t=-x; return -(t-y); } \
  inline R sub(const X& x, const Y& y, RoundApprox) { return x-y; } \
  inline R sub(const X& x, const Y& y, RoundUp)     { return (volatile X&)x-(volatile Y&)y; } \
  inline R sub(const X& x, const Y& y, RoundDown)   { volatile X t=-x; return -(t+y); } \
  inline R mul(const X& x, const Y& y, RoundApprox) { return x*y; } \
  inline R mul(const X& x, const Y& y, RoundUp)     { return (volatile X&)x*(volatile Y&)y; } \
  inline R mul(const X& x, const Y& y, RoundDown)   { volatile X t=-x; return -(t*y); } \
  inline R div(const X& x, const Y& y, RoundApprox) { return x/y; } \
  inline R div(const X& x, const Y& y, RoundUp)     { return (volatile X&)x/(volatile Y&)y; } \
  inline R div(const X& x, const Y& y, RoundDown)   { volatile X t=-x; return -(t/y); } \
  
#else

#error "No arithmetic mode defined"

#endif



#if defined ARIADNE_STD_TRANSCENDENTAL

#include <cmath>
#define ARIADNE_ROUNDED_FUNCTION(F) \
  template<class RM> inline void F##_(double& r, const double& x, Round<RM> rnd) { \
    set_rounding_mode(rnd); \
    (volatile double&)r = std::F( (volatile double&)x ); \
  } \

#elif defined ARIADNE_MPFR_TRANSCENDENTAL

#include <mpfr.h>
#define ARIADNE_ROUNDED_FUNCTION(F) \
  template<class RM> inline void F##_(double& r, const double& x, Round<RM> rnd) { \
    mp_rnd_t rm=mpfr_rounding_mode(rnd); \
    mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_##F(t,t,rm); r=mpfr_get_d(t,rm); \
  } \

#elif defined ARIADNE_INTERNAL_TRANSCENDENTAL

#include "numeric/transcendental.h"

#else

#error "No transcendental functions defined"

#endif




namespace Ariadne {


inline double med(double x, double y) { return (x+y)/2; }
inline double rad(double x, double y) { return (y-x)/2; }

#if defined ARIADNE_STANDARD_ROUNDING 

inline bool initialise() { return true; }

inline void incr_(double& r) {
  set_rounding_mode(round_up); (volatile double&)r += std::numeric_limits<double>::min(); }
inline void decr_(double& r) {
  set_rounding_mode(round_down); (volatile double&)r -= std::numeric_limits<double>::min(); }

inline void prec_(double& r, const double& x) {
  set_rounding_mode(round_down); (volatile double&)r = x - std::numeric_limits<double>::min(); }
inline void succ_(double& r, const double& x) {
  set_rounding_mode(round_up); (volatile double&)r = x + std::numeric_limits<double>::min(); }

inline void med_(double& r, const double& x, const double& y, RoundApprox rnd) {
  set_rounding_mode(rnd); (volatile double&)r = ( (volatile double&)x + (volatile double&)y ) / 2; }
inline  void rad_(double& r, const double& x, const double& y, RoundUp rnd) {
  set_rounding_mode(rnd); (volatile double&)r = ( (volatile double&)y - (volatile double&)x) / 2; }


#elif defined ARIADNE_OPPOSITE_ROUNDING 

inline bool initialise() { 
  //std::cerr << "Initialising rounding mode\n"; 
  fesetround(FE_UPWARD); return true; }

inline void prec_(double& r, const double& x) {
  r = -(-x+std::numeric_limits<double>::min()); }
inline void succ_(double& r, const double& x) {
  r = x + std::numeric_limits<double>::min(); }

inline void decr_(double& r) {
  r = -(-r+std::numeric_limits<double>::min()); }
inline void incr_(double& r) {
  r = r+std::numeric_limits<double>::min(); }

inline void med_(double& r, const double& x, const double& y, RoundApprox rnd) { 
  r = (x+y)/2; }

inline void rad_(double& r, const double& x, const double& y, RoundUp rnd) {
  r = (y-x)/2; }

#endif

ARIADNE_ROUNDED_ARITHMETIC(double,double,double);
ARIADNE_ROUNDED_ARITHMETIC(double,double,int);
ARIADNE_ROUNDED_ARITHMETIC(double,int,double);
ARIADNE_ROUNDED_ARITHMETIC(double,double,unsigned int);
ARIADNE_ROUNDED_ARITHMETIC(double,unsigned int,double);









// Transcendental functions

template<class RM> inline void pi_(double& r, Round<RM>) { r=3.1415926535897931; }
template<> inline void pi_(double& r, Round<Up>) { r=3.1415926535897936; }

ARIADNE_ROUNDED_FUNCTION(sqrt)
ARIADNE_ROUNDED_FUNCTION(exp)
ARIADNE_ROUNDED_FUNCTION(log)
ARIADNE_ROUNDED_FUNCTION(sin)
ARIADNE_ROUNDED_FUNCTION(cos)
ARIADNE_ROUNDED_FUNCTION(tan)
ARIADNE_ROUNDED_FUNCTION(asin)
ARIADNE_ROUNDED_FUNCTION(acos)
ARIADNE_ROUNDED_FUNCTION(atan)
ARIADNE_ROUNDED_FUNCTION(sinh)
ARIADNE_ROUNDED_FUNCTION(cosh)
ARIADNE_ROUNDED_FUNCTION(tanh)


#if defined ARIADNE_MPFR_TRANSCENDENTAL

template<class Rnd> inline void pow_(double& r, const double& x, const unsigned int& m, Rnd rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);                                \
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_pow_ui(t,t,m,rm); r=mpfr_get_d(t,rm);
}

template<class Rnd> inline void pow_(double& r, const double& x, const int& n, Rnd rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);                                \
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_pow_si(t,t,n,rm); r=mpfr_get_d(t,rm);
}

template<class Rnd> inline void hypot_(double& r, const double& x, const double& y, Rnd rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);                                \
  mpfr_t s,t; mpfr_init_set_d(s,x,GMP_RNDN); mpfr_init_set_d(t,y,GMP_RNDN); mpfr_hypot(s,s,t,rm); r=mpfr_get_d(s,rm);
}


#elif defined ARIADNE_STD_TRANSCENDENTAL

template<class Rnd> inline void pow_(double& r, const double& x, const int& n, Rnd rnd) {
  set_rounding_mode(rnd); r=std::pow(x,n); }

template<class Rnd> inline void pow_(double& r, const double& x, const unsigned int& m, Rnd rnd) {
  set_rounding_mode(rnd); int n=m; r=std::pow(x,n); }

template<class Rnd> inline void hypot_(double& r, const double& x, const double& y, Rnd rnd) {
  set_rounding_mode(rnd); double s=std::fabs(x); double t=std::fabs(y); if(s>t) { std::swap(s,t); } 
  s/=t; s*=s; s=std::sqrt(s); r=s*t; }

#elif defined ARIADNE_INTERNAL_TRANSCENDENTAL


#endif


}



#undef ARIADNE_ROUNDED_ARITHMETIC
#undef ARIADNE_ROUNDED_FUNCTION
