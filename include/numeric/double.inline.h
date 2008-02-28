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


#include <cmath>
#include <limits>
#include "numeric/rounding.h"
#include "numeric/double.h"

#define MPFR_TRANSCENDENTAL 1
#undef STD_TRANSCENDENTAL
#undef ARIADNE_TRANSCENDENTAL

#ifdef MPFR_TRANSCENDENTAL
#include <mpfr.h>
#endif

#ifdef STD_TRANSCENDENTAL
#include <cmath>
#endif

#ifdef ARIADNE_TRANSCENDENTAL
#include "numeric/transcendental.h"
#endif


namespace Ariadne {
namespace Numeric {


inline void incr_(double& r) {
  set_rounding_mode(hardware_round_up); (volatile double&)r += std::numeric_limits<double>::min(); }
inline void decr_(double& r) {
  set_rounding_mode(hardware_round_down); (volatile double&)r -= std::numeric_limits<double>::min(); }

inline void prec_(double& r, const double& x) {
  set_rounding_mode(hardware_round_down); (volatile double&)r = x - std::numeric_limits<double>::min(); }
inline void succ_(double& r, const double& x) {
  set_rounding_mode(hardware_round_up); (volatile double&)r = x + std::numeric_limits<double>::min(); }

inline void add_(double& r, const double& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x + (volatile double&)y; }
inline void sub_(double& r, const double& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x - (volatile double&)y; }
inline void mul_(double& r, const double& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x * (volatile double&)y; }
inline void div_(double& r, const double& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x / (volatile double&)y; }

inline void add_(double& r, const int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile int&)x + (volatile double&)y; }
inline void sub_(double& r, const int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile int&)x - (volatile double&)y; }
inline void mul_(double& r, const int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile int&)x * (volatile double&)y; }
inline void div_(double& r, const int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile int&)x / (volatile double&)y; }

inline void add_(double& r, const double& x, const int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x + (volatile int&)y; }
inline void sub_(double& r, const double& x, const int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x - (volatile int&)y; }
inline void mul_(double& r, const double& x, const int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x * (volatile int&)y; }
inline void div_(double& r, const double& x, const int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x / (volatile int&)y; }

inline void add_(double& r, const unsigned int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile unsigned int&)x + (volatile double&)y; }
inline void sub_(double& r, const unsigned int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile unsigned int&)x - (volatile double&)y; }
inline void mul_(double& r, const unsigned int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile unsigned int&)x * (volatile double&)y; }
inline void div_(double& r, const unsigned int& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile unsigned int&)x / (volatile double&)y; }

inline void add_(double& r, const double& x, const unsigned int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x + (volatile unsigned int&)y; }
inline void sub_(double& r, const double& x, const unsigned int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x - (volatile unsigned int&)y; }
inline void mul_(double& r, const double& x, const unsigned int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x * (volatile unsigned int&)y; }
inline void div_(double& r, const double& x, const unsigned int& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = (volatile double&)x / (volatile unsigned int&)y; }

inline void med_(double& r, const double& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = ( (volatile double&)x + (volatile double&)y ) / 2; }
inline  void rad_(double& r, const double& x, const double& y, rounding_mode rnd) {
  set_rounding_mode(rnd); (volatile double&)r = ( (volatile double&)y - (volatile double&)x) / 2; }


inline void pow_(double& r, const double& x, const unsigned int& n, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  volatile double p(x); unsigned int m=n; r=1.0; 
  while(m) { if(m%2) { r*=p; } p*=p; m/=2; }
}

inline void pow_(double& r, const double& x, const int& n, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  unsigned int m=(n>=0) ? n : -n;
  volatile double p=(n>=0) ? x : 1/x;
  r=1.0; while(m) { if(m%2) { r*=p; } p*=p; m/=2; }
}





#if defined MPFR_TRANSCENDENTAL

inline mp_rnd_t mpfr_rounding_mode(rounding_mode rnd) {
  switch(rnd) {
    case hardware_round_near: return GMP_RNDN;
    case hardware_round_down: return GMP_RNDD;
    case hardware_round_up  : return GMP_RNDU;
    case hardware_round_chop: return GMP_RNDZ;
  }
  return GMP_RNDN;
}

inline void sqrt_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,rm); mpfr_sqrt(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void hypot_(double& r, const double& x, const double& y, rounding_mode rnd) {
  if(std::fabs(x)>std::fabs(y)) { hypot_(r,y,x,rnd); return; }
  set_rounding_mode(rnd);
  double t=std::fabs(x)/std::fabs(y); t*=t; sqrt_(t,t,rnd); (volatile double&)r=t*std::fabs(y);
}

inline void exp_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,rm); mpfr_exp(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void log_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,rm); mpfr_log(t,t,rm); r=mpfr_get_d(t,rm); 
}


inline void pi_(double& r, rounding_mode rnd) {
  r = (rnd==hardware_round_up ? 3.1415926535897936 : 3.1415926535897931); }

inline void sin_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_sin(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void cos_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_cos(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void tan_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_tan(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void asin_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_asin(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void acos_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_acos(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void atan_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_atan(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void sinh_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_sinh(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void cosh_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_cosh(t,t,rm); r=mpfr_get_d(t,rm); 
}

inline void tanh_(double& r, const double& x, rounding_mode rnd) {
  mp_rnd_t rm=mpfr_rounding_mode(rnd);
  mpfr_t t; mpfr_init_set_d(t,x,GMP_RNDN); mpfr_tanh(t,t,rm); r=mpfr_get_d(t,rm); 
}

#elif defined STD_TRANSCENDENTAL

inline void sqrt_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::sqrt( (volatile double&)x ); 
}

inline void hypot_(double& r, const double& x, const double& y, rounding_mode rnd) {
  if(std::fabs(x)>std::fabs(y)) { hypot_(r,y,x,rnd); return; }
  set_rounding_mode(r-+*nd);
  volatile double t=std::fabs(x)/std::fabs(y); t*=t; t=std::sqrt(t); (volatile double&)r=t*fabs(y);
}

 


inline void exp_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::exp( (volatile double&)x ); 
}

inline void log_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::log( (volatile double&)x ); 
}


inline void pi_(double& r, rounding_mode rnd) {
  r = (rnd==hardware_round_up ? 3.1415926535897936 : 3.1415926535897931); }

inline void sin_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::sin( (volatile double&)x ); 
}

inline void cos_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::cos( (volatile double&)x ); 
}

inline void tan_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::tan( (volatile double&)x ); 
}

inline void asin_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::asin( (volatile double&)x ); 
}

inline void acos_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::acos( (volatile double&)x ); 
}

inline void atan_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::atan( (volatile double&)x ); 
}

inline void sinh_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::sinh( (volatile double&)x ); 
}

inline void cosh_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::cosh( (volatile double&)x ); 
}

inline void tanh_(double& r, const double& x, rounding_mode rnd) {
  set_rounding_mode(rnd); 
  (volatile double&)r = std::tanh( (volatile double&)x ); 
}

#endif


}
}
