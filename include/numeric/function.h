/***************************************************************************
 *            function.h
 *
 *  Wed 18 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
  
/*! \file function.h
 *  \brief Standard functions on double precision and dyadic number types.
 */

#ifndef _ARIADNE_FUNCTION_H
#define _ARIADNE_FUNCTION_H

#include <cmath>

#include <gmpxx.h>
#include <mpfr.h>

#include <boost/numeric/interval.hpp>

#include "../numeric/numerical_types.h"
#include "../numeric/approximation.h"

namespace Ariadne {
  template<typename R> class rounding;
  
  template<>
  struct rounding<Float64> :
    boost::numeric::interval_lib::rounded_arith_opp<Float64>
  {
   private:
    typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
  
    Float64 invoke_mpfr(Float64 x, mpfr_func f, mp_rnd_t r) {
      mpfr_t xx;
      mpfr_init_set_d(xx, x, GMP_RNDN);
      f(xx, xx, r);
      Float64 res = mpfr_get_d(xx, r);
      mpfr_clear(xx);
      return res;
    }
   public:
    Float64 exp_down(Float64 x) { return invoke_mpfr(x, mpfr_exp, GMP_RNDD); } 
    Float64 log_down(Float64 x) { return invoke_mpfr(x, mpfr_log, GMP_RNDD); } 
    Float64 sin_down(Float64 x) { return invoke_mpfr(x, mpfr_sin, GMP_RNDD); } 
    Float64 cos_down(Float64 x) { return invoke_mpfr(x, mpfr_cos, GMP_RNDD); } 
    Float64 tan_down(Float64 x) { return invoke_mpfr(x, mpfr_tan, GMP_RNDD); } 
    Float64 asin_down(Float64 x) { return invoke_mpfr(x, mpfr_asin, GMP_RNDD); } 
    Float64 acos_down(Float64 x) { return invoke_mpfr(x, mpfr_acos, GMP_RNDD); } 
    Float64 atan_down(Float64 x) { return invoke_mpfr(x, mpfr_atan, GMP_RNDD); } 
    Float64 sinh_down(Float64 x) { return invoke_mpfr(x, mpfr_sinh, GMP_RNDD); } 
    Float64 cosh_down(Float64 x) { return invoke_mpfr(x, mpfr_cosh, GMP_RNDD); } 
    Float64 tanh_down(Float64 x) { return invoke_mpfr(x, mpfr_tanh, GMP_RNDD); }
    Float64 asinh_down(Float64 x) { return invoke_mpfr(x, mpfr_asinh, GMP_RNDD); } 
    Float64 acosh_down(Float64 x) { return invoke_mpfr(x, mpfr_acosh, GMP_RNDD); } 
    Float64 atanh_down(Float64 x) { return invoke_mpfr(x, mpfr_atanh, GMP_RNDD); }
  
    Float64 exp_up (Float64 x) { return invoke_mpfr(x, mpfr_exp, GMP_RNDU); }
    Float64 log_up (Float64 x) { return invoke_mpfr(x, mpfr_log, GMP_RNDU); }
    Float64 sin_up (Float64 x) { return invoke_mpfr(x, mpfr_sin, GMP_RNDU); }
    Float64 cos_up (Float64 x) { return invoke_mpfr(x, mpfr_cos, GMP_RNDU); }
    Float64 tan_up (Float64 x) { return invoke_mpfr(x, mpfr_tan, GMP_RNDU); }
    Float64 asin_up (Float64 x) { return invoke_mpfr(x, mpfr_asin, GMP_RNDU); }
    Float64 acos_up (Float64 x) { return invoke_mpfr(x, mpfr_acos, GMP_RNDU); }
    Float64 atan_up (Float64 x) { return invoke_mpfr(x, mpfr_atan, GMP_RNDU); }
    Float64 sinh_up (Float64 x) { return invoke_mpfr(x, mpfr_sinh, GMP_RNDU); }
    Float64 cosh_up (Float64 x) { return invoke_mpfr(x, mpfr_cosh, GMP_RNDU); }
    Float64 tanh_up (Float64 x) { return invoke_mpfr(x, mpfr_tanh, GMP_RNDU); }
    Float64 asinh_up (Float64 x) { return invoke_mpfr(x, mpfr_asinh, GMP_RNDU); }
    Float64 acosh_up (Float64 x) { return invoke_mpfr(x, mpfr_acosh, GMP_RNDU); }
    Float64 atanh_up (Float64 x) { return invoke_mpfr(x, mpfr_atanh, GMP_RNDU); }
  };
  
  template<>
  struct rounding<MPFloat> :
    boost::numeric::interval_lib::rounded_arith_opp<Float64>
  {
   public:
    MPFloat exp_down(MPFloat x); 
    MPFloat log_down(MPFloat x);
    MPFloat sin_down(MPFloat x);
    MPFloat cos_down(MPFloat x);
    MPFloat tan_down(MPFloat x);
    MPFloat asin_down(MPFloat x);
    MPFloat atan_down(MPFloat x);
    MPFloat sinh_down(MPFloat x);
    MPFloat cosh_down(MPFloat x);
    MPFloat tanh_down(MPFloat x);
    MPFloat asinh_down(MPFloat x);
    MPFloat acosh_down(MPFloat x);
    MPFloat atanh_down(MPFloat x);
  
    MPFloat exp_up (MPFloat x);
    MPFloat log_up (MPFloat x);
    MPFloat sin_up (MPFloat x);
    MPFloat cos_up (MPFloat x);
    MPFloat tan_up (MPFloat x);
    MPFloat asin_up (MPFloat x); 
    MPFloat acos_up (MPFloat x);
    MPFloat atan_up (MPFloat x); 
    MPFloat sinh_up (MPFloat x);
    MPFloat cosh_up (MPFloat x);
    MPFloat tanh_up (MPFloat x); 
    MPFloat asinh_up (MPFloat x);
    MPFloat acosh_up (MPFloat x);
    MPFloat atanh_up (MPFloat x);
  };

  inline MPFloat div_approx(const MPFloat& x1, const MPFloat& x2, const MPFloat& e) {
    return div_approx(x1,x2,precision(e));
  }
    
  inline MPFloat div_approx(const MPFloat& x1, const MPFloat& x2, const uint& n) {
    MPFloat r(0,n);
    mpf_div(r.get_mpf_t(),x1.get_mpf_t(),x2.get_mpf_t());
    return r;
  }
 
  inline Dyadic div_approx(const Dyadic& x1, const Dyadic& x2, const Dyadic& e) {
    Rational q=Rational(x1)/Rational(x2);
    return Numeric::approximate<Dyadic>(q,e);
  }
    
  template<typename R> 
  R sqrt_approx(R x, R e) {
    assert(false);
  }

  template<typename R> 
  R exp_approx(R x, R e) {
    R result=0;
    R term=1;
    R error=term;
    R n=0;    
    while(error>e) {
      result+=term;
      n+=1;
      term*=x/n;
      error=term;
    }
    return result;
  }
  
  template<typename R> 
  R exp_down(R x, R e) {
    R result=0;
    R term=1;
    R error=term;
    R n=0;    
    while(error>e && x>n) {
      result+=term;
      n+=1;
      term*=x/n;
      error=term;
    }
    if(term<0) {
      result+=term;
    }
    return result;
  }
  
  template<typename R> 
  R exp_up(R x, R e) {
    R result=0;
    R term=1;
    R error=term;
    R n=0;    
    while(error>e && x>n) {
      result+=term;
      n+=1;
      term*=x/n;
      error=term;
    }
    if(term>0) {
      result+=term;
    }
    return result;
  }

  template<typename R> 
  R cos_approx(R x, R e) {
    R result=0;
    R term=1;
    R error=term;
    R n=0;
    
    while(error>e) {
      result+=term;
      n+=2;
      term*=-x*x/(n*(n-1));
      error=term;
    }
    
    return result;
  }
  
  template<typename R> 
  R sin_approx(R x, R e) {
    R result=0;
    R term=x;
    R error=term;
    R n=1;
    
    while(error>e) {
      result+=term;
      n+=2;
      term*=-x*x/(n*(n-1));
      error=term;
    }
    
    return result;
  }

}

#endif /* _ARIADNE_FUNCTION_H */
