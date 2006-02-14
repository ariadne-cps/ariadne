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
#include <boost/numeric/interval/io.hpp>

#include "numerical_type.h"

namespace Ariadne {
  template<typename R> class rounding;
  
  template<>
  struct rounding<double> :
    boost::numeric::interval_lib::rounded_arith_opp<double>
  {
   private:
    typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
  
    double invoke_mpfr(double x, mpfr_func f, mp_rnd_t r) {
      mpfr_t xx;
      mpfr_init_set_d(xx, x, GMP_RNDN);
      f(xx, xx, r);
      double res = mpfr_get_d(xx, r);
      mpfr_clear(xx);
      return res;
    }
   public:
    double exp_down(double x) { return invoke_mpfr(x, mpfr_exp, GMP_RNDD); } 
    double log_down(double x) { return invoke_mpfr(x, mpfr_log, GMP_RNDD); } 
    double sin_down(double x) { return invoke_mpfr(x, mpfr_sin, GMP_RNDD); } 
    double cos_down(double x) { return invoke_mpfr(x, mpfr_cos, GMP_RNDD); } 
    double tan_down(double x) { return invoke_mpfr(x, mpfr_tan, GMP_RNDD); } 
    double asin_down(double x) { return invoke_mpfr(x, mpfr_asin, GMP_RNDD); } 
    double acos_down(double x) { return invoke_mpfr(x, mpfr_acos, GMP_RNDD); } 
    double atan_down(double x) { return invoke_mpfr(x, mpfr_atan, GMP_RNDD); } 
    double sinh_down(double x) { return invoke_mpfr(x, mpfr_sinh, GMP_RNDD); } 
    double cosh_down(double x) { return invoke_mpfr(x, mpfr_cosh, GMP_RNDD); } 
    double tanh_down(double x) { return invoke_mpfr(x, mpfr_tanh, GMP_RNDD); }
    double asinh_down(double x) { return invoke_mpfr(x, mpfr_asinh, GMP_RNDD); } 
    double acosh_down(double x) { return invoke_mpfr(x, mpfr_acosh, GMP_RNDD); } 
    double atanh_down(double x) { return invoke_mpfr(x, mpfr_atanh, GMP_RNDD); }
  
    double exp_up (double x) { return invoke_mpfr(x, mpfr_exp, GMP_RNDU); }
    double log_up (double x) { return invoke_mpfr(x, mpfr_log, GMP_RNDU); }
    double sin_up (double x) { return invoke_mpfr(x, mpfr_sin, GMP_RNDU); }
    double cos_up (double x) { return invoke_mpfr(x, mpfr_cos, GMP_RNDU); }
    double tan_up (double x) { return invoke_mpfr(x, mpfr_tan, GMP_RNDU); }
    double asin_up (double x) { return invoke_mpfr(x, mpfr_asin, GMP_RNDU); }
    double acos_up (double x) { return invoke_mpfr(x, mpfr_acos, GMP_RNDU); }
    double atan_up (double x) { return invoke_mpfr(x, mpfr_atan, GMP_RNDU); }
    double sinh_up (double x) { return invoke_mpfr(x, mpfr_sinh, GMP_RNDU); }
    double cosh_up (double x) { return invoke_mpfr(x, mpfr_cosh, GMP_RNDU); }
    double tanh_up (double x) { return invoke_mpfr(x, mpfr_tanh, GMP_RNDU); }
    double asinh_up (double x) { return invoke_mpfr(x, mpfr_asinh, GMP_RNDU); }
    double acosh_up (double x) { return invoke_mpfr(x, mpfr_acosh, GMP_RNDU); }
    double atanh_up (double x) { return invoke_mpfr(x, mpfr_atanh, GMP_RNDU); }
  };
  
  template<>
  struct rounding<Dyadic> :
    boost::numeric::interval_lib::rounded_arith_opp<double>
  {
   public:
    Dyadic exp_down(Dyadic x); 
    Dyadic log_down(Dyadic x);
    Dyadic sin_down(Dyadic x);
    Dyadic cos_down(Dyadic x);
    Dyadic tan_down(Dyadic x);
    Dyadic asin_down(Dyadic x);
    Dyadic atan_down(Dyadic x);
    Dyadic sinh_down(Dyadic x);
    Dyadic cosh_down(Dyadic x);
    Dyadic tanh_down(Dyadic x);
    Dyadic asinh_down(Dyadic x);
    Dyadic acosh_down(Dyadic x);
    Dyadic atanh_down(Dyadic x);
  
    Dyadic exp_up (Dyadic x);
    Dyadic log_up (Dyadic x);
    Dyadic sin_up (Dyadic x);
    Dyadic cos_up (Dyadic x);
    Dyadic tan_up (Dyadic x);
    Dyadic asin_up (Dyadic x); 
    Dyadic acos_up (Dyadic x);
    Dyadic atan_up (Dyadic x); 
    Dyadic sinh_up (Dyadic x);
    Dyadic cosh_up (Dyadic x);
    Dyadic tanh_up (Dyadic x); 
    Dyadic asinh_up (Dyadic x);
    Dyadic acosh_up (Dyadic x);
    Dyadic atanh_up (Dyadic x);
  };

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
