/***************************************************************************
 *            mofloat.h
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
 
/*! \file mpfloat.h
 *  \brief Type definitions and conversion operators for multiple-precision floating-point numbers.
 */

#ifndef _ARIADNE_MPFLOAT_H
#define _ARIADNE_MPFLOAT_H

#include <mpfr.h>
// Don't use mpf2mpfr, since it ignores rounding modes
// #include <mpf2mpfr.h>

#include <gmpxx.h>

#include "../declarations.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/function.h"
#include "../numeric/rational.h"

namespace Ariadne {
  namespace Numeric {
    #ifdef DOXYGEN
    /*! \brief A multiple-precision floating-point type.
     *  \ingroup Numeric
     * 
     * Currently implemented using mpf_class from the GNU Multiple Precision library.
     */
    class MPFloat { };
#else
    typedef mpf_class MPFloat;
#endif

/*
    class MPFloat {
     private:
      mpfr_t _value;
      typedef __mpfr_struct* pointer;
      typedef const __mpfr_struct* const_pointer;
     public:
      MPFloat() { mpfr_init(_value); }
      MPFloat(int x) { mpfr_init_set_si(_value,x,GMP_RNDN); }
      MPFloat(long int x) { mpfr_init_set_si(_value,x,GMP_RNDN); }
      MPFloat(long unsigned int x) { mpfr_init_set_ui(_value,x,GMP_RNDN); }
      MPFloat(double x) { mpfr_init_set_d(_value,x,GMP_RNDN); }
      MPFloat(const MPFloat& x) { mpfr_init_set(_value,x._value,GMP_RNDN); }
      MPFloat(const mpf_class& x) { mpfr_init_set_f(_value,x.get_mpf_t(),GMP_RNDN); }
      MPFloat(const mpf_t x) { mpfr_init_set_f(_value,x,GMP_RNDN); }
      MPFloat& operator=(const MPFloat& x) { if(this != &x) { mpfr_set(_value,x._value,GMP_RNDN); } return *this; }
      mpfr_srcptr get_mpfr_t() const { return _value; }
      mpfr_ptr get_mpfr_t() { return _value; }
      mpf_class get_mpf_class() const { mpf_class result; mpfr_get_f(result.get_mpf_t(),_value,GMP_RNDD); return result; }
      double get_d() const { return mpfr_get_d(_value,GMP_RNDN); }
      operator mpf_class () const { mpf_class result; mpfr_get_f(result.get_mpf_t(),_value,GMP_RNDD); return result; }
      operator const_pointer () const { return _value; }
      operator pointer () { return _value; }
    };
    
    inline std::ostream& operator<<(std::ostream& os, const MPFloat& x) {
      return os<<x.get_d(); }
    inline std::istream& operator>>(std::istream& is, MPFloat& x) {
      //FIXME: Improve this!
      double dx; is >> dx; x=MPFloat(dx); return is; }
*/
    
    inline int precision(const MPFloat& num) {
      return mpf_get_prec(num.get_mpf_t());
    }
  
    inline int exponent(const MPFloat& num) {
      long int res; 
      mpf_get_d_2exp(&res,num.get_mpf_t()); 
      return res-1; 
    }
  
    inline MPFloat mantissa(const MPFloat& num) {
      long int exp; 
      MPFloat res; 
      mpf_get_d_2exp(&exp,num.get_mpf_t()); 
      exp=exp-1;
      mpf_div_2exp(res.get_mpf_t(),num.get_mpf_t(),exp); 
      return res; 
    }
  
    template<> class numerical_traits<MPFloat> {
     public:
      typedef ring_tag algebraic_category;
      typedef Rational field_extension_type;
    };
      
    template<> inline std::string name<Numeric::MPFloat>() { return "MPFloat"; }
    template<> inline std::string name<Numeric::Interval<Numeric::MPFloat> >() { return "Interval<MPFloat>"; }
    

    int mpfr_hypot(mpfr_t y, const __mpfr_struct* x1, const __mpfr_struct* x2, mpfr_rnd_t r);

    typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
    typedef int mpfr_bin_func(mpfr_t, const __mpfr_struct*, const __mpfr_struct*, mp_rnd_t);
    inline void invoke_mpfr(MPFloat& y, const MPFloat& x, mpfr_func f, mp_rnd_t r) {
      mpfr_t xx;
      mpfr_init_set_f(xx, x.get_mpf_t(), r);
      f(xx, xx, r);
      mpfr_get_f(y.get_mpf_t(), xx, r);
      mpfr_clear(xx);
    }     
    
    inline 
    void invoke_mpfr(MPFloat& y, const MPFloat& x1, const MPFloat& x2, mpfr_bin_func f, mp_rnd_t r) {
      mpfr_t xx1;
      mpfr_t xx2;
      mpfr_init_set_f(xx1, x1.get_mpf_t(), r);
      mpfr_init_set_f(xx2, x2.get_mpf_t(), r);
      f(xx1, xx1, xx2, r);
      mpfr_get_f(y.get_mpf_t(), xx1, r);
      mpfr_clear(xx1);
      mpfr_clear(xx2);
    }

    inline 
    MPFloat invoke_mpfr(const MPFloat& x, mpfr_func f, mp_rnd_t r) {
      MPFloat y; invoke_mpfr(y,x,f,r); return y; 
    }
    
    inline 
    MPFloat invoke_mpfr(const MPFloat& x1, const MPFloat& x2, mpfr_bin_func f, mpfr_rnd_t r) 
    {
      MPFloat y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
    /*! \brief Conversion to a double. */
    template<> inline double conv_approx(const MPFloat& x) { return x.get_d(); }
 
    /*! \brief Conversion from a double. */
    template<> inline MPFloat conv_exact(const double& x) { return MPFloat(x); }
    template<> inline MPFloat conv_down(const double& x) { return conv_exact<MPFloat>(x); }
    template<> inline MPFloat conv_up(const double& x) { return conv_exact<MPFloat>(x); }
    template<> inline MPFloat conv_approx(const double& x) { return conv_exact<MPFloat>(x); }
  
    /*! \brief Conversion to a rational. */
    template<> inline mpq_class conv_exact(const MPFloat& x) { return x; }
    template<> inline mpq_class conv_down(const MPFloat& x) { return conv_exact<mpq_class>(x); }
    template<> inline mpq_class conv_up(const MPFloat& x) { return conv_exact<mpq_class>(x); }
 
    /*! \brief Conversion from a rational. */
    template<> inline MPFloat conv_down(const mpq_class& x) { return MPFloat(x); }
    template<> inline MPFloat conv_up(const mpq_class& x) { return MPFloat(x); }
    template<> inline MPFloat conv_approx(const mpq_class& x) { return MPFloat(x); }
  
    /*! \brief Unary negation. */
    template<> inline MPFloat neg_exact(const MPFloat& x) { return invoke_mpfr(x,mpfr_neg,GMP_RNDN); }
    inline MPFloat neg_down(const MPFloat& x) { return neg_exact(x); }
    template<> inline MPFloat neg_up(const MPFloat& x) { return neg_exact(x); }
   
    /*! \brief Absolute value. */
    template<> inline MPFloat abs_exact(const MPFloat& x) { return invoke_mpfr(x,mpfr_abs,GMP_RNDN); }
    template<> inline MPFloat abs_down(const MPFloat& x) { return abs_exact(x);  }
    template<> inline MPFloat abs_up(const MPFloat& x) { return abs_exact(x);  }
   
    /*! \brief Addition. */
    template<> inline MPFloat add_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDD); }
    template<> inline MPFloat add_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDU); }
    template<> inline MPFloat add_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDN); }
    
    /*! \brief Subtraction. */
    template<> inline MPFloat sub_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDD); }
    template<> inline MPFloat sub_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDU); }
    template<> inline MPFloat sub_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDN); }
    
     /*! \brief Multiplication. */
    template<> inline MPFloat mul_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDD); }
    template<> inline MPFloat mul_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDU); }
    template<> inline MPFloat mul_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDN); }
    
    /*! \brief Division. */
    template<> inline MPFloat div_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDD); }
    template<> inline MPFloat div_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDU); }
    template<> inline MPFloat div_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDN); }
    
    /*! \brief Median. */
    template<> inline MPFloat med_approx(const MPFloat& x1, const MPFloat& x2) {
      return div_approx(add_approx(x1,x2),MPFloat(2)); }
    
    /*! \brief Square root. */
    template<> inline MPFloat sqrt_down(const MPFloat& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDD); }
    template<> inline MPFloat sqrt_up(const MPFloat& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDU); }
    template<> inline MPFloat sqrt_approx(const MPFloat& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDN); }

    /*! \brief Hypoteneuse. */
    template<> inline MPFloat hypot_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDD); }
    template<> inline MPFloat hypot_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDU); }
    template<> inline MPFloat hypot_approx(const MPFloat& x1,const MPFloat& x2) { 
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDN); }
    
    /*! \brief Integer bounds. */
    template<> inline Integer int_down(const MPFloat& x);
    template<> inline Integer int_up(const MPFloat& x);

    /*! \brief Exponential. */
    template<> inline MPFloat exp_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDD); }
    template<> inline MPFloat exp_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDU); };

    /*! \brief Natural logarithm. */
    template<> inline MPFloat log_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDD); }
    template<> inline MPFloat log_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDU); };

    /*! \brief Sine function. */
    template<> inline MPFloat sin_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDD); }
    template<> inline MPFloat sin_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDU); };

    /*! \brief Cosine function. */
    template<> inline MPFloat cos_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDD); }
    template<> inline MPFloat cos_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDU); };

    /*! \brief Tangent function. */
    template<> inline MPFloat tan_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDD); }
    template<> inline MPFloat tan_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDU); };

    /*! \brief Hyperbolic sine function. */
    template<> inline MPFloat sinh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDD); }
    template<> inline MPFloat sinh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDU); };

    /*! \brief Hyperbolic cosine function. */
    template<> inline MPFloat cosh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDD); }
    template<> inline MPFloat cosh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDU); };

    /*! \brief Hyperbolic tangent function. */
    template<> inline MPFloat tanh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDD); }
    template<> inline MPFloat tanh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDU); };

    /*! \brief Inverse sine function. */
    template<> inline MPFloat asin_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDD); }
    template<> inline MPFloat asin_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDU); };

    /*! \brief Inverse cosine function. */
    template<> inline MPFloat acos_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDD); }
    template<> inline MPFloat acos_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDU); };

    /*! \brief Inverse tangent function. */
    template<> inline MPFloat atan_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDD); }
    template<> inline MPFloat atan_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDU); };

    /*! \brief Inverse hyperbolic sine function. */
    template<> inline MPFloat asinh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDD); }
    template<> inline MPFloat asinh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDU); };

    /*! \brief Inverse hyperbolic cosine function. */
    template<> inline MPFloat acosh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDD); }
    template<> inline MPFloat acosh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDU); };

    /*! \brief Inverse hyperbolic tangent function. */
    template<> inline MPFloat atanh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDD); }
    template<> inline MPFloat atanh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDU); };
    

  }
}

#endif /* _ARIADNE_MPFLOAT_H */
