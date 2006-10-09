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
    /*!\ingroup Numeric
     * \brief A multiple-precision floating-point type.
     *  
     * Currently implemented using mpf_class from the GNU Multiple Precision library.
     */
    class MPFloat { };
#else
    //typedef mpf_class MPFloat;
#endif


    class MPFloat {
     private:
      mpfr_t _value;
      typedef __mpfr_struct* pointer;
      typedef const __mpfr_struct* const_pointer;
     public:
      MPFloat() { mpfr_init_set_si(_value,0,GMP_RNDN); }
      MPFloat(int x) { mpfr_init_set_si(_value,x,GMP_RNDN); }
      MPFloat(long int x) { mpfr_init_set_si(_value,x,GMP_RNDN); }
      MPFloat(unsigned int x) { mpfr_init_set_ui(_value,x,GMP_RNDN); }
      MPFloat(long unsigned int x) { mpfr_init_set_ui(_value,x,GMP_RNDN); }
      MPFloat(double x) { mpfr_init_set_d(_value,x,GMP_RNDN); }
      MPFloat(const MPFloat& x) { mpfr_init_set(_value,x._value,GMP_RNDN); }
      MPFloat(const mpf_class& x) { mpfr_init_set_f(_value,x.get_mpf_t(),GMP_RNDN); }
      MPFloat(const mpf_t x) { mpfr_init_set_f(_value,x,GMP_RNDN); }
      MPFloat& operator=(const MPFloat& x) { if(this != &x) { mpfr_set(_value,x._value,GMP_RNDN); } return *this; }
      mpfr_srcptr get_mpfr_t() const { return _value; }
      mpfr_ptr get_mpfr_t() { return _value; }
      mpf_class get_mpf_class() const { mpf_class result; mpfr_get_f(result.get_mpf_t(),_value,GMP_RNDD); return result; }
      operator mpf_class () const { mpf_class result; mpfr_get_f(result.get_mpf_t(),_value,GMP_RNDD); return result; }
      operator const_pointer () const { return _value; }
      operator pointer () { return _value; }
      operator Rational () const { 
        mpf_class f; 
        mpfr_get_f(f.get_mpf_t(),this->get_mpfr_t(),GMP_RNDN); 
        return mpq_class(f);
      }
      
      int precision() const { return mpfr_get_prec(this->get_mpfr_t()); }
    };
    
    inline std::ostream& operator<<(std::ostream& os, const MPFloat& x) {
      return os<<x.get_mpf_class().get_d(); }
    inline std::istream& operator>>(std::istream& is, MPFloat& x) {
      //FIXME: Improve this!
      double dx; is >> dx; x=MPFloat(dx); return is; }

    
    inline int precision(const MPFloat& num) {
      return mpfr_get_prec(num.get_mpfr_t());
    }
  
    template<> class numerical_traits<MPFloat> {
     public:
      typedef ring_tag algebraic_category;
      typedef Rational field_extension_type;
      typedef Interval<MPFloat> arithmetic_type;
    };
      
    template<> inline std::string name<Numeric::MPFloat>() { return "MPFloat"; }
    template<> inline std::string name<Numeric::Interval<Numeric::MPFloat> >() { return "Interval<MPFloat>"; }
    

    int mpfr_hypot(mpfr_t y, const __mpfr_struct* x1, const __mpfr_struct* x2, mpfr_rnd_t r);

    typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
    typedef int mpfr_bin_func(mpfr_t, const __mpfr_struct*, const __mpfr_struct*, mp_rnd_t);
    typedef int mpfr_bin_int_func(mpfr_t, const __mpfr_struct*, long int, mp_rnd_t);
    typedef int mpfr_bin_uint_func(mpfr_t, const __mpfr_struct*, unsigned long int, mp_rnd_t);


    inline void invoke_mpfr(MPFloat& y, const MPFloat& x, mpfr_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x.get_mpfr_t(), r);
    }     
    
    inline 
    void invoke_mpfr(MPFloat& y, const MPFloat& x1, const MPFloat& x2, mpfr_bin_func f, mp_rnd_t r) {
      //std::cerr << "invoke_mpfr(MPFloat& y, const MPFloat& x1, const MPFloat& x2, mpfr_bin_func f, mp_rnd_t r)\n";
      f(y.get_mpfr_t(), x1.get_mpfr_t(), x2.get_mpfr_t(), r);
    }

    inline 
    void invoke_mpfr(MPFloat& y, const MPFloat& x1, const long int& x2, mpfr_bin_int_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x1.get_mpfr_t(), x2, r);
    }

    inline 
    void invoke_mpfr(MPFloat& y, const MPFloat& x1, const unsigned long int& x2, mpfr_bin_uint_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x1.get_mpfr_t(), x2, r);
    }

    
    inline 
    MPFloat invoke_mpfr(const MPFloat& x, mpfr_func f, mp_rnd_t r) {
      MPFloat y; invoke_mpfr(y,x,f,r); return y; 
    }
    
    inline 
    MPFloat invoke_mpfr(const MPFloat& x1, const MPFloat& x2, mpfr_bin_func f, mpfr_rnd_t r) 
    {
      //std::cerr << "MPFloat invoke_mpfr(const MPFloat& x1, const MPFloat& x2, mpfr_bin_func f, mpfr_rnd_t r)\n";
      MPFloat y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
    inline 
    MPFloat invoke_mpfr(const MPFloat& x1, const long int& x2, mpfr_bin_int_func f, mpfr_rnd_t r) 
    {
      MPFloat y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
    inline 
    MPFloat invoke_mpfr(const MPFloat& x1, const long int& x2, mpfr_bin_uint_func f, mpfr_rnd_t r) 
    {
      MPFloat y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
  
    inline bool operator==(const MPFloat& x1, const MPFloat& x2) {
      return mpfr_equal_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator!=(const MPFloat& x1, const MPFloat& x2) {
      return mpfr_lessgreater_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator<=(const MPFloat& x1, const MPFloat& x2) {
      return mpfr_lessequal_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator>=(const MPFloat& x1, const MPFloat& x2) {
      return mpfr_greaterequal_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator< (const MPFloat& x1, const MPFloat& x2) {
      return mpfr_less_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator> (const MPFloat& x1, const MPFloat& x2) {
      return mpfr_greater_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
     
    inline bool operator==(const MPFloat& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)==0; }
    inline bool operator!=(const MPFloat& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)!=0; }
    inline bool operator<=(const MPFloat& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)<=0; }
    inline bool operator>=(const MPFloat& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)>=0; }
     inline bool operator< (const MPFloat& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)< 0; }
    inline bool operator> (const MPFloat& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)> 0; }
  
    inline bool operator==(const MPFloat& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)==0; }
    inline bool operator!=(const MPFloat& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)!=0; }
    inline bool operator<=(const MPFloat& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)<=0; }
    inline bool operator>=(const MPFloat& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)>=0; }
     inline bool operator< (const MPFloat& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)< 0; }
    inline bool operator> (const MPFloat& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)> 0; }
  
      
      
    template<> inline MPFloat min(const MPFloat& x1, const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_min,GMP_RNDN); }
    template<> inline MPFloat max(const MPFloat& x1, const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_max,GMP_RNDN); }
    template<> inline MPFloat neg(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_neg,GMP_RNDN); }
    template<> inline MPFloat abs(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_abs,GMP_RNDN); }




    template<> inline MPFloat conv_exact(const MPFloat& x) { return x; }
    template<> inline MPFloat conv_approx(const MPFloat& x) { return conv_exact<MPFloat>(x); }
    template<> inline MPFloat conv_down(const MPFloat& x) { return conv_exact<MPFloat>(x); }
    template<> inline MPFloat conv_up(const MPFloat& x) { return conv_exact<MPFloat>(x); }
 

    template<> inline double conv_approx(const MPFloat& x) { return mpfr_get_d(x.get_mpfr_t(),GMP_RNDN); }
 
    template<> inline MPFloat conv_exact(const int& n) { return MPFloat(n); }
    template<> inline MPFloat conv_down(const int& n) { return conv_exact<MPFloat>(n); }
    template<> inline MPFloat conv_up(const int& n) { return conv_exact<MPFloat>(n); }
    template<> inline MPFloat conv_approx(const int& n) { return conv_exact<MPFloat>(n); }

    template<> inline MPFloat conv_exact(const double& x) { return MPFloat(x); }
    template<> inline MPFloat conv_down(const double& x) { return conv_exact<MPFloat>(x); }
    template<> inline MPFloat conv_up(const double& x) { return conv_exact<MPFloat>(x); }
    template<> inline MPFloat conv_approx(const double& x) { return conv_exact<MPFloat>(x); }
  
    template<> inline int conv_down(const MPFloat& x) { return mpfr_get_si(x.get_mpfr_t(),GMP_RNDD); }
    template<> inline int conv_up(const MPFloat& x) { return mpfr_get_si(x.get_mpfr_t(),GMP_RNDU); }
 
    template<> inline Rational conv_exact(const MPFloat& x) { return Rational(x); }
    template<> inline Rational conv_approx(const MPFloat& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_down(const MPFloat& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_up(const MPFloat& x) { return conv_exact<Rational>(x); }
 

    template<> inline MPFloat min_exact(const MPFloat& x1, const MPFloat& x2) {
      return (x1<=x2) ? x1 : x2; }
    template<> inline MPFloat min_approx(const MPFloat& x1, const MPFloat& x2) {
      return min_exact(x1,x2); }
    template<> inline MPFloat min_down(const MPFloat& x1, const MPFloat& x2) {
      return min_exact(x1,x2); }
    template<> inline MPFloat min_up(const MPFloat& x1, const MPFloat& x2) {
      return min_exact(x1,x2); }
  
    template<> inline MPFloat max_exact(const MPFloat& x1, const MPFloat& x2) {
      return (x1>=x2) ? x1 : x2; }
    template<> inline MPFloat max_approx(const MPFloat& x1, const MPFloat& x2) {
      return max_exact(x1,x2); }
    template<> inline MPFloat max_down(const MPFloat& x1, const MPFloat& x2) {
      return max_exact(x1,x2); }
    template<> inline MPFloat max_up(const MPFloat& x1, const MPFloat& x2) {
      return max_exact(x1,x2); }
  
  
  
    template<> inline MPFloat neg_exact(const MPFloat& x) { return invoke_mpfr(x,mpfr_neg,GMP_RNDN); }
    template<> inline MPFloat neg_approx(const MPFloat& x) { return neg_exact(x); }
    template<> inline MPFloat neg_down(const MPFloat& x) { return neg_exact(x); }
    template<> inline MPFloat neg_up(const MPFloat& x) { return neg_exact(x); }


    template<> inline MPFloat abs_exact(const MPFloat& x) { return invoke_mpfr(x,mpfr_abs,GMP_RNDN); }
    template<> inline MPFloat abs_approx(const MPFloat& x) { return abs_exact(x);  }
    template<> inline MPFloat abs_down(const MPFloat& x) { return abs_exact(x);  }
    template<> inline MPFloat abs_up(const MPFloat& x) { return abs_exact(x);  }

    template<> inline MPFloat add_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDN); }
    template<> inline MPFloat add_down(const MPFloat& x1,const MPFloat& x2) {
      //std::cerr << "add_down(MPFloat,MPFloat)" << std::endl;
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDD); }
    template<> inline MPFloat add_up(const MPFloat& x1,const MPFloat& x2) {
      //std::cerr << "add_up(MPFloat,MPFloat)" << std::endl;
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDU); }
    
    template<> inline MPFloat sub_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDN); }
    template<> inline MPFloat sub_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDD); }
    template<> inline MPFloat sub_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDU); }
    
    template<> inline MPFloat mul_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDN); }
    template<> inline MPFloat mul_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDD); }
    template<> inline MPFloat mul_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDU); }
    
    template<> inline MPFloat div_approx(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDN); }
    template<> inline MPFloat div_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDD); }
    template<> inline MPFloat div_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDU); }
    
    template<> inline MPFloat pow_approx(const MPFloat& x1, const unsigned int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_ui,GMP_RNDN); }
    template<> inline MPFloat pow_down(const MPFloat& x1,const unsigned int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_ui,GMP_RNDD); }
    template<> inline MPFloat pow_up(const MPFloat& x1,const unsigned int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_ui,GMP_RNDU); }
    
    template<> inline MPFloat pow_approx(const MPFloat& x1, const int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_si,GMP_RNDN); }
    template<> inline MPFloat pow_down(const MPFloat& x1,const int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_si,GMP_RNDD); }
    template<> inline MPFloat pow_up(const MPFloat& x1,const int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_si,GMP_RNDU); }
    
    template<> inline MPFloat med_approx(const MPFloat& x1, const MPFloat& x2) {
      return div_approx(add_approx(x1,x2),MPFloat(2)); }
    
    inline MPFloat mul_approx(const int& n, const MPFloat& x) {
      return invoke_mpfr(x,n,mpfr_mul_si,GMP_RNDN); }
    inline MPFloat mul_approx(const MPFloat& x, const int& n) {
      return invoke_mpfr(x,n,mpfr_mul_si,GMP_RNDN); }
    inline MPFloat div_approx(const MPFloat& x, const int& n) {
      return invoke_mpfr(x,n,mpfr_div_si,GMP_RNDN); }
      
    inline MPFloat mul_approx(const uint& n, const MPFloat& x) {
      return invoke_mpfr(x,n,mpfr_mul_ui,GMP_RNDN); }
    inline MPFloat mul_approx(const MPFloat& x, const uint& n) {
      return invoke_mpfr(x,n,mpfr_mul_ui,GMP_RNDN); }
    inline MPFloat div_approx(const MPFloat& x, const uint& n) {
      return invoke_mpfr(x,n,mpfr_div_ui,GMP_RNDN); }
     
    inline MPFloat mul_approx(const double& d, const MPFloat& x) {
      return mul_approx(x,MPFloat(d)); }
    inline MPFloat mul_approx(const MPFloat& x, const double& d) {
      return mul_approx(x,MPFloat(d)); }
      
      
    template<> inline MPFloat sqrt_approx(const MPFloat& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDN); }
    template<> inline MPFloat sqrt_down(const MPFloat& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDD); }
    template<> inline MPFloat sqrt_up(const MPFloat& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDU); }

    template<> inline MPFloat hypot_approx(const MPFloat& x1,const MPFloat& x2) { 
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDN); }
    template<> inline MPFloat hypot_down(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDD); }
    template<> inline MPFloat hypot_up(const MPFloat& x1,const MPFloat& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDU); }
    
    template<> inline Integer int_down(const MPFloat& x);
    template<> inline Integer int_up(const MPFloat& x);

    template<> inline MPFloat exp_approx(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDN); }
    template<> inline MPFloat exp_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDD); }
    template<> inline MPFloat exp_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDU); };

    template<> inline MPFloat log_approx(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDN); }
    template<> inline MPFloat log_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDD); }
    template<> inline MPFloat log_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDU); };

    template<> inline MPFloat sin_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDD); }
    template<> inline MPFloat sin_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDU); };

    template<> inline MPFloat cos_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDD); }
    template<> inline MPFloat cos_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDU); };

    template<> inline MPFloat tan_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDD); }
    template<> inline MPFloat tan_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDU); };

    template<> inline MPFloat sinh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDD); }
    template<> inline MPFloat sinh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDU); };

    template<> inline MPFloat cosh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDD); }
    template<> inline MPFloat cosh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDU); };

    template<> inline MPFloat tanh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDD); }
    template<> inline MPFloat tanh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDU); };

    template<> inline MPFloat asin_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDD); }
    template<> inline MPFloat asin_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDU); };

    template<> inline MPFloat acos_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDD); }
    template<> inline MPFloat acos_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDU); };

    template<> inline MPFloat atan_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDD); }
    template<> inline MPFloat atan_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDU); };

    template<> inline MPFloat asinh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDD); }
    template<> inline MPFloat asinh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDU); };

    template<> inline MPFloat acosh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDD); }
    template<> inline MPFloat acosh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDU); };

    template<> inline MPFloat atanh_down(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDD); }
    template<> inline MPFloat atanh_up(const MPFloat& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDU); };
    
    inline MPFloat operator-(const MPFloat& x) {
      return neg(x); }

    Interval<MPFloat> operator+(const MPFloat&, const MPFloat&);
    Interval<MPFloat> operator-(const MPFloat&, const MPFloat&);
    Interval<MPFloat> operator*(const MPFloat&, const MPFloat&);
    Interval<MPFloat> operator/(const MPFloat&, const MPFloat&);
      
  }
}

#endif /* _ARIADNE_MPFLOAT_H */
