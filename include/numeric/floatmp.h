/***************************************************************************
 *            floatmp.h
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
 
/*! \file numeric/floatmp.h
 *  \brief Type definitions and conversion operators for multiple-precision floating-point numbers.
 */

#ifndef ARIADNE_FLOATMP_H
#define ARIADNE_FLOATMP_H

#include <mpfr.h>
// Don't use mpf2mpfr, since it ignores rounding modes
// #include <mpf2mpfr.h>

#include <gmpxx.h>

#include "../numeric/numerical_traits.h"
#include "../numeric/function.h"
#include "../numeric/rational.h"

namespace Ariadne {
  namespace Numeric {

    /*!\ingroup Numeric
     * \brief A multiple-precision floating-point type.
     *  
     * Currently implemented using mpf_class from the GNU Multiple Precision library.
     */
    class FloatMP {
     private:
      mpfr_t _value;
      typedef __mpfr_struct* pointer;
      typedef const __mpfr_struct* const_pointer;
     public:
      FloatMP() { mpfr_init_set_si(_value,0,GMP_RNDN); }
      FloatMP(int x) { mpfr_init_set_si(_value,x,GMP_RNDN); }
      FloatMP(long int x) { mpfr_init_set_si(_value,x,GMP_RNDN); }
      FloatMP(unsigned int x) { mpfr_init_set_ui(_value,x,GMP_RNDN); }
      FloatMP(long unsigned int x) { mpfr_init_set_ui(_value,x,GMP_RNDN); }
      FloatMP(double x) { mpfr_init_set_d(_value,x,GMP_RNDN); }
      FloatMP(const FloatMP& x) { mpfr_init_set(_value,x._value,GMP_RNDN); }
      FloatMP(const mpf_class& x) { mpfr_init_set_f(_value,x.get_mpf_t(),GMP_RNDN); }
      FloatMP(const mpf_t x) { mpfr_init_set_f(_value,x,GMP_RNDN); }
      FloatMP& operator=(const FloatMP& x) { if(this != &x) { mpfr_set(_value,x._value,GMP_RNDN); } return *this; }
      // FIXME: Use mpfr_class instead
      mpf_class get_base() const { mpf_class result; mpfr_get_f(result.get_mpf_t(),this->get_mpfr_t(),GMP_RNDN); return result; }
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
      static uint default_precision() { return mpfr_get_default_prec(); }
      static void set_default_precision(uint p) { mpfr_set_default_prec(p); }
      int precision() const { return mpfr_get_prec(this->get_mpfr_t()); }
    };
    
    inline std::ostream& operator<<(std::ostream& os, const FloatMP& x) {
      return os<<x.get_mpf_class().get_d(); }
    inline std::istream& operator>>(std::istream& is, FloatMP& x) {
      //FIXME: Improve this!
      double dx; is >> dx; x=FloatMP(dx); return is; }

    
    inline int precision(const FloatMP& num) {
      return mpfr_get_prec(num.get_mpfr_t());
    }
  
    
    template<> inline std::string name<Numeric::FloatMP>() { return "FloatMP"; }
    template<> inline std::string name<Numeric::Interval<Numeric::FloatMP> >() { return "Interval<FloatMP>"; }
    

    int mpfr_hypot(mpfr_t y, const __mpfr_struct* x1, const __mpfr_struct* x2, mpfr_rnd_t r);

    typedef int mpfr_const(mpfr_t, mp_rnd_t);
    typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
    typedef int mpfr_func_rnd(mpfr_t, const __mpfr_struct*);
    typedef int mpfr_bin_func(mpfr_t, const __mpfr_struct*, const __mpfr_struct*, mp_rnd_t);
    typedef int mpfr_bin_int_func(mpfr_t, const __mpfr_struct*, long int, mp_rnd_t);
    typedef int mpfr_bin_uint_func(mpfr_t, const __mpfr_struct*, unsigned long int, mp_rnd_t);


    inline void invoke_mpfr(FloatMP& y, mpfr_const f, mp_rnd_t r) {
      f(y.get_mpfr_t(), r);
    }     
    
    inline void invoke_mpfr(FloatMP& y, const FloatMP& x, mpfr_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x.get_mpfr_t(), r);
    }     
    
    inline 
    void invoke_mpfr(FloatMP& y, const FloatMP& x1, const FloatMP& x2, mpfr_bin_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x1.get_mpfr_t(), x2.get_mpfr_t(), r);
    }

    inline 
    void invoke_mpfr(FloatMP& y, const FloatMP& x1, const long int& x2, mpfr_bin_int_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x1.get_mpfr_t(), x2, r);
    }

    inline 
    void invoke_mpfr(FloatMP& y, const FloatMP& x1, const unsigned long int& x2, mpfr_bin_uint_func f, mp_rnd_t r) {
      f(y.get_mpfr_t(), x1.get_mpfr_t(), x2, r);
    }

    inline 
    void invoke_mpfr(FloatMP& y, const FloatMP& x1, mpfr_func_rnd f) {
      f(y.get_mpfr_t(), x1.get_mpfr_t());
    }

    
    inline 
    FloatMP invoke_mpfr(mpfr_const f, mp_rnd_t r) {
      FloatMP y; invoke_mpfr(y,f,r); return y; 
    }
    
    inline 
    FloatMP invoke_mpfr(const FloatMP& x, mpfr_func f, mp_rnd_t r) {
      FloatMP y; invoke_mpfr(y,x,f,r); return y; 
    }
    
    inline 
    FloatMP invoke_mpfr(const FloatMP& x1, const FloatMP& x2, mpfr_bin_func f, mpfr_rnd_t r) 
    {
      FloatMP y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
    inline 
    FloatMP invoke_mpfr(const FloatMP& x1, const long int& x2, mpfr_bin_int_func f, mpfr_rnd_t r) 
    {
      FloatMP y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
    inline 
    FloatMP invoke_mpfr(const FloatMP& x1, const long int& x2, mpfr_bin_uint_func f, mpfr_rnd_t r) 
    {
      FloatMP y; invoke_mpfr(y,x1,x2,f,r); return y; 
    }
    
    inline 
    FloatMP invoke_mpfr(const FloatMP& x, mpfr_func_rnd f) 
    {
      FloatMP y; invoke_mpfr(y,x,f); return y; 
    }
    
  
    inline bool operator==(const FloatMP& x1, const FloatMP& x2) {
      return mpfr_equal_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator!=(const FloatMP& x1, const FloatMP& x2) {
      return mpfr_lessgreater_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator<=(const FloatMP& x1, const FloatMP& x2) {
      return mpfr_lessequal_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator>=(const FloatMP& x1, const FloatMP& x2) {
      return mpfr_greaterequal_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator< (const FloatMP& x1, const FloatMP& x2) {
      return mpfr_less_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
    inline bool operator> (const FloatMP& x1, const FloatMP& x2) {
      return mpfr_greater_p(x1.get_mpfr_t(),x2.get_mpfr_t()); }
     
    inline bool operator==(const FloatMP& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)==0; }
    inline bool operator!=(const FloatMP& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)!=0; }
    inline bool operator<=(const FloatMP& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)<=0; }
    inline bool operator>=(const FloatMP& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)>=0; }
     inline bool operator< (const FloatMP& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)< 0; }
    inline bool operator> (const FloatMP& x1, const int& x2) {
      return mpfr_cmp_si(x1.get_mpfr_t(),x2)> 0; }
  
    inline bool operator==(const FloatMP& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)==0; }
    inline bool operator!=(const FloatMP& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)!=0; }
    inline bool operator<=(const FloatMP& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)<=0; }
    inline bool operator>=(const FloatMP& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)>=0; }
     inline bool operator< (const FloatMP& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)< 0; }
    inline bool operator> (const FloatMP& x1, const double& x2) {
      return mpfr_cmp_d(x1.get_mpfr_t(),x2)> 0; }
  
    inline bool operator==(const FloatMP& x1, const mpq_class& x2) {
      return mpfr_cmp_q(x1.get_mpfr_t(),x2.get_mpq_t())==0; }
    inline bool operator!=(const FloatMP& x1, const mpq_class& x2) {
      return mpfr_cmp_q(x1.get_mpfr_t(),x2.get_mpq_t())!=0; }
    inline bool operator<=(const FloatMP& x1, const mpq_class& x2) {
      return mpfr_cmp_q(x1.get_mpfr_t(),x2.get_mpq_t())<=0; }
    inline bool operator>=(const FloatMP& x1, const mpq_class& x2) {
      return mpfr_cmp_q(x1.get_mpfr_t(),x2.get_mpq_t())>=0; }
     inline bool operator< (const FloatMP& x1, const mpq_class& x2) {
      return mpfr_cmp_q(x1.get_mpfr_t(),x2.get_mpq_t())< 0; }
    inline bool operator> (const FloatMP& x1, const mpq_class& x2) {
      return mpfr_cmp_q(x1.get_mpfr_t(),x2.get_mpq_t())> 0; }
  
  
  template<> inline FloatMP next_up(const FloatMP& x) { FloatMP y(x); mpfr_nextabove(y.get_mpfr_t()); return y; }
  template<> inline FloatMP next_down(const FloatMP& x) { FloatMP y(x); mpfr_nextbelow(y.get_mpfr_t()); return y; }
      
      
    template<> inline FloatMP min(const FloatMP& x1, const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_min,GMP_RNDN); }
    template<> inline FloatMP max(const FloatMP& x1, const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_max,GMP_RNDN); }
    template<> inline FloatMP neg(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_neg,GMP_RNDN); }
    template<> inline FloatMP abs(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_abs,GMP_RNDN); }



    template<> inline FloatMP conv_exact(const FloatMP& x) { return x; }
    template<> inline FloatMP conv_approx(const FloatMP& x) { return conv_exact<FloatMP>(x); }
    template<> inline FloatMP conv_down(const FloatMP& x) { return conv_exact<FloatMP>(x); }
    template<> inline FloatMP conv_up(const FloatMP& x) { return conv_exact<FloatMP>(x); }
 

    template<> inline double conv_approx(const FloatMP& x) { return mpfr_get_d(x.get_mpfr_t(),GMP_RNDN); }
 
    template<> inline FloatMP conv_exact(const int& n) { return FloatMP(n); }
    template<> inline FloatMP conv_down(const int& n) { return conv_exact<FloatMP>(n); }
    template<> inline FloatMP conv_up(const int& n) { return conv_exact<FloatMP>(n); }
    template<> inline FloatMP conv_approx(const int& n) { return conv_exact<FloatMP>(n); }

    template<> inline FloatMP conv_exact(const double& x) { return FloatMP(x); }
    template<> inline FloatMP conv_down(const double& x) { return conv_exact<FloatMP>(x); }
    template<> inline FloatMP conv_up(const double& x) { return conv_exact<FloatMP>(x); }
    template<> inline FloatMP conv_approx(const double& x) { return conv_exact<FloatMP>(x); }
   
    template<> inline FloatMP conv_down(const Rational& x) { 
      FloatMP r; mpfr_init_set_q(r.get_mpfr_t(),x.get_mpq_t(),GMP_RNDD); return r; }
    template<> inline FloatMP conv_up(const Rational& x) {
      FloatMP r; mpfr_init_set_q(r.get_mpfr_t(),x.get_mpq_t(),GMP_RNDU); return r; }
    template<> inline FloatMP conv_approx(const Rational& x) {
      FloatMP r; mpfr_init_set_q(r.get_mpfr_t(),x.get_mpq_t(),GMP_RNDD); return r; }
   
    template<> inline mpf_class conv_exact(const FloatMP& x) { 
      mpf_class r; mpfr_get_f(r.get_mpf_t(),x.get_mpfr_t(),GMP_RNDN); return r; }
    template<> inline mpf_class conv_approx(const FloatMP& x) { 
      mpf_class r; mpfr_get_f(r.get_mpf_t(),x.get_mpfr_t(),GMP_RNDN); return r; }
    template<> inline mpf_class conv_down(const FloatMP& x) { 
      mpf_class r; mpfr_get_f(r.get_mpf_t(),x.get_mpfr_t(),GMP_RNDD); return r; }
    template<> inline mpf_class conv_up(const FloatMP& x) { 
      mpf_class r; mpfr_get_f(r.get_mpf_t(),x.get_mpfr_t(),GMP_RNDU); return r; }
  
    template<> inline Rational conv_exact(const FloatMP& x) { return Rational(x); }
    template<> inline Rational conv_approx(const FloatMP& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_down(const FloatMP& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_up(const FloatMP& x) { return conv_exact<Rational>(x); }
 
    template<> inline Integer int_down(const FloatMP& x) { 
      mpz_class z; mpfr_get_z(z.get_mpz_t(),x.get_mpfr_t(),GMP_RNDD); return z; }
    template<> inline Integer int_up(const FloatMP& x){ 
      mpz_class z; mpfr_get_z(z.get_mpz_t(),x.get_mpfr_t(),GMP_RNDD); return z; }

    template<> inline int int_down(const FloatMP& x) { return mpfr_get_si(x.get_mpfr_t(),GMP_RNDD); }
    template<> inline int int_up(const FloatMP& x){ return mpfr_get_si(x.get_mpfr_t(),GMP_RNDU); }


    template<> inline FloatMP min_exact(const FloatMP& x1, const FloatMP& x2) {
      return (x1<=x2) ? x1 : x2; }
    template<> inline FloatMP min_approx(const FloatMP& x1, const FloatMP& x2) {
      return min_exact(x1,x2); }
    template<> inline FloatMP min_down(const FloatMP& x1, const FloatMP& x2) {
      return min_exact(x1,x2); }
    template<> inline FloatMP min_up(const FloatMP& x1, const FloatMP& x2) {
      return min_exact(x1,x2); }
  
    template<> inline FloatMP max_exact(const FloatMP& x1, const FloatMP& x2) {
      return (x1>=x2) ? x1 : x2; }
    template<> inline FloatMP max_approx(const FloatMP& x1, const FloatMP& x2) {
      return max_exact(x1,x2); }
    template<> inline FloatMP max_down(const FloatMP& x1, const FloatMP& x2) {
      return max_exact(x1,x2); }
    template<> inline FloatMP max_up(const FloatMP& x1, const FloatMP& x2) {
      return max_exact(x1,x2); }
  
  
  
    template<> inline FloatMP neg_exact(const FloatMP& x) { return invoke_mpfr(x,mpfr_neg,GMP_RNDN); }
    template<> inline FloatMP neg_approx(const FloatMP& x) { return neg_exact(x); }
    template<> inline FloatMP neg_down(const FloatMP& x) { return neg_exact(x); }
    template<> inline FloatMP neg_up(const FloatMP& x) { return neg_exact(x); }


    template<> inline FloatMP abs_exact(const FloatMP& x) { return invoke_mpfr(x,mpfr_abs,GMP_RNDN); }
    template<> inline FloatMP abs_approx(const FloatMP& x) { return abs_exact(x);  }
    template<> inline FloatMP abs_down(const FloatMP& x) { return abs_exact(x);  }
    template<> inline FloatMP abs_up(const FloatMP& x) { return abs_exact(x);  }

    template<> inline FloatMP add_approx(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDN); }
    template<> inline FloatMP add_down(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDD); }
    template<> inline FloatMP add_up(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDU); }
    
    template<> inline FloatMP sub_approx(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDN); }
    template<> inline FloatMP sub_down(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDD); }
    template<> inline FloatMP sub_up(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDU); }
    
    template<> inline FloatMP mul_approx(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDN); }
    template<> inline FloatMP mul_down(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDD); }
    template<> inline FloatMP mul_up(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDU); }
    
    template<> inline FloatMP mul_approx(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_mul_si,GMP_RNDN); }
    template<> inline FloatMP mul_down(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_mul_si,GMP_RNDD); }
    template<> inline FloatMP mul_up(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_mul_si,GMP_RNDU); }
    
    template<> inline FloatMP div_approx(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDN); }
    template<> inline FloatMP div_down(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDD); }
    template<> inline FloatMP div_up(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDU); }

    template<> inline FloatMP div_approx(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_div_si,GMP_RNDN); }
    template<> inline FloatMP div_down(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_div_si,GMP_RNDD); }
    template<> inline FloatMP div_up(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_div_si,GMP_RNDU); }

    template<> inline FloatMP pow_approx(const FloatMP& x1, const unsigned int& n2) {
      return invoke_mpfr(x1,(unsigned long int)n2,mpfr_pow_ui,GMP_RNDN); }
    template<> inline FloatMP pow_down(const FloatMP& x1,const unsigned int& n2) {
      return invoke_mpfr(x1,(unsigned long int)n2,mpfr_pow_ui,GMP_RNDD); }
    template<> inline FloatMP pow_up(const FloatMP& x1,const unsigned int& n2) {
      return invoke_mpfr(x1,(unsigned long int)n2,mpfr_pow_ui,GMP_RNDU); }
    
    template<> inline FloatMP pow_approx(const FloatMP& x1, const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_pow_si,GMP_RNDN); }
    template<> inline FloatMP pow_down(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_pow_si,GMP_RNDD); }
    template<> inline FloatMP pow_up(const FloatMP& x1,const int& n2) {
      return invoke_mpfr(x1,(long int)n2,mpfr_pow_si,GMP_RNDU); }
    
    template<> inline FloatMP pow_approx(const FloatMP& x1, const unsigned long int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_ui,GMP_RNDN); }
    template<> inline FloatMP pow_down(const FloatMP& x1,const unsigned long int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_ui,GMP_RNDD); }
    template<> inline FloatMP pow_up(const FloatMP& x1,const unsigned long int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_ui,GMP_RNDU); }
    
    template<> inline FloatMP pow_approx(const FloatMP& x1, const long int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_si,GMP_RNDN); }
    template<> inline FloatMP pow_down(const FloatMP& x1,const long int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_si,GMP_RNDD); }
    template<> inline FloatMP pow_up(const FloatMP& x1,const long int& n2) {
      return invoke_mpfr(x1,n2,mpfr_pow_si,GMP_RNDU); }
     
    template<> inline FloatMP med_approx(const FloatMP& x1, const FloatMP& x2) {
      return div_approx(add_approx(x1,x2),FloatMP(2)); }
    
    inline FloatMP mul_approx(const int& n, const FloatMP& x) {
      return invoke_mpfr(x,n,mpfr_mul_si,GMP_RNDN); }
    inline FloatMP mul_approx(const FloatMP& x, const int& n) {
      return invoke_mpfr(x,n,mpfr_mul_si,GMP_RNDN); }
    inline FloatMP div_approx(const FloatMP& x, const int& n) {
      return invoke_mpfr(x,n,mpfr_div_si,GMP_RNDN); }

    inline FloatMP mul_approx(const long int& n, const FloatMP& x) {
      return invoke_mpfr(x,n,mpfr_mul_si,GMP_RNDN); }
    inline FloatMP mul_approx(const FloatMP& x, const long int& n) {
      return invoke_mpfr(x,n,mpfr_mul_si,GMP_RNDN); }
    inline FloatMP div_approx(const FloatMP& x, const long int& n) {
      return invoke_mpfr(x,n,mpfr_div_si,GMP_RNDN); }

    inline FloatMP mul_approx(const uint& n, const FloatMP& x) {
      return invoke_mpfr(x,n,mpfr_mul_ui,GMP_RNDN); }
    inline FloatMP mul_approx(const FloatMP& x, const unsigned int& n) {
      return invoke_mpfr(x,n,mpfr_mul_ui,GMP_RNDN); }
    inline FloatMP div_approx(const FloatMP& x, const unsigned int& n) {
      return invoke_mpfr(x,n,mpfr_div_ui,GMP_RNDN); }
    
    inline FloatMP mul_approx(const long unsigned int& n, const FloatMP& x) {
      return invoke_mpfr(x,n,mpfr_mul_ui,GMP_RNDN); }
    inline FloatMP mul_approx(const FloatMP& x, const long unsigned int& n) {
      return invoke_mpfr(x,n,mpfr_mul_ui,GMP_RNDN); }
    inline FloatMP div_approx(const FloatMP& x, const long unsigned int& n) {
      return invoke_mpfr(x,n,mpfr_div_ui,GMP_RNDN); }
    
    inline FloatMP mul_approx(const double& d, const FloatMP& x) {
      return mul_approx(x,FloatMP(d)); }
    inline FloatMP mul_approx(const FloatMP& x, const double& d) {
      return mul_approx(x,FloatMP(d)); }
      
      
    template<> inline FloatMP sqrt_approx(const FloatMP& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDN); }
    template<> inline FloatMP sqrt_down(const FloatMP& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDD); }
    template<> inline FloatMP sqrt_up(const FloatMP& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDU); }

    template<> inline FloatMP hypot_approx(const FloatMP& x1,const FloatMP& x2) { 
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDN); }
    template<> inline FloatMP hypot_down(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDD); }
    template<> inline FloatMP hypot_up(const FloatMP& x1,const FloatMP& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDU); }
    
    template<> inline FloatMP floor(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_floor); }
    template<> inline FloatMP ceil(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_ceil); }
    
    template<> inline FloatMP exp_approx(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDN); }
    template<> inline FloatMP exp_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDD); }
    template<> inline FloatMP exp_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDU); };

    template<> inline FloatMP log_approx(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDN); }
    template<> inline FloatMP log_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDD); }
    template<> inline FloatMP log_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDU); };

    template<> inline FloatMP pi_approx() {
      return invoke_mpfr(mpfr_const_pi,GMP_RNDN); }
    template<> inline FloatMP pi_down() {
      return invoke_mpfr(mpfr_const_pi,GMP_RNDD); }
    template<> inline FloatMP pi_up() {
      return invoke_mpfr(mpfr_const_pi,GMP_RNDU); };

    template<> inline FloatMP sin_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDD); }
    template<> inline FloatMP sin_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDU); };

    template<> inline FloatMP cos_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDD); }
    template<> inline FloatMP cos_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDU); };

    template<> inline FloatMP tan_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDD); }
    template<> inline FloatMP tan_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDU); };

    template<> inline FloatMP sinh_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDD); }
    template<> inline FloatMP sinh_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDU); };

    template<> inline FloatMP cosh_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDD); }
    template<> inline FloatMP cosh_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDU); };

    template<> inline FloatMP tanh_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDD); }
    template<> inline FloatMP tanh_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDU); };

    template<> inline FloatMP asin_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDD); }
    template<> inline FloatMP asin_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDU); };

    template<> inline FloatMP acos_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDD); }
    template<> inline FloatMP acos_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDU); };

    template<> inline FloatMP atan_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDD); }
    template<> inline FloatMP atan_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDU); };

    template<> inline FloatMP asinh_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDD); }
    template<> inline FloatMP asinh_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDU); };

    template<> inline FloatMP acosh_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDD); }
    template<> inline FloatMP acosh_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDU); };

    template<> inline FloatMP atanh_down(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDD); }
    template<> inline FloatMP atanh_up(const FloatMP& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDU); };
    
    inline FloatMP operator+(const FloatMP& x) { return x; }
    inline FloatMP operator-(const FloatMP& x) { return neg(x); }

    Interval<FloatMP> operator+(const FloatMP&, const FloatMP&);
    Interval<FloatMP> operator-(const FloatMP&, const FloatMP&);
    Interval<FloatMP> operator*(const FloatMP&, const FloatMP&);
    Interval<FloatMP> operator/(const FloatMP&, const FloatMP&);
      
    
    
      
/*  Expression-templated functions for possible future use.
      
    template<class Rnd> inline mpfr_rnd_t mpfr_rounding_mode();
    template<> inline mpfr_rnd_t mpfr_rounding_mode<approx>() { return GMP_RNDN; }
    template<> inline mpfr_rnd_t mpfr_rounding_mode<down>() { return GMP_RNDD; }
    template<> inline mpfr_rnd_t mpfr_rounding_mode<up>() { return GMP_RNDU; }

    class MpfrConstantFunction
      : public ConstantFunction<FloatMP>
    {
     public:
      MpfrConstantFunction(mpfr_const_func ff, const mpfr_rnd_t rr) : f(ff), r(rr) { }
      void operator() (FloatMP& r) const;
      FloatMP operator() () const;
      mpfr_const_func* f; mpfr_rnd_t r;
    };

    class MpfrUnaryFunction
      : public UnaryFunction<FloatMP,FloatMP>
    {
     public:
      MpfrUnaryFunction(mpfr_func ff, const mpfr_rnd_t rr) : f(ff), r(rr) { }
      void operator() (FloatMP& r, const FloatMP& a) const;
      FloatMP operator()(const FloatMP& a) const;
      mpfr_func* f; mpfr_rnd_t r;
    };

    class MpfrBinaryFunction
      : public BinaryFunction<FloatMP,FloatMP,FloatMP>
    {
     public:
      MpfrBinaryFunction(mpfr_bin_func ff, const mpfr_rnd_t rr) : f(ff), r(rr) { }
      void operator() (FloatMP& r, const FloatMP& a1, const FloatMP& a2) const;
      FloatMP operator()(const FloatMP& a1, const FloatMP& a2) const;
      mpfr_bin_func* f; mpfr_rnd_t r;
    };

    template<typename A>
    inline
    UnaryExpression<A,MpfrUnaryFunction> 
    mpfr_unary_expression(const A& a, mpfr_func f, mpfr_rnd_t r)
    {
      return UnaryExpression<A,MpfrUnaryFunction>(a,MpfrUnaryFunction(f,r));
    }

    template<typename A1, typename A2>
    inline
    BinaryExpression<A1,A2,MpfrBinaryFunction> 
    mpfr_binary_expression(const A1& a1, const A2& a2, mpfr_bin_func f, mpfr_rnd_t r)
    {
      return BinaryExpression<A1,A2,MpfrBinaryFunction>(a1,a2,MpfrBinaryFunction(f,r));
    }

    template<typename A1, typename A2, typename Rnd>
    inline
    BinaryExpression<A1,A2,MpfrBinaryFunction> 
    mpfr_binary_expression(const A1& a1, const A2& a2, mpfr_bin_func f, Rnd)
    {
      return BinaryExpression<A1,A2,MpfrBinaryFunction>(a1,a2,MpfrBinaryFunction(f,mpfr_rounding_mode<Rnd>()));
    }

    template<class E1, class E2> 
    inline 
    BinaryExpression<E1,E2,MpfrBinaryFunction> 
    add_approx(const Expression<FloatMP,E1>& x1, const Expression<FloatMP,E2>&x2) {
      return mpfr_binary_expression(x1.promote(),x2.promote(),mpfr_add,GMP_RNDN);
    }

*/
      
      
  }
}

#endif /* ARIADNE_MPFLOAT_H */
