/***************************************************************************
 *            rounding.h
 *
 *  Copyright 2006  Alberto Casagrande, Pieter Collins
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
  
/*! \file rounding.h
 *  \brief Specify rounding modes for Boost interval.
 */

#ifndef _ARIADNE_ROUNDING_H
#define _ARIADNE_ROUNDING_H

#include <boost/numeric/interval.hpp>

#include "../numeric/numerical_traits.h"
#include "../numeric/function.h"

namespace Ariadne {
  namespace Numeric {

    template<typename R> struct rounded_math;
      
    template<typename R>
    struct rounded_arith_exact
    {
      R neg_down(const R& x) { return Numeric::neg_exact(x); } 
      R abs_down(const R& x) { return Numeric::abs_exact(x); } 
      R add_down(const R& x1, const R& x2) { return Numeric::add_exact(x1,x2); } 
      R sub_down(const R& x1, const R& x2) { return Numeric::sub_exact(x1,x2); } 
      R mul_down(const R& x1, const R& x2) { return Numeric::mul_exact(x1,x2); } 
      R div_down(const R& x1, const R& x2) { return Numeric::div_exact(x1,x2); } 
 
      R neg_up(const R& x) { return Numeric::neg_exact(x); } 
      R abs_up(const R& x) { return Numeric::abs_exact(x); } 
      R add_up(const R& x1, const R& x2) { return Numeric::add_exact(x1,x2); } 
      R sub_up(const R& x1, const R& x2) { return Numeric::sub_exact(x1,x2); } 
      R mul_up(const R& x1, const R& x2) { return Numeric::mul_exact(x1,x2); } 
      R div_up(const R& x1, const R& x2) { return Numeric::div_exact(x1,x2); } 
    };


    template<typename R>
    struct rounded_math
    {
     public:
      template<typename X> R conv_down(const X& x) { return Numeric::conv_down<R>(x); }
      template<typename X> R conv_up(const X& x) { return Numeric::conv_up<R>(x); }
      
      R neg_down(const R& x) { return neg_down(x); } 
      R abs_down(const R& x) { return Numeric::abs_down(x); } 
      R add_down(const R& x1, const R& x2) { return Numeric::add_down(x1,x2); } 
      R sub_down(const R& x1, const R& x2) { return Numeric::sub_down(x1,x2); } 
      R mul_down(const R& x1, const R& x2) { return Numeric::mul_down(x1,x2); } 
      R div_down(const R& x1, const R& x2) { return Numeric::div_down(x1,x2); } 
      R sqrt_down(R x) { return Numeric::sqrt_down(x); } 
      R hypot_down(const R& x1, const R& x2) { return Numeric::hypot_down(x1,x2); } 
      
      R exp_down(const R& x) { return Numeric::exp_down(x); } 
      R log_down(const R& x) { return Numeric::log_down(x); }
      R sin_down(const R& x) { return Numeric::sin_down(x); }
      R cos_down(const R& x) { return Numeric::cos_down(x); }
      R tan_down(const R& x) { return Numeric::tan_down(x); }
      R asin_down(const R& x) { return Numeric::asin_down(x); }
      R acos_down(const R& x) { return Numeric::acos_down(x); }
      R atan_down(const R& x) { return Numeric::atan_down(x); }
      R sinh_down(const R& x) { return Numeric::sinh_down(x); }
      R cosh_down(const R& x) { return Numeric::cosh_down(x); }
      R tanh_down(const R& x) { return Numeric::tanh_down(x); }
      R asinh_down(const R& x) { return Numeric::asinh_down(x); }
      R acosh_down(const R& x) { return Numeric::acosh_down(x); }
      R atanh_down(const R& x) { return Numeric::atanh_down(x); }
    
      R neg_up(const R& x) { return Numeric::neg_up(x); } 
      R abs_up(const R& x) { return Numeric::abs_up(x); } 
      R add_up(const R& x1, const R& x2) { return Numeric::add_up(x1,x2); } 
      R sub_up(const R& x1, const R& x2) { return Numeric::sub_up(x1,x2); } 
      R mul_up(const R& x1, const R& x2) { return Numeric::mul_up(x1,x2); } 
      R div_up(const R& x1, const R& x2) { return Numeric::div_up(x1,x2); } 
      R sqrt_up(const R& x) { return Numeric::sqrt_up(x); } 
      R hypot_up(const R& x1, const R& x2) { return Numeric::hypot_up(x1,x2); } 
      
      R exp_up(const R& x) { return Numeric::exp_up(x); }
      R log_up(const R& x) { return Numeric::log_up(x); }
      R sin_up(const R& x) { return Numeric::sin_up(x); }
      R cos_up(const R& x) { return Numeric::cos_up(x); }
      R tan_up(const R& x) { return Numeric::tan_up(x); }
      R asin_up(const R& x) { return Numeric::asin_up(x); }
      R acos_up(const R& x) { return Numeric::acos_up(x); }
      R atan_up(const R& x) { return Numeric::atan_up(x); }
      R sinh_up(const R& x) { return Numeric::sinh_up(x); }
      R cosh_up(const R& x) { return Numeric::cosh_up(x); }
      R tanh_up(const R& x) { return Numeric::tanh_up(x); }
      R asinh_up(const R& x) { return Numeric::asinh_up(x); }
      R acosh_up(const R& x) { return Numeric::acosh_up(x); }
      R atanh_up(const R& x) { return Numeric::atanh_up(x); }
    };
    
  }
}

#endif /* _ARIADNE_ROUNDING_H */
