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
  
/*! \file numeric/function.h
 *  \brief Standard functions on real number types.
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

namespace Ariadne {
  namespace Numeric {

 #ifdef DOXYGEN

   //! \name Approximate arithmetical operations.
    //@{ 
    //! \ingroup Numeric

    /*! \brief Infinity and NotAnyNumber. */
    template<class R> inline R nan();
    template<class R> inline R inf();
    template<class R> inline R infinity();

    /*! \brief Rounding. */
    template<class R> inline R next_down(const R& x);
    template<class R> inline R next_up(const R& x);

    /*! \brief Minumum. */
    template<class R> inline R min_exact(const R& x1,const R& x2);
    template<class R> inline R min_approx(const R& x1,const R& x2);
    template<class R> inline R min_down(const R& x1,const R& x2);
    template<class R> inline R min_up(const R& x1,const R& x2);

    /*! \brief Maximum . */
    template<class R> inline R max_exact(const R& x1,const R& x2);
    template<class R> inline R max_approx(const R& x1,const R& x2);
    template<class R> inline R max_down(const R& x1,const R& x2);
    template<class R> inline R max_up(const R& x1,const R& x2);
   
    /*! \brief Median . */
    template<class R> inline R med_approx(const R& x1, const R& x2);
    template<class R> inline R med_exact(const R& x1, const R& x2);
   
    /*! \brief Unary neg_ation. */
    template<class R> inline R neg_exact(const R& x);
    template<class R> inline R neg_approx(const R& x);
    template<class R> inline R neg_down(const R& x);
    template<class R> inline R neg_up(const R& x);
   
    /*! \brief Reciprocal. */
    template<class R> inline R rec_approx(const R& x);
    template<class R> inline R rec_down(const R& x);
    template<class R> inline R rec_up(const R& x);
   
    /*! \brief Absolute value. */
    template<class R> inline R abs_exact(const R& x);
    template<class R> inline R abs_approx(const R& x);
    template<class R> inline R abs_down(const R& x);
    template<class R> inline R abs_up(const R& x);
   
    /*! \brief Addition. */
    template<class R> inline R add_approx(const R& x1,const R& x2);
    template<class R> inline R add_exact(const R& x1,const R& x2);
    template<class R> inline R add_down(const R& x1,const R& x2);
    template<class R> inline R add_up(const R& x1,const R& x2);
    
    /*! \brief Subtraction. */
    template<class R> inline R sub_approx(const R& x1,const R& x2);
    template<class R> inline R sub_exact(const R& x1,const R& x2);
    template<class R> inline R sub_down(const R& x1,const R& x2);
    template<class R> inline R sub_up(const R& x1,const R& x2);
    
     /*! \brief Multiplication. */
    template<class R> inline R mul_approx(const R& x1,const R& x2);
    template<class R> inline R mul_exact(const R& x1,const R& x2);
    template<class R> inline R mul_down(const R& x1,const R& x2);
    template<class R> inline R mul_up(const R& x1,const R& x2);
    
     /*! \brief Multiplication. */
    template<class R> inline R mul_approx(const R& x1,const int& x2);
    template<class R> inline R mul_exact(const R& x1,const int& n2);
    template<class R> inline R mul_down(const R& x1,const int& n2);
    template<class R> inline R mul_up(const R& x1,const int& n2);
    
    /*! \brief Division. */
    template<class R> inline R div_approx(const R& x1,const R& x2);
    template<class R> inline R div_exact(const R& x1,const R& x2);
    template<class R> inline R div_down(const R& x1,const R& x2);
    template<class R> inline R div_up(const R& x1,const R& x2);
    
    template<class R> inline R div_approx(const R& x1,const int& n2);
    template<class R> inline R div_exact(const R& x1,const int& n2);
    template<class R> inline R div_down(const R& x1,const int& n2);
    template<class R> inline R div_up(const R& x1,const int& n2);
    
    /*! \brief Power. */
    template<class R, class N> inline R pow_approx(const R& x1,const N& x2);
    template<class R, class N> inline R pow_exact(const R& x1,const N& x2);
    template<class R, class N> inline R pow_down(const R& x1,const N& x2);
    template<class R, class N> inline R pow_up(const R& x1,const N& x2);
    
    /*! \brief %Integer lower bound. */
    template<class N, class R> inline N int_down(const R& x);
    /*! \brief %Integer upper bound. */
    template<class N, class R> inline N int_up(const R& x);
    //@}

    //! \name Algebraic and transcendental functions.
    //@{
    //! \ingroup Numeric
    /*! \brief Square root. */
    template<class R> inline R sqrt_approx(const R& x);
    template<class R> inline R sqrt_down(const R& x);
    template<class R> inline R sqrt_up(const R& x);

    /*! \brief Hypoteneuse. */
    template<class R> inline R hypot_approx(const R& x1,const R& x2);
    template<class R> inline R hypot_down(const R& x1,const R& x2);
    template<class R> inline R hypot_up(const R& x1,const R& x2);
   
    /*! \brief Exponential. */
    template<class R> inline R exp_approx(const R& x);
    template<class R> inline R exp_down(const R& x);
    template<class R> inline R exp_up(const R& x);

    /*! \brief Natural logarithm. */
    template<class R> inline R log_approx(const R& x);
    template<class R> inline R log_down(const R& x);
    template<class R> inline R log_up(const R& x);

    /*! \brief The constant \f$\pi\f$. */
    template<typename R> R pi_approx();
    template<typename R> R pi_down();
    template<typename R> R pi_up();

    /*! \brief Sine function. */
    template<class R> inline R sin_approx(const R& x);
    template<class R> inline R sin_down(const R& x);
    template<class R> inline R sin_up(const R& x);

    /*! \brief Cosine function. */
    template<class R> inline R cos_approx(const R& x);
    template<class R> inline R cos_down(const R& x);
    template<class R> inline R cos_up(const R& x);

    /*! \brief Tangent function. */
    template<class R> inline R tan_approx(const R& x);
    template<class R> inline R tan_down(const R& x);
    template<class R> inline R tan_up(const R& x);

    /*! \brief Hyperbolic sine function. */
    template<class R> inline R sinh_approx(const R& x);
    template<class R> inline R sinh_down(const R& x);
    template<class R> inline R sinh_up(const R& x);

    /*! \brief Hyperbolic cosine function. */
    template<class R> inline R cosh_approx(const R& x);
    template<class R> inline R cosh_down(const R& x);
    template<class R> inline R cosh_up(const R& x);

    /*! \brief Hyperbolic tangent function. */
    template<class R> inline R tanh_approx(const R& x);
    template<class R> inline R tanh_down(const R& x);
    template<class R> inline R tanh_up(const R& x);

    /*! \brief Inverse sine function. */
    template<class R> inline R asin_approx(const R& x);
    template<class R> inline R asin_down(const R& x);
    template<class R> inline R asin_up(const R& x);

    /*! \brief Inverse cosine function. */
    template<class R> inline R acos_approx(const R& x);
    template<class R> inline R acos_down(const R& x);
    template<class R> inline R acos_up(const R& x);

    /*! \brief Inverse tangent function. */
    template<class R> inline R atan_approx(const R& x);
    template<class R> inline R atan_down(const R& x);
    template<class R> inline R atan_up(const R& x);

    /*! \brief Inverse hyperbolic sine function. */
    template<class R> inline R asinh_approx(const R& x);
    template<class R> inline R asinh_down(const R& x);
    template<class R> inline R asinh_up(const R& x);

    /*! \brief Inverse hyperbolic cosine function. */
    template<class R> inline R acosh_approx(const R& x);
    template<class R> inline R acosh_down(const R& x);
    template<class R> inline R acosh_up(const R& x);

    /*! \brief Inverse hyperbolic tangent function. */
    template<class R> inline R atanh_approx(const R& x);
    template<class R> inline R atanh_down(const R& x);
    template<class R> inline R atanh_up(const R& x);
  //@}

#endif

  }
}

#endif /* ARIADNE_FUNCTION_H */
