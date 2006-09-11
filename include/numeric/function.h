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

namespace Ariadne {
  namespace Numeric {

    //! \name Approximate arithmetical operations.
    //@{ 
    //! \ingroup Numeric
    /*! \brief Minumum. */
    template<typename R> R min_exact(const R& x1,const R& x2);

    /*! \brief Maximum . */
    template<typename R> R max_exact(const R& x1,const R& x2);
   
     /*! \brief Median . */
    template<typename R> R med_approx(const R& x1, const R& x2);
    template<typename R> R med_exact(const R& x1, const R& x2);
   
   /*! \brief Unary negation. */
    template<typename R> R neg_exact(const R& x);
    template<typename R> R neg_down(const R& x);
    template<typename R> R neg_up(const R& x);
   
    /*! \brief Absolute value. */
    template<typename R> R abs_exact(const R& x);
    template<typename R> R abs_down(const R& x);
    template<typename R> R abs_up(const R& x);
   
    /*! \brief Addition. */
    template<typename R> R add_approx(const R& x1,const R& x2);
    template<typename R> R add_exact(const R& x1,const R& x2);
    template<typename R> R add_down(const R& x1,const R& x2);
    template<typename R> R add_up(const R& x1,const R& x2);
    
    /*! \brief Subtraction. */
    template<typename R> R sub_approx(const R& x1,const R& x2);
    template<typename R> R sub_exact(const R& x1,const R& x2);
    template<typename R> R sub_down(const R& x1,const R& x2);
    template<typename R> R sub_up(const R& x1,const R& x2);
    
     /*! \brief Multiplication. */
    template<typename R> R mul_approx(const R& x1,const R& x2);
    template<typename R> R mul_exact(const R& x1,const R& x2);
    template<typename R> R mul_down(const R& x1,const R& x2);
    template<typename R> R mul_up(const R& x1,const R& x2);
    
    /*! \brief Division. */
    template<typename R> R div_approx(const R& x1,const R& x2);
    template<typename R> R div_exact(const R& x1,const R& x2);
    template<typename R> R div_down(const R& x1,const R& x2);
    template<typename R> R div_up(const R& x1,const R& x2);
    
    /*! \brief %Integer lower bound. */
    template<typename N, typename R> N int_down(const R& x);
    /*! \brief %Integer upper bound. */
    template<typename N, typename R> N int_up(const R& x);
    //@}

    //! \name Algebraic and transcendental functions.
    //@{
    //! \ingroup Numeric
    /*! \brief Square root. */
    template<typename R> R sqrt_approx(const R& x);
    template<typename R> R sqrt_down(const R& x);
    template<typename R> R sqrt_up(const R& x);

    /*! \brief Hypoteneuse. */
    template<typename R> R hypot_approx(const R& x1,const R& x2);
    template<typename R> R hypot_down(const R& x1,const R& x2);
    template<typename R> R hypot_up(const R& x1,const R& x2);
   
    /*! \brief Exponential. */
    template<typename R> R exp_approx(const R& x);
    template<typename R> R exp_down(const R& x);
    template<typename R> R exp_up(const R& x);

    /*! \brief Natural logarithm. */
    template<typename R> R log_approx(const R& x);
    template<typename R> R log_down(const R& x);
    template<typename R> R log_up(const R& x);

    /*! \brief Sine function. */
    template<typename R> R sin_approx(const R& x);
    template<typename R> R sin_down(const R& x);
    template<typename R> R sin_up(const R& x);

    /*! \brief Cosine function. */
    template<typename R> R cos_approx(const R& x);
    template<typename R> R cos_down(const R& x);
    template<typename R> R cos_up(const R& x);

    /*! \brief Tangent function. */
    template<typename R> R tan_approx(const R& x);
    template<typename R> R tan_down(const R& x);
    template<typename R> R tan_up(const R& x);

    /*! \brief Hyperbolic sine function. */
    template<typename R> R sinh_approx(const R& x);
    template<typename R> R sinh_down(const R& x);
    template<typename R> R sinh_up(const R& x);

    /*! \brief Hyperbolic cosine function. */
    template<typename R> R cosh_approx(const R& x);
    template<typename R> R cosh_down(const R& x);
    template<typename R> R cosh_up(const R& x);

    /*! \brief Hyperbolic tangent function. */
    template<typename R> R tanh_approx(const R& x);
    template<typename R> R tanh_down(const R& x);
    template<typename R> R tanh_up(const R& x);

    /*! \brief Inverse sine function. */
    template<typename R> R asin_approx(const R& x);
    template<typename R> R asin_down(const R& x);
    template<typename R> R asin_up(const R& x);

    /*! \brief Inverse cosine function. */
    template<typename R> R acos_approx(const R& x);
    template<typename R> R acos_down(const R& x);
    template<typename R> R acos_up(const R& x);

    /*! \brief Inverse tangent function. */
    template<typename R> R atan_approx(const R& x);
    template<typename R> R atan_down(const R& x);
    template<typename R> R atan_up(const R& x);

    /*! \brief Inverse hyperbolic sine function. */
    template<typename R> R asinh_approx(const R& x);
    template<typename R> R asinh_down(const R& x);
    template<typename R> R asinh_up(const R& x);

    /*! \brief Inverse hyperbolic cosine function. */
    template<typename R> R acosh_approx(const R& x);
    template<typename R> R acosh_down(const R& x);
    template<typename R> R acosh_up(const R& x);

    /*! \brief Inverse hyperbolic tangent function. */
    template<typename R> R atanh_approx(const R& x);
    template<typename R> R atanh_down(const R& x);
    template<typename R> R atanh_up(const R& x);
  //@}

  }
}

#endif /* _ARIADNE_FUNCTION_H */
