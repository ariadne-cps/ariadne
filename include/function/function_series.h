/***************************************************************************
 *            function_series.h
 *
 *  Copyright 2007  Pieter Collins
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
 
/*! \file function_series.h
 *  \brief Univariate Taylor series for univariate functions.
 */
 
#ifndef ARIADNE_FUNCTION_SERIES_H
#define ARIADNE_FUNCTION_SERIES_H

#include "numeric/declarations.h"

namespace Ariadne {
  namespace Function {
  
    template<class X> class TaylorSeries;

    template<class X> 
    class ArithmeticSeries
    {
     public:
      /*! \brief Construct the Taylor series of the reciprocal function. */
      static TaylorSeries<X> rec(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the \a n<sup>th</sup> power function. */
      static TaylorSeries<X> pow(smoothness_type d, const X& c, const uint& k); 
    };

    template<class X> 
    class TranscendentalSeries
    {
     public:
      /*! \brief Construct the Taylor series of the square-root function. */
      static TaylorSeries<X> sqrt(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the exponential function. */
      static TaylorSeries<X> exp(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the logarithm function. */
      static TaylorSeries<X> log(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the sine function. */
      static TaylorSeries<X> sin(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the cosine function. */
      static TaylorSeries<X> cos(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the tangent function. */
      static TaylorSeries<X> tan(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inverse sine function. */
      static TaylorSeries<X> asin(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inversecosine function. */
      static TaylorSeries<X> acos(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inversetangent function. */
      static TaylorSeries<X> atan(smoothness_type d, const X& c); 
    };


    template<class X>
    class FunctionSeries
      : public ArithmeticSeries<X>,
        public TranscendentalSeries<X>
    { };

    template<>
    class FunctionSeries<Numeric::Rational>
      : public ArithmeticSeries<Numeric::Rational>
    { };


    
  }
}

#endif /* ARIADNE_FUNCTION_SERIES_H */

