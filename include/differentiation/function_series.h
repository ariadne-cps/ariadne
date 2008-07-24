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


namespace Ariadne {

   

  
    template<class X> class PowerSeries;

    template<class X> 
    class ArithmeticSeries
    {
     public:
      /*! \brief Construct the Taylor series of the reciprocal function. */
      static PowerSeries<X> rec(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the \a n<sup>th</sup> power function. */
      static PowerSeries<X> pow(smoothness_type d, const X& c, const uint& k); 
    };

    template<class X> 
    class TranscendentalSeries
    {
     public:
      /*! \brief Construct the Taylor series of the square-root function. */
      static PowerSeries<X> sqrt(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the exponential function. */
      static PowerSeries<X> exp(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the logarithm function. */
      static PowerSeries<X> log(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the sine function. */
      static PowerSeries<X> sin(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the cosine function. */
      static PowerSeries<X> cos(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the tangent function. */
      static PowerSeries<X> tan(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inverse sine function. */
      static PowerSeries<X> asin(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inversecosine function. */
      static PowerSeries<X> acos(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inversetangent function. */
      static PowerSeries<X> atan(smoothness_type d, const X& c); 
    };


    template<class X>
    class FunctionSeries
      : public ArithmeticSeries<X>,
        public TranscendentalSeries<X>
    { };

    template<>
    class FunctionSeries<Rational>
      : public ArithmeticSeries<Rational>
    { };


    
  
} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_SERIES_H */

