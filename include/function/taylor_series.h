/***************************************************************************
 *            taylor_series.h
 *
 *  Copyright 2007  Pieter Collins
 *  Email Pieter.Collins@cwi.nl
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
 
/*! \file taylor_series.h
 *  \brief Univariate Taylor series.
 */
 
#ifndef ARIADNE_TAYLOR_SERIES_H
#define ARIADNE_TAYLOR_SERIES_H

#include <iostream>
#include <stdexcept>
#include <cassert>

#include "base/tribool.h"
#include "base/exceptions.h"
#include "base/stlio.h"

#include "numeric/exceptions.h"
#include "numeric/traits.h"
#include "numeric/conversion.h"
#include "numeric/arithmetic.h"
#include "numeric/function.h"

namespace Ariadne {
  namespace Function {
  
    /*!\ingroup Function
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a multiple arguments.
     */
    template<class X>
    class TaylorSeries
    {
     private:
      array<X> _data;
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      TaylorSeries();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      TaylorSeries(smoothness_type d);
      /*! \brief A taylor variable of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> TaylorSeries(smoothness_type d, const XX* ptr);

      /*! \brief Copy constructor. */
      template<class XX> TaylorSeries(const TaylorSeries<XX>& ts); 
      /*! \brief Copy assignment operator. */
      template<class XX> TaylorSeries<X>& operator=(const TaylorSeries<XX>& ts);

      /*! \brief Assign a constant \a c. */
      template<class XX> TaylorSeries<X>& operator=(const XX& c);

      /*! \brief Construct a constant variable of degree \a d with respect to \a as variables and value \a c. */
      template<class XX> static TaylorSeries<X> constant(smoothness_type d, const XX& c); 
      /*! \brief Construct the variable of degree \a d at value \a x with respect to the \a i<sup>th</sup> variable of \a as. */
      template<class XX> static TaylorSeries<X> variable(smoothness_type d, const XX& x);

      /*! \brief The degree (number of derivatives computed). */
      smoothness_type degree() const; 
      /*! \brief The value of the quantity. */
      const X& value() const;
      /*! \brief A reference to the value of the quantity. */
      X& value();
      /*! \brief The array of derivative values. */
      const array<X>& data() const;
      /*! \brief A reference to the array of derivative values. */
      array<X>& data();
      /*! \brief A reference to the \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$. */
      X& operator[](const smoothness_type& d); 
      /*! \brief The \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$. */
      const X& operator[](const smoothness_type& a) const; 

      /*! \brief Construct the Taylor series of the reciprocal function. */
      static TaylorSeries<X> rec(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the \a n<sup>th</sup> power function. */
      static TaylorSeries<X> pow(smoothness_type d, const X& c, const uint& k); 
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
     private:
      static void instantiate();
    };


  
    template<>
    class TaylorSeries<Numeric::Rational>
    {
      typedef Numeric::Rational X;
     private:
      array<X> _data;
     public:
      typedef X value_type;

      TaylorSeries();
      TaylorSeries(smoothness_type d);
      template<class XX> TaylorSeries(smoothness_type d, const XX* ptr);

      template<class XX> TaylorSeries(const TaylorSeries<XX>& ts); 
      template<class XX> TaylorSeries<X>& operator=(const TaylorSeries<XX>& ts);

      template<class XX> TaylorSeries<X>& operator=(const XX& c);

      static TaylorSeries<X> constant(smoothness_type d, const X& c); 
      static TaylorSeries<X> variable(smoothness_type d, const X& x);

      smoothness_type degree() const; 
      const X& value() const;
      X& value();
      const array<X>& data() const;
      array<X>& data();
      X& operator[](const smoothness_type& d); 
      const X& operator[](const smoothness_type& a) const; 

      static TaylorSeries<X> rec(smoothness_type d, const X& c); 
      static TaylorSeries<X> pow(smoothness_type d, const X& c, const uint& k); 
     private:
      static void instantiate();
    };


    template<class X> TaylorSeries<X> compose(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> inverse(const TaylorSeries<X>& ts, const X& c);

    template<class X> TaylorSeries<X> min(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> max(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> abs(const TaylorSeries<X>& ts);

    template<class X> TaylorSeries<X> pos(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> neg(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> rec(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> add(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> sub(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> mul(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> div(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> pow(const TaylorSeries<X>& ts, const uint& n);
    template<class X> TaylorSeries<X> pow(const TaylorSeries<X>& ts, const int& n);

    template<class X> TaylorSeries<X> sqrt(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> exp(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> log(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> sin(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> cos(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> tan(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> asin(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> acos(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> atan(const TaylorSeries<X>& ts);
    
    template<class X> TaylorSeries<X> operator+(const TaylorSeries<X>& ts);
    template<class X> TaylorSeries<X> operator-(const TaylorSeries<X>& ts);

    template<class X> TaylorSeries<X> operator+(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> operator-(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> operator*(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);
    template<class X> TaylorSeries<X> operator/(const TaylorSeries<X>& ts1, const TaylorSeries<X>& ts2);

    template<class X, class R> TaylorSeries<X> operator+(const TaylorSeries<X>& ts, const R& c);
    template<class X, class R> TaylorSeries<X> operator-(const TaylorSeries<X>& ts, const R& c);
    template<class X, class R> TaylorSeries<X> operator*(const TaylorSeries<X>& ts, const R& c);
    template<class X, class R> TaylorSeries<X> operator/(const TaylorSeries<X>& ts, const R& c);
    template<class X, class R> TaylorSeries<X>& operator/=(TaylorSeries<X>& ts, const R& c);

    template<class X, class R> TaylorSeries<X> operator+(const R& c, const TaylorSeries<X>& ts);
    template<class X, class R> TaylorSeries<X> operator-(const R& c, const TaylorSeries<X>& ts);
    template<class X, class R> TaylorSeries<X> operator*(const R& c, const TaylorSeries<X>& ts);
    template<class X, class R> TaylorSeries<X> operator/(const R& c, const TaylorSeries<X>& ts);

    template<class X> std::ostream& operator<<(std::ostream& os, const TaylorSeries<X>& ts);
  }
}



#include "taylor_series.inline.h"
#include "taylor_series.template.h"


#endif /* ARIADNE_TAYLOR_SERIES_H */

