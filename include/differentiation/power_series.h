/***************************************************************************
 *            power_series.h
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
 
/*! \file power_series.h
 *  \brief Univariate Taylor series.
 */
 
#ifndef ARIADNE_POWER_SERIES_H
#define ARIADNE_POWER_SERIES_H

#include "base/array.h"

namespace Ariadne {
  
    template<class X> class Differential;

    /*!\ingroup Differentiation
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a single argument.
     */
    template<class X>
    class PowerSeries
    {
     private:
      array<X> _data;
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      PowerSeries();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      PowerSeries(smoothness_type d);
      /*! \brief A taylor variable of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> PowerSeries(smoothness_type d, const XX* ptr);

      /*! \brief Copy constructor. */
      template<class XX> PowerSeries(const PowerSeries<XX>& ts); 
      /*! \brief Copy assignment operator. */
      template<class XX> PowerSeries<X>& operator=(const PowerSeries<XX>& ts);

      /*! \brief Assign a constant \a c. */
      template<class XX> PowerSeries<X>& operator=(const XX& c);

      /*! \brief Construct a constant series of degree \a d with value zero. */
      static PowerSeries<X> zero(smoothness_type d); 
      /*! \brief Construct a constant serues of degree \a d with value \a c. */
      static PowerSeries<X> constant(smoothness_type d, const X& c); 
      /*! \brief Construct the variable of degree \a d at value \a x with respect to a single scalar variable. */
      static PowerSeries<X> variable(smoothness_type d, const X& x);

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
     private:
      static void instantiate();
    };


  
    template<class X> PowerSeries<X> compose(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> inverse(const PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X> derivative(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> antiderivative(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> antiderivative(const PowerSeries<X>& ts, const X& c);

    template<class X> PowerSeries<X> min(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> max(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> abs(const PowerSeries<X>& ts);

    template<class X> PowerSeries<X> pos(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> neg(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> rec(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> add(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> sub(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> mul(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> div(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> pow(const PowerSeries<X>& ts, const uint& n);
    template<class X> PowerSeries<X> pow(const PowerSeries<X>& ts, const int& n);

    template<class X> PowerSeries<X> sqrt(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> exp(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> log(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> sin(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> cos(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> tan(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> asin(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> acos(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> atan(const PowerSeries<X>& ts);
    
    template<class X> PowerSeries<X> operator+(const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator-(const PowerSeries<X>& ts);

    template<class X> PowerSeries<X> operator+(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> operator-(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> operator*(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);
    template<class X> PowerSeries<X> operator/(const PowerSeries<X>& ts1, const PowerSeries<X>& ts2);

    template<class X> PowerSeries<X>& operator+=(PowerSeries<X>& ts1, const PowerSeries<X>& ts2);

    template<class X> PowerSeries<X> operator+(const PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X> operator-(const PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X> operator*(const PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X> operator/(const PowerSeries<X>& ts, const X& c);

    template<class X> PowerSeries<X> operator+(const X& c, const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator-(const X& c, const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator*(const X& c, const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator/(const X& c, const PowerSeries<X>& ts);

    template<class X> PowerSeries<X> operator+(const PowerSeries<X>& ts, const double& c);
    template<class X> PowerSeries<X> operator-(const PowerSeries<X>& ts, const double& c);
    template<class X> PowerSeries<X> operator*(const PowerSeries<X>& ts, const double& c);
    template<class X> PowerSeries<X> operator/(const PowerSeries<X>& ts, const double& c);

    template<class X> PowerSeries<X> operator+(const double& c, const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator-(const double& c, const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator*(const double& c, const PowerSeries<X>& ts);
    template<class X> PowerSeries<X> operator/(const double& c, const PowerSeries<X>& ts);

    template<class X> PowerSeries<X>& operator+=(PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X>& operator-=(PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X>& operator*=(PowerSeries<X>& ts, const X& c);
    template<class X> PowerSeries<X>& operator/=(PowerSeries<X>& ts, const X& c);

    template<class X> PowerSeries<X>& operator+=(PowerSeries<X>& ts, const double& c);
    template<class X> PowerSeries<X>& operator-=(PowerSeries<X>& ts, const double& c);
    template<class X> PowerSeries<X>& operator*=(PowerSeries<X>& ts, const double& c);
    template<class X> PowerSeries<X>& operator/=(PowerSeries<X>& ts, const double& c);

    template<class X, class R> PowerSeries<X>& operator*=(PowerSeries<X>& ts, const R& c);


    template<class X> std::ostream& operator<<(std::ostream& os, const PowerSeries<X>& ts);
  
} // namespace Ariadne



#include "power_series.inline.h"


#endif /* ARIADNE_POWER_SERIES_H */

