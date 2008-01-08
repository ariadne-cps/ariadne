/***************************************************************************
 *            taylor_variable.h
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
 
/*! \file taylor_variable.h
 *  \brief Derivatives of scalar functions of many variables.
 */
 
#ifndef ARIADNE_TAYLOR_VARIABLE_H
#define ARIADNE_TAYLOR_VARIABLE_H


namespace Ariadne {
  namespace Function {
  
    class MultiIndex;
    template<class X> class TaylorSeries;
    template<class X> class TaylorVariable;
    template<class X> class TaylorDerivative;
  
    /*!\ingroup Function
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a multiple arguments.
     */
    template<class X>
    class TaylorVariable
    {
      friend class TaylorDerivative<X>;
     public:
      /*! The type used to represent numbers. */
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      TaylorVariable();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      TaylorVariable(size_type a, smoothness_type d);
      /*! \brief A taylor variable of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> TaylorVariable(size_type a, smoothness_type d, const XX* ptr);

      /*! \brief Construct from a univariate Taylor series. */
      template<class XX> TaylorVariable(const TaylorSeries<XX>& ts); 
      /*! \brief Copy constructor. */
      template<class XX> TaylorVariable(const TaylorVariable<XX>& tv); 
      /*! \brief Copy assignment operator. */
      template<class XX> TaylorVariable<X>& operator=(const TaylorVariable<XX>& tv);

      /*! \brief Assign a constant \a c. */
      template<class XX> TaylorVariable<X>& operator=(const XX& c);

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const TaylorVariable<XX>& other) const;
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const TaylorVariable<XX>& other) const;

      /*! \brief Construct a constant variable of degree \a d with respect to \a as variables and value \a c. */
      template<class XX> static TaylorVariable<X> constant(size_type as, smoothness_type d, const XX& c); 
      /*! \brief Construct the variable of degree \a d at value \a value with respect to the \a i<sup>th</sup> variable of \a as. */
      template<class XX> static TaylorVariable<X> variable(size_type as, smoothness_type d, const XX& value, size_type i);

      /*! \brief The number of variables of the argument. */
      size_type argument_size() const; 
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
      X& operator[](const MultiIndex& a); 
      /*! \brief The \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$. */
      const X& operator[](const MultiIndex& a) const; 

      /*! \brief Inplace addition. */
      TaylorVariable<X>& operator+=(const TaylorVariable<X>& x);
      /*! \brief Inplace multiplication. */
      TaylorVariable<X>& operator*=(const X& x);
      /*! \brief Inplace division. */
      TaylorVariable<X>& operator/=(const X& x);
#ifdef DOXYGEN
      //@{ 
      //! \name Friend operations
      /*! \brief The composition of two variables computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
       *  The composition inductively by
       *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
       *  \f$ y = a_0 + x ( a_1 + x ( a_2/2 + x ( a_3/3! + \cdots)))\f$.
       */
      friend TaylorVariable<X> compose(const TaylorVariable<X>& y, const TaylorVariable<X>& x);
      /*! \brief The derivatives of the inverse of \f$y\f$ evaluated at \f$x\f$. (Not currently implemented.) */
      friend TaylorVariable<X> inverse(const TaylorVariable<X>& y, const X& x);
      /*! \brief The derivative of \f$x\f$ with respect to the variable \a k .*/
      friend TaylorVariable<X> derivative(const TaylorVariable<x>& x, const size_type& k);
      /*! \brief The minimum of two variables. Returns the variable whose zero-th order value is minimal. */
      friend TaylorVariable<X> min(const TaylorVariable<X>& x1, const TaylorVariable<X>& x2);
      /*! \brief The maximum of two variables. Returns the variable whose zero-th order value is maximal. */
      friend TaylorVariable<X> max(const TaylorVariable<X>& x1, const TaylorVariable<X>& x2);
      /*! \brief The derivatives of \f$+x\f$. Returns a copy. */
      friend TaylorVariable<X> pos(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$-x\f$. */
      friend TaylorVariable<X> neg(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$x+y\f$. */
      friend TaylorVariable<X> add(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x-y\f$. */
      friend TaylorVariable<X> sub(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x*y\f$. */
      friend TaylorVariable<X> mul(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x/y\f$. */
      friend TaylorVariable<X> div(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x^n\f$. */
      friend TaylorVariable<X> pow(const TaylorVariable<X>& x, const Integer& n);
      /*! \brief The derivatives of \f$\sqrt{x}\f$. */
      friend TaylorVariable<X> sqrt(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\exp(x)\f$. */
      friend TaylorVariable<X> exp(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\log(x)\f$. */
      friend TaylorVariable<X> log(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\sin(x)\f$. */
      friend TaylorVariable<X> sin(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\cos(x)\f$. */
      friend TaylorVariable<X> cos(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\tan(x)\f$. */
      friend TaylorVariable<X> tan(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\sin^{-1}(x)\f$. (Not currently implemented.) */
      friend TaylorVariable<X> asin(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\cos^{-1}(x)\f$. (Not currently implemented.) */
      friend TaylorVariable<X> acos(const TaylorVariable<X>& x);
      /*! \brief The derivatives of \f$\tan^{-1}(x)\f$. (Not currently implemented.) */
      friend TaylorVariable<X> atan(const TaylorVariable<X>& x);

      /*! \brief The derivatives of \f$x+y\f$. */
      friend TaylorVariable<X> operator+(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x-y\f$. */
      friend TaylorVariable<X> operator-(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x*y\f$. */
      friend TaylorVariable<X> operator*(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
      /*! \brief The derivatives of \f$x/y\f$. */
      friend TaylorVariable<X> operator/(const TaylorVariable<X>& x, const TaylorVariable<X>& y);

      /*! \brief The derivatives of \f$c+x\f$ for a constant \f$c\f$. (Other mixed-mode arithmetic is also supported.) */
      friend TaylorVariable<X> operator+(const R& c, const TaylorVariable<X>& x);

      /*! \brief Stream output operator. */
      friend std::ostream& operator<<(std::ostream& os, const TaylorVariable<X>& x);
      //@}
#endif 
     private:
      static void instantiate();
     private:
      size_type _argument_size;
      smoothness_type _degree;
      array<X> _data;
    };


  template<class X> TaylorVariable<X> compose(const TaylorSeries<X>& y, const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> compose(const TaylorVariable<X>& y, const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> reduce(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> derivative(const TaylorVariable<X>& x, const size_type& k);
  template<class X> TaylorVariable<X> min(const TaylorVariable<X>& x1, const TaylorVariable<X>& x2); 
  template<class X> TaylorVariable<X> max(const TaylorVariable<X>& x1,const TaylorVariable<X>& x2); 
  template<class X> TaylorVariable<X> pos(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> neg(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> abs(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> rec(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> add(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X> TaylorVariable<X> sub(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X> TaylorVariable<X> mul(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X> TaylorVariable<X> div(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X, class N> TaylorVariable<X> pow(const TaylorVariable<X>& x, N k);

  template<class X> TaylorVariable<X> sqrt(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> exp(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> log(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> sin(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> cos(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> tan(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> asin(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> acos(const TaylorVariable<X>& x); 
  template<class X> TaylorVariable<X> atan(const TaylorVariable<X>& x); 

  template<class X> TaylorVariable<X> operator+(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> operator-(const TaylorVariable<X>& x);
  template<class X> TaylorVariable<X> operator+(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X> TaylorVariable<X> operator-(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X> TaylorVariable<X> operator*(const TaylorVariable<X>& x, const TaylorVariable<X>& y);
  template<class X> TaylorVariable<X> operator/(const TaylorVariable<X>& x, const TaylorVariable<X>& y);

  template<class X, class R> TaylorVariable<X> operator+(const TaylorVariable<X>& x, const R& c);
  template<class X, class R> TaylorVariable<X> operator+(const R& c, const TaylorVariable<X>& x);
  template<class X, class R> TaylorVariable<X> operator-(const TaylorVariable<X>& x, const R& c);
  template<class X, class R> TaylorVariable<X> operator-(const R& c, const TaylorVariable<X>& x);
  template<class X, class R> TaylorVariable<X> operator*(const TaylorVariable<X>& x, const R& c);
  template<class X, class R> TaylorVariable<X> operator*(const R& c, const TaylorVariable<X>& x);
  template<class X, class R> TaylorVariable<X> operator/(const TaylorVariable<X>& x, const R& c);
  template<class X, class R> TaylorVariable<X> operator/(const R& c, const TaylorVariable<X>& x);

  template<class X, class R> TaylorVariable<X> operator/=(const TaylorVariable<X>& x, const R& c);

  template<class X> std::ostream& operator<<(std::ostream& os, const TaylorVariable<X>& x);


  }
}


#include "taylor_variable.inline.h"
#include "taylor_variable.template.h"


#endif /* ARIADNE_TAYLOR_VARIABLE_H */

