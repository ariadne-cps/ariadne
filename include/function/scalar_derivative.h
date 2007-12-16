/***************************************************************************
 *            scalar_derivative.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file scalar_derivative.h
 *  \brief Derivatives of scalar functions of a single variable.
 */
 
#ifndef ARIADNE_SCALAR_DERIVATIVE_H
#define ARIADNE_SCALAR_DERIVATIVE_H

#include <iostream>
#include <stdexcept>
#include <cassert>

#include "../base/tribool.h"
#include "../base/exceptions.h"
#include "../base/stlio.h"

#include "../numeric/exceptions.h"
#include "../numeric/traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

namespace Ariadne {
  namespace Function {


    /*!\ingroup Function
     * \brief \deprecated A templated class representing a the derivatives of a scalar quantity with respect to a scalar argument. (%Deprecated; use TaylorVariable instead)
     *
     * A scalar derivative is an array \f$(x,\dot{x},\ddot{x},\ldots)\f$ where \f$x\f$ is the value of the quantity,
     * and \f$\dot{x}\f$ is the variation with respect to some independent variable, \f$\ddot{x}\f$ is the second variation etc. 
     * A scalar derivative of the form \f$(c,0,0,\ldots)\f$ represents a constant. 
     * A derivative of the form \f$(x,1,0,0,\ldots)\f$ represents the independent variable with value \f$x\f$.
     *
     * The operation of a scalar function \f$f\f$ on the derivative \f$(x,\dot{x})\f$ is defined \f[f(x,\dot{x},\ddot{x}):=(f(x),f'(x)\dot{x},f''(x)\dot{x}^2+f'(x)\ddot{x},\ldots) . \f]
     *
     * To compute the \f$n^\mathrm{th}\f$ derivative of a scalar function \f$f\f$ at a point \f$x\f$, evaluate <tt>f(%Derivative<%Float>::variable(n,x)).%derivative()</tt>.
     * To compute range of derivatives of a scalar function \f$f\f$ over an interval \f$I\f$, evaluate <tt>f(%Derivative<%Interval>(I,1)).%derivative()</tt>.
     *
     * To construct a constant \f$c\f$, use <tt>ScalarDerivative(n,c) or ScalarDerivative::constant(n,c)</tt>.
     * To construct the derivative of a variable \f$x\f$, use <tt>ScalarDerivative(n,x,1)</tt>. or <tt>ScalarDerivative::variable(n,x)</tt>.
     * 
     */
    template<class X>
    class ScalarDerivative
    {
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      ScalarDerivative()
        : _data(1u) { }
      /*! \brief The constant zero of degree \a degree. */
      ScalarDerivative(smoothness_type degree)
        : _data(degree+1u) { }
      /*! \brief The constant \a constant of degree \a degree. */
      ScalarDerivative(smoothness_type degree, const X& constant)
        : _data(degree+1u) { this->_data[0]=constant; }
      /*! \brief A scalar derivative of degree \a degree, with value \a value and first derivative \a first_derivative. Higher derivatives are set to zero. */
      ScalarDerivative(smoothness_type degree, const X& value, const X& first_derivative)
        : _data(degree+1u) 
      { this->_data[0]=value; this->_data[1]=first_derivative; }
      /*! \brief A scalar derivative of degree \a degree, with value \a value and first derivative \a first_derivative. Higher derivatives are set to zero. */
      template<class XX> ScalarDerivative(smoothness_type degree, const XX* values)
        : _data(values,values+(degree+1u)) 
      { }
      /*! \brief A scalar derivative of degree \a degree, constructed from a list of values. */
      template<class XX> ScalarDerivative(smoothness_type degree, const XX* values, const XX* values2, const XX* values3)
        : _data(values,values+(degree+1u)) 
      { }

      /*! \brief A scalar derivative with values given by \a ary. */
      template<class XX> ScalarDerivative(const array<XX>& ary)
        : _data(ary) { }
      /*! \brief A scalar derivative with values given by array literal \a str. */
      ScalarDerivative(const std::string& str) { 
        std::stringstream ss(str); read_array(ss,_data); }
    
      /*! \brief Copy constructor. */
      template<class XX> ScalarDerivative(const ScalarDerivative<XX>& other) 
        : _data(other._data) { }
      /*! \brief Copy assignment operator. */
      template<class XX> ScalarDerivative<X>& operator=(const ScalarDerivative<XX>& other) {
        this->_data=other._data; return *this; }

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const ScalarDerivative<XX>& other) {
        return this->_data==other._data; }
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const ScalarDerivative<XX>& other) {
        return !(*this==other); }

      /*! \brief Construct a constant derivative of degree \a degree and value \a constant. */
      static ScalarDerivative<X> constant(smoothness_type degree, const X& constant) { 
        return ScalarDerivative<X>(degree,constant); }
      /*! \brief Construct the derivative of degree \a degree for the indepentent variable at value \a value. */
      static ScalarDerivative<X> variable(smoothness_type degree, const X& value) {
        return ScalarDerivative<X>(degree,value,1); }

      /*! \brief The degree (number of derivatives computed). */
      smoothness_type degree() const { 
        return this->_data.size()-1; }
      /*! \brief The array of derivative values. */
      const array<X>& data() const {
        return this->_data; }
      /*! \brief A reference to the array of derivative values. */
      array<X>& data() {
        return this->_data; }
      /*! \brief The \a i<sup> th</sup> derivative \f$d^if/dx^i\f$. */
      const X& derivative(smoothness_type i) const { 
        return this->_data[i]; }
      /*! \brief The value of the derived variable. */
      const X& value() const { 
        return this->_data[0]; }
      /*! \brief A reference to the \a i<sup> th</sup> derivative \f$d^if/dx^i\f$. */
      X& operator[](smoothness_type i) { 
        return this->_data[i]; }
      /*! \brief The \a i<sup> th</sup> derivative \f$d^if/dx^i\f$. */
      const X& operator[](smoothness_type i) const { 
        return this->_data[i]; }
#ifdef DOXYGEN
    //@{ 
    //! \name Friend operations
    /*! \brief The composition of two derivatives computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
     *  The composition inductively by
     *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
     */
    friend ScalarDerivative<X> compose(const ScalarDerivative<X>& y, const ScalarDerivative<X>& x);
    /*! \brief The derivative of the inverse of \f$y\f$ evaluated at \f$x\f$. (Not currently implemented.) */
    friend ScalarDerivative<X> inverse(const ScalarDerivative<X>& y, const X& x);
    /*! \brief The minimum of two derivatives. Returns the derivative whose zero-th order value is minimal. */
    friend ScalarDerivative<X> min(const ScalarDerivative<X>& x1, const ScalarDerivative<X>& x2);
    /*! \brief The maximum of two derivatives. Returns the derivative whose zero-th order value is maximal. */
    friend ScalarDerivative<X> max(const ScalarDerivative<X>& x1, const ScalarDerivative<X>& x2);
    /*! \brief The derivatives of \f$+x\f$. Returns a copy. */
    friend ScalarDerivative<X> pos(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$-x\f$. */
    friend ScalarDerivative<X> neg(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$x+y\f$. */
    friend ScalarDerivative<X> add(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend ScalarDerivative<X> sub(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x*y\f$. */
    friend ScalarDerivative<X> mul(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x/y\f$. */
    friend ScalarDerivative<X> div(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x^n\f$. */
    friend ScalarDerivative<X> pow(const ScalarDerivative<X>& x, const Integer& n);
    /*! \brief The derivatives of \f$\sqrt{x}\f$. */
    friend ScalarDerivative<X> sqrt(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\exp(x)\f$. */
    friend ScalarDerivative<X> exp(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\log(x)\f$. */
    friend ScalarDerivative<X> log(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\sin(x)\f$. */
    friend ScalarDerivative<X> sin(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\cos(x)\f$. */
    friend ScalarDerivative<X> cos(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\tan(x)\f$. */
    friend ScalarDerivative<X> tan(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\sin^{-1}(x)\f$. (Not currently implemented.) */
    friend ScalarDerivative<X> asin(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\cos^{-1}(x)\f$. (Not currently implemented.) */
    friend ScalarDerivative<X> acos(const ScalarDerivative<X>& x);
    /*! \brief The derivatives of \f$\tan^{-1}(x)\f$. (Not currently implemented.) */
    friend ScalarDerivative<X> atan(const ScalarDerivative<X>& x);

    /*! \brief The derivatives of \f$x+y\f$. */
    friend ScalarDerivative<X> operator+(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend ScalarDerivative<X> operator-(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x*y\f$. */
    friend ScalarDerivative<X> operator*(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
    /*! \brief The derivatives of \f$x/y\f$. */
    friend ScalarDerivative<X> operator/(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);

    /*! \brief The derivatives of \f$c+x\f$ for a constant \f$c\f$. (Other mixed-mode arithmetic is also supported.) */
    friend ScalarDerivative<X> operator+(const R& c, const ScalarDerivative<X>& x);

    /*! \brief Stream output operator. */
    friend std::ostream& operator<<(std::ostream& os, const ScalarDerivative<X>& x);
    //@}
#endif 
     private:
      array<X> _data;
    };


  template<class X> void compute_composition(ScalarDerivative<X>& y, const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> compose(const ScalarDerivative<X>& y, const ScalarDerivative<X>& x);

  template<class X> ScalarDerivative<X> inverse(const ScalarDerivative<X>& x, const X&);

  template<class X> ScalarDerivative<X> min(const ScalarDerivative<X>& x1, const ScalarDerivative<X>& x2); 
  template<class X> ScalarDerivative<X> max(const ScalarDerivative<X>& x1,const ScalarDerivative<X>& x2); 
  template<class X> ScalarDerivative<X> pos(const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> neg(const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> abs(const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> inv(const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> add(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X> ScalarDerivative<X> sub(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X> ScalarDerivative<X> mul(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X> ScalarDerivative<X> div(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X, class N> ScalarDerivative<X> pow(const ScalarDerivative<X>& x, N k);

  template<class X> ScalarDerivative<X> sqrt(const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> exp(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> log(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> sin(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> cos(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> tan(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> asin(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> acos(const ScalarDerivative<X>& x); 
  template<class X> ScalarDerivative<X> atan(const ScalarDerivative<X>& x); 

  template<class X> ScalarDerivative<X> operator-(const ScalarDerivative<X>& x);
  template<class X> ScalarDerivative<X> operator+(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X> ScalarDerivative<X> operator-(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X> ScalarDerivative<X> operator*(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);
  template<class X> ScalarDerivative<X> operator/(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y);

  template<class X, class R> ScalarDerivative<X> operator+(const ScalarDerivative<X>& x, const R& c);
  template<class X, class R> ScalarDerivative<X> operator+(const R& c, const ScalarDerivative<X>& x);
  template<class X, class R> ScalarDerivative<X> operator-(const ScalarDerivative<X>& x, const R& c);
  template<class X, class R> ScalarDerivative<X> operator-(const R& c, const ScalarDerivative<X>& x);
  template<class X, class R> ScalarDerivative<X> operator*(const ScalarDerivative<X>& x, const R& c);
  template<class X, class R> ScalarDerivative<X> operator*(const R& c, const ScalarDerivative<X>& x);
  template<class X, class R> ScalarDerivative<X> operator/(const ScalarDerivative<X>& x, const R& c);
  template<class X, class R> ScalarDerivative<X> operator/(const R& c, const ScalarDerivative<X>& x);

  template<class X> std::ostream& operator<<(std::ostream& os, const ScalarDerivative<X>& x);


  }
}


#include "scalar_derivative.inline.h"
#include "scalar_derivative.template.h"


#endif /* ARIADNE_SCALAR_DERIVATIVE_H */

