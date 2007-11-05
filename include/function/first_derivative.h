/***************************************************************************
 *            first_derivative.h
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
 
/*! \file first_derivative.h
 *  \brief First derivatives with respect to multiple arguments.
 */
 
#ifndef ARIADNE_FIRST_DERIVATIVE_H
#define ARIADNE_FIRST_DERIVATIVE_H

#include <iostream>
#include <stdexcept>
#include <cassert>

#include "../base/tribool.h"
#include "../base/exceptions.h"
#include "../base/stlio.h"

#include "../numeric/exceptions.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

namespace Ariadne {
  namespace Function {

    /*!\ingroup Function
     * \brief A templated class representing a quantity and its first variation or differential with respect to one or more independent variables.
     *
     * A first derivative \f$(x,v)\f$ or \f$(x,\dot{x})\f$ represents a quantity whose value is \f$x\f$ and whose derivative with respect to the ith independent
     * variable is \f$\partial x/\partial t_i=v_i\f$. A
     * A first derivative of the form \f$(c,0)\f$ represents a constant. 
     * A first derivative of the form \f$(x,e_i)\f$, where \f$e_i\f$ is the ith unit vector, represents the independent variable \f$x_i\f$.
     *
     * Addition of first derivatives is given by \f[(x,v)+(y,w):=(x+y,v+w) . \f]
     * By the Leibnitz rule, \f[(x,v)\times(y,w):=(xy,\;yv+xw) . \f]
     * The action of a scalar function \f$f\f$ on the differential \f$(x,v)\f$ is defined \f[f(x,v):=(f(x),\;f'(x)v) . \f]
     * The chain rule is automatically obtained: \f[(g\circ f)(x,v)=\bigl(g(f(x)),\;g'(f(x))f'(x)v\bigr) . \f]
     *
     * The constructor <tt>FirstDerivative(x,v)</tt> constructs a variable whose value is \f$x\f$ and whose partial derivative with respect to the
     * \f$i\f$<sup>th</sup> independent variable is <tt>v[i]</tt> .
     * To construct a "constant" quantity with respect to \f$n\f$ arguments, use <tt>%FirstDerivative(c, zero_vector(n))</tt>.
     * To construct the \f$i\f$<sup>th</sup> "variable" \f$x\f$ of \a n, use <tt>%FirstDerivative(x,unit_vector(n,i))</tt>.
     *
     *
     * \internal 
     * Maybe a better name for the class would be "Variation" or "Derivative". 
     * The name "Differential" could then be used for backward differentiation. 
     * Alternatively, "DifferentialForm" could be used for backward differentiation.
     *
     * This class needs to be specialised or extended to provide better support for vector-valued functions.
     */
    template<class X, class V>
    class FirstDerivative
    {
     private:
      X _x; V _dx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs the differential of the constant zero with respect to no independent variable. */
      FirstDerivative() : _x(), _dx() { }
      /*! \brief Constuct a differential based on variable x and derivative vector dx. */
      template<class XX,class VV> FirstDerivative(const XX& x, const VV& dx) : _x(x), _dx(dx) { }
      /*! \brief Construct the differential of a constant with respect to \a n independent variables. (Deprecated) */
      template<class XX> FirstDerivative(const XX& x, size_type n) : _x(x), _dx(n) { }
      /*! \brief Construct the differential of the \a i th variable of \a n independent variables. (Deprecated) */
      template<class XX> FirstDerivative(const XX& x, size_type n, size_type i) : _x(x), _dx(n) { _dx[i]=1; }
      /*! \brief Assign the differential of a constant zero with respect to the current independent variables (useful as an alternative to a constructor).*/
      template<class XX> FirstDerivative<X,V> operator=(const XX& c) { 
        this->_x=c; this->_dx=this->_dx*0; return *this; }
      //@}

      //@{
      //! \name Data access
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The differential of the variable with respect to the indepenent variables. */
      const V& derivative() const { return this->_dx; }
      /*! \brief The differential of the variable with respect to the \a j th indepenent variable. (Deprecated)  */
      const typename V::value_type& derivative(size_type j) const { return this->_dx[j]; }
      //@}
      
 #ifdef DOXYGEN
      //@{
      //! \name Comparison operators.
      /*! \brief Equality operator. */
      friend template<class X1, class V1, class X2, class V2> tribool operator==(const FirstDerivative<X1,V1>& x1, const FirstDerivative<X2,V2>& x2); 
      /*! \brief Inequality operator. */
      friend template<class X1, class V1, class X2, class V2> tribool operator!=(const FirstDerivative<X1,V1>& x1, const FirstDerivative<X2,V2>& x2); 
      //@}
     
      //@{
      //! \name Arithmetic operations
      
      /*! \brief The differential with the minimal value from \a x1 or \a x2. */
      friend FirstDerivative<X,V> min(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2);
      /*! \brief The differential with the maximal value from \a x1 or \a x2. */
      friend FirstDerivative<X,V> max(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2);
      /*! \brief Absolute value function. Returns \f$(-x,-\dot{x})\f$ if \f$x<0\f$ and \f$(x,\dot{x})\f$ if \f$x>0\f$. 
       *  If interval arithmetic is used, returns \f$(0,[-1,1]\cdot\dot{x})\f$ at zero. */
      friend FirstDerivative<X,V> abs(const FirstDerivative<X,V>& x);
      
      /*! \brief The positive value of \a x. Returns a copy \f$(x,\; \dot{x})\f$. */
      friend FirstDerivative<X,V> operator+(const FirstDerivative<X,V>& x);
      /*! \brief Negation of a differential. Returns \f$(-x,\; -\dot{x})\f$. */
      friend FirstDerivative<X,V> operator-(const FirstDerivative<X,V>& x);
      /*! \brief  Addition of differentials. Returns \f$(x_1+x_2,\; \dot{x}_1+\dot{x_2})\f$. Addition of a constant is also available. */
      friend FirstDerivative<X,V> operator+(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2);
      /*! \brief  Subtraction of differentials. Returns \f$(x_1-x_2,\; \dot{x}_1-\dot{x_2})\f$. Subtraction of/from a constant is also available. */
      friend FirstDerivative<X,V> operator-(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2);
      /*! \brief Multiplication of differentials. Returns \f$(x_1x_2,\; x_2\dot{x}_1+x_1\dot{x_2})\f$. Multiplication by a constant is also available. */
      friend FirstDerivative<X,V> operator*(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2);
      /*! \brief Division of differentials. Returns \f$(x_1/x_2,\; \dot{x}_1/x2-x_1\dot{x_2}/x_2^2)\f$. Division of/by a constant is also available. */
      friend FirstDerivative<X,V> operator/(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2);
      /*! \brief Integer power of a differential. Returns \f$(x^n,\; nx^{n-1}\dot{x})\f$. */
      friend template<class N> FirstDerivative<X,V> pow(const FirstDerivative<X,V>& x, const N& n);
      //@}
      
      //@{
      //! \name Algebraic and transcendental operations
      /*! \brief Square root. Returns \f$(\sqrt{x},\; \dot{x}/2\sqrt{x})\f$. */
      friend FirstDerivative<X,V> sqrt(const FirstDerivative<X,V>& x); 
      /*! \brief Exponential. Returns \f$(e^x,\; \dot{x}e^x)\f$. */
      friend FirstDerivative<X,V> exp(const FirstDerivative<X,V>& x); 
      /*! \brief Natural logarithm. Returns \f$(\log x,\; \dot{x}/x)\f$. */
      friend FirstDerivative<X,V> log(const FirstDerivative<X,V>& x); 
      /*! \brief Sine function. Returns \f$(\sin x,\; \dot{x}\cos x)\f$. */
      friend FirstDerivative<X,V> sin(const FirstDerivative<X,V>& x); 
      /*! \brief Cosine function. Returns \f$(\cos x,\; -\dot{x}\sin x)\f$. */
      friend FirstDerivative<X,V> cos(const FirstDerivative<X,V>& x); 
      /*! \brief Tangent function. Returns \f$(\tan x,\; \dot{x}\sec^2 x)\f$. */
      friend FirstDerivative<X,V> tan(const FirstDerivative<X,V>& x); 
      /*! \brief Inverse sine function. Returns \f$(\sin^{-1} x,\; \dot{x}/\sqrt{1-x^2})\f$. */
      friend FirstDerivative<X,V> asin(const FirstDerivative<X,V>& x); 
      /*! \brief Inverse cosine function. Returns \f$(\cos^{-1} x,\; \dot{x}/\sqrt{1-x^2})\f$. */
      friend FirstDerivative<X,V> acos(const FirstDerivative<X,V>& x); 
      /*! \brief Inverse tangent function. Returns \f$(\tan^{-1} x,\; \dot{x}/(1+x^2)\f$. */
      friend FirstDerivative<X,V> atan(const FirstDerivative<X,V>& x); 
      //@}
      
      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const FirstDerivative<X,V>& x); 
     //@}
#endif
      
    };

  template<class X1, class V1, class X2, class V2> bool operator==(const FirstDerivative<X1,V1>& x1, const FirstDerivative<X2,V2>& x2); 
  template<class X1, class V1, class X2, class V2> bool operator!=(const FirstDerivative<X1,V1>& x1, const FirstDerivative<X2,V2>& x2); 

  template<class X, class V> std::ostream& operator<<(std::ostream& os, const FirstDerivative<X,V>& x);

  template<class X, class V> FirstDerivative<X,V> operator+(const FirstDerivative<X,V>& x);
  template<class X, class V> FirstDerivative<X,V> operator-(const FirstDerivative<X,V>& x);
  template<class X, class V> FirstDerivative<X,V> operator+(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> operator-(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> operator*(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> operator/(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 

  template<class C, class X, class V> FirstDerivative<X,V> operator+(const C& c, const FirstDerivative<X,V>& x); 
  template<class C, class X, class V> FirstDerivative<X,V> operator+(const FirstDerivative<X,V>& x, const C& c); 
  template<class C, class X, class V> FirstDerivative<X,V> operator-(const C& c, const FirstDerivative<X,V>& x); 
  template<class C, class X, class V> FirstDerivative<X,V> operator-(const FirstDerivative<X,V>& x, const C& c); 
  template<class C, class X, class V> FirstDerivative<X,V> operator*(const C& c, const FirstDerivative<X,V>& x); 
  template<class C, class X, class V> FirstDerivative<X,V> operator*(const FirstDerivative<X,V>& x, const C& c); 
  template<class C, class X, class V> FirstDerivative<X,V> operator/(const C& c, const FirstDerivative<X,V>& x); 
  template<class C, class X, class V> FirstDerivative<X,V> operator/(const FirstDerivative<X,V>& x, const C& c); 

  template<class X, class V> FirstDerivative<X,V> min(const FirstDerivative<X,V>& x1,const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> max(const FirstDerivative<X,V>& x1,const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> abs(const FirstDerivative<X,V>& x); 

  template<class X, class V> FirstDerivative<X,V> pos(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> neg(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> add(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> sub(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> mul(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> div(const FirstDerivative<X,V>& x1, const FirstDerivative<X,V>& x2); 
  template<class X, class V> FirstDerivative<X,V> pow(const FirstDerivative<X,V>& x, int n); 
  template<class X, class V> FirstDerivative<X,V> sqrt(const FirstDerivative<X,V>& x);
  template<class X, class V> FirstDerivative<X,V> exp(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> log(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> sin(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> cos(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> tan(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> asin(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> acos(const FirstDerivative<X,V>& x); 
  template<class X, class V> FirstDerivative<X,V> atan(const FirstDerivative<X,V>& x); 


  }

}

#include "first_derivative.inline.h"

#endif /* ARIADNE_FIRST_DERIVATIVE_H */
