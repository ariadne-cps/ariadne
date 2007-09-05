/***************************************************************************
 *            differential.h
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
 
/*! \file differential.h
 *  \brief Differentials for automatic differentiation.
 */
 
#ifndef ARIADNE_DIFFERENTIAL_H
#define ARIADNE_DIFFERENTIAL_H

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
  namespace Numeric {

    /*!\ingroup Numeric
     * \brief A templated class representing a quantity and its first variation or differential with respect to one or more independent variables.
     *
     * A differential \f$(x,v)\f$ or \f$(x,\dot{x})\f$ represents a quantity whose value is \f$x\f$ and whose derivative with respect to the ith independent
     * variable is \f$\partial x/\partial t_i=v_i\f$. A
     * A differential of the form \f$(c,0)\f$ represents a constant. 
     * A differential of the form \f$(x,e_i)\f$, where \f$e_i\f$ is the ith unit vector, represents the independent variable \f$x_i\f$.
     *
     * Addition of differentials is given by \f[(x,v)+(y,w):=(x+y,v+w) . \f]
     * By the Leibnitz rule, \f[(x,v)\times(y,w):=(xy,\;yv+xw) . \f]
     * The action of a scalar function \f$f\f$ on the differential \f$(x,v)\f$ is defined \f[f(x,v):=(f(x),\;f'(x)v) . \f]
     * The chain rule is automatically obtained: \f[(g\circ f)(x,v)=\bigl(g(f(x)),\;g'(f(x))f'(x)v\bigr) . \f]
     *
     * The constructor <tt>Differential(x,v)</tt> constructs a variable whose value is \f$x\f$ and whose partial derivative with respect to the
     * \f$i\f$<sup>th</sup> independent variable is <tt>v[i]</tt> .
     * To construct a "constant" quantity with respect to \f$n\f$ arguments, use <tt>%Differential(c, zero_vector(n))</tt>.
     * To construct the \f$i\f$<sup>th</sup> "variable" \f$x\f$ of \a n, use <tt>%Differential(x,unit_vector(n,i))</tt>.
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
    class Differential
    {
     private:
      X _x; V _dx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs the differential of the constant zero with respect to no independent variable. */
      Differential() : _x(), _dx() { }
      /*! \brief Constuct a differential based on variable x and derivative vector dx. */
      template<class XX,class VV> Differential(const XX& x, const VV& dx) : _x(x), _dx(dx) { }
      /*! \brief Construct the differential of a constant with respect to \a n independent variables. (Deprecated) */
      template<class XX> Differential(const XX& x, size_type n) : _x(x), _dx(n) { }
      /*! \brief Construct the differential of the \a i th variable of \a n independent variables. (Deprecated) */
      template<class XX> Differential(const XX& x, size_type n, size_type i) : _x(x), _dx(n) { _dx[i]=1; }
      /*! \brief Assign the differential of a constant zero with respect to the current independent variables (useful as an alternative to a constructor).*/
      template<class XX> Differential<X,V> operator=(const XX& c) { 
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
      friend template<class X1, class V1, class X2, class V2> tribool operator==(const Differential<X1,V1>& x1, const Differential<X2,V2>& x2); 
      /*! \brief Inequality operator. */
      friend template<class X1, class V1, class X2, class V2> tribool operator!=(const Differential<X1,V1>& x1, const Differential<X2,V2>& x2); 
      //@}
     
      //@{
      //! \name Arithmetic operations
      
      /*! \brief The differential with the minimal value from \a x1 or \a x2. */
      friend Differential<X,V> min(const Differential<X,V>& x1, const Differential<X,V>& x2);
      /*! \brief The differential with the maximal value from \a x1 or \a x2. */
      friend Differential<X,V> max(const Differential<X,V>& x1, const Differential<X,V>& x2);
      /*! \brief Absolute value function. Returns \f$(-x,-\dot{x})\f$ if \f$x<0\f$ and \f$(x,\dot{x})\f$ if \f$x>0\f$. 
       *  If interval arithmetic is used, returns \f$(0,[-1,1]\cdot\dot{x})\f$ at zero. */
      friend Differential<X,V> abs(const Differential<X,V>& x);
      
      /*! \brief The positive value of \a x. Returns a copy \f$(x,\; \dot{x})\f$. */
      friend Differential<X,V> operator+(const Differential<X,V>& x);
      /*! \brief Negation of a differential. Returns \f$(-x,\; -\dot{x})\f$. */
      friend Differential<X,V> operator-(const Differential<X,V>& x);
      /*! \brief  Addition of differentials. Returns \f$(x_1+x_2,\; \dot{x}_1+\dot{x_2})\f$. Addition of a constant is also available. */
      friend Differential<X,V> operator+(const Differential<X,V>& x1, const Differential<X,V>& x2);
      /*! \brief  Subtraction of differentials. Returns \f$(x_1-x_2,\; \dot{x}_1-\dot{x_2})\f$. Subtraction of/from a constant is also available. */
      friend Differential<X,V> operator-(const Differential<X,V>& x1, const Differential<X,V>& x2);
      /*! \brief Multiplication of differentials. Returns \f$(x_1x_2,\; x_2\dot{x}_1+x_1\dot{x_2})\f$. Multiplication by a constant is also available. */
      friend Differential<X,V> operator*(const Differential<X,V>& x1, const Differential<X,V>& x2);
      /*! \brief Division of differentials. Returns \f$(x_1/x_2,\; \dot{x}_1/x2-x_1\dot{x_2}/x_2^2)\f$. Division of/by a constant is also available. */
      friend Differential<X,V> operator/(const Differential<X,V>& x1, const Differential<X,V>& x2);
      /*! \brief Integer power of a differential. Returns \f$(x^n,\; nx^{n-1}\dot{x})\f$. */
      friend template<class N> Differential<X,V> pow(const Differential<X,V>& x, const N& n);
      //@}
      
      //@{
      //! \name Algebraic and transcendental operations
      /*! \brief Square root. Returns \f$(\sqrt{x},\; \dot{x}/2\sqrt{x})\f$. */
      friend Differential<X,V> sqrt(const Differential<X,V>& x); 
      /*! \brief Exponential. Returns \f$(e^x,\; \dot{x}e^x)\f$. */
      friend Differential<X,V> exp(const Differential<X,V>& x); 
      /*! \brief Natural logarithm. Returns \f$(\log x,\; \dot{x}/x)\f$. */
      friend Differential<X,V> log(const Differential<X,V>& x); 
      /*! \brief Sine function. Returns \f$(\sin x,\; \dot{x}\cos x)\f$. */
      friend Differential<X,V> sin(const Differential<X,V>& x); 
      /*! \brief Cosine function. Returns \f$(\cos x,\; -\dot{x}\sin x)\f$. */
      friend Differential<X,V> cos(const Differential<X,V>& x); 
      /*! \brief Tangent function. Returns \f$(\tan x,\; \dot{x}\sec^2 x)\f$. */
      friend Differential<X,V> tan(const Differential<X,V>& x); 
      /*! \brief Inverse sine function. Returns \f$(\sin^{-1} x,\; \dot{x}/\sqrt{1-x^2})\f$. */
      friend Differential<X,V> asin(const Differential<X,V>& x); 
      /*! \brief Inverse cosine function. Returns \f$(\cos^{-1} x,\; \dot{x}/\sqrt{1-x^2})\f$. */
      friend Differential<X,V> acos(const Differential<X,V>& x); 
      /*! \brief Inverse tangent function. Returns \f$(\tan^{-1} x,\; \dot{x}/(1+x^2)\f$. */
      friend Differential<X,V> atan(const Differential<X,V>& x); 
      //@}
      
      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const Differential<X,V>& x); 
     //@}
#endif
      
    };
    
    template<class X1, class V1, class X2, class V2> inline
    bool operator==(const Differential<X1,V1>& x1, const Differential<X2,V2>& x2) {
      return x1.value()==x2.value() && x1.derivative()==x2.derivative();
    }

    template<class X1, class V1, class X2, class V2> inline
    bool operator!=(const Differential<X1,V1>& x1, const Differential<X2,V2>& x2) {
      return !(x1==x2);
    }

    template<class X, class V> inline
    std::ostream& operator<<(std::ostream& os, const Differential<X,V>& x) {
      return os << "("<<x.value()<<","<<x.derivative()<<")";
    }


    template<class C, class X, class V> inline 
    Differential<X,V> operator+(const C& c, const Differential<X,V>& x) {
      return Differential<X,V>(c+x.value(),x.derivative());
    }

    template<class C, class X, class V> inline 
    Differential<X,V> operator+(const Differential<X,V>& x, const C& c) {
      return Differential<X,V>(x.value()+c,x.derivative());
    }

    template<class C, class X, class V> inline 
    Differential<X,V> operator-(const C& c, const Differential<X,V>& x) {
      return Differential<X,V>(c-x.value(),-x.derivative());
    }

    template<class C, class X, class V> inline 
    Differential<X,V> operator-(const Differential<X,V>& x, const C& c) {
      return Differential<X,V>(x.value()-c,x.derivative());
    }

    template<class C, class X, class V> inline 
    Differential<X,V> operator*(const C& c, const Differential<X,V>& x) {
      return Differential<X,V>(c*x.value(),c*x.derivative());
    }
    
    template<class C, class X, class V> inline 
    Differential<X,V> operator*(const Differential<X,V>& x, const C& c) {
      return Differential<X,V>(c*x.value(),c*x.derivative());
    }
    
    template<class C, class X, class V> inline 
    Differential<X,V> operator/(const C& c, const Differential<X,V>& x) {
      // Use this form to get right dimension of constant.
	  return (c+0*x)/x;
    }
    
    template<class C, class X, class V> inline 
    Differential<X,V> operator/(const Differential<X,V>& x, const C& c) {
      return Differential<X,V>(x.value()/c,x.derivative()/c);
    }
    


    template<class X, class V> inline 
    Differential<X,V> min(const Differential<X,V>& x1,const Differential<X,V>& x2) {
      if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(Differntial x1, Differential x2)","x1.value()==x2.value()");
      }
      return x1.value()<x2.value() ? x1 : x2;
    }
    
    template<class X, class V> inline 
    Differential<X,V> max(const Differential<X,V>& x1,const Differential<X,V>& x2) {
      if(x1.value()==x2.value()) { 
        ARIADNE_THROW(std::runtime_error,"max(Differntial x1, Differential x2)","x1.value()==x2.value()"); 
      }
      return x1.value()>x2.value() ? x1 : x2;
    }
    
    template<class X, class V> inline 
    Differential<X,V> abs(const Differential<X,V>& x) {
      if(x.value()==0) { ARIADNE_THROW(std::runtime_error,"abs(Differntial x)","x.value()==0"); }
      if(x.value()>0) { return x; } else { return -x; }
    }
    
    template<class X, class V> inline 
    Differential<X,V> pos(const Differential<X,V>& x) {
      return Differential<X,V>(x.value(),x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> neg(const Differential<X,V>& x) {
      return Differential<X,V>(-x.value(),-x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> add(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return Differential<X,V>(x1.value()+x2.value(),x1.derivative()+x2.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> sub(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return Differential<X,V>(x1.value()-x2.value(),x1.derivative()-x2.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> mul(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return Differential<X,V>(x1.value()*x2.value(),x1.derivative()*x2.value()+x1.value()*x2.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> div(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      X y=x1.value()/x2.value();
      V v=x1.derivative()-y*x2.derivative();
      return Differential<X,V>(y,v/x2.value());
    }
    
    template<class X, class V> inline 
    Differential<X,V> pow(const Differential<X,V>& x, int n) {
      //std::cerr << "pow(Differential x, int n)"<<std::endl;
      if(n==0) {
        return Differential<X,V>(pow(x.value(),0),static_cast<X>(0)*x.derivative());
      } else {
        X y=pow(x.value(),n-1);
        return Differential<X,V>((x.value()*y),static_cast<X>(n*y)*x.derivative());
      }
    }
    
    template<class X, class V> inline 
    Differential<X,V> sqrt(const Differential<X,V>& x) {
      X y=sqrt(x.value());
      return Differential<X,V>(y,x.derivative()/(2*y));
    }
    
    template<class X, class V> inline 
    Differential<X,V> exp(const Differential<X,V>& x) {
      X y=exp(x.value());
      return Differential<X,V>(y,y*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> log(const Differential<X,V>& x) {
      return Differential<X,V>(log(x.value()),x.derivative()/x.value());
    }
    
    template<class X, class V> inline 
    Differential<X,V> sin(const Differential<X,V>& x) {
      return Differential<X,V>(sin(x.value()),cos(x.value())*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> cos(const Differential<X,V>& x) {
      return Differential<X,V>(cos(x.value()),-sin(x.value())*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> tan(const Differential<X,V>& x) {
      X y=tan(x.value());
      return Differential<X,V>(y,(X(1)-y*y)*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> asin(const Differential<X,V>& x) {
      return Differential<X,V>(asin(x.value()),x.derivative()/sqrt(X(1)-x.value()*x.value()));
    }
    
    template<class X, class V> inline 
    Differential<X,V> acos(const Differential<X,V>& x) {
      return Differential<X,V>(acos(x.value()),-x.derivative()/sqrt(X(1)-x.value()*x.value()));
    }
    
    template<class X, class V> inline 
    Differential<X,V> atan(const Differential<X,V>& x) {
      return Differential<X,V>(atan(x.value()),x.derivative()/(X(1)+x.value()*x.value()));
    }
 


    template<class X, class V> inline 
    Differential<X,V> operator+(const Differential<X,V>& x) {
      return pos(x);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator-(const Differential<X,V>& x) {
      return neg(x);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator+(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return add(x1,x2);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator-(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return sub(x1,x2);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator*(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return mul(x1,x2);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator/(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return div(x1,x2);
    }

  }
}

#endif /* ARIADNE_DIFFERENTIAL_H */
