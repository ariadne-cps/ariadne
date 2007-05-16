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

#include "../base/tribool.h"
#include "../base/exceptions.h"

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
     * The constructor Differential(x,v) constructs a variable whose value is \a x and whose partial derivative with respect to the
     * \a i th independent variable is \a v[i] .
     * To construct a "constant" quantity with respect to \a n arguments, use %Differential(c, zero_vector(n)) .
     * To construct the \a i th "variable" x of \a n, use %Differential(x,unit_vector(n,i)) .
     *
     * The specialization Differential<X,X> represents a quantiy and its first variation with respect to a single independent variable.
     * The SecondDifferential<X,X,X> class represents a quantity and its first two variations with respect to an independent variable.
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
    
    
    /*!\ingroup Numeric
     * \brief A templated class representing a the differential of a quantity with respect to a scalar argument, 
     * which can be used to compute the derivative of a scalar-valued function.
     *
     * A scalar differential is a pair \f$(x,\dot{x})\f$ where \f$x\f$ is the value of the quantity,
     * and \f$\dot{x}\f$ is the variation with respect to some independent variable. 
     * A differential of the form \f$(c,0)\f$ represents a constant. 
     * A differential of the form \f$(x,1)\f$ represents the independent variable with value \f$x\f$.
     *
     * The operation of a scalar function \f$f\f$ on the differential \f$(x,\dot{x})\f$ is defined \f[f(x,\dot{x}):=(f(x),f'(x)\dot{x}) . \f]
     *
     * To compute the derivative of a scalar function \f$f\f$ at a point \f$x\f$, evaluate f(%Differential<%Float,%Float>(x,1)).%derivative().
     * To compute range of derivatives of a scalar function \f$f\f$ over an interval \f$I\f$, evaluate f(%Differential<%Interval,%Interval>(I,1)).%derivative().
     *
     * To construct a constant differential \a c, use Differential(c,0).
     * To construct the differential of a variable \a x, use Differential(x,1).
     * 
     * See the documentation for the Differential template for a list of supported operations. See the SecondDifferential<X,X,X> template for computing second derivatives.
     */
    template<class X>
    class Differential< X, X >
    {
     private:
      X _x; X _dx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs a differential based on the constant zero. */
      Differential() : _x(), _dx() { }
      /*! \brief Constuct a differential based on variable x and derivative dx. */
      template<class XX0,class XX1> Differential(const XX0& x, const XX1& dx) : _x(x), _dx(dx) { }
      /*! \brief Assign a differential from a constant quantity. */
      template<class XX> Differential<X,X>& operator=(const XX& c) { this->_x=c; this->_dx=0; return *this; }
      //@}

      //@{
      //! \name Data access
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The derivative of the variable with respect to the independent variable. */
      const X& derivative() const { return this->_dx; }
      //@}
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


    template<class X> inline
    std::ostream& operator<<(std::ostream& os, const Differential< X, Differential<X,X> >& x) {
      return os << "("<<x.value()<<","<<x.first_derivative()<<","<<x.second_derivative()<<")";
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







    template<class X,class V, class H> class SecondDifferential;

    /*!\ingroup Numeric
     * \brief A templated class representing a the second differential of a quantity with respect to a scalar argument,
     * which can be used to compute the second derivative of a scalar-valued function.
     *
     * The image of the second differential \f$(x,\dot{x},\ddot{x})\f$ under a scalar function \f$f\f$ is 
     *   \f[f(x,\dot{x},\ddot{x}):=(f(x),\;f'(x)\dot{x},\;f'(x)\ddot{x}+f''(x)\dot{x}^2) . \f]
     *
     * To compute the second derivative of a scalar function \f$f\f$ at a point \f$x\f$, evaluate f(%SecondDifferential<%Float,%Float,%FLoat>(x,1,0)).%second_derivative().
     * To compute range of second of a scalar function \f$f\f$ over an interval \f$I\f$, evaluate f(%Differential<%Interval,%Interval,%Interval>(I,1,0)).%derivative().
     *
     * To construct a constant second differential \a c, use SecondDifferential(c,0,0).
     * To construct the differential of a variable \a x, use SecondDifferential(x,1,0).
     */
    template<class X>
    class SecondDifferential< X, X, X >
    {
      private:
      X _x; X _dx; X _ddx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor. */
      SecondDifferential() : _x(), _dx(), _ddx() { }
      /*! \brief Constuct a second differential based on variable \a x, derivative \a dx and second derivative \a ddx. */
      template<class XX0,class XX1,class XX2> SecondDifferential(const XX0& x, const XX1& dx, const XX2& ddx)
        : _x(x), _dx(dx), _ddx(ddx) { }
      /*! \brief Assign a second differential from a constant quantity. */
      template<class XX> SecondDifferential<X,X,X>& operator=(const XX& c) { 
        this->_x=c; this->_dx=0; this->_ddx=0; return *this; }
      //@}

      //@{
      //! \name Data access
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The derivative of the variable with respect to the indepenent variable. */
      const X& derivative() const { return this->_dx; }
      /*! \brief The first derivative of the variable with respect to the indepenent variable. */
      const X& first_derivative() const { return this->_dx; }
      /*! \brief The second derivative of the variable with respect to the indepenent variable. */
      const X& second_derivative() const { return this->_ddx; }
      //@}
      
    };
    
      
      
    template<class X1, class V1, class H1, class X2, class V2, class H2> 
    bool operator==(const SecondDifferential<X1,V1,H1>& x1, const SecondDifferential<X2,V2,H2>& x2) {
      return x1.value()==x2.value() 
        && x1.derivative()==x2.derivative()
        && x1.second_derivative()==x2.second_derivative();
    }

    template<class X1, class V1, class H1, class X2, class V2, class H2> 
    bool operator!=(const SecondDifferential<X1,V1,H1>& x1, const SecondDifferential<X2,V2,H2>& x2) {
      return !(x1==x2);
    }

    template<class X, class V, class H> 
    std::ostream& operator<<(std::ostream& os, const SecondDifferential<X,V,H>& x) {
    return os << "("<<x.value()<<","<<x.derivative()<<","<<x.second_derivative()<<")";
    }


    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator+(const C& c, const SecondDifferential<X,V,H>& x) {
      return SecondDifferential<X,V,H>(c+x.value(),x.first_derivative(),x.second_derivative());
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator+(const SecondDifferential<X,V,H>& x, const C& c) {
      return SecondDifferential<X,V,H>(x.value()+c,x.first_derivative(),x.second_derivative());
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator-(const C& c, const SecondDifferential<X,V,H>& x) {
      return SecondDifferential<X,V,H>(c-x.value(),-x.first_derivative(),-x.second_derivative());
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator-(const SecondDifferential<X,V,H>& x, const C& c) {
      return SecondDifferential<X,V,H>(x.value()-c,x.first_derivative(),x.second_derivative());
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator*(const C& c, const SecondDifferential<X,V,H>& x) {
      return SecondDifferential<X,V,H>(c*x.value(),c*x.first_derivative(),c*x.second_derivative());
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator*(const SecondDifferential<X,V,H>& x, const C& c) {
      return SecondDifferential<X,V,H>(x.value()*c,x.first_derivative()*c,x.second_derivative()*c);
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator/(const C& c, const SecondDifferential<X,V,H>& x) {
      X r=1/x.value();
      X dr=-r*r;
      X ddr=-2*r*dr;
      //std::cerr << "y="<<y<<" s="<<s<<" dy="<<dy<<" ds="<<ds<<" ddy="<<ddy<<std::endl;
      return SecondDifferential<X,X,X>(c*r,c*dr*x.first_derivative(),
		                               c*(dr*x.second_derivative()+ddr*pow(x.first_derivative(),2)));
    }
 
    template<class X, class V, class H, class C> inline 
    SecondDifferential<X,V,H> operator/(const SecondDifferential<X,V,H>& x, const C& c) {
      return SecondDifferential<X,V,H>(x.value()/c,x.first_derivative()/c,x.second_derivative()/c);
    }



    template<class X>  
    SecondDifferential<X,X,X> abs(const SecondDifferential<X,X,X>& x) {
      if(x.value()==0) { ARIADNE_THROW(std::runtime_error,"abs(Differntial x)","x.value()==0"); }
      if(x.value()>0) { return pos(x); } else { return neg(x); }
    }
    
    template<class X>  
    SecondDifferential<X,X,X> pos(const SecondDifferential<X,X,X>& x) {
      return SecondDifferential<X,X,X>(x.value(),x.derivative(),x.second_derivative());
    }
    
    template<class X>  
    SecondDifferential<X,X,X> neg(const SecondDifferential<X,X,X>& x) {
      return SecondDifferential<X,X,X>(-x.value(),-x.derivative(),-x.second_derivative());
    }
    
    template<class X>  
    SecondDifferential<X,X,X> add(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()+x2.value();
      X dy=x1.derivative()+x2.derivative();
      X ddy=x1.second_derivative()+x2.second_derivative();
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> sub(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()-x2.value();
      X dy=x1.derivative()-x2.derivative();
      X ddy=x1.second_derivative()-x2.second_derivative();
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> mul(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()*x2.value();
      X dy=x1.derivative()*x2.value()+x1.value()*x2.derivative();
      X ddy=x1.second_derivative()*x2.value()+2*x1.derivative()*x2.derivative()+x1.value()*x2.second_derivative();
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    // (ddx1*x2-2*x2*dx1*dx2-x1*ddx2)/x2^2+2*x1*dx2^2/x2^3
    template<class X>  
    SecondDifferential<X,X,X> div(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()/x2.value();
      X s=x1.derivative()-y*x2.derivative();
      X dy=s/x2.value();
      X ds=x1.second_derivative()-dy*x2.derivative()-y*x2.second_derivative();
      X ddy=(ds-dy*x2.derivative())/x2.value();
      //std::cerr << "y="<<y<<" s="<<s<<" dy="<<dy<<" ds="<<ds<<" ddy="<<ddy<<std::endl;
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> pow(const SecondDifferential<X,X,X>& x, int n) {
      //std::cerr << "pow(SecondDifferential x, int n)"<<std::endl;
      if(n==0) {
        return SecondDifferential<X,X,X>(1,0,0);
      } else if(n==1) {
        return SecondDifferential<X,X,X>(x.value(),1,0);
      } else {
        X z2=(n*(n-1))*pow(x.value(),n-2);
        X z1=n*pow(x.value(),n-1);
        X z0=pow(x.value(),n);
        X w=pow(x.derivative(),2);
        return SecondDifferential<X,X,X>(z0,z1*x.derivative(),z1*x.second_derivative()+z2*w);
      }
    }
    
    template<class X>  
    SecondDifferential<X,X,X> sqrt(const SecondDifferential<X,X,X>& x) {
      X y=sqrt(x.value());
      X r=static_cast<X>(0.5)/y;
      X dy=x.derivative()*r;
      X ddy=(x.second_derivative()-pow(x.derivative(),2)/(2*x.value()))*r;
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> exp(const SecondDifferential<X,X,X>& x) {
      X y=exp(x.value());
      X dy=y*x.derivative();
      X ddy=y*(pow(x.derivative(),2)+x.second_derivative());
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> log(const SecondDifferential<X,X,X>& x) {
      X y=log(x.value());
      X dy=x.derivative()/x.value();
      X ddy=x.second_derivative()/x.value()-pow(dy,2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> sin(const SecondDifferential<X,X,X>& x) {
      X s=sin(x.value());
      X c=cos(x.value());
      const X& y=s;
      X dy=c*x.derivative();
      X ddy=c*x.second_derivative()-s*pow(x.derivative(),2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> cos(const SecondDifferential<X,X,X>& x) {
      X c=cos(x.value());
      X s=sin(x.value());
      const X& y=c;
      X dy=-s*x.derivative();
      X ddy=-s*x.second_derivative()-c*pow(x.derivative(),2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> tan(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      X y=tan(x.value());
      X z=1+pow(y,2);
      X dy=z*x.derivative();
      X ddy=dy*x.second_derivative()+2*y*z*pow(x.derivative(),2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> asin(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      // TODO: Check these formulae
      X r=1/sqrt(1-pow(x.value(),2));
      X y=asin(x.value());
      X dy=r*x.derivative();
      X ddy=r*(x.second_derivative()+pow(dy,2));
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }


    template<class X>  
    SecondDifferential<X,X,X> acos(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      X r=1/sqrt(1-pow(x.value(),2));
      X y=acos(x.value());
      X dy=r*x.derivative();
      X ddy=r*(x.second_derivative()+pow(dy,2));
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }

    template<class X>  
    SecondDifferential<X,X,X> atan(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      X r=1/(1+pow(x.value(),2));
      X y=atan(x.value());
      X dy=x.derivative()*r;
      X ddy=x.second_derivative()*r-2*x.value()*pow(dy,2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }


    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator+(const SecondDifferential<X,V,H>& x) {
      return pos(x);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator-(const SecondDifferential<X,V,H>& x) {
      return neg(x);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator+(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return add(x1,x2);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator-(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return sub(x1,x2);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator*(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return mul(x1,x2);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator/(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return div(x1,x2);
    }



  }
}

#endif /* ARIADNE_DIFFERENTIAL_H */
