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





    


    /*!\ingroup Numeric
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a scalar argument.
     *
     * A scalar derivative is an array \f$(x,\dot{x},\ddot{x},\ldots)\f$ where \f$x\f$ is the value of the quantity,
     * and \f$\dot{x}\f$ is the variation with respect to some independent variable, \f$\ddot{x}\f$ is the second variation etc. 
     * A scalar derivative of the form \f$(c,0,0,\ldots)\f$ represents a constant. 
     * A differential of the form \f$(x,1,0,0,\ldots)\f$ represents the independent variable with value \f$x\f$.
     *
     * The operation of a scalar function \f$f\f$ on the differential \f$(x,\dot{x})\f$ is defined \f[f(x,\dot{x},\ddot{x}):=(f(x),f'(x)\dot{x},f''(x)\dot{x}^2+f'(x)\ddot{x},\ldots) . \f]
     *
     * To compute the \f$n^\th\f$ derivative of a scalar function \f$f\f$ at a point \f$x\f$, evaluate f(%Differential<%Float>::variable(n,x)).%derivative().
     * To compute range of derivatives of a scalar function \f$f\f$ over an interval \f$I\f$, evaluate f(%Differential<%Interval>(I,1)).%derivative().
     *
     * To construct a constant \a c, use ScalarDerivative(n,c) or ScalarDerivative::constant(n,c).
     * To construct the derivative of a variable \a x, use ScalarDerivative(n,x,1) or ScalarDerivative::variable(n,x)..
     * 
     */
    template<class X>
    class ScalarDerivative
    {
     public:
      typedef X value_type;

      ScalarDerivative()
        : _values(1u) { }
      ScalarDerivative(uint degree)
        : _values(degree+1u) { }
      template<class X0> ScalarDerivative(uint degree, const X0& constant)
        : _values(degree+1u) { this->_values[0]=constant; }
      template<class X0, class X1> ScalarDerivative(uint degree, const X0& constant, const X1& first_derivative)
        : _values(degree+1u) 
      { this->_values[0]=constant; this->_values[1]=first_derivative; }

      template<class XX> ScalarDerivative(const array<XX>& ary)
        : _values(ary) { }
      ScalarDerivative(const std::string& str) { 
        std::stringstream ss(str); read_array(ss,_values); }
    
      template<class XX> ScalarDerivative(const ScalarDerivative<XX>& other) 
        : _values(other._values) { }
      template<class XX> ScalarDerivative<X>& operator=(const ScalarDerivative<XX>& other) {
        this->_values=other._values; return *this; }

      template<class XX> bool operator==(const ScalarDerivative<XX>& other) {
        return this->_values==other._values; }
      template<class XX> bool operator!=(const ScalarDerivative<XX>& other) {
        return !(*this==other); }

      static ScalarDerivative<X> constant(uint degree, const X& constant) { 
        return ScalarDerivative<X>(degree,constant); }
      static ScalarDerivative<X> variable(uint degree, const X& value) {
        return ScalarDerivative<X>(degree,value,1); }

      uint degree() const { 
        return this->_values.size()-1; }
      const array<X>& values() const {
        return this->_values; }
      const X& derivative(uint i) const { 
        return this->_values[i]; }
      X& operator[](uint i) { 
        return this->_values[i]; }
      const X& operator[](uint i) const { 
        return this->_values[i]; }
     private:
      array<X> _values;
    };


    /*! The composition inductively by
     *  \f[ y^\[n\] = \sum_{i=0}^{n-1} \choose{n}{i} {\dot{y}}^{[i]} x^{(n-i)}
     */
    template<class X> inline
    void compute_composition(ScalarDerivative<X>& y, const ScalarDerivative<X>& x)
    {
      assert(y.degree()==x.degree());
      int d=y.degree();
      for(int n=0; n<d; ++n) {
        for(int i=0; i<=n; ++i) {
          //std::cout<<"y["<<d-i<<"] = 1*y["<<d-i<<"]*x[1]"<<std::flush;
          y[d-i]*=x[1];
          for(int j=1; j<=n-i; ++j) {
            //std::cout<<" + "<<choose(n-i,j)<<"*y["<<d-i-j<<"]*x["<<j+1<<"]"<<std::flush;
            y[d-i] += choose<int>(n-i,j) * y[d-i-j] * x[j+1];
          }
          //std::cout << std::endl;
        }
      }
      return;
    }


    template<class X> inline
    void compose(ScalarDerivative<X>& y, const ScalarDerivative<X>& x)
    {
      ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
      for(uint n=0; n!=result.degree(); ++n) { result[n]=y[n]; }
      compute_composition(result,x);
      return result;
    }


    template<class X> inline 
    ScalarDerivative<X> 
    min(const ScalarDerivative<X>& x1, const ScalarDerivative<X>& x2) 
    {
      if(x1[0]==x2[0]) {
        ARIADNE_THROW(std::runtime_error,"min(ScalarDerivative x1, ScalarDerivative x2)","x1[0]==x2[0]");
      }
      return x1[0]<x2[0] ? x1 : x2;
    }
    
    template<class X> inline 
    ScalarDerivative<X> 
    max(const ScalarDerivative<X>& x1,const ScalarDerivative<X>& x2) 
    {
      if(x1[0]==x2[0]) { 
        ARIADNE_THROW(std::runtime_error,"max(ScalarDerivative x1, ScalarDerivative x2)","x1[0]==x2[0]"); 
      }
      return x1[0]>x2[0] ? x1 : x2;
    }
    
    template<class X> inline
    ScalarDerivative<X> 
    pos(const ScalarDerivative<X>& x)
    {
      return x;
    }

    template<class X> inline
    ScalarDerivative<X> 
    neg(const ScalarDerivative<X>& x)
    {
      ScalarDerivative<X> result(x.degree());
      for(uint n=0; n<=result.degree(); ++n) {
          result[n] = -x[n];
      }
      return result;
    }

    template<class X>  
    ScalarDerivative<X> 
    abs(const ScalarDerivative<X>& x) 
    {
      if(x[0]==0) { 
        ARIADNE_THROW(std::runtime_error,"abs(ScalarDerivative x)","x[0]==0"); 
      }
      return x[0]>0 ? pos(x) : neg(x); 
    }
    
    template<class X> inline
    ScalarDerivative<X> 
    inv(const ScalarDerivative<X>& x)
    {
      ScalarDerivative<X> y(x.degree());
      X mr = X(-1)/x[0];
      for(uint i=0; i<=y.degree(); ++i) {
        y[i]=(-factorial<int>(i))*pow(mr,i+1);
      }
      //std::cerr << y << std::endl;
      compute_composition(y,x);
      //std::cerr << y << std::endl;
      return y;
    }

    template<class X> inline
    ScalarDerivative<X> 
    add(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
      for(uint n=0; n<=result.degree(); ++n) {
          result[n] = x[n]+y[n];
      }
      return result;
    }

    template<class X> inline
    ScalarDerivative<X> 
    sub(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
      for(uint n=0; n<=result.degree(); ++n) {
          result[n] = x[n]-y[n];
      }
      return result;
    }

    template<class X> inline
    ScalarDerivative<X> 
    mul(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
      for(uint n=0; n<=result.degree(); ++n) {
        for(uint i=0; i<=n; ++i) {
          result[n] += choose<int>(n,i)*x[i]*y[n-i];
        }
      }
      return result;
    }

    template<class X> inline
    ScalarDerivative<X> 
    div(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      return x*inv(y);
    }

    template<class X, class N> inline
    ScalarDerivative<X> 
    pow(const ScalarDerivative<X>& x, N k)
    {
      uint n=k;
      ScalarDerivative<X> result(x.degree());
      for(uint i=0; i<=std::min(result.degree(),n); ++i) {
        int j=n-i;
        result[i]=(factorial<int>(n)/factorial<int>(j))*pow(x[0],j);
      }
      compute_composition(result,x);
      return result;
    }

    template<class X>  
    ScalarDerivative<X> sqrt(const ScalarDerivative<X>& x) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class X>  
    ScalarDerivative<X> exp(const ScalarDerivative<X>& x) {
      ScalarDerivative<X> y(x.degree());
      y[0]=exp(x[0]);
      for(uint i=1; i<=y.degree(); ++i) {
        y[i]=y[0];
      }
      compute_composition(y,x);
      return y;
    }
    
    template<class X>  
    ScalarDerivative<X> log(const ScalarDerivative<X>& x) {
      ScalarDerivative<X> y(x.degree());
      y[0]=log(x[0]);
      X mr=(-1)/x[0];
      for(uint i=1; i!=y.degree();++i) {
        y[i]=(-factorial<int>(i-1))*pow(x[0],i);
      }
      compute_composition(y,x);
      return y;
    }
    
    template<class X>  
    ScalarDerivative<X> sin(const ScalarDerivative<X>& x) {
      uint d=x.degree();
      ScalarDerivative<X> y(d);
      y[0]=sin(x[0]);
      y[1]=cos(x[0]);
      for(uint i=2; i!=d; ++i) {
        y[i]=-y[i-2];
      }
      compute_composition(y,x);
      return y;
    }
    
    template<class X>  
    ScalarDerivative<X> cos(const ScalarDerivative<X>& x) {
      uint d=x.degree();
      ScalarDerivative<X> y(d);
      y[0]=cos(x[0]);
      y[1]=-sin(x[0]);
      for(uint i=2; i!=d; ++i) {
        y[i]=-y[i-2];
      }
      compute_composition(y,x);
      return y;
    }
    
    template<class X>  
    ScalarDerivative<X> tan(const ScalarDerivative<X>& x) {
      return sin(x)/cos(x);
    }
    
    template<class X>  
    ScalarDerivative<X> asin(const ScalarDerivative<X>& x) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class X>  
    ScalarDerivative<X> acos(const ScalarDerivative<X>& x) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class X>  
    ScalarDerivative<X> atan(const ScalarDerivative<X>& x) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }


    template<class X> inline
    ScalarDerivative<X> 
    operator-(const ScalarDerivative<X>& x)
    {
      return neg(x);
    }

    template<class X> inline
    ScalarDerivative<X> 
    operator+(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      return add(x,y);
    }

    template<class X> inline
    ScalarDerivative<X> 
    operator-(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      return sub(x,y);
    }

    template<class X> inline
    ScalarDerivative<X> 
    operator*(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      return mul(x,y);
    }

    template<class X> inline
    ScalarDerivative<X> 
    operator/(const ScalarDerivative<X>& x, const ScalarDerivative<X>& y)
    {
      return div(x,y);
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator+(const ScalarDerivative<X>& x, const R& c)
    {
      return add(x,ScalarDerivative<X>::constant(x.degree(),c));
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator+(const R& c, const ScalarDerivative<X>& x)
    {
      return add(ScalarDerivative<X>::constant(x.degree(),c),x);
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator-(const ScalarDerivative<X>& x, const R& c)
    {
      return sub(x,ScalarDerivative<X>::constant(x.degree(),c));
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator-(const R& c, const ScalarDerivative<X>& x)
    {
      return sub(ScalarDerivative<X>::constant(x.degree(),c),x);
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator*(const ScalarDerivative<X>& x, const R& c)
    {
      return mul(x,ScalarDerivative<X>::constant(x.degree(),c));
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator*(const R& c, const ScalarDerivative<X>& x)
    {
      return mul(ScalarDerivative<X>::constant(x.degree(),c),x);
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator/(const ScalarDerivative<X>& x, const R& c)
    {
      return div(x,ScalarDerivative<X>::constant(x.degree(),c));
    }

    template<class X, class R> inline
    ScalarDerivative<X> 
    operator/(const R& c, const ScalarDerivative<X>& x)
    {
      return div(ScalarDerivative<X>::constant(x.degree(),c),x);
    }

    template<class X> inline
    std::ostream& operator<<(std::ostream& os, const ScalarDerivative<X>& x) {
      os << "ScalarDerivative";
      for(uint i=0; i<=x.degree(); ++i) {
        os << (i==0 ? '(' : ',') << x[i]; 
      }
      os << ")";
      return os;
    }


  }

}

#endif /* ARIADNE_DIFFERENTIAL_H */
