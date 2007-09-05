/***************************************************************************
 *            scalar_differential.h
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
 
/*! \file scalar_differential.h
 *  \brief Differentials of scalar functions of a single variable.
 */
 
#ifndef ARIADNE_SCALAR_DIFFERENTIAL_H
#define ARIADNE_SCALAR_DIFFERENTIAL_H

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
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a scalar argument.
     *
     * A scalar derivative is an array \f$(x,\dot{x},\ddot{x},\ldots)\f$ where \f$x\f$ is the value of the quantity,
     * and \f$\dot{x}\f$ is the variation with respect to some independent variable, \f$\ddot{x}\f$ is the second variation etc. 
     * A scalar derivative of the form \f$(c,0,0,\ldots)\f$ represents a constant. 
     * A differential of the form \f$(x,1,0,0,\ldots)\f$ represents the independent variable with value \f$x\f$.
     *
     * The operation of a scalar function \f$f\f$ on the differential \f$(x,\dot{x})\f$ is defined \f[f(x,\dot{x},\ddot{x}):=(f(x),f'(x)\dot{x},f''(x)\dot{x}^2+f'(x)\ddot{x},\ldots) . \f]
     *
     * To compute the \f$n^\mathrm{th}\f$ derivative of a scalar function \f$f\f$ at a point \f$x\f$, evaluate <tt>f(%Differential<%Float>::variable(n,x)).%derivative()</tt>.
     * To compute range of derivatives of a scalar function \f$f\f$ over an interval \f$I\f$, evaluate <tt>f(%Differential<%Interval>(I,1)).%derivative()</tt>.
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
        : _values(1u) { }
      /*! \brief The constant zero of degree \a degree. */
      ScalarDerivative(uint degree)
        : _values(degree+1u) { }
      /*! \brief The constant \a constant of degree \a degree. */
      template<class X0> ScalarDerivative(uint degree, const X0& constant)
        : _values(degree+1u) { this->_values[0]=constant; }
      /*! \brief A scalar derivative of degree \a degree, with value \a value and first derivative \a first_derivative. Higher derivatives are set to zero. */
      template<class X0, class X1> ScalarDerivative(uint degree, const X0& value, const X1& first_derivative)
        : _values(degree+1u) 
      { this->_values[0]=value; this->_values[1]=first_derivative; }

      /*! \brief A scalar derivative with values given by \a ary. */
      template<class XX> ScalarDerivative(const array<XX>& ary)
        : _values(ary) { }
      /*! \brief A scalar derivative with values given by array literal \a str. */
      ScalarDerivative(const std::string& str) { 
        std::stringstream ss(str); read_array(ss,_values); }
    
      /*! \brief Copy constructor. */
      template<class XX> ScalarDerivative(const ScalarDerivative<XX>& other) 
        : _values(other._values) { }
      /*! \brief Copy assignment operator. */
      template<class XX> ScalarDerivative<X>& operator=(const ScalarDerivative<XX>& other) {
        this->_values=other._values; return *this; }

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const ScalarDerivative<XX>& other) {
        return this->_values==other._values; }
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const ScalarDerivative<XX>& other) {
        return !(*this==other); }

      /*! \brief Construct a constant derivative of degree \a degree and value \a constant. */
      static ScalarDerivative<X> constant(uint degree, const X& constant) { 
        return ScalarDerivative<X>(degree,constant); }
      /*! \brief Construct the derivative of degree \a degree for the indepentent variable at value \a value. */
      static ScalarDerivative<X> variable(uint degree, const X& value) {
        return ScalarDerivative<X>(degree,value,1); }

      /*! \brief The degree (number of derivatives computed). */
      uint degree() const { 
        return this->_values.size()-1; }
      /*! \brief The array of derivative values. */
      const array<X>& values() const {
        return this->_values; }
      /*! \brief The \a i<sup> th</sup> derivative \f$d^if/dx^i\f$. */
      const X& derivative(uint i) const { 
        return this->_values[i]; }
      /*! \brief A reference to the \a i<sup> th</sup> derivative \f$d^if/dx^i\f$. */
      X& operator[](uint i) { 
        return this->_values[i]; }
      /*! \brief The \a i<sup> th</sup> derivative \f$d^if/dx^i\f$. */
      const X& operator[](uint i) const { 
        return this->_values[i]; }
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
      array<X> _values;
    };


    /* The composition of two derivatives, computed in-place. 
     *
     *  The composition inductively by
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
    ScalarDerivative<X> compose(const ScalarDerivative<X>& y, const ScalarDerivative<X>& x)
    {
      ScalarDerivative<X> result(std::min(x.degree(),y.degree()));
      for(uint n=0; n!=result.degree(); ++n) { result[n]=y[n]; }
      compute_composition(result,x);
      return result;
    }


    template<class X> inline
    ScalarDerivative<X> inverse(const ScalarDerivative<X>& y, const X& x)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
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
      ScalarDerivative<X> y(x.degree());
      y[0]=sqrt(x[0]);
      X mhr=(-0.5)/x[0];
      for(uint i=1; i<=y.degree(); ++i) {
        y[i]=(2*int(i)-3)*mhr*y[i-1];
      }
      compute_composition(y,x);
      return y;
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

#endif /* ARIADNE_SCALAR_DIFFERENTIAL_H */
