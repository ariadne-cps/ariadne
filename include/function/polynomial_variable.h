/***************************************************************************
 *            polynomial_variable.h
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
 
/*! \file polynomial_variable.h
 *  \brief A scalar quantity represented as a polynomial in several variables.
 */
 
#ifndef ARIADNE_POLYNOMIAL_VARIABLE_H
#define ARIADNE_POLYNOMIAL_VARIABLE_H

#include <iostream>
#include <stdexcept>
#include <cassert>

#include "base/tribool.h"
#include "base/exceptions.h"
#include "base/stlio.h"

namespace Ariadne {
  
    class MultiIndex;
  
    /*!\ingroup Function
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a multiple arguments.
     */
    template<class X>
    class PolynomialVariable
    {
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      PolynomialVariable();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      PolynomialVariable(size_type a, smoothness_type d);
      /*! \brief A polynomial variable of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> PolynomialVariable(size_type a, smoothness_type d, const XX* ptr);

      /*! \brief Copy constructor. */
      PolynomialVariable(const PolynomialVariable<X>& tv); 
      /*! \brief Copy assignment operator. */
      PolynomialVariable<X>& operator=(const PolynomialVariable<X>& tv);

      /*! \brief Assign a constant \a c. */
      template<class XX> PolynomialVariable<X>& operator=(const XX& c);

      /*! \brief Equality operator. */
      bool operator==(const PolynomialVariable<X>& other);
      /*! \brief Inequality operator. */
      bool operator!=(const PolynomialVariable<X>& other);

      /*! \brief Construct a constant variable of degree \a d with respect to \a as variables and value \a c. */
      template<class XX> static PolynomialVariable<X> constant(size_type as, smoothness_type d, const XX& c); 
      /*! \brief Construct the variable of degree \a d at value \a value with respect to the \a i<sup>th</sup> variable of \a as. */
      template<class XX> static PolynomialVariable<X> variable(size_type as, smoothness_type d, const XX& value, size_type i);

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


      /*! \brief Construct the Taylor series of the reciprocal function. */
      static PolynomialVariable<X> rec(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the \a n<sup>th</sup> power function. */
      static PolynomialVariable<X> pow(smoothness_type d, const X& c, const uint& k); 
      /*! \brief Construct the Taylor series of the square-root function. */
      static PolynomialVariable<X> sqrt(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the exponential function. */
      static PolynomialVariable<X> exp(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the logarithm function. */
      static PolynomialVariable<X> log(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the sine function. */
      static PolynomialVariable<X> sin(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the cosine function. */
      static PolynomialVariable<X> cos(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the tangent function. */
      static PolynomialVariable<X> tan(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inverse sine function. */
      static PolynomialVariable<X> asin(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inversecosine function. */
      static PolynomialVariable<X> acos(smoothness_type d, const X& c); 
      /*! \brief Construct the Taylor series of the inversetangent function. */
      static PolynomialVariable<X> atan(smoothness_type d, const X& c); 

#ifdef DOXYGEN
    //@{ 
    //! \name Friend operations
    /*! \brief The composition of two variables computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
     *  The composition inductively by
     *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
     */
    friend PolynomialVariable<X> compose(const PolynomialVariable<X>& y, const PolynomialVariable<X>& x);
    /*! \brief The derivative of the inverse of \f$y\f$ evaluated at \f$x\f$. (Not currently implemented.) */
    friend PolynomialVariable<X> inverse(const PolynomialVariable<X>& y, const X& x);
    /*! \brief The minimum of two variables. Returns the variable whose zero-th order value is minimal. */
    friend PolynomialVariable<X> min(const PolynomialVariable<X>& x1, const PolynomialVariable<X>& x2);
    /*! \brief The maximum of two variables. Returns the variable whose zero-th order value is maximal. */
    friend PolynomialVariable<X> max(const PolynomialVariable<X>& x1, const PolynomialVariable<X>& x2);
    /*! \brief The derivatives of \f$+x\f$. Returns a copy. */
    friend PolynomialVariable<X> pos(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$-x\f$. */
    friend PolynomialVariable<X> neg(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$x+y\f$. */
    friend PolynomialVariable<X> add(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend PolynomialVariable<X> sub(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x*y\f$. */
    friend PolynomialVariable<X> mul(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x/y\f$. */
    friend PolynomialVariable<X> div(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x^n\f$. */
    friend PolynomialVariable<X> pow(const PolynomialVariable<X>& x, const Integer& n);
    /*! \brief The derivatives of \f$\sqrt{x}\f$. */
    friend PolynomialVariable<X> sqrt(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\exp(x)\f$. */
    friend PolynomialVariable<X> exp(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\log(x)\f$. */
    friend PolynomialVariable<X> log(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\sin(x)\f$. */
    friend PolynomialVariable<X> sin(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\cos(x)\f$. */
    friend PolynomialVariable<X> cos(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\tan(x)\f$. */
    friend PolynomialVariable<X> tan(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\sin^{-1}(x)\f$. (Not currently implemented.) */
    friend PolynomialVariable<X> asin(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\cos^{-1}(x)\f$. (Not currently implemented.) */
    friend PolynomialVariable<X> acos(const PolynomialVariable<X>& x);
    /*! \brief The derivatives of \f$\tan^{-1}(x)\f$. (Not currently implemented.) */
    friend PolynomialVariable<X> atan(const PolynomialVariable<X>& x);

    /*! \brief The derivatives of \f$x+y\f$. */
    friend PolynomialVariable<X> operator+(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend PolynomialVariable<X> operator-(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x*y\f$. */
    friend PolynomialVariable<X> operator*(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
    /*! \brief The derivatives of \f$x/y\f$. */
    friend PolynomialVariable<X> operator/(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);

    /*! \brief The derivatives of \f$c+x\f$ for a constant \f$c\f$. (Other mixed-mode arithmetic is also supported.) */
    friend PolynomialVariable<X> operator+(const R& c, const PolynomialVariable<X>& x);

    /*! \brief Stream output operator. */
    friend std::ostream& operator<<(std::ostream& os, const PolynomialVariable<X>& x);
    //@}
#endif 
     private:
      size_type _argument_size;
      smoothness_type _degree;
      array<X> _data;
    };


  template<class X> PolynomialVariable<X> compose(const PolynomialVariable<X>& y, const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> reduce(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> derivative(const PolynomialVariable<X>& x, const size_type& k);
  template<class X> PolynomialVariable<X> min(const PolynomialVariable<X>& x1, const PolynomialVariable<X>& x2); 
  template<class X> PolynomialVariable<X> max(const PolynomialVariable<X>& x1,const PolynomialVariable<X>& x2); 
  template<class X> PolynomialVariable<X> pos(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> neg(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> abs(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> rec(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> add(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X> PolynomialVariable<X> sub(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X> PolynomialVariable<X> mul(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X> PolynomialVariable<X> div(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X, class N> PolynomialVariable<X> pow(const PolynomialVariable<X>& x, N k);

  template<class X> PolynomialVariable<X> sqrt(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> exp(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> log(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> sin(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> cos(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> tan(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> asin(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> acos(const PolynomialVariable<X>& x); 
  template<class X> PolynomialVariable<X> atan(const PolynomialVariable<X>& x); 

  template<class X> PolynomialVariable<X> operator-(const PolynomialVariable<X>& x);
  template<class X> PolynomialVariable<X> operator+(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X> PolynomialVariable<X> operator-(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X> PolynomialVariable<X> operator*(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);
  template<class X> PolynomialVariable<X> operator/(const PolynomialVariable<X>& x, const PolynomialVariable<X>& y);

  template<class X, class R> PolynomialVariable<X> operator+(const PolynomialVariable<X>& x, const R& c);
  template<class X, class R> PolynomialVariable<X> operator+(const R& c, const PolynomialVariable<X>& x);
  template<class X, class R> PolynomialVariable<X> operator-(const PolynomialVariable<X>& x, const R& c);
  template<class X, class R> PolynomialVariable<X> operator-(const R& c, const PolynomialVariable<X>& x);
  template<class X, class R> PolynomialVariable<X> operator*(const PolynomialVariable<X>& x, const R& c);
  template<class X, class R> PolynomialVariable<X> operator*(const R& c, const PolynomialVariable<X>& x);
  template<class X, class R> PolynomialVariable<X> operator/(const PolynomialVariable<X>& x, const R& c);
  template<class X, class R> PolynomialVariable<X> operator/(const R& c, const PolynomialVariable<X>& x);

  template<class X, class R> PolynomialVariable<X> operator/=(const PolynomialVariable<X>& x, const R& c);

  template<class X> std::ostream& operator<<(std::ostream& os, const PolynomialVariable<X>& x);


  }
}


#include "linear_algebra/vector.h"

namespace Ariadne {

inline size_type compute_polynomial_data_size(size_type as, smoothness_type d) {
  return bin(d+as,as);
}

template<class X> template<class XX> 
PolynomialVariable<X>::PolynomialVariable(size_type a, smoothness_type d, const XX* ptr)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d)) 
{
  for(size_type i=0; i!=this->_data.size(); ++i) {
    this->_data[i]=ptr[i];
  }
}

template<class X> template<class XX> 
PolynomialVariable<X>& 
PolynomialVariable<X>::operator=(const XX& c) 
{
  this->_data[0]=c;
  for(size_type i=1; i!=this->_data.size(); ++i) {
    this->_data[i]=0;
  }
  return *this;
}

template<class X> template<class XX> 
PolynomialVariable<X> 
PolynomialVariable<X>::constant(size_type a, smoothness_type d, const XX& c) 
{
  PolynomialVariable<X> result(a,d);
  result._data[0]=c; 
  return result;
}

template<class X> template<class XX> 
PolynomialVariable<X> 
PolynomialVariable<X>::variable(size_type a, smoothness_type d, const XX& x, size_type i) 
{
  PolynomialVariable<X> result(a,d);
  result._data[0]=x; 
  result._data[i+1u]=1; 
  return result;
}




template<class X, class R> inline
PolynomialVariable<X> 
operator+(const R& c, const PolynomialVariable<X>& x)
{
  PolynomialVariable<X> r=x; r.value()+=c; return r;
}

template<class X, class R> inline
PolynomialVariable<X> 
operator+(const PolynomialVariable<X>& x, const R& c)
{
  PolynomialVariable<X> r=x; r.value()+=c; return r;
}

template<class X, class R> inline
PolynomialVariable<X> 
operator-(const R& c, const PolynomialVariable<X>& x)
{
  PolynomialVariable<X> r=x; 
  reinterpret_cast<Vector<X>&>(r.data())*=X(-1);
  r.value()+=c; 
  return r;
}

template<class X, class R> inline
PolynomialVariable<X> 
operator-(const PolynomialVariable<X>& x, const R& c)
{
  PolynomialVariable<X> r=x; r.value()-=c; return r;
}

template<class X, class R> inline
PolynomialVariable<X> 
operator*(const R& c, const PolynomialVariable<X>& x)
{
  PolynomialVariable<X> r=x; 
  reinterpret_cast<Vector<X>&>(r.data())*=X(c);
  return r;
}

template<class X, class R> inline
PolynomialVariable<X> 
operator*(const PolynomialVariable<X>& x, const R& c)
{
  return c*x;
}

template<class X, class R> inline
PolynomialVariable<X> 
operator/(const R& c, const PolynomialVariable<X>& x)
{
  return c*rec(x);
}

template<class X, class R> inline
PolynomialVariable<X> 
operator/(const PolynomialVariable<X>& x, const R& c)
{
  return X(1/c)*x;
}


} // namespace Ariadne

#include "polynomial_variable.template.h"




#endif /* ARIADNE_POLYNOMIAL_VARIABLE_H */

