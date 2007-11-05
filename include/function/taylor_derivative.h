/***************************************************************************
 *            taylor_derivative.h
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
 
/*! \file taylor_derivative.h
 *  \brief Derivatives of scalar functions of many variables.
 */
 
#ifndef ARIADNE_TAYLOR_DERIVATIVE_H
#define ARIADNE_TAYLOR_DERIVATIVE_H

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
  
    class MultiIndex;
    template<class X> class ScalarDerivative;
  
    /*!\ingroup Function
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a multiple arguments.
     */
    template<class X>
    class TaylorDerivative
    {
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      TaylorDerivative();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      TaylorDerivative(uint a, uint d);
      /*! \brief The constant \a x of degree \a d in \a arguments. */
      TaylorDerivative(uint a, uint d, const X& c);
      /*! \brief The  \a i<sup>th</sup> variable of \a a, with value \a x and degree \a d. */
      TaylorDerivative(uint a, uint d, uint i, const X& x);
      /*! \brief A taylor derivative of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> TaylorDerivative(uint a, uint d, const XX* ptr);

      /*! \brief Copy constructor. */
      template<class XX> TaylorDerivative(const TaylorDerivative<XX>& other); 
      /*! \brief Copy assignment operator. */
      template<class XX> TaylorDerivative<X>& operator=(const TaylorDerivative<XX>& other);

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const TaylorDerivative<XX>& other);
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const TaylorDerivative<XX>& other);

      /*! \brief Construct a constant derivative of degree \a d with respect to \a as variables and value \a c. */
      static TaylorDerivative<X> constant(uint as, uint d, const X& c); 
      /*! \brief Construct the derivative of degree \a d at value \a value with respect to the \a i<sup>th</sup> variable of \a as. */
      static TaylorDerivative<X> variable(uint as, uint d, uint i, const X& value);

      /*! \brief The number of variables of the argument. */
      uint argument_size() const; 
      /*! \brief The degree (number of derivatives computed). */
      uint degree() const; 
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
#ifdef DOXYGEN
    //@{ 
    //! \name Friend operations
    /*! \brief The composition of two derivatives computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
     *  The composition inductively by
     *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
     */
    friend TaylorDerivative<X> compose(const TaylorDerivative<X>& y, const TaylorDerivative<X>& x);
    /*! \brief The derivative of the inverse of \f$y\f$ evaluated at \f$x\f$. (Not currently implemented.) */
    friend TaylorDerivative<X> inverse(const TaylorDerivative<X>& y, const X& x);
    /*! \brief The minimum of two derivatives. Returns the derivative whose zero-th order value is minimal. */
    friend TaylorDerivative<X> min(const TaylorDerivative<X>& x1, const TaylorDerivative<X>& x2);
    /*! \brief The maximum of two derivatives. Returns the derivative whose zero-th order value is maximal. */
    friend TaylorDerivative<X> max(const TaylorDerivative<X>& x1, const TaylorDerivative<X>& x2);
    /*! \brief The derivatives of \f$+x\f$. Returns a copy. */
    friend TaylorDerivative<X> pos(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$-x\f$. */
    friend TaylorDerivative<X> neg(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$x+y\f$. */
    friend TaylorDerivative<X> add(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend TaylorDerivative<X> sub(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x*y\f$. */
    friend TaylorDerivative<X> mul(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x/y\f$. */
    friend TaylorDerivative<X> div(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x^n\f$. */
    friend TaylorDerivative<X> pow(const TaylorDerivative<X>& x, const Integer& n);
    /*! \brief The derivatives of \f$\sqrt{x}\f$. */
    friend TaylorDerivative<X> sqrt(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\exp(x)\f$. */
    friend TaylorDerivative<X> exp(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\log(x)\f$. */
    friend TaylorDerivative<X> log(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\sin(x)\f$. */
    friend TaylorDerivative<X> sin(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\cos(x)\f$. */
    friend TaylorDerivative<X> cos(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\tan(x)\f$. */
    friend TaylorDerivative<X> tan(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\sin^{-1}(x)\f$. (Not currently implemented.) */
    friend TaylorDerivative<X> asin(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\cos^{-1}(x)\f$. (Not currently implemented.) */
    friend TaylorDerivative<X> acos(const TaylorDerivative<X>& x);
    /*! \brief The derivatives of \f$\tan^{-1}(x)\f$. (Not currently implemented.) */
    friend TaylorDerivative<X> atan(const TaylorDerivative<X>& x);

    /*! \brief The derivatives of \f$x+y\f$. */
    friend TaylorDerivative<X> operator+(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend TaylorDerivative<X> operator-(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x*y\f$. */
    friend TaylorDerivative<X> operator*(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
    /*! \brief The derivatives of \f$x/y\f$. */
    friend TaylorDerivative<X> operator/(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);

    /*! \brief The derivatives of \f$c+x\f$ for a constant \f$c\f$. (Other mixed-mode arithmetic is also supported.) */
    friend TaylorDerivative<X> operator+(const R& c, const TaylorDerivative<X>& x);

    /*! \brief Stream output operator. */
    friend std::ostream& operator<<(std::ostream& os, const TaylorDerivative<X>& x);
    //@}
#endif 
     private:
      uint _argument_size;
      uint _degree;
      array<X> _data;
    };


  template<class X0, class X1, class X2> void compute_product(TaylorDerivative<X0>& x0, const TaylorDerivative<X1>& x1, const TaylorDerivative<X2>& x2);
  template<class X0, class X1, class X2> void compute_composition(TaylorDerivative<X0>& z, const ScalarDerivative<X1>& y, const TaylorDerivative<X2>& x);

  template<class X> TaylorDerivative<X> compose(const ScalarDerivative<X>& y, const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> reduce(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> derivative(const TaylorDerivative<X>& x, const size_type& k);
  template<class X> TaylorDerivative<X> min(const TaylorDerivative<X>& x1, const TaylorDerivative<X>& x2); 
  template<class X> TaylorDerivative<X> max(const TaylorDerivative<X>& x1,const TaylorDerivative<X>& x2); 
  template<class X> TaylorDerivative<X> pos(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> neg(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> abs(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> inv(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> add(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X> TaylorDerivative<X> sub(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X> TaylorDerivative<X> mul(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X> TaylorDerivative<X> div(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X, class N> TaylorDerivative<X> pow(const TaylorDerivative<X>& x, N k);

  template<class X> TaylorDerivative<X> sqrt(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> exp(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> log(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> sin(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> cos(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> tan(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> asin(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> acos(const TaylorDerivative<X>& x); 
  template<class X> TaylorDerivative<X> atan(const TaylorDerivative<X>& x); 

  template<class X> TaylorDerivative<X> operator-(const TaylorDerivative<X>& x);
  template<class X> TaylorDerivative<X> operator+(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X> TaylorDerivative<X> operator-(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X> TaylorDerivative<X> operator*(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);
  template<class X> TaylorDerivative<X> operator/(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y);

  template<class X, class R> TaylorDerivative<X> operator+(const TaylorDerivative<X>& x, const R& c);
  template<class X, class R> TaylorDerivative<X> operator+(const R& c, const TaylorDerivative<X>& x);
  template<class X, class R> TaylorDerivative<X> operator-(const TaylorDerivative<X>& x, const R& c);
  template<class X, class R> TaylorDerivative<X> operator-(const R& c, const TaylorDerivative<X>& x);
  template<class X, class R> TaylorDerivative<X> operator*(const TaylorDerivative<X>& x, const R& c);
  template<class X, class R> TaylorDerivative<X> operator*(const R& c, const TaylorDerivative<X>& x);
  template<class X, class R> TaylorDerivative<X> operator/(const TaylorDerivative<X>& x, const R& c);
  template<class X, class R> TaylorDerivative<X> operator/(const R& c, const TaylorDerivative<X>& x);

  template<class X> std::ostream& operator<<(std::ostream& os, const TaylorDerivative<X>& x);


  }
}


#include "taylor_derivative.inline.h"
#include "taylor_derivative.template.h"


#endif /* ARIADNE_TAYLOR_DERIVATIVE_H */

