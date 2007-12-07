/***************************************************************************
 *            multivariable_derivative.h
 *
 *  Copyright 2007  Pieter Collins <Pieter.Collins@cwi.nl>
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
 
/*! \file multivariable_derivative.h
 *  \brief Derivatives of scalar functions of many variables.
 */
 
#ifndef ARIADNE_MULTIVARIABLE_DERIVATIVE_H
#define ARIADNE_MULTIVARIABLE_DERIVATIVE_H

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

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace Function {
  
    class MultiIndex;
    template<class X> class ScalarDerivative;
    template<class X> class TaylorDerivative;
  
    /*!\ingroup Function
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to a multiple arguments.
     */
    template<class X>
    class MultivariableDerivative
    {
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      MultivariableDerivative();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      MultivariableDerivative(uint r, uint a, uint d);
      /*! \brief The constant \a x of degree \a d in \a arguments. */
      MultivariableDerivative(uint r, uint a, uint d, const X& c);
      /*! \brief A multivariable derivative of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> MultivariableDerivative(uint r, uint a, uint d, const XX* ptr);

      /*! \brief Copy constructor. */
      template<class XX> MultivariableDerivative(const MultivariableDerivative<XX>& other); 
      /*! \brief Copy assignment operator. */
      template<class XX> MultivariableDerivative<X>& operator=(const MultivariableDerivative<XX>& other);

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const MultivariableDerivative<XX>& other);
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const MultivariableDerivative<XX>& other);

      /*! \brief The number of variables of the result. */
      uint result_size() const; 
      /*! \brief The number of variables of the argument. */
      uint argument_size() const; 
      /*! \brief The degree (number of derivatives computed). */
      uint degree() const; 
      /*! \brief The value of the vector derivative. */
      const LinearAlgebra::Vector<X> value() const;
      /*! \brief The Jacobian derivative of the vector quantity. */
      const LinearAlgebra::Matrix<X> jacobian() const;
      /*! \brief The array of derivative values. */
      const array<X>& data() const;
      /*! \brief A reference to the array of derivative values. */
      array<X>& data();

      /*! \brief Get the \a j<sup>th</sup> derivative of the \a i<sup>th</sup> component. */
      const X& get(const uint& i, const MultiIndex& j) const; 
      /*! \brief Set the \a j<sup>th</sup> derivative of the \a i<sup>th</sup> component. */
      template<class XX> void set(const uint& i, const MultiIndex& j, const XX& x); 
      /*! \brief Set the \a i<sup>th</sup> component to the scalar quantity \a x. */
      template<class XX> void set(const uint& i, const TaylorDerivative<XX>& x); 

      /*! \brief The \a i<sup>th</sup> component of the multivariable derivative. */
      const TaylorDerivative<X> operator[](const uint& i) const; 
      /*! \brief A reference to the \a i<sup>th</sup> component. */
      TaylorDerivative<X>& operator[](const uint& i); 
      /*! \brief The \a j<sup>th</sup> derivative of the \a i<sup>th</sup> component. */
      const X& operator()(const uint& i, const MultiIndex& j) const; 
      /*! \brief A reference to the \a j<sup>th</sup> derivative of the \a i<sup>th</sup> component. */
      X& operator()(const uint& i, const MultiIndex& j); 
#ifdef DOXYGEN
    //@{ 
    //! \name Friend operations
    /*! \brief The composition of two derivatives computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
     *  The composition inductively by
     *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
     */
    friend MultivariableDerivative<X> compose(const MultivariableDerivative<X>& y, const MultivariableDerivative<X>& x);
    /*! \brief The derivative of the inverse of \f$y\f$ evaluated at \f$x\f$. */
    friend MultivariableDerivative<X> inverse(const MultivariableDerivative<X>& x, const Vector<X>& c);
    /*! \brief The implicit function defined by \f$z(x,y)=\mathrm{const}\f$. */
    friend MultivariableDerivative<X> implicit(const MultivariableDerivative<X>& z, const MultivariableDerivative<X>& x);

    /*! \brief The derivatives of \f$-x\f$. */
    friend MultivariableDerivative<X> operator-(const MultivariableDerivative<X>& x);
    /*! \brief The derivatives of \f$x+y\f$. */
    friend MultivariableDerivative<X> operator+(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend MultivariableDerivative<X> operator-(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);
    /*! \brief The derivatives of \f$s*x\f$. */
    friend MultivariableDerivative<X> operator*(const TaylorDerivative<X>& s, const MultivariableDerivative<X>& x);
    /*! \brief The derivatives of \f$x*s\f$. */
    friend MultivariableDerivative<X> operator*(const MultivariableDerivative<X>& x, const TaylorDerivative<X>& s);
    /*! \brief The derivatives of \f$x/s\f$. */
    friend MultivariableDerivative<X> operator/(const MultivariableDerivative<X>& x, const TaylorDerivative<X>& s);

    /*! \brief Stream output operator. */
    friend std::ostream& operator<<(std::ostream& os, const MultivariableDerivative<X>& x);
    //@}
#endif 
     private:
      uint _result_size;
      uint _argument_size;
      uint _degree;
      array<X> _data;
    };


  template<class X0, class X1, class X2> void compute_composition(MultivariableDerivative<X0>& z, const MultivariableDerivative<X1>& y, const MultivariableDerivative<X2>& x);

  template<class X> MultivariableDerivative<X> compose(const ScalarDerivative<X>& y, const MultivariableDerivative<X>& x);
  template<class X> MultivariableDerivative<X> reduce(const MultivariableDerivative<X>& x, const uint& k);

  template<class X> MultivariableDerivative<X> neg(const MultivariableDerivative<X>& x);
  template<class X> MultivariableDerivative<X> add(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);
  template<class X> MultivariableDerivative<X> sub(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);

  template<class X> MultivariableDerivative<X> operator-(const MultivariableDerivative<X>& x);
  template<class X> MultivariableDerivative<X> operator+(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);
  template<class X> MultivariableDerivative<X> operator-(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);
  template<class X> MultivariableDerivative<X> operator*(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);
  template<class X> MultivariableDerivative<X> operator/(const MultivariableDerivative<X>& x, const MultivariableDerivative<X>& y);

  template<class X, class R> MultivariableDerivative<X> operator+(const MultivariableDerivative<X>& x, const R& c);
  template<class X, class R> MultivariableDerivative<X> operator+(const R& c, const MultivariableDerivative<X>& x);
  template<class X, class R> MultivariableDerivative<X> operator-(const MultivariableDerivative<X>& x, const R& c);
  template<class X, class R> MultivariableDerivative<X> operator-(const R& c, const MultivariableDerivative<X>& x);
  template<class X, class R> MultivariableDerivative<X> operator*(const MultivariableDerivative<X>& x, const R& c);
  template<class X, class R> MultivariableDerivative<X> operator*(const R& c, const MultivariableDerivative<X>& x);
  template<class X, class R> MultivariableDerivative<X> operator/(const MultivariableDerivative<X>& x, const R& c);
  template<class X, class R> MultivariableDerivative<X> operator/(const R& c, const MultivariableDerivative<X>& x);

  template<class X> std::ostream& operator<<(std::ostream& os, const MultivariableDerivative<X>& x);


  }
}


#include "multivariable_derivative.inline.h"
#include "multivariable_derivative.template.h"


#endif /* ARIADNE_MULTIVARIABLE_DERIVATIVE_H */

