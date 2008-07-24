/***************************************************************************
 *            differential.h
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
 
/*! \file differential.h
 *  \brief Derivatives of scalar functions of many variables.
 */
 
#ifndef ARIADNE_DIFFERENTIAL_H
#define ARIADNE_DIFFERENTIAL_H


namespace Ariadne {
  
    class MultiIndex;
    template<class X> class PowerSeries;
    template<class X> class Differential;
    template<class X> class DifferentialVector;
  
    /*!\ingroup Differentiation
     * \brief A templated class representing a the derivatives of a scalar quantity with respect to multiple arguments. 
     */
    template<class X>
    class Differential
    {
      friend class DifferentialVector<X>;
     public:
      /*! The type used to represent indices. */
      typedef MultiIndex index_type;
      /*! The type used to represent numbers. */
      typedef X value_type;
      /*! The type used for a vector of variables. */
      typedef DifferentialVector<X> vector_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      Differential();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      Differential(size_type a, smoothness_type d);
      /*! \brief A taylor variable of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> Differential(size_type a, smoothness_type d, const XX* ptr);

      /*! \brief Construct from a univariate Taylor series. */
      template<class XX> Differential(const PowerSeries<XX>& ts); 
      /*! \brief Copy constructor. */
      template<class XX> Differential(const Differential<XX>& tv); 
      /*! \brief Copy assignment operator. */
      template<class XX> Differential<X>& operator=(const Differential<XX>& tv);

      /*! \brief Assign a constant \a c. */
      template<class XX> Differential<X>& operator=(const XX& c);

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const Differential<XX>& other) const;
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const Differential<XX>& other) const;

      /*! \brief Construct a constant variable of degree \a d with respect to \a as variables and value \a c. */
      template<class XX> static Differential<X> constant(size_type as, smoothness_type d, const XX& c); 
      /*! \brief Construct the variable of degree \a d at value \a value with respect to the \a i<sup>th</sup> variable of \a as. */
      template<class XX> static Differential<X> variable(size_type as, smoothness_type d, const XX& value, size_type i);

      /*! \brief The number of variables of the argument. */
      size_type argument_size() const; 
      /*! \brief The degree (number of derivatives computed). */
      smoothness_type degree() const; 
      /*! \brief The value of the quantity. */
      const X& value() const;
      /*! \brief A reference to the value of the quantity. */
      X& value();
      /*! \brief The variation of the quantity with respect to the \a j<sup>th</sup> argument. */
      const X& gradient(size_type j) const;
      /*! \brief The array of derivative values. */
      const array<X>& data() const;
      /*! \brief A reference to the array of derivative values. */
      array<X>& data();
      /*! \brief A reference to the \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$. */
      X& operator[](const MultiIndex& a); 
      /*! \brief The \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$. */
      const X& operator[](const MultiIndex& a) const; 

      /*! \brief Assign all elements of degree less than the degree of \a x to those of \a x . */
      Differential<X>& assign(const Differential<X>& x);
#ifdef DOXYGEN
      //@{ 
      //! \name Friend operations
      /*! \brief The composition of two variables computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
       *  The composition inductively by
       *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
       *  \f$ y = a_0 + x ( a_1 + x ( a_2/2 + x ( a_3/3! + \cdots)))\f$.
       */
      friend Differential<X> compose(const Differential<X>& y, const Differential<X>& x);
      /*! \brief Embed \a x in a space of independent variables of dimension \a size, with the first independent variable of \a x becoming the \a start independen variable of the result. */
      friend Differential<X> embed(const Differential<X>& x, uint size, uint start);
      /*! \brief The derivatives of the inverse of \f$y\f$ evaluated at \f$x\f$. (Not currently implemented.) */
      friend Differential<X> inverse(const Differential<X>& y, const X& x);
      /*! \brief The derivative of \f$x\f$ with respect to the variable \a k .*/
      friend Differential<X> derivative(const Differential<x>& x, size_type k);
      /*! \brief The antiderivative of \f$x\f$ with respect to the variable \a k .*/
      friend Differential<X> antiderivative(const Differential<x>& x, size_type k);
      /*! \brief The minimum of two variables. Returns the variable whose zero-th order value is minimal. */
      friend Differential<X> min(const Differential<X>& x1, const Differential<X>& x2);
      /*! \brief The maximum of two variables. Returns the variable whose zero-th order value is maximal. */
      friend Differential<X> max(const Differential<X>& x1, const Differential<X>& x2);
      /*! \brief The derivatives of \f$+x\f$. Returns a copy. */
      friend Differential<X> pos(const Differential<X>& x);
      /*! \brief The derivatives of \f$-x\f$. */
      friend Differential<X> neg(const Differential<X>& x);
      /*! \brief The derivatives of \f$x+y\f$. */
      friend Differential<X> add(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x-y\f$. */
      friend Differential<X> sub(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x*y\f$. */
      friend Differential<X> mul(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x/y\f$. */
      friend Differential<X> div(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x^n\f$. */
      friend Differential<X> pow(const Differential<X>& x, const Integer& n);
      /*! \brief The derivatives of \f$\sqrt{x}\f$. */
      friend Differential<X> sqrt(const Differential<X>& x);
      /*! \brief The derivatives of \f$\exp(x)\f$. */
      friend Differential<X> exp(const Differential<X>& x);
      /*! \brief The derivatives of \f$\log(x)\f$. */
      friend Differential<X> log(const Differential<X>& x);
      /*! \brief The derivatives of \f$\sin(x)\f$. */
      friend Differential<X> sin(const Differential<X>& x);
      /*! \brief The derivatives of \f$\cos(x)\f$. */
      friend Differential<X> cos(const Differential<X>& x);
      /*! \brief The derivatives of \f$\tan(x)\f$. */
      friend Differential<X> tan(const Differential<X>& x);
      /*! \brief The derivatives of \f$\sin^{-1}(x)\f$. (Not currently implemented.) */
      friend Differential<X> asin(const Differential<X>& x);
      /*! \brief The derivatives of \f$\cos^{-1}(x)\f$. (Not currently implemented.) */
      friend Differential<X> acos(const Differential<X>& x);
      /*! \brief The derivatives of \f$\tan^{-1}(x)\f$. (Not currently implemented.) */
      friend Differential<X> atan(const Differential<X>& x);

      /*! \brief The derivatives of \f$x+y\f$. */
      friend Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x-y\f$. */
      friend Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x*y\f$. */
      friend Differential<X> operator*(const Differential<X>& x, const Differential<X>& y);
      /*! \brief The derivatives of \f$x/y\f$. */
      friend Differential<X> operator/(const Differential<X>& x, const Differential<X>& y);

      /*! \brief The derivatives of \f$c+x\f$ for a constant \f$c\f$. (Other mixed-mode arithmetic is also supported.) */
      friend Differential<X> operator+(const R& c, const Differential<X>& x);

      /*! \brief Stream output operator. */
      friend std::ostream& operator<<(std::ostream& os, const Differential<X>& x);
      //@}
#endif 
     private:
      static void instantiate();
     private:
      size_type _argument_size;
      smoothness_type _degree;
      array<X> _data;
    };

  template<class X> bool operator<(const Differential<X>& x1, const Differential<X>& x2);

  template<class X, class R> R evaluate(const Differential<X>& y, const array<R>& x);

  template<class X> Differential<X> compose(const PowerSeries<X>& y, const Differential<X>& x);
  template<class X> Differential<X> compose(const Differential<X>& y, const Differential<X>& x);
  template<class X> Differential<X> reduce(const Differential<X>& x);
  template<class X> Differential<X> derivative(const Differential<X>& x, size_type k);
  template<class X> Differential<X> antiderivative(const Differential<X>& x, size_type k);

  template<class X> Differential<X> min(const Differential<X>& x1, const Differential<X>& x2); 
  template<class X> Differential<X> max(const Differential<X>& x1,const Differential<X>& x2); 
  template<class X> Differential<X> pos(const Differential<X>& x);
  template<class X> Differential<X> neg(const Differential<X>& x);
  template<class X> Differential<X> abs(const Differential<X>& x);
  template<class X> Differential<X> rec(const Differential<X>& x);
  template<class X> Differential<X> add(const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X> sub(const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X> mul(const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X> div(const Differential<X>& x, const Differential<X>& y);
  template<class X, class N> Differential<X> pow(const Differential<X>& x, N k);

  template<class X> Differential<X>& acc(Differential<X>& r, const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X>& acc(Differential<X>& r, const X& c, const Differential<X>& x);

  template<class X> Differential<X> sqrt(const Differential<X>& x);
  template<class X> Differential<X> exp(const Differential<X>& x); 
  template<class X> Differential<X> log(const Differential<X>& x); 
  template<class X> Differential<X> sin(const Differential<X>& x); 
  template<class X> Differential<X> cos(const Differential<X>& x); 
  template<class X> Differential<X> tan(const Differential<X>& x); 
  template<class X> Differential<X> asin(const Differential<X>& x); 
  template<class X> Differential<X> acos(const Differential<X>& x); 
  template<class X> Differential<X> atan(const Differential<X>& x); 

  template<class X> Differential<X> operator+(const Differential<X>& x);
  template<class X> Differential<X> operator-(const Differential<X>& x);
  template<class X> Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X> operator*(const Differential<X>& x, const Differential<X>& y);
  template<class X> Differential<X> operator/(const Differential<X>& x, const Differential<X>& y);

  template<class X, class R> Differential<X> operator+(const Differential<X>& x, const R& c);
  template<class X, class R> Differential<X> operator+(const R& c, const Differential<X>& x);
  template<class X, class R> Differential<X> operator-(const Differential<X>& x, const R& c);
  template<class X, class R> Differential<X> operator-(const R& c, const Differential<X>& x);
  template<class X, class R> Differential<X> operator*(const Differential<X>& x, const R& c);
  template<class X, class R> Differential<X> operator*(const R& c, const Differential<X>& x);
  template<class X, class R> Differential<X> operator/(const Differential<X>& x, const R& c);
  template<class X, class R> Differential<X> operator/(const R& c, const Differential<X>& x);

  /*
  template<class X> Differential<X> operator+(const Differential<X>& x, const X& c);
  template<class X> Differential<X> operator+(const X& c, const Differential<X>& x);
  template<class X> Differential<X> operator-(const Differential<X>& x, const X& c);
  template<class X> Differential<X> operator-(const X& c, const Differential<X>& x);
  template<class X> Differential<X> operator*(const Differential<X>& x, const X& c);
  template<class X> Differential<X> operator*(const X& c, const Differential<X>& x);
  template<class X> Differential<X> operator/(const Differential<X>& x, const X& c);
  template<class X> Differential<X> operator/(const X& c, const Differential<X>& x);

  template<class X> Differential<X> operator+(const Differential<X>& x, const double& c);
  template<class X> Differential<X> operator+(const double& c, const Differential<X>& x);
  template<class X> Differential<X> operator-(const Differential<X>& x, const double& c);
  template<class X> Differential<X> operator-(const double& c, const Differential<X>& x);
  template<class X> Differential<X> operator*(const Differential<X>& x, const double& c);
  template<class X> Differential<X> operator*(const double& c, const Differential<X>& x);
  template<class X> Differential<X> operator/(const Differential<X>& x, const double& c);
  template<class X> Differential<X> operator/(const double& c, const Differential<X>& x);
  */

  template<class X> Differential<X>& operator+=(Differential<X>& r, const Differential<X>& x);
  template<class X> Differential<X>& operator-=(Differential<X>& r, const Differential<X>& x);

  template<class X, class R> Differential<X>& operator+=(Differential<X>& x, const R& c);
  template<class X, class R> Differential<X>& operator-=(Differential<X>& x, const R& c);
  template<class X, class R> Differential<X>& operator*=(Differential<X>& x, const R& c);
  template<class X, class R> Differential<X>& operator/=(Differential<X>& x, const R& c);

  template<class X> std::ostream& operator<<(std::ostream& os, const Differential<X>& x);


  
} // namespace Ariadne


#include "differential.inline.h"
#include "differential.template.h"


#endif /* ARIADNE_DIFFERENTIAL_H */

