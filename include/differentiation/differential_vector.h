/***************************************************************************
 *            differential_vector.h
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
 
/*! \file differential_vector.h
 *  \brief Derivatives of scalar functions of many variables.
 */
 
#ifndef ARIADNE_DIFFERENTIAL_VECTOR_H
#define ARIADNE_DIFFERENTIAL_VECTOR_H

#include "linear_algebra/declarations.h"

namespace Ariadne {
  
    class MultiIndex;
    template<class X> class PowerSeries;
    template<class X> class Differential;
  
    // Base template algorithm for evaluating a polynomial
    template<class P, class M> void evaluate_polynomial(M& r, const P& p, const M& x);

    /*!\ingroup Differentiation
     * \brief A templated class representing a the derivatives of a vector quantity with respect to multiple arguments.
     *
     * Formally, a Taylor variable \f$y\f$ is a member of a differential algebra over some field, denoted by the template variable X. 
     * A Taylor variable therefore supports addition, subtraction and multiplication, and also division, as long as the division has a
     * non-zero \em value. 
     * We can think of \f$y\f$ as denoting a quiantity and it's first \f$d\f$ derivatives 
     * with respect to independent variables \f$x_1,\ldots,x_n\f$, 
     * 
     * In the current implementation, we store the coefficients of the power series expansion \f$y=\sum a_\alpha x^\alpha\f$ rather than
     * the derivative values \f$d^\alpha y/dx^\alpha\f$, as the former are more convenient to work with.
     *
     */
    template<class X>
    class DifferentialVector
    {
     public:
      typedef X value_type;

      /*! \brief Default constructor constructs a constant of degree zero. */
      DifferentialVector();
      /*! \brief The constant zero of degree \a d in \a a arguments. */
      DifferentialVector(size_type r, size_type a, smoothness_type d);
      /*! \brief A taylor derivative of degree \a d in \a arguments, with values given by the array based at \a ptr. */
      template<class XX> DifferentialVector(size_type r, size_type a, smoothness_type d, const XX* ptr);

      /*! \brief Construct from a univariate Taylor series. */
      template<class XX> DifferentialVector(const PowerSeries<XX>& ts); 
      /*! \brief Construct from a scalar Taylor variable. */
      template<class XX> DifferentialVector(const Differential<XX>& tv); 
      /*! \brief Copy constructor. */
      template<class XX> DifferentialVector(const DifferentialVector<XX>& td); 
      /*! \brief Copy assignment operator. */
      template<class XX> DifferentialVector<X>& operator=(const DifferentialVector<XX>& td);

      /*! \brief Equality operator. */
      template<class XX> bool operator==(const DifferentialVector<XX>& other) const;
      /*! \brief Inequality operator. */
      template<class XX> bool operator!=(const DifferentialVector<XX>& other) const;

      /*! \brief Construct a constant derivative of degree \a d with respect to \a as variables and value \a c. */
      template<class V> static DifferentialVector<X> constant(size_type rs, size_type as, smoothness_type d, const V& c); 
      /*! \brief Construct the derivative of degree \a d at values \a x. Requires rs==as. */
      template<class V> static DifferentialVector<X> variable(size_type rs, size_type as, smoothness_type d, const V& x);
      template<class V> static DifferentialVector<X> variable(const V& x, smoothness_type d);

      /*! \brief The number of variables of the argument. */
      size_type result_size() const; 
      // Synonym for result size; needed by template algorithm
      size_type size() const; 
      /*! \brief The number of variables of the argument. */
      size_type argument_size() const; 
      /*! \brief The degree (number of derivatives computed). */
      smoothness_type degree() const; 
      /*! \brief The value of the vector quantity. */
      Vector<X> value() const; 
      /*! \brief The value of the vector quantity. */
      Matrix<X> jacobian() const; 

      /*! \brief Set the value of the quantity. */
      void set_value(const Vector<X>&); 
      /*! \brief Set the Jacobian matrix of the quantity. */
      void set_jacobian(const Matrix<X>&); 

      /*! \brief The derivative of the \a i<sup>th</sup> variable with respect to multi-index \a j. */
      const X& get(const size_type& i, const MultiIndex& j) const;
      /*! \brief Set the derivative of the \a i<sup>th</sup> variable with respect to multi-index \a j. */
      template<class XX> void set(const size_type& i, const MultiIndex& j, const XX& x);
      /*! \brief The derivative values of the \a i<sup>th</sup> variable. */
      const Differential<X>& get(const size_type& i) const;
      /*! \brief Set the derivative of the \a i<sup>th</sup> variable. */
      template<class XX> void set(const size_type& i, const Differential<XX>& tv);
      /*! \brief The array of Taylor variables. */
      const array< Differential<X> >& variables() const;
      /*! \brief The array of derivative values. */
      //const array<X>& data() const;
      /*! \brief A reference to the array of derivative values. */
      //array<X>& data();
      /*! \brief A reference to the \a i<sup> th</sup> component. */
      Differential<X>& operator[](const size_type& i); 
      /*! \brief The \a i<sup> th</sup> component. */
      const Differential<X>& operator[](const size_type& i) const; 

#ifdef DOXYGEN
    //@{ 
    //! \name Friend operations
    /*! \brief The composition of two derivatives computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
     *  The composition inductively by
     *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
     */
    friend Differential<X> compose(const Differential<X>& y, const DifferentialVector<X>& x);
    /*! \brief The composition of two derivatives computes \f$d^iy/dt^i\f$ from \f$d^iy/dx^i\f$ and \f$d^ix/dt^i\f$. 
     *  The composition inductively by
     *  \f$ y^{[n]} = \sum_{i=0}^{n-1} \Bigl(\!\begin{array}{c}n\\i\end{array}\!\Bigr) {\dot{y}}^{[i]} x^{(n-i)} \f$
     */
    friend DifferentialVector<X> compose(const DifferentialVector<X>& y, const DifferentialVector<X>& x);
    /*! \brief The derivative of the inverse of \f$y\f$, assuming \f$c\f$ is the centre of the approximation. */
    friend DifferentialVector<X> inverse(const DifferentialVector<X>& y, const Vector<X>& c);
    /*! \brief The taylor derivatives of the function \f$z:\R^m\rightarrow\R^n\f$ defined by \f$y(x,z(x))=\text{const}\f$. The value of \f$z(x)\f$ is set to be \f$c\f$. */
    friend DifferentialVector<X> implicit(const DifferentialVector<X>& y, const Vector<X>& c);
    /*! \brief The derivatives of \f$x+y\f$. */
    friend DifferentialVector<X> add(const DifferentialVector<X>& x, const DifferentialVector<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend DifferentialVector<X> sub(const DifferentialVector<X>& x, const DifferentialVector<X>& y);

    /*! \brief The derivatives of \f$x+y\f$. */
    friend DifferentialVector<X> operator+(const DifferentialVector<X>& x, const DifferentialVector<X>& y);
    /*! \brief The derivatives of \f$x-y\f$. */
    friend DifferentialVector<X> operator-(const DifferentialVector<X>& x, const DifferentialVector<X>& y);

    /*! \brief Stream output operator. */
    friend std::ostream& operator<<(std::ostream& os, const DifferentialVector<X>& x);
    //@}
#endif 
     private:
      size_type _increment() const;
      static void instantiate();
     private:
      size_type _result_size;
      size_type _argument_size;
      smoothness_type _degree;
      array< Differential<X> > _variables;
    };

  template<class X, class R> array<R> evaluate(const DifferentialVector<X>& y, const array<R>& x);

  template<class X0, class X1, class X2> void compute_composition(Differential<X0>& z, const Differential<X1>& y, const DifferentialVector<X2>& x);
  template<class X0, class X1, class X2> void compute_composition(DifferentialVector<X0>& z, const DifferentialVector<X1>& y, const DifferentialVector<X2>& x);

  template<class X> Differential<X> evaluate(const Differential<X>& y, const DifferentialVector<X>& x);
  template<class X> Differential<X> compose(const Differential<X>& y, const DifferentialVector<X>& x);

  template<class X> Vector<X> evaluate(const DifferentialVector<X>& y, const Vector<X>& x);

  template<class X> DifferentialVector<X> translate(const DifferentialVector<X>& x, const Vector<X>& c);
  template<class X> DifferentialVector<X> compose(const DifferentialVector<X>& y, const DifferentialVector<X>& x);
  template<class X> DifferentialVector<X> inverse(const DifferentialVector<X>& x, const Vector<X>& c);
  template<class X> DifferentialVector<X> implicit(const DifferentialVector<X>& x, const Vector<X>& c);
  template<class X> DifferentialVector<X> concatenate(const DifferentialVector<X>& x, const DifferentialVector<X>& y);
  template<class X> DifferentialVector<X> reduce(const DifferentialVector<X>& x);
  template<class X> DifferentialVector<X> derivative(const DifferentialVector<X>& x, const size_type& k);
  template<class X> DifferentialVector<X> antiderivative(const DifferentialVector<X>& x, const size_type& k);
  
  template<class X>
  array< PowerSeries< Differential<X> > > 
  integrate(const array<Differential<X> >& y, const array<Differential<X> >& x);


  template<class X> DifferentialVector<X> neg(const DifferentialVector<X>& x);
  template<class X> DifferentialVector<X> add(const DifferentialVector<X>& x, const DifferentialVector<X>& y);
  template<class X> DifferentialVector<X> sub(const DifferentialVector<X>& x, const DifferentialVector<X>& y);

  template<class X> DifferentialVector<X> operator+(const DifferentialVector<X>& x);
  template<class X> DifferentialVector<X> operator-(const DifferentialVector<X>& x);
  template<class X> DifferentialVector<X> operator+(const DifferentialVector<X>& x, const DifferentialVector<X>& y);
  template<class X> DifferentialVector<X> operator-(const DifferentialVector<X>& x, const DifferentialVector<X>& y);

  template<class X> DifferentialVector<X> operator*(const Matrix<X>& A, const DifferentialVector<X>& x);

  template<class X> DifferentialVector<X>& operator-=(DifferentialVector<X>& x, const Vector<X>& v);
  template<class X> DifferentialVector<X> operator-(const DifferentialVector<X>& x, const Vector<X>& v);

  template<class X> DifferentialVector<X>& operator*=(DifferentialVector<X>& x, const double& c);
  template<class X> DifferentialVector<X>& operator*=(DifferentialVector<X>& x, const X& c);
  template<class X> DifferentialVector<X> operator*(const X& c, const DifferentialVector<X>& x);
  template<class X> DifferentialVector<X> operator*(const DifferentialVector<X>& x, const X& c);


  template<class X> std::ostream& operator<<(std::ostream& os, const DifferentialVector<X>& x);


  
} // namespace Ariadne


#include "differential_vector.inline.h"
#include "differential_vector.template.h"


#endif /* ARIADNE_DIFFERENTIAL_VECTOR_H */

