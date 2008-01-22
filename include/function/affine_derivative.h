/***************************************************************************
 *            affine_derivative.h
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
 
/*! \file affine_derivative.h
 *  \brief First derivatives with respect to multiple arguments.
 */
 
#ifndef ARIADNE_AFFINE_DERIVATIVE_H
#define ARIADNE_AFFINE_DERIVATIVE_H

#include <iostream>
#include <stdexcept>
#include <cassert>

#include "base/tribool.h"
#include "base/exceptions.h"
#include "base/stlio.h"

#include "numeric/exceptions.h"
#include "numeric/traits.h"
#include "numeric/conversion.h"
#include "numeric/arithmetic.h"
#include "numeric/function.h"

#include "linear_algebra/covector.h"
#include "function/affine_variable.h"

namespace Ariadne {
  namespace Function {

    template<class X> class AffineVariableReference;

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
    template<class X>
    class AffineDerivative
    {
     private:
      uint _result_size;
      uint _argument_size;
      array< AffineVariable<X> > _variables;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs the differential of the constant zero with respect to no independent variable. */
      AffineDerivative();
      /*! \brief Construct the zero constant in \a r variables with respect to \a a arguments. */
      AffineDerivative(const uint& r, const uint& a);
      /*! \brief Construct a first derivative in \a r variables with respect to \a a arguments from the data starting at \a ptr. */
      template<class XX> AffineDerivative(const uint& r, const uint& a, const XX* ptr);

      /*! \brief Type conversion copy constructor. */
      template<class XX> AffineDerivative(const AffineDerivative<XX>& fd);
      /*! \brief Type conversion copy assignment. */
      template<class XX> AffineDerivative<X>& operator=(const AffineDerivative<XX>& fd);

      /*! \brief Construct a constant with the value \a x */
      static AffineDerivative<X> constant(size_type as, const LinearAlgebra::Vector<X>& x);
      
      /*! \brief Construct a variable with the value \a x */
      static AffineDerivative<X> variable(const LinearAlgebra::Vector<X>& x);
      
      //@}

      //@{
      //! \name Data access
      /*! \brief The number of variables of the result. */
      size_type result_size() const;
      /*! \brief The number of variables of the argument. */
      size_type argument_size() const;
      /*! \brief The degree of the derivative map. Returns the constant one. */
      smoothness_type degree() const;
      /*! \brief The variables. */
      const array<AffineVariable<X> >& variables() const;
      /*! \brief Resize to hold the derivatives of a function with result size \a rs and argument size \a as. */
      void resize(const size_type& rs, const size_type& as);

      /*! \brief A reference to the \a i<sup>th</sup> component. */
      AffineVariable<X>& operator[](size_type i);
      /*! \brief A constant reference to the \a i<sup>th</sup> component. */
      const AffineVariable<X>& operator[](size_type i) const;

      /*! \brief The value of the variable. */
      LinearAlgebra::Vector<X> value() const;
      /*! \brief The jacobian derivative. */
      LinearAlgebra::Matrix<X> jacobian() const;
      //@}
      
      
    };

    template<class X1, class X2> bool operator==(const AffineDerivative<X1>& x1, const AffineDerivative<X2>& x2); 
    template<class X1, class X2> bool operator!=(const AffineDerivative<X1>& x1, const AffineDerivative<X2>& x2); 

    template<class X> std::ostream& operator<<(std::ostream& os, const AffineDerivative<X>& x);

    template<class X> AffineDerivative<X> operator+(const AffineDerivative<X>& x);
    template<class X> AffineDerivative<X> operator-(const AffineDerivative<X>& x);
    template<class X> AffineDerivative<X> operator+(const AffineDerivative<X>& x1, const AffineDerivative<X>& x2); 
    template<class X> AffineDerivative<X> operator-(const AffineDerivative<X>& x1, const AffineDerivative<X>& x2); 



  }

}

#include "affine_derivative.inline.h"

#endif /* ARIADNE_AFFINE_DERIVATIVE_H */
