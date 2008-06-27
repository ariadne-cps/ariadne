/***************************************************************************
 *            function_interface.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
/*! \file function_interface.h
 *  \brief General function interface.
 */
 
#ifndef ARIADNE_FUNCTION_INTERFACE_H
#define ARIADNE_FUNCTION_INTERFACE_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "base/types.h"
#include "base/array.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"

namespace Ariadne {

    template<class X> class Vector;
    template<class X> class Matrix;
    class MultiIndex;
    template<class X> class SparseDifferentialVector;

    /*!\ingroup Function
     * \brief Abstract base class for continuous functionss.
     * 
     * The function is specified by the method operator()(const Vector< Interval<R> >& A) const,
     * This method should compute an interval vector \f$\overline{f}(A)\f$ with the
     * following properties:
     *   -# \f$f(A)\subset\overline{f}(A)\f$,
     *   -# If \f$A_1\subset A_0\f$, then \f$\overline{f}(A_1)\subset 
     *       \overline{f}(A_0)\f$, and
     *   -# If \f$\bigcap_{n\in\mathbb{N}}A_n=\{x\}\f$, then 
     *       \f$\bigcap_{n\in\mathbb{N}}\overline{f}(A_n)=\{f(x)\}\f$.
     *
     * More succinctly, we say that \f$\overline{f}(A_n)\f$ converges monotonically 
     * as \f$A_n\f$ tends to a point.
     *
     * Additional accuracy can be obtained be using derivatives.
     * The method derivative(const Vector< Interval<R> >& A) const computes the \a i th component of the derivative over the set \a A 
     * with respect to the variables in the multi-index \a j.
     */
    template<class R>
    class FunctionInterface {
      typedef typename traits<R>::approximate_arithmetic_type A;
      typedef typename traits<R>::arithmetic_type F;
      typedef typename traits<R>::interval_type I;
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      
      /*! \brief Virtual destructor. */
      virtual ~FunctionInterface();
     
      /*! \brief Make a copy (clone) of the function. */
      virtual FunctionInterface<R>* clone() const = 0;
     
      /*! \brief Evaluate the function. */
      Vector<F> operator() (const Vector<F>& x) const;

      /*! \brief Evaluate the function. */
      virtual Vector<F> evaluate(const Vector<F>& x) const = 0;
      /*! \brief Evaluate the Jacobian derivative matrix at the point \a x. */
      virtual Matrix<F> jacobian(const Vector<F>& x) const = 0;
      /*! \brief Evaluate the derivative of the function. */
      virtual TaylorDerivative<F> derivative(const Vector<F>& x, const smoothness_type& s) const = 0;
      /*! \brief Evaluate the derivative of the function using approximate arithmetic. */
      virtual SparseDifferentialVector<A> expansion(const Vector<A>& x, const smoothness_type& s) const = 0;

      /*! \brief The degree of differentiability of the function. */
      virtual smoothness_type smoothness() const = 0;
      /*! \brief The size of the function argument. */
      virtual size_type argument_size() const = 0;
      /*! \brief The size of the function result. */
      virtual size_type result_size() const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
    };
  



    template<class R> inline 
    FunctionInterface<R>::~FunctionInterface() {
    }; 

    template<class R> inline 
    Vector<typename FunctionInterface<R>::F>
    FunctionInterface<R>::operator() (const Vector<F>& x) const {
      return this->evaluate(x);
    }; 

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const FunctionInterface<R>& f) {
      return f.write(os);
    }; 



  
} // namespace Ariadne

#endif /* ARIADNE_SYSTEM_FUNCTION_H */
