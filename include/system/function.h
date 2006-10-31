/***************************************************************************
 *            function.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file function.h
 *  \brief General function interface.
 */
 
#ifndef _ARIADNE_FUNCTION_H
#define _ARIADNE_FUNCTION_H

#include "../declarations.h"

namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Abstract base class for (differentiable) functionss.
     * 
     * The function is specified by the method operator()(const LinearAlgebra::Vector< Interval<R> >& A) const,
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
     * The method derivative(const LinearAlgebra::Vector< Interval<R> >& A) const computes the \a i th component of the derivative over the set \a A 
     * with respect to the variables in the multi-index \a j.
     */
    template<class R>
    class Function {
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      
      /*! \brief Virtual destructor. */
      virtual ~Function() { }
     
      /*! \brief Make a copy (clone) of the vector field. */
      virtual Function<R>* clone() const = 0;
     
      /*! \brief A bound for the function over a set of vectors. */
      LinearAlgebra::Vector<I> operator() (const LinearAlgebra::Vector<I>& A) const { return this->image(A); };

      /*! \brief A bound for the vector field over aa set of vectors. */
      virtual LinearAlgebra::Vector<I> image(const LinearAlgebra::Vector<I>& A) const = 0;

      /*! \brief A bound for the vector field over a set of vectors. */
      virtual I derivative(const LinearAlgebra::Vector<I>& A, const size_type& i, const multi_index_type& j) const = 0;

    
      /*! \brief The degree of differentiability of the function. */
      virtual size_type smoothness() const = 0;
      /*! \brief The dimension of the function argument. */
      virtual dimension_type argument_dimension() const = 0;
      /*! \brief The dimension of the function result. */
      virtual dimension_type result_dimension() const = 0;
    };
   
  }
}

#endif /* _ARIADNE_FUNCTION_H */
