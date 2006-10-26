/***************************************************************************
 *            vector_field.h
 *
 *  Thu Feb  3 21:06:54 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
/*! \file vector_field.h
 *  \brief Vector_type field interface.
 */
 
#ifndef _ARIADNE_VECTOR_FIELD_H
#define _ARIADNE_VECTOR_FIELD_H

#include <limits>

#include "../declarations.h"


namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup ContinuousTime
     * \brief Abstract base class for (differentiable) vector fields.
     * 
     * The system is specified by the method operator()(const Geometry::Rectangle<R>& A) const,
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
     * The method derivative(const Geometry::Rectangle<R>& A) const computes the \a i th component of the derivative over the set \a A 
     * with respect to the variables in the multi-index \a j.
     */
    template<class R>
    class VectorField {
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Virtual destructor. */
      virtual ~VectorField();
     
      /*! \brief An approximation to the vector field at a point. */
      LinearAlgebra::Vector<F> operator() (const Geometry::Point<R>& x) const { return this->image(x); }
      /*! \brief A bound for the vector field over a rectangle. */
      LinearAlgebra::Vector<I> operator() (const Geometry::Rectangle<R>& A) const { return this->image(A); };

      /*! \brief An approximation to the vector field at a point. */
      virtual LinearAlgebra::Vector<F> image(const Geometry::Point<R>& x) const;
      /*! \brief A bound for the vector field over a rectangle. */
      virtual LinearAlgebra::Vector<I> image(const Geometry::Rectangle<R>& A) const;

      /*! \brief An approximation to the vector field at a point. */
      virtual F derivative(const Geometry::Point<R>& x, const size_type& i, const multi_index_type& j) const;
      /*! \brief A bound for the vector field over a rectangle. */
      virtual I derivative(const Geometry::Rectangle<R>& A, const size_type& i, const multi_index_type& j) const;

      /*! \brief An approximation to the Jacobian derivative at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<R>& x) const;
      /*! \brief A bound for the Jacobian derivative over a rectangle. */
      virtual LinearAlgebra::Matrix<I> jacobian(const Geometry::Rectangle<R>& A) const;
    
      /*! \brief The degree of differentiability of the map. */
      virtual size_type smoothness() const = 0;
      /*! \brief The dimension of the space the vector field lives in. */
      virtual dimension_type dimension() const = 0;

      /*! \brief The name of the system. */
      virtual std::string name() const = 0;
    };
   
  }
}

#endif /* _ARIADNE_VECTOR_FIELD_H */
