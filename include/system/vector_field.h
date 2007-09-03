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
 
#ifndef ARIADNE_VECTOR_FIELD_H
#define ARIADNE_VECTOR_FIELD_H

#include <limits>

#include "../base/types.h"
#include "../numeric/declarations.h"
#include "../numeric/numerical_traits.h"
#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"


namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup ContinuousTime
     * \brief Abstract base class for (differentiable) vector fields.
     * 
     * The system is specified by the method operator()(const Geometry::Point<F>& pt) const,
     * This method should compute an interval vector \f$v=\overline{f}(A)\f$ with the
     * following properties:
     *   -# \f$f(p)\subset\overline{f}(A)\f$,
     *   -# If \f$A_1\subset A_0\f$, then \f$\overline{f}(A_1)\subset 
     *       \overline{f}(A_0)\f$, and
     *   -# If \f$\bigcap_{n\in\mathbb{N}}A_n=\{x\}\f$, then 
     *       \f$\bigcap_{n\in\mathbb{N}}\overline{f}(A_n)=\{f(x)\}\f$.
     *
     * More succinctly, we say that \f$\overline{f}(A_n)\f$ converges monotonically 
     * as \f$A_n\f$ tends to a point.
     *
     * Additional accuracy can be obtained be using derivatives.
     * The method jacobian(const Geometry::Point<F>& pt) const computes the derivative matrix at/over the point \a pt.
     */
    template<class R>
    class VectorFieldInterface {
     protected:
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Virtual destructor. */
      virtual ~VectorFieldInterface();
     
      /*! \brief Make a copy (clone) of the vector field. */
      virtual VectorFieldInterface<R>* clone() const = 0;
     
      /*! \brief An approximation to the vector field at a point. */
      LinearAlgebra::Vector<F> operator() (const Geometry::Point<F>& x) const { return this->image(x); }

      /*! \brief An approximation to the vector field at a point. */
      virtual LinearAlgebra::Vector<F> image(const Geometry::Point<F>& x) const = 0;

      /*! \brief An approximation to the vector field at a point. */
      virtual F derivative(const Geometry::Point<F>& x, const size_type& i, const LinearAlgebra::MultiIndex& j) const;

      /*! \brief An approximation to the Jacobian derivative at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& x) const;
    
      /*! \brief The degree of differentiability of the map. */
      virtual smoothness_type smoothness() const = 0;
      /*! \brief The dimension of the space the vector field lives in. */
      virtual dimension_type dimension() const = 0;

      /*! \brief The name of the system. */
      virtual std::string name() const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
    };
   
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const VectorFieldInterface<R>& vf) {
      return vf.write(os);
    }
    
  }
}

#endif /* ARIADNE_VECTOR_FIELD_H */
