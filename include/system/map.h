/***************************************************************************
 *            map.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_MAP_H
#define ARIADNE_MAP_H

#include <string>

#include "base/types.h"
#include "numeric/declarations.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"

#include <boost/shared_ptr.hpp>
#include "geometry/euclidean_space.h"

/*! \file map.h
 *  \brief MapInterface interface.
 */



namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief A discrete-time dynamical system, defined by a function.
     * 
     * The system is specified by the method operator()(const Geometry::Point<F>& A) const.
     * This method should compute a basic set \f$\overline{f}(A)\f$ with the
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
     * The method derivative(const Geometry::Point<F>& A) const computes the \a i th component of the derivative over the set \a A 
     * with respect to the variables in the multi-index \a j.
     */
    template<class R>
    class Map
    {
     protected:
      typedef typename Numeric::traits<R>::arithmetic_type A; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The type used to represent time. */
      typedef Numeric::Integer time_type;
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type used to describe the state space. */
      typedef Geometry::EuclideanSpace state_space_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      /*! \brief The type obtained by applying the map to a state. */
      typedef Geometry::Point<A> result_type;
      
      /*! \brief Construct from a function interface. */
      Map(const Function::FunctionInterface<R>& f);
      /*! \brief Construct a dynamically-allocated copy. */
      Map<R>* clone() const { return new Map<R>(*this); }
      /*! \brief The defining function. */
      const Function::FunctionInterface<R>& function() const;

      /*! \brief An over-approximation to the image of a point. */
      Geometry::Point<A> operator() (const Geometry::Point<A>& pt) const;
      /*! \brief An over-approximation to the image of a point. */
      Geometry::Point<A> image(const Geometry::Point<A>& pt) const;
      /*! \brief The Jacobian derivative matrix over a rectangle. */
      LinearAlgebra::Matrix<A> jacobian(const Geometry::Point<A>& x) const;
      /*! \brief The derivatives up to order \a s. */
      Function::TaylorDerivative<A> derivative(const Geometry::Point<A>& x, const smoothness_type& s) const;
        
      /*! \brief The state space of a self-map. */
      Geometry::EuclideanSpace state_space() const;
      /*! \brief The dimension of the space of a self-map. */
      dimension_type dimension() const;
      /*! \brief The degree of differentiability of the map. */
      smoothness_type smoothness() const;
      /*! \brief The dimension of the domain space. */
      dimension_type argument_dimension() const;
      /*! \brief The dimension of the range space. */
      dimension_type result_dimension() const;
    
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< Function::FunctionInterface<R> > _function_ptr;
    };
   
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Map<R>& f) {
      return f.write(os); }

  }
}

#endif /* ARIADNE_MAP_H */
