/***************************************************************************
 *            map.h
 *
 *  Wed Feb  2 18:33:10 2005
 *  Copyright  2005, 2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file map.h
 *  \brief MapInterface interface.
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


namespace Ariadne {
  namespace System {

    /*!\ingroup System
     * \ingroup DiscreteTime
     * \brief Abstract base class for (differentiable) functions.
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
    class MapInterface {
     protected:
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      /*! \brief The type obtained by applying the map to a state. */
      typedef Geometry::Point<F> result_type;
      
      /*! \brief Virtual destructor. */
      virtual ~MapInterface();
      
      /*! \brief Make a copy (clone) of the map. */
      virtual MapInterface<R>* clone() const = 0;

      /*! \brief An over-approximation to the image of a point. */
      Geometry::Point<F> operator() (const Geometry::Point<F>& pt) const {
        return this->image(pt); }
        
      /*! \brief An over-approximation to the image of a point. */
      virtual Geometry::Point<F> image(const Geometry::Point<F>& pt) const = 0;
      /*! \brief The derivative of the \a i th component with respect to the multi-index j. */
      virtual F derivative(const Geometry::Point<F>& r, const size_type& i, const Function::MultiIndex& j) const;
      /*! \brief The Jacobian derivative matrix over a rectangle. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& r) const;
        
      /*! \brief The degree of differentiability of the map. */
      virtual smoothness_type smoothness() const = 0;
      /*! \brief The dimension of the domain space. */
      virtual dimension_type argument_dimension() const = 0;
      /*! \brief The dimension of the range space. */
      virtual dimension_type result_dimension() const = 0;
    
      /*! \brief The name of the map. */
      virtual std::string name() const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
    };
   
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const MapInterface<R>& f) {
      return f.write(os);
    }  
    
  }
}

#endif /* ARIADNE_MAP_H */
