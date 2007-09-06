/***************************************************************************
 *            curve.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file curve.h
 *  \brief A arbitraty curve in Euclidean space.
  */

#ifndef ARIADNE_CURVE_H
#define ARIADNE_CURVE_H

#include "../function/function_interface.h"
#include "curve_interface.h"

namespace Ariadne {
  namespace Geometry {
    
    //! \ingroup ExactSet
    /*! \brief A curve in Euclidean space
     */
    template<class R>
    class Curve
      : public DifferentiableCurveInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Destructor. */
      virtual ~Curve();
      /*! \brief Constructor. */
      Curve(const Function::DifferentiableFunctionInterface<R>& f);
      /*! \brief Copy constructor. */
      Curve(const Curve<R>& c);
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual Curve<R>* clone() const;
      /*! \brief The dimension of the set. */
      virtual dimension_type dimension() const;
      /*! \brief The smoothness of the curve. */
      virtual smoothness_type smoothness() const;

      /*! \brief The value at a point. */
      virtual Point<A> value(const A& s) const;
      /*! \brief The tangent at a point. */
      virtual LinearAlgebra::Vector<A> tangent(const A& s) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      Function::DifferentiableFunctionInterface<R>* _function_ptr;
    };
    
  }
}

#endif /* ARIADNE_CURVE_H */
