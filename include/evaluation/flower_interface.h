/***************************************************************************
 *            flower_interface.h
 *
 *  Copyright  2007-8 Pieter Collins
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
 
/*! \file flower_interface.h
 *  \brief Class for computing the image of points under a vector field. 
 */

#ifndef ARIADNE_FLOWER_INTERFACE_H
#define ARIADNE_FLOWER_INTERFACE_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

namespace Ariadne {
  namespace Evaluation {


     /*! \brief Interface for computing the flow of a vector field and its spacial variations.
     *   \ingroup EvaluatorInterfaces \ingroup Integrators
     */
    template<class R>
    class FlowerInterface
    {
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Virtual destructor. */
      virtual ~FlowerInterface() { }

      /*! \brief Make a dynamically-allocated copy. */
      virtual FlowerInterface<R>* clone() const = 0;

      virtual 
      /*! \brief A model for the flow at time \a t with centre \a c, given that the flow remains in \a bb. */
      Function::AffineModel<R> 
      affine_flow_model(const System::VectorField<R>& vf,
                        const Geometry::Point<R>& c,
                        const Geometry::Box<R>& d,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const = 0;
     
      /*! \brief A model for the flow at time \a t with centre \a c, given that the flow remains in \a bb. */
      virtual 
      Function::TaylorModel<R> 
      taylor_flow_model(const System::VectorField<R>& vf,
                        const Geometry::Point<R>& c,
                        const Geometry::Box<R>& d,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const = 0;
      
       /*! \brief Compute the flow of a point. */
      virtual 
      Geometry::Point<I> 
      flow_step(const System::VectorField<R>& f, 
                const Geometry::Point<I>& s, 
                const Numeric::Rational& t, 
                const Geometry::Box<R>& bb) const = 0;

      /*! \brief Compute the spacial jacobian over a flow step of time \a t starting at \a p assuming that the flow remains within \a bb. */
      virtual LinearAlgebra::Matrix<I> flow_step_jacobian(const System::VectorField<R>& vf,
                                                          const Geometry::Point<I>& p,
                                                          const Numeric::Rational& t,
                                                          const Geometry::Box<R>& bb) const = 0;
    };


  }
}

#endif /* ARIADNE_FLOWER_INTERFACE_H */
