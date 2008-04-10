/***************************************************************************
 *            standard_flower.h
 *
 *  Copyright  2006-8  Pieter Collins
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
 

#ifndef ARIADNE_STANDARD_FLOWER_H
#define ARIADNE_STANDARD_FLOWER_H

#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "system/vector_field.h"

#include "evaluation/flower_interface.h"



namespace Ariadne { 
  namespace Evaluation { 

    /*!\ingroup Integrators
     * \brief A class for integrating a vector field to obtain a higher-order models for the flow.
     *
     * To compute an enclosure for the solution \f$\xi\f$, 
     * we use the criterion that if \f$B\f$ is the initial condition and \f$D\f$ is a bound, then
     *   \f[ B + [0,h] f(D) \subset D . \f]
     * To compute an enclosure for the solution \f$\xi\f$ and its spacial derivatives up to order \a n, we use the same formula, 
     * but on the variation of the flow. Explicitly, if \f$W\f$ is to be a bound on the first variation \a V, we have
     *   \f[ I + [0,h] Df(D) W \subset W . \f$
     *
     */
    template<class R>
    class StandardFlower
      : public FlowerInterface<R>
    {
      typedef Numeric::Interval<R> I;
     public:

      /*! \brief Virtual destructor. */
      virtual ~StandardFlower() { }

      /*! \brief Constructor. */
      StandardFlower(smoothness_type temporal_order, smoothness_type spacial_order)
        : _temporal_order(temporal_order), _spacial_order(spacial_order) { }
      /*! \brief Cloning operator. */
      StandardFlower<R>* clone() const { return new StandardFlower<R>(*this); }

      /*! \brief The order of the temporal model used. */
      smoothness_type temporal_order() const { return this->_temporal_order; }
      /*! \brief The order of the spacial model used. */
      smoothness_type spacial_order() const { return this->_spacial_order; }

    public:
      /*! \brief A model for the flow at time \a t with centre \a c, given that the flow remains in \a bb. */
      virtual 
      Function::AffineModel<R> 
      affine_flow_model(const System::VectorField<R>& vf,
                        const Geometry::Point<R>& c,
                        const Geometry::Box<R>& d,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const;
     
      /*! \brief A model for the flow at time \a t with centre \a c, given that the flow remains in \a bb. */
      virtual 
      Function::TaylorModel<R> 
      taylor_flow_model(const System::VectorField<R>& vf,
                        const Geometry::Point<R>& c,
                        const Geometry::Box<R>& d,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const;
      
      /*! \brief Evolve the flow for a step of time \a step_size, given that the flow remains in \a bounding_box. */
      virtual 
      Geometry::Point<I>
      flow_step(const System::VectorField<R>& vector_field,
                const Geometry::Point<I>& initial_point,
                const Numeric::Rational& step_size,
                const Geometry::Box<R>& bounding_box) const;
      
      /*! \brief Compute the Jacobian derivative (first variation) of the flow. */
      virtual
      LinearAlgebra::Matrix< Numeric::Interval<R> >
      flow_step_jacobian(const System::VectorField<R>& vector_field, 
                         const Geometry::Point<I>& initial_point, 
                         const Numeric::Rational& step_size, 
                         const Geometry::Box<R>& bounding_box) const;
     private:
      smoothness_type _temporal_order;
      smoothness_type _spacial_order;
    };
    



  }
}

#endif /* STANDARD_FLOWER_H */
