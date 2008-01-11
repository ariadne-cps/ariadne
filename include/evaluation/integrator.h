/***************************************************************************
 *            integrator.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file integrator.h
 *  \brief An integrator based on computations using Taylor derivatives.
 */

#ifndef ARIADNE_INTEGRATOR_H
#define ARIADNE_INTEGRATOR_H

#include "base/types.h"
#include "numeric/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*!\ingroup Integrate
     * \brief A first order in space integration scheme.
     */
    template<class R>
    class IntegratorBase
    {
      typedef Numeric::Interval<R> I;
     public:

      /*! \brief Virtual destructor. */
      virtual ~IntegratorBase() { }

      /*! \brief Constructor. */
      IntegratorBase(smoothness_type temporal_order, smoothness_type spacial_order)
        : _temporal_order(temporal_order), _spacial_order() { }

      /*! \brief The order of the temporal model used. */
      smoothness_type temporal_order() const { return this->_temporal_order; }
      /*! \brief The order of the spacial model used. */
      smoothness_type spacial_order() const { return this->_spacial_order; }

    public:
      /*! \brief Compute an integration time and a bounding box, given a bounding box for the intitial set, and a maximum allowable flow time. */
      virtual 
      std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& f, 
                  const Geometry::Box<R>& bx,
                  const Numeric::Rational& t) const; 

      /*! \brief A model for the flow at time \a t with centre \a c, given that the flow remains in \a bb. */
      virtual 
      Function::AffineModel<R> 
      affine_flow_model(const System::VectorField<R>& vf,
                        const Geometry::Point<R>& c,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const;
     
      /*! \brief A model for the flow at time \a t with centre \a c, given that the flow remains in \a bb. */
      virtual 
      Function::TaylorModel<R> 
      taylor_flow_model(const System::VectorField<R>& vf,
                        const Geometry::Point<R>& c,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const;
      
     private:
      smoothness_type _temporal_order;
      smoothness_type _spacial_order;
    };
    
      
    
  }
}

#endif /* ARIADNE_INTEGRATOR_H */
