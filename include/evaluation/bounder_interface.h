/***************************************************************************
 *            bounder_interface.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file bounder_interface.h
 *  \brief Class for bounding the flow of a vector field.
 */

#ifndef ARIADNE_BOUNDER_INTERFACE_H
#define ARIADNE_BOUNDER_INTERFACE_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for bounding the flow of a vector field.
     *  \ingroup VectorFieldEvolver
     */
    template<class R>
    class BounderInterface
    {
      typedef Numeric::Interval<R> I;
     public:
      //@{ 
      //! \name Destructors, constructors and cloning operations.
      /*! \brief Destructor. */
      virtual ~BounderInterface() { }

      /*! \brief Make a dynamically-allocated copy. */
      virtual BounderInterface<R>* clone() const = 0;
     
      //@}


      //@{ 
      //! \name Methods for computing bounding boxes for a flow.

      /*! \brief Computes a step size and bounding box for the flow of \a vector_field starting in \a initial_set remains in the bound for times up to the step size. The step size must be less than \a maximum_step_size. */
      virtual std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& vector_field,
                  const Geometry::Box<R>& initial_set,
                  const Numeric::Rational& maximum_step_size) const = 0;

      /*! \brief Computes a flow time and bounding box for the flow of \a vector_field starting in \a initial_set.  */
      virtual 
      std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& vector_field,
                  const Geometry::Box<R>& initial_set) const = 0;

      



      /*! \brief Compute a bound for the Jacobian of the flow over the time interval [-h,h], assuming that the flow remains inside the set \a b. */
      virtual LinearAlgebra::Matrix<I> estimate_flow_jacobian_bounds(const System::VectorField<R>& vf,
                                                                     const Geometry::Box<R>& b,
                                                                     const Numeric::Rational& h) const = 0;


      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Box<R> estimate_interval_flow_bounds(const System::VectorField<R>& vector_field,
                                                                   const Geometry::Box<R>& initial_set,
                                                                   Numeric::Interval<R>& integration_time) const = 0;

      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Box<R> refine_interval_flow_bounds(const System::VectorField<R>& vector_field,
                                                                 const Geometry::Box<R>& initial_set,
                                                                 const Geometry::Box<R>& estimated_bound,
                                                                 const Numeric::Interval<R>& integration_time) const = 0;

    };

  }
}

#endif /* ARIADNE_BOUNDER_INTERFACE_H */
