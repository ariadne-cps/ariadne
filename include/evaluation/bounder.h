/***************************************************************************
 *            bounder.h
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
 
/*! \file bounder.h
 *  \brief Class for bounding the flow of a vector field.
 */

#ifndef ARIADNE_BOUNDER_H
#define ARIADNE_BOUNDER_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "bounder_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for bounding the flow of a vector field.
     *  \ingroup VectorFieldEvolver
     */
    template<class R>
    class Bounder
      : public BounderInterface<R>
    {
      typedef Numeric::Interval<R> I;
     public:
      //@{ 
      //! \name Destructors, constructors and cloning operations.
 
      /*! \brief Default constructor. */
      Bounder();

      /*! \brief Make a dynamically-allocated copy. */
      Bounder<R>* clone() const;
     
      //@}


      //@{ 
      //! \name Methods for applying a system to a basic set.

      /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. 
       *
       *  This method may return \a false even if the flow remains in \a bound. 
       */
      virtual bool check_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                     const Geometry::Box<R>& initial_set,
                                     const Geometry::Box<R>& bound,
                                     const Numeric::Rational& integration_time) const;
      
      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Box<R> flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                 const Geometry::Box<R>& initial_set,
                                                 Numeric::Rational& integration_time) const;

      
      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Box<R> flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Box<R>& initial_set,
                                                          const Numeric::Rational& integration_time) const;

      /*! \brief Gives an inital estimated bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Box<R> estimate_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Box<R>& initial_set,
                                                          Numeric::Rational& integration_time) const;

      
      /*! \brief Computes an inital estimated bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Box<R> estimate_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Box<R>& initial_set,
                                                          const Numeric::Rational& integration_time) const;

      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Box<R> estimate_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Box<R>& initial_set,
                                                          const Numeric::Rational& integration_time,
                                                          const unsigned int& maximum_iterations) const;

      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Geometry::Box<R> refine_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                        const Geometry::Box<R>& initial_set,
                                                        const Geometry::Box<R>& estimated_bound,
                                                        const Numeric::Rational& integration_time) const;

      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting at \a initial_point remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Geometry::Box<R> refine_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                        const Geometry::Point<I>& initial_point,
                                                        const Geometry::Box<R>& estimated_bound,
                                                        const Numeric::Rational& integration_time) const;





      /*! \brief Compute a bound for the Jacobian of the flow over the time interval [-h,h], assuming that the flow remains inside the set \a b. */
      virtual LinearAlgebra::Matrix<I> estimate_flow_jacobian_bounds(const System::VectorFieldInterface<R>& vf,
                                                                     const Geometry::Box<R>& b,
                                                                     const Numeric::Rational& h) const;


      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Box<R> estimate_interval_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                                   const Geometry::Box<R>& initial_set,
                                                                   Numeric::Interval<R>& integration_time) const;

      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Box<R> refine_interval_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                                 const Geometry::Box<R>& initial_set,
                                                                 const Geometry::Box<R>& estimated_bound,
                                                                 const Numeric::Interval<R>& integration_time) const;

    };

  }
}

#endif /* ARIADNE_BOUNDER_H */
