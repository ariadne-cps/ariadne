/***************************************************************************
 *            standard_bounder.h
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
 
/*! \file standard_bounder.h
 *  \brief Class for bounding the flow of a vector field.
 */

#ifndef ARIADNE_STANDARD_BOUNDER_H
#define ARIADNE_STANDARD_BOUNDER_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "bounder_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \ingroup Integrators
     *  \brief A class for bounding the flow of a vector field.
     */
    template<class R>
    class StandardBounder
      : public BounderInterface<R>
    {
      typedef Numeric::Interval<R> I;
     public:
      //@{ 
      //! \name Destructors, constructors and cloning operations.
 
      /*! \brief Default constructor. */
      StandardBounder();

      /*! \brief Construct with a given maximum step size. */
      StandardBounder(const Numeric::Rational& max_step_size);

      /*! \brief Make a dynamically-allocated copy. */
      StandardBounder<R>* clone() const;
     
      //@}

      //@{ 
      //! \name Data access.
      Numeric::Rational maximum_step_size() const;
      //@}

      //@{ 
      //! \name Methods for computing bounding set for a flow. */

      /*! \brief Computes an integration time and a bounding box for the flow of \a vector_field starting in \a initial_set. */
      virtual std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& vector_field,
                  const Geometry::Box<R>& initial_set) const;

      /*! \brief Computes an integration time and a bounding box for the flow of \a vector_field starting in \a initial_set. */
      virtual std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& vector_field,
                  const Geometry::Box<R>& initial_set,
                  const Numeric::Rational& maximum_step_size) const;

      /*! \brief Compute an integration time and bounds for both the flow and its derivatives up to order \a o. */
      virtual 
      std::pair< Numeric::Rational, Function::TaylorDerivative<I> >
      variation_flow_bounds(const System::VectorField<R>& f, 
                            const Geometry::Box<R>& bx,
                            const Numeric::Rational& t,
                            smoothness_type o) const; 

     /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. 
       *
       *  This method may return \a false even if the flow remains in \a bound. 
       */
      virtual bool check_flow_bounds(const System::VectorField<R>& vector_field,
                                     const Geometry::Box<R>& initial_set,
                                     const Numeric::Rational& integration_time,
                                     const Geometry::Box<R>& bound) const;
      
      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Geometry::Box<R> refine_flow_bounds(const System::VectorField<R>& vector_field,
                                                        const Geometry::Box<R>& initial_set,
                                                        const Geometry::Box<R>& estimated_bound,
                                                        const Numeric::Rational& integration_time) const;



      /*! \brief Compute a bound for the Jacobian of the flow over the time interval [-h,h], assuming that the flow remains inside the set \a b. */
      virtual LinearAlgebra::Matrix<I> estimate_flow_jacobian_bounds(const System::VectorField<R>& vf,
                                                                     const Geometry::Box<R>& b,
                                                                     const Numeric::Rational& h) const;


     private:
      Numeric::Rational _maximum_step_size;
    };





  }
}

#endif /* ARIADNE_BOUNDER_H */
