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
  

    /*! \ingroup Integrators
     *  \brief A class for bounding the flow of a vector field.
     */
    template<class R>
    class StandardBounder
      : public BounderInterface<R>
    {
      typedef Interval<R> I;
     public:
      //@{ 
      //! \name Destructors, constructors and cloning operations.
 
      /*! \brief Default constructor. */
      StandardBounder();

      /*! \brief Construct with a given maximum step size. */
      StandardBounder(const Rational& max_step_size);

      /*! \brief Make a dynamically-allocated copy. */
      StandardBounder<R>* clone() const;
     
      //@}

      //@{ 
      //! \name Data access.
      Rational maximum_step_size() const;
      //@}

      //@{ 
      //! \name Methods for computing bounding set for a flow. */

      /*! \brief Computes an integration time and a bounding box for the flow of \a vector_field starting in \a initial_set. */
      virtual std::pair< Rational, Box<R> >
      flow_bounds(const VectorField<R>& vector_field,
                  const Box<R>& initial_set) const;

      /*! \brief Computes an integration time and a bounding box for the flow of \a vector_field starting in \a initial_set. */
      virtual std::pair< Rational, Box<R> >
      flow_bounds(const VectorField<R>& vector_field,
                  const Box<R>& initial_set,
                  const Rational& maximum_step_size) const;

      /*! \brief Compute an integration time and bounds for both the flow and its derivatives up to order \a o. */
      virtual 
      std::pair< Rational, TaylorDerivative<I> >
      variation_flow_bounds(const VectorField<R>& f, 
                            const Box<R>& bx,
                            const Rational& t,
                            smoothness_type o) const; 

     /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. 
       *
       *  This method may return \a false even if the flow remains in \a bound. 
       */
      virtual bool check_flow_bounds(const VectorField<R>& vector_field,
                                     const Box<R>& initial_set,
                                     const Rational& integration_time,
                                     const Box<R>& bound) const;
      
      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Box<R> refine_flow_bounds(const VectorField<R>& vector_field,
                                                        const Box<R>& initial_set,
                                                        const Box<R>& estimated_bound,
                                                        const Rational& integration_time) const;



      /*! \brief Compute a bound for the Jacobian of the flow over the time interval [-h,h], assuming that the flow remains inside the set \a b. */
      virtual Matrix<I> estimate_flow_jacobian_bounds(const VectorField<R>& vf,
                                                                     const Box<R>& b,
                                                                     const Rational& h) const;


     private:
      Rational _maximum_step_size;
    };





  
} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_H */
