/***************************************************************************
 *            integrate.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file integrate.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef _ARIADNE_INTEGRATE_H
#define _ARIADNE_INTEGRATE_H

#include "../geometry/geometry_declarations.h"
#include "../linear_algebra/linear_algebra_declarations.h"

#include "../evaluation/evaluation_declarations.h"
#include "../evaluation/vector_field.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
    template<typename R>
    bool
    check_flow_bounds(const Evaluation::VectorField<R>& vector_field,
                      const Geometry::Rectangle<R>& initial_set,
                      const Geometry::Rectangle<R>& bound,
                      const R& integration_time);

    /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to \a integration_time. 
     *
     *  The algorithm is guarenteed to terminate if the flow has linear growth at infinity.
     */
    template<typename R>
    Geometry::Rectangle<R>
    estimate_flow_bounds(const Evaluation::VectorField<R>& vector_field,
                         const Geometry::Rectangle<R>& initial_set,
                         R& integration_time);

    /*! \brief Compute a set \a bound such that the flow of \a vector_field starting at \a initial_point remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
    template<typename R>
    Geometry::Rectangle<R>
    refine_flow_bounds(const Evaluation::VectorField<R>& vector_field,
                       const Geometry::Point<R>& initial_point,
                       const Geometry::Rectangle<R>& estimated_bound,
                       const R& integration_time);

    
    /*! \brief An inaccurate C0 algorithm for integrating forward a rectangle. */
    template<typename R>
    Geometry::Rectangle<R> 
    integration_step(const VectorField<R>& vector_field, const Geometry::Rectangle<R>& initial_set, R& step_size);

    /*! \brief A C1 algorithm for integrating forward a parallelotope. */
    template<typename R>
    Geometry::Parallelotope<R> 
    integration_step(const VectorField<R>& vector_field, const Geometry::Parallelotope<R>& initial_set, R& step_size);

    /*! \brief A specialized algorithm for integrating forward a parallelotope. */
    template<typename R>
    Geometry::Parallelotope<R> 
    integration_step(const AffineVectorField<R>& vector_field, const Geometry::Parallelotope<R>& initial_set, R& step_size);

    
    /*! \brief An inaccurate C0 algorithm for integrating forward a rectangle up to a certain time. */
    template<typename R>
    Geometry::Rectangle<R> 
    reach_step(const VectorField<R>& vector_field, const Geometry::Rectangle<R>& initial_set, R& step_size);

    /*! \brief A C1 algorithm for integrating forward a parallelotope up to a certain time. */
    template<typename R>
    Geometry::Parallelotope<R> 
    reach_step(const VectorField<R>& vector_field, const Geometry::Parallelotope<R>& initial_set, R& step_size);

    /*! \brief A C1 algorithm for integrating forward a zonotope up to a certain time. */
    template<typename R>
    Geometry::Zonotope<R> 
    reach_step(const VectorField<R>& vector_field, const Geometry::Zonotope<R>& initial_set, R& step_size);



    /*! \brief Integrate a rectangle. */
    template<typename R>
    Geometry::Rectangle<R> 
    integrate(const VectorField<R>& vector_field, const Geometry::Rectangle<R>& initial_set, const R& time, const R& step_size);

    /*! \brief Integrate a parallelotope. */
    template<typename R>
    Geometry::Parallelotope<R> 
    integrate(const VectorField<R>& vector_field, const Geometry::Parallelotope<R>& initial_set, const R& time, const R& step_size);

    /*! \brief Integrate a set. */
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    integrate(const VectorField<R>& vector_field, const Geometry::ListSet<R,BS>& initial_set, const R& time, const R& step_size);

    /*! \brief Integrate a set using a C1 algorithm. */
    template<typename R>
    Geometry::GridMaskSet<R> 
    integrate(const VectorField<R>& vector_field, const Geometry::GridMaskSet<R>& initial_set, const R& time, const R& step_size);

    
    /*! \brief Compute the reachable set. */
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    reach(const VectorField<R>& vector_field, const Geometry::ListSet<R,BS>& initial_set, const R& time, const R& step_size);

    /*! \brief Integrate a set using a C1 algorithm. */
    template<typename R>
    Geometry::GridMaskSet<R> 
    reach(const VectorField<R>& vector_field, const Geometry::GridMaskSet<R>& initial_set, const R& time, const R& step_size);

    
/*
    template<typename R>
    Ariadne::Geometry::GridMaskSet<R> 
    chainreach(const VectorField<R>& f, 
               const Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>& is, 
               const Ariadne::Geometry::FiniteGrid<R>& g, 
               const Ariadne::Geometry::Rectangle<R>& bb);
*/
    
  }
}

#endif /* _ARIADNE_INTEGRATE_H */
