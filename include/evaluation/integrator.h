/***************************************************************************
 *            integrator.h
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
 
/*! \file integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef _ARIADNE_INTEGRATOR_H
#define _ARIADNE_INTEGRATOR_H

#include "../declarations.h"
#include "../numeric/rational.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*!\brief The type used to denote time. */
    typedef Numeric::Rational time_type;
    
    /*! \brief %Base class for integration schemes. 
     *  \ingroup Integrate
     */
    template<typename R>
    class Integrator {
     public:
      /*! \brief Virtual destructor. */
      virtual ~Integrator();

      /*! \brief Constructor. */
      Integrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. */
      virtual Geometry::GridMaskSet<R> integrate(const System::VectorField<R>& vector_field,
                                                 const Geometry::GridMaskSet<R>& initial_set,
                                                 const Geometry::GridMaskSet<R>& bounding_set,
                                                 const time_type& time) const = 0;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field, while remaining in \a bounding_set. */
      virtual Geometry::GridMaskSet<R> reach(const System::VectorField<R>& vector_field,
                                             const Geometry::GridMaskSet<R>& initial_set,
                                             const Geometry::GridMaskSet<R>& bounding_set,
                                             const time_type& time) const = 0;

      /*! \brief Integrate \a intial_set for all times under \a vector_field, while remaining in \a bounding_set. 
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual Geometry::GridMaskSet<R> chainreach(const System::VectorField<R>& vector_field,
                                                  const Geometry::GridMaskSet<R>& initial_set,
                                                  const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. 
       *
       *  This method may return \a false even if the flow remains in \a bound. 
       */
      virtual bool check_flow_bounds(const System::VectorField<R>& vector_field,
                                     const Geometry::Rectangle<R>& initial_set,
                                     const Geometry::Rectangle<R>& bound,
                                     const time_type& integration_time) const;
      
      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Rectangle<R> estimate_flow_bounds(const System::VectorField<R>& vector_field,
                                                          const Geometry::Rectangle<R>& initial_set,
                                                          time_type& integration_time) const;

      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Rectangle<R> estimate_flow_bounds(const System::VectorField<R>& vector_field,
                                                          const Geometry::Rectangle<R>& initial_set,
                                                          const time_type& integration_time,
                                                          const unsigned int& maximum_iterations) const;

      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_point remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Geometry::Rectangle<R> refine_flow_bounds(const System::VectorField<R>& vector_field,
                                                        const Geometry::Rectangle<R>& initial_set,
                                                        const Geometry::Rectangle<R>& estimated_bound,
                                                        const time_type& integration_time) const;

      /*! \brief A suggested minimum step size for integration. */
      virtual time_type minimum_step_size() const;
      /*! \brief The maximum allowable step size for integration. */
      virtual time_type maximum_step_size() const;
      /*! \brief A suggested minimum radius of a basic set after a subdivision (not a strict bound). */
      virtual R minimum_basic_set_radius() const;
      /*! T\brief he maximum allowable radius of a basic set during integration. */
      virtual R maximum_basic_set_radius() const;
      /*! \brief The time after which an integrator may approximate computed sets on a grid,
       *  in order to use previously caches integration results for the grid. */
      virtual time_type lock_to_grid_time() const;

     protected:
      /*! \brief Template for integrating a list set. */
      template<template<typename> class BS>
      Geometry::ListSet<R,BS> 
      integrate_list_set(const System::VectorField<R>& vector_field, 
                         const Geometry::ListSet<R,BS>& initial_set, 
                         const time_type& time) const;

      /*! \brief Template for integrating a basic set. */
      template<template<typename> class BS>
      BS<R>
      integrate_basic_set(const System::VectorField<R>& vector_field, 
                          const BS<R>& initial_set, 
                          const time_type& time) const;

      /*! \brief An algorithm for integrating forward a rectangle. */
      virtual Geometry::Rectangle<R> integration_step(const System::VectorField<R>&,
                                                      const Geometry::Rectangle<R>&,
                                                      time_type&) const = 0;

      /*! \brief An algorithm for integrating forward a parallelotope. */
      virtual Geometry::Parallelotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Parallelotope<R>&,
                                                          time_type&) const = 0;

       /*! \brief An algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<R> integration_step(const System::VectorField<R>&,
                                                     const Geometry::Zonotope<R>&,
                                                     time_type&) const = 0;
     protected:
      /*! \brief Constructs a zonotopic generator for the convex hull of \a iv and \a -iv. */ 
      static LinearAlgebra::Matrix<R> symmetrize(const LinearAlgebra::Vector< Interval<R> >& r);
     private:
      time_type _minimum_step_size;
      time_type _maximum_step_size;
      time_type _lock_to_grid_time;
      R _minimum_basic_set_radius;
      R _maximum_basic_set_radius;

    };

    /*! \brief %Base class for integration schemes which do not use the derivative of the vector field. 
     *  \ingroup Integrate
     */
    template<typename R>
    class C0Integrator : public Integrator<R> {
     public:
      /*! \brief Constructor. */
      C0Integrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      /*! \brief A C0 algorithm for integrating forward a rectangle. 
       */
      virtual Geometry::Rectangle<R> integrate(const System::VectorField<R>& vector_field,
                                               const Geometry::Rectangle<R>& initial_set,
                                               const time_type& time) const;

      /*! \brief A C0 algorithm for integrating forward a set of rectangles. 
       *
       * The algorithm first finds \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes \f$ R_{n+1} = R_{n}+ h_{n} Df(B_{n+1})\f$.
       *
       * C0 algorithms are typically too inaccurate for even the simplest systems.
       */
      virtual Geometry::ListSet<R,Geometry::Rectangle> integrate(const System::VectorField<R>& vector_field,
                                                                 const Geometry::ListSet<R,Geometry::Rectangle>& initial_set,
                                                                 const time_type& time) const;

      /*! \brief A C0 algorithm for integrating forward a set of rectangles up to a certain time. */
      virtual Geometry::ListSet<R,Geometry::Rectangle> reach(const System::VectorField<R>& vector_field,
                                                             const Geometry::ListSet<R,Geometry::Rectangle>& initial_set,
                                                             const time_type& time) const;
     public:
      /*! \brief A C0 algorithm for integrating forward a rectangle. */
      virtual Geometry::Rectangle<R> integration_step(const System::VectorField<R>& vector_field, 
                                                      const Geometry::Rectangle<R>& initial_set, 
                                                      time_type& step_size) const = 0;

      /*! \brief A C0 algorithm for integrating forward a parallelotope. */
      virtual Geometry::Parallelotope<R> integration_step(const System::VectorField<R>& vector_field, 
                                                          const Geometry::Parallelotope<R>& initial_set, 
                                                          time_type& step_size) const;

      /*! \brief A C0 algorithm for integrating forward a rectangle up to a certain time. */
      virtual Geometry::Rectangle<R> reachability_step(const System::VectorField<R>& vector_field, 
                                                       const Geometry::Rectangle<R>& initial_set,
                                                       time_type& step_size) const = 0;
    };

 
    /*! \brief %Base class for integration schemes which require at least a \f$C^1\f$ vector field.
     *  \ingroup Integrate
     */
    template<typename R>
    class C1Integrator : public Integrator<R> {
     public:
      /*! \brief Constructor. */
      C1Integrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      /*! \brief A C1 algorithm for integrating forward a rectangle.
       */
      virtual Geometry::Rectangle<R> integrate(const System::VectorField<R>&,
                                                   const Geometry::Rectangle<R>&,
                                                   const time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a parallelotope.
       */
      virtual Geometry::Parallelotope<R> integrate(const System::VectorField<R>&,
                                                   const Geometry::Parallelotope<R>&,
                                                   const time_type&) const;
      
      /*! \brief A C1 algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<R> integrate(const System::VectorField<R>&,
                                                   const Geometry::Zonotope<R>&,
                                                   const time_type&) const;
      
      /*! \brief A C1 algorithm for integrating forward a parallelotope.
       *
       * The algorithm first finds \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes an interval matrix \f$ \mathcal{A}_{n} \f$ such that \f$ Df(B_{n+1}) \in \mathcal{A}_{n} \f$.
       * It then computes a rectangle \f$ \mathcal{c}_{n+1} \f$ such that \f$ \Phi(t,c_{n})\in \mathcal{c}_{n+1} \f$.
       * We then compute \f$ \mathcal{P}_{n} \f$ such that \f$ D\Phi(h,R_{n}) \subset \mathcal{P}_{n} \f$.
       * We then compute \f$ A_{n+1} \f$ such that \f$ A_{n+1} \mathcal{e} \supset \mathcal{P}_{n} \mathcal{e} \f$.
       */
      virtual Geometry::ListSet<R,Geometry::Parallelotope> integrate(const System::VectorField<R>&,
                                                                     const Geometry::ListSet<R,Geometry::Parallelotope>&,
                                                                     const time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope.
       *
       * The algorithm first finds \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes an interval matrix \f$ \mathcal{A}_{n} \f$ such that \f$ Df(B_{n+1}) \in \mathcal{A}_{n} \f$.
       * It then computes a rectangle \f$ \mathcal{c}_{n+1} \f$ such that \f$ \Phi(t,c_{n})\in \mathcal{c}_{n+1} \f$.
       * We then compute \f$ \mathcal{P}_{n} \f$ such that \f$ D\Phi(h,R_{n}) \subset \mathcal{P}_{n} \f$.
       * We then compute \f$ A_{n+1} \f$ such that \f$ A_{n+1} \mathcal{e} \supset \mathcal{P}_{n} \mathcal{e} \f$.
       */
      virtual Geometry::ListSet<R,Geometry::Zonotope> integrate(const System::VectorField<R>&,
                                                                     const Geometry::ListSet<R,Geometry::Zonotope>&,
                                                                     const time_type&) const;

      
      /*! \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. */
      virtual Geometry::GridMaskSet<R> integrate(const System::VectorField<R>& vector_field,
                                                 const Geometry::GridMaskSet<R>& initial_set,
                                                 const Geometry::GridMaskSet<R>& bounding_set,
                                                 const time_type& time) const;


      /*! \brief A C1 algorithm for integrating forward a set of parallelotopes up to a certain time. */
      virtual Geometry::ListSet<R,Geometry::Parallelotope> reach(const System::VectorField<R>&,
                                                            const Geometry::ListSet<R,Geometry::Parallelotope>&,
                                                            const time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a set of zonotopes up to a certain time. */
      virtual Geometry::ListSet<R,Geometry::Zonotope> reach(const System::VectorField<R>&,
                                                            const Geometry::ListSet<R,Geometry::Zonotope>&,
                                                            const time_type&) const;

      /*! \brief A C1 algorithm for computing the reachable set of in initial set of Parallelotopes up to a certain time. */
      inline Geometry::ListSet<R,Geometry::Parallelotope> reach(
                      const System::VectorField<R>& vector_field, 
                      const Geometry::Parallelotope<R>& initial_set, 
                      const time_type& time) const
      {
        return this->reach(vector_field, initial_set.subdivide(), time);
      }

      /*! \brief A C1 algorithm for computing the reachable set of in initial set of Zonotopes up to a certain time. */
      inline Geometry::ListSet<R,Geometry::Zonotope> reach(
                      const System::VectorField<R>& vector_field,
                      const Geometry::Zonotope<R>& initial_set, 
                      const time_type& time) const
      {
        return this->reach(vector_field, initial_set.subdivide(), time);
      }
    
      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field, while remaining in \a bounding_set. */
      virtual Geometry::GridMaskSet<R> reach(const System::VectorField<R>& vector_field,
                                             const Geometry::GridMaskSet<R>& initial_set,
                                             const Geometry::GridMaskSet<R>& bounding_set,
                                             const time_type& time) const;

     public:
      /*! \brief A C1 algorithm for integrating forward a rectangle. */
      virtual Geometry::Rectangle<R> integration_step(const System::VectorField<R>&,
                                                      const Geometry::Rectangle<R>&,
                                                      time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a parallelotope. */
      virtual Geometry::Parallelotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Parallelotope<R>&,
                                                          time_type&) const = 0;

      /*! \brief A C1 algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Zonotope<R>&,
                                                          time_type&) const = 0;

      /*! \brief A C1 algorithm for integrating forward a rectangle up to a certain time. 
       */
      virtual Geometry::Rectangle<R> reachability_step(const System::VectorField<R>&,
                                                       const Geometry::Rectangle<R>&,
                                                       time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a parallelotope up to a certain time. 
       *
       * The result type is a zonotope, since the map is from \f$\mathbb{R}^n\times\mathbb{R}\rightarrow\mathbb{R}^n\f$.
       * and hence is not invertible. 
       *
       * This method forwards its argument to 
       * reachability_step(const System::VectorField<R>&, const Geometry::Zonotope<R>&, const time_type&)
       */
      virtual Geometry::Zonotope<R> reachability_step(const System::VectorField<R>&,
                                                      const Geometry::Parallelotope<R>&,
                                                      time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope up to a certain time. 
        */ 
      virtual Geometry::Zonotope<R> reachability_step(const System::VectorField<R>&,
                                                      const Geometry::Zonotope<R>&,
                                                      time_type&) const = 0;
    };






    
  }
}

#endif /* _ARIADNE_INTEGRATE_H */
