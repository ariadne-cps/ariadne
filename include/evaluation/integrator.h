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

#ifndef ARIADNE_INTEGRATOR_H
#define ARIADNE_INTEGRATOR_H

#include "../base/types.h"
#include "../base/declarations.h"
#include "../numeric/declarations.h"
#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

namespace Ariadne {

  namespace Geometry {
    template<class R> class SetInterface;
    template<class R> class ArbitrarySet;
    template<class BS> class ListSet;
    template<class R> class GridMaskSet;
  }

  namespace Evaluation {

   
    /*! \brief The identity function. */
    template<class BS> inline
    BS identity(const BS& bs) 
    { 
      return bs;
    }

  

   



    /*! \brief %Base class for integration schemes. 
     *  \ingroup Integrate
     */
    template<class R>
    class Integrator {
     private:
      typedef Numeric::Interval<R> I;
      time_type _minimum_step_size;
      time_type _maximum_step_size;
      time_type _lock_to_grid_time;
      R _grid_size;
      R _minimum_basic_set_radius;
      R _maximum_basic_set_radius;
     public:
      /*! The type used to denote real numbers. */
      typedef R real_type;

      /*! The type used to denote an interval of real numbers. */
      typedef Numeric::Interval<R> interval_type;
     
      /*! \brief Virtual destructor. */
      virtual ~Integrator();

      /*! \brief Make a dynamically-allocated copy. */
      virtual Integrator<R>* clone() const = 0;

      /*! \brief Constructor. */
      Integrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      //@{
      //! \name Parameters controlling the accuracy

      /*! \brief A suggested minimum step size for integration. */
      virtual time_type minimum_step_size() const;
      /*! \brief The maximum allowable step size for integration. */
      virtual time_type maximum_step_size() const;
      /*! \brief A suggested minimum radius of a basic set after a subdivision (not a strict bound). */
      virtual R minimum_basic_set_radius() const;
      /*! \brief The maximum allowable radius of a basic set during integration. */
      virtual R maximum_basic_set_radius() const;
      /*! \brief The time after which an integrator may approximate computed sets on a grid,
       *  in order to use previously caches integration results for the grid. */
      virtual time_type lock_to_grid_time() const;
      /*! \brief The time after which an integrator may approximate computed sets on a grid,
       *  in order to use previously caches integration results for the grid. */
      virtual R grid_size() const;
      /*! \brief Set the grid size. */
      virtual void set_grid_size(const R&);
      //@}

      //@{
      //! \name Integration routines

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* integrate(const System::VectorFieldInterface<R>& vector_field,
                                                   const Geometry::SetInterface<R>& initial_set,
                                                   const time_type& time) const = 0;

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* integrate(const System::VectorFieldInterface<R>& vector_field,
                                                   const Geometry::SetInterface<R>& initial_set,
                                                   const Geometry::SetInterface<R>& bounding_set,
                                                   const time_type& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const time_type& time) const = 0;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Geometry::SetInterface<R>& bounding_set,
                                               const time_type& time) const;

      /*! \brief Integrate \a intial_set for all times under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set) const = 0;

      /*! \brief Integrate \a intial_set for all times under \a vector_field, while remaining in \a bounding_set. 
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual Geometry::SetInterface<R>* chainreach(const System::VectorFieldInterface<R>& vector_field,
                                                   const Geometry::SetInterface<R>& initial_set,
                                                   const Geometry::SetInterface<R>& bounding_set) const;

      /*! \brief Computes the set of points which remain in \a bounding_set under evolution of \a vector_field.
       */
      virtual Geometry::SetInterface<R>* viable(const System::VectorFieldInterface<R>& vector_field,
                                                const Geometry::SetInterface<R>& bounding_set) const;

      /*! \brief  Verifies that the flow of \a vector_field starting in \a initial_set remains in \a safe_set all times.
       */
      virtual tribool verify(const System::VectorFieldInterface<R>& vector_field,
                             const Geometry::SetInterface<R>& initial_set,
                             const Geometry::SetInterface<R>& safe_set) const;

      //@}

      //@{
      //! \name Integration of list and grid sets. (Deprecated)


      /*!  \deprecated \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. (Deprecated) 
       */
      virtual Geometry::ListSet< Geometry::Rectangle<R> > integrate(const System::VectorFieldInterface<R>& vector_field,
                                                                    const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
                                                                    const time_type& time) const = 0;

      /*!  \deprecated \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. (Deprecated) 
       */
      virtual Geometry::ListSet< Geometry::Rectangle<R> > reach(const System::VectorFieldInterface<R>& vector_field,
                                                                const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
                                                                const time_type& time) const = 0;

      /*!  \deprecated \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. (Deprecated) 
       */
      virtual Geometry::GridMaskSet<R> integrate(const System::VectorFieldInterface<R>& vector_field,
                                                  const Geometry::GridMaskSet<R>& initial_set,
                                                  const Geometry::GridMaskSet<R>& bounding_set,
                                                  const time_type& time) const = 0;

      /*! \deprecated \brief Integrate \a intial_set for times up to \a time under \a vector_field, while remaining in \a bounding_set. (Deprecated) 
       */
      virtual Geometry::GridMaskSet<R> reach(const System::VectorFieldInterface<R>& vector_field,
                                             const Geometry::GridMaskSet<R>& initial_set,
                                             const Geometry::GridMaskSet<R>& bounding_set,
                                             const time_type& time) const = 0;

      /*! \deprecated \brief Integrate \a intial_set for all times under \a vector_field, while remaining in \a bounding_set. (Deprecated)
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual Geometry::GridMaskSet<R> chainreach(const System::VectorFieldInterface<R>& vector_field,
                                                   const Geometry::GridMaskSet<R>& initial_set,
                                                   const Geometry::GridMaskSet<R>& bounding_set) const = 0;

      /*! \brief Computes the set of points which remain in \a bounding_set under evolution of \a vector_field.
       */
      virtual Geometry::GridMaskSet<R> viable(const System::VectorFieldInterface<R>& vector_field,
                                              const Geometry::GridMaskSet<R>& bounding_set) const = 0;

      /*! \brief  Verifies that the flow of \a vector_field starting in \a initial_set remains in \a safe_set all times. (Deprecated) \deprecated
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual tribool verify(const System::VectorFieldInterface<R>& vector_field,
                             const Geometry::GridMaskSet<R>& initial_set,
                             const Geometry::GridMaskSet<R>& safe_set) const = 0;


      //@}


      /*! \brief Verifies that the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. 
       *
       *  This method may return \a false even if the flow remains in \a bound. 
       */
      virtual bool check_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                     const Geometry::Rectangle<R>& initial_set,
                                     const Geometry::Rectangle<R>& bound,
                                     const time_type& integration_time) const;
      
      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. The integration time may be dynamically varied to allow the bounding box to be computed. */
      virtual Geometry::Rectangle<R> estimate_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Rectangle<R>& initial_set,
                                                          time_type& integration_time) const;

      
      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Rectangle<R> estimate_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Rectangle<R>& initial_set,
                                                          const time_type& integration_time) const;

      /*! \brief Computes a bounding box for the flow of \a vector_field starting in \a initial_set remains in \a bound for times up to time \a integration_time. */
      virtual Geometry::Rectangle<R> estimate_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                          const Geometry::Rectangle<R>& initial_set,
                                                          const time_type& integration_time,
                                                          const unsigned int& maximum_iterations) const;

      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_point remains in \a bound for times up to \a integration_time, given a bound \a estimated_bound. */
      virtual Geometry::Rectangle<R> refine_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                        const Geometry::Rectangle<R>& initial_set,
                                                        const Geometry::Rectangle<R>& estimated_bound,
                                                        const time_type& integration_time) const;


      /*! \brief Compute a set \a bound such that the flow of \a vector_field starting in \a initial_point remains in \a bound for times within \a integration_time, given a bound \a estimated_bound. */
      virtual Geometry::Rectangle<R> refine_flow_bounds(const System::VectorFieldInterface<R>& vector_field,
                                                        const Geometry::Rectangle<R>& initial_set,
                                                        const Geometry::Rectangle<R>& estimated_bound,
                                                        const Numeric::Interval<time_type>& integration_times) const;


      /*! \brief Compute a bound for the Jacobian of the flow over the time interval [-h,h], assuming that the flow remains inside the set \a b. */
      virtual LinearAlgebra::Matrix<I> estimate_flow_jacobian_bounds(const System::VectorFieldInterface<R>& vf,
                                                                     const Geometry::Rectangle<R>& b,
                                                                     const time_type& h) const;
    };






    /*! \brief %Base class for integration schemes. 
     *  \ingroup Integrate
     */
    template<class R, class VF, class BS>
    class IntegratorBase
      : public Integrator<R> 
    {
     public:
      typedef Numeric::Interval<R> I;
      /*! \brief  */
      typedef VF vector_field_type;
      /*! \brief  */
      typedef BS basic_set_type;
      /*! \brief  */
      typedef Geometry::ListSet<BS> list_set_type;
      /*! \brief  */
      typedef Geometry::Point<I> point_type;
      /*! \brief  */
      typedef Geometry::Rectangle<R> bounding_set_type;
     protected:
      IntegratorBase(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
        : Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius) { }
     public:
      //@{
      //! \name Supporting geometric routines
      /*! \brief Subdivide the basic set into smaller pieces whose radius tends to zero with repeated subdivisions. */
      virtual list_set_type subdivide(const basic_set_type&) const;
      //@}

      //@{
      //! \name Integration of points and basic sets given bounds for the flow. */
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual point_type bounded_flow(const vector_field_type&,
                                      const point_type&,
                                      const bounding_set_type&,
                                      const time_type&) const = 0;
     
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual LinearAlgebra::Matrix<I> bounded_flow_jacobian(const vector_field_type&,
                                                             const point_type&,
                                                             const bounding_set_type&,
                                                             const time_type&) const = 0;
     
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual basic_set_type bounded_integration_step(const vector_field_type&,
                                              const basic_set_type&,
                                              const bounding_set_type&,
                                              const time_type&) const = 0;
     
      /*! \brief A reachability step for a basic set within a bounding set. */
      virtual basic_set_type bounded_reachability_step(const vector_field_type&,
                                               const basic_set_type&,
                                               const bounding_set_type&,
                                               const time_type&) const = 0;
     
      //@}

     public:

      //@{
      //! \name Single-step integration of points and basic sets, with a suggested step size. */
      /*! \brief Integrate a basic set. */
      virtual basic_set_type integration_step(const vector_field_type&,
                                  const basic_set_type&,
                                  time_type&) const;

      /*! \brief Integrate a basic set up to a given time. */
      virtual basic_set_type reachability_step(const vector_field_type&,
                                               const basic_set_type&,
                                               time_type&) const;

      //@}

      /*! \brief Integrate a basic set. */
      virtual
      Geometry::SetInterface<R>*
      integrate(const System::VectorFieldInterface<R>& vector_field, 
                const Geometry::SetInterface<R>& initial_set, 
                const time_type& time) const;

      virtual
      Geometry::SetInterface<R>*
      reach(const System::VectorFieldInterface<R>& vector_field, 
            const Geometry::SetInterface<R>& initial_set, 
            const time_type& time) const;

      virtual
      Geometry::SetInterface<R>*
      reach(const System::VectorFieldInterface<R>& vector_field, 
                const Geometry::SetInterface<R>& initial_set) const;

     public:
      
      /*! \brief Template for integrating a basic set. */
      BS 
      integrate(const VF& vector_field, 
                const BS& initial_set, 
                const time_type& time) const;

      /*! \brief Template for computing the reachable set from a basic. */
      Geometry::ListSet<BS> 
      reach(const VF& vector_field, 
            const BS& initial_set, 
            const time_type& time) const;



     public:
      
      /*! \brief Template for integrating a list set. */
      Geometry::ListSet<BS> 
      lower_integrate(const VF& vector_field, 
                      const Geometry::ListSet<BS>& initial_set, 
                      const time_type& time) const;

      
      /*! \brief Template for computing the reachable set from a list set. */
      Geometry::ListSet<BS> 
      lower_reach(const VF& vector_field, 
                  const Geometry::ListSet<BS>& initial_set, 
                  const time_type& time) const;

       /*! \brief Template for integrating a list set. */
      Geometry::ListSet<BS>
      upper_integrate(const VF& vector_field, 
                      const Geometry::ListSet<BS>& initial_set, 
                      const time_type& time) const;

      
      /*! \brief Template for computing the reachable set from a list set. */
      Geometry::ListSet<BS>
      upper_reach(const VF& vector_field, 
                  const Geometry::ListSet<BS>& initial_set, 
                  const time_type& time) const;

     public:


      /*! \brief Template for integrating a list set (computes a lower-approximation). */
      Geometry::ListSet< Geometry::Rectangle<R> >
      integrate(const System::VectorFieldInterface<R>& vector_field, 
                const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set, 
                const time_type& time) const;

      
      /*! \brief Template for computing the reachable set from a list set (computes a lower-approximation). */
       
      Geometry::ListSet< Geometry::Rectangle<R> >
      reach(const System::VectorFieldInterface<R>& vector_field, 
            const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set, 
            const time_type& time) const;


      /*! \brief Template for computing the integral set from a list set. */
      Geometry::GridMaskSet<R>
      integrate(const System::VectorFieldInterface<R>& vector_field, 
                const Geometry::GridMaskSet<R>& initial_set, 
                const Geometry::GridMaskSet<R>& bounding_set,
                const time_type& time) const;

      /*! \brief Template for computing the reachable set. */
      Geometry::GridMaskSet<R>
      reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::GridMaskSet<R>& initial_set,
                  const Geometry::GridMaskSet<R>& bounding_set,
                  const time_type& time) const;

      /*! \brief Template for computing the chain-reachable set. */
      Geometry::GridMaskSet<R>
      chainreach(const System::VectorFieldInterface<R>& vector_field, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Template for computing the viability kernel. */
      Geometry::GridMaskSet<R>
      viable(const System::VectorFieldInterface<R>& vector_field, 
             const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Template for verifying system evolution. */
      tribool
      verify(const System::VectorFieldInterface<R>& vector_field, 
             const Geometry::GridMaskSet<R>& initial_set, 
             const Geometry::GridMaskSet<R>& bounding_set) const;


     protected:
      /*! \brief The identity function. */ 
      BS identity(const BS& bs) const;
    };




    
  }
}

#endif /* ARIADNE_INTEGRATE_H */
