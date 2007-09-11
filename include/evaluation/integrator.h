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

   
    template<class R> class EvolutionParameters;
    template<class BS> class IntegratorPluginInterface;
    template<class BS> class BounderPluginInterface;



    /*! \brief %Base class for integration schemes. 
     *  \ingroup Integrate
     */
    template<class BS>
    class Integrator {
     private:
      typedef typename BS::real_type R;
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef Numeric::Interval<R> I;

      EvolutionParameters<R>* _parameters;
      BounderPluginInterface<R>* _bounder_plugin;
      IntegratorPluginInterface<BS>* _integrator_plugin;
     public:
      /*! The type used to denote real numbers. */
      typedef R real_type;

      /*! The type used to denote an interval of real numbers. */
      typedef Numeric::Interval<R> interval_type;
     
      /*! The type used to denote real numbers. */
      typedef Geometry::Point<I> point_type;

      /*! The type used to denote real numbers. */
      typedef Geometry::Rectangle<R> bounding_set_type;

      /*! The type used to denote real numbers. */
      typedef BS basic_set_type;

      /*! The type used to denote real numbers. */
      typedef Geometry::ListSet<BS> list_set_type;

      /*! The type used to denote real numbers. */
      typedef System::VectorFieldInterface<R> vector_field_type;

      /*! \brief Virtual destructor. */
      virtual ~Integrator();

      /*! \brief Constructor. */
      Integrator(const EvolutionParameters<R>& parameters, const IntegratorPluginInterface<BS>& plugin);

      /*! \brief Copy constructor. */
      Integrator(const Integrator<BS>& i);

      /*! \brief Make a dynamically-allocated copy. */
      virtual Integrator<BS>* clone() const;

      //@{
      //! \name Parameters controlling the accuracy

      /*! \brief The parameters controlling the accuracy. */
      const EvolutionParameters<R>& parameters() const;
      /*! \brief A reference to the parameters controlling the accuracy. */
      EvolutionParameters<R>& parameters();

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
      virtual R grid_length() const;
      //@}

      //@{
      //! \name Integration routines


     /*! \brief Integrate \a intial_set for time \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* integrate(const System::VectorFieldInterface<R>& vector_field,
                                                   const Geometry::SetInterface<R>& initial_set,
                                                   const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* integrate(const System::VectorFieldInterface<R>& vector_field,
                                                   const Geometry::SetInterface<R>& initial_set,
                                                   const Geometry::SetInterface<R>& bounding_set,
                                                   const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Geometry::SetInterface<R>& bounding_set,
                                               const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for all times under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set) const;

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


     public:
      //@{
      //! \name Supporting geometric routines
      /*! \brief Subdivide the basic set into smaller pieces whose radius tends to zero with repeated subdivisions. */
      virtual list_set_type subdivide(const basic_set_type&) const;
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

     public:
      
      /*! \brief Template for integrating a basic set. */
      BS
      integrate(const System::VectorFieldInterface<R>& vector_field, 
                const BS& initial_set, 
                const time_type& time) const;

      /*! \brief Template for computing the reachable set from a basic. */
      Geometry::ListSet<BS> 
      reach(const System::VectorFieldInterface<R>& vector_field, 
            const BS& initial_set, 
            const time_type& time) const;



     public:
      
      /*! \brief Template for integrating a list set. */
      Geometry::ListSet<BS> 
      lower_integrate(const System::VectorFieldInterface<R>& vector_field, 
                      const Geometry::ListSet<BS>& initial_set, 
                      const time_type& time) const;

      
      /*! \brief Template for computing the reachable set from a list set. */
      Geometry::ListSet<BS> 
      lower_reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::ListSet<BS>& initial_set, 
                  const time_type& time) const;

       /*! \brief Template for integrating a list set. */
      Geometry::ListSet<BS>
      upper_integrate(const System::VectorFieldInterface<R>& vector_field, 
                      const Geometry::ListSet<BS>& initial_set, 
                      const time_type& time) const;

      
      /*! \brief Template for computing the reachable set from a list set. */
      Geometry::ListSet<BS>
      upper_reach(const System::VectorFieldInterface<R>& vector_field, 
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
