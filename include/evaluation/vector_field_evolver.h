/***************************************************************************
 *            vector_field_evolver.h
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
 
/*! \file vector_field_evolver.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_H
#define ARIADNE_VECTOR_FIELD_EVOLVER_H

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
    template<class BS> class IntegratorInterface;
    template<class R> class BounderInterface;



    /*! \brief Class for computing the evolution of a continuous-time autonomous system.
     *  \ingroup Integrate
     */
    template<class R>
    class VectorFieldEvolver {
     private:
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I,I> BS;

      EvolutionParameters<R>* _parameters;
      BounderInterface<R>* _bounder;
      IntegratorInterface<BS>* _integrator;
     public:
      typedef R real_type;

      typedef Numeric::Interval<R> interval_type;
     
      typedef Geometry::Point<I> point_type;

      typedef Geometry::Rectangle<R> bounding_set_type;

      typedef BS basic_set_type;

      typedef Geometry::ListSet<BS> list_set_type;

      typedef System::VectorFieldInterface<R> vector_field_type;

      //@{
      //! \name Constructors and destructors

      /*! \brief Virtual destructor. */
      virtual ~VectorFieldEvolver();

      /*! \brief Construct from evolution paramters. */
      VectorFieldEvolver(const EvolutionParameters<R>& parameters);

      /*! \brief Construct from evolution parameters and an integration method. */
      template<class BST>
      VectorFieldEvolver(const EvolutionParameters<R>& parameters, const IntegratorInterface<BST>& plugin);

      /*! \brief Copy constructor. */
      VectorFieldEvolver(const VectorFieldEvolver<R>& i);

      /*! \brief Make a dynamically-allocated copy. */
      virtual VectorFieldEvolver<R>* clone() const;
      //@}

      //@{
      //! \name Parameters controlling the accuracy

      /*! \brief The parameters controlling the accuracy. */
      const EvolutionParameters<R>& parameters() const;
      /*! \brief A reference to the parameters controlling the accuracy. */
      EvolutionParameters<R>& parameters();
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
      virtual Geometry::SetInterface<R>* lower_reach(const System::VectorFieldInterface<R>& vector_field,
                                                     const Geometry::SetInterface<R>& initial_set) const;

      /*! \brief Integrate \a intial_set for all times under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* lower_reach(const System::VectorFieldInterface<R>& vector_field,
                                                     const Geometry::SetInterface<R>& initial_set,
                                                     const Geometry::SetInterface<R>& bounding_set) const;

      /*! \brief Integrate \a intial_set for all times under \a vector_field. 
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual Geometry::SetInterface<R>* chainreach(const System::VectorFieldInterface<R>& vector_field,
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
      //! \name Single-step integration of points and basic sets, with a suggested step size. */
      /*! \brief Integrate a basic set. */
      virtual bounding_set_type flow_bounds(const vector_field_type&,
                                            const basic_set_type&,
                                            time_type&) const;

      /*! \brief Integrate a basic set, assuming the flow remains in the bounding set. */
      virtual basic_set_type integration_step(const vector_field_type&,
                                              const basic_set_type&,
                                              const time_type&,
                                              const bounding_set_type&) const;

      /*! \brief Integrate a basic set up to a given time assuming the flow remains in the bounding set. */
      virtual basic_set_type reachability_step(const vector_field_type&,
                                               const basic_set_type&,
                                               const time_type&,
                                               const bounding_set_type&) const;


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
      //@{ 
      //! \name Integration of concrete sets. 

      /*! \brief Integrating a list set (computes a lower-approximation). */
      Geometry::ListSet< Geometry::Rectangle<R> >
      lower_integrate(const System::VectorFieldInterface<R>& vector_field, 
                      const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set, 
                      const time_type& time) const;
      
      
      /*! \brief Computing the timed reachable set from a list set (computes a lower-approximation). */
      Geometry::ListSet< Geometry::Rectangle<R> >
      lower_reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set, 
                  const time_type& time) const;

      /*! \brief Computing the timed reachable set from a list set (computes a lower-approximation). */
      Geometry::ListSet< Geometry::Rectangle<R> >
      lower_reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set, 
                  const Geometry::SetInterface<R>& bounding_set, 
                  const time_type& time) const;

      /*! \brief Integrating a list set. (%Deprecated) */
      Geometry::ListSet<BS> 
      lower_integrate(const System::VectorFieldInterface<R>& vector_field, 
                      const Geometry::ListSet<BS>& initial_set, 
                      const time_type& time) const;

      
      /*! \brief Compute the reachable set from a list set. (%Deprecated) */
      Geometry::ListSet<BS> 
      lower_reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::ListSet<BS>& initial_set, 
                  const time_type& time) const;

      /*! \brief Compute the reachable set from a list set. (%Deprecated) */
      Geometry::ListSet<BS> 
      lower_reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::ListSet<BS>& initial_set, 
                  const Geometry::SetInterface<R>& bounding_set, 
                  const time_type& time) const;

       /*! \brief Integrate a list set. (%Deprecated) */
      Geometry::ListSet<BS>
      upper_integrate(const System::VectorFieldInterface<R>& vector_field, 
                      const Geometry::ListSet<BS>& initial_set, 
                      const time_type& time) const;

      
      /*! \brief Template for computing the reachable set from a list set. (%Deprecated) */
      Geometry::ListSet<BS>
      upper_reach(const System::VectorFieldInterface<R>& vector_field, 
                  const Geometry::ListSet<BS>& initial_set, 
                  const time_type& time) const;


      /*! \brief Computing the integral set from a grid mask set while remaining in a bounding set. */
      Geometry::GridMaskSet<R>
      bounded_integrate(const System::VectorFieldInterface<R>& vector_field, 
                        const Geometry::GridMaskSet<R>& initial_set, 
                        const Geometry::GridMaskSet<R>& bounding_set,
                        const time_type& time) const;

      /*! \brief Template for computing the reachable set. */
      Geometry::GridMaskSet<R>
      bounded_reach(const System::VectorFieldInterface<R>& vector_field, 
                    const Geometry::GridMaskSet<R>& initial_set,
                    const Geometry::GridMaskSet<R>& bounding_set,
                    const time_type& time) const;

      /*! \brief Compute the chain-reachable set. */
      Geometry::GridMaskSet<R>
      chainreach(const System::VectorFieldInterface<R>& vector_field, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Compute the viability kernel. */
      Geometry::GridMaskSet<R>
      viable(const System::VectorFieldInterface<R>& vector_field, 
             const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Verify system evolution. */
      tribool
      verify(const System::VectorFieldInterface<R>& vector_field, 
             const Geometry::GridMaskSet<R>& initial_set, 
             const Geometry::GridMaskSet<R>& bounding_set) const;

      //@}

     private:

      /*! \brief Integrate a basic set. */
      BS
      integrate(const System::VectorFieldInterface<R>& vector_field, 
                const BS& initial_set, 
                const time_type& time) const;

      /*! \brief Compute the reachable set from a basic. */
      Geometry::ListSet<BS> 
      reach(const System::VectorFieldInterface<R>& vector_field, 
            const BS& initial_set, 
            const time_type& time) const;


    private:
      //@{
      //! \name Supporting geometric routines
      /*! \brief Subdivide the basic set into smaller pieces whose radius tends to zero with repeated subdivisions. */
      virtual list_set_type subdivide(const basic_set_type&) const;
      //@}

    };
 

    
  }
}

#include "vector_field_evolver.inline.h"

#endif /* ARIADNE_INTEGRATE_H */
