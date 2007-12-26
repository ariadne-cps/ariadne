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

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

namespace Ariadne {

  namespace Geometry {
    template<class R> class SetInterface;
    template<class R> class ArbitrarySet;
    template<class R> class GridMaskSet;
  }

  namespace Evaluation {

   
    template<class R> class EvolutionParameters;
    template<class R> class VectorFieldOrbiterInterface;



    /*! \brief Class for computing the evolution of a continuous-time autonomous system.
     *  \ingroup Integrate
     */
    template<class R>
    class VectorFieldEvolver {
     private:
      boost::shared_ptr< EvolutionParameters<R> >  _parameters;
      boost::shared_ptr< VectorFieldOrbiterInterface<R> >  _orbiter;
     public:
      typedef R real_type;

      //@{
      //! \name Constructors and destructors

      /*! \brief Virtual destructor. */
      virtual ~VectorFieldEvolver();

      /*! \brief Construct from evolution paramters. */
      VectorFieldEvolver(const EvolutionParameters<R>& parameters);

      /*! \brief Construct from evolution parameters and an integration method. */
      template<class BS>
      VectorFieldEvolver(const EvolutionParameters<R>& parameters, 
			 const IntegratorInterface<BS>& integrator,
			 const ApproximatorInterface<BS>& approximator);

      /*! \brief Copy constructor. */
      VectorFieldEvolver(const VectorFieldEvolver<R>& i);

      /*! \brief Cloning operator. */
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
      virtual Geometry::SetInterface<R>* evolve(const System::VectorFieldInterface<R>& vector_field,
                                                const Geometry::SetInterface<R>& initial_set,
                                                const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* bounded_evolve(const System::VectorFieldInterface<R>& vector_field,
                                                        const Geometry::SetInterface<R>& initial_set,
                                                        const Geometry::SetInterface<R>& bounding_set,
                                                        const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field, while remaining in \a bounding_set. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* bounded_reach(const System::VectorFieldInterface<R>& vector_field,
                                                       const Geometry::SetInterface<R>& initial_set,
                                                       const Geometry::SetInterface<R>& bounding_set,
                                                       const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for all times under \a vector_field, while remaining in \a bounding_box. 
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual Geometry::SetInterface<R>* chainreach(const System::VectorFieldInterface<R>& vector_field,
                                                    const Geometry::SetInterface<R>& initial_set,
                                                    const Geometry::Box<R>& bounding_box) const;

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

     private:
      //  Generate a grid
      Geometry::Grid<R> grid() const;

      //  Integrate a box.
      Geometry::GridCellListSet<R>
      evolve(const System::VectorFieldInterface<R>& vector_field, 
             const Geometry::GridCell<R>& initial_set, 
             const Numeric::Rational& time) const;

      //  Compute the reachable set from a basic.
      Geometry::GridCellListSet<R>
      reach(const System::VectorFieldInterface<R>& vector_field, 
            const Geometry::GridCell<R>& initial_set, 
            const Numeric::Rational& time) const;
    };
 

    
  }
}

#include "vector_field_evolver.inline.h"

#endif /* ARIADNE_VECTOR_FIELD_EVOLVER_H */
