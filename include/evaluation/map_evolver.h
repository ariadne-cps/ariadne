/***************************************************************************
 *            map_evolver.h
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
 
/*! \file map_evolver.h
 *  \brief Methods for computing the images of sets under maps.
 *
 * This class works by approximating sets based on grids. 
 * A System::DiscreteMap class is made using a MapOrbiter which is fed into a ModelChecker which runs the compuations.
 */

#ifndef ARIADNE_MAP_EVOLVER_H
#define ARIADNE_MAP_EVOLVER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"


// For templated constructor
#include "evaluation/applicator_interface.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/map_orbiter.h"


namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class R>
    class MapEvolver {
      typedef Numeric::Integer T;
      typedef Numeric::Interval<R> I;
     private:
      EvolutionParameters<R>* _parameters;
      MapOrbiterInterface<R>* _orbiter;
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Default constructor chooses appropriate parameter values for maximum basic set radius and grid size. */
      MapEvolver();
      
      /*! \brief Construct from evolution parameters. */
      MapEvolver(const EvolutionParameters<R>& parameters);
      
      /*! \brief Construct from evolution parameters and a method for iterating basic sets. */
      template<class BS>
      MapEvolver(const EvolutionParameters<R>& parameters, 
                 const ApplicatorInterface<BS>& applicator, 
                 const ApproximatorInterface<BS>& approximator);
      
      /*! \brief Copy constructor. */
      MapEvolver(const MapEvolver<R>& other);
      
      /*! \brief Destructor. */
      virtual ~MapEvolver();

      /*! \brief Make a dynamically-allocated copy. */
      MapEvolver<R>* clone() const;
      
      //@}


      //@{ 
      //! \name Methods to set and get the parameters controlling the accuracy.

      /*! \brief The parameters controlling the accuracy. */
      virtual const EvolutionParameters<R>& parameters() const;

      /*! \brief A reference to the parameters controlling the accuracy. */
      virtual EvolutionParameters<R>& parameters();

      //@}


      //@{
      //! \name Evaluation of maps on abstract sets

    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      lower_evolve(const System::MapInterface<R>& map, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      lower_reach(const System::MapInterface<R>& map, 
                  const Geometry::SetInterface<R>& initial_set,
                  const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      upper_evolve(const System::MapInterface<R>& map, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      upper_reach(const System::MapInterface<R>& map, 
                  const Geometry::SetInterface<R>& initial_set,
                  const Numeric::Integer& steps) const;
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      chainreach(const System::MapInterface<R>& map, 
                 const Geometry::SetInterface<R>& initial_set) const;
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::SetInterface<R>*
      chainreach(const System::MapInterface<R>& map, 
                 const Geometry::SetInterface<R>& initial_set, 
                 const Geometry::Box<R>& bounding_set) const;
    
      /*! \brief Compute an outer-approximation to the viability kernel of \a map within \a bounding_set. */
      virtual
      Geometry::SetInterface<R>* 
      viable(const System::MapInterface<R>& map, 
             const Geometry::SetInterface<R>& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      virtual
      tribool
      verify(const System::MapInterface<R>& map, 
             const Geometry::SetInterface<R>& initial_set, 
             const Geometry::SetInterface<R>& safe_set) const;
      //@}


      //@{
      //! \name Methods for computing discretizations of maps on grids

      /*! \brief Discretize  \a map for cells in \a domain_set with image discretized in \a range_grid. */
      virtual 
      System::GridMultiMap<R> 
      discretize(const System::MapInterface<R>& map, 
                 const Geometry::GridMaskSet<R>& domain_set,
                 const Geometry::Grid<R>& range_grid) const;

      //@}

 
      //@{ 
      //! \name Methods for computing orbits.

      /*! \brief Compute the orbit of a rectangle under steps of continuous function. */
      virtual 
      Geometry::OrbitInterface< Numeric::Integer>*
      orbit(const System::MapInterface<R>& f, const Geometry::Box<R>& r, const Numeric::Integer& n) const;

      //@}

     private:
      // Helper functions
      Evaluation::EvolutionParameters<R>* default_parameters();
      Evaluation::MapOrbiterInterface<R>* default_orbiter();
      Evaluation::ModelChecker<R> model_checker() const;
      System::TransitionSystem<R> discrete_map(const System::MapInterface<R>& f) const;
      Geometry::GridMaskSet<R> outer_approximation(const Geometry::SetInterface<R>& f) const;
      Geometry::ListSet< Geometry::Box<R> > lower_approximation(const Geometry::SetInterface<R>& f) const;
      Geometry::GridMaskSet<R> inner_approximation(const Geometry::SetInterface<R>& f) const;
    };



  }
}


// Inline functions

namespace Ariadne {

template<class R> template<class BS> inline
Evaluation::MapEvolver<R>::MapEvolver(const EvolutionParameters<R>& parameters,
                                      const ApplicatorInterface<BS>& applicator,
                                      const ApproximatorInterface<BS>& approximator)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _orbiter(new MapOrbiter<BS>(parameters,applicator,approximator))
{
}

} 


#endif /* ARIADNE_MAP_EVOLVER_H */
