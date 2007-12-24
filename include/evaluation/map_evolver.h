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
      MapEvolver(const EvolutionParameters<R>& parameters, const ApplicatorInterface<BS>& applicator);
      
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

      /*! \brief Compute the image of \a set under \a map. */
      virtual
      Geometry::SetInterface<R>* 
      image(const System::MapInterface<R>& map, 
            const Geometry::SetInterface<R>& set) const;
    
      /*! \brief Compute the preimage of \a set under \a map. */
      virtual
      Geometry::SetInterface<R>*
      preimage(const System::MapInterface<R>& map, 
               const Geometry::SetInterface<R>& set, 
               const Geometry::SetInterface<R>& bound) const;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      iterate(const System::MapInterface<R>& map, 
              const Geometry::SetInterface<R>& initial_set,
              const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      reach(const System::MapInterface<R>& map, 
            const Geometry::SetInterface<R>& initial_set,
            const Numeric::Integer& steps) const;
    
      /*! \brief Compute a lower-approximation to the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      lower_reach(const System::MapInterface<R>& map, 
                  const Geometry::SetInterface<R>& initial_set) const;
    
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
      //! \name Methods for applying a system to a basic set and computing orbits.

     public:
      /*! \brief Compute the image of a rectangle under a continuous function. */
      virtual 
      Geometry::Box<R> 
      apply(const System::MapInterface<R>& f, const Geometry::Box<R>& r) const;

      /*! \brief Compute the image of a grid cell under a continuous self-map. */
      virtual 
      Geometry::GridCellListSet<R> 
      apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& r) const;

      /*! \brief Compute the image of a grid cell under a continuous function, approximating on grid \a g which may lie in a different space. */
      virtual 
      Geometry::GridCellListSet<R> 
      apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& r, const Geometry::Grid<R>& g) const;
      //@}


      //@{ 
      //! \name Methods for computing orbits.

      /*! \brief Compute the orbit of a rectangle under steps of continuous function. */
      virtual 
      Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Box<R> >
      orbit(const System::MapInterface<R>& f, const Geometry::Box<R>& r, const Numeric::Integer& n) const;

      /*! \brief Compute the orbit of a rectangle under \a n steps of continuous function, until the size reaches \a s. */
      virtual 
      Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Box<R> >
      orbit(const System::MapInterface<R>& f, const Geometry::Box<R>& r, const Numeric::Integer& n, const R& s) const;

      /*! \brief Compute the orbit of a grid cell under steps of continuous function. */
      virtual 
      Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::GridCellListSet<R> >
      orbit(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const;
      //@}

     public:
      //@{ 
      //! \name Evaluation of maps on concrete sets

      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< Geometry::Box<R> > 
      image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Box<R> >& ds) const;
       
      
      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridCellListSet<R> 
      image(const System::MapInterface<R>& map, 
            const Geometry::GridCellListSet<R>& initial_set,
            const Geometry::Grid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set while remaining in \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::MapInterface<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::GridMaskSet<R>& bounding_set) const;


      /*! \brief Compute the preimage of \a set under \a map contained in \a bound. */
      virtual
      Geometry::GridMaskSet<R> 
      preimage(const System::MapInterface<R>& map, 
               const Geometry::GridMaskSet<R>& set,
               const Geometry::GridMaskSet<R>& bound) const;

      /*! \brief Compute the preimage of \a set under \a map contained in \a bound. */
      virtual
      Geometry::PartitionTreeSet<R> 
      preimage(const System::MapInterface<R>& map, 
               const Geometry::PartitionTreeSet<R>& set,
               const Geometry::Box<R>& bound) const;


      /*! \brief Compute an approximation to the iterated set of \a map starting in \a initial_set, iterating \a steps times. */
      virtual
      Geometry::ListSet< Geometry::Box<R> > 
      iterate(const System::MapInterface<R>& map, 
              const Geometry::ListSet< Geometry::Box<R> >& initial_set,
              const Numeric::Integer& steps) const;

      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set, iterating at most \a steps times. */
      virtual
      Geometry::ListSet< Geometry::Box<R> > 
      reach(const System::MapInterface<R>& map, 
            const Geometry::ListSet< Geometry::Box<R> >& initial_set,
            const Numeric::Integer& steps) const;

      /*! \brief Compute a lower-approximation to the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::ListSet< Geometry::Box<R> > 
      lower_reach(const System::MapInterface<R>& map, 
                  const Geometry::ListSet< Geometry::Box<R> >& initial_set) const;

      
     
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      chainreach(const System::MapInterface<R>& map, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;
    
      /*! \brief Compute the viability kernel of \a map within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      viable(const System::MapInterface<R>& map, 
             const Geometry::GridMaskSet<R>& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      virtual
      tribool
      verify(const System::MapInterface<R>& map, 
             const Geometry::GridMaskSet<R>& initial_set, 
             const Geometry::GridMaskSet<R>& safe_set) const;
      //@}

      

      //@{
      //! \name Methods for control systems

      /*! \brief Compute a controller for a control-to-target problem. */ 
      virtual 
      System::GridMultiMap<R> 
      control_synthesis(const System::DiscreteTimeSystem<R>& f, 
                        const Geometry::SetInterface<R>& initial_set,
                        const Geometry::SetInterface<R>& target_set,
                        const Geometry::GridMaskSet<R>& state_bounding_set,
                        const Geometry::GridMaskSet<R>& control_set,
                        const Geometry::GridMaskSet<R>& noise_set) const;
     
      //@}


     private:
      // Helper functions
      Evaluation::EvolutionParameters<R>* default_parameters();
      Evaluation::MapOrbiterInterface<R>* default_orbiter();
      Evaluation::ModelChecker<R> model_checker() const;
      System::DiscreteMap<R> discrete_map(const System::MapInterface<R>& f) const;
      
    };



  }
}


// Inline functions

namespace Ariadne {

template<class R> template<class BS> inline
Evaluation::MapEvolver<R>::MapEvolver(const EvolutionParameters<R>& parameters,
                                      const ApplicatorInterface<BS>& applicator)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _orbiter(new MapOrbiter<BS>(parameters,applicator))
{
}

}


#endif /* ARIADNE_APPLY_H */
