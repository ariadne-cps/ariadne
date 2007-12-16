/***************************************************************************
 *            model_checker.h
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
 
/*! \file model_checker.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef ARIADNE_MODEL_CHECKER_H
#define ARIADNE_MODEL_CHECKER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "system/discrete_map.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/orbiter_interface.h"

namespace Ariadne {
  namespace Evaluation {



    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class R>
    class ModelChecker {
      typedef Numeric::Interval<R> I;
     private:
      EvolutionParameters<R>* _parameters;
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Default constructor chooses appropriate parameter values for maximum basic set radius and grid size. */
      ModelChecker();
      
      /*! \brief Construct from evolution parameters. */
      ModelChecker(const EvolutionParameters<R>& parameters);
      
      /*! \brief Copy constructor. */
      ModelChecker(const ModelChecker<R>& other);
      
      /*! \brief Destructor. */
      virtual ~ModelChecker();

      /*! \brief Make a dynamically-allocated copy. */
      ModelChecker<R>* clone() const;
      
      //@}


 
      public:
      //@{ 
      //! \name Evaluation of discretized maps on concrete sets

      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< Geometry::Rectangle<R> > 
      image(const System::DiscreteMapInterface<R>& map, const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set) const;
       
      
      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridCellListSet<R> 
      image(const System::DiscreteMapInterface<R>& map, 
            const Geometry::GridCellListSet<R>& initial_set,
            const Geometry::Grid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set while remaining in \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::DiscreteMapInterface<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::GridMaskSet<R>& bounding_set) const;


      /*! \brief Compute the preimage of \a set under \a map contained in \a bound. */
      virtual
      Geometry::GridMaskSet<R> 
      preimage(const System::DiscreteMapInterface<R>& map, 
               const Geometry::GridMaskSet<R>& set,
               const Geometry::GridMaskSet<R>& bound) const;

      /*! \brief Compute the preimage of \a set under \a map contained in \a bound. */
      virtual
      Geometry::PartitionTreeSet<R> 
      preimage(const System::DiscreteMapInterface<R>& map, 
               const Geometry::PartitionTreeSet<R>& set,
               const Geometry::Rectangle<R>& bound) const;


      /*! \brief Compute an approximation to the iterated set of \a map starting in \a initial_set, iterating \a steps times. */
      virtual
      Geometry::ListSet< Geometry::Rectangle<R> > 
      iterate(const System::DiscreteMapInterface<R>& map, 
              const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
              const Numeric::Integer& steps) const;

      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set, iterating at most \a steps times. */
      virtual
      Geometry::ListSet< Geometry::Rectangle<R> > 
      reach(const System::DiscreteMapInterface<R>& map, 
            const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
            const Numeric::Integer& steps) const;

      /*! \brief Compute a lower-approximation to the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::ListSet< Geometry::Rectangle<R> > 
      lower_reach(const System::DiscreteMapInterface<R>& map, 
                  const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set) const;

      
     
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      chainreach(const System::DiscreteMapInterface<R>& map, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;
    
      /*! \brief Compute the viability kernel of \a map within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      viable(const System::DiscreteMapInterface<R>& map, 
             const Geometry::GridMaskSet<R>& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      virtual
      tribool
      verify(const System::DiscreteMapInterface<R>& map, 
             const Geometry::GridMaskSet<R>& initial_set, 
             const Geometry::GridMaskSet<R>& safe_set) const;
      //@}

      

      //@{
      //! \brief Methods for control systems

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

      //@{
      //! \name Methods for computing discretizations

      /*! \brief Discretize a system on a grid. */ 
      virtual 
      System::GridMultiMap<R> 
      discretize(const System::DiscreteMapInterface<R>& f, 
                 const Geometry::GridMaskSet<R>& dom,
                 const Geometry::Grid<R>& range_grid) const;

      //@}

     private:
      const EvolutionParameters<R>& parameters() const;

      Geometry::Rectangle<R> apply(const System::DiscreteMapInterface<R>&, const Geometry::Rectangle<R>&) const;

      Geometry::GridCellListSet<R> apply(const System::DiscreteMapInterface<R>&, const Geometry::GridCell<R>&) const;     
      
      Geometry::GridCellListSet<R> apply(const System::DiscreteMapInterface<R>&, const Geometry::GridCell<R>&, const Geometry::Grid<R>& g) const;     
      
      Geometry::DiscreteTimeOrbit< Numeric::Integer,Geometry::Rectangle<R> > 
      orbit(const System::DiscreteMapInterface<R>&, const Geometry::Rectangle<R>&, const Numeric::Integer&) const;
      
      Geometry::DiscreteTimeOrbit< Numeric::Integer,Geometry::Rectangle<R> > 
      orbit(const System::DiscreteMapInterface<R>&, const Geometry::Rectangle<R>&, const Numeric::Integer&, const R&) const;
      
      Geometry::DiscreteTimeOrbit< Numeric::Integer,Geometry::GridCellListSet<R> > 
      orbit(const System::DiscreteMapInterface<R>&, const Geometry::GridCell<R>&, const Numeric::Integer&) const;

    };



  }
}

#endif /* ARIADNE_APPLY_H */
