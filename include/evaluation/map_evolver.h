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
 */

#ifndef ARIADNE_APPLICATOR_H
#define ARIADNE_APPLICATOR_H

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

#include "../evaluation/evolution_parameters.h"
#include "../evaluation/applicator_plugin_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a set under a map. 
     *  \ingroup Applicators
     */
    template<class R>
    class MapEvolver {
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I,R> BS;
     private:
      EvolutionParameters<R>* _parameters;
      ApplicatorPluginInterface<R>* _plugin;
     public:
      /*! \brief Default constructor chooses appropriate parameter values for maximum basic set radius and grid size. */
      MapEvolver();
      
      /*! \brief Construct from evolution parameters. */
      MapEvolver(const EvolutionParameters<R>& parameters);
      
      /*! \brief Copy constructor. */
      MapEvolver(const MapEvolver<R>& other);
      
      /*! \brief Compute the image of a basic set under a continuous function. */
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
      //! \name Supporting geometric routines.
      /*! \brief Subdivide a basic set into smaller pieces. */
      virtual Geometry::ListSet< BS > subdivide(const BS& bs) const;
      //@}


      //@{ 
      //! \name Methods for applying a system to a basic set.

      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual 
      BS 
      evaluate(const System::MapInterface<R>& f, const BS& bs) const;

      //@}

     public:
      //@{ 
      //! \name Evaluation of maps on concrete sets

      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< Geometry::Rectangle<R> > 
      image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Rectangle<R> >& ds) const;
       
      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< BS > 
      image(const System::MapInterface<R>& f, const Geometry::ListSet< BS >& ds) const;
       
      
      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridCellListSet<R> 
      image(const System::MapInterface<R>& map, 
            const Geometry::GridCellListSet<R>& initial_set,
            const Geometry::Grid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      /*
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::MapInterface<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::FiniteGrid<R>& grid) const;
      */

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


      /*! \brief Compute the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::ListSet< Geometry::Rectangle<R> > 
      reach(const System::MapInterface<R>& map, 
            const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set) const;

           
      /*! \brief Compute the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::ListSet<BS> 
      reach(const System::MapInterface<R>& map, 
            const Geometry::ListSet<BS>& initial_set) const;

           
      /*! \brief Compute the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::GridMaskSet<R> 
      reach(const System::MapInterface<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set) const;

            
      /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
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
    
      /*! \brief Compute the reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::SetInterface<R>*
      reach(const System::MapInterface<R>& map, 
            const Geometry::SetInterface<R>& initial_set) const;
    
      /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::SetInterface<R>*
      chainreach(const System::MapInterface<R>& map, 
                 const Geometry::SetInterface<R>& initial_set, 
                 const Geometry::SetInterface<R>& bounding_set) const;
    
      /*! \brief Compute the viability kernel of \a map within \a bounding_set. */
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
      //! \brief Methods for computing discretizations

      /*! \brief Discretize a system on a grid. */ 
      virtual 
      System::GridMultiMap<R> 
      discretize(const System::MapInterface<R>& f, 
                 const Geometry::GridMaskSet<R>& dom,
                 const Geometry::Grid<R>& range_grid) const;

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




       
    };



  }
}

#endif /* ARIADNE_APPLY_H */
