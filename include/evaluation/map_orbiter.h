/***************************************************************************
 *            map_orbiter.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file map_orbiter.h
 *  \brief Methods for computing orbits of boxes.
 */

#ifndef ARIADNE_MAP_ORBITER_H
#define ARIADNE_MAP_ORBITER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/orbiter_interface.h"

namespace Ariadne {
  namespace Evaluation {

      
    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class BS>
    class MapOrbiter 
      : public MapOrbiterInterface<typename BS::real_type>
    {
      typedef typename BS::real_type R;
     public:
      /*! \brief Default construnctor. */
      //MapOrbiter();

      /*! \brief Construct from the evolution parameters and an applicator. */
      MapOrbiter(const EvolutionParameters<R>&, const ApplicatorInterface<BS>&);

      /*! \brief Copy constructor. */
      MapOrbiter(const MapOrbiter<BS>&);

      /*! \brief Make a dynamically-allocated copy. */
      MapOrbiter<BS>* clone() const;


      /*! \brief Compute the image of a rectangle under a continuous function. */
      virtual 
      Geometry::Rectangle<R> 
      apply(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r) const;

      /*! \brief Compute the image of a grid cell under a continuous self-map. */
      virtual 
      Geometry::GridCellListSet<R> 
      apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& r) const;

      /*! \brief Compute the image of a grid cell under a continuous function, approximating on grid \a g which may lie in a different space. */
      virtual 
      Geometry::GridCellListSet<R> 
      apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& r, const Geometry::Grid<R>& g) const;



      /*! \brief Compute the orbit of a rectangle under \a n steps of continuous function. */
      virtual 
      Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Rectangle<R> >
      orbit(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r, const Numeric::Integer& n) const;

      /*! \brief Compute the orbit of a rectangle under at most \a n steps of continuous function, until the size reaches \a s. */
      virtual 
      Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Rectangle<R> >
      orbit(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r, const Numeric::Integer& n, const R& s) const;

      /*! \brief Compute the orbit of a grid cell under steps of continuous function. */
      virtual 
      Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::GridCellListSet<R> >
      orbit(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const;

     private:
      // Convenience wrapper for applicator services
      BS apply(const System::MapInterface<R>& f, const BS& bs) const;
     private:
      boost::shared_ptr< ApplicatorInterface<BS> > _applicator;
    };

  }

}

#endif /* ARIADNE_MAP_ORBITER_H */
