/***************************************************************************
 *            map_orbiter_interface.h
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
 
/*! \file map_orbiter_interface.h
 *  \brief Abstract class for computing for computing orbits of boxes.
 */

#ifndef ARIADNE_MAP_ORBITER_INTERFACE_H
#define ARIADNE_MAP_ORBITER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"


namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class R>
    class MapOrbiterInterface 
    {
     public:
      /*! \brief Destructor. */
     virtual ~MapOrbiterInterface();

      /*! \brief Make a dynamically-allocated copy. */
      virtual MapOrbiterInterface<R>* clone() const = 0;

      /*! \brief Compute the evolved set under a continuous function. */
      virtual 
      Geometry::GridCellListSet<R> 
      upper_evolve(const System::MapInterface<R>& f, 
                   const Geometry::Box<R>& s, 
                   const Numeric::Integer& n) const = 0;

      /*! \brief Compute the reach set under a continuous function. */
      virtual 
      Geometry::GridCellListSet<R> 
      upper_reach(const System::MapInterface<R>& f, const 
                  Geometry::Box<R>& s, 
                  const Numeric::Integer& n) const = 0;

      /*! \brief Compute the evolved set under a continuous function. */
      virtual 
      Geometry::ListSet< Geometry::Box<R> > 
      lower_evolve(const System::MapInterface<R>& f, 
                   const Geometry::Box<R>& s, 
                   const Numeric::Integer& n) const = 0;

      /*! \brief Compute the reach set under a continuous function. */
      virtual  
      Geometry::ListSet< Geometry::Box<R> > 
      lower_reach(const System::MapInterface<R>& f, 
                  const Geometry::Box<R>& s, 
                  const Numeric::Integer& n) const = 0;

      /*! \brief Compute an orbit under a continuous function. */
      virtual
      Geometry::OrbitInterface<Numeric::Integer>*
      orbit(const System::MapInterface<R>& f, 
            const Geometry::Box<R>& s, 
            const Numeric::Integer& n) const = 0;
  
    };

    template<class R> MapOrbiterInterface<R>::~MapOrbiterInterface() { }


  }

}

#endif /* ARIADNE_ORBITER_INTERFACE_H */
