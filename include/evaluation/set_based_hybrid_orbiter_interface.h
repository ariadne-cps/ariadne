/***************************************************************************
 *            set_based_hybrid_orbiter_interface.h
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
 
/*! \file set_based_hybrid_orbiter_interface.h
 *  \brief Methods for computing orbits of boxes under a continuous-time system.
 */

#ifndef ARIADNE_SET_BASED_HYBRID_ORBITER_INTERFACE_H
#define ARIADNE_SET_BASED_HYBRID_ORBITER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"


namespace Ariadne {
  namespace Evaluation {

 

    /*! \brief A class for computing the evolution of a continuous-time autonomous system.
     *  \ingroup Integrators
     */
    template<class R>
    class SetBasedHybridOrbiterInterface 
    {
     public:
      /*! \brief The type used to denote time. */
      typedef Numeric::Rational time_type;
      /*! \brief The type used to represent space. */
      typedef R real_type;
     public:
      /*! \brief Destructor. */
      virtual ~SetBasedHybridOrbiterInterface();

      /*! \brief Make a dynamically-allocated copy. */
      virtual SetBasedHybridOrbiterInterface<R>* clone() const = 0;

      /*! \brief Compute an over-approximation to the time \a t evolution of the grid cell \a gc. */
      virtual 
      Geometry::HybridGridCellListSet<R>
      upper_evolve(const System::HybridAutomaton<R>& ha, 
                   const Geometry::HybridGridCell<R>& gc, 
                   const Numeric::Rational& t) const = 0;

      /*! \brief Compute an over-approximation to the time \a t reachable set of the grid cell \a gc. */
      virtual 
      Geometry::HybridGridCellListSet<R>
      upper_reach(const System::HybridAutomaton<R>& vf, 
                  const Geometry::HybridGridCell<R>& bx, 
                  const Numeric::Rational& t) const = 0;

      /*! \brief Compute a lower-approximation to the time \a t evolution of the box \a bx. */
      virtual 
      Geometry::HybridBox<R>
      lower_evolve(const System::HybridAutomaton<R>& vf, 
                   const Geometry::HybridBox<R>& bx, 
                   const Numeric::Rational& t) const = 0;

      /*! \brief Compute a lower-approximation to the time \a t reachable set of the box \a bx. */
      virtual 
      Geometry::HybridBoxListSet<R>
      lower_reach(const System::HybridAutomaton<R>& vf, 
                  const Geometry::HybridBox<R>& bx, 
                  const Numeric::Rational& t) const = 0;

         
    };

    template<class R> SetBasedHybridOrbiterInterface<R>::~SetBasedHybridOrbiterInterface() { }

  }

}

#endif /* ARIADNE_SET_BASED_HYBRID_ORBITER_INTERFACE_H */
