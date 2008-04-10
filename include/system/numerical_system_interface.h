/***************************************************************************
 *            numerical_system_interface.h
 *
 *  Copyright  2007  Pieter Collins
 *
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
 
/*! \file numerical_system_interface.h
 *  \brief Interface for discrete numerical systems with evolve and reach steps.
 */

#ifndef ARIADNE_NUMERICAL_SYSTEM_INTERFACE_H
#define ARIADNE_NUMERICAL_SYSTEM_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

namespace Ariadne {

  namespace System {

    /*!\ingroup System
     * \brief Interface for numerical systems defined by reach and evolve operators. 
     */
    template<class T, class ES>
    class NumericalSystemInterface 
    {
      typedef Geometry::ListSet<ES> ESL;
     public:
      /*! \brief The type used for a single enclosure set. */
      typedef ES enclosure_set_type;
      /*! \brief The type used for a single enclosure set. */
      typedef ESL enclosure_set_list_type;
      /*! \brief The type used to describe the state space. */
      typedef typename ES::space_type state_space_type;
      /*! \brief The type used to represent real numbers. */
      typedef typename ES::real_type real_type;
      /*! \brief The type used to represent time. */
      typedef T time_type;

      /*! \brief Virtual destructor. */
      virtual ~NumericalSystemInterface() { }
      /*! \brief Cloning operation. */
      virtual NumericalSystemInterface<T,ES>* clone() const = 0;

      /*! \brief The underlying state space of the system. */
      virtual Geometry::EuclideanSpace state_space() const = 0;

      /*! \brief Compute a lower approximation to the evolution for time \a t starting in set \a s. */
      virtual ESL evolve(const ES& initial, const T& time, Evaluation::Semantics semantics) const = 0;
      /*! \brief Compute a lower approximation to the evolution for time \a t starting in set \a s. */
      virtual ESL reach(const ES& initial, const T& time, Evaluation::Semantics semantics) const = 0;
      /*! \brief Compute a lower approximation to the evolution for time \a t starting in set \a s. */
      virtual std::pair<ESL,ESL> reach_evolve(const ES& initial, const T& t, Evaluation::Semantics semantics) const = 0;
    };

  }
}

#endif /* ARIADNE_NUMERICAL_SYSTEM_INTERFACE_H */
