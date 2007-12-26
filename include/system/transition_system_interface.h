/***************************************************************************
 *            transition_system_interface.h
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
 
/*! \file transition_system_interface.h
 *  \brief Interface for discrete transition systems with evolve and reach steps.
 */

#ifndef ARIADNE_TRANSITION_SYSTEM_INTERFACE_H
#define ARIADNE_TRANSITION_SYSTEM_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

namespace Ariadne {
  namespace System {

    template<class R>
    class TransitionSystemInterface 
    {
     public:
      virtual ~TransitionSystemInterface<R>() { }
      virtual TransitionSystemInterface<R>* clone() const = 0;
      virtual size_type dimension() const = 0;
      virtual Geometry::Box<R> lower_evolve(const Geometry::Box<R>&) const = 0;
      virtual Geometry::ListSet< Geometry::Box<R> > lower_reach(const Geometry::Box<R>&) const = 0;
      virtual Geometry::GridCellListSet<R> upper_evolve(const Geometry::GridCell<R>&) const = 0;
      virtual Geometry::GridCellListSet<R> upper_reach(const Geometry::GridCell<R>&) const = 0;
    };

  }
}

#endif /* ARIADNE_TRANSITION_SYSTEM_INTERFACE_H */
