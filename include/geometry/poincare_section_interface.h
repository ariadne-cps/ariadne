/***************************************************************************
 *            poincare_section_interface.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file poincare_section_interface.h
 *  \brief A hypersurface suitable for defining a Poincare' section.
  */

#ifndef ARIADNE_POINCARE_SECTION_INTERFACE_H
#define ARIADNE_POINCARE_SECTION_INTERFACE_H

#include <iosfwd>

#include "base/types.h"

#include "function/declarations.h"

namespace Ariadne {
  namespace Geometry {
    
    // Forward declarations for friends
    template<class R> class ConstraintInterface;
    
    //! \ingroup SetInterface
    /*! \brief A hypersurface suitable for use as a Poincare' section.
     */
    template<class R>
    class PoincareSectionInterface
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Destructor. */
      virtual ~PoincareSectionInterface() { };
      /*! \brief Return a new dynamically-allocated copy of the section. */
      virtual PoincareSectionInterface<R>* clone() const = 0;
      /*! \brief The dimension of the section. */
      virtual dimension_type dimension() const = 0;
      /*! \brief The smoothness of the section. */
      virtual smoothness_type smoothness() const = 0;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;

      /*! \brief The map used to define the inclusion of the section in state space. */
      virtual const Function::FunctionInterface<R>& inclusion_map() const = 0;
      /*! \brief The map used to define the (local) projection of state space to the section. */
      virtual const Function::FunctionInterface<R>& projection_map() const = 0;
      /*! \brief A function whose zero set is the section. */
      virtual const Function::FunctionInterface<R>& crossing_condition() const = 0;
    };
    
    template<class R> inline std::ostream& operator<<(std::ostream& os, const PoincareSectionInterface<R>& ps) {
      return ps.write(os);
    }
    
  }
}

#endif /* ARIADNE_POINCARE_SECTION_INTERFACE_H */
