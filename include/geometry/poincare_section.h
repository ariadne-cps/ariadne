/***************************************************************************
 *            poincare_section.h
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
 
/*! \file poincare_section.h
 *  \brief A hypersurface suitable for defining a Poincare' section.
  */

#ifndef ARIADNE_POINCARE_SECTION_H
#define ARIADNE_POINCARE_SECTION_H

#include "poincare_section_interface.h"

namespace Ariadne {
  namespace Geometry {
    

    /*! \brief A hypersurface suitable for use as a Poincare' section.
     */
    template<class R>
    class PoincareSection
      : public PoincareSectionInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief Destructor. */
      virtual ~PoincareSection();
      /*! \brief Constructor. */
      PoincareSection(const System::MapInterface<R>& im, const System::MapInterface<R>& pm, const ConstraintInterface<R>& c);
      /*! \brief Copy constructor. */
      PoincareSection(const PoincareSection<R>& s);
      /*! \brief Return a new dynamically-allocated copy of the constraint. */
      virtual PoincareSection<R>* clone() const;
      /*! \brief The dimension of the section. */
      virtual dimension_type dimension() const;
      /*! \brief The smoothness of the constraint function. */
      virtual smoothness_type smoothness() const;
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;

      virtual const System::MapInterface<R>& inclusion_map() const;
      virtual const System::MapInterface<R>& projection_map() const;
      virtual const Geometry::ConstraintInterface<R>& crossing_condition() const;
     private:
      System::MapInterface<R>* _inclusion_map_ptr;
      System::MapInterface<R>* _projection_map_ptr;
      Geometry::ConstraintInterface<R>* _constraint_ptr;
    };
    

  }
}

#endif /* ARIADNE_POINCARE_SECTION_H */
