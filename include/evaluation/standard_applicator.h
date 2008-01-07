/***************************************************************************
 *            standard_applicator.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file standard_applicator.h
 *  \brief Class for computing the image of a basic set under a map.
 */

#ifndef ARIADNE_STANDARD_APPLICATOR_H
#define ARIADNE_STANDARD_APPLICATOR_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "applicator_interface.h"

namespace Ariadne {
  namespace Evaluation {


  
    /*! \brief A class for computing the image of a basic set under a map. 
     *  \ingroup Applicators
     */
    template<class R>
    class StandardApplicator
      : public ApplicatorInterface< Geometry::Rectangle<R> >,
        public ApplicatorInterface< Geometry::Zonotope<R> >
    {
      typedef Numeric::Interval<R> I;
     public:
      //@{ 
      //! \name Constructors and cloning operations.
      /*! \brief Default constructor. */
      StandardApplicator() { }
      /*! \brief Make a dynamically-allocated copy. */
      StandardApplicator<R>* clone() const { return new StandardApplicator<R>(*this); }
      //@}

      //@{ 
      //! \name Methods for applying a system to a basic set.
      
      /*! \brief Compute the image of a rectangle under a continuous function. */
      virtual Geometry::Rectangle<R> apply(const System::Map<R>& f, const Geometry::Rectangle<R>& bs) const;
      /*! \brief Compute the image of a zonotope under a continuous function. */
      virtual Geometry::Zonotope<R> apply(const System::Map<R>& f, const Geometry::Zonotope<R>& bs) const;
      //@}
    };

  }
}

#endif /* ARIADNE_STANDARD_APPLICATOR_H */
