/***************************************************************************
 *            kuhn_applicator.h
 *
 *  Copyright  2006-7  Pieter Collins
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
 
/*! \file kuhn_applicator.h
 *  \brief Methods for computing the image of a zonotope under a continuous map.
 */

#ifndef ARIADNE_KUHN_APPLICATOR_H
#define ARIADNE_KUHN_APPLICATOR_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "evaluation/applicator_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a basic set under a map. 
     *  \ingroup Applicators
     */
    template<class R>
    class KuhnApplicator
      : public ApplicatorInterface< Geometry::Zonotope<R> >
    {
      typedef Geometry::Zonotope<R> BS;
     public:
      /*! \brief Construct with a parameter specifying the cascade size. */
      KuhnApplicator(const size_type& cascade_size) : _cascade_size(cascade_size) { }

      /*! \brief Make a dynamically-allocated copy. */
      virtual KuhnApplicator<R>* clone() const { return new KuhnApplicator<R>(*this); }

      /*! \brief Compute the image of a zonotope set under a differentiable function. */
      virtual Geometry::Zonotope<R> apply(const System::MapInterface<R>& f, const Geometry::Zonotope<R>& s) const;
     private:
      size_type _cascade_size;
    };


  }
}

#endif /* ARIADNE_KUHN_APPLICATOR_H */
