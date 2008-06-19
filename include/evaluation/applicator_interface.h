/***************************************************************************
 *            applicator_interface.h
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
 
/*! \file applicator_interface.h
 *  \brief Interface for computing the image of a basic set under a map.
 */

#ifndef ARIADNE_APPLICATOR_INTERFACE_H
#define ARIADNE_APPLICATOR_INTERFACE_H

#include <boost/shared_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

namespace Ariadne {
  

    /*! \ingroup EvaluatorInterfaces \ingroup Applicators
     *  \brief Interface for computing the image of an enclosure set under a map. 
     *  
     */
    template<class ES>
    class ApplicatorInterface
    {
      typedef typename ES::real_type R;
      typedef Interval<R> I;
     public:
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ~ApplicatorInterface() { }

      /*! \brief Make a dynamically-allocated copy. */
      virtual ApplicatorInterface<ES>* clone() const = 0;
      
      //@}


      //@{ 
      //! \name Methods for applying a system to a basic set.

      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ES apply(const Map<R>& f, const ES& bs) const = 0;

      /*! \brief Compute the image of a basic set under a continuous function. */
      ES operator() (const Map<R>& f, const ES& bs) const {
        return this->apply(f,bs); }

      //@}
      
    };


  
} // namespace Ariadne

#endif /* ARIADNE_APPLICATOR_INTERFACE_H */
