/***************************************************************************
 *            applicator_plugin_interface.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file applicator_plugin_interface.h
 *  \brief Interface for computing the image of a basic set under a map.
 */

#ifndef ARIADNE_APPLICATOR_PLUGIN_INTERFACE_H
#define ARIADNE_APPLICATOR_PLUGIN_INTERFACE_H

#include <boost/shared_ptr.hpp>

#include "../base/types.h"
#include "../base/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a basic set under a map. 
     *  \ingroup Applicators
     */
    template<class BS>
    class ApplicatorPluginInterface
    {
      typedef typename BS::real_type R;
     public:
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ~ApplicatorPluginInterface() { }

      /*! \brief Make a dynamically-allocated copy. */
      virtual ApplicatorPluginInterface<BS>* clone() const = 0;
      
      //@}


      //@{ 
      //! \name Methods for applying a system to a basic set.

      /*! \brief Compute the image of a basic set under a continuous function. Returns a dynamically allocated set. */
      virtual 
      BS
      evaluate(const System::MapInterface<R>& f, const BS& s) const = 0;

      //@}
      
    };


  }
}

#endif /* ARIADNE_APPLICATOR_PLUGIN_INTERFACE_H */
