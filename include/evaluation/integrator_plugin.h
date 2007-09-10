/***************************************************************************
 *            integrator_plugin.h
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
 
/*! \file integrator_plugin.h
 *  \brief Class for computing the image of a basic set under a map.
 */

#ifndef ARIADNE_INTEGRATOR_PLUGIN_H
#define ARIADNE_INTEGRATOR_PLUGIN_H

#include <boost/shared_ptr.hpp>

#include "../base/types.h"
#include "../base/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

#include "continuous_time_evolver.h"
#include "integrator_plugin_interface.h"
#include "integrator_plugin_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a basic set under a map. 
     *  \ingroup Integrator
     */
    template<class R>
    class IntegratorPlugin
      : public IntegratorPluginInterface<R>
    {
     public:
      //@{ 
      //! \name Constructors and cloning operations.
      /*! \brief Default constructor. */
      IntegratorPlugin();

      /*! \brief Make a dynamically-allocated copy. */
      IntegratorPlugin<R>* clone() const;
     
      //@}


      //@{ 
      //! \name Methods for applying a system to a basic set.

      /*! \brief Compute the image of a basic set under a continuous function. Returns a dynamically allocated set. */
      virtual 
      Geometry::BasicSetInterface<R>*
      integration_step(const System::VectorFieldInterface<R>& f, 
                       const Geometry::BasicSetInterface<R>& s,
                       const time_type& h,
                       const Geometry::Rectangle<R>& bb) const;

      /*! \brief Compute the image of a basic set under a continuous function. Returns a dynamically allocated set. */
      virtual 
      Geometry::BasicSetInterface<R>*
      reachability_step(const System::VectorFieldInterface<R>& f, 
                        const Geometry::BasicSetInterface<R>& s,
                        const time_type& h,
                        const Geometry::Rectangle<R>& bb) const;
    };

  }
}

#endif /* ARIADNE_INTEGRATOR_PLUGIN_H */
