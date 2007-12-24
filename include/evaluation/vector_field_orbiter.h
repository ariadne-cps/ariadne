/***************************************************************************
 *            vector_field_orbiter.h
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
 
/*! \file vector_field_orbiter.h
 *  \brief Methods for computing orbits of boxes.
 */

#ifndef ARIADNE_VECTOR_FIELD_ORBITER_H
#define ARIADNE_VECTOR_FIELD_ORBITER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/bounder_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/orbiter_interface.h"

namespace Ariadne {
  namespace Evaluation {

      
    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class BS>
    class VectorFieldOrbiter 
      : public VectorFieldOrbiterInterface<typename BS::real_type>
    {
      typedef typename BS::real_type R;
     public:
      /*! \brief Construct from evolution parameters and an integrator. */
      VectorFieldOrbiter<BS>(const EvolutionParameters<R>&, const IntegratorInterface<BS>&);

      /*! \brief Copy constructor. */
      VectorFieldOrbiter<BS>(const VectorFieldOrbiter<BS>&);

      /*! \brief Make a dynamically-allocated copy. */
      VectorFieldOrbiter<BS>* clone() const;



      /*! \brief Compute the orbit of a rectangle for at most time \a t. */
      virtual 
      Geometry::ContinuousTimeOrbit< Numeric::Rational, Geometry::Box<R>, Geometry::Box<R>  >
      orbit(const System::VectorFieldInterface<R>& f, const Geometry::Box<R>& r, const Numeric::Rational& t) const;

      /*! \brief Compute the orbit of a grid cell for time \a t under a vector field. */
      virtual 
      Geometry::ContinuousTimeOrbit< Numeric::Rational, Geometry::GridCellListSet<R>, Geometry::GridCellListSet<R> >
      orbit(const System::VectorFieldInterface<R>& vf, const Geometry::GridCell<R>& gc, const Numeric::Rational& t) const;

     private:
      // Convenience wrapper for integrator services
     private:
      boost::shared_ptr< BounderInterface<R> > _bounder;
      boost::shared_ptr< IntegratorInterface<BS> > _integrator;
      boost::shared_ptr< ApproximatorInterface<BS> > _approximator;
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
    };

  }

}

#endif /* ARIADNE_VECTOR_FIELD_ORBITER_H */
