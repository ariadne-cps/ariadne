/***************************************************************************
 *            euler_integrator.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file euler_integrator.h
 *  \brief Simple methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_EULER_INTEGRATOR_H
#define ARIADNE_EULER_INTEGRATOR_H

#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/integrator_interface.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief An integrator based on the Euler method.
     */
    template<class R>
    class EulerIntegrator
      : public IntegratorInterface< Geometry::Rectangle<R> >
    {
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Constructor. */
      EulerIntegrator();


      /*! \brief Cloning operator. */
      virtual EulerIntegrator<R>* clone() const;

      /*! \brief Compute an integration time and a bounding box, given a bounding box for the intitial set, and a maximum allowable flow time. */
      virtual 
      std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& f, 
                  const Geometry::Box<R>& bx,
                  const Numeric::Rational& t) const; 

      /*! \brief A C0 algorithm for integrating forward a rectangle. */
      virtual Geometry::Rectangle<R> 
      integration_step(const System::VectorField<R>&,
                       const Geometry::Rectangle<R>&,
                       const Numeric::Interval<R>&,
                       const Geometry::Box<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a rectangle up to a certain time. */
      virtual Geometry::Rectangle<R> 
      reachability_step(const System::VectorField<R>&,
                        const Geometry::Rectangle<R>&,
                        const Numeric::Interval<R>&,
                        const Geometry::Box<R>&) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
    };

    
  }
}


#endif /* ARIADNE_EULER_INTEGRATOR_H */
