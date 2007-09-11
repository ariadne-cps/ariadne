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

#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/integrator_interface.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief An integrator based on the Euler method.
     */
    template<class R>
    class EulerIntegrator
      : public IntegratorInterface<R>
    {
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief Constructor. */
      EulerIntegrator();


      /*! \brief Cloning operator. */
      virtual EulerIntegrator<R>* clone() const;

      /*! \brief A C0 algorithm for integrating forward a point. */
      virtual Geometry::Point<I> 
      flow_step(const System::VectorFieldInterface<R>&,
                const Geometry::Point<I>&,
                const Numeric::Interval<R>&,
                const Geometry::Rectangle<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a rectangle. */
      virtual Geometry::Rectangle<R> 
      integration_step(const System::VectorFieldInterface<R>&,
                       const Geometry::Rectangle<R>&,
                       const Numeric::Interval<R>&,
                       const Geometry::Rectangle<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a zonotope up to a certain time. */
      virtual Geometry::Rectangle<R> 
      reachability_step(const System::VectorFieldInterface<R>&,
                        const Geometry::Rectangle<R>&,
                        const Numeric::Interval<R>&,
                        const Geometry::Rectangle<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<I,R> 
      integration_step(const System::VectorFieldInterface<R>&,
                       const Geometry::Zonotope<I,R>&,
                       const Numeric::Interval<R>&,
                       const Geometry::Rectangle<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a zonotope up to a certain time. */
      virtual Geometry::Zonotope<I,R> 
      reachability_step(const System::VectorFieldInterface<R>&,
                        const Geometry::Zonotope<I,R>&,
                        const Numeric::Interval<R>&,
                        const Geometry::Rectangle<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<I,I> 
      integration_step(const System::VectorFieldInterface<R>&,
                       const Geometry::Zonotope<I,I>&,
                       const Numeric::Interval<R>&,
                       const Geometry::Rectangle<R>&) const;

      /*! \brief A C0 algorithm for integrating forward a zonotope up to a certain time. */
      virtual Geometry::Zonotope<I,I> 
      reachability_step(const System::VectorFieldInterface<R>&,
                        const Geometry::Zonotope<I,I>&,
                        const Numeric::Interval<R>&,
                        const Geometry::Rectangle<R>&) const;

    };

    
  }
}

#endif /* ARIADNE_EULER_INTEGRATOR_H */
