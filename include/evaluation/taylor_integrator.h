/***************************************************************************
 *            taylor_integrator.h
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
 
/*! \file taylor_integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_TAYLOR_INTEGRATOR_H
#define ARIADNE_TAYLOR_INTEGRATOR_H

#include "numeric/declarations.h"
#include "geometry/declarations.h"
#include "system/taylor_flow.h"
#include "evaluation/integrator_interface.h"

namespace Ariadne {
  namespace Evaluation {
   

    /*!\ingroup Integrate
     * \brief An integrator based on the Taylor algorithm.
     */
    template<class R>
    class TaylorIntegrator
      : public IntegratorInterface< Geometry::Zonotope< Numeric::Interval<R> > >
    {
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I,I> BS;
     public:
      
      /*! \brief Constructor. */
      TaylorIntegrator();

      /*! \brief Cloning operator. */
      virtual TaylorIntegrator<R>* clone() const;

     public:
      
      
      /*! \brief Compute the flow map of a vector field. */
      virtual System::TaylorFlow<I> flow(const System::VectorField<R>& vf, 
                                         const Geometry::Box<R>& bb) const;

      /*! \brief Integrate a basic set for time \a t within a bounding set. */
      virtual Geometry::Point<I> flow_step(const System::VectorField<R>& vf,
                                           const Geometry::Point<I>& p,
                                           const Numeric::Interval<R>& t,
                                           const Geometry::Box<R>& bb) const;
     
      /*! \brief An algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<I,I> 
      integration_step(const System::VectorField<R>& vf,
                       const Geometry::Zonotope<I,I>& s,
                       const Numeric::Interval<R>& t,
                       const Geometry::Box<R>& bb) const;

      /*! \brief An algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Geometry::Zonotope<I,I> 
      reachability_step(const System::VectorField<R>& vf,
                        const Geometry::Zonotope<I,I>& s,
                        const Numeric::Interval<R>& t,
                        const Geometry::Box<R>& bb) const;


      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
    };
    
      
    

    
  }
}

#endif /* ARIADNE_TAYLOR_INTEGRATOR_H */
