/***************************************************************************
 *            kuhn_integrator.h
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
 
/*! \file kuhn_integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_KUHN_INTEGRATOR_H
#define ARIADNE_KUHN_INTEGRATOR_H

#include "numeric/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "evaluation/integrator_interface.h"
#include "evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
   

    /*!\ingroup Integrate
     * \brief An integrator based on the Kuhn cascade reduction. 
     */
    template<class R>
    class KuhnIntegrator
      : public IntegratorInterface< Geometry::Zonotope<R> >
    {
      typedef Numeric::Interval<R> I;
     public:
      
      /*! \brief Constructor. */
      KuhnIntegrator(smoothness_type temporal_order, uint cascade_size)
        : _integrator(new IntegratorBase<R>(temporal_order, 1u))
        , _cascade_size(cascade_size) { }

      /*! \brief Cloning operator. */
      virtual KuhnIntegrator<R>* clone() const { return new KuhnIntegrator<R>(*this); }

     public:
      
      /*! \brief Compute an integration time and a bounding box, given a bounding box for the intitial set, and a maximum allowable flow time. */
      virtual 
      std::pair< Numeric::Rational, Geometry::Box<R> >
      flow_bounds(const System::VectorField<R>& f, 
                  const Geometry::Box<R>& bx,
                  const Numeric::Rational& t) const; 


      /*! \brief A C1 algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<R> 
      integration_step(const System::VectorField<R>& vf,
                       const Geometry::Zonotope<R>& s,
                       const Numeric::Rational& t,
                       const Geometry::Box<R>& bb) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Geometry::Zonotope<R> 
      reachability_step(const System::VectorField<R>& vf,
                        const Geometry::Zonotope<R>& s,
                        const Numeric::Rational& t,
                        const Geometry::Box<R>& bb) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
     private:
      boost::shared_ptr< IntegratorBase<R> > _integrator;
      size_type _cascade_size;
    };
    
      
    
  }
}

#endif /* ARIADNE_KUHN_INTEGRATOR_H */
