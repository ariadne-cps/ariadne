/***************************************************************************
 *            standard_integrator.h
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
 

#ifndef ARIADNE_STANDARD_INTEGRATOR_H
#define ARIADNE_STANDARD_INTEGRATOR_H

#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "system/vector_field.h"

#include "evaluation/bounder_interface.h"
#include "evaluation/flower_interface.h"
#include "evaluation/reducer_interface.h"
#include "evaluation/integrator_interface.h"

#include "evaluation/integrator_base.h"


namespace Ariadne { 
   

    template<class R> class Flower;

    /*!\ingroup Integrators
     * \brief A template for computing a step of the evolution of an enclosure set under a vector field.
     */
    template<class ES> class StandardIntegrator;
  
    /*!\ingroup Integrators
     * \brief A template for computing a step of the evolution of a zonotope under a vector field.
     *
     * The update rule for an integration step on a zonotope is
     * \f[  \Phi(X,t) \subset c+tf(c)+\frac{t^2}{2} Df(B)f(B) + \bigl(I+t\,Df(X)\bigr)\cdot(x-c) . \f]
     * where \f$B\f$ is a bound for \f$\Phi(c,[0,t])\f$.
     * Subdivision is performed using orthogonal over-approximation.
     *
     * The update rule for an integration step for an interval zonotope is 
     * \f[  \Phi(x,t) \subset c+tf(c)+\frac{t^2}{2}\,Df(B_c)\,f(B_c) + \bigl(I+t\,Df(B)\,W\bigr)\cdot(x-c) \f]
     * where \f$B_c\f$ is a bound for \f$\Phi(c,[0,t])\f$, \f$B\f$ is a bound for \f$\Phi(X,[0,t])\f$ and \f$W\f$ is a bound for \f$D\Phi(X,[0,t])\f$.
     * Subdivision is performed using orthogonal over-approximation.
     * See the section on the \ref c1lohnerintegrator for details.
     */
    template<class R> class StandardIntegrator< Zonotope<R> >
      : public IntegratorBase< Zonotope<R> >
    {
      typedef Interval<R> I;
      typedef Zonotope<R> ES;
     private:
      boost::shared_ptr< BounderInterface<R> > _bounder;
      boost::shared_ptr< FlowerInterface<R> > _flower;
     public:
      
      /*! \brief Construct and integrator with a given temporal order. */
      StandardIntegrator();

      /*! \brief Construct an integrator from a bounding box method, an method for computing a flow, and an enclosure reduction method. */
      StandardIntegrator(const BounderInterface<R>& bounder, 
                         const FlowerInterface<R>& flower);

      /*! \brief Cloning operator. */
      virtual StandardIntegrator< Zonotope<R> >* clone() const;

     public:
      
      /*! \brief Compute an integration time and a bounding box, given a bounding box for the intitial set, and a maximum allowable flow time. */
      virtual 
      std::pair< Rational, Box<R> >
      flow_bounds(const VectorField<R>& vector_field, 
                  const Zonotope<R>& initial_set,
                  const Rational& maximum_step_size) const; 

      /*! \brief An algorithm for integrating forward a zonotope.  */
      virtual Zonotope<R> 
      integration_step(const VectorField<R>& vector_field,
                       const Zonotope<R>& initial_set,
                       const Rational& step_size,
                       const Box<R>& bounding_set) const;

      /*! \brief An algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Zonotope<R> 
      reachability_step(const VectorField<R>& vector_field,
                        const Zonotope<R>& initial_set,
                        const Rational& step_size,
                        const Box<R>& bounding_set) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
    };
    


  
} // namespace Ariadne

#endif /* STANDARD_INTEGRATOR_H */
