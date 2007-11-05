/***************************************************************************
 *            lohner_integrator.h
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
 
/*! \file lohner_integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_LOHNER_INTEGRATOR_H
#define ARIADNE_LOHNER_INTEGRATOR_H

#include "../numeric/declarations.h"
#include "../geometry/declarations.h"
#include "../evaluation/integrator_interface.h"

namespace Ariadne {
  namespace Evaluation {
   

    /*!\ingroup Integrate
     * \brief An integrator based on the Lohner algorithm.
     *
     * The update rule for an integration step on a zonotope is
     * \f[  \Phi(X,t) \subset c+tf(c)+\frac{t^2}{2} Df(B)f(B) + \bigl(I+t\,Df(X)\bigr)\cdot(x-c) . \f]
     * where \f$B\f$ is a bound for \f$\Phi(c,[0,t])\f$.
     * Subdivision is performed using orthogonal over-approximation.
     *
     * See the section on the \ref c1lohnerintegrator for details.
     */
    template<class R>
    class LohnerIntegrator
      : public IntegratorInterface< Geometry::Zonotope<Numeric::Interval<R>,R> >
    {
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I,R> BS;
     public:
      
      /*! \brief Constructor. */
      LohnerIntegrator();

      /*! \brief Cloning operator. */
      virtual LohnerIntegrator<R>* clone() const;

     public:
      
      /*! \brief Integrate a basic set for time \a t within a bounding set. */
      virtual Geometry::Point<I> flow_step(const System::VectorFieldInterface<R>& vf,
                                           const Geometry::Point<I>& p,
                                           const Numeric::Interval<R>& t,
                                           const Geometry::Rectangle<R>& bb) const;
     
      /*! \brief A C1 algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<I,R> 
      integration_step(const System::VectorFieldInterface<R>& vf,
                       const Geometry::Zonotope<I,R>& s,
                       const Numeric::Interval<R>& t,
                       const Geometry::Rectangle<R>& bb) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Geometry::Zonotope<I,R> 
      reachability_step(const System::VectorFieldInterface<R>& vf,
                        const Geometry::Zonotope<I,R>& s,
                        const Numeric::Interval<R>& t,
                        const Geometry::Rectangle<R>& bb) const;


      /*! \brief A C1 algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<I,I> 
      integration_step(const System::VectorFieldInterface<R>& vf,
                       const Geometry::Zonotope<I,I>& s,
                       const Numeric::Interval<R>& t,
                       const Geometry::Rectangle<R>& bb) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<I,I> 
      reachability_step(const System::VectorFieldInterface<R>& vf,
                       const Geometry::Zonotope<I,I>& s,
                       const Numeric::Interval<R>& t,
                       const Geometry::Rectangle<R>& bb) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
    };
    
      
    
    /*!\ingroup Integrate
     * \brief An integrator based on the C<sup>1</sup>-Lohner algorithm on zonotopes. 
     *
     * The update rule for an integration step is 
     * \f[  \Phi(x,t) \subset c+tf(c)+\frac{t^2}{2}\,Df(B_c)\,f(B_c) + \bigl(I+t\,Df(B)\,W\bigr)\cdot(x-c) \f]
     * where \f$B_c\f$ is a bound for \f$\Phi(c,[0,t])\f$, \f$B\f$ is a bound for \f$\Phi(X,[0,t])\f$ and \f$W\f$ is a bound for \f$D\Phi(X,[0,t])\f$.
     * Subdivision is performed using orthogonal over-approximation.
     *
     * See the section \ref c1lohnerintegrator for details.
     */
    template<class R>
    class C1LohnerIntegrator
      : public IntegratorInterface< Geometry::Zonotope<Numeric::Interval<R> > > 
    {
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I,I> BS;
     public:
      /*! \brief Constructor. */
      C1LohnerIntegrator();

      /*! \brief Cloning operator. */
      virtual C1LohnerIntegrator<R>* clone() const;

     public:
      
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Geometry::Point<I> flow_step(const System::VectorFieldInterface<R>& vf,
                                           const Geometry::Point<I>& p,
                                           const Numeric::Interval<R>& t,
                                           const Geometry::Rectangle<R>& bb) const;
     

      /*! \brief Integrate a basic set for within a bounding set. */
      virtual LinearAlgebra::Matrix<I> flow_step_jacobian(const System::VectorFieldInterface<R>& vf,
                                                          const Geometry::Point<I>& p,
                                                          const Numeric::Interval<R>& t,
                                                          const Geometry::Rectangle<R>& bb) const;
     
      /*! \brief A C1 algorithm for integrating forward a zonotope.
       */
      virtual Geometry::Zonotope<I,I> 
      integration_step(const System::VectorFieldInterface<R>& vf,
                       const Geometry::Zonotope<I,I>& s,
                       const Numeric::Interval<R>& t,
                       const Geometry::Rectangle<R>& bb) const;


      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Geometry::Zonotope<I,I> 
      reachability_step(const System::VectorFieldInterface<R>& vf,
                        const Geometry::Zonotope<I,I>& s,
                        const Numeric::Interval<R>& t,
                        const Geometry::Rectangle<R>& bb) const;


      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;
    };
    
      
    
  }
}

#endif /* ARIADNE_LOHNER_INTEGRATOR_H */
