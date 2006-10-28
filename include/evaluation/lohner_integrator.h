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

#ifndef _ARIADNE_LOHNER_INTEGRATOR_H
#define _ARIADNE_LOHNER_INTEGRATOR_H

#include "../declarations.h"
#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
   

    /*!\ingroup Integrate
     * \brief An integrator based on the \f$C^1\f$-Lohner algorithm. 
     *  
     * The \f$C^1\f$-Lohner algorithm is a Taylor method.
     */
    template<class R>
    class LohnerIntegrator : public Integrator<R> {
      typedef Interval<R> I;
     public:
      /*! \brief Constructor. */
      LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

     public:
      /*! \brief A C1 algorithm for integrating forward a parallelotope.
       *
       * The algorithm first finds \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes an interval Matrix \f$ \mathcal{A}_{n} \f$ such that \f$ Df(B_{n+1}) \in \mathcal{A}_{n} \f$.
       * It then computes a rectangle \f$ C_{n+1} \f$ such that \f$ \Phi(t,C_{n})\in C_{n+1} \f$.
       * We then compute \f$ \mathcal{P}_{n} \f$ such that \f$ D\Phi(h,R_{n}) \subset \mathcal{P}_{n} \f$.
       * We then compute \f$ A_{n+1} \f$ such that \f$ A_{n+1} e \supset \mathcal{P}_{n} e \f$.
       */
      virtual Geometry::Parallelotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Parallelotope<R>&,
                                                          time_type&) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope.
       *
       * The algorithm first finds \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes an interval Matrix \f$ \mathcal{A}_{n} \f$ such that \f$ Df(B_{n+1}) \in \mathcal{A}_{n} \f$.
       * It then computes a rectangle \f$ C_{n+1} \f$ such that \f$ \Phi(t,C_{n})\in C_{n+1} \f$.
       * We then compute \f$ \mathcal{P}_{n} \f$ such that \f$ D\Phi(h,R_{n}) \subset \mathcal{P}_{n} \f$.
       * We then compute \f$ A_{n+1} \f$ such that \f$ A_{n+1} e \supset \mathcal{P}_{n} e \f$.
       */
      virtual Geometry::Zonotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Zonotope<R>&,
                                                          time_type&) const;

      /*! \brief A C1 algorithm for integrating forward an interval parallelotope. */
      virtual Geometry::Parallelotope< Interval<R> > 
      integration_step(const System::VectorField<R>&,
                       const Geometry::Parallelotope< Interval<R> >&,
                       time_type&) const;

      /*! \brief A C1 algorithm for integrating forward an interval zonotope. */
      virtual Geometry::Zonotope< Interval<R> > 
      integration_step(const System::VectorField<R>&,
                       const Geometry::Zonotope< Interval<R> >&,
                       time_type&) const;


      
      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      virtual Geometry::Zonotope<R> reachability_step(const System::VectorField<R>&,
                                                      const Geometry::Zonotope<R>&,
                                                      time_type& step_size) const;

      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      virtual Geometry::Zonotope< Interval<R> > reachability_step(const System::VectorField<R>&,
                                                      const Geometry::Zonotope< Interval<R> >&,
                                                      time_type& step_size) const;


      /*! \brief A specialized algorithm for integrating forward a parallelotope under an affine vector field. */
      virtual Geometry::Parallelotope<R> integration_step(const System::AffineVectorField<R>&,
                                                          const Geometry::Parallelotope<R>&,
                                                          time_type&) const;

      /*! \brief A specialized algorithm for integrating forward a zonotope under an affine vector field. */
      virtual Geometry::Zonotope<R> integration_step(const System::AffineVectorField<R>&,
                                                     const Geometry::Zonotope<R>&,
                                                     time_type&) const;

      /*! \brief A specialized algorithm for computing the reachable set of a zonotope under an affine vector field. */
      virtual Geometry::Zonotope<R> reachability_step(const System::AffineVectorField<R>&,
                                                      const Geometry::Zonotope<R>&,
                                                      time_type&) const;
                                                          
    };
    
    
  }
}

#endif /* _ARIADNE_LOHNER_INTEGRATOR_H */
