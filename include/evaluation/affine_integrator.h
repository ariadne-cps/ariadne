/***************************************************************************
 *            affine_integrator.h
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
 
/*! \file affine_integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef _ARIADNE_AFFINE_INTEGRATOR_H
#define _ARIADNE_AFFINE_INTEGRATOR_H

#include "../declarations.h"
#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief Compute an over-approximation to \f$\sum_{n=0}^{\infty} \frac{x^n}{(k+n)!}\f$. */
    template<class R> R gexp_up(const R& x, uint k);
      
    
    /*! \brief Compute \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}b\f$ with an error of at most \a epsilon. */
    template<class R> 
    LinearAlgebra::Vector< Numeric::Interval<R> > 
    gexp(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b, const time_type& t, const uint& k, const R& err);
      
    
    /*!\ingroup Integrate
     * \brief An integrator based on using the exponential formula to integrate an affine vector field. 
     *  
     * The \f$C^1\f$-Affine algorithm is a Taylor method.
     */
    template<class R>
    class AffineIntegrator : public Integrator<R> {
      typedef Interval<R> I;
     public:
      /*! \brief Constructor. */
      AffineIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

     public:
      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a parallelotope. */
      virtual Geometry::Parallelotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Parallelotope<R>&,
                                                          time_type&) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<R> integration_step(const System::VectorField<R>&,
                                                          const Geometry::Zonotope<R>&,
                                                          time_type&) const;

      
      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      virtual Geometry::Zonotope<R> reachability_step(const System::VectorField<R>&,
                                                      const Geometry::Zonotope<R>&,
                                                      time_type& step_size) const;
     public:
      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope< Interval<R> > integration_step(const System::AffineVectorField<R>&,
                                                                 const Geometry::Zonotope< Interval<R> >&,
                                                                 time_type&) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<R> integration_step(const System::AffineVectorField<R>&,
                                                          const Geometry::Zonotope<R>&,
                                                          time_type&) const;

      
      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. 
       *
       * We use the formula that the reached set is 
       * \f\[ c + Ge + t(Ac+b) + 
       */
      virtual Geometry::Zonotope<R> reachability_step(const System::AffineVectorField<R>&,
                                                      const Geometry::Zonotope<R>&,
                                                      time_type& step_size) const;

     };

     
  }
}

#endif /* _ARIADNE_LOHNER_INTEGRATOR_H */
