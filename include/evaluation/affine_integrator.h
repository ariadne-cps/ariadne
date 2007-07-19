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

#ifndef ARIADNE_AFFINE_INTEGRATOR_H
#define ARIADNE_AFFINE_INTEGRATOR_H

#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief Compute an over-approximation to \f$\sum_{n=0}^{\infty} \frac{x^n}{(k+n)!}\f$. */
    template<class R> R gexp_up(const R& x, uint k);
  
    /*! \brief Compute an interval approximation to \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}b\f$. */
    template<class R> 
    LinearAlgebra::Vector< Numeric::Interval<R> > 
    gexp(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b, const time_type& t, const uint& k);
      
    /*! \brief Compute an interval approximation to \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}\f$. */
    template<class R> 
    LinearAlgebra::Matrix< Numeric::Interval<R> > 
    gexp(const LinearAlgebra::Matrix<R>& Ab, const time_type& t, const uint& k);
      
    
    /*!\ingroup Integrate
     * \brief An integrator based on using the exponential formula to integrate an affine vector field. 
     *  
     * The \f$C^1\f$-Affine algorithm is a Taylor method.
     */
    template<class R>
    class AffineIntegrator
      : public IntegratorBase< R, System::AffineVectorField<R>, Geometry::Zonotope< Numeric::Interval<R> > > 
    {
      typedef Numeric::Interval<R> I;
      typedef IntegratorBase< R, System::AffineVectorField<R>, Geometry::Zonotope<I> > Base_;
     public:
      /*! \brief Constructor. */
      AffineIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      /*! \brief Cloning operator. */
      virtual AffineIntegrator<R>* clone() const;

     public:

      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Geometry::Point<I> bounded_flow(const System::AffineVectorField<R>& vf,
                                              const Geometry::Point<I>& p,
                                              const Geometry::Rectangle<R>& bb,
                                              const time_type&) const;
     
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual LinearAlgebra::Matrix<I> bounded_flow_jacobian(const System::AffineVectorField<R>& vf,
                                                             const Geometry::Point<I>& p,
                                                             const Geometry::Rectangle<R>& bb,
                                                             const time_type&) const;
     


      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<I> bounded_integration_step(const System::AffineVectorField<R>& vector_field,
                                                             const Geometry::Zonotope<I>& initial_set,
                                                             const Geometry::Rectangle<R>& bounding_set,
                                                             const time_type& step_size) const;

      
      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. 
       */
      virtual Geometry::Zonotope<I> bounded_reachability_step(const System::AffineVectorField<R>& vector_field,
                                                              const Geometry::Zonotope<I>& initial_set, 
                                                              const Geometry::Rectangle<R>& bounding_set,
                                                              const time_type& step_size) const;


      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope. 
      *
      *  Overrides method in Integrator since we don't use the bounding set. 
      */
      virtual Geometry::Zonotope<I> integration_step(const System::AffineVectorField<R>& vector_field,
                                                     const Geometry::Zonotope<I>& initial_set,
                                                     time_type& step_size) const;


      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size.
       *
       *  Overrides method in Integrator since we don't use the bounding set. 
       */
      virtual Geometry::Zonotope<I> reachability_step(const System::AffineVectorField<R>& vector_field,
                                                      const Geometry::Zonotope<I>& initial_set,
                                                      time_type& step_size) const;



     };

     
  }
}

#endif /* ARIADNE_AFFINE_INTEGRATOR_H */
