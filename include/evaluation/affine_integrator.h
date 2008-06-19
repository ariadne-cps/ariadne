/***************************************************************************
 *            affine_integrator.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file affine_integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_AFFINE_INTEGRATOR_H
#define ARIADNE_AFFINE_INTEGRATOR_H

#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/integrator_interface.h"

#include "system/affine_vector_field.h"

namespace Ariadne {
  
   
    /*! \brief Compute an over-approximation to \f$\sum_{n=0}^{\infty} \frac{x^n}{(k+n)!}\f$. */
    template<class R> R gexp_up(const R& x, uint k);
  
    /*! \brief Compute an interval approximation to \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}b\f$. */
    template<class R> 
    Vector< Interval<R> > 
    gexp(const Matrix< Interval<R> >& A, 
         const Vector< Interval<R> >& b, 
         const Interval<R>& t, const uint& k);
      
    /*! \brief Compute an interval approximation to \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}\f$. */
    template<class R> 
    Matrix< Interval<R> > 
    gexp(const Matrix< Interval<R> >& Ab, const Interval<R>& t, const uint& k);
      
    template<class BS> class AffineIntegrator;
    
    /*!\ingroup Integrators
     * \brief An integrator based on using the exponential formula to integrate an affine vector field. 
     *  
     * The \f$C^1\f$-Affine algorithm is a Taylor method.
     */
    template<class R>
    class AffineIntegrator< Zonotope<R> >
      : public IntegratorInterface< Zonotope<R> >      
    {
      typedef Interval<R> I;
     public:
      /*! \brief Constructor. */
      AffineIntegrator();

      /*! \brief Cloning operator. */
      virtual AffineIntegrator< Zonotope<R> >* clone() const;

     public:

      /*! \brief Compute an integration time and a bounding box, given a bounding box for the intitial set, and a maximum allowable flow time. */
      virtual 
      std::pair< Rational, Box<R> >
      flow_bounds(const VectorField<R>& vector_field, 
                  const Zonotope<R>& initial_set,
                  const Rational& maximum_step_size) const; 

      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Point<I> flow_step(const VectorField<R>& vector_field,
                                           const Point<I>& initial_point,
                                           const Rational& step_size,
                                           const Box<R>& bounding_set) const;
     
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Matrix<I> flow_step_jacobian(const VectorField<R>& vector_field,
                                                          const Point<I>& initial_point,
                                                          const Rational& step_size,
                                                          const Box<R>& bounding_set) const;
     

      virtual Zonotope<R> 
      integration_step(const VectorField<R>& vector_field,
                       const Zonotope<R>& initial_set,
                       const Rational& step_size,
                       const Box<R>& bounding_set) const;
      
      virtual Zonotope<R> 
      reachability_step(const VectorField<R>& affine_vector_field,
                        const Zonotope<R>& initial_set,
                        const Rational& step_size,
                        const Box<R>& bounding_set) const;

     public:
      /*! \brief Integrate \a initial point for time \a step_size. */
      Point<I> 
      flow_step(const AffineVectorField<R>& affine_vector_field,
                const Point<I>& initial_point,
                const Rational& step_size) const;

      /*! \brief Comput the spacial Jacobian derivative at \a initial point for time \a step_size. */
      Matrix<I> 
      flow_step_jacobian(const AffineVectorField<R>& vector_field,
                         const Point<I>& affine_initial_point,
                         const Rational& step_size) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      Zonotope<R> 
      integration_step(const AffineVectorField<R>& affine_vector_field,
                       const Zonotope<R>& initial_set,
                       const Rational& step_size) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      Zonotope<R> 
      reachability_step(const AffineVectorField<R>& affine_vector_field,
                        const Zonotope<R>& initial_set,
                        const Rational& step_size) const;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream&) const;

     private:
      const AffineVectorField<R>* cast(const VectorField<R>*) const;
    };

     
  
} // namespace Ariadne

#endif /* ARIADNE_AFFINE_INTEGRATOR_H */
