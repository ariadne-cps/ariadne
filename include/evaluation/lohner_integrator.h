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
#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
   

    /*!\ingroup Integrate
     * \brief An integrator based on the \f$C^1\f$-Lohner algorithm. 
     *  
     * The \f$C^1\f$-Lohner algorithm is a Taylor method.
     */
    template<class R>
    class LohnerIntegrator
      : public IntegratorBase<R, System::VectorFieldInterface<R>, Geometry::Zonotope< Numeric::Interval<R> > > 
    {
      typedef IntegratorBase<R, System::VectorFieldInterface<R>, Geometry::Zonotope< Numeric::Interval<R> > > Base_;
     public:
      typedef Numeric::Interval<R> I;
      
      /*! \brief Constructor. */
      LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      /*! \brief Cloning operator. */
      virtual LohnerIntegrator<R>* clone() const;

     public:
      
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Geometry::Point<I> bounded_flow(const System::VectorFieldInterface<R>& vf,
                                              const Geometry::Point<I>& p,
                                              const Geometry::Rectangle<R>& bb,
                                              const time_type&) const;
     
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual LinearAlgebra::Matrix<I> bounded_flow_jacobian(const System::VectorFieldInterface<R>& vf,
                                                             const Geometry::Point<I>& p,
                                                             const Geometry::Rectangle<R>& bb,
                                                             const time_type&) const;
     
       /*! \brief A C1 algorithm for integrating forward a zonotope.
       *
       * The algorithm uses a set \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes an interval Matrix \f$ \mathcal{A}_{n} \f$ such that \f$ Df(B_{n+1}) \in \mathcal{A}_{n} \f$.
       * It then computes a rectangle \f$ C_{n+1} \f$ such that \f$ \Phi(t,C_{n})\in C_{n+1} \f$.
       * We then compute \f$ \mathcal{P}_{n} \f$ such that \f$ D\Phi(h,R_{n}) \subset \mathcal{P}_{n} \f$.
       * We then compute \f$ A_{n+1} \f$ such that \f$ A_{n+1} e \supset \mathcal{P}_{n} e \f$.
       */
      virtual Geometry::Zonotope<I> 
      bounded_integration_step(const System::VectorFieldInterface<R>& vf,
                               const Geometry::Zonotope<I>& s,
                               const Geometry::Rectangle<R>& bb,
                               const time_type&) const;



      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Geometry::Zonotope<I> 
      bounded_reachability_step(const System::VectorFieldInterface<R>& vf,
                       const Geometry::Zonotope<I>& s,
                       const Geometry::Rectangle<R>& bb,
                       const time_type&) const;


    };
    
      
    
    /*!\ingroup Integrate
     * \brief An integrator based on the \f$C^1\f$-Lohner algorithm. 
     *  
     * The \f$C^1\f$-Lohner algorithm is a Taylor method.
     */
    template<class R>
    class C1LohnerIntegrator
      : public IntegratorBase<R, System::VectorFieldInterface<R>, Geometry::Zonotope< Numeric::Interval<R> > > 
    {
      typedef IntegratorBase<R, System::VectorFieldInterface<R>, Geometry::Zonotope< Numeric::Interval<R> > > Base_;
     public:
      typedef Numeric::Interval<R> I;
      
      /*! \brief Constructor. */
      C1LohnerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_set_radius);

      /*! \brief Cloning operator. */
      virtual C1LohnerIntegrator<R>* clone() const;

     public:
      
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Geometry::Point<I> bounded_flow(const System::VectorFieldInterface<R>& vf,
                                              const Geometry::Point<I>& p,
                                              const Geometry::Rectangle<R>& bb,
                                              const time_type&) const;
     

      /*! \brief Integrate a basic set for within a bounding set. */
      virtual LinearAlgebra::Matrix<I> bounded_flow_jacobian(const System::VectorFieldInterface<R>& vf,
                                                             const Geometry::Point<I>& p,
                                                             const Geometry::Rectangle<R>& bb,
                                                             const time_type&) const;
     
       /*! \brief A C1 algorithm for integrating forward a zonotope.
       *
       * The algorithm uses a set \f$B_{n+1}\f$ such that \f$R_{n+1}\subset B_{n+1}\f$. 
       * It then computes an interval Matrix \f$ \mathcal{A}_{n} \f$ such that \f$ Df(B_{n+1}) \in \mathcal{A}_{n} \f$.
       * It then computes a rectangle \f$ C_{n+1} \f$ such that \f$ \Phi(t,C_{n})\in C_{n+1} \f$.
       * We then compute \f$ \mathcal{P}_{n} \f$ such that \f$ D\Phi(h,R_{n}) \subset \mathcal{P}_{n} \f$.
       * We then compute \f$ A_{n+1} \f$ such that \f$ A_{n+1} e \supset \mathcal{P}_{n} e \f$.
       */
      virtual Geometry::Zonotope<I> 
      bounded_integration_step(const System::VectorFieldInterface<R>& vf,
                               const Geometry::Zonotope<I>& s,
                               const Geometry::Rectangle<R>& bb,
                               const time_type&) const;


      /*! \brief A C1 algorithm for integrating forward a zonotope for a time up to time \a step_size, assuming the set \a bb is a bounding box for the integration. */
      virtual Geometry::Zonotope<I> 
      bounded_reachability_step(const System::VectorFieldInterface<R>& vf,
                       const Geometry::Zonotope<I>& s,
                       const Geometry::Rectangle<R>& bb,
                       const time_type&) const;


    };
    
      
    
  }
}

#endif /* ARIADNE_LOHNER_INTEGRATOR_H */
