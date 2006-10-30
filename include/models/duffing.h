/***************************************************************************
 *            duffing.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.itm, Pieter.Collins@cwi.nl
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

/*! \file duffing.h
 *  \brief The Duffing equation \f$(\ddot{x}+\delta\dot{x}+(\beta x^3\pm\omega_0^2)=\gamma\cos(\omega t+\phi)\f$.
 */

#ifndef _ARIADNE_DUFFING_EQUATION_H
#define _ARIADNE_DUFFING_EQUATION_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"

#include "../system/vector_field.h"

namespace Ariadne {
  namespace System {

    /*! \brief The Duffing equation. */
    template<class R>
    class DuffingEquation : public System::VectorField<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief Construct the unforced Duffing system with parameter values \a delta,
       * \a beta and \a alpha.
       */
      explicit DuffingEquation(R delta=0.0, 
                               R beta=1.0, 
                               R alpha=1.0)
        : _delta(delta), _beta(beta), _alpha(alpha),
          _gamma(0.0), _omega(1.0), _phi(0.0) { }
      
      
      DuffingEquation<R>* clone() const { return new DuffingEquation<R>(this->_delta, this->_beta, this->_alpha); }
       
       
      /*! \brief  The vector field applied to a state. */
      virtual LinearAlgebra::Vector<F> operator() (const Geometry::Point<R>& x) const;
      /*! \brief  The map applied to a rectangle basic set. */
      virtual LinearAlgebra::Vector< Interval<R> > operator() (const Geometry::Rectangle<R>& r) const;
     
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<R>& x) const;
      /*! \brief  The derivative of the map over a rectangular basic set. */
      virtual LinearAlgebra::Matrix< Interval<R> > jacobian(const Geometry::Rectangle<R>& r) const;
            
      /*! \brief  The parameter \f$\delta\f$. */
      const R& delta() const { return _delta; }
      /*! \brief  The parameter \f$\beta\f$. */
      const R& beta() const { return _beta; }
      /*! \brief  The parameter \f$\alpha\f$. */
      const R& alpha() const { return _alpha; }
      
      
      /*! \brief  The dimension of the space. */
      dimension_type dimension() const { return 2; }
      
      /*! \brief  The smoothness of the vector field. */
      size_type smoothness() const { return std::numeric_limits<size_type>::max(); }
      
       /*! \brief  The name of the system. */
      std::string name() const { return "DuffingEquations"; }

     private:
      R _delta, _beta, _alpha, _gamma, _omega,_phi;
    };
      
    template<class R>
    LinearAlgebra::Vector<typename DuffingEquation<R>::F>
    DuffingEquation<R>::operator() (const Geometry::Point<R>& x) const
    {
      LinearAlgebra::Vector<F> result(2); 
      result(0)=x[1];
      result(1)=-_delta*x[1]-x[0]*(_alpha+_beta*x[0]*x[0]);
      return result;
    }
     
    template<class R>
    LinearAlgebra::Vector< Interval<R> >
    DuffingEquation<R>::operator() (const Geometry::Rectangle<R>& x) const
    {
      LinearAlgebra::Vector< Interval<R> > result(2); 
      result(0)=x[1];
      result(1)=-_delta*x[1]-x[0]*(_alpha+_beta*x[0]*x[0]);
      return result;
    }
     
    template<class R>
    LinearAlgebra::Matrix<typename DuffingEquation<R>::F>
    DuffingEquation<R>::jacobian(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::Matrix<F> result(2,2); 
      result(0,0) = 0;
      result(0,1) = 1;
      result(1,0) = -(_alpha+3.0*_beta*x[0]*x[0]);
      result(1,1) = -_delta;
      return result;
    }
     
    template<class R>
    LinearAlgebra::Matrix< Interval<R> >
    DuffingEquation<R>::jacobian(const Geometry::Rectangle<R>& x) const
    {
      LinearAlgebra::Matrix< Interval<R> > result(2,2); 
      result(0,0) = 0;
      result(0,1) = 1;
      result(1,0) = -(_alpha+3.0*_beta*x[0]*x[0]);
      result(1,1) = -_delta;
      return result;
    }
     
     
    template<class R>
    std::ostream& operator<<(std::ostream& os, const DuffingEquation<R>& de) {
      os << "DuffingEquation( delta=" << de.delta() << ", beta=" << de.beta() << ", alpha=" << de.alpha() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* _ARIADNE_DUFFING_EQUATION_H */
