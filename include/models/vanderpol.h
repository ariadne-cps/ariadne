/***************************************************************************
 *            vanderpol.h
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

/*! \file vanderpol.h
 *  \brief The van der Pol equation \f$\ddot{y}-\mu(1-y^2)\dot{y}+y=0\f$.
 */

#ifndef _ARIADNE_VANDERPOL_EQUATION_H
#define _ARIADNE_VANDERPOL_EQUATION_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"

#include "../system/vector_field.h"

namespace Ariadne {
  namespace System {

    /*! \brief The van der Pol equation. */
    template<class R>
    class VanDerPolEquation : public System::VectorField<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef Interval<R> I;
     public:
      /*! \brief Constructor. */
      explicit VanDerPolEquation(R mu=0.0)
        : _mu(mu) { }
      
      VanDerPolEquation<R>* clone() const { return new VanDerPolEquation<R>(this->_mu); }
       
      /*! \brief  The vector field applied to a state. */
      virtual LinearAlgebra::Vector<F> operator() (const Geometry::Point<R>& x) const;
      /*! \brief  The map applied to a rectangle basic set. */
      virtual LinearAlgebra::Vector<I> operator() (const Geometry::Rectangle<R>& r) const;
     
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<R>& x) const;
      /*! \brief  The derivative of the map over a rectangular basic set. */
      virtual LinearAlgebra::Matrix<I> jacobian(const Geometry::Rectangle<R>& r) const;
            
      /*! \brief  The parameter \f$\mu\f$. */
      const R& mu() const { return _mu; }
      
      
      /*! \brief  The dimension of the space. */
      dimension_type dimension() const { return 2; }
      
      /*! \brief  The smoothness of the vector field. */
      size_type smoothness() const { return std::numeric_limits<size_type>::max(); }
      
       /*! \brief  The name of the system. */
      std::string name() const { return "VanDerPolEquation"; }

     private:
      R _mu;
    };
      
    template<class R>
    LinearAlgebra::Vector<typename VanDerPolEquation<R>::F>
    VanDerPolEquation<R>::operator() (const Geometry::Point<R>& x) const
    {
      LinearAlgebra::Vector<F> result(2); 
      result(0)=x[1];
      result(1)=_mu*(1-x[0]*x[0])*x[1]-x[0];
      return result;
    }
     
    template<class R>
    LinearAlgebra::Vector< Interval<R> >
    VanDerPolEquation<R>::operator() (const Geometry::Rectangle<R>& x) const
    {
      LinearAlgebra::Vector< Interval<R> > result(2); 
      result(0)=x[1];
      result(1)=_mu*(1.0-x[0]*x[0])*x[1]-x[0];
      return result;
    }
     
    template<class R>
    LinearAlgebra::Matrix<typename VanDerPolEquation<R>::F>
    VanDerPolEquation<R>::jacobian(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::Matrix<F> result(2,2); 
      result(0,0) = 0.0;
      result(0,1) = 1.0;
      result(1,0) = -2.0*_mu*x[0]*x[1]-1.0;
      result(1,1) = _mu*(1.0-x[0]*x[0]);
      return result;
    }
     
    template<class R>
    LinearAlgebra::Matrix< Interval<R> >
    VanDerPolEquation<R>::jacobian(const Geometry::Rectangle<R>& x) const
    {
      LinearAlgebra::Matrix< Interval<R> > result(2,2); 
      result(0,0) = 0.0;
      result(0,1) = 1.0;
      result(1,0) = -2.0*_mu*x[0]*x[1]-1.0;
      result(1,1) = _mu*(1.0-x[0]*x[0]);
      return result;
    }
     
     
    template<class R>
    std::ostream& operator<<(std::ostream& os, const VanDerPolEquation<R>& vdp) {
      os << "VanDerPolEquation( mu=" << vdp.mu() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* _ARIADNE_DUFFING_EQUATION_H */
