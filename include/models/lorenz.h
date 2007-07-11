/***************************************************************************
 *            lorenz_system.h
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
 
/*! \file lorenz_system.h
 *  \brief The Lorenz system \f$(\dot{x},\dot{y},\dot{z}) = (\sigma(y-x),\rho x-y-xz,-\beta z+xy)\f$.
 */

#ifndef ARIADNE_LORENZ_SYSTEM_H
#define ARIADNE_LORENZ_SYSTEM_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"

#include "../system/vector_field.h"

namespace Ariadne {
  namespace System {

    /*! \brief The Lorenz system. */
    template<class R>
    class LorenzSystem : public System::VectorFieldInterface<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief Construct the Lorenz system with parameter values \a beta,
       * \a rho and \a sigma.
       */
      explicit LorenzSystem(R beta=R(8.0/3.0), 
                            R rho=R(28.0), 
                            R sigma=R(10.0))
       : _b(beta), _p(rho), _s(sigma) { }
      
      
      LorenzSystem<R>* clone() const { return new LorenzSystem<R>(this->_b, this->_p, this->_s); }
       
       
      /*! \brief  The vector field applied to a state. */
      virtual LinearAlgebra::Vector<F> image(const Geometry::Point<F>& x) const;
     
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& x) const;
            
      /*! \brief  The parameter \f$\beta\f$. */
      const R& beta() const { return _b; }
      /*! \brief  The parameter \f$\rho\f$. */
      const R& rho() const { return _p; }
      /*! \brief  The parameter \f$\sigma\f$. */
      const R& sigma() const { return _s; }
      
      
      /*! \brief  The dimension of the space. */
      dimension_type dimension() const { return 3; }
      
      /*! \brief  The smoothness of the vector field. */
      size_type smoothness() const { return std::numeric_limits<size_type>::max(); }
      
       /*! \brief  The name of the system. */
      std::string name() const { return "LorenzSystem"; }

     private:
      R _b;
      R _p;
      R _s;
    };
      
    template<class R>
    LinearAlgebra::Vector<typename LorenzSystem<R>::F>
    LorenzSystem<R>::image(const Geometry::Point<F>& x) const
    {
      LinearAlgebra::Vector<F> result(3); 
      result(0)=_s*(x[1]-x[0]);
      result(0)=_p*x[0]-x[1]-x[0]*x[2];
      result(0)=-_b*x[2]+x[0]*x[1];
      return result;
    }
     
     
    template<class R>
    LinearAlgebra::Matrix<typename LorenzSystem<R>::F>
    LorenzSystem<R>::jacobian(const Geometry::Point<F>& x) const
    {
      LinearAlgebra::Matrix<F> result(3,3); 
      result(0,0) = -_s;
      result(0,1) = _s;
      result(0,2) = 0;
      result(1,0) = _p-x[2];
      result(1,1) = -1;
      result(1,2) = -x[0];
      result(2,0) = x[1];
      result(2,1) = x[0];
      result(2,2) = -_b;
      return result;
    }
     
  
     
    template<class R>
    std::ostream& operator<<(std::ostream& os, const LorenzSystem<R>& ls) {
      os << "LorenzSystem( beta=" << ls.beta() << ", rho=" << ls.rho() << ", sigma=" << ls.sigma() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* ARIADNE_LORENZ_SYSTEM_H */
