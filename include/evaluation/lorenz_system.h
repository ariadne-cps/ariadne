/***************************************************************************
 *            lorenz_flow.h
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

#ifndef _ARIADNE_LORENZ_FLOW_H
#define _ARIADNE_LORENZ_FLOW_H

#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"

#include "../evaluation/vector_field.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief The Lorenz system. */
    template <typename R>
    class LorenzSystem : public VectorField<R> 
    {
      typedef R Real;
      typedef Geometry::Point<R> State;
     public:
      typedef Real real_type;
      typedef State state_type;
      
      explicit LorenzSystem(Real beta=Real(8.0/3.0), 
                            Real rho=Real(28.0), 
                            Real sigma=Real(10.0))
       : _b(beta), _p(rho), _s(sigma) { }
      
      /*! \brief  The vector field applied to a state. */
      virtual LinearAlgebra::vector<R> apply(const State& x) const;
      /*! \brief  The map applied to a rectangle basic set. */
      virtual LinearAlgebra::vector< Interval<R> > apply(const Geometry::Rectangle<R>& r) const;
     
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::matrix<R> derivative(const State& x) const;
      /*! \brief  The derivative of the map over a rectangular basic set. */
      virtual LinearAlgebra::matrix< Interval<R> > derivative(const Geometry::Rectangle<R>& r) const;
            
      /*! \brief  The parameter \f$\beta$\f. */
      inline const Real& beta() const { return _b; }
      /*! \brief  The parameter \f$\rho$\f. */
      inline const Real& rho() const { return _p; }
      /*! \brief  The parameter \f$\sigma$\f. */
      inline const Real& sigma() const { return _s; }
      
      
      /*! \brief  The dimension of the space. */
      inline dimension_type dimension() const { return 3; }
      
     private:
      Real _b;
      Real _p;
      Real _s;
    };
      
    template <typename R>
    LinearAlgebra::vector<R>
    LorenzSystem<R>::apply(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::vector<R> result(3); 
      result[0]=_s*(x[1]-x[0]);
      result[1]=_p*x[0]-x[1]-x[0]*x[2];
      result[2]=-_b*x[2]+x[0]*x[1];
      return result;
    }
     
    template <typename R>
    LinearAlgebra::vector< Interval<R> >
    LorenzSystem<R>::apply(const Geometry::Rectangle<R>& X) const
    {
      LinearAlgebra::vector< Interval<R> > result(3); 
      result[0]=_s*(X[1]-X[0]);
      result[1]=_p*X[0]-X[1]-X[0]*X[2];
      result[2]=X[0]*X[1]-_b*X[2];
      return result;
    }
     
    template <typename R>
    LinearAlgebra::matrix<R>
    LorenzSystem<R>::derivative(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::matrix<R> result(3,3); 
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
     
    template <typename R>
    LinearAlgebra::matrix< Interval<R> >
    LorenzSystem<R>::derivative(const Geometry::Rectangle<R>& X) const
    {
      LinearAlgebra::matrix< Interval<R> > result(3,3); 
      result(0,0) = -_s;
      result(0,1) = _s;
      result(0,2) = 0;
      result(1,0) = _p-X[2];
      result(1,1) = -1;
      result(1,2) = -X[0];
      result(2,0) = X[1];
      result(2,1) = X[0];
      result(2,2) = -_b;
      return result;
    }
     
     
    template <typename R>
    std::ostream& operator<<(std::ostream& os, const LorenzSystem<R>& ls) {
      os << "LorenzSystem( beta=" << ls.beta() << ", rho=" << ls.rho() << ", sigma=" << ls.sigma() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* _ARIADNE_LORENZ_FLOW_H */
