/***************************************************************************
 *            henon_map.h
 *
 *  Wed Feb  2 18:52:36 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
/*! \file henon_map.h
 *  \brief The Henon map \f$(x,y) \hookrightarrow (a-x^2-by,x)\f$.
 */

#ifndef _ARIADNE_HENON_MAP_H
#define _ARIADNE_HENON_MAP_H

#include "linear_algebra.h"

#include "point.h"
#include "rectangle.h"
#include "parallelopiped.h"

#include "map.h"
#include "apply.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief An affine map on Euclidean space. */
    template <typename R>
    class HenonMap : public Map<R> 
    {
     public:
      typedef R Real;
      typedef Geometry::Point<Real> State;
      
      inline explicit HenonMap(Real a=Real(1.5),Real b=Real(0.3)) : _a(a), _b(b) { }
      
      /*! \brief  The map applied to a state. */
      virtual State apply(const State& x) const;
      /*! \brief  The map applied to a rectangle basic set. */
      virtual Geometry::Rectangle<R> apply(const Geometry::Rectangle<R>& r) const;
      /*! \brief  The map applied to a parallelopiped basic set. */
      virtual Geometry::Parallelopiped<R> apply(const Geometry::Parallelopiped<R>& p) const;
      
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::matrix<R> derivative(const State& x) const;
      /*! \brief  The derivative of the map over a rectangular basic set. */
      virtual LinearAlgebra::matrix< Interval<R> > derivative(const Geometry::Rectangle<R>& r) const;
            
      /*! \brief  The parameter a. */
      inline const Real& a() const { return _a; }
      /*! \brief  The parameter b. */
      inline const Real& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      inline dimension_type argument_dimension() const { return 2; }
      /*! \brief The dimension of the result. */
      inline dimension_type result_dimension() const { return 2; }
      
      inline bool invertible() const {
        return _b!=0; }
     private:
      Real _a;
      Real _b;
    };
      
    template <typename R>
    Geometry::Point<R>
    HenonMap<R>::apply(const Geometry::Point<R>& x) const
    {
      State result(2); 
      result[0]=_a-x[0]*x[0]-_b*x[1]; 
      result[1]=x[0]; 
      return result;
    }
     
    template <typename R>
    Geometry::Rectangle<R>
    HenonMap<R>::apply(const Geometry::Rectangle<R>& A) const
    {
      Geometry::Rectangle<R> result(2); 
      result[0]=_a-A[0]*A[0]-_b*A[1]; 
      result[1]=A[0]; 
      return result;
    }
     
    template <typename R>
    Geometry::Parallelopiped<R>
    HenonMap<R>::apply(const Geometry::Parallelopiped<R>& p) const
    {
      std::cerr << "HenonMap::apply(const Parallelopiped) const" << std::endl;
      return Evaluation::apply(*this,p);
      std::cerr << "Done HenonMap::apply(const Parallelopiped) const" << std::endl;
    }

    template <typename R>
    LinearAlgebra::matrix<R>
    HenonMap<R>::derivative(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::matrix<R> result(2,2); 
      result(0,0) = -2*x[0];
      result(0,1) = -_b;
      result(1,0) = 1;
      result(1,1) = 0;
      return result;
    }
     
    template <typename R>
    LinearAlgebra::matrix< Interval<R> >
    HenonMap<R>::derivative(const Geometry::Rectangle<R>& r) const
    {
      LinearAlgebra::matrix< Interval<R> > result(2,2); 
      result(0,0) = R(-2)*r[0];
      result(0,1) = Interval<R>(-_b);
      result(1,0) = Interval<R>(1);
      result(1,1) = Interval<R>(0);
      return result;
    }
     
     
    template <typename R>
    std::ostream& operator<<(std::ostream& os, const HenonMap<R>& hm) {
      os << "HenonMap( a=" << hm.a() << ", b=" << hm.b() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* _ARIADNE_AFFINE_MAP_H */
