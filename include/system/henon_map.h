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
 *  \brief The Henon map \f$(x,y) \rightarrow (a-x^2-by,x)\f$.
 */

#ifndef _ARIADNE_HENON_MAP_H
#define _ARIADNE_HENON_MAP_H

#include <limits>

#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"

#include "../system/map.h"

namespace Ariadne {
  namespace System {

    /*! \brief The Henon map \f$(x,y)\mapsto(a-x^2-by,x)\f$. */
    template<class R>
    class HenonMap : public Map<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Construct the Henon map with parameters \a a and \a b. */
      explicit HenonMap(R a=R(1.5), R b=R(0.3)) : _a(a), _b(b) { }
      
      /*! \brief Returns a pointer to a dynamically-allocated copy of the map. */
      HenonMap<R>* clone() const { return new HenonMap<R>(this->_a,this->_b); }

      /*! \brief  The map applied to a state. */
      virtual Geometry::Point<F> image(const Geometry::Point<R>& x) const;
      /*! \brief  The map applied to a rectangle basic set. */
      virtual Geometry::Rectangle<R> image(const Geometry::Rectangle<R>& r) const;
      
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<R>& x) const;
      /*! \brief  The derivative of the map over a rectangular basic set. */
      virtual LinearAlgebra::Matrix< Interval<R> > jacobian(const Geometry::Rectangle<R>& r) const;
            
      /*! \brief  The parameter a. */
      const R& a() const { return _a; }
      /*! \brief  The parameter b. */
      const R& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      size_type smoothness() const { return std::numeric_limits<size_type>::max(); }
      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const { return 2; }
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const { return 2; }
      
      /*! \brief The Henon map is invertible unless \f$b=0\f$. */
      bool invertible() const {
        return _b!=0; }
        
      /*! \brief  The name of the system. */
      virtual std::string name() const { return "HenonMap"; }
     private:
      R _a;
      R _b;
    };
      
    template<class R>
    Geometry::Point<typename HenonMap<R>::F>
    HenonMap<R>::image(const Geometry::Point<R>& p) const
    {
      Geometry::Point<F> result(2); 
      const F& x=p[0];
      const F& y=p[0];
      result[0]=_a-x*x-_b*y; 
      result[1]=x; 
      return result;
    }
     
    template<class R>
    Geometry::Rectangle<R>
    HenonMap<R>::image(const Geometry::Rectangle<R>& A) const
    {
      Geometry::Rectangle<R> result(2); 
      result[0]=_a-A[0]*A[0]-_b*A[1]; 
      result[1]=A[0]; 
      return result;
    }
     
    template<class R>
    LinearAlgebra::Matrix<typename HenonMap<R>::F>
    HenonMap<R>::jacobian(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::Matrix<F> result(2,2); 
      result(0,0) = F(-2)*x[0];
      result(0,1) = R(-_b);
      result(1,0) = R(1);
      result(1,1) = R(0);
      return result;
    }
     
    template<class R>
    LinearAlgebra::Matrix< Interval<R> >
    HenonMap<R>::jacobian(const Geometry::Rectangle<R>& r) const
    {
      LinearAlgebra::Matrix< Interval<R> > result(2,2); 
      result(0,0) = R(-2)*r[0];
      result(0,1) = Interval<R>(R(-_b));
      result(1,0) = Interval<R>(R(1));
      result(1,1) = Interval<R>(R(0));
      return result;
    }
     
     
    template<class R>
    std::ostream& operator<<(std::ostream& os, const HenonMap<R>& hm) {
      os << "HenonMap( a=" << hm.a() << ", b=" << hm.b() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* _ARIADNE_HENON_MAP_H */
