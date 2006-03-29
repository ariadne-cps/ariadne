/***************************************************************************
 *            affine_map.h
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
 
/*! \file affine_map.h
 *  \brief Affine maps of the form \f$x\rightarrow Ax+b\f$.
 */

#ifndef _ARIADNE_AFFINE_MAP_H
#define _ARIADNE_AFFINE_MAP_H

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../evaluation/map.h"

#include "../geometry/point.h"
#include "../geometry/simplex.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief An affine map on Euclidean space. */
    template <typename R>
    class AffineMap // : public Map<R,Geometry::Point> 
    {
     public:
      typedef R Real;
      typedef Geometry::Point<Real> State;
      
      typedef Ariadne::LinearAlgebra::matrix<Real> Matrix;
      typedef Ariadne::LinearAlgebra::vector<Real> Vector;
      
      inline explicit AffineMap() {}
      inline explicit AffineMap(const Matrix& A, const Vector& b) : _A(A), _b(b) { }
      inline explicit AffineMap(const Matrix& A) : _A(A), _b(A.columns()) { }
      inline explicit AffineMap(const Vector& b) : _A(b.size(),b.size()), _b(b) { }
      
      inline AffineMap(const AffineMap<Real>& T) : _A(T._A), _b(T._b) { }
      inline AffineMap<Real>& operator=(const AffineMap<Real>& T) {
        this->_A=T._A; this->_b=T._b; return *this; }
      
      /*! \brief  The map applied to a state. */
      State operator() (const State& x) const;
        
      /*! \brief  The map applied to a simplex basic set. */
      Geometry::Simplex<R> operator() (const Geometry::Simplex<R>& A) const;
      
      /*! \brief  The map applied to a rectangle basic set. */
      Geometry::Rectangle<R> operator() (const Geometry::Rectangle<R>& A) const;
      
      /*! \brief  The map applied to a parallelotope basic set. */
      Geometry::Parallelotope<R> operator() (const Geometry::Parallelotope<R>& A) const;
      
      /*! \brief  The map applied to a zonotope basic set. */
      Geometry::Zonotope<R> operator() (const Geometry::Zonotope<R>& A) const;
      
      /*! \brief  The map applied to a polyhedron basic set. */
      Geometry::Polyhedron<R> operator() (const Geometry::Polyhedron<R>& A) const;
      
      /*! \brief  The linear transformation of the map. */
      inline const Matrix& A() const { return _A; }
      /*! \brief  The offset vector of the map. */
      inline const Vector& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      inline size_type argument_dimension() const {
        return _A.size2();
      }
      
      /*! \brief The dimension of the result. */
      inline size_type result_dimension() const {
        return _b.size();
      }
      
      inline bool invertible() const {
        return _A.invertible(); }
     private:
      Matrix _A;
      Vector _b;
    };
      
    template <typename R>
    typename AffineMap<R>::State
    AffineMap<R>::operator() (const State& s) const
    {
      const Matrix& A=this->A();
      if (LinearAlgebra::number_of_columns(A)!=s.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Point& s): the map does not have the same dimension of the point.");
      }
      const Vector& pv=s.position_vector();
      Vector v=prod(A,pv);
      v+= this->b();
      return State(v);
    }
    
    template <typename R>
    Geometry::Simplex<R>
    AffineMap<R>::operator() (const Geometry::Simplex<R>& s) const
    {
      const Matrix& A=this->A();
      if (LinearAlgebra::number_of_columns(A)!=s.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Simplex<R>& s): the map does not have the same dimension of the simplex.");
      }
      
      const array<State>&  v=s.vertices();
      array<State> new_v(s.dimension());
      
      for(size_t i=0; i<s.dimension(); i++) {
        new_v[i]= (*this)(v[i]);
      }
     
      return Geometry::Simplex<R>(new_v);
    }

    template <typename R>
    Geometry::Rectangle<R>
    AffineMap<R>::operator() (const Geometry::Rectangle<R>& r) const
    {
      const Matrix& A=this->A();
      const Vector& b=this->b();
      if (LinearAlgebra::number_of_columns(A)!=r.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Rectangle<R>& r): the map does not have the same dimension of the rectangle.");
      }
      
      Base::array<Interval<R> > imv(r.dimension());
      
      for(size_t j=0; j<LinearAlgebra::number_of_rows(A); j++) {
        imv[j]=b[j];
        for(size_t i=0; i!=r.dimension(); ++i) {
           imv[j]+=A(j,i)*r[i];
        }
      }
      
      return Geometry::Rectangle<R>(imv);
    }
     
    template <typename R>
    Geometry::Parallelotope<R>
    AffineMap<R>::operator() (const Geometry::Parallelotope<R>& p) const
    {
      const Matrix& A=this->A();
      if (LinearAlgebra::number_of_columns(A)!=p.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Parallelotope<R>& p): the map does not have the same dimension of the parallelotope.");
      }
      State new_centre=(*this)(p.centre());
      return Geometry::Parallelotope<R>(new_centre,A*p.generators());
    }

    template <typename R>
    Geometry::Zonotope<R>
    AffineMap<R>::operator() (const Geometry::Zonotope<R>& z) const
    {
      const Matrix& A=this->A();
      if (LinearAlgebra::number_of_columns(A)!=z.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Zonotope<R>& z): the map does not have the same dimension of the zonotope."); 
      }
      
      State new_centre=(*this)(z.centre());
      return Geometry::Zonotope<R>(new_centre,A*z.principle_directions());
    }    
     
    template <typename R>
    Geometry::Polyhedron<R>
    AffineMap<R>::operator() (const Geometry::Polyhedron<R>& p) const
    {
      const Matrix& A=this->A();
      if (LinearAlgebra::number_of_columns(A)!=p.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Polyhedron<R>& p): the map does not have the same dimension of the polyhedron.");
      }
      throw std::domain_error("AffineMap<R>::operator() (const Geometry::Polyhedron<R>& p) not implemented.");
    }    

  }
}


#endif /* _ARIADNE_AFFINE_MAP_H */
