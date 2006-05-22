/***************************************************************************
 *            affine_map.tpl
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 

#include "../evaluation/affine_map.h"


namespace Ariadne {
  namespace Evaluation {

  
    template <typename R>
    Geometry::Point<R>
    AffineMap<R>::operator() (const Geometry::Point<R>& s) const
    {
      const Matrix_type& A=this->A();
      if (A.size2()!=s.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Point& s): the map does not have the same dimension of the point.");
      }
      const Vector_type& pv=s.position_vector();
      Vector_type v=prod(A,pv);
      v+= this->b();
      return Geometry::Point<R>(v);
    }
    
    template <typename R>
    Geometry::Simplex<R>
    AffineMap<R>::operator() (const Geometry::Simplex<R>& s) const
    {
      const Matrix_type& A=this->A();
      if (A.size2()!=s.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Simplex<R>& s): the map does not have the same dimension of the simplex.");
      }
      
      const array<state_type>&  v=s.vertices();
      array<state_type> new_v(s.dimension());
      
      for(size_t i=0; i<s.dimension(); i++) {
        new_v[i]= (*this)(v[i]);
      }
     
      return Geometry::Simplex<R>(new_v);
    }

    template <typename R>
    Geometry::Rectangle<R>
    AffineMap<R>::operator() (const Geometry::Rectangle<R>& r) const
    {
      const Matrix_type& A=this->A();
      const Vector_type& b=this->b();
      if (A.size2()!=r.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Rectangle<R>& r): the map does not have the same dimension of the rectangle.");
      }
      
      Base::array<Interval<R> > imv(r.dimension());
      
      for(size_t j=0; j<A.size1(); j++) {
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
      const Matrix_type& A=this->A();
      if (A.size2()!=p.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Parallelotope<R>& p): the map does not have the same dimension of the parallelotope.");
      }
      state_type new_centre=(*this)(p.centre());
      return Geometry::Parallelotope<R>(new_centre,A*p.generators());
    }

    template <typename R>
    Geometry::Zonotope<R>
    AffineMap<R>::operator() (const Geometry::Zonotope<R>& z) const
    {
      const Matrix_type& A=this->A();
      if (A.size2()!=z.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Zonotope<R>& z): the map does not have the same dimension of the zonotope."); 
      }
      
      state_type new_centre=(*this)(z.centre());
      return Geometry::Zonotope<R>(new_centre,A*z.generators());
    }    
    
    /* TO IMPROVE */
    template <typename R>
    Geometry::Polyhedron<R>
    AffineMap<R>::operator() (const Geometry::Polyhedron<R>& p) const
    {
      const Matrix_type& A=this->A();
      if (A.size2()!=p.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Geometry::Polyhedron<R>& p): the map does not have the same dimension of the polyhedron.");
      }

      std::vector< Geometry::Point<R> > vert=p.vertices();
      LinearAlgebra::Vector<R> new_point_pos;

      for (size_t i=0; i< vert.size(); i++) {
       new_point_pos=this->b()+prod(A,vert[i].position_vector());
       vert[i]=Geometry::Point<R>(new_point_pos); 
      }

      return Geometry::Polyhedron<R>(vert);
    }    

    /* TO IMPROVE */
    template <typename R>
    Geometry::ListSet<R,Geometry::Parallelotope>
    AffineMap<R>::operator() (const Geometry::GridMaskSet<R>& cms) const
    {
      using namespace Ariadne::Geometry;

      ListSet<R,Rectangle> lrs=ListSet<R,Rectangle>(cms);

      ListSet<R,Geometry::Parallelotope> output(cms.dimension());

      for (size_t i=0; i< lrs.size(); i++) {
	Parallelotope<R> p=(*this)(Parallelotope<R>(lrs[i]));

	output.push_back(p);
      }
      
      return output;
    }  
  }
}
