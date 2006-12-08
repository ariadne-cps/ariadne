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
 

#include "affine_map.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"




namespace Ariadne {
  namespace System {

    template<class R>
    Geometry::Point< typename AffineMap<R>::F >
    AffineMap<R>::image(const Geometry::Point<R>& pt) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;

      check_argument_dimension(*this,pt,__PRETTY_FUNCTION__);
      LinearAlgebra::Vector<F> image(this->A()*LinearAlgebra::Vector<F>(pt.position_vector())+this->b());
      return Geometry::Point<F>(image);
    }
    

    template<class R>
    Geometry::Rectangle<R>
    AffineMap<R>::image(const Geometry::Rectangle<R>& r) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;

      check_argument_dimension(*this,r,__PRETTY_FUNCTION__);
      return Geometry::Rectangle<R>(this->A()*r.position_vectors()+this->b());
    }
    
    template<class R>
    Geometry::Parallelotope<R>
    AffineMap<R>::image(const Geometry::Parallelotope<R>& p) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;

      check_argument_dimension(*this,p,__PRETTY_FUNCTION__);
      Geometry::Point<F> nc=this->image(p.centre());
      LinearAlgebra::Matrix<F> ng=this->A()*LinearAlgebra::Matrix<F>(p.generators());
      return Geometry::over_approximation(Geometry::Parallelotope<F>(nc,ng));
    }

    template<class R>
    Geometry::Zonotope<R>
    AffineMap<R>::image(const Geometry::Zonotope<R>& z) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;

      check_argument_dimension(*this,z,__PRETTY_FUNCTION__);
      return Geometry::over_approximation(Geometry::Zonotope< Interval<R> >(
        this->image(z.centre()),
        this->A()*LinearAlgebra::Matrix< Interval<R> >(z.generators())
      ));
    }   
    
    template<class R>
    Geometry::Polytope<R>
    AffineMap<R>::image(const Geometry::Polytope<R>& p) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }   
    
    template<class R>
    std::ostream& 
    AffineMap<R>::write(std::ostream& os) const
    {
      return os << "AffineMap( A=" << this->A()
                << ", b=" << this->b() << " )";
    }
    
/*
    template<class R>
    Geometry::ListSet<R,Geometry::Parallelotope>
    AffineMap<R>::operator() (const Geometry::GridMaskSet<R>& gms) const
    {
      Geometry::ListSet<R,Geometry::Parallelotope> result(gms.dimension());

      for(typename Geometry::GridMaskSet<R>::const_iterator iter=gms.begin();
          iter!=gms.end(); ++iter) {
        result.push_back(this->operator()(Geometry::Parallelotope<R>(*iter)));
      }
      
      return result;
    }  

    Geometry::Point<Rational>
    AffineMap<Rational>::operator() (const Geometry::Point<Rational>& pt) const
    {
      if (this->argument_dimension()!=pt.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Point&): "
          "the map does not have the same dimension of the point.");
      }
      return Geometry::Point<Rational>(this->A()*pt.position_vector()+this->b());
    }
    
    Geometry::Rectangle<Rational>
    AffineMap<Rational>::operator() (const Geometry::Rectangle<Rational>& r) const
    {
      if (this->argument_dimension()!=r.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Rectangle&): "
          "the map does not have the same dimension of the rectangle.");
      }
      return Geometry::Rectangle<Rational>(this->A()*r.position_vectors()+this->b());
    }
    
    Geometry::Parallelotope<Rational>
    AffineMap<Rational>::operator() (const Geometry::Parallelotope<Rational>& p) const
    {
      if (this->argument_dimension()!=p.dimension()) {
        throw std::domain_error("AffineMap<R>::operator() (const Parallelotope&): "
          "the map does not have the same dimension of the rectangle.");
      }
      return Geometry::Parallelotope<Rational>(
        this->operator() (p.centre()),A()*p.generators());
    }
    
    Geometry::Zonotope<Rational>
    AffineMap<Rational>::operator() (const Geometry::Zonotope<Rational>& z) const
    {
      if (this->argument_dimension()!=z.dimension()) {
        throw std::domain_error("AffineMap<Rational>::operator() (const Zonotope<Rational>&): "
          "the map does not have the same dimension of the simplex.");
      }

      return Geometry::Zonotope<Rational>(
          this->operator() (z.central_block()),this->A()*z.generators());
    }

    Geometry::Simplex<Rational>
    AffineMap<Rational>::operator() (const Geometry::Simplex<Rational>& s) const
    {
      if (this->argument_dimension()!=s.dimension()) {
        throw std::domain_error("AffineMap<Rational>::operator() (const Simplex<Rational>&): "
          "the map does not have the same dimension of the simplex.");
      }
      
      array< Geometry::Point<Rational> > new_vertices=s.vertices();
      for (size_type i=0; i< new_vertices.size(); ++i) {
        new_vertices[i]=this->operator() (new_vertices[i]);
      }
      return Geometry::Simplex<Rational>(new_vertices);
    }

    
    Geometry::Polytope<Rational>
    AffineMap<Rational>::operator() (const Geometry::Polytope<Rational>& p) const
    {
      if (this->argument_dimension()!=p.dimension()) {
        throw std::domain_error("AffineMap<Rational>::operator() (const Polytope<Rational>&): "
          "the map does not have the same dimension of the polytope.");
      }

      std::vector< Geometry::Point<Rational> > new_vertices=p.vertices();
      for (size_type i=0; i< new_vertices.size(); ++i) {
        new_vertices[i]=this->operator() (new_vertices[i]);
      }
      return Geometry::Polytope<Rational>(new_vertices);
    }    
*/

  }
}
