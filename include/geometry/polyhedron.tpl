/***************************************************************************
 *            polyhedron.tpl
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include <ppl.hh>

#include "polyhedron.h"

#include "../utility/stlio.h"
#include "../numeric/interval.h"

#include "../numeric/arithmetic.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/ppl_polyhedron.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {

    using Parma_Polyhedra_Library::C_Polyhedron;
    
    template <typename R> inline 
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const Point<R>& pt) 
    {
      return ppl_polyhedron(LinearAlgebra::Vector<Rational>(pt.position_vector()));
    }
    


    template <typename R>
    Polyhedron<R>::Polyhedron(const LinearAlgebra::Matrix<R>& A)
      : _vertices(A)
    {
    }
   
    template <typename R>
    Polyhedron<R>::Polyhedron(const PointList<R>& pts)
      : _vertices(pts)
    {
    }
   
    template <typename R>
    Polyhedron<R>::Polyhedron(const Rectangle<R>& r)
      : _vertices(r.vertices())
    {
    }
   
    template <typename R>
    Polyhedron<R>::Polyhedron(const Polyhedron<R>& p)
      : _vertices(p._vertices)
    {
    }
   
    
    template <typename R>
    Polyhedron<R>::operator Parma_Polyhedra_Library::C_Polyhedron() const
    {
      return ppl_polyhedron(this->_vertices);
    }
   
    template <typename R>
    dimension_type 
    Polyhedron<R>::dimension() const
    {
      return this->_vertices.size(1);
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::empty() const
    {
      return Geometry::empty(ppl_polyhedron(this->_vertices));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::empty_interior() const
    {
      return Geometry::empty_interior(ppl_polyhedron(this->_vertices));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::contains(const Point<R>& pt) const
    {
      return Geometry::subset(ppl_polyhedron(pt),ppl_polyhedron(this->_vertices));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::interior_contains(const Point<R>& pt) const
    {
      return Geometry::inner_subset(ppl_polyhedron(pt),ppl_polyhedron(this->_vertices));
    }

    template <typename R>
    PointList<R> 
    Polyhedron<R>::vertices() const 
    {
      return PointList<R>(this->_vertices);
    }

    template <typename R>
    Rectangle<R> 
    Polyhedron<R>::bounding_box() const 
    {
      Rectangle<R> result(this->dimension());
      PointList<R> v=this->_vertices;
      for(typename PointList<R>::const_iterator pt_iter=v.begin(); pt_iter!=v.end(); ++pt_iter) {
        result=rectangular_hull(result,Rectangle<R>(*pt_iter));
      }
      return result;
    }
      

    template <typename R>
    bool 
    Polyhedron<R>::equal(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::equal(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::disjoint(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::disjoint(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::interiors_intersect(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::interiors_intersect(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::subset(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::subset(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polyhedron<R>::inner_subset(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::inner_subset(C_Polyhedron(A),C_Polyhedron(B));
    }
          
    template <typename R>
    Polyhedron<R>
    Polyhedron<R>::convex_hull(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      throw std::runtime_error("Polyhedron<R>::convex_hull(const Polyhedron& A, const Polyhedron& B) not implemented");
    }
    
    template<typename R>  
    std::ostream& operator<<(std::ostream& os, const Polyhedron<R>& p) {
      return os << "Polyhedron( vertices=" << p.vertices() << " )";
    }
    
    template<typename R>  
    std::istream& operator>>(std::istream& is, Polyhedron<R>& p) {
      throw std::runtime_error("std::istream& operator>>(std::istream&, Polyhedron<R>&) not implemented");
    }

    
    
    
    
    template <typename R>
    Polytope<R>::operator Parma_Polyhedra_Library::C_Polyhedron() const
    {
      return ppl_polyhedron(this->_A,this->_b);
    }
    
    template <typename R>
    PointList<Rational>
    Polytope<R>::vertices() const
    {
      return PointList<Rational>(Geometry::generators(ppl_polyhedron(this->_A,this->_b)));
    }

    template <typename R>
    bool 
    Polytope<R>::contains(const Point<R>& pt) const
    {
      return Geometry::subset(ppl_polyhedron(pt),
                              ppl_polyhedron(this->_A,this->_b));
    }
    
    template <typename R>
    bool 
    Polytope<R>::interior_contains(const Point<R>& pt) const
    {
      return Geometry::inner_subset(ppl_polyhedron(pt),
                                    ppl_polyhedron(this->_A,this->_b));
    }
    
    template <typename R>
    bool 
    Polytope<R>::equal(const Polytope<R>& A, const Polytope<R>& B)
    {
      return Geometry::equal(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polytope<R>::disjoint(const Polytope<R>& A, const Polytope<R>& B)
    {
      return Geometry::disjoint(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polytope<R>::interiors_intersect(const Polytope<R>& A, const Polytope<R>& B)
    {
      return Geometry::interiors_intersect(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polytope<R>::subset(const Polytope<R>& A, const Polytope<R>& B)
    {
      return Geometry::subset(C_Polyhedron(A),C_Polyhedron(B));
    }
    
    template <typename R>
    bool 
    Polytope<R>::inner_subset(const Polytope<R>& A, const Polytope<R>& B)
    {
      return Geometry::inner_subset(C_Polyhedron(A),C_Polyhedron(B));
    }
          
    template<typename R> inline 
    std::ostream& operator<<(std::ostream& os, const Polytope<R>& p) {
      return os << "Polytope( A=" << p.A() << ", b=" << p.b() << " )";
    }
    
    template<typename R> inline 
    std::istream& operator>>(std::istream& is, Polytope<R>& p) {
      throw std::runtime_error("std::istream& operator>>(std::istream&, Polytope<R>&) not implemented");
    }


  }
}
