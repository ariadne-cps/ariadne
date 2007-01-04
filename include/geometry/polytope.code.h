/***************************************************************************
 *            polytope.code.h
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

#include "polytope.h"

#include "../base/stlio.h"

#include "../numeric/interval.h"
#include "../numeric/arithmetic.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../geometry/ddconv.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {

    template<class R>
    void
    Polytope<R>::_instantiate_geometry_operators()
    {
      Rectangle<R> r;
      Polytope<R> p;
      
      Geometry::disjoint(r,p);
      Geometry::disjoint(p,r);
      Geometry::disjoint(p,p);
      Geometry::subset(r,p);
      Geometry::subset(p,r);
      Geometry::subset(p,p);
      Geometry::convex_hull(p,p);
    }


    template<class R>
    Polytope<R>::Polytope(dimension_type d)
      : _generators(d+1,0)
    {
    }
   
    template<class R>
    Polytope<R>::Polytope(dimension_type d, size_type nv, const R* data)
      : _generators(d+1,nv,data,1,d+1)
    {
    }
   
   
    template<class R>
    Polytope<R>::Polytope(const LinearAlgebra::Matrix<R>& G)
      : _generators(G.number_of_rows()+1,G.number_of_columns())
    {
      for(size_type j=0; j!=G.number_of_columns(); ++j) {
        for(size_type i=0; i!=G.number_of_rows(); ++i) {
          this->_generators(i,j)=G(i,j);
        }
        this->_generators(this->dimension(),j)=static_cast<R>(1);
      }
    }

   
    template<class R>
    Polytope<R>::Polytope(const PointList<R>& pts)
      : _generators(pts.dimension()+1,pts.size())
    {
      for(size_type j=0; j!=pts.size(); ++j) {
        for(size_type i=0; i!=pts.dimension(); ++i) {
          this->_generators(i,j)=pts[j][i];
        }
        this->_generators(this->dimension(),j)=1;
      }
    }
   
    template<class R>
    Polytope<R>::Polytope(const Rectangle<R>& r)
      : _generators(r.dimension()+1,1<<r.dimension())
    {
      dimension_type d=r.dimension();
      size_type j=0;
      for(typename Rectangle<R>::vertices_const_iterator v=r.vertices_begin();
          v!=r.vertices_end(); ++v)
      {
        for(size_type i=0; i!=this->dimension(); ++i) {
          this->_generators(i,j)=(*v)[i];
        }
        this->_generators(d,j)=1;
      }
    }
   
    template<class R>
    Polytope<R>::Polytope(const Polyhedron<R>& pltp)
      : _generators()
    {   
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<>
    Polytope<Rational>::Polytope(const Polyhedron<Rational>& plhd)
      : _generators()
    {
      //std::cerr << "Polytope<Rational>::Polytope(const Polyhedron<Rational>&)" << std::endl;
      typedef Rational R;
      typedef Rational F;
      
      dimension_type d=plhd.dimension();
      size_type nc=plhd.number_of_constraints();
      
      const LinearAlgebra::Matrix<R>& A=plhd.A();
      const LinearAlgebra::Vector<R>& b=plhd.b();
      LinearAlgebra::Matrix<R>& G=this->_generators;
      
      std::vector< LinearAlgebra::Vector<F> > result;
      std::vector< LinearAlgebra::Vector<F> > argument;
      
      LinearAlgebra::Vector<F> tmp(d+1);
      
      for(size_type i=0; i!=nc; ++i) {
        for(size_type j=0; j!=d; ++j) {
          tmp(j)=-A(i,j);
        }
        tmp(d)=b(i);
        argument.push_back(tmp);
      }
      ddconv(result,argument);     
      
      //std::cerr << "argument=" << argument << std::endl;
      //std::cerr << "result=" << result << std::endl;
      
      size_type nv=result.size();
      G.resize(d+1,nv);
      for(size_type j=0; j!=nv; ++j) {
        F s=result[j](d);
        if(s!=F(0)) {
          for(size_type i=0; i!=d; ++i) {
            G(i,j)=result[j](i)/s;
          }
          G(d,j)=1;
        } else {
          throw std::runtime_error("Unbounded polytope");
        }
      }
      //std::cerr << "G=" << G << std::endl;
    }
    
    template<class R>
    Polytope<R>&
    Polytope<R>::operator=(const Polytope<R>& p)
      
    {
      if(this!=&p) { 
        this->_generators=p._generators;
      }
      return *this;
    }
   
    
    template<class R>
    dimension_type 
    Polytope<R>::dimension() const
    {
      return this->_generators.number_of_rows()-1;
    }
    
    template<class R>
    const LinearAlgebra::Matrix<R>&
    Polytope<R>::generators() const
    {
      return this->_generators;
    }
    
    
    template<class R>
    size_type 
    Polytope<R>::number_of_vertices() const
    {
      return this->_generators.number_of_columns();
    }
    
    template<class R>
    PointList<R> 
    Polytope<R>::vertices() const 
    {
      return PointList<R>(this->generators());
    }

    template<class R>
    Point<R>
    Polytope<R>::vertex(const size_type& j) const 
    {
      Point<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=this->_generators(i,j);
      }
      return result;
    }

    template<class R>
    typename Polytope<R>::vertices_const_iterator
    Polytope<R>::vertices_begin() const 
    {
      return PolytopeVerticesIterator<R>(*this,0);
    }

    template<class R>
    typename Polytope<R>::vertices_const_iterator
    Polytope<R>::vertices_end() const 
    {
      return PolytopeVerticesIterator<R>(*this,this->number_of_vertices());
    }

    template<class R>
    Rectangle<R> 
    Polytope<R>::bounding_box() const 
    {
      //std::cerr << "Polytope<R>::bounding_box()" << std::endl;
      vertices_const_iterator pt_iter=this->vertices_begin();
      Rectangle<R> result(*pt_iter);
      ++pt_iter;
      for( ; pt_iter!=this->vertices_end(); ++pt_iter) {
        result=rectangular_hull(result,Rectangle<R>(*pt_iter));
      }
      return result;
    }
      

    template<class R>
    tribool 
    Polytope<R>::empty() const
    {
      return indeterminate;
    }
    
    template<class R>
    tribool 
    Polytope<R>::bounded() const
    {
      return true;
    }
    
    /*!Set up linear programming problem
     * Try to simultaneously solve A*x=p where A is the extended vertex matrix
     * d+1 auxiliary variables, d+1 equations
     */
    template<class R>
    tribool 
    contains(const Polytope<R>& ply, const Point<R>& pt)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    Polytope<R>::contains(const Point<R>& pt) const
    {
      return Geometry::contains(*this,pt);
    }

    template<class R>
    tribool 
    disjoint(const Polytope<R>& ply, const Rectangle<R>& rect)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    disjoint(const Rectangle<R>& A, const Polytope<R>& B)
    { 
      return disjoint(B,A);
    }

    template<class R>
    tribool 
    disjoint(const Polytope<R>& A, const Polytope<R>& B)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
     
    template<class R>
    tribool 
    subset(const Polytope<R>& A, const Polytope<R>& B)
    {
      //typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef Rational F;
      return subset(Polytope<F>(A),Polyhedron<F>(B));
    }
      
    
    
    template<class R>
    tribool 
    subset(const Polytope<R>& A, const Rectangle<R>& B)
    {
      tribool result=true;
      for(typename Polytope<R>::vertices_const_iterator v=A.vertices_begin(); v!=A.vertices_end(); ++v) {
        result = result && B.contains(*v);
        if(!result) { return result; }
      }
      return result;
    }
    
    
    template<class R>
    tribool 
    subset(const Rectangle<R>& A, const Polytope<R>& B)
    {
      return Geometry::subset(Rectangle<Rational>(A),Polyhedron<Rational>(B));
    }
    
    template<class R>
    Polytope<R>
    convex_hull(const Polytope<R>& A, const Polytope<R>& B)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>  
    std::ostream& 
    Polytope<R>::write(std::ostream& os) const 
    {
      return os << "Polytope( vertices=" << this->vertices() << " )";
    }
    
    template<class R>  
    std::istream& 
    Polytope<R>::read(std::istream& is)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    
    
    
    template<class R>
    dimension_type 
    Polytope< Interval<R> >::dimension() const 
    { 
      return _vertices.dimension(); 
    }
    
    
    template<class R>
    size_type 
    Polytope< Interval<R> >::number_of_vertices() const 
    { 
      return _vertices.size(); 
    }
    
    
    template<class R>
    const PointList< Interval<R> >&
    Polytope< Interval<R> >::vertices() const 
    {
      return this->_vertices;
    }
    
    

  }
}
