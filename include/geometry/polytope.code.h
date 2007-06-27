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

#include "../linear_programming/linear_program.h"

#include "../geometry/ddconv.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/polyhedron.h"

#include "../output/logging.h"

namespace Ariadne {
  namespace Geometry {

    

    template<class R>
    Polytope<R>::Polytope(dimension_type d)
      : _dimension(d), _number_of_vertices(0), _data()
    {
    }
   
    template<class R>
    Polytope<R>::Polytope(dimension_type d, size_type nv, const R* data)
      : _dimension(d), _number_of_vertices(nv), _data(data,data+(d+1)*nv)
    {
    }
   
   
    template<class R>
    Polytope<R>::Polytope(const LinearAlgebra::Matrix<R>& G)
      : _dimension(G.number_of_rows()-1), 
        _number_of_vertices(G.number_of_columns()), 
        _data(G.number_of_rows()*G.number_of_columns())
    {
      dimension_type& d=this->_dimension;
      size_type& nv=this->_number_of_vertices;
      R* ptr=this->_data.begin();
      LinearAlgebra::MatrixSlice<R>(d+1,nv,ptr,1u,d+1)=G;
    }

   
    template<class R>
    Polytope<R>::Polytope(const PointList<R>& pts)
      : _dimension(pts.dimension()),
        _number_of_vertices(pts.size()),
        _data((pts.dimension()+1)*pts.size())
    {
      dimension_type& d=this->_dimension;
      size_type& nv=this->_number_of_vertices;
      R* ptr=this->_data.begin();
      LinearAlgebra::MatrixSlice<R> g(d+1,nv,ptr,1,d+1);
      for(size_type j=0; j!=nv; ++j) {
        for(size_type i=0; i!=d; ++i) {
          g(i,j)=pts[j][i];
        }
        g(d,j)=1;
      }
    }
   
    template<class R>
    Polytope<R>::Polytope(const Rectangle<R>& r)
      : _dimension(r.dimension()),
        _number_of_vertices(r.number_of_vertices()),
        _data((r.dimension()+1)*r.number_of_vertices())
    {
      dimension_type& d=this->_dimension;
      size_type& nv=this->_number_of_vertices;
      R* ptr=this->_data.begin();
      LinearAlgebra::MatrixSlice<R> g(d+1,nv,ptr,1,d+1);

      size_type j=0;
      for(typename Rectangle<R>::vertices_const_iterator v=r.vertices_begin();
          v!=r.vertices_end(); ++v)
      {
        for(size_type i=0; i!=this->dimension(); ++i) {
          this->_generators_()(i,j)=(*v)[i];
        }
        this->_generators_()(d,j)=1;
      }
    }
   
    
    template<class R>
    Polytope<R>
    polytope(const Rectangle<R>& r)
    {
      return Polytope<R>(r);
    }
    
    
    template<class R>
    Polytope<typename Numeric::traits<R>::arithmetic_type>
    polytope(const Polyhedron<R>& plhd)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      
      dimension_type d=plhd.dimension();
      size_type nc=plhd.number_of_constraints();
      
      const LinearAlgebra::Matrix<R> A=plhd.A();
      const LinearAlgebra::Vector<R> b=plhd.b();
      
      std::vector< LinearAlgebra::Vector<F> > constraints;
      std::vector< LinearAlgebra::Vector<F> > generators;
      
      LinearAlgebra::Vector<F> tmp(d+1);
      // Add positivity constraint
      tmp(d)=1;
      constraints.push_back(tmp);
      for(size_type i=0; i!=nc; ++i) {
        for(size_type j=0; j!=d; ++j) {
          tmp(j)=-A(i,j);
        }
        tmp(d)=b(i);
        constraints.push_back(tmp);
      }

      ddconv(generators,constraints);
      
      size_type nv=generators.size();
      for(size_type j=0; j!=nv; ++j) {
        if(generators[j](d)==0) {
          //std::cerr << "Warning: Unbounded polytope\n";
          ARIADNE_THROW(UnboundedSet,"Polytope::Polytope(Polyhedron plhd)","plhd="<<plhd<<", generators="<<generators);
        }
      }
      array<F> data(nv*(d+1));
      for(size_type j=0; j!=nv; ++j) {
        for(size_type i=0; i!=d+1u; ++i) {
          data[j*(d+1)+i]=generators[j](i);
        }
      }
      Polytope<F> pltp(d,nv,data.begin());
      return pltp;
    }
    
    template<class R>
    Polytope<R>&
    Polytope<R>::operator=(const Polytope<R>& p)
      
    {
      if(this!=&p) { 
        this->_dimension=p._dimension;
        this->_number_of_vertices=p._number_of_vertices;
        this->_data=p._data;
      }
      return *this;
    }
   
    
    template<class R>
    dimension_type 
    Polytope<R>::dimension() const
    {
      return this->_dimension;
    }
    
    template<class R>
    const LinearAlgebra::Matrix<R>
    Polytope<R>::generators() const
    {
      return const_cast<Polytope<R>*>(this)->_generators_();
    }
    
    
    template<class R>
    size_type 
    Polytope<R>::number_of_vertices() const
    {
      return this->_number_of_vertices;
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
      LinearAlgebra::Matrix<R> g=this->generators();
      Point<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=g(i,j);
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


    template<class R> inline
    Rectangle<R> 
    bounding_box(const Polytope<R>& p)  
    {
      //std::cerr << "Polytope<R>::bounding_box()" << std::endl;
      typename Polytope<R>::vertices_const_iterator pt_iter=p.vertices_begin();
      Rectangle<R> result(*pt_iter);
      ++pt_iter;
      for( ; pt_iter!=p.vertices_end(); ++pt_iter) {
        result=rectangular_hull(result,Rectangle<R>(*pt_iter));
      }
      return result;
    }
      
    template<class R> inline
    Rectangle< Numeric::Interval<R> > 
    bounding_box(const Polytope< Numeric::Interval<R> >& p)  
    {
      typename Polytope< Numeric::Interval<R> >::vertices_const_iterator pt_iter=p.vertices_begin();
      Rectangle<R> result(*pt_iter);
      ++pt_iter;
      for( ; pt_iter!=p.vertices_end(); ++pt_iter) {
        result=rectangular_hull(result,Rectangle<R>(*pt_iter));
      }
      return result;
    }
      
    template<class R>
    Rectangle<R> 
    Polytope<R>::bounding_box() const 
    {
      return Geometry::bounding_box(*this);
    }


    template<class R>
    tribool 
    Polytope<R>::empty() const
    {
      return this->number_of_vertices()==0;
    }
    
    template<class R>
    tribool 
    Polytope<R>::bounded() const
    {
      tribool result=true;
      dimension_type d=this->dimension();
      size_type nv=this->number_of_vertices();
      R zero=0;
      for(size_type i=0; i!=nv; ++i) {
        result = result && (this->_data[i*(d+1u)+d]!=zero);
      }
      return result;
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
    contains(const Polytope< Numeric::Interval<R> >& ply, const Point< Numeric::Interval<R> >& pt)
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
      typedef Numeric::Rational F;
      return subset(Polytope<F>(A),Polyhedron<F>(Polytope<F>(B)));
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
      return Geometry::subset(Rectangle<Numeric::Rational>(A),
                              Polyhedron<Numeric::Rational>(Polytope<Numeric::Rational>(B)));
    }
    
    
    template<class R>
    tribool 
    disjoint(const Polytope< Numeric::Interval<R> >& ply, const Rectangle< Numeric::Interval<R> >& rect) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    disjoint(const Polytope< Numeric::Interval<R> >& ply1, const Polytope< Numeric::Interval<R> >& ply2) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    subset(const Rectangle< Numeric::Interval<R> >& ply1, const Polytope< Numeric::Interval<R> >& ply2) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    subset(const Polytope< Numeric::Interval<R> >& ply1, const Polytope< Numeric::Interval<R> >& ply2) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Polytope<R>
    convex_hull(const Polytope<R>& A, const Polytope<R>& B)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    

    template<class R>
    std::string
    Polytope<R>::name()
    {
      return std::string("Polytope")+"<"+Numeric::name<R>()+">";
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
    void
    Polytope<R>::_instantiate_geometry_operators()
    {   
      Rectangle<R> r;
      Polytope<R> p;
      Polyhedron<R> h;
      
      Geometry::disjoint(r,p);
      Geometry::disjoint(p,r);
      Geometry::disjoint(p,p);
      Geometry::subset(r,p);
      Geometry::subset(p,r);
      Geometry::subset(p,p);
      Geometry::convex_hull(p,p);
      
      polytope(r);
      polytope(h);
      polyhedron(p);
      
    }

    
  /*
    
    template<class R>
    dimension_type 
    Polytope< Numeric::Interval<R> >::dimension() const 
    { 
      return _vertices.dimension(); 
    }
    
    
    template<class R>
    size_type 
    Polytope< Numeric::Interval<R> >::number_of_vertices() const 
    { 
      return _vertices.size(); 
    }
    
    
    template<class R>
    const PointList< Numeric::Interval<R> >&
    Polytope< Numeric::Interval<R> >::vertices() const 
    {
      return this->_vertices;
    }
    
  */

  }
}
