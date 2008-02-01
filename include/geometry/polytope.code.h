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

#include "base/stlio.h"

#include "numeric/interval.h"
#include "numeric/arithmetic.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "linear_programming/linear_program.h"

#include "geometry/ddconv.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/box.h"
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"

#include "output/logging.h"



namespace {

using namespace Ariadne;

inline
tribool 
disjoint(const Geometry::Polytope<Numeric::Rational>& ply, const Geometry::Box<Numeric::Rational>& bx) 
{
  return Geometry::disjoint(Geometry::Polyhedron<Numeric::Rational>(ply),bx);
}

template<class T>
inline
tribool 
disjoint(const Geometry::Polytope< Numeric::Float<T> >& ply, 
         const Geometry::Box< Numeric::Float<T> >& bx) 
{
  return ::disjoint(Geometry::Polytope<Numeric::Rational>(ply),Geometry::Box<Numeric::Rational>(bx));
}

template<class R>
inline
tribool 
disjoint(const Geometry::Polytope< Numeric::Interval<R> >& ply, const Geometry::Box<R>& bx) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


inline
tribool 
disjoint(const Geometry::Polytope<Numeric::Rational>& ply1, const Geometry::Polytope<Numeric::Rational>& ply2) 
{
  return Geometry::disjoint(Geometry::Polyhedron<Numeric::Rational>(ply1),Geometry::Polyhedron<Numeric::Rational>(ply1));
}

template<class R>
inline
tribool 
disjoint(const Geometry::Polytope<R>& ply1, const Geometry::Polytope<R>& ply2) 
{
  return Geometry::disjoint(Geometry::Polytope<Numeric::Rational>(ply1),Geometry::Polytope<Numeric::Rational>(ply1));
}

template<class R>
tribool 
disjoint(const Geometry::Polytope< Numeric::Interval<R> >& ply1, const Geometry::Polytope< Numeric::Interval<R> >& ply2) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




inline
tribool 
subset(const Geometry::Box<Numeric::Rational>& bx, const Geometry::Polytope<Numeric::Rational>& ply) 
{
  return Geometry::subset(bx,Geometry::Polyhedron<Numeric::Rational>(ply));
}

template<class R>
inline
tribool 
subset(const Geometry::Box<R>& bx, const Geometry::Polytope<R>& ply) 
{
  return Geometry::subset(Geometry::Box<Numeric::Rational>(bx),Geometry::Polytope<Numeric::Rational>(ply));
}

template<class R>
inline
tribool 
subset(const Geometry::Box<R>& bx, const Geometry::Polytope< Numeric::Interval<R> >& ply) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



inline
tribool 
subset(const Geometry::Polytope<Numeric::Rational>& ply1, const Geometry::Polytope<Numeric::Rational>& ply2) 
{
  return Geometry::subset(ply1,Geometry::Polyhedron<Numeric::Rational>(ply2));
}

template<class R>
inline
tribool 
subset(const Geometry::Polytope<R>& ply1, const Geometry::Polytope<R>& ply2) 
{
  return Geometry::subset(Geometry::Polytope<Numeric::Rational>(ply1),Geometry::Polytope<Numeric::Rational>(ply2));
}

template<class R>
inline
tribool 
subset(const Geometry::Polytope< Numeric::Interval<R> >& ply1, const Geometry::Polytope< Numeric::Interval<R> >& ply2) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




} // namespace 





namespace Ariadne {


    

template<class X>
Geometry::Polytope<X>::Polytope(dimension_type d)
  : _dimension(d), _number_of_vertices(0), _data()
{
}

template<class X>
Geometry::Polytope<X>::Polytope(dimension_type d, size_type nv, const X* data)
  : _dimension(d), _number_of_vertices(nv), _data(data,data+(d+1)*nv)
{
}


template<class X>
Geometry::Polytope<X>::Polytope(const LinearAlgebra::Matrix<X>& G)
  : _dimension(G.number_of_rows()-1), 
    _number_of_vertices(G.number_of_columns()), 
    _data(G.number_of_rows()*G.number_of_columns())
{
  dimension_type& d=this->_dimension;
  size_type& nv=this->_number_of_vertices;
  X* ptr=this->_data.begin();
  LinearAlgebra::MatrixSlice<X>(d+1,nv,ptr,1u,d+1)=G;
}


template<class X>
Geometry::Polytope<X>::Polytope(const PointList<X>& pts)
  : _dimension(pts.dimension()),
    _number_of_vertices(pts.size()),
    _data((pts.dimension()+1)*pts.size())
{
  dimension_type& d=this->_dimension;
  size_type& nv=this->_number_of_vertices;
  X* ptr=this->_data.begin();
  LinearAlgebra::MatrixSlice<X> g(d+1,nv,ptr,1,d+1);
  for(size_type j=0; j!=nv; ++j) {
    for(size_type i=0; i!=d; ++i) {
      g(i,j)=pts[j][i];
    }
    g(d,j)=1;
  }
}

template<class X> 
Geometry::Polytope<X>::Polytope(const Box<R>& bx)
  : _dimension(bx.dimension()),
    _number_of_vertices(bx.number_of_vertices()),
    _data((bx.dimension()+1)*bx.number_of_vertices())
{
  dimension_type& d=this->_dimension;
  size_type& nv=this->_number_of_vertices;
  X* ptr=this->_data.begin();
  LinearAlgebra::MatrixSlice<X> g(d+1,nv,ptr,1,d+1);
  
  size_type j=0;
  for(typename Box<R>::vertices_const_iterator v=bx.vertices_begin();
      v!=bx.vertices_end(); ++v)
  {
    for(size_type i=0; i!=this->dimension(); ++i) {
      this->_generators_()(i,j)=(*v)[i];
    }
    this->_generators_()(d,j)=1;
  }
}


template<class X>
Geometry::Polytope<X>
Geometry::polytope(const Rectangle<X>& r)
{
  dimension_type d=r.dimension();
  size_type nv=r.number_of_vertices();
  array<X> data((d+1)*nv);
  size_type j=0;
  for(typename Rectangle<X>::vertices_const_iterator v=r.vertices_begin();
      v!=r.vertices_end(); ++v)
  {
    for(size_type i=0; i!=d; ++i) {
      data[j]=(*v)[i];
      ++j;
    }
    data[d]=1;
    ++j;
  }
  return Polytope<X>(d,nv,data.begin());
}  


template<class X>
Geometry::Polytope<typename Numeric::traits<X>::arithmetic_type>
Geometry::polytope(const Polyhedron<X>& plhd)
{
  typedef typename Numeric::traits<X>::arithmetic_type F;
  
  dimension_type d=plhd.dimension();
  size_type nc=plhd.number_of_constraints();
  
  const LinearAlgebra::Matrix<X> A=plhd.A();
  const LinearAlgebra::Vector<X> b=plhd.b();
  
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


template<class X>
dimension_type 
Geometry::Polytope<X>::dimension() const
{
  return this->_dimension;
}

template<class X>
const LinearAlgebra::Matrix<X>
Geometry::Polytope<X>::generators() const
{
  return const_cast<Polytope<X>*>(this)->_generators_();
}


template<class X>
size_type 
Geometry::Polytope<X>::number_of_vertices() const
{
  return this->_number_of_vertices;
}

template<class X>
Geometry::PointList<X> 
Geometry::Polytope<X>::vertices() const 
{
  return PointList<X>(this->generators());
}

template<class X>
Geometry::Point<X>
Geometry::Polytope<X>::vertex(const size_type& j) const 
{
  LinearAlgebra::Matrix<X> g=this->generators();
  Point<X> result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result[i]=g(i,j);
  }
  return result;
}

template<class X>
typename Geometry::Polytope<X>::vertices_const_iterator
Geometry::Polytope<X>::vertices_begin() const 
{
  return PolytopeVerticesIterator<X>(*this,0);
}

template<class X>
typename Geometry::Polytope<X>::vertices_const_iterator
Geometry::Polytope<X>::vertices_end() const 
{
  return PolytopeVerticesIterator<X>(*this,this->number_of_vertices());
}


template<class X> 
Geometry::Box<typename Geometry::Polytope<X>::R> 
Geometry::Polytope<X>::bounding_box() const
{
  //std::cerr << "Polytope<X>::bounding_box()" << std::endl;
  typename Polytope<X>::vertices_const_iterator pt_iter=this->vertices_begin();
  Box<R> result(*pt_iter);
  ++pt_iter;
  for( ; pt_iter!=this->vertices_end(); ++pt_iter) {
    result=rectangular_hull(result,Box<R>(*pt_iter));
  }
  return result;
}

template<class X>
tribool 
Geometry::Polytope<X>::empty() const
{
  return this->number_of_vertices()==0;
}

template<class X>
tribool 
Geometry::Polytope<X>::bounded() const
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
template<class X>
tribool 
Geometry::contains(const Polytope<X>& ply, const Point<X>& pt)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>
tribool 
Geometry::contains(const Polytope< Numeric::Interval<X> >& ply, const Point< Numeric::Interval<X> >& pt)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
tribool 
Geometry::Polytope<X>::contains(const Point<X>& pt) const
{
  return Geometry::contains(*this,pt);
}

template<class X, class R>
tribool 
Geometry::disjoint(const Polytope<X>& ply, const Box<R>& bx)
{
  return ::disjoint(ply,bx);
}


template<class X, class R>
tribool 
Geometry::disjoint(const Box<R>& bx, const Polytope<X>& ply)
{ 
  return ::disjoint(ply,bx);
}

template<class X>
tribool 
Geometry::disjoint(const Polytope<X>& ply1, const Polytope<X>& ply2)
{
  return ::disjoint(ply1,ply2);
}

template<class X>
tribool 
Geometry::subset(const Polytope<X>& ply1, const Polytope<X>& ply2)
{
  return ::subset(ply1,ply2);
}

template<class X, class R>
tribool 
Geometry::subset(const Box<R>& bx, const Polytope<X>& ply)
{
  return ::subset(bx,ply);
}

template<class X, class R>
tribool 
Geometry::subset(const Polytope<X>& ply, const Box<R>& bx)
{
  tribool result=true;
  for(typename Polytope<X>::vertices_const_iterator v=ply.vertices_begin(); v!=ply.vertices_end(); ++v) {
    result = result && bx.contains(*v);
    if(!result) { return result; }
  }
  return result;
}





template<class X>
Geometry::Polytope<X>
Geometry::convex_hull(const Polytope<X>& A, const Polytope<X>& B)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
std::string
Geometry::Polytope<X>::name()
{
  return std::string("Polytope")+"<"+Numeric::name<X>()+">";
}


template<class X>  
std::ostream& 
Geometry::Polytope<X>::write(std::ostream& os) const 
{
  return os << "Polytope( vertices=" << this->vertices() << " )";
}

template<class X>  
std::istream& 
Geometry::Polytope<X>::read(std::istream& is)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
void
Geometry::Polytope<X>::_instantiate()
{   
  typedef typename Numeric::traits<X>::number_type R;
  typedef typename Numeric::traits<X>::arithmetic_type I;

  tribool tb;
  Box<R>* bx=0;
  Rectangle<X>* r=0;
  Polytope<X>* p=0;
  Polyhedron<X>* h=0;
  //Rectangle<I>* ir=0;
  Polytope<I>* ip=0;
  Polyhedron<I>* ih=0;

  tb=Geometry::disjoint(*bx,*p);
  Geometry::disjoint(*p,*bx);
  Geometry::disjoint(*p,*p);
  Geometry::subset(*bx,*p);
  Geometry::subset(*p,*bx);
  Geometry::subset(*p,*p);
  Geometry::convex_hull(*p,*p);
  
  *p=polytope(*r);
  *ip=polytope(*h);
  *ih=polyhedron(*p);
  
}



} // namespace Ariadne
