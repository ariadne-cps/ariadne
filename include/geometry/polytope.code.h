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

#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "linear_programming/linear_program.h"

#include "geometry/ddconv.h"
#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/box.h"
#include "geometry/polyhedron.h"

#include "output/logging.h"



namespace {

using namespace Ariadne;

inline
tribool 
contains(const Polytope<Rational>& ply, const Point<Rational>& pt) 
{
  return contains(Polyhedron<Rational>(ply),pt);
}

template<class T>
inline
tribool 
contains(const Polytope< Float<T> >& ply, 
         const Point< Float<T> >& pt) 
{
  return ::contains(Polytope<Rational>(ply),Point<Rational>(pt));
}

template<class R>
inline
tribool 
contains(const Polytope< Interval<R> >& ply, const Point<R>& pt) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



inline
tribool 
disjoint(const Polytope<Rational>& ply, const Box<Rational>& bx) 
{
  return disjoint(Polyhedron<Rational>(ply),bx);
}

template<class T>
inline
tribool 
disjoint(const Polytope< Float<T> >& ply, 
         const Box< Float<T> >& bx) 
{
  return ::disjoint(Polytope<Rational>(ply),Box<Rational>(bx));
}

template<class R>
inline
tribool 
disjoint(const Polytope< Interval<R> >& ply, const Box<R>& bx) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


inline
tribool 
disjoint(const Polytope<Rational>& ply1, const Polytope<Rational>& ply2) 
{
  return disjoint(Polyhedron<Rational>(ply1),Polyhedron<Rational>(ply1));
}

template<class R>
inline
tribool 
disjoint(const Polytope<R>& ply1, const Polytope<R>& ply2) 
{
  return disjoint(Polytope<Rational>(ply1),Polytope<Rational>(ply1));
}

template<class R>
tribool 
disjoint(const Polytope< Interval<R> >& ply1, const Polytope< Interval<R> >& ply2) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




inline
tribool 
subset(const Box<Rational>& bx, const Polytope<Rational>& ply) 
{
  return subset(bx,Polyhedron<Rational>(ply));
}

template<class R>
inline
tribool 
subset(const Box<R>& bx, const Polytope<R>& ply) 
{
  return subset(Box<Rational>(bx),Polytope<Rational>(ply));
}

template<class R>
inline
tribool 
subset(const Box<R>& bx, const Polytope< Interval<R> >& ply) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



inline
tribool 
subset(const Polytope<Rational>& ply1, const Polytope<Rational>& ply2) 
{
  return subset(ply1,Polyhedron<Rational>(ply2));
}

template<class R>
inline
tribool 
subset(const Polytope<R>& ply1, const Polytope<R>& ply2) 
{
  return subset(Polytope<Rational>(ply1),Polytope<Rational>(ply2));
}

template<class R>
inline
tribool 
subset(const Polytope< Interval<R> >& ply1, const Polytope< Interval<R> >& ply2) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}




} // namespace 





namespace Ariadne {


    

template<class X>
Polytope<X>::Polytope(dimension_type d)
  : _dimension(d), _number_of_vertices(0), _data()
{
}

template<class X>
Polytope<X>::Polytope(dimension_type d, size_type nv, const X* data)
  : _dimension(d), _number_of_vertices(nv), _data(data,data+(d+1)*nv)
{
}


template<class X>
Polytope<X>::Polytope(const Matrix<X>& G)
  : _dimension(G.number_of_rows()-1), 
    _number_of_vertices(G.number_of_columns()), 
    _data(G.number_of_rows()*G.number_of_columns())
{
  dimension_type& d=this->_dimension;
  size_type& nv=this->_number_of_vertices;
  X* ptr=this->_data.begin();
  MatrixSlice<X>(d+1,nv,ptr,1u,d+1)=G;
}


template<class X>
Polytope<X>::Polytope(const PointList<X>& pts)
  : _dimension(pts.dimension()),
    _number_of_vertices(pts.size()),
    _data((pts.dimension()+1)*pts.size())
{
  dimension_type& d=this->_dimension;
  size_type& nv=this->_number_of_vertices;
  X* ptr=this->_data.begin();
  MatrixSlice<X> g(d+1,nv,ptr,1,d+1);
  for(size_type j=0; j!=nv; ++j) {
    for(size_type i=0; i!=d; ++i) {
      g(i,j)=pts[j][i];
    }
    g(d,j)=1;
  }
}

template<class X> 
Polytope<X>::Polytope(const Box<R>& bx)
  : _dimension(bx.dimension()),
    _number_of_vertices(bx.number_of_vertices()),
    _data((bx.dimension()+1)*bx.number_of_vertices())
{
  dimension_type& d=this->_dimension;
  size_type& nv=this->_number_of_vertices;
  X* ptr=this->_data.begin();
  MatrixSlice<X> g(d+1,nv,ptr,1,d+1);
  
  size_type j=0;
  for(typename Box<R>::vertices_const_iterator v=bx.vertices_begin();
      v!=bx.vertices_end(); ++v)
  {
    for(size_type i=0; i!=this->dimension(); ++i) {
      this->_generators_()(i,j)=(*v)[i];
    }
    this->_generators_()(d,j)=1;
    ++j;
  }
}



template<class R>
Polytope<R>
polytope(const Box<R>& r)
{
  dimension_type d=r.dimension();
  size_type nv=r.number_of_vertices();
  array<R> data((d+1)*nv);
  size_type j=0;
  for(typename Box<R>::vertices_const_iterator v=r.vertices_begin();
      v!=r.vertices_end(); ++v)
  {
    for(size_type i=0; i!=d; ++i) {
      data[j]=(*v)[i];
      ++j;
    }
    data[j]=1;
    ++j;
  }
  return Polytope<R>(d,nv,data.begin());
}  


template<class X>
Polytope<typename traits<X>::arithmetic_type>
polytope(const Polyhedron<X>& plhd)
{
  typedef typename traits<X>::arithmetic_type F;
  
  dimension_type d=plhd.dimension();
  size_type nc=plhd.number_of_constraints();
  
  const Matrix<X> A=plhd.A();
  const Vector<X> b=plhd.b();
  
  std::vector< Vector<F> > constraints;
  std::vector< Vector<F> > generators;
  
  Vector<F> tmp(d+1);
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
Polytope<R>
approx_polytope(const Polyhedron<R>& plhd)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  
  dimension_type d=plhd.dimension();
  size_type nc=plhd.number_of_constraints();
  
  const Matrix<A> a=plhd.A();
  const Vector<A> b=plhd.b();
  
  std::vector< Vector<A> > constraints;
  std::vector< Vector<A> > generators;
  
  Vector<A> tmp(d+1);
  // Add positivity constraint
  tmp(d)=1;
  constraints.push_back(tmp);
  for(size_type i=0; i!=nc; ++i) {
    for(size_type j=0; j!=d; ++j) {
      tmp(j)=-a(i,j);
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
  array<R> data(nv*(d+1));
  for(size_type j=0; j!=nv; ++j) {
    for(size_type i=0; i!=d+1u; ++i) {
      data[j*(d+1)+i]=R(generators[j](i));
    }
  }
  Polytope<R> pltp(d,nv,data.begin());
  return pltp;
}


template<class X>
dimension_type 
Polytope<X>::dimension() const
{
  return this->_dimension;
}

template<class X>
const Matrix<X>
Polytope<X>::generators() const
{
  return const_cast<Polytope<X>*>(this)->_generators_();
}


template<class X>
size_type 
Polytope<X>::number_of_vertices() const
{
  return this->_number_of_vertices;
}

template<class X>
PointList<X> 
Polytope<X>::vertices() const 
{
  return PointList<X>(this->generators());
}

template<class X>
Point<X>
Polytope<X>::vertex(const size_type& j) const 
{
  Matrix<X> g=this->generators();
  Point<X> result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result[i]=g(i,j);
  }
  return result;
}

template<class X>
typename Polytope<X>::vertices_const_iterator
Polytope<X>::vertices_begin() const 
{
  return PolytopeVerticesIterator<X>(*this,0);
}

template<class X>
typename Polytope<X>::vertices_const_iterator
Polytope<X>::vertices_end() const 
{
  return PolytopeVerticesIterator<X>(*this,this->number_of_vertices());
}


template<class X> 
Box<typename Polytope<X>::R> 
Polytope<X>::bounding_box() const
{
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
Polytope<X>::empty() const
{
  return this->number_of_vertices()==0;
}

template<class X>
tribool 
Polytope<X>::bounded() const
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



template<class X>
tribool 
Polytope<X>::contains(const Point<R>& pt) const
{
  return ::contains(*this,pt);
}

template<class X>
tribool 
contains(const Polytope<X>& ply, const Point<typename Polytope<X>::real_type>& pt) 
{
  return ::contains(ply,pt);
}


template<class X>
tribool 
disjoint(const Polytope<X>& ply, const Box<typename Polytope<X>::real_type>& bx)
{
  return ::disjoint(ply,bx);
}


template<class X>
tribool 
disjoint(const Box<typename Polytope<X>::real_type>& bx, const Polytope<X>& ply)
{ 
  return ::disjoint(ply,bx);
}

template<class X>
tribool 
disjoint(const Polytope<X>& ply1, const Polytope<X>& ply2)
{
  return ::disjoint(ply1,ply2);
}

template<class X>
tribool 
subset(const Polytope<X>& ply1, const Polytope<X>& ply2)
{
  return ::subset(ply1,ply2);
}

template<class X>
tribool 
subset(const Box<typename Polytope<X>::real_type>& bx, const Polytope<X>& ply)
{
  return ::subset(bx,ply);
}

template<class X>
tribool 
subset(const Polytope<X>& ply, const Box<typename Polytope<X>::real_type>& bx)
{
  tribool result=true;
  for(typename Polytope<X>::vertices_const_iterator v=ply.vertices_begin(); v!=ply.vertices_end(); ++v) {
    result = result && bx.contains(*v);
    if(!result) { return result; }
  }
  return result;
}





template<class X>
Polytope<X>
convex_hull(const Polytope<X>& A, const Polytope<X>& B)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
std::string
Polytope<X>::name()
{
  return std::string("Polytope")+"<"+Ariadne::name<X>()+">";
}


template<class X>  
std::ostream& 
Polytope<X>::write(std::ostream& os) const 
{
  return os << "Polytope( vertices=" << this->vertices() << " )";
}

template<class X>  
std::istream& 
Polytope<X>::read(std::istream& is)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class X>
void
Polytope<X>::_instantiate()
{   
  typedef typename traits<X>::number_type R;
  typedef typename traits<X>::arithmetic_type I;

  tribool tb;
  Point<R>* pt=0;
  Box<R>* bx=0;
  Polyhedron<R>* ph=0;
  Polytope<X>* p=0;
  Polyhedron<X>* h=0;
  Polytope<I>* ip=0;
  Polyhedron<I>* ih=0;

  tb=Ariadne::contains(*p,*pt);
  tb=disjoint(*bx,*p);
  disjoint(*p,*bx);
  disjoint(*p,*p);
  subset(*bx,*p);
  subset(*p,*bx);
  subset(*p,*p);
  convex_hull(*p,*p);
  
  *p=polytope(*bx);
  *ip=polytope(*h);
  *ih=polyhedron(*p);
  
  *p=approx_polytope(*ph);
}



} // namespace Ariadne
