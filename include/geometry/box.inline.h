/***************************************************************************
 *            box.inline.h
 *
 *  Copyright 2005-6  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 


#include "macros/throw.h"
#include "base/array.h"
#include "base/iterator.h"
#include "base/tribool.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/box_expression.h"

namespace Ariadne {

  
template<class R> inline
Box<R>::Box(size_type d)
  : _data(2*d)
{ 
  if(d!=0) { this->_data[0]=1; this->_data[1]=0; }
}

template<class R> template<class ForwardIterator> inline
Box<R>::Box(ForwardIterator b, ForwardIterator e)
  : _data(2*std::distance(b,e))
{
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    this->set_lower_bound(i,b->lower());
    this->set_upper_bound(i,b->upper());
    ++b;
  }
}

template<class R> template<class RR> inline
Box<R>::Box(const dimension_type& d, const RR* ptr )
  : _data(2*d)
{
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    this->set_lower_bound(i,ptr[2*i]);
    this->set_upper_bound(i,ptr[2*i+1]);
  }
}

template<class R> template<class RR> inline
Box<R>::Box(const dimension_type& d, const RR ary[][2])
  : _data(2*d)
{
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    this->set_lower_bound(i,ary[i][0]);
    this->set_upper_bound(i,ary[i][1]);
  }
}

template<class R> inline
Box<R>::Box(const array< Interval<R> >& a)
  : _data(2*a.size())
{
  for(dimension_type i=0; i!=a.size(); ++i) {
    this->set_lower_bound(i,a[i].lower());
    this->set_upper_bound(i,a[i].upper());
  }
}

template<class R> inline
Box<R>::Box(const std::vector< Interval<R> >& v)
  : _data(2*v.size())
{
  for(dimension_type i=0; i!=v.size(); ++i) {
    this->set_lower_bound(i,v[i].lower());
    this->set_upper_bound(i,v[i].upper());
  }
}

template<class R> inline
Box<R>::Box(const Point<R>& pt)
  : _data(2*pt.dimension())
{
  for(dimension_type i=0; i!=pt.dimension(); ++i) {
    this->set_lower_bound(i,pt[i]);
    this->set_upper_bound(i,pt[i]);
  }
}

template<class R> inline
Box<R>::Box(const Point< Interval<R> >& pt)
  : _data(2*pt.dimension())
{
  for(dimension_type i=0; i!=pt.dimension(); ++i) {
    this->set_lower_bound(i,pt[i].lower());
    this->set_upper_bound(i,pt[i].upper());
  }
}

template<class R> inline
Box<R>::Box(const Point<R>& pt1, const Point<R>& pt2) 
  : _data(2*pt1.dimension())
{
  if(pt1.dimension()!=pt2.dimension()) {
    ARIADNE_THROW(IncompatibleDimensions,"Box(Point pt1, Point pt2)","pt1=" << pt1 << ", pt2=" << pt2);
  }
  for (size_type i=0; i!=this->dimension(); ++i) {
    this->set_lower_bound(i,min(pt1[i],pt2[i]));
    this->set_upper_bound(i,max(pt1[i],pt2[i]));
  }
}


template<class R> inline
Box<R>::Box(const Vector< Interval<R> >& iv)
  : _data(2*iv.size())
{
  for (size_type i=0; i!=this->dimension(); ++i) {
    this->set_lower_bound(i,iv(i).lower());
    this->set_upper_bound(i,iv(i).upper());
  }
}

template<class R> template<class E> inline
Box<R>::Box(const BoxExpression<E>& original)
  : _data(2*original().dimension())
{         
  const E& expression=original();
  for (size_type i=0; i!=this->dimension(); ++i) {
    this->_data[2*i]=expression.lower_bound(i);
    this->_data[2*i+1]=expression.upper_bound(i);
  }
}

template<class R> inline
Box<R>::Box(const Box<R>& original)
  : _data(original._data)
{ 
}

template<class R> inline
Box<R>& 
Box<R>::operator=(const Box<R>& A) 
{
  if(this != &A) {
    this->_data = A._data;
  }
  return *this;
}

template<class R> template<class E> inline
Box<R>& 
Box<R>::operator=(const BoxExpression<E>& original)
{         
  const E& expression=original();
  this->_data.resize(2*expression.dimension());
  for (size_type i=0; i!=this->dimension(); ++i) {
    this->set_lower_bound(i,expression.lower_bound(i));
    this->set_upper_bound(i,expression.upper_bound(i));
  }
  return *this;
}


template<class R> inline
Box<R> 
Box<R>::empty_box(dimension_type d)
{
  Box<R> r(d);
  r.set_lower_bound(0,+1);
  r.set_upper_bound(0,0);
  return r;
}


template<class R> inline
Box<R> 
Box<R>::unit_box(dimension_type d)
{
  Box<R> r(d);
  for(dimension_type i=0; i!=d; ++i) {
    r.set_lower_bound(i,-1);
    r.set_upper_bound(i,+1);
  }
  return r;
}

template<class R> inline
Box<R> 
Box<R>::positive_orthant(dimension_type d)
{
  Box<R> r(d);
  R inf=Ariadne::inf<R>();
  for(dimension_type i=0; i!=d; ++i) {
    r.set_lower_bound(i,0);
    r.set_upper_bound(i,inf);
  }
  return r;
}

template<class R> inline
Box<R> 
Box<R>::entire_space(dimension_type d)
{
  Box<R> r(d);
  R pinf=inf();
  R minf=-pinf;
  for(dimension_type i=0; i!=d; ++i) {
    r.set_lower_bound(i,minf);
    r.set_upper_bound(i,pinf);
  }
  return r;
}


// Conversion operators
template<class R> inline
Box<R>::operator Point< Interval<R> >() const 
{
  return Point< Interval<R> >(this->dimension(),reinterpret_cast<const Interval<R>*>(this->_data.begin()));
}    


// Comparison operators
template<class R> inline
bool 
Box<R>::operator==(const Box<R>& A) const
{
  if(this->dimension()!=A.dimension()) { return false; }
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    if (this->lower_bound(i)!=A.lower_bound(i)) { return false; }
    if (this->upper_bound(i)!=A.upper_bound(i)) { return false; }
  }
  return true;
}

template<class R> inline
bool 
Box<R>::operator!=(const Box<R>& A) const 
{
  return !(*this == A);
}


// Data access
template<class R> inline
array<R>& 
Box<R>::data()
{
  return this->_data;
}

template<class R> inline
const array<R>& 
Box<R>::data() const
{
  return this->_data;
}



template<class R> inline
const R& 
Box<R>::lower_bound(dimension_type i) const 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::lower_bound(dimension_type i) const","*this=" << *this << ", i=" << i);
  }
  return this->_data[2*i];
}

template<class R> inline
R& 
Box<R>::lower_bound(dimension_type i) 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::lower_bound(dimension_type i)","self=" << *this << ", i=" << i);
  }
  return this->_data[2*i];
}

template<class R> inline
const R& 
Box<R>::upper_bound(dimension_type i) const 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::upper_bound(dimension_type i) const","self=" << *this << ", i=" << i);
  }
  return this->_data[2*i+1];
}

template<class R> inline
R& 
Box<R>::upper_bound(dimension_type i) 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::upper_bound(dimension_type i)","self=" << *this << ", i=" << i);
  }
  return this->_data[2*i+1];
}


template<class R> inline
Interval<R>& 
Box<R>::operator[] (dimension_type i) 
{
  return reinterpret_cast<Interval<R>&>(this->_data[2*i]);
}

template<class R> inline
const Interval<R>& 
Box<R>::operator[] (dimension_type i) const 
{
  return reinterpret_cast<const Interval<R>&>(this->_data[2*i]);
}


template<class R> inline
const Interval<R>& 
Box<R>::interval(dimension_type i) const 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::interval(dimension_type i) const","self=" << *this << ", i=" << i);
  }
  return reinterpret_cast<const Interval<R>&>(this->_data[2*i]);
}


template<class R> inline
Point<R> 
Box<R>::lower_corner() const 
{
  Point<R> result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result[i]=this->lower_bound(i);
  }
  return result;
}

template<class R> inline
Point<R> 
Box<R>::upper_corner() const 
{
  Point<R> result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result[i]=this->upper_bound(i);
  }
  return result;
}


template<class R> inline
Vector< Interval<R> > 
Box<R>::position_vectors() const 
{
  Vector< Interval<R> > result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result(i)=this->interval(i);
  }
  return result;
}


template<class R> inline
tribool 
Box<R>::empty() const
{
  tribool result=false;
  if(this->dimension()==0) {
    return true;
  }
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    if(this->lower_bound(i) > this->upper_bound(i)) {
      return true;
    }
    if(this->lower_bound(i)== this->upper_bound(i)) {
      result=indeterminate;
    }
  }
  return result;
}

template<class R> inline
tribool 
Box<R>::bounded() const
{
  R pos_inf=inf();
  R neg_inf=-pos_inf;
  
  tribool result=true;
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    if(this->lower_bound(i) == neg_inf || this->upper_bound(i)== pos_inf) {
      return false;
    }
  }
  return result;
}

template<class R> inline
const Box<R>&
Box<R>::bounding_box() const
{
  return *this;
}

// Modifying operations
template<class R> inline
void 
Box<R>::clear()
{
  if(this->_data.size()!=0) {
    this->_data[0]=1;
    this->_data[1]=0;
  }
}

template<class R> inline
void Box<R>::set_interval(dimension_type i, Interval<R> x)
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::set_interval(dimension_type i, Interval x) const","self=" << *this << ", i=" << i);
  }
  this->set_lower_bound(i,x.lower());
  this->set_upper_bound(i,x.upper());
}

template<class R> inline
void Box<R>::set_lower_bound(dimension_type i, const R& l) 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::set_lower_bound(dimension_type i, Real l) const","self=" << *this << ", i=" << i);
  }
  this->_data[2*i]=l;
}

template<class R> inline
void Box<R>::set_upper_bound(dimension_type i, const R& u) 
{
  if(i>=this->dimension()) {
    ARIADNE_THROW(InvalidCoordinate,"Box::set_upper_bound(dimension_type i, Real u) const","self=" << *this << ", i=" << i);
  }
  this->_data[2*i+1]=u;
}

template<class R> inline
Box<R> 
Box<R>::neighbourhood(const R& delta) const 
{
  Box<R> result(this->dimension());
  Interval<R> expand(-delta,delta);
  for (size_type j=0; j< this->dimension(); ++j) {
    result[j]=(*this)[j]+expand;
  }
  return result;
}




// Box geometric operations
template<class R> inline
dimension_type 
Box<R>::dimension() const 
{
  return this->_data.size()/2;
}


template<class R> inline
Point<R> 
Box<R>::centre() const
{
  Point<R> result(this->dimension());
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result[i]=this->interval(i).midpoint();
  }
  return result;
}

template<class R> inline
R 
Box<R>::radius() const 
{
  R result=0;
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result=max(this->interval(i).radius(),result);
  }
  return result;
}



template<class R> inline
tribool 
Box<R>::intersects(const Box<R>& r) const
{
  return !this->disjoint(r);
}

template<class R> inline
tribool 
Box<R>::superset(const Box<R>& r) const
{
  return r.subset(*this);
}

template<class R, class X> inline
tribool 
contains(const Box<R>& bx, const Point<X>& pt)
{
  return bx.contains(pt);
}

template<class R> inline
tribool 
disjoint(const Box<R>& r1, const Box<R>& r2)
{
  return r1.disjoint(r2);
}

template<class R> inline
tribool 
intersect(const Box<R>& r1, const Box<R>& r2)
{
  return !r1.disjoint(r2);
}


template<class R> inline
tribool 
subset(const Box<R>& r1, const Box<R>& r2)
{
  return r1.subset(r2);
}

template<class R> inline
tribool 
superset(const Box<R>& r1, const Box<R>& r2)
{
  return r2.subset(r1);
}


template<class R> inline
Box<R> 
bounding_box(const Box<R>& r)
{
  return r.bounding_box();
}


template<class R> inline
ListSet< Box<R> >
split(const Box<R>& bx) 
{
  return bx.split();
}


template<class R> inline
Box<R> 
closed_intersection(const Box<R>& r1, const Box<R>& r2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Box closed_intersection(Box r1, Box r2)");
  Box<R> r3(r1.dimension());
  for(size_type i=0; i != r3.dimension(); ++i) {
    r3[i]=intersection(r1[i],r2[i]);
  }
  return r3;
}

template<class R> inline
Box<R> 
open_intersection(const Box<R>& r1, const Box<R>& r2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Box closed_intersection(Box r1, Box r2)");
  Box<R> r3(r1.dimension());
  for(size_type i=0; i != r3.dimension(); ++i) {
    r3[i]=intersection(r1[i],r2[i]);
    if(r3[i].lower()>=r3[i].upper()) {
      r3[i]=Interval<R>();
    }
  }
  return r3;
}

template<class R> inline
Box<R>
rectangular_hull(const Box<R>& r1, const Box<R>& r2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Box rectangular_hull(Box r1, Box r2)");
  Box<R> r3(r1.dimension());
  for(size_type i=0; i != r3.dimension(); ++i) {
    r3[i]=hull(r1[i],r2[i]);
  }
  return r3;
}

template<class R> inline
Box<R>
rectangular_hull(const Box<R>& r, const Point<R>& pt)
{
  return rectangular_hull(r,Box<R>(pt));
}

template<class R> inline
Box<R>
rectangular_hull(const Point<R>& pt, const Box<R>& r)
{
  return rectangular_hull(Box<R>(pt),r);
}

template<class R> inline
Box<R>
rectangular_hull(const Point<R>& pt1, const Point<R>& pt2)
{
  return Box<R>(pt1,pt2);
}


template<class R1, class R2> inline
Box<typename traits<R1,R2>::arithmetic_type> 
minkowski_sum(const Box<R1>& r1, const Box<R2>& r2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Box minkowski_sum(Box r1, Box r2)");
  Box<typename traits<R1,R2>::arithmetic_type> r3(r1.dimension());
  for(dimension_type i=0; i!=r3.dimension(); ++i) {
    r3.set_lower_bound(i,r1.lower_bound(i)+r2.lower_bound(i));
    r3.set_lower_bound(i,r1.upper_bound(i)+r2.upper_bound(i));
  }
  return r3;
}

template<class R1, class R2> inline
Box<typename traits<R1,R2>::arithmetic_type> 
minkowski_difference(const Box<R1>& r1, const Box<R2>& r2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Box minkowski_difference(Box r1, Box r2)");
  Box<typename traits<R1,R2>::arithmetic_type> r3(r1.dimension());
  for(dimension_type i=0; i!=r3.dimension(); ++i) {
    r3.set_lower_bound(i,r1.lower_bound(i)-r2.lower_bound(i));
    r3.set_upper_bound(i,r1.upper_bound(i)-r2.upper_bound(i));
  }
  return r3;
}



template<class R> inline
Box<R> 
operator+(const Box<R>& r, 
                    const Vector< Interval<R> >& iv)
{
  Box<R> result(r.dimension());
  ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(r,iv,"Vector operator+(Box r, IntervalVector iv)");
  
  for(size_type i=0; i!=result.dimension(); ++i) {
    result.set_interval(i,r[i]+iv(i));
  }
  return result;
}





template<class R> inline
BoxVerticesIterator<R>::BoxVerticesIterator(const Box<R>& r, const bool end)
  : _bx(&r), _i(end==true ? (1<<(r.dimension()-1))*3 : 0), _parity(0), _pt(r.lower_corner()) { }

template<class R> inline
bool BoxVerticesIterator<R>::equal(const BoxVerticesIterator<R>& other) const {
  return this->_i==other._i && this->_bx==other._bx; }

template<class R> inline
const Point<R>& BoxVerticesIterator<R>::dereference() const { 
  return this->_pt; }

template<class R> inline
void BoxVerticesIterator<R>::increment() { 
  uint j=0; uint m=1; if(this->_parity) { while(!(m&(this->_i))) { ++j; m*=2u; } ++j; m*=2u; }
  this->_parity=!this->_parity;
  if(j==this->_bx->dimension()) { this->_i+=m; return; }
  if(m&(this->_i)) { this->_pt[j]=this->_bx->lower_bound(j); this->_i-=m; }
  else { this->_pt[j]=this->_bx->upper_bound(j); this->_i+=m; }
}




} // namespace Ariadne
