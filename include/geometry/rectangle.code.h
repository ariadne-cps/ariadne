/***************************************************************************
 *            rectangle.code.h
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
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
 
#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "rectangle.h"

#include "combinatoric/binary_word.h" 
#include "base/array.h" 
#include "geometry/point.h" 
#include "geometry/point_list.h" 
#include "geometry/list_set.h" 

namespace {

using namespace Ariadne;


Geometry::Rectangle<Numeric::Rational> 
neighbourhood(const Geometry::Rectangle<Numeric::Rational>& r, 
              const Numeric::Rational& delta) 
{
  Geometry::Rectangle<Numeric::Rational> result(r.dimension());
  for (size_type j=0; j!=r.dimension(); ++j) {
    result.set_lower_bound(j,r.lower_bound(j)-delta);
    result.set_upper_bound(j,r.upper_bound(j)+delta);
  }
  return result;
}

template<class T>
Geometry::Rectangle< Numeric::Float<T> > 
neighbourhood(const Geometry::Rectangle< Numeric::Float<T> >&r, 
              const Numeric::Float<T>& delta) 
{
  Geometry::Rectangle< Numeric::Float<T> > result(r.dimension());
  for (size_type j=0; j!=r.dimension(); ++j) {
    result.set_lower_bound(j,sub_down(r.lower_bound(j),delta));
    result.set_upper_bound(j,add_up(r.upper_bound(j),delta));
  }
  return result;
}

} // namespace



namespace Ariadne { namespace Numeric {


}}


namespace Ariadne {


template<class R>
Geometry::Rectangle<R>::Rectangle(const std::string& s)
  : _data()
{
  std::stringstream ss(s);
  ss >> *this;
}



template<class R> Geometry::Rectangle<R> 
Geometry::Rectangle<R>::neighbourhood(const R& delta) const
{
  return ::neighbourhood(*this,delta);
}

template<class R> inline
R
Geometry::Rectangle<R>::volume() const 
{
  R result=1;
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result=mul_approx(result,sub_approx(this->upper_bound(i),this->lower_bound(i)));
  }
  return result;
}

template<class R>
Geometry::Rectangle<R>
Geometry::Rectangle<R>::quadrant(const Combinatoric::BinaryWord& w) const 
{
  ARIADNE_CHECK_DIMENSION(*this,w.size(),"Rectangle Rectangle::quadrant(BinaryWord w)");
  Rectangle<R> quadrant(this->dimension());
  
  for (size_type i=0; i!=this->dimension(); ++i) {
    if(w[i]) {
      quadrant[i]=Numeric::Interval<R>(this->interval(i).midpoint(),this->upper_bound(i));
    } 
    else {
      quadrant[i]=Numeric::Interval<R>(this->lower_bound(i),this->interval(i).midpoint());
    }
  }
  return quadrant;
}

template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Geometry::Rectangle<R>::subdivide() const 
{
  ListSet< Rectangle<R> > result(this->dimension());
  size_type n=this->dimension();
  
  Point<R> lwr_crnr=this->lower_corner();
  Point<R> cntr=this->centre();
  Point<R> upr_crnr=this->upper_corner();
  Point<R> new_lwr_crnr(n);
  Point<R> new_upr_crnr(n);
  for(size_type i=0; i!=1u<<n; ++i) {
    for(size_type j=0; j!=n; ++j) {
      if(i&(1u<<j)) {
        new_lwr_crnr[j]=lwr_crnr[j];
        new_upr_crnr[j]=cntr[j];
      }
      else {
        new_lwr_crnr[j]=cntr[j];
        new_upr_crnr[j]=upr_crnr[j];
      }
    }
    Rectangle<R> new_rect(new_lwr_crnr,new_upr_crnr);
    result.adjoin(new_rect);
  }
  return result;
}


template<class R>
Geometry::PointList<R>
Geometry::Rectangle<R>::vertices() const
{
  size_type number_of_vertices=(1<<this->dimension());
  PointList<R> result(this->dimension(),number_of_vertices);
  Point<R> vertex(this->dimension());
  
  for (size_type i=0; i<number_of_vertices; ++i) {
    for (size_type j=0; j<this->dimension(); ++j) {
      if ((1<<j)&(i)) {
        vertex[j]=this->upper_bound(j);
      } 
      else {
        vertex[j]=this->lower_bound(j);
      }
    }
    result[i]=vertex;
  }
  return result;   
}


template<class R>
size_type 
Geometry::Rectangle<R>::number_of_vertices() const 
{
  return 1<<this->dimension();
}


template<class R>
Geometry::Point<R> 
Geometry::Rectangle<R>::vertex(size_type i) const 
{
  size_type d=this->dimension();
  state_type result(d); 
  
  ARIADNE_CHECK_VERTEX_INDEX(*this,i,"Point Rectangle::vertex(size_type i)");
  
  for (size_type j=0; j<d; ++j) {
    if (i%2) {
      result[j]=this->lower_bound(j);
    }
    else {
      result[j]=this->upper_bound(j);
    }
    i=i/2;
  }
  
  return result;
}


template<class R>
typename Geometry::Rectangle<R>::vertices_const_iterator
Geometry::Rectangle<R>::vertices_begin() const
{
  return RectangleVerticesIterator<R>(*this,false);
}

template<class R>
typename Geometry::Rectangle<R>::vertices_const_iterator
Geometry::Rectangle<R>::vertices_end() const
{
  return RectangleVerticesIterator<R>(*this,true);
}





template<class R> inline
tribool 
Geometry::Rectangle< Numeric::Interval<R> >::empty() const 
{
  tribool result=false;
  if(this->dimension()==0) {
    return true;
  }
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    if(this->lower_bound(i) > this->upper_bound(i)) {
      return true;
    }
    if(possibly(this->lower_bound(i) >= this->upper_bound(i))) {
      result=indeterminate;
    }
  }
  return result;
}


template<class R>
inline
tribool 
Geometry::Rectangle< Numeric::Interval<R> >::contains(const Point< Numeric::Interval<R> >& pt) const 
{
  tribool result=true;
  ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,pt,"tribool Rectangle<Interval>::contains(Point<Interval> pt)");
  const Rectangle< Numeric::Interval<R> >& self=*this;
  for (size_type i=0; i!=self.dimension(); ++i) {
    if(self.lower_bound(i)>pt[i] || pt[i]>self.upper_bound(i)) {
      return false;
    }
    if(self.lower_bound(i)==pt[i] || pt[i]==self.upper_bound(i)) { 
      result=indeterminate;
    }
  }
  return result;
}


template<class R>
Geometry::Point< Numeric::Interval<R> > 
Geometry::Rectangle< Numeric::Interval<R> >::lower_corner() const 
{
  return Point< Numeric::Interval<R> >(this->dimension(),this->_data.begin(),2u);
}


template<class R>
Geometry::Point< Numeric::Interval<R> > 
Geometry::Rectangle< Numeric::Interval<R> >::upper_corner() const 
{
  return Point< Numeric::Interval<R> >(this->dimension(),this->_data.begin()+1u,2u);
}


template<class R>
size_type 
Geometry::Rectangle< Numeric::Interval<R> >::number_of_vertices() const 
{
  return 1<<this->dimension();
}


template<class R>
typename Geometry::Rectangle< Numeric::Interval<R> >::vertices_const_iterator
Geometry::Rectangle< Numeric::Interval<R> >::vertices_begin() const
{
  return RectangleVerticesIterator< Numeric::Interval<R> >(*this,false);
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
typename Geometry::Rectangle< Numeric::Interval<R> >::vertices_const_iterator
Geometry::Rectangle< Numeric::Interval<R> >::vertices_end() const
{
  return RectangleVerticesIterator< Numeric::Interval<R> >(*this,true);
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
std::string
Geometry::Rectangle<R>::name()
{
  return std::string("Rectangle")+"<"+Numeric::name<R>()+">";
}


template<class R>
std::ostream&
Geometry::Rectangle<R>::write(std::ostream& os) const 
{
  const Rectangle<R>& self=*this;
  if(self.dimension()==0) {
    os << "EmptyRectangle";
  }
  else {
    os << "[" << self.lower_bound(0) << "," << self.upper_bound(0) << "]";
    for(dimension_type i=1; i!=self.dimension(); ++i) {
      os << "x[" << self.lower_bound(i) << "," << self.upper_bound(i) << "]";
    }
  }
  return os;
}

template<class R>
std::string
Geometry::Rectangle< Numeric::Interval<R> >::name()
{
  return std::string("Rectangle")+"<"+Numeric::name< Numeric::Interval<R> >()+">";
}


template<class R>
std::ostream&
Geometry::Rectangle< Numeric::Interval<R> >::write(std::ostream& os) const 
{
  const Rectangle< Numeric::Interval<R> >& self=*this;
  if(self.dimension()==0) {
    os << "EmptyRectangle";
  }
  else {
    os << "[" << self.lower_bound(0) << "," << self.upper_bound(0) << "]";
    for(dimension_type i=1; i!=self.dimension(); ++i) {
      os << "x[" << self.lower_bound(i) << "," << self.upper_bound(i) << "]";
    }
  }
  return os;
}

template<class R>
std::istream& 
Geometry::Rectangle<R>::read(std::istream& is)
{
  
  char c;
  is >> c;
  is.putback(c);
  if(c=='[') {
    /* Representation as a literal [a1,b1]x[a2,b2]x...x[an,bn] */
    std::vector< Numeric::Interval<R> > v;
    Numeric::Interval<R> i;
    c='x';
    while(c=='x') {
      is >> i;
      v.push_back(i);
      c=' ';
      while( is && c==' ') {
        is >> c;
      }
    }
    if(is) {
      is.putback(c);
    }
    (*this)=Rectangle<R>(v.begin(),v.end());
  }
  else {
    /* representation as lower and upper corners */
    /* FIXME */
    ARIADNE_THROW(InvalidInput,"Rectangle::read(istream&)","");
  }
  return is;
}


} //namespace Ariadne
