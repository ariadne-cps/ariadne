/***************************************************************************
 *            box.code.h
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

#include "box.h"

#include "combinatoric/binary_word.h" 
#include "base/array.h" 
#include "geometry/point.h" 


namespace Ariadne {

template<class R>
Geometry::Box<R>::Box(const std::string& s)
  : _data()
{
  std::stringstream ss(s);
  ss >> *this;
}



template<class R> 
R
Geometry::Box<R>::volume() const 
{
  R result=1;
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result=mul_approx(result,sub_approx(this->upper_bound(i),this->lower_bound(i)));
  }
  return result;
}

template<class R>
Geometry::Box<R>
Geometry::Box<R>::quadrant(const Combinatoric::BinaryWord& w) const 
{
  ARIADNE_CHECK_DIMENSION(*this,w.size(),"Box Box::quadrant(BinaryWord w)");
  Box<R> quadrant(this->dimension());
  
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
size_type 
Geometry::Box<R>::number_of_vertices() const 
{
  return 1<<this->dimension();
}


template<class R>
Geometry::Point<R> 
Geometry::Box<R>::vertex(size_type i) const 
{
  size_type d=this->dimension();
  state_type result(d); 
  
  ARIADNE_CHECK_VERTEX_INDEX(*this,i,"Point Box::vertex(size_type i)");
  
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
typename Geometry::Box<R>::vertices_const_iterator
Geometry::Box<R>::vertices_begin() const
{
  return BoxVerticesIterator<R>(*this,false);
}

template<class R>
typename Geometry::Box<R>::vertices_const_iterator
Geometry::Box<R>::vertices_end() const
{
  return BoxVerticesIterator<R>(*this,true);
}







template<class R>
std::ostream&
Geometry::operator<<(std::ostream& os, const Box<R>& bx)
{
  if(bx.dimension()==0) {
    os << "EmptyBox";
  }
  else {
    os << "[" << bx.lower_bound(0) << "," << bx.upper_bound(0) << "]";
    for(dimension_type i=1; i!=bx.dimension(); ++i) {
      os << "x[" << bx.lower_bound(i) << "," << bx.upper_bound(i) << "]";
    }
  }
  return os;
}


template<class R>
std::istream& 
Geometry::operator>>(std::istream& is, Box<R>& bx)
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
    bx=Box<R>(v.begin(),v.end());
  }
  else {
    /* representation as lower and upper corners */
    /* FIXME */
    ARIADNE_THROW(InvalidInput,"Box::read(istream&)","");
  }
  return is;
}


} //namespace Ariadne
