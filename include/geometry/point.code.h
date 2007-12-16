/***************************************************************************
 *            point.code.h
 *
 *  Sun Jan 23 18:00:21 2005
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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



#include "point.h"

#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "base/stlio.h"

namespace Ariadne {


template<class X>
Geometry::Point<X>::Point(const std::string& s) : _vector(1)
{
  std::stringstream ss(s);
  ss >> *this;
}


template<class X>
std::ostream& 
Geometry::Point<X>::write(std::ostream& os) const
{
  const Point<X>& pt=*this;
  os << "(";
  if(pt.dimension() > 0) {
    os << pt[0] ;
    for (size_type i=1; i<pt.dimension(); i++) {
      os << "," << pt[i];
    }
  }
  os << ")" ;
  
  return os;
}

template<class X>
std::istream& 
Geometry::Point<X>::read(std::istream& is)
{
  Point<X>& pt=*this;
  static size_type last_size;
  
  std::vector<X> v;
  v.reserve(last_size);
  char c;
  is >> c;
  is.putback(c);
  if(c=='(') {
    Base::read_vector(is, v, '(', ')');
  } else if(c=='[') {
    Base::read_vector(is, v, '[', ']');
  } else {
    throw InvalidInput("Invalid point input");
  }
  last_size = v.size();
  
  pt._vector=LinearAlgebra::Vector<X>(v.size());
  for(size_t i=0; i!=v.size(); ++i) {
    pt._vector(i)=v[i];
  }
  return is;
}


} // namespace Ariadne
