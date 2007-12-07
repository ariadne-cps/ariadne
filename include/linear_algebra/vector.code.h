/***************************************************************************
 *            vector.code.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
#include "vector.h"

#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "base/stlio.h"

#include "numeric/interval.code.h"

namespace Ariadne {

template<class R>
LinearAlgebra::Vector<R>::Vector(const std::string& str)
  : _array(1)
{  
  std::istringstream ss(str); 
  ss >> *this; 
}      

template<class R>
std::ostream&
LinearAlgebra::Vector<R>::write(std::ostream& os) const
{  
  os << "[";
  if(this->size()>0) {
    os << (*this)(0);
    for(uint i=1; i!=this->size(); ++i) {
      os << "," << (*this)(i);
    }
  }
  os << "]";
  return os;
}      

template<class R>
std::istream&
LinearAlgebra::Vector<R>::read(std::istream& is)
{  
  std::vector<R> stdvec;
  is >> stdvec;
  *this=Vector<R>(stdvec.size());
  for(size_type i=0; i!=this->size(); ++i) {
    (*this)(i)=stdvec[i];
  }
  return is;
}      


}
