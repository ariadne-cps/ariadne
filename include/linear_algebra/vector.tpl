/***************************************************************************
 *            vector.tpl
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

#include "../utility/stlio.h"

namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    Vector<R>::Vector(const std::string& s)
      : _Base(1)
    {  
      std::istringstream ss(s); 
      ss >> *this; 
    }      

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Vector<R>& v)
    {
      os << "[";
      if(v.size()>0) {
        os << v[0];
      }
      for(uint i=1; i!=v.size(); ++i) {
        os << "," << v[i];
      }
      os << "]";
      return os;
    }
    
    template <typename R>
    std::istream&
    operator>>(std::istream& is, Vector<R>& v)
    {
      std::vector<R> stdvec;
      is >> stdvec;
      v=Vector<R>(stdvec.size());
      for(size_type i=0; i!=v.size(); ++i) {
        v(i)=stdvec[i];
      }
      return is;
    }

  }
}
