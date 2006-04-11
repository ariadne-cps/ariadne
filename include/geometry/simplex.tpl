/***************************************************************************
 *            simplex.tpl
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include "simplex.h"

#include <iostream>
#include <vector>

#include "../utility/stlio.h"
#include "../geometry/parallelotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    Simplex<R>::Simplex(size_type n)
      : _vertices(n+1,state_type(n)) 
    {
      for(size_type i=0; i!=n; ++i) {
        _vertices[i][i]=1;
      }
    }
  
    template<typename R>
    Simplex<R>::Simplex(const std::vector<state_type>& v)
      : _vertices(v.begin(),v.end())
    {
      size_type d=_vertices.size()-1;
      for(size_type i=0; i!=d+1; ++i) {
        if(_vertices[i].dimension()!=d) {
          throw std::domain_error("The the list of vertices is invalid");
        }
      }
    }
    
    template<typename R>
    Simplex<R>::Simplex(const array<state_type>& v)
      : _vertices(v.begin(),v.end())
    {
      size_type d=_vertices.size()-1;
      for(size_type i=0; i!=d+1; ++i) {
        if(_vertices[i].dimension()!=d) {
          throw std::domain_error("The the list of vertices is invalid");
        }
      }
    }
    
    template<typename R>
    Simplex<R>::Simplex(const std::string& s)
      : _vertices()
    {
      std::stringstream ss(s);
      ss >> *this;
    }
    
    template<typename R>
    Simplex<R>::operator Polyhedron<R> () const 
    {
      std::vector<state_type> vert_vec(_vertices.begin(),_vertices.end());
      return Polyhedron<R>(vert_vec);
    }
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Simplex<R>& s) 
    {
//      if(s.empty()) {
//        os << "Empty";
//     }
//      else 
      if(s.dimension() > 0) {
        os << "Simplex( vertices=";
        Utility::write_sequence(os,s.vertices().begin(),s.vertices().end());
        os << ")";
      }
      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Simplex<R>& s)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}
